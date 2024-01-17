(* library of functiones needed for the evidence calculation based on MH steps
   and stepwise tightening *)

open Dagger
module Gen = Stats.Gen.Make (RNG)
module Pdfs = Stats.Pdfs
open Containers

module Gsl_dist = Gsl_dist.Make (Gsl)


let rng_state = ref (RNG.make [| 0x1337; 0x533D |])

(* drawing d-dimensional gaussian random variables *)

(** attention: inplace; the input is modified *)
let eigendecomp m_sym =
  let open Gsl.Eigen in
  symmv ~protect:false (`M m_sym)
  |> Fun.tap (Fun.flip symmv_sort VAL_DESC)

let vmap f v =
  let dim = Gsl.Vector.length v in
  Bigarray.(Array1.init float64 c_layout dim
              (fun i -> Gsl.Vector.get v i |> f))

let ev_mvnd ~eigenbasis ~magnitudes ~mean =
  let open Gsl in
  (* manual protect: *)
  let indepsamples =
    magnitudes |> vmap (fun sigma ->
        Randist.gaussian (Gsl_dist.rng !rng_state) ~sigma) in
  let y = Vector.copy mean in
  Blas.(gemv NoTrans ~alpha:1. ~a:eigenbasis ~x:indepsamples ~beta:1. ~y); y

let gen_mvnd ~mean ~covar =
  (* manual protect: *)
  let m = Gsl.Matrix.copy covar in
  let eigenvalues, eigenbasis = eigendecomp m in
  let magnitudes = eigenvalues |> vmap (fun f -> sqrt (abs_float f)) in
  fun () -> ev_mvnd ~eigenbasis ~magnitudes ~mean


(* fake data on demand. params are not read yet! *)

let data_gaussian_1d params =
  Gen.run
    Gen.(iid (gaussian ~mean:params#data_mean ~std:params#data_std)) !rng_state
  |> Seq.take params#data_N

let data_gaussian_2d params =
  let mean = Gsl.Vector.create ~init:params#data_mean 2 in
  let covar =
    let s = params#data_std ** 2. in
    let rho = params#data_2d_corr in
    Gsl.Matrix.of_arrays [|[|s; s *. rho|]; [|s *. rho; s|]|] in
  let sample = gen_mvnd ~mean ~covar in
  Seq.of_dispenser (fun () -> Some (sample ()))
  |> Seq.map (Gsl.Vector.to_array)
  |> Seq.take params#data_N


(* now we have samples (abbreviated d). we want to get an evidence calculation
   in the step wise scheme. what do we need?

   1. a prior distribution for μ,ρ and a metropolis sampler from it.

   2. a generative model. we choose a model that yields 2d samples {x,y} from a 2d Gaussian
   with 2 free parameters: the correlation coefficient ρ and the common mean μ --
   as an example.

   3. a collection of particles (walkers) which hold a current value for ρ,μ
   and produce a set of generated points (abbreviated g) upon request,
   initialized from the prior

   4. summary statistics to judge the accuracy of a parametrized model. given a
   set of generated samples, and the set of data samples, we compute the
   empirical means m_{g,d} and empirical correlation coefficients r_{g,d}. then
   we calculate the (squared) distance d = (m_g - m_d)^2/m_d^2 + (r_g - r_d)^2

   5. MC loop:
    a) start with a generous distance allowance.
    b) let all walkers generate points and throw out the bottom fraction (1-X)
   with the worst distance
    c) resample from the remaining ones and perform an independent proposal
   change in ρ,μ for each. we need an update that preserves the prior
   distribution here! e.g. symmetric proposal and metropolis acceptance with
   prior weights; also, constrain the dynamics to the attained value of d by
   rejecting steps with higher d unconditionally.
    d) repeat the loop from b) until a stop criterion is met. e.g. fraction of
   rejected moves to forbidden d too high, desired d reached, etc

   6. during the loop, keep track of the fraction φ of walkers that were
   retained at step b. φ is always around X. by throwing out the other fraction
   1-X, we are conditioning on distance below the threshold distance which
   yields X. we can renormalize the weights of the retained fraction to sum to
   1 (or collection size) and keep φ on record.

   6. tunable parameters: jump distance protocol for the parameter update, discarded
   fraction X, number of generated samples G, collection size M

*)

(* to score the model fit we need the distance. we make one from mean and
     variance, normalizing to the data for comparability. *)
let sq_dist data pred =
  let m_data, m_pred, v_data, v_pred =
    Stats.Emp.Float.(
      empirical_mean data, empirical_mean pred,
      empirical_variance data, empirical_variance pred) in
  ((m_pred -. m_data) /. m_data) ** 2.
  +. ((v_pred -. v_data) /. v_data) ** 2.


(* TODO: types do not feel right. we will want to generalize to other particle
 * outputs. basically the types here are suitable for thresholding but the
 * particle output is too restrictive. *)

module type Thresholding_Param_Space = sig type param end

module Make_thresholding_types (Param:Thresholding_Param_Space) = struct
  type distance = float
  type acceptance = Accept | Reject
  type param = Param.param
  type particle_output =
    {acceptance:acceptance; distance:distance; param:param}
  type resampling_state =
    {distance_threshold:distance;                       (* threshold reached *)
     iteration:int;                                  (* resampling iteration *)
     log_evidence:Log_space.t;             (* reduction in log weight so far *)
     mobile_fraction:float;            (* % accepted moves inside the region *)
    }                    (* possibly more data needed for stopping criterion *)
end

module type Thresholding_types =                         (* abstract version *)
  module type of (Make_thresholding_types (struct type param end))

module Thresholding (Types:Thresholding_types) = struct
  module Smc = Smc_inference.Make (Types)

  (* basic insights:
     - internally, Smc has a particles module (call it P) which operates on an
       internal population data type with a double buffer
     - P.iter iterates on the front buffer, as does P.fold
     - P.append appends to the back buffer.
     - after doing resampling, the buffers are swapped by Smc.
     - thus we can use P.append below as adding to an initially empty bag of
       particles for the next round. *)

  (*   there is confusion as to the correct type for the weights
       to use. in the Smc_inference, the presupplied resampling strategies use
       float (!) although the particles internally store weights as log. this
       surely is in need of improvement. for now stick with whatever Smc does
       to avoid pitfalls. this is to just use float, but somehow internally the
       particle weights are still log_space... messy. *)

  type resampling_strategy =
    (Types.particle_output, float, Types.resampling_state) Resampling.strategy

  let evidence_resampling
      ~state_tracker    (* ugly: a reference to the running resampling state *)
      ~cutoff_fraction     (* fraction of particles dropped at high distance *)
    : resampling_strategy = fun
    ~target_size                                  (* desired population size *)
    (module P)
    (* P contains the mutable population as well as the operations on it. it
       is analyzed and updated during the resampling step. *)
    state                   (* necessary for resampling and propagated along *)
    rng_state ->
    let retained_pop =
      let a = Vector.create () in
      P.iter (fun p w -> Vector.push a (w, p));
      Vector.sort'
        (fun (w, p) (w', p') ->
           let o, o' =
             let disto p =
               P.get_output p |> Option.map (fun o -> o.Types.distance) in
             disto p, disto p' in
           Pair.compare (Option.compare Float.compare) Float.compare
             (o, w) (o', w')) a; a in
    (* truncate population to lower effective distance threshold *)
    let retained = truncate
        (float (Vector.length retained_pop) *. (1. -. cutoff_fraction)) in
    let () = Vector.truncate retained_pop retained in
    let new_threshold =
      let (_, p) = Vector.get retained_pop (pred retained) in
      (Option.get_exn_or "internal error: cannot get particle output"
         (P.get_output p)).distance in
    Format.printf "@\niteration: %d@,threshold: old %g, new %g@."
      state.Types.iteration state.distance_threshold new_threshold;
    (* the new threshold can get bigger than the old one due to random
       fluctuations in the distance estimate!  *)
    let retained_mass =            (* FIXME: susceptible to underflow etc. *)
      Vector.fold (fun w (w', _) -> w +. w') 0. retained_pop in
    let retained_mass_fraction = retained_mass /. P.total () in
    Format.printf "retained length mass, mass fraction: %d, %f, %f@."
      retained retained_mass retained_mass_fraction;
    let log_evidence =
      Log_space.(mul state.log_evidence (of_float retained_mass_fraction)) in
    Format.printf "accumulated evidence: %g@." (Log_space.to_float log_evidence);
    let mobile_fraction =
      let accepted = Vector.fold (fun acc (_, p) ->
          match P.get_output p with
          | Some {Types.acceptance=Accept; _ } -> succ acc | Some _ | None -> acc)
          0 retained_pop in
      float accepted /. float retained in
    Format.printf "mobile fraction %f@." mobile_fraction;
    (* now we have a new threshold and reduced number of particles. next we
       need to resample them. ie. go through the truncated population vector
       and append to the population for the next iteration.

       copy-and-paste from Resampling.Make because otherwise not usable with
       the truncated pop *)
    let resample_generic_from_vec pop f =
      let cumulative = ref 0. in
      let partition_index = ref 1 in
      let last = ref (f !partition_index rng_state) in
      let res_w = retained_mass /. (float_of_int retained) in
      Vector.iter (fun (w, particle) ->
          cumulative := !cumulative +. w ;
          while (!last <. !cumulative) do
            P.append particle res_w ;
            last := f !partition_index rng_state ;
            incr partition_index done) pop in
    let _resample_stratified =
      if retained < 2 then invalid_arg "< 2 retained particles";
      if target_size < 2 then invalid_arg "target_size < 2";
      let inv = retained_mass /. float_of_int target_size in
      resample_generic_from_vec retained_pop
        (fun i rng_state ->
           retained_mass *. float_of_int i /. float_of_int target_size
           +. RNG.float rng_state inv) in
    let new_state = { Types.distance_threshold=new_threshold;
                      iteration=succ state.iteration;
                      log_evidence;
                      mobile_fraction } in
    let () = Option.iter ((Fun.flip (:=)) new_state) state_tracker in
    new_state


  (* now we have to implement mc updates in the new population which respect
     the new distance threshold. this should be the model. can we drive all of
     this using Smc.run? *)

  (* answer, citing from the docs:

     In a nutshell:
     particles give their current value as argument to yield, and the outcome
     of resampling (of type resampling_state) is broadcast to each particle as
     the value to which yield evaluates to. Resampling steps take the outcome
     of the previous resampling step as argument.
  *)

end
