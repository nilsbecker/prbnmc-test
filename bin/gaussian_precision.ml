(* evidence calculation based on MH steps and stepwise tightening *)

(* naive runtime *)
let t_start = Sys.time ()

open Containers
open Dagger
open Evidence

let () as _init_random =
  rng_state := RNG.make [| 0x13374; 0x533D99 |]

(* parameter object. keep this out of evidence lib which is supposed to be
   compatible with different models and cases. *)

let params = object (self)
  method data_mean = 3.
  method data_std = 1.
  method data_N = 100
  method data_2d_corr = 0.9
  method data_2d_angle_rad = Gsl.Math.pi /. 4.
  method prior_alpha = 1
  method prior_beta = 2.
  method prior_kernel_step = 0.01
  method kernel_walker_steps = 3
  method model_N = self#data_N
end


(* gaussian walk kernel and helper functions *)

(*
start with a simplified version of the scheme above where data are 1d, and the
model is a 1d gaussian where mean μ is known and the precitions τ=1/std is
to be found. this has a known conjugate prior for the precision, which is a
Gamma distribution
*)

(* TODO bug report: Gamma dist should take a non-integer shape *)

let prior_precision ~alpha ~beta =
  Gen.run (Gen.gamma ~shape:alpha ~scale:(1./. beta))

(* first, make a kernel that generates prior_precision as it's stationary
   distribution. to do that we need to provide
   1. an initial condition. here we take the mean of the prior.
   2. a sampler. to do that we need a step width and some way to accept/reject.
   3. a log likelihood ratio proposed / current; this we can just take from the
   explicit prior *)

(** prior_kernel ~alpha ~beta ~step is a kernel that generates proposals that
    converge to the gamma distribution; the step width is the step in
    cumulative-probability space used by the transformation method; should be
    much less than 1. *)
let prior_kernel ~start ~alpha ~beta ~step =
  let _ = assert (step <. 1.) in
  let dist = Pdfs.gamma_ln ~shape:alpha ~scale:(1./.beta) in
  let a, b = float alpha, beta in
  let sampler tau_current rng_state =
    (* FIXME: in case there are numerical problems, switch to the complementary
       cdf close to 1 and add float error handling *)
    let cdf x = Gsl.Cdf.gamma_P ~x ~a ~b in
    let inv_cdf p = Gsl.Cdf.gamma_Pinv ~p ~a ~b in
    let p_current = cdf tau_current in
    let p_proposed =
      p_current +. Gen.run (Gen.range {min=(-.step); max=step}) rng_state in
    let p_proposed =
      if (0. <. p_proposed) && (p_proposed <. 1.) then p_proposed
      else p_current in
    let tau_proposed = inv_cdf p_proposed in
    tau_proposed in
  let log_pdf cur prop = Log_space.unsafe_cast @@ dist prop -. dist cur in
  Dist.kernel start sampler log_pdf

(* model definition *)

module Types = Make_thresholding_types (struct type param=float end)

module T = Thresholding (Types)

module Model = struct
  module Smc = T.Smc

  (* we describe one walker as it proceeds from one reweighting to the next. the
     driver should become a stop criterion which looks at the current population
     and the current distance thresholds and decides to stop. what's missing is
     the rejection for escaping the distance threshold. unclear where we should
     add this: at each individual walker step? also, the current distance
     threshold is data that we need to get into the algorithm, but which comes
     only from the last reweighting.

     in other words we need to pass around resampling_state! this is returned by
     yield (!). *)

  open Types

  let gaussian_distance params data tau =
    (* sample synthetic data. *)
    let synthetic_gaussian_data params ~tau =
      let dist = Gen.gaussian ~mean:params#data_mean ~std:(1. /. tau) in
      Array.init params#model_N (fun _ -> Gen.run dist !rng_state) in
    let synthetic_data = synthetic_gaussian_data params ~tau in
    sq_dist data synthetic_data

  let model
      stopping      (* stopping criterion, evaluated on the resampling state *)
      threshold_rejection
      (* threshold rejection for exceeding previous distance level from data *)
    =
    let kernel start =
      (* the kernel for this walker; carries current tau as mutable state! *)
      prior_kernel
        ~start
        ~alpha:params#prior_alpha
        ~beta:params#prior_beta
        ~step:params#prior_kernel_step in
    let walk_to_new_tau state current_output =
      (* generate a new tau proposal in the vicinity and with recording of weights
         in each step; this should effectively be one compound step in the prior.
         arguably, we could just adjust the step size if we do nothing in between *)
      let open Lmh_inference in
      let walking_model = sample (kernel current_output.param) in
      let proposed_tau = stream_samples walking_model !rng_state
      (* this is not strictly correct. would need to check threshold after each
         step in order to really implement update with the kernel constrained
         to the distance threshold. *)
                         |> Seq.take params#kernel_walker_steps
                         |> Seq.fold (fun _ a -> a) 0. in
      match threshold_rejection state proposed_tau with
      | (Reject, _distance) ->                                     (* rejected *)
        ( (* Format.printf "staying  at tau=%f@." current_tau;*)
          {current_output with acceptance=Reject})
      | (Accept, distance) ->
        ( (* Format.printf "accepted to tau=%f@." proposed_tau; *)
          {acceptance=Accept; distance; param=proposed_tau}) in

    let rec loop =
      let open Smc in let open Infix in
      fun output ->
        let* state = yield output in
        if stopping state then return output.param else
          loop (walk_to_new_tau state output) in
    loop
      (* start at indep. draw from prior for better mixing (cheating) *)
      {acceptance=Accept;
       distance=Float.nan;
       param=(prior_precision
               ~alpha:(float params#prior_alpha) ~beta:params#prior_beta)
           !rng_state}

end


(* main *)

let _1d_gaussian_test = Format.printf "gaussian samples with μ=%.2f and σ=%.2f:@\n%a@."
    params#data_mean params#data_std
    (Seq.pp Float.pp) (data_gaussian_1d params)

let _2d_gaussian_test =
let f_pp fmt f = Format.fprintf fmt "%.2f" f in
  let pp_sep fmt () =  Format.pp_print_string fmt ":" in
  let a_pp fmt a = (Array.pp ~pp_sep f_pp) fmt a  in
  Format.printf
    "2d gaussian samples with equal μ=%.2f, σ=%.2f and corr coeff ρ = %.2f@\n%a@."
    params#data_mean params#data_std params#data_2d_corr
    (Seq.pp a_pp) (data_gaussian_2d params)

let _gamma_sampling_kernel =
  if true then () else  (* disabled. manual way; see below for idiomatic way *)
  let n, alpha, beta, step = 100_000, 2, 3., 0.1 in
  let start = let mean = float alpha /. beta in mean in
  let samples =
    match prior_kernel ~start ~alpha ~beta ~step with
    | Dist.Stateless _ ->
      failwith "bad idea to pack stateless and kernel like that"
    | Dist.Kernel k -> (* mc by hand using private type *)
      k.start
      |> Seq.iterate (fun state -> k.sample state !rng_state)
      |> Seq.take n in
  Format.printf "gamma prior samples with α=%d, β=%.2f, η=%.*f:@\n%a@."
    alpha beta 4 step
    (Seq.pp Float.pp) (Seq.take 50 samples);
  (* tested: this seems to yield the correct pdf! nice. *)
  Csv.save "gamma_samples.csv"
    (samples |> Seq.map Float.to_string |> Seq.to_list |> List.pure)

let _gamma_sampling_kernel =
  let n, alpha, beta, step = 1_000_000, 2, 3., 0.1 in
  let start = let mean = float alpha /. beta in mean in
  let kernel = prior_kernel ~start ~alpha ~beta ~step in
  (* automatic metropolis-hastings sampling, where no actual inference
     is done. checked, works! *)
  let model : float Lmh_inference.t = Lmh_inference.sample kernel in
  let samples =
    (* the stream_samples function generates new proposals from the transition
       kernel, which we internally implemented by the inversion method. it also
       tracks the log-likelihood of the proposed steps. in an actual model with
       a longer evaluation trace, an elementary resampling would be introduced
       at a random position, and the entire model reevaluated to the end; then
       a MH step would take a step in the space of exectution traces. this is
       not what we do here: the model is just a single sample draw! note that
       the current state of the kernel is correctly propagated to the next
       draw. *)
    Lmh_inference.stream_samples model !rng_state
    |> Seq.take n in
  Format.printf "lmh gamma prior samples with α=%d, β=%.2f, η=%.*f:@\n%a@."
    alpha beta 4 step
    (Seq.pp Float.pp) (Seq.take 50 samples);
  (* tested: this seems to yield the correct pdf! nice. *)
  samples
  |> Seq.map Float.to_string
  |> Seq.to_list
  |> List.return
  |> Csv.save "gamma_samples.csv"


(* ok, sampling works. not instantiate a model
   and test evidence calculation! *)

module M = Model

(* comparison to data *)
let the_data = params |> data_gaussian_1d |> Seq.to_array
let compute_distance tau = M.gaussian_distance params the_data tau


let initial_state =
  { Types.distance_threshold=max_float; iteration=0;
    log_evidence=Log_space.one; mobile_fraction=1.}

(* kludgy *)
let state_tracker = Some (ref initial_state)

let threshold_resampling =
  T.evidence_resampling
    ~cutoff_fraction:0.3
    ~state_tracker

(* desired epsilon *)
let final_distance_threshold = 0.01

(* model definition including stopping criterion *)
let model =
  let open Types in
  let stopping state =
    state.distance_threshold <. final_distance_threshold
    || state.iteration > 125
  in
  let threshold_rejection state proposed_tau =
    let distance = compute_distance proposed_tau in
    let acceptance = if distance >. state.distance_threshold
      then Reject else Accept in
    acceptance, distance in
  M.model stopping threshold_rejection

let sequence_of_populations =
  T.Smc.run ~nthreads:8 ~npart:1_000
    threshold_resampling
    initial_state
    model
    !rng_state

(* problem: we don't get to keep the evidence; the population state does not
   contain it. as the population state is inherited from Smc_inference.Make, we
   cannot easily add it either. *)

(* only here do we calculate (!) *)
let pops = sequence_of_populations |> Seq.to_list

let _l =
  pops |> List.map (fun pop -> pop.T.Smc.active |> Array.length)

let _ = Option.iter (fun tr ->
    Format.printf "@,total evidence after reaching threshold sq_dist %g: %g@."
      final_distance_threshold (Log_space.to_float !tr.Types.log_evidence))
    state_tracker

(*let _ = Format.printf "@.active: %a@." (List.pp Int.pp) l*)

let _ =
  match List.rev pops with
  | [] -> ()
  | h :: _ ->
    h.T.Smc.terminated
    |> Array.map (fun (tau, _) -> [Float.to_string tau])
    |> Array.to_list
    |> Csv.save "posterior_tau.csv"

let _ =
  Format.printf "@\ntotal time: %.2fs@." (Sys.time () -. t_start);
  exit 0
