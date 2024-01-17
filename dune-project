(lang dune 3.8)

(name timecount)

(generate_opam_files true)

(authors "Nils")

(maintainers "Nils")

(package
 (name timecount)
 (allow_empty)
 (synopsis "Test of probabilistic progamming")
 (description "test")
 (depends ocaml dune prbncm-dagger containers)
 (tags
  (topics "probabilistic programming")))

; See the complete stanza docs at https://dune.readthedocs.io/en/stable/dune-files.html#dune-project
