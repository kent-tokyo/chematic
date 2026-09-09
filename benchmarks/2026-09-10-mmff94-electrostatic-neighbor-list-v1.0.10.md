# v1.0.10 MMFF94 electrostatic neighbor-list slice

The prepared MMFF94 model now exposes an opt-in 10 Å coordinate-dependent
neighbor-list path for electrostatic energy and analytic gradient evaluation.
The compatibility-default evaluator remains all-pair for electrostatics.

The implementation uses the same strict `distance <= 10.0 Å` boundary as the
prepared vdW list, handles negative coordinates, and falls back to the
prepared all-pair list for non-finite coordinates. Tests verify the candidate
set and energy against an independent all-pair cutoff calculation, then verify
the combined cutoff gradient against central differences on a nine-atom
ethanol fixture. The prepared analytic objective also has a separate butane
finite-difference gate, and the opt-in L-BFGS entry point switches the cutoff
energy and gradient together.

Reproduce:

```sh
cargo test -p chematic-ff --offline --lib \
  coordinate_dependent_electrostatic_neighbor_list_matches_cutoff_reference \
  --quiet
cargo test -p chematic-ff --offline --lib \
  cutoff_nonbonded_gradient_matches_cutoff_energy_finite_difference \
  --quiet
cargo test -p chematic-ff --offline --lib \
  cutoff_bounded_analytic_gradient_matches_cutoff_total_energy \
  --quiet
```

Machine-readable evidence: [`mmff94-electrostatic-neighbor-list-v1.0.10.json`](../validation/results/mmff94-electrostatic-neighbor-list-v1.0.10.json).

This remains a bounded opt-in slice. It does not claim that 10 Å is the
compatibility-default MMFF94 electrostatic cutoff, production-default
minimizer behavior, or periodic electrostatics correctness.
