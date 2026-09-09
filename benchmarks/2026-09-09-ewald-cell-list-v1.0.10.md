# Ewald real-space cell-list parity — 2026-09-09

This local gate verifies the bounded implementation slice for the PME
real-space term. `spme_energy` now calls
`direct_coulomb_damped_cutoff`, which uses a uniform cell list with cell edge
`r_cut` and visits only the 27 neighboring cells. The strict `r < r_cut`
boundary and `r > 1e-6` overlap guard are retained.

The implementation is checked against an independent nested-loop reference on
fixtures containing negative coordinates, coordinates on cell boundaries, a
pair exactly at the cutoff, and pairs outside the cutoff. The package test
command passed:

```text
CARGO_TARGET_DIR=/private/tmp/chematic-target-ewald-cell \
  cargo test -p chematic-ewald --offline --lib
```

This is correctness and boundary evidence, not a throughput claim. The
periodic `chematic-crystal` neighbor enumerator still uses its validated
bounded all-pairs/image path and remains open for a separate cell-list design
and numerical gate.
