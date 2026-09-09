# 3D force-field energy sanity: v1.0.10 candidate

The bounded `pipeline_v2` integration corpus was rerun on 2026-09-09 with
release optimizations. It contains 63 molecules: the frozen 58-molecule corpus
plus five specification stress cases.

```text
energy_sanity arm=E_mmff94_strict attempted=63 successful=61 finite=61 non_increasing=61 tolerance=1e-6
energy_sanity arm=F_mmff94_widened attempted=63 successful=61 finite=61 non_increasing=61 tolerance=1e-6
energy_sanity arm=G_mmff94_uff_fb attempted=63 successful=61 finite=61 non_increasing=61 tolerance=1e-6
energy_sanity arm=H_dreiding attempted=63 successful=61 finite=61 non_increasing=61 tolerance=1e-6
```

All successful force-field outcomes reported finite `before` and `after`
energies, and every accepted result was non-increasing within `1e-6` force-field
energy units. The four force-field arms each had 61 successful outcomes; the
two typed failures in each arm remain part of the denominator and are not
silently counted as passes.

Reproduce with:

```sh
cargo run --release -p chematic-3d --example pipeline_v2_integration_gate
```

The machine-readable record is
[`2026-09-09-3d-energy-sanity-v1.0.10.json`](2026-09-09-3d-energy-sanity-v1.0.10.json).
This is a local soundness slice, not an independent MMFF94/UFF oracle or a
completion claim for typing, parameters, convergence, stereo, analytic
gradients, or the remaining symmetry-aware 3D quality gates.
