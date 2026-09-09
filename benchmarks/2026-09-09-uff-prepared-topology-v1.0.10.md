# UFF prepared-energy topology — 2026-09-09

The UFF minimizer now prepares its bond, angle, and graph-based van der Waals
exclusion terms once per minimization. Each energy and analytic-gradient
evaluation then reuses those immutable terms instead of rebuilding the
molecule traversal and parameter combinations. The public `uff_total_energy`
API remains a compatibility wrapper.

Validation:

```text
CARGO_TARGET_DIR=/private/tmp/chematic-target-uff-prepared \
  cargo test -p chematic-ff --offline --lib
```

The command passed the current `chematic-ff` library suite, including
analytic-gradient comparisons with finite-difference references. A separate
dispatch probe covers 18 public metal/halogen element mappings, while the
prepared UFF parameter fixture covers all 45 declared `UffType` variants.
This is a topology and gradient correctness gate, not a speed claim. The
remaining independent force-field oracle parity, convergence, stereo, and
broader quality gates remain open roadmap work.
