# MMFF94 prepared nonbonded gradient — 2026-09-09

`Mmff94EnergyModel` now computes analytic gradients for prepared bond, angle,
buffered 14-7 van der Waals, and buffered Coulomb terms. These remain
separately verified building blocks while the existing L-BFGS finite-difference
path remains unchanged until a full minimization soundness gate is available.

Validation:

```text
CARGO_TARGET_DIR=/private/tmp/chematic-target-mmff94-gradient \
  cargo test -p chematic-ff mmff94_minimizer --offline --lib --quiet
CARGO_TARGET_DIR=/private/tmp/chematic-target-mmff94-gradient \
  cargo clippy -p chematic-ff --offline --all-targets -- -D warnings
```

The focused suite passed 47 tests and clippy passed. The bond/angle test uses
methane coordinates and the nonbonded test uses propane coordinates; both use
central differences at `delta = 1e-5 Å`. This is a bounded correctness slice,
not a full MMFF94 analytic-gradient or performance claim.
