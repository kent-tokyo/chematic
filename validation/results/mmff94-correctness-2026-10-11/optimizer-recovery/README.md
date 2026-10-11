# MMFF94 optimizer recovery and exact coordinate JSON, 2026-10-11

This uncommitted CheMatic 1.0.20 source candidate fixes the L-BFGS
recovery step and the final convergence report. It also enables exact
binary64 JSON decoding after an external stereo check exposed coordinate
rounding at an E/Z boundary. No public API signature or force-field
energy/gradient kernel changed in this phase.

## Implementation

- Both the quasi-Newton line search and its steepest-descent recovery
  require finite coordinates/energy, Armijo sufficient decrease and the
  caller's geometry predicate. The previous fixed 0.01 angstrom recovery
  displacement could increase energy near a stationary point.
- A rejected or invalid quasi-Newton direction clears curvature history.
  Recovery restarts with a normalized displacement capped at 0.01 angstrom
  and backtracks up to 40 trials. If neither search finds an admissible
  step, the optimizer retains the last accepted geometry.
- Non-finite gradient components prevent convergence. The convergence
  flag is checked on the final accepted state, including the last allowed
  iteration and a zero-iteration budget.
- The 3D crate enables serde_json/float_roundtrip. The external common
  scorer had decoded some Python JSON numbers to an adjacent binary64
  value. For input 191, this changed the declared alkene angle from
  -90.0 to -89.99999999999997 degrees and produced a false E/Z violation.
  The production stereo verifier and its boundary are unchanged. The
  regression checks coordinate bits through both typed decoding and
  serde_json::Value, then verifies the declared alkene.

The [incremental patch](change.patch) is relative to the preceding
[cation/water source candidate](../cation-water/README.md), whose cumulative
patch and recorded HEAD define the previous optimizer. It contains only
this phase's implementation/configuration and tests; it is not a snapshot
of unrelated concurrent work.

## Fixed-coordinate MMFF94 validation

The original 265-input manifests, explicit-H reference coordinates,
protocol and RDKit 2026.03.6 oracle are unchanged. All candidate JSONL rows
are identical to the preceding candidate: 262 comparable, two RDKit
embedding failures and one declared unsupported. Maximum total-energy
difference remains 9.094947017729282e-13 kcal/mol.

The direct analytic-gradient comparison covers all 37,407 Cartesian
components in the same 262 molecules. Maximum absolute component error
is 3.811030246936298e-7 kcal/mol/angstrom; maximum scaled error is
7.609023152183738e-8. Both pass the pre-existing 1e-6 / 1e-7 gates.
The separate 20-input cation/water development regressions retain matching
numeric atom types and their prior energy/charge/gradient accuracy.
These are exposed regression inputs, not an unseen test set.

## Pipeline and convergence

Both candidates were evaluated with the same rebuilt common scorer,
using exact coordinate decoding. Every one of the fixed 265 inputs
succeeds and passes the existing independent geometry check; there are
zero gross clashes, violated or unevaluable declared stereoelements.

The optimizer reports convergence for 89/265, compared with 88/265
previously, at the unchanged 200-iteration limit and 1e-4 kcal/mol/angstrom
gradient threshold. The extra report is input 225: coordinates and energy
are identical, and its actual final analytic maximum gradient is
1.9082225835909128e-5. This is a reporting correction, not a newly reached
minimum or a speedup. The convergence audit reconstructs the actual
force-field graph (177 input graphs and 88 explicit-H expansions), checks
exact energy equality on force_field.coords, and verifies every flag
against that graph's analytic gradient. Candidate flag mismatches: zero.

Four force-field coordinate sets change. Input 170's energy decreases by
0.24598673344900135 kcal/mol, and its finite-difference residual force
falls from 17.356900428922017 to 4.281430790342711 kcal/mol/angstrom.
It remains nonconverged. The other three changes are tiny accepted
decreases at geometry constraints. Of 176 nonconverged inputs, 172 reach
the iteration limit and four exhaust admissible steps earlier.

[pre-fix-scorer/](pre-fix-scorer/) preserves the intermediate scorer result
and its diagnosis. Its single apparent violation was a JSON decoding
error, not an accepted violation in the original binary64 geometry.
The boundary fixture still sits at 90 degrees: zero verifier violations
does not establish a chemically realistic planar alkene. Constrained
stationarity and conformer quality remain development work.

## Tests and evidence

Rust: 866 passed, 15 ignored across chematic-ff and chematic-3d.
Python: 92 passed across fixed-coordinate diagnostic, pipeline v2 and
both ensemble suites. Strict Clippy for both crates, scoped formatting
and diff checks pass. Six new Rust test functions cover recovery
overshoot/descent, rejected/non-finite trials, final-step and zero-budget
convergence, and exact JSON coordinate decoding.

[comparison.json](comparison.json) records artifact/source SHA-256 values,
fixed-input measurements, per-input changes and test totals.
[convergence-audit.json](convergence-audit.json) records all 265 gradient
checks. [scorer-build.json](scorer-build.json) pins the rebuilt scorer and
the baseline rescore. The extension and scorer hashes are checked before
and after validation; source files are checked against build-source.json.
The installed rdBase.so is pinned independently, not the whole RDKit
distribution.

From the repository root, in the recorded development environment:

~~~sh
cargo test -p chematic-ff -p chematic-3d --offline -j 2 --quiet
cargo clippy -p chematic-ff -p chematic-3d --all-targets --offline -j 2 -- -D warnings
cargo fmt -p chematic-ff -p chematic-3d --check
cargo build -p chematic-3d --example pipeline_v2_vs_rdkit_common_scorer --offline -j 2
VIRTUAL_ENV=/private/tmp/chematic-mmff-development-20261011 .venv-v1gate/bin/maturin develop --release --offline -j 2 -m crates/chematic-py/Cargo.toml -- -C strip=none -C debuginfo=1
/private/tmp/chematic-mmff-development-20261011/bin/python -B validation/results/mmff94-correctness-2026-10-11/optimizer-recovery/revalidate.py
/private/tmp/chematic-mmff-development-20261011/bin/python -B validation/results/mmff94-correctness-2026-10-11/optimizer-recovery/compare.py
~~~

The exact frozen-energy command is also in candidate-summary.json.
Its underlying exit status 2 correctly records the dirty source tree;
all other provenance/accounting gates pass. The wrapper checks that
specific condition rather than claiming a clean release gate.
No package was released. These results do not establish superior physical
accuracy, conformer search coverage, conformer quality or runtime speed
versus RDKit.
