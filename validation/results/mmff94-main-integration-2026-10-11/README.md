# MMFF94 main integration, 2026-10-11

Source candidate based on main 8f5a49a32dfbbce39109e5feaaa04ea4d7e1a3f5 (version 1.1.0).
Original commit: 0bae2e0bdec5ef299d1c4dd6e16e1c06cb164c61, developed on a version 1.0.20 checkout.
Historical measurements in ../mmff94-correctness-2026-10-11/ remain unchanged;
they are not measurements of this main candidate.

## Integration scope

Retain main's newer numeric typing, angle, stretch-bend, van der Waals,
out-of-plane, molecule-aware torsion, MMFF94S and BFGS implementations.
Transplant bounded L-BFGS recovery with Armijo backtracking, finite-state checks
and final-state convergence reporting, preserving main's termination enum and
bindings. Keep the bond-aware empirical torsion helper as an additive API.
Preserve coordinate binary64 bits with serde_json float_roundtrip.

New tests exposed triangle torsions in coverage accounting and invalid outer
numeric types resolving through an out-of-plane wildcard. Exclude repeated
torsion endpoints and reject unregistered nonzero outer types, retaining explicit
zero table queries and main's larger valid eligibility list.
The shared-coordinate diagnostic checks indexed bond graphs and element order,
and retains main's per-term energy output.

Unrelated local IO and WASM commits are excluded. Main's versions and dependencies
are retained.

## Validation

- Rust: 959 passed, 4 ignored. Command: CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0 CARGO_INCREMENTAL=0 cargo test -p chematic-ff -p chematic-3d --offline -j 2 --lib --tests --no-fail-fast --quiet.
- Clippy exit 0: CARGO_PROFILE_DEV_DEBUG=0 CARGO_INCREMENTAL=0 cargo clippy -p chematic-ff -p chematic-3d --all-targets --offline -j 2 -- -D warnings.
- Format check exit 0: cargo fmt -p chematic-ff -p chematic-3d --check.
- Python diagnostic, pipeline and ensemble suites: 95 passed.
- Frozen 265-input pipeline: 265 successful, 265 independently sound, 0 stereo violations, 0 unevaluable checks; reported convergence 135/265.
- Shared RDKit-produced explicit-H coordinates: 262 comparable, two RDKit embedding failures, one declared unsupported. Maximum total energy difference 9.6633812063373625e-13 kcal/mol. Seven term results are retained in energy-summary.json. Full-corpus gradients were not rerun.

RDKit runtime 2026.03.6 and its pinned binary hash are recorded in validation
records. Python uses an unoptimized dev source wheel extracted into this worktree,
separate from the original editable installation. This is not a timing comparison
or a published artifact. The source patch and file hashes identify the tested code
before its final commit exists.

validate_python.py and validate_energy.py retain exact commands and artifact
hashes. The energy diagnostic's exit 2 is expected solely because the source-tree
clean check runs before committing; its other gates pass.

These bounded regression checks do not establish broad RDKit superiority,
performance, physical conformer accuracy, or a release.
