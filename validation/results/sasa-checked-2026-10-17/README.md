# Explicit SASA radius contract, 2026-10-17

Base: `bef0d008d49f3d432fc35a900ccbf6acab2cd957`. Source hashes are in
`source-sha256.json`; host, native hash and measured counts are in
`summary.json`. The build uses Python 3.13.6 on macOS arm64.

## Upstream and scope

[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) and
[PR #9611](https://github.com/rdkit/rdkit/pull/9611) were OPEN when inspected.
The issue includes checking one FreeSASA radius per atom. Snapshots are
`upstream-issue.json` and `upstream-pr.json`.

Chematic previously exposed SASA with built-in radii. This change adds an
opt-in Rust/Python function for supplied coordinates/radii with strict counts,
finite values, numerical range and work bounds. It uses independent midpoint
Fibonacci Shrake-Rupley sampling. It is not a FreeSASA port or general RDKit
API parity, and does not implement other work in the tracking issue.
Existing SASA APIs are unchanged.

## Results

- New Rust geometry/error tests: **8 passed**. Three collinear spheres test
  multiple occluders against independent nonoverlapping spherical-cap areas.
- Existing Rust SASA unit tests: **27 passed**, 674 unrelated tests filtered.
- New Python tests: **43 passed**, including explicit H cardinality, iterables,
  translated/reordered geometry, numerical and resource limits.
- Strict Clippy for both changed crates, all targets: passed.
- Analytic probe: **96 geometries, 288 native runs**, zero analytic failures.
  Largest per-atom error / full isolated sphere area: **0.0015869140625**;
  tolerance **0.005**. Translation allows one sample's area because translating
  decimal coordinates changes their floating-point representation. Atom-order
  reversal must exactly reverse the per-atom areas.
- Extreme-scale controls: **5 passed**, `1e-100` through `1e100`.
- RDKit 2026.09.1 default Lee-Richards values are diagnostic only: **95 within
  the separate 0.025 full-sphere tolerance, 1 nonfinite**. At radii `[1,1.5]`,
  centers `[0,0,0]` and `[0,0,1]`, probe zero, its total is NaN and its atom
  values are `[5.654866776461628, NaN]`. JSON uses `null` for nonfinite values
  and records the status separately. No native failure is hidden by this.
- Isolated upstream short-radii probe: two atoms / one radius returned
  `72.38229473870885` without an exception or crash on this host. This is an
  observation of missing validation, not a stable expected numerical value.
  Core dumps were disabled and the subprocess had a 15-second timeout.
- Absence baseline: a retained native from PR #972, built on the same base
  plus an unrelated CML writer, fails the single-sphere new-API test with
  `AttributeError`. This is not represented as a fresh unmodified-main build.
  The base source also lacks the API. Its native SHA256 is
  `6fada006000ecb1ad861a3846f61b06b6efc04d41ea4a4c7b587ef6342fd3443`.

Two-sphere analytic cap checks, plus the three-sphere unit control, do not
establish general molecular, protein or high-density accuracy. Finite sampling
has orientation dependence. The work cap is conservative and is not a timeout.

## Reproduction

Run from this checkout. Set a writable target and these build parameters:

```sh
export CARGO_TARGET_DIR=/private/tmp/chematic-sasa-checked-target
export CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1
export CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0
export CARGO_PROFILE_DEV_STRIP=none
export CARGO_PROFILE_DEV_OPT_LEVEL=1 CARGO_PROFILE_TEST_OPT_LEVEL=1
export CARGO_PROFILE_DEV_CODEGEN_UNITS=1 CARGO_PROFILE_TEST_CODEGEN_UNITS=1
export PYO3_PYTHON=/Library/Frameworks/Python.framework/Versions/3.13/bin/python3
cargo test -p chematic-3d --offline --test sasa_checked
cargo test -p chematic-3d --offline --lib sasa::tests
cargo clippy -p chematic-3d -p chematic-py --all-targets --offline --no-deps -- -D warnings
cargo rustc -p chematic-py --offline -- -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup

mkdir -p /private/tmp/chematic-sasa-checked-python/chematic
cp "$CARGO_TARGET_DIR/debug/libchematic.dylib" /private/tmp/chematic-sasa-checked-python/chematic/chematic.so
cp crates/chematic-py/python/chematic/{__init__.py,__init__.pyi,rdkit_compat.py,py.typed} /private/tmp/chematic-sasa-checked-python/chematic/
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-sasa-checked-python python3.13 -m pytest -q crates/chematic-py/tests/test_sasa_checked.py

# The second path contains independently installed RDKit 2026.09.1.
PYTHONPATH=/private/tmp/chematic-sasa-checked-python:/private/tmp/chematic-path-rdkit python3.13 validation/rdkit_issues/sasa_checked/probe.py --output /private/tmp/sasa-comparison.json
```

`comparison.json.gz` retains every geometry, analytic/native/reference value,
diagnostic status, numerical threshold and the isolated short-radii result.
Compressed logs retain source tests, baseline failure, native build, Clippy,
documentation build and the probe summary. `gates.log.gz` records source
consistency checks. Python uses the built native package, not a pure Python
replacement. No upstream comments or changes were posted.
