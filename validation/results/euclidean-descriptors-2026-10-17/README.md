# Checked descriptor distances, 2026-10-17

Base: `de41dadd8f123a4b1368937c7dd86df823434588`.
`source-sha256.json` records source, tests, manifests, probe and binding inventory.
`summary.json` records the host, compiled native hash and measured counts.

## Scope and upstream

[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) remained OPEN at
inspection; `upstream-issue.json` retains the API snapshot. Its remaining
wrapper work explicitly includes array/sequence inputs to `GetEuclideanDistMat`.
The new native Rust/Python API is an explicit descriptor-space capability
gap, with typed stubs and bounded iterable handling. It does not implement the
whole tracking issue, nanobind or a compatibility-shim alias.

Independent valid-input probes on RDKit 2026.09.1 return infinity for the
one-dimensional distance `1e200` and zero for `1e-200`. These observations
are retained in `upstream-numeric-baseline.json`. They are numeric behavior
discovered locally, not the bug described in the upstream wrapper-typing issue.
No upstream test data or implementation was copied.

## Results

- New Rust tests: **7 passed**, covering triangle layout, duplicates, zero
  dimensions, validation, exact caps, extreme norms and source correspondence.
- Existing native graph-distance tests: **9 passed**.
- New Python tests: **33 passed**, including numeric NumPy float64/float32/int64
  rows, strides, generators, bounded consumption and original exceptions.
- Existing screening-enrichment registration tests: **25 passed**.
- All **74 tests** passed with no failures or skips/ignored tests.
- Live ordinary RDKit comparisons: **300 cases**, `n = 1/2/3/8/16`,
  `d = 0/1/3/16/64`, three scales and four fixed seeds. The native runs each
  case as lists, tuples, inner/outer generators, strided NumPy rows and reversed
  row order: **1500 native runs**. Python `math.dist` is an independent oracle.
  Both engines pass relative tolerance `2e-12` or absolute tolerance `1e-12`.
  Maximum ordinary absolute error to `math.dist` for both is
  `2.7284841053187847e-12`, at distances where the relative tolerance applies.
- One descriptor workflow: MW/logP/TPSA of six native molecules are fed to
  both distance engines. This is shared descriptor input, not a comparison of
  how chematic/RDKit generate descriptors or a calibrated similarity model.
- Extreme controls: **5 passed** against `math.dist`, using relative tolerance
  `2e-15` and **zero absolute tolerance**. Native output retains `1e200`,
  `1e-200`, the minimum positive subnormal, the finite two-axis `1e308` norm
  and the maximum finite float. The RDKit diagnostics return infinity or zero;
  JSON uses the string `"inf"` for nonfinite reference values.
- Strict all-target Clippy for both changed crates, formatting, source
  consistency gates and strict MkDocs pass. Compressed logs retain outputs.
- Source-only binding inventory updates Python exports/registrations from
  213 to 214 and keeps WASM at 317; this bookkeeping is not semantic parity.
- Absence baseline: retained PR #974 native from older `bef0d008` plus unrelated
  MMFF custom-charge work lacks the new function (`AttributeError`). This is
  labeled as an older-base retained native, not a fresh build of this base.
  Its SHA256 is `6f194706eb5ed4846d097b72f79ea29484e7f0bf90145bc153bdf8729212e2ad`.
  The current base source also lacks the new API.

Outputs are condensed Python lists in lower-triangle row order, not NumPy
arrays. Empty outer iteration returns empty output; NumPy shape metadata is
not interpreted. Finite-precision error remains. No feature normalization,
clustering, descriptor generation, chemical-equivalence or whole-RDKit parity
claim is made. Existing graph and Cartesian distance APIs remain separate.

## Reproduction

Run from this checkout with the matching macOS arm64 Python 3.13 runtime:

```sh
export CARGO_TARGET_DIR=/private/tmp/chematic-euclidean-target
export CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=1
export CARGO_PROFILE_DEV_DEBUG=0 CARGO_PROFILE_TEST_DEBUG=0
export CARGO_PROFILE_DEV_STRIP=none
export CARGO_PROFILE_DEV_OPT_LEVEL=1 CARGO_PROFILE_TEST_OPT_LEVEL=1
export CARGO_PROFILE_DEV_CODEGEN_UNITS=1 CARGO_PROFILE_TEST_CODEGEN_UNITS=1
export PYO3_PYTHON=/Library/Frameworks/Python.framework/Versions/3.13/bin/python3
cargo test -p chematic-chem --offline --test descriptor_distances
cargo test -p chematic-chem --offline --test native_distance_matrix
cargo clippy -p chematic-chem -p chematic-py --all-targets --offline --no-deps -- -D warnings
cargo rustc -p chematic-py --offline -- -C strip=none -C link-arg=-undefined -C link-arg=dynamic_lookup
mkdir -p /private/tmp/chematic-euclidean-python/chematic
cp "$CARGO_TARGET_DIR/debug/libchematic.dylib" /private/tmp/chematic-euclidean-python/chematic/chematic.so
cp crates/chematic-py/python/chematic/{__init__.py,__init__.pyi,rdkit_compat.py,py.typed} /private/tmp/chematic-euclidean-python/chematic/
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-euclidean-python python3.13 -m pytest -q crates/chematic-py/tests/test_descriptor_distances.py crates/chematic-py/tests/test_screening_enrichment.py
# The second path contains independently installed RDKit 2026.09.1 / NumPy.
PYTHONPATH=/private/tmp/chematic-euclidean-python:/private/tmp/chematic-path-rdkit python3.13 validation/rdkit_issues/descriptor_distances/probe.py --output /private/tmp/descriptor-distance-comparison.json
```

`comparison.json.gz` retains all ordinary rows, expected/native/reference
values, errors and extreme diagnostics. The compiled extension is used by
the Python tests and probe. No upstream changes or comments were posted.
