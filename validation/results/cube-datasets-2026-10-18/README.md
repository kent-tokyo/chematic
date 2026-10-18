# Explicit CUBE dataset selection — 2026-10-18

Related open issue: https://github.com/openbabel/openbabel/issues/2579
Related upstream PR: https://github.com/openbabel/openbabel/pull/2580 (open).

Native chematic scalar output already uses nonnegative NATOMS and no orbital
ID record; this change fills the missing multi-dataset input/selection boundary.
It does not change existing Rust/Python/WASM scalar readers or writers.

- Base: `f99c7b8ab33098a2e1fd9fd568aee370bf79ab38`.
- `pytest.log`: 106 tests pass, 56 new and 50 existing format regressions.
  Covers both dataset forms, asymmetric nonorthogonal grid order, ECP nuclear
  charge, immutable ownership, malformed/trailing data, byte/size limits and
  actual 64-dataset/262144-value/50000-atom boundaries.
- `oracle.py`, `oracle-summary.json`, `oracle.log`: ASE 3.26.0 checks 192 inputs,
  240 original positive-NVAL fields and 480 scalar exports (12724 scalar values).
  Values and geometry agree exactly; source IDs/order and authored exact voxel
  values also check the negative form. ASE accepts 24 negative-form originals
  and refuses 72 multi-field originals. These failures are reported separately,
  not counted as original-input cross-read successes.
- `orbital-input.cube`, `single-input.cube`: original synthetic controls, with
  ECP nuclear charge 4 for carbon and IDs 7/9 in the orbital form. No upstream
  structure/attachment or third-party source implementation is copied.
- `openbabel-*.cube`, `openbabel-*.stderr`, `openbabel-summary.json`: local
  Open Babel 3.2.1 still writes negative NATOMS for single-field input; the
  two-field rewrite renumbers IDs 7/9 to 1/2, and both rewrite nuclear charge
  4 to 6. Scalar values survive this control. These observations are distinct
  from field/geometry cross-read successes and are not all-version/API claims.
- `typing-valid.log`, `typing-invalid.log`: strict helper/public consumer pass;
  four expected errors for the negative consumer. Imported modules/stubs are
  silent, so package-wide strict typing is not claimed.
- `docs-and-typing.log`, `mkdocs.log`: Python 3.9 grammar, executed guide, strict
  docs build and three rendered API anchors pass.
- `source-gates.log`: source stubs, manifest, binding inventory, unsafe surface,
  refactor budget and canonical names pass. No wheel supplied.
- `native-source-comparison.log`: Cube reader, shared volumetric type and Python
  grid binding match the verified native fixture's source revision. No fresh
  current-main native build or whole-extension equivalence claim.
- `provenance.json`: source/reference/evidence hashes and checked issue status.

The independent implementation uses h5cube's public format grammar:
https://h5cube-spec.readthedocs.io/en/latest/cubeformat.html

Hard caps: 8 MiB UTF-8 input, 50000 atoms, 64 datasets, 262144 grid points and
total scalar values, 16384 characters per input line and 16 MiB scalar export.
The API materializes bounded data and requires explicit positional selection.
It keeps the complete source and omitted positions alongside the output.
Lengths remain in bohr; negative axis counts, ghost atoms, nonfinite scalars,
duplicate IDs, extra/missing values and Fortran D notation refuse. Values are
float64; no physical orbital interpretation, scalar-unit conversion, resampling,
GUI execution, performance or published-wheel validation is claimed.

## Reproduce

Put this revision's Python package with a compatible native extension on
PYTHONPATH; add an independent ASE 3.26.0 installation for the oracle.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_cube_datasets.py \
  crates/chematic-py/tests/test_new_formats.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 python3.13 validation/results/cube-datasets-2026-10-18/oracle.py
obabel -icube validation/results/cube-datasets-2026-10-18/orbital-input.cube -ocube
obabel -icube validation/results/cube-datasets-2026-10-18/single-input.cube -ocube
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-cube-site
```

Strict mypy uses strict=True, follow_imports=silent,
follow_imports_for_stubs=True, ignore_missing_imports=True and MYPYPATH pointing
to the candidate package. Check the helper with typing_valid.py, then separately
check typing_invalid.py and expect four errors. Native fixture identity and
relevant-source comparisons are recorded separately from the Python-source and
reference-corpus hashes.
