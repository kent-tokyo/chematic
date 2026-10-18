# Typed metric matrix facade — 2026-10-18

Related open issue: https://github.com/rdkit/rdkit/issues/9601

Existing native descriptor distances and sparse-set Tanimoto already work.
This change adds the missing typed `rdkit_compat.rdMetricMatrixCalc` namespace
and owned 1D float64 NumPy return contract. Binary sources provide declared
widths and on-bit positions; staged validation checks equal widths and position
ranges before native calculation. The facade adopts RDKit's empty/empty
similarity 0 and distance 1; the native sparse-set API remains unchanged.

- Base: `f99c7b8ab33098a2e1fd9fd568aee370bf79ab38`.
- `pytest.log`: 281 tests pass, including 46 new, 40 native sparse-set,
  33 descriptor-distance and 162 compatibility cases. Exact small binary sets,
  asymmetric pair order, source lifetime/error propagation, owned results,
  numeric extremes and actual row/raw-ID/aggregate boundaries are covered.
- `oracle.py`, `oracle-summary.json`, `oracle.log`: live RDKit 2026.09.1 and
  independent math.dist/set formulas. 512 Euclidean plus 1024 binary calls,
  262128 pair values; binary matches exactly. Maximum Euclidean RDKit error
  4.263256414560601e-14, maximum math.dist error 2.842170943040401e-14.
  Sources include nested generators, strided descriptor arrays, explicit/sparse
  vectors, mixed-source generators and chematic's own explicit vectors.
- `reference-diagnostics.json`: eight separately recorded numeric/input-profile
  controls: huge/tiny distance, empty and mixed-type collections, high signed
  sparse views and unequal declared widths. They are excluded from ordinary
  parity successes; no ordinary failures were discarded.
- `typing-valid.log`, `typing-invalid.log`, `mypy.ini`: strict helper/public
  consumer pass and four expected negative errors, including wrong result type.
  Return annotations describe ndarray[tuple[int], dtype[float64]] directly.
  Imported legacy modules/stubs are silent; NumPy imports are followed normally.
  No package-wide strict typing claim.
- `docs-and-typing.log`, `mkdocs.log`: Python 3.9 grammar, executed guide, strict
  docs build and five rendered API anchors pass.
- `source-gates.log`: source stubs, manifest, binding inventory, unsafe surface,
  refactor budget and canonical names pass. No wheel supplied.
- `native-source-comparison.log`: both native metric implementations and Python
  bindings match the hash-verified extension fixture revision. No fresh native
  current-main build or whole-extension equivalence claim.
- `provenance.json`: source/runtime/reference/evidence hashes and issue snapshot.

Binary caps: 1024 sources, 65536 raw IDs per source, one million raw IDs and
50 million native pair-bit work units. Positions stay sparse with no dense
expansion; duplicate IDs are idempotent but count against raw limits. Widths
must match in 0..4294967295, and IDs must be unsigned-32-bit representable and
less than width. Bool widths/IDs, negative signed views and count/raw-row input
refuse. Empty/singleton collections return empty arrays. These profiles can
differ from RDKit and do not assert all-input compatibility.

Descriptor caps/numerics are inherited: 1024 rows, 256 features, 50 million
pair-feature operations, rectangular finite data and representable distances.
No descriptor generation, scaling, physical interpretation, new native
algorithm, performance or published-wheel validation is claimed.

## Reproduce

Put this revision's Python package with a compatible native extension on
PYTHONPATH; add independent RDKit 2026.09.1 for the oracle. Tests do not need
RDKit. The runtime fixture and source-match scope are recorded in provenance.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_metric_matrix_facade.py \
  crates/chematic-py/tests/test_sparse_tanimoto.py \
  crates/chematic-py/tests/test_descriptor_distances.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 python3.13 validation/results/metric-matrix-facade-2026-10-18/oracle.py
python3.13 -m mypy --config-file validation/results/metric-matrix-facade-2026-10-18/mypy.ini \
  crates/chematic-py/python/chematic/_metric_matrix.py \
  validation/results/metric-matrix-facade-2026-10-18/typing_valid.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-metric-site
```

For typing, set MYPYPATH to the candidate package parent and check the negative
consumer separately, expecting exit 1 with four errors. The checked-in config
keeps legacy imports silent but follows NumPy's typing normally.
