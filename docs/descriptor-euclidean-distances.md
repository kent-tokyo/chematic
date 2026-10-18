# Checked Euclidean distances between descriptor vectors

[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) lists
`GetEuclideanDistMat` among wrappers whose array/sequence inputs still need
typed handling. Chematic adds a distinct Rust/Python API for supplied numeric
descriptor vectors. It does not generate or normalize descriptors.

```python
import chematic

rows = [[0.0, 0.0], [3.0, 4.0], [3.0, 0.0]]
distances = chematic.euclidean_distance_matrix_checked(rows)
assert distances == [5.0, 3.0, 4.0]

molecules = [chematic.from_smiles(s) for s in ["CO", "CCO", "CC(=O)O"]]
rows = [[m.mw, m.logp, m.tpsa] for m in molecules]
distances = chematic.euclidean_distance_matrix_checked(rows)
```

The result is the strict lower triangle in row-major order:
`(1,0), (2,0), (2,1), (3,0), ...`. For `i > j`, its index is `i*(i-1)/2+j`. The diagonal is
omitted and mirrored entries are equal. This matches the ordering documented
for [RDKit's metric calculator](https://www.rdkit.org/docs/source/rdkit.DataManip.Metric.rdMetricMatrixCalc.html),
but the new Python API returns a `list[float]` rather than a NumPy array.

Rust exposes `chematic_chem::descriptor_distances::euclidean_distance_matrix_checked`
for a borrowed slice of vectors implementing `AsRef<[f64]>`, with a separate
`DescriptorDistanceError`. Inputs must be rectangular and finite. Empty
input and a singleton return empty output after validation; zero-feature
vectors have zero distances. The maximum output is 523776 values.

Python accepts nested numeric iterables, including tuples, generators and
NumPy rows with strides. Text/bytes are refused. Outer iteration consumes at
most 1025 rows, and each row consumes at most 257 feature items, to detect the
hard limits below. Source exceptions propagate.
Invalid counts/shapes, nonfinite data, limits and unrepresentable distances
raise `ValueError`; wrong types raise `TypeError`; allocation errors raise
`MemoryError`. NumPy shape metadata is not inspected: empty outer iteration
is empty input regardless of a declared second dimension.

## Computation boundary

- At most 1024 rows and 256 features per row.
- `n*(n-1)/2 * features` must not exceed 50 million. This bounds numeric work;
  it is not a timeout for arbitrary Python iterator callbacks.
- Native output allocation is fallible. Input collection is bounded and uses
  fallible Rust reservations; Python still allocates its own result objects.
- The norm uses scaled sums of squares, avoiding squaring absolute differences.
  For example, one-dimensional distances `1e200` and `1e-200` remain finite and
  positive. Overflowing subtraction or an overflowing final norm refuses the
  whole result, without partial output.

On RDKit 2026.09.1 on this host, valid one-dimensional cases returned `inf`
for `1e200` and zero for `1e-200`. These are independently observed numeric
boundaries, not defects described by the upstream wrapper-typing issue.
Ordinary live-RDKit comparisons and extreme-scale checks against Python's
independent `math.dist` are recorded separately.

Finite-precision error and sensitivity to feature scaling remain. This is
unweighted Euclidean distance in the caller's feature space, without missing-
value imputation, feature weighting, clustering, molecule correspondence or
descriptor agreement claims. Existing molecular graph and 3D distance APIs
are unchanged. A later [typed metric facade](metric-matrix-facade.md) supplies
NumPy-array results through the compatibility namespace. This native API retains
its list return, and no WASM binding is supplied.

The independent probe and dated source/native evidence are under
`validation/rdkit_issues/descriptor_distances/` and
`validation/results/euclidean-descriptors-2026-10-17/`.
