# Typed metric matrix namespace

[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) includes typed array,
sequence and explicit/sparse bit-vector inputs for metric matrices. Native
chematic already computes checked [descriptor distances](descriptor-euclidean-distances.md)
and [sparse Tanimoto matrices](sparse-tanimoto-matrix.md), but has no matching
`rdMetricMatrixCalc` namespace or NumPy-array results.

Python now exposes three typed adapters through `chematic.rdkit_compat`.
They reuse those native algorithms and bounds; no new native metric or RDKit
runtime dependency is added. This covers one wrapper boundary from the tracking
issue, without claiming a complete metric library or whole-issue solution.

```python
from chematic.rdkit_compat import ExplicitBitVect, rdMetricMatrixCalc as Metric

rows = [[0.0, 0.0], [3.0, 4.0], [3.0, 0.0]]
assert Metric.GetEuclideanDistMat(descripMat=rows).tolist() == [5.0, 3.0, 4.0]

a, b, empty = (ExplicitBitVect(8) for _ in range(3))
a.SetBit(1)
a.SetBit(3)
b.SetBit(3)
assert Metric.GetTanimotoSimMat(bitVectList=[a, b, empty]).tolist() == [0.5, 0.0, 0.0]
assert Metric.GetTanimotoDistMat([empty, empty]).tolist() == [1.0]
```

Every result is a fresh owned one-dimensional `numpy.float64` array. Entries
are the strict lower triangle `(1,0), (2,0), (2,1), (3,0), ...`, without a
diagonal. This follows the
[RDKit metric calculator's documented layout](https://www.rdkit.org/docs/source/rdkit.DataManip.Metric.rdMetricMatrixCalc.html).
Inputs are neither mutated nor retained by results. Returned arrays can be
changed independently. The namespace is imported from `chematic.rdkit_compat`;
no new `rdkit.DataManip` package is installed.

## Descriptor rows

`GetEuclideanDistMat(descripMat)` accepts rectangular numeric nested iterables,
including lists/tuples, generators and strided NumPy rows. Existing native
validation applies: at most 1,024 rows, 256 features, 50 million pair-feature
operations and finite float64 inputs/results. Text containers/rows refuse.
Empty or singleton inputs give an empty array; zero-feature vectors give zero
distances. Iteration consumes at most 1,025 rows or 257 items per row to detect
overflow. Array shape metadata is not inspected: no outer rows means empty
input regardless of a declared second dimension.

The existing scaled norm retains representable distances such as `1e200` or
`1e-200`, while the inspected RDKit 2026.09.1 returns infinity or zero on those
controls. These numeric differences already existed in the native API; they
are separate from the upstream wrapper-typing issue. Descriptors, scaling,
weights and physical meaning remain caller supplied.

## Binary vector sources

`GetTanimotoSimMat(bitVectList)` and `GetTanimotoDistMat(bitVectList)` accept
iterables of sources satisfying the exported `MetricBitVector` protocol:
`GetNumBits()` returns an integer width and `GetOnBits()` returns integer
positions. chematic explicit vectors and RDKit explicit/sparse vectors work;
mixed source types and lazy sources are accepted. Raw position rows and count
vectors are outside this profile. No dense expansion occurs.

All widths must be equal and in `0..4294967295`; every ID must be nonnegative,
unsigned-32-bit representable and less than the declared width. Python/NumPy
bool widths/IDs refuse. IDs can be unordered or duplicated; duplicates are
idempotent but count toward raw limits. High unsigned IDs remain sparse.
The inspected RDKit sparse `GetOnBits()` can expose high IDs as negative signed
values; such views refuse rather than silently reinterpreting them. A caller
may instead supply a protocol source with checked unsigned positions.

Hard caps are 1,024 sources, 65,536 raw IDs per source, one million raw IDs
overall and native pair-bit work `(n-1)*sum(unique ID counts) <= 50 million`.
The largest output has 523,776 values. All sources and positions are copied
before native calculation; a lazy source stays alive while its iterator is
consumed. At most one excess row/ID is consumed to detect a cap. These bounds
do not impose a timeout on arbitrary caller methods.

Similarity is intersection size divided by union size, with empty/empty
similarity **0** and distance **1**, matching RDKit. This facade-specific rule
differs from the existing native sparse-set API's similarity 1/distance 0;
that native API is unchanged. Empty/singleton collections return an empty
array; the inspected RDKit raises on them. Unequal declared widths refuse,
even where a particular RDKit input might be accepted. These are explicit
input profiles, without all-input or fingerprint-generation parity claims.

Wrong types raise `TypeError`; invalid widths/positions, limits, shapes and
numeric ranges raise `ValueError`; IDs outside u32 raise `OverflowError`.
Allocation failures propagate as `MemoryError`. Source exceptions propagate;
`AttributeError` during method lookup is treated as an absent method. Results
are never exposed partially if source collection or native calculation fails.

## Validation

281 Python tests pass: 46 new cases and existing descriptor, sparse-set and
compatibility regressions. They cover exact set expectations over all small
bit patterns, lower-triangle ordering, NumPy/generator inputs, empty conventions,
high IDs, ownership, source lifetimes/errors and actual raw-ID/row/total limits.

RDKit 2026.09.1 comparisons cover 512 Euclidean and 1,024 binary calls, totaling
262,128 pair values. Tanimoto results agree exactly; the maximum Euclidean
difference is `4.27e-14`, and maximum difference from independent `math.dist`
is `2.85e-14`. Binary set formulas provide a separate independent reference.
Eight numeric/empty/source-profile diagnostic cases are recorded separately from
ordinary comparisons; no ordinary failures are discarded.

Evidence and the executable oracle are in
`validation/results/metric-matrix-facade-2026-10-18/`. Tests reuse a hash-verified
native extension with unchanged relevant metric implementations and bindings.
No fresh current-main native build, published wheel or performance claim is made.

## Python reference

::: chematic._metric_matrix
