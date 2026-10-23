# Checked sparse binary Tanimoto matrices

[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) includes typed explicit
and sparse bit-vector matrix inputs in its remaining wrapper work. RDKit already
calculates these matrices. This addition supplies a missing native chematic API,
with typed stubs, bounded input consumption and owned source extraction.
It implements a narrow slice of that input contract, not the whole issue.

## Use

```python
import chematic

rows = [[1, 9, 1], [9], [3]]
similarity = chematic.sparse_tanimoto_matrix_checked(rows, metric="similarity")
assert similarity == [0.5, 0.0, 0.0]
distance = chematic.sparse_tanimoto_matrix_checked(iter(rows))
assert distance == [0.5, 1.0, 1.0]
```

The list contains the strict lower triangle `(1,0), (2,0), (2,1), ...`.
Rows may be iterables of integer indices, including NumPy integer arrays,
or objects providing `GetOnBits()`, including RDKit `ExplicitBitVect` and
`SparseBitVect`. A source row stays alive while its returned iterator is consumed.
Input order and duplicate indices do not affect the result; sources are not mutated.
On the inspected RDKit build, `SparseBitVect.GetOnBits()` exposes high IDs above
`2**31-1` as negative signed values; those views are refused with `OverflowError`.
Unsigned caller-supplied index rows support the full u32 range.

Rust exports `chematic_fp::sparse_tanimoto_matrix_checked(&[Vec<u32>],
SparseTanimotoMetric) -> Result<Vec<f64>, SparseMatrixError>`. Python returns
`list[float]`; `metric` is exactly `"similarity"` or `"distance"` (default).
A later [typed metric facade](metric-matrix-facade.md) adds the compatibility
namespace with checked declared widths, NumPy output and RDKit's empty-pair
convention. This native API retains its set/empty policy and list return;
no WASM binding is supplied.

## Contract and limits

Indices span the entire unsigned 32-bit range, including `2**32-1`, without
a dense expansion. They describe binary membership in a caller-defined shared
universe. Source bit-length metadata is ignored, so this does not reproduce
RDKit's handling of mismatched declared dimensions. Count fingerprints,
MinHash signatures and fingerprint generation are outside this API.

Similarity is intersection size divided by union size; distance is `1-similarity`.
Empty input and singleton input produce an empty list after validation. Empty
versus nonempty similarity is zero. Empty versus empty similarity is **one**,
matching chematic's existing `BitVec2048` convention. On this host RDKit
2026.09.1 returns **zero** for empty/empty; that deliberate policy difference is
recorded separately from ordinary comparisons.

| Bound | Maximum |
| --- | ---: |
| Fingerprints | 1,024 |
| Raw IDs per fingerprint | 65,536 |
| Total raw IDs | 1,000,000 |
| Pair-bit work `(n-1)*sum(unique counts)` | 50,000,000 |
| Output values | 523,776 |

Raw caps apply before sorting or duplicate removal. The pair-work bound is
conservative and checked before calculating the matrix. Native numeric output
occupies at most about 4 MiB, excluding normalized input and Python objects.
Fallible reservations return errors without exposing partial output. Python
consumes at most one excess item to detect a cap. These bounds do not time out
arbitrary caller iterator/getter callbacks.

Wrong container/index types, including Python `bool` indices, raise `TypeError`;
negative or out-of-range indices raise `OverflowError`; invalid metrics/limits
raise `ValueError`; allocation failure raises `MemoryError`. Iteration and
`GetOnBits()` call errors propagate. Attribute lookup treats `AttributeError`
as absence of `GetOnBits()` and falls back to iterating the row.

## Local evidence

The original deterministic probe at
`validation/rdkit_issues/sparse-tanimoto/probe.py` compares five input forms
against an independent Python set calculation and live RDKit 2026.09.1 explicit
and sparse matrices. At most one row is empty in each ordinary dataset;
empty/empty policy, singleton refusal upstream and high sparse IDs are separate
diagnostics.
Results, artifact hash and upstream issue snapshot are retained in
`validation/results/sparse-tanimoto-2026-10-17/`. These are local candidate
checks, not a published package or whole-wrapper compatibility claim.
