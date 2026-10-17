# Typed binary feature ranking

`chematic.rdkit_compat.rdInfoTheory.InfoBitRanker` ranks binary fingerprint
features by information gain or Pearson chi-square against caller-supplied
integer class labels. It implements part of the typed explicit/sparse bit-vector
input work in [RDKit #9601](https://github.com/rdkit/rdkit/issues/9601),
verified open on 2026-10-17. This is an independently implemented Python/NumPy
compatibility API; RDKit is not required at runtime.

```python
from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory

ranker = rdInfoTheory.InfoBitRanker(4, 2)
for on_bits, label in [([0, 1], 0), ([1, 2], 1)]:
    fp = ExplicitBitVect(4)
    for bit in on_bits:
        fp.SetBit(bit)
    ranker.AccumulateVotes(fp, label)

rows = ranker.GetTopN(4)
# Columns: bit ID, information gain (bits), class-0 on-count, class-1 on-count.
```

Input vectors implement `GetNumBits()` and `GetOnBits()`. This structural
protocol accepts chematic's `ExplicitBitVect`, actual RDKit explicit/sparse
binary vectors, or other objects with those methods. No new sparse-vector
storage class is added. Width and bit IDs implement the integer index protocol,
as do dimensions, class labels and top counts; Python and NumPy integers are
supported. Booleans and floats are rejected. Raw strings, packed bytes, count
fingerprints and plain lists of bit IDs are not input vector alternatives.

The width must match the ranker. IDs must be unique and in `[0, nBits)`;
class labels must be in `[0, nClasses)`. The complete current vector is
validated before any counters change. Invalid labels are checked before
reading the vector, and failures during `GetOnBits()` iteration propagate
without committing a partial vote. Vectors are not retained; later source
mutation has no effect. Empty vectors count as observations with all bits off.
This guarantee covers the ranker's own updates, not side effects performed
by user-defined vector methods.

Supported dimensions are 1–65,536 bits, 1–256 classes, and at most 1,048,576
bit/class count entries, checked before allocation. The count table is uint64;
at most 1,048,576 observations can be accumulated. These bounds also keep bit
IDs and returned counts exactly representable in float64. Work and temporary
storage for ranking scale with the bounded count table, plus sorting by bit.

`GetTopN(num)` returns a new float64 array with `num` rows and `nClasses + 2`
columns. The default ENTROPY information gain is label entropy minus weighted conditional label
entropy for bit-on/off subsets, using log2. Zero-probability terms contribute
zero. Floating-point results are clipped to `[0, label entropy]`. Ranking is
by decreasing gain and then ascending bit ID for exactly equal computed
scores. RDKit's heap gives a different, unstable order among equal scores;
exact ranking-order parity is not promised. Empty rankers raise `ValueError`.
After at least one observation, `GetTopN(0)` returns an empty array with the
correct column count. Changing a returned array does not alter later results.

`InfoType.ENTROPY` and [`InfoType.CHISQUARE`](chi-square-bit-ranking.md) are
implemented. CHISQUARE returns an uncorrected Pearson statistic, not information
gain or a p-value; it uses the same count table and deterministic tie policy.
Biased enum metrics raise `ValueError` at construction. `SetBiasList` and
`SetMaskBits` explicitly raise `NotImplementedError`. Biased scores, masks, correlation
matrices and file output are outside this profile. This is supervised binary
feature ranking; it does not generate chemical fingerprints or infer activity
labels. No speed advantage over RDKit is claimed.

The output contract follows the official
[RDKit InfoTheory documentation](https://www.rdkit.org/docs/source/rdkit.ML.InfoTheory.rdInfoTheory.html).
Local evidence in `validation/results/info-bit-ranker-2026-10-17/` includes
source and runtime-fixture hashes, input-state tests, public type-check examples
and a live RDKit 2026.09.1 comparison. The implementation uses Python/NumPy
only; native code was reused as an import fixture and was not rebuilt here.

The original entropy validation passed 207 Python tests, including 45 new ranker tests, and
768 live explicit/sparse ranker comparisons (10,560 paired votes). Per-class
counts matched exactly; maximum information-gain difference was 1.23e-15.
All 1,536 order-reversal/data-duplication checks passed. Ranking order differed
in 462 comparisons, with the sorted score sequences agreeing within 1e-12;
this is consistent with the documented tie/roundoff policy, not an exact-rank
compatibility claim. The isolated implementation passes strict mypy, public
valid examples type-check, and four deliberately wrong inputs are rejected.
This is not a package-wide strict typing or published-wheel validation claim.
