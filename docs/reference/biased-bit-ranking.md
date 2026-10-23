# Class-biased binary feature ranking

`chematic.rdkit_compat.rdInfoTheory.InfoBitRanker` accepts
`InfoType.BIASENTROPY`, `InfoType.BIASCHISQUARE`, and typed
`SetBiasList(classList: Iterable[SupportsIndex])` inputs. RDKit already provides
these features; chematic previously refused them. This extends the
[explicit/sparse input profile](info-bit-ranker.md) associated with the open
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) tracking issue.
It does not implement the full nanobind migration or resolve that issue.

```python
from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I

ranker = I.InfoBitRanker(3, 2, I.InfoType.BIASCHISQUARE)
ranker.SetBiasList(iter([0]))
for bits, label in [([0], 0), ([0], 0), ([1], 1), ([1], 1)]:
    vector = ExplicitBitVect(3)
    for bit in bits:
        vector.SetBit(bit)
    ranker.AccumulateVotes(vector, label)
rows = ranker.GetTopN(3)
assert rows.tolist() == [[0, 4, 2, 0], [1, 0, 0, 2], [2, 0, 0, 0]]
ranker.SetBiasList([1])  # rerank the same votes for the other class
assert ranker.GetTopN(1).tolist() == [[1, 4, 0, 2]]
```

For each bit, the on-fraction for class `c` is its on-count divided by the
number of observations in that class. Unobserved classes have fraction zero.
A bit keeps its ordinary information-gain or Pearson score when **any**
preferred class has a fraction **at least as large** as the maximum fraction
among nonpreferred classes. Otherwise its score becomes zero. This is the
[RDKit bias rule](https://github.com/rdkit/rdkit/blob/7f3d37d16b71ff8c8a6b0de82c12a8b371e0a8df/Code/ML/InfoTheory/InfoBitRanker.cpp),
not class weighting, a statistical correction, or a p-value. Equality passes;
raw count comparisons would be incorrect when class sizes differ.

An empty preference, including the default, keeps every base score. As in
RDKit's public setter, the list must contain fewer than `nClasses` unique
IDs; a list containing all classes is refused. A single-class ranker therefore
allows only an empty list. Unbiased ENTROPY/CHISQUARE modes validate and store
preferences but ignore them when scoring.

IDs must support the integer index protocol and lie in `[0, nClasses)`;
Python/NumPy booleans, floats, duplicates and out-of-range IDs are refused.
Lists, tuples, NumPy integer sequences and generators are accepted. Validation
is bounded by the configured class count and completes before replacing the
preference; failed iteration or invalid input leaves the old preference and
votes intact. The ranker owns the resulting preference and retains no caller
container or generator. Changing preferences after votes reranks existing
counts without resetting them.

The existing resource limits and float64 output columns
`[bit_id, score, class_0_on_count, ...]` remain unchanged. Rejected features keep
their counts and participate in ranking with zero scores. Exactly equal
computed scores sort by ascending bit ID; this may differ from RDKit's heap
ordering. Masks remain unsupported. No Rust/WASM API, chemical fingerprint
generator, RDKit runtime dependency or performance claim is added.

## Validation

The focused and existing compatibility tests cover preference direction,
inclusive ties, unequal class sizes, multiple preferred classes, empty classes,
empty/default preferences, unbiased modes, state preservation, iterable
ownership, NumPy integers and explicit/sparse vectors. Exhaustive small 2-by-3
contingency tables check the direction rule and Pearson score independently
using exact rational arithmetic.

The seeded live oracle compares 1,152 complete rankers against RDKit 2026.09.1
across 72 datasets, two metrics, both vector representations and four preference
choices. All 34,720 paired votes produce identical per-bit class counts; the
maximum score difference is `1.42e-14`. An independent calculation checks 10,368
scores, with maximum difference `3.56e-15`: exact rational Pearson statistics
and fractions for gating, plus direct mutual information for entropy. Ranking
order differs in 916 cases; per-bit comparisons pass, and no reference ranker
returns duplicate/missing bit IDs in this corpus. This is bounded numerical
evidence, not all-input or exact-order parity.

Reproduction and provenance are in
`validation/results/biased-bit-ranking-2026-10-18/`. The pure Python/NumPy
implementation was tested with a verified existing extension as an import
fixture. No current-main native build or published-wheel validation was run.
