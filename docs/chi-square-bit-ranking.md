# Chi-square binary feature ranking

`chematic.rdkit_compat.rdInfoTheory.InfoBitRanker` now accepts
`InfoType.CHISQUARE`. This extends the typed explicit/sparse binary-vector
input work associated with [RDKit #9601](https://github.com/rdkit/rdkit/issues/9601),
which remains open. RDKit already provides this metric; chematic previously
refused it. This is one part of the tracking issue, not a complete wrapper port.

```python
from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I

ranker = I.InfoBitRanker(3, 2, I.InfoType.CHISQUARE)
for bits, label in [([0], 0), ([0], 0), ([1], 1), ([1], 1)]:
    vector = ExplicitBitVect(3)
    for bit in bits:
        vector.SetBit(bit)
    ranker.AccumulateVotes(vector, label)
rows = ranker.GetTopN(3)
assert rows[:, 0].tolist() == [0.0, 1.0, 2.0]
assert rows[:, 1].tolist() == [4.0, 4.0, 0.0]
```

Each feature defines a 2-by-class contingency table: bit on/off observations
against caller-supplied class labels. Its score is the Pearson statistic
`sum((observed - expected)**2 / expected)`, where
`expected = row_total * class_total / observation_count`. Zero expected cells
contribute zero: unobserved classes and constant features are supported, and
a single observed class yields zero. No Yates correction is applied. The
result is a statistic, not a p-value or evidence of predictive performance.

The float64 output columns remain `[bit_id, score, class_0_on_count, ...]`.
Scores sort descending, with ascending bit ID for exactly equal computed scores.
Every selected bit appears once. This deterministic policy can differ from
RDKit's heap and floating-point ordering. Duplicating all observations doubles
the chi-square statistic and counts; entropy information gain instead stays
unchanged. Reversing vote order leaves results unchanged.

The existing [input and resource contract](info-bit-ranker.md) applies:
structural `GetNumBits()` / `GetOnBits()` vectors, bounded dimensions and record
count, uint64 counters, staged validated votes, independent output arrays,
and explicit errors. `GetTopN(0)` works after a vote; empty rankers still refuse.
Biasing, masks, p-values, chemical fingerprint generation and Rust/WASM APIs
remain outside this Python/NumPy addition. It has no RDKit runtime dependency.

## Validation

The new tests cover hand-calculated tables, empty rows/classes, complement and
duplication laws, failed input iteration, source/output ownership, record limits,
and counts beyond 65535. Existing entropy and compatibility tests also pass.

The seeded oracle attempts 768 live ranker comparisons against RDKit 2026.09.1
using explicit and sparse vectors, with 32,896 paired votes. Of those, 764 return
one row per bit upstream and match per-class counts exactly and scores within
1.96e-14 (chi-square) / 1.34e-15 (entropy). Four CHISQUARE cases return duplicate
bit IDs in RDKit and are recorded separately, with complete inputs and outputs;
they are not counted as parity successes. A reduced table is a regression test.
This is a locally observed numerical boundary, separate from #9601's typing work.

An independent exact-rational calculation checks all 2,064 distinct chi-square
scores, including those four reference anomalies, with maximum difference
1.78e-15. All 384 vote-reversal and 384 data-duplication checks pass, as do 18
separate degenerate-table checks. These establish the bounded statistic contract,
not exact ranking order or all-input RDKit parity.

Reproduction, input/source/native-fixture hashes and logs:
`validation/results/chi-square-bit-ranker-2026-10-18/`. No native code changed;
the verified existing local extension is only an import fixture. Published
wheel behavior and performance comparisons were not tested.
