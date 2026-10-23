# Typed feature masks for information ranking

`InfoBitRanker.SetMaskBits(maskBits)` now selects the eligible binary features
for the ENTROPY and CHISQUARE modes in `chematic.rdkit_compat.rdInfoTheory`.
This is an adjacent typed InfoTheory capability motivated by the still-open
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) tracking issue.
RDKit already implements masks; this adds the missing chematic implementation
and does not implement nanobind or resolve the whole upstream issue.

```python
from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory

ranker = rdInfoTheory.InfoBitRanker(4, 2)
ranker.SetMaskBits(bit for bit in [0, 2])
for bits, label in [([0, 1], 0), ([2, 3], 1)]:
    fp = ExplicitBitVect(4)
    for bit in bits:
        fp.SetBit(bit)
    ranker.AccumulateVotes(fp, label)

rows = ranker.GetTopN(2)
assert rows[:, 0].tolist() == [0.0, 2.0]
ranker.SetMaskBits(range(4))  # Restores eligibility, not previously skipped votes.
all_rows = ranker.GetTopN(4)
by_bit = {int(row[0]): row[2:].tolist() for row in all_rows}
assert by_bit[1] == [0.0, 0.0]
```

## Mask and vote semantics

The initial state includes every feature. A mask is a replacement set, not a
union with the previous one. Duplicate IDs collapse to one feature. Input
order does not change the existing deterministic tie order (ascending bit ID).
An empty iterable selects no features; pass `range(nBits)` to restore them all.
`None` is not a reset value.

The mask applies both to ranking and to **future on-votes**. Each valid record
still increments its class total. Votes for excluded bits are skipped, not
saved for later. Replacing the mask does not erase existing counts, restore
skipped votes, or recompute observations. Changing masks during accumulation
therefore measures different feature observation windows. To compare all
features over one common dataset, choose the mask before collecting votes or
start a fresh ranker and replay the dataset.

This follows the reference behavior in
[RDKit's InfoBitRanker source](https://github.com/rdkit/rdkit/blob/7f3d37d16b71ff8c8a6b0de82c12a8b371e0a8df/Code/ML/InfoTheory/InfoBitRanker.cpp).
The complete current vector is still checked before counts change; a mask
does not hide invalid, duplicate or out-of-range on-bit IDs.

## Input, errors and limits

`maskBits` is `Iterable[SupportsIndex]`, accepting generators and Python/NumPy
integer IDs in `[0, nBits)`. Booleans, floats, text and packed bytes are refused.
At most `nBits` raw entries are consumed successfully; a further entry raises
`ValueError`, including when entries are duplicates. This bounds infinite
iterators. The existing [ranker resource limits](info-bit-ranker.md) still apply.

Validation completes before replacing the mask. Invalid IDs, iteration errors
and oversized input leave the previous mask and counts unchanged. Exceptions
from caller iteration propagate. Side effects in caller code are outside this
guarantee. Only an owned immutable set is retained.

`GetTopN(num)` requires `0 <= num <= number_of_unique_selected_bits` and at
least one accumulated record. It returns a new float64 array with `nClasses+2`
columns. With an empty mask, only `GetTopN(0)` is allowed. Selecting a small
mask does not reduce the existing count table's memory allocation; ranking
still computes scores over that bounded table. No performance claim is made.

## Verification and reference differences

Tests cover both metrics, mask timing/replacement, duplicate IDs, source and
output ownership, complete-vector validation, failed iterators and bounded
consumption. Python/NumPy implements this feature; RDKit is not a runtime
dependency. The existing native extension is an import/legacy regression
fixture, not a newly built native implementation.

The seeded reference run records successful full-mask comparisons separately
from RDKit results with duplicate/missing eligible IDs. All candidate results
also pass independent per-bit count and scalar entropy/exact-rational
chi-square checks. Two isolated diagnostics record RDKit's zero-top signal
exit with a nonempty mask and partial mask replacement after an invalid ID;
chematic returns the existing empty-array result and preserves its prior mask.
These diagnostics are not counted as successful parity tests. Tie ordering,
raw-item bounds and stricter input errors remain explicit compatibility
limits. Details and reproducible commands are in
`validation/results/info-ranker-masks-2026-10-18/`.
