# InfoTheory ranking table export

`InfoBitRanker.WriteTopBitsToFile(fileName)` writes the last successful
`GetTopN` result as an ASCII table. This fills a missing chematic capability
alongside the typed InfoTheory work tracked in
[RDKit #9601](https://github.com/rdkit/rdkit/issues/9601), verified open on
2026-10-19. RDKit already has this method; this is not an upstream feature
invention or a resolution of the entire wrapper issue.

```python
from pathlib import Path
from chematic.rdkit_compat import ExplicitBitVect, rdInfoTheory as I

ranker = I.InfoBitRanker(4, 2)
for label, on_bits in [(0, [0, 1]), (1, [1, 2])]:
    fp = ExplicitBitVect(4)
    for bit in on_bits:
        fp.SetBit(bit)
    ranker.AccumulateVotes(fp, label)

rows = ranker.GetTopN(4)
ranker.WriteTopBitsToFile(Path("ranked-bits.txt"))
```

The columns are `Bit`, `InfoContent`, then `class0`, `class1`, etc. IDs and
counts are decimal integers; scores use five significant digits. The
`InfoContent` heading follows the upstream format even for CHISQUARE, where
the score is a Pearson statistic. The existing
[score definitions and tie policy](info-bit-ranker.md) apply. In particular,
equal scores use ascending bit IDs in chematic; RDKit's heap order differs.
Use `GetTopN` arrays when full numeric precision is required.

## Snapshot lifecycle

Every successful `GetTopN`, including `GetTopN(0)`, replaces an independent
snapshot for export. Editing its returned array, accumulating more votes,
or changing masks/bias preferences does not change that snapshot. Call
`GetTopN` again to export an updated ranking. A failed ranking preserves the
previous successful snapshot. A zero-row snapshot produces the header alone.

Before the first successful ranking, export raises `ValueError` without
opening the destination. Ranker state is unchanged by export, including on
failure. The snapshot adds at most approximately 8.6 MiB of float64 storage
under the existing ranker dimension limits; text formatting also uses bounded
temporary storage proportional to the selected rows and class count.

## Paths and file errors

`fileName` accepts `str`, `bytes`, `os.PathLike[str]` and
`os.PathLike[bytes]`, including `pathlib.Path`. It is a filesystem path,
not a stream or descriptor. Wrong types raise `TypeError`; invalid paths
such as embedded NUL raise `ValueError`; open/write failures raise `OSError`.
Exceptions from caller-defined `__fspath__` methods propagate.

The complete table is formatted before the file is opened. Existing files
are overwritten; missing parent directories are not created. Newlines are
always LF. As with an ordinary text file write, an I/O failure after opening
can leave a partial destination. The method does not provide atomic file
replacement or synchronization with simultaneous ranker mutations.

## Validation boundary

The implementation is Python/NumPy and requires no RDKit runtime. Tests cover
literal file bytes, typed/non-ASCII paths, small-score scientific notation,
multiclass headers, snapshot ownership and failed ranking/I/O behavior.

`scripts/validation/validate_info_ranker_export.py` runs 64 datasets through four metrics
and explicit/sparse source variants: 512 paired exports. All 1,024 per-runtime
format checks and all 1,024 saved-snapshot checks pass against RDKit 2026.09.1.
Whole files match in 98 cases. Other files can differ because their ranking
orders or score values differ, not because the formatting contract differs.
Among 508 cases with complete reference bit-ID permutations, per-class counts
match exactly and scores differ by at most `1.78e-14`; eight formatted score
fields differ. Four reference results contain duplicate/missing bit IDs and
are recorded separately with their rows, not counted as matching rankings.
No case is dropped from the report.

Raw evidence and the local package-import boundary are in
`validation/results/info-ranker-export-2026-10-19/`. These measurements do not
establish exact ranking parity, a native rebuild, or published-wheel support.
