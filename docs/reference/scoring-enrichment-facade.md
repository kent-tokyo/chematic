# Ranked-record enrichment scoring

`chematic.rdkit_compat.Scoring.CalcEnrichment(scores, col, fractions)` now
accepts the ranked-record calling pattern used by RDKit's scoring module.
The existing native [binary-label enrichment API](screening-enrichment.md)
already implements the corrected ceiling/full-library behavior associated
with open [RDKit #7981](https://github.com/rdkit/rdkit/issues/7981).
This adds the missing typed record adapter and Scoring namespace, without
claiming a new native algorithm or a fix for every scoring function.

```python
from chematic.rdkit_compat import Scoring

# Caller-ranked records: score, binary activity. No automatic sorting.
scores = [(0.95, 1), (0.80, 1), (0.70, 0), (0.60, 0), (0.50, 0)]
assert Scoring.CalcEnrichment(scores, col=1, fractions=[0.2, 0.9, 1.0]) == [2.5, 1.0, 1.0]
assert Scoring.CalcEnrichment((row for row in scores), -1, [1.0, 0.2, 0.2]) == [1.0, 2.5, 2.5]
```

This follows the [RDKit ranked scoring convention](https://www.rdkit.org/docs/source/rdkit.ML.Scoring.Scoring.html):
records are already ordered best first, and `col` selects their activity label.
For `N` records, `A` total actives and fraction `f`, select `k = ceil(N*f)`
records and calculate `(actives_in_top_k / k) / (A / N)`. The actual selected
count supplies the denominator. All records therefore give EF 1, including a
single active record, while a selected prefix with no actives gives EF 0.
The existing native calculation and its floating-point cutoff convention are
unchanged. Fractions are converted to Python floats before native calculation.

`Scoring` is imported from `chematic.rdkit_compat`, not a new
`rdkit.ML.Scoring` module. Only `CalcEnrichment` is supplied; ROC/AUC, RIE and
BEDROC are not added. No RDKit runtime dependency is introduced.

## Typed input and bounds

`scores` can be an iterable/generator of indexed records, lists/tuples or a
2D NumPy array. The exported `RankedScoreRow` protocol requires length and
integer indexing. Only the selected column is read; metadata and score values
in other columns are ignored, not converted or used to reorder records.
Python-style negative column indices are supported independently for each row.
Rows must be indexed records, not text, scalar labels or row generators.

Selected labels must be Python/NumPy bools or real numeric values exactly
0 or 1. Integer/float NumPy activity columns are accepted, including when a
homogeneous score array stores binary labels as floats. Strings, None,
complex values, NaN, infinity, other numeric labels and arbitrary truthy objects
are refused. This deliberately narrows RDKit's arbitrary truthiness to binary
activity; it does not interpret weights or multiclass labels.

Fractions are a nonempty iterable of finite real numbers in `(0, 1]`, with
Python/NumPy bools refused. Output preserves their order and duplicates,
including multiple fractions selecting the same record count. Empty rankings,
empty fractions and no-active rankings raise `ValueError`; EF for no actives
is undefined under the native profile. Wrong types raise `TypeError`.

Hard caps are 100,000 records, 4,096 requested fractions and 1,024 columns per
record, with column indices in `-1024..1023`. These are not overrideable.
Fraction/column validation precedes reading score rows. Iteration consumes at
most the configured cap plus one item to detect overflow; infinite valid-input
iterables therefore terminate with errors. Caller iteration/indexing exceptions
propagate unchanged. All binary labels/fractions are copied before native
calculation, and no source container is retained or changed. Results are fresh
independent lists. Arbitrary caller iterator methods are not preempted.

## Validation

Tests cover the reported 90% endpoint, full-library/single-record cases,
selected-count rounding, duplicate/unsorted fractions, caller ranking, binary
label types, generators/NumPy/negative indexing, exact resource boundaries,
iteration/index failures before native calculation, independent outputs and
public typing. Every nonzero binary ranking up to eight records is checked
against an exact rational count definition, including boundary cutoffs.

The live RDKit 2026.09.1 oracle compares 1,024 ordinary rankings through four
representations: 4,096 calls and 137,216 paired records. Maximum difference
from RDKit and the independent count formula is `4.45e-16`. Six separate
endpoint/duplicate/order cases reproduce upstream differences and validate
the corrected native results against the exact count formula. They are not
counted as parity successes. The new tests do not claim broad scoring-module
or all-input RDKit compatibility.

Evidence, reproducible oracle and provenance are under
`validation/results/scoring-enrichment-facade-2026-10-18/`. Local Python tests
reuse a verified existing extension; the native screening source and Python
binding are unchanged from that fixture. No fresh current-main native build,
published-wheel or performance validation was performed.

## Python reference

::: chematic._scoring_enrichment
