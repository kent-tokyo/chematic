# Ranked screening enrichment factors

RDKit [issue #7981](https://github.com/rdkit/rdkit/issues/7981) remains open
as checked on 2026-10-14. Its example asks for the top 90% of five ranked
compounds, so the selected count is five. Locally installed RDKit 2026.09.1
returns an empty list instead of one enrichment factor. The same defect occurs
at 100%, including a one-compound library.

chematic now provides `enrichment_factors` in Rust (`chematic_chem`) and Python.
The input is a binary activity label for every ranked compound, best first.
The function does not sort scores or break ties. The caller must use a ranking
and tie policy appropriate to its screening experiment.

```python
import chematic

ranked_actives = [True, True, False, False, False]
chematic.enrichment_factors(ranked_actives, [0.2, 0.9, 1.0])
# [2.5, 1.0, 1.0]
```

For `N` ranked compounds, `A` total actives and fraction `f`, the calculation is:

```text
k = ceil(N * f)
EF = (actives among the first k compounds / k) / (A / N)
```

The actual selected count supplies the denominator. For example, fraction
0.21 of five compounds selects two; with both actives at the top, EF is 2.5.
It is not divided by the nominal 0.21 fraction. The `ceil` operation uses the
floating-point product `N*f`, as RDKit does. A positive fraction smaller than
one rank selects one entry. A fraction of 1 returns 1 for any nonempty library
with at least one active.

Output has one value per requested fraction, preserving order and duplicates.
Fractions need not be sorted, and multiple fractions may select the same count.
The implementation scans the ranked labels in cutoff order, with
`O(N + F log F)` work and `O(F)` auxiliary memory for `F` fractions.

## Error and compatibility policy

Rust returns `Result<Vec<f64>, ScreeningError>`; Python maps input errors to
`ValueError`:

- Empty rankings or empty fraction lists are errors.
- Fractions must be finite and strictly greater than zero and at most one.
- A library with no actives is an error because EF has a zero denominator.
  A library with some actives but none in the selected prefix has EF zero.

Python requires boolean labels (including NumPy booleans). Integer labels,
weights and score values are not converted to activity truth values; convert
validated binary labels to booleans before calling.

This is a native API using the enrichment definition and RDKit's ceiling
convention. It intentionally differs from `CalcEnrichment` at broken boundary
conditions and undefined inputs: RDKit accepts zero fractions and returns
zeros for no actives. No general RDKit scoring compatibility is claimed.
ROC, RIE, BEDROC, score sorting and statistical uncertainty are outside this API.

## Validation

The original ranking is retained in
`validation/rdkit_issues/screening/enrichment-7981.json`.
`crates/chematic-chem/tests/screening_enrichment.rs` adds eight Rust tests,
including exhaustive binary rankings up to eight compounds against an
integer-count oracle, as well as boundaries, repeated and unsorted fractions,
rounding, invalid inputs, and the smallest positive float. The complete local
`chematic-chem` suite passed 1,149 tests with 5 ignored and no failures.

The Python suite `crates/chematic-py/tests/test_screening_enrichment.py` adds
25 cases. All passed alongside seven existing bulk/similarity tests using the
locally compiled source extension. Clippy for `chematic-chem` and `chematic-py`
passed with warnings denied.

`scripts/validate_rdkit_enrichment.py` ran 3,048 ordinary single-cutoff
comparisons against both live RDKit 2026.09.1 and an independent count formula.
The maximum absolute error was `2.220446049250313e-16`. Four separate boundary
cases validate the corrected behavior, and five undefined-input probes record
explicit errors. The source, fixture, RDKit function and extension hashes are
in `validation/results/rdkit-enrichment-2026-10-14.json`.

These checks concern local uncommitted source, not a published package or an
upstream RDKit change.
