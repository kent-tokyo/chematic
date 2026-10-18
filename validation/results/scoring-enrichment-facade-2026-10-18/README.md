# Ranked-record enrichment facade — 2026-10-18

Related open issue: https://github.com/rdkit/rdkit/issues/7981

The native enrichment calculator already handles the reported endpoint defect.
This change adds the missing typed `rdkit_compat.Scoring.CalcEnrichment` record
adapter: checked activity-column selection, binary labels, generators and NumPy
arrays. It preserves caller ranking and fraction order/duplicates.

- Base: `f99c7b8ab33098a2e1fd9fd568aee370bf79ab38`.
- `pytest.log`: 242 tests pass (55 new, 25 native enrichment, 162 compatibility).
  Includes every nonzero binary ranking of length 1..8 against exact rational
  counts, actual hard limits, caller exceptions and ownership.
- `oracle.py`, `oracle-summary.json`, `oracle.log`: RDKit 2026.09.1 reference;
  1,024 ordinary rankings through four input representations, 4,096 comparisons
  and 137,216 paired records. Maximum RDKit/independent-formula error is
  4.440892098500626e-16.
- `boundary-comparisons.json`: six endpoint/order/duplicate controls differ
  from upstream and agree with independent exact counts. They are excluded from
  ordinary parity successes; no failed ordinary comparisons were discarded.
- `typing-valid.log`, `typing-invalid.log`: strict helper/public consumer pass;
  negative consumer exits 1 with four expected errors. Imported legacy modules
  and stubs are silent; package-wide strict typing is not claimed.
- `docs-and-typing.log`, `mkdocs.log`: Python 3.9 source grammar, executed guide,
  strict docs build and three rendered API anchors pass.
- `source-gates.log`: source stubs, cross-binding manifest, binding inventory,
  unsafe surface, refactor budget and canonical naming pass. No wheel supplied.
- `native-source-comparison.log`: native screening implementation and EF Python
  binding match the verified fixture revision. This does not establish complete
  current-main extension equivalence or a fresh native build.
- `provenance.json`: source/evidence hashes, native fixture and issue snapshot.

Hard caps are 100000 records, 4096 fractions and 1024 columns per record.
Column indices support -1024..1023 and must exist in every row. Fractions must
be finite real numbers in (0, 1], excluding bool. Activity labels are bool or
real numeric 0/1, deliberately stricter than RDKit's arbitrary truthiness.
Empty rankings/fractions and no-active rankings refuse. Caller iterators are
bounded by consumed items, but individual caller methods are not preempted.
No sorting, tie policy, ROC/AUC/RIE/BEDROC, performance or wheel claim is made.

## Reproduce

Put this revision's Python package plus a compatible native extension on
PYTHONPATH. The oracle additionally needs an independent RDKit 2026.09.1 install;
the tests do not require RDKit. The original extension fixture SHA256 is in
provenance.json; the relevant native source comparison is recorded separately.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_scoring_enrichment.py \
  crates/chematic-py/tests/test_screening_enrichment.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 python3.13 validation/results/scoring-enrichment-facade-2026-10-18/oracle.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-scoring-site
```

Strict mypy uses strict=True, follow_imports=silent,
follow_imports_for_stubs=True, ignore_missing_imports=True and MYPYPATH pointing
to the candidate package. Check the helper with typing_valid.py; separately
check typing_invalid.py and expect four errors. Python/NumPy/reference versions
and exact oracle corpus/function hashes are captured in provenance and summary.
