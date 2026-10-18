# Typed supplied-candidate score callbacks — 2026-10-18

Adds `rdkit_pick_tautomer_scored` alongside the unchanged default native picker.
Related to RDKit #9601's typed iterable/callback and ownership work; RDKit
already supports custom PickCanonical scoring. This is a missing chematic
capability, not an upstream algorithm defect or whole-issue resolution.

## Results

- `baseline-gap.log`: base-source capability inspection; the native picker
  only accepts candidates and max_candidates. No fabricated before-crash test.
- `pytest.log`: 221 tests pass, including 41 new score-callback cases and 180
  existing pickle/tautomer/compatibility cases. Inputs are snapshotted before
  advancing iterators; callbacks get separate views after all validation.
  Tests cover property types, reused generators, GC, external source edits,
  iterator/callback exceptions, count/atom/bond/payload limits, integer score
  bounds, empty molecules, and stereo/isotope/map preservation.
- `oracle.py` / `oracle-summary.json`: 128 seeded trials, 8 acyclic supplied
  families, 4 pure scoring policies and 3 candidate source forms. All 1,536
  selected-index/graph comparisons match RDKit 2026.09.1. Both sides make
  exactly 3,456 score calls in candidate order. The reference always retains
  its input list; there is no RDKit temporary-generator lifetime claim.
- `separate-diagnostics.json`: one deliberately different lexical/stereo tie
  policy, excluded from parity. Newer RDKit prefers specified stereochemistry;
  this API retains the existing native picker's lexical policy.
- Strict helper/public-consumer mypy, exactly four negative type errors,
  Python 3.9 AST, executed guide code, source gates/Black and strict MkDocs pass.

`native-source-comparison.log` verifies the native codec, properties,
canonical wrapper, bond-table and tautomer-score files against the fixture.
Current-main macrocycle aromaticity and unrelated scalar-cache additions are
also recorded; the new reference packet is acyclic. The extension is an older
hash-verified backend, not a freshly built whole-current-main artifact.
No native chemistry algorithm was changed. The callback can choose arbitrary
supplied graphs; it does not check tautomer-family membership or perform final
stereo cleanup. Snapshot copies follow the versioned native state's contract,
not a promise to preserve input source spelling or derived caches.

## Reproduce

Prepare a candidate package from `crates/chematic-py/python/chematic` with a
compatible native extension; exact host paths/hashes are in provenance.json.
Then run from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-tautomer-score-python python3.13 -m pytest crates/chematic-py/tests/test_tautomer_scored.py crates/chematic-py/tests/test_pickle.py crates/chematic-py/tests/test_rdkit_tautomer.py crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/chematic-tautomer-score-python:/private/tmp/chematic-sparse-tanimoto-rdkit python3.13 validation/results/tautomer-score-callback-2026-10-18/oracle.py
MYPYPATH=/private/tmp/chematic-tautomer-score-python python3.13 -m mypy --config-file validation/results/tautomer-score-callback-2026-10-18/mypy.ini crates/chematic-py/python/chematic/_scored_tautomer.py validation/results/tautomer-score-callback-2026-10-18/typing_valid.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-tautomer-score-site
```

The negative typing example intentionally exits 1. Host-specific package
paths may be replaced with compatible local installs. Caller callbacks can
retain views or perform arbitrary work; no deadline or total-memory guarantee
is made. Hosted CI and released artifacts are independent integration gates.
