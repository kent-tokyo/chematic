# Class-biased InfoTheory evidence — 2026-10-18

This is a Python/NumPy compatibility feature addition associated with RDKit
#9601's typed explicit/sparse ranker inputs. RDKit already implements both
biased metrics; this does not fix or close its tracking issue.

- Candidate base: `540086f978612df1e21baa731c0884dea7441f46`.
- `pytest.log`: 250 tests pass, including 27 new test cases. One case checks
  all 2,916 small-table/preference combinations with exact Fraction arithmetic.
- `oracle.py` / `oracle-summary.json` / `oracle.log`: seeded live RDKit
  2026.09.1 explicit/sparse per-bit score/count comparisons. All 1,152 rankers
  return complete permutations, counts agree exactly, maximum score difference
  1.42e-14. Independent Fraction/direct-mutual-information calculation checks
  10,368 scores, maximum error 3.56e-15. No reference anomalies in this corpus.
- `typing-valid.log`: strict helper and public consumer pass. Unrelated legacy
  imported stubs use silent imports; this is not package-wide strict typing.
- `typing-invalid.log`: four deliberate public usage errors; exit 1 expected.
- `docs-and-typing.log`: Python 3.9 grammar, guide example, rendered guide and
  positive/negative checker status verified.
- `mkdocs.log`: strict documentation build passes.
- `source-gates.log`: source stub, cross-binding manifest, binding inventory,
  unsafe surface, refactor budget and product spelling checks pass. No wheel
  was supplied to the stub gate.
- `provenance.json`: hashes and environment; old native extension is reused
  only as an import fixture. No current-main native build, published wheel,
  Rust/WASM feature or performance validation is claimed.

The bias rule is verified against RDKit source revision
`7f3d37d16b71ff8c8a6b0de82c12a8b371e0a8df`,
`Code/ML/InfoTheory/InfoBitRanker.cpp` blob
`81037caa45399213eaf665a405d67268e88c5523`. Its public setter requires strictly
fewer than nClasses unique IDs, including refusing a full-class preference;
empty/default preference leaves all base scores intact. chematic uses typed
inputs, bounded extraction, explicit Python errors and staged replacement.
Rejected scores are zeroed without removing their counts or bit IDs. Empty
classes use zero fractions. ENTROPY/CHISQUARE ignore valid stored preferences.
Masks remain unsupported.

The deterministic chematic ranking differs from the upstream heap order in
916 comparisons (ties and floating point); only per-bit scores/counts are
compared. No exact-order or all-input compatibility claim is made.

## Reproduction

Use a copy of this source Python package with a compatible native extension
for package import; add that candidate directory to PYTHONPATH. For the
oracle also add an independent RDKit 2026.09.1 installation to PYTHONPATH.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_biased_bit_ranker.py \
  crates/chematic-py/tests/test_info_bit_ranker.py \
  crates/chematic-py/tests/test_chi_square_bit_ranker.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 python3.13 validation/results/biased-bit-ranking-2026-10-18/oracle.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-biased-site
```

Use strict mypy with `follow_imports=silent`, `follow_imports_for_stubs=True`,
`ignore_missing_imports=True`, and MYPYPATH pointing to the candidate package.
Check `_info_theory.py` plus `typing_valid.py`, then `typing_invalid.py` separately
and expect exactly four errors. No RDKit import is needed for the unit tests.
