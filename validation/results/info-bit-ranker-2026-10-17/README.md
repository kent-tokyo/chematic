# Typed entropy bit ranking — 2026-10-17

This Python/NumPy-only candidate implements the ENTROPY InfoBitRanker input
item in RDKit #9601. Explicit/sparse binary vectors share a structural
GetNumBits/GetOnBits protocol. It does not implement biased entropy,
chi-squared, feature masks, correlation matrices or file output. Invalid
votes are staged before count updates; counts do not retain vector objects.

Evidence:

- Clean baseline Python source at f20304f4 has no rdInfoTheory namespace.
- 207 selected tests pass, including 45 new ranker tests, covering entropy,
  multiclass counts, vector storage/order, state preservation after invalid
  input or failing iteration, result isolation, numeric protocols and limits.
- 768 live RDKit 2026.09.1 ranker comparisons use actual RDKit ExplicitBitVect
  and SparseBitVect objects and 10,560 paired votes. All per-bit/class counts
  match; maximum absolute information-gain difference is
  1.2212453270876722e-15. Cases span 1/4/16/64 bits, 2/3/5 classes and
  1/3/11/40 observations. This is a synthetic feature corpus, not a chemistry
  accuracy or performance benchmark.
- 1,536 order-reversal/data-duplication checks pass. Candidate exact-score
  ties use ascending bit IDs; 462 whole rank sequences differ from RDKit's
  heap ordering. Sorted score sequences still agree within 1e-12. No exact
  tie/rank-order parity is claimed.
- Isolated helper strict mypy passes, public valid examples (including NumPy
  integer scalars and vector IDs) type-check, and four invalid type examples
  produce expected arg-type errors. No package-wide strict typing claim.
- Refactor/source-stub/native-binding-inventory/cross-binding/unsafe/naming
  gates and strict MkDocs pass. Both production Python files parse with Python
  3.9 grammar; runtime validation used Python 3.13.

The reused native extension from source commit
`a56bc1353bd70eb5381dce38f818c0b78755d027` is import scaffolding. The new
ranker calls Python/NumPy only. It was not rebuilt, and no wheel build,
publication or native Rust ranking API is claimed. Hashes, versions and
upstream issue state are recorded separately. Native binding inventory is
unchanged.

Reproduce with candidate Python package source and a compatible extension:

```sh
export PYTHONPATH=/private/tmp/chematic-info-bit-ranker-python:/private/tmp/chematic-sparse-tanimoto-rdkit
export PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1
python3.13 -m pytest crates/chematic-py/tests/test_info_bit_ranker.py crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
python3.13 validation/validate_info_bit_ranker.py
MYPYPATH=crates/chematic-py/python python3.13 -m mypy --follow-imports=silent --ignore-missing-imports --no-incremental validation/typing/info_bit_ranker_valid.py
# Expected exit 1, exactly four arg-type errors:
MYPYPATH=crates/chematic-py/python python3.13 -m mypy --follow-imports=silent --ignore-missing-imports --no-incremental validation/typing/info_bit_ranker_invalid.py
# Isolate the helper from legacy package stubs when testing strictness:
mkdir -p /private/tmp/chematic-info-bit-ranker-typecheck
cp crates/chematic-py/python/chematic/_info_theory.py /private/tmp/chematic-info-bit-ranker-typecheck/info_theory.py
python3.13 -m mypy --strict --follow-imports=silent --ignore-missing-imports --no-incremental /private/tmp/chematic-info-bit-ranker-typecheck/info_theory.py
```
