# Typed binary bit-pair co-occurrence — 2026-10-18

Related open tracking issue: https://github.com/rdkit/rdkit/issues/9601

RDKit already supplies BitCorrMatGenerator. chematic has a typed InfoTheory
ranker, but lacks this adjacent operation. This adds a checked Python/NumPy
implementation accepting explicit/sparse structural binary sources, with
selected-order pair counts, reset semantics, staged validation and owned arrays.
It is related to the issue's typed InfoTheory input family; no upstream
co-occurrence algorithm defect or whole-issue completion is asserted.

- Base: `f99c7b8ab33098a2e1fd9fd568aee370bf79ab38`.
- `pytest.log`: 262 tests pass (37 new, 44 entropy, 19 chi-square, 162
  compatibility). Covers exact small binary patterns and selected permutations,
  state failures, source lifetime, immutable ownership and actual 1024-selected,
  65536-raw-ID and production pair-work boundaries. The record cap/reset test
  lowers its budget; one million voting calls are not executed.
- `oracle.py`, `oracle-summary.json`, `oracle.log`: RDKit 2026.09.1 and independent
  integer membership sums on 128 datasets, four source representations, 512
  final matrices / 22240 pair counts / 12056 paired votes. Count error is zero.
  All 1536 prefix checks and 1536 reversal/duplication/reset checks pass.
- `reference-diagnostics.json`: two source-profile diagnostics, width changes
  without reset and a high signed sparse view, recorded separately from ordinary
  parity. No ordinary failures were discarded.
- `typing-valid.log`, `typing-invalid.log`, `mypy.ini`: strict helper/public
  consumer pass and four expected negative errors, including wrong result type.
  Legacy imports/stubs are silent, NumPy typing is followed normally; no
  package-wide strict typing claim.
- `docs-and-typing.log`, `mkdocs.log`: Python 3.9 grammar, executed guide, strict
  docs build and four rendered API anchors pass.
- `source-gates.log`: source stubs, manifest, binding inventory, unsafe surface,
  refactor budget and canonical names pass. No wheel supplied.
- `provenance.json`: source/runtime/reference/evidence hashes and issue snapshot.

The official behavioral reference is
https://www.rdkit.org/docs/source/rdkit.ML.InfoTheory.rdInfoTheory.html
No RDKit source implementation was copied. The algorithm counts each selected
pair's joint membership directly; output is not normalized Pearson correlation,
a probability, a p-value or a chemical fingerprint.

Hard caps per configuration: 2..1024 unique selected IDs, 65536 raw IDs per
record, one million records and 50 million actual selected-pair increments.
Widths are consistent in 1..4294967295 and must contain all selected IDs.
Selected/raw IDs exclude bool, are unique and fit the stated unsigned range.
No allocation scales with declared bit width. Each successful setter resets
counts, width and budgets. Validation/limit failures preserve previous state;
caller-method side effects and concurrent mutations are outside that guarantee.

The new calculation is entirely Python/NumPy. A hash-verified old extension is
reused only to import the package and run unchanged legacy regressions. No fresh
current-main native build, whole-extension equivalence, published wheel,
performance, all-input RDKit compatibility or normalized correlation claim.

## Reproduce

Put this revision's Python package with a compatible extension on PYTHONPATH;
add independent RDKit 2026.09.1 for the oracle. Tests do not need RDKit.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_bit_cooccurrence.py \
  crates/chematic-py/tests/test_info_bit_ranker.py \
  crates/chematic-py/tests/test_chi_square_bit_ranker.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 python3.13 validation/results/bit-cooccurrence-2026-10-18/oracle.py
python3.13 -m mypy --config-file validation/results/bit-cooccurrence-2026-10-18/mypy.ini \
  crates/chematic-py/python/chematic/_bit_correlation.py \
  validation/results/bit-cooccurrence-2026-10-18/typing_valid.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-bit-correlation-site
```

For typing set MYPYPATH to the candidate package parent; check the negative
consumer separately, expecting exit 1 with four errors.
