# Typed InfoBitRanker masks — 2026-10-18

Adds the missing Python/NumPy `SetMaskBits` implementation for ENTROPY and
CHISQUARE. This is an adjacent typed InfoTheory capability motivated by
RDKit #9601, not a fix for a reported upstream mask defect or the whole
nanobind tracking issue. RDKit already has this API. The new code is original;
reference semantics were inspected in the official repository and confirmed
against live RDKit 2026.09.1.

## Behavior

Masks replace the previous selection. Duplicate IDs collapse; empty selects
nothing. Masked-out future on-votes are skipped but class totals always grow.
Past counts remain, and later selection cannot restore skipped votes.
`GetTopN` ranks only the currently eligible features. Source iteration is
staged and bounded to nBits raw items; failures preserve the previous mask.
The public method takes `Iterable[SupportsIndex]`, with stricter runtime
rejection of bool/text/packed bytes and a deterministic bit-ID tie rule.

## Evidence

- `before.log`: all 26 new mask cases fail on the former unsupported method.
- `pytest.log`: 251 passing tests, including the 26 new mask cases and existing
  entropy, chi-square and compatibility regressions.
- `oracle.py` / `oracle-summary.json`: seeded 96 datasets × 2 metrics × 4 input
  sources, with four changing mask phases. 2,064 valid full-mask RDKit array
  comparisons pass (25,192 rows), with exact counts/IDs and maximum score
  difference 7.105427357601002e-15. Another 768 empty-mask shape checks pass.
- All 2,304 nonempty candidate arrays pass independent Python integer-count,
  scalar log2 entropy and exact-Fraction Pearson checks.
- `reference-anomalies.json`: 240 additional reference runs return duplicate or
  missing eligible IDs (60 distinct dataset/metric/phase cases across four
  sources). They are excluded from successful comparison counts and have
  committed seed/context/output evidence. Replay the script to regenerate
  full phase inputs. No full-ranking parity claim is made.
- `separate-diagnostics.json`: two isolated RDKit controls, with core dumps
  disabled: a zero-top query with a nonempty mask exits by signal (-11);
  an invalid replacement mask raises IndexError after replacing the prior mask
  with a partial one. The candidate keeps its existing safe zero-top behavior
  and atomically preserves its previous mask. These are not parity cases.
- `typing-valid.log`, `typing-invalid.log`, `docs-and-typing.log`: strict checking
  of the modified helper and a public consumer, exactly four expected negative
  errors, Python 3.9 AST compatibility and executed guide code. Legacy imports
  are silent; this is not package-wide strict verification.
- `source-gates.log`, `mkdocs.log`: source gates, Black and strict docs passed.

The verified older native extension is an import/legacy fixture only. All new
mask behavior is Python/NumPy; no rebuilt native binary, released wheel or
Rust/WASM extension claim. No RDKit dependency at runtime.

## Reproduce

Prepare a candidate package from `crates/chematic-py/python/chematic` with a
compatible local native extension; original exact paths/hashes are in
`provenance.json`. Then run from the repository root:

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 PYTHONPATH=/private/tmp/chematic-info-ranker-masks-python python3.13 -m pytest crates/chematic-py/tests/test_info_ranker_masks.py crates/chematic-py/tests/test_info_bit_ranker.py crates/chematic-py/tests/test_chi_square_bit_ranker.py crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/private/tmp/chematic-info-ranker-masks-python:/private/tmp/chematic-sparse-tanimoto-rdkit python3.13 validation/results/info-ranker-masks-2026-10-18/oracle.py
MYPYPATH=/private/tmp/chematic-info-ranker-masks-python python3.13 -m mypy --config-file validation/results/info-ranker-masks-2026-10-18/mypy.ini crates/chematic-py/python/chematic/_info_theory.py validation/results/info-ranker-masks-2026-10-18/typing_valid.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-info-ranker-masks-site
```

The negative typing consumer intentionally exits 1. Host-specific package paths
may be replaced with compatible local installs; all comparator versions and
module hashes are recorded. Hosted CI remains a separate integration gate.
