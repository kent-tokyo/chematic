# mmCIF chain record evidence — 2026-10-18

Related open issue: https://github.com/openbabel/openbabel/issues/2458

The existing native reader already retains full residue identities. This adds
an explicit eager Python view keyed by model number and label_asym_id; it
neither changes the default aggregate reader nor equates chains with connected
molecules. Native first-atom_site-loop scope and numeric defaults remain.

- Base: `f99c7b8ab33098a2e1fd9fd568aee370bf79ab38`.
- `pytest.log`: 238 tests pass, including 51 new cases. One new case checks
  all 720 permutations of six atoms for complete, ordered partitions and exact
  correspondence to source rows. Existing native residue/compatibility tests
  also pass.
- `oracle.py`, `oracle.log`, `oracle-summary.json`: 96 seeded synthetic cases
  with shuffled columns/rows, quoted identifiers and wrapped rows, independently
  parsed with Gemmi 0.7.3. All 384 chain/model partitions and 4,992 atom/residue
  records agree; maximum observed coordinate difference is zero.
- `chains.cif`: original eight-atom fixture with four model/label-chain pairs,
  shared author aliases, long component names, alternate locations, distinct
  waters and an opaque extra column. No upstream attachment was copied.
- `openbabel.stdout`, `openbabel.stderr`: local Open Babel 3.2.1 reports one
  aggregate record with all eight atoms. This demonstrates default granularity,
  not an assertion that aggregation is always wrong or that Open Babel lacks
  connected-component splitting. Command is recorded in the oracle.
- `typing-valid.log`, `typing-invalid.log`: strict helper/public consumer passes;
  the negative consumer exits 1 with four expected argument errors. Unrelated
  imported legacy stubs are silent; package-wide strict typing is not claimed.
- `docs-and-typing.log`, `mkdocs.log`: Python 3.9 implementation grammar,
  executed guide, strict build and five rendered API anchors verified.
- `source-gates.log`: source stub, manifest, native binding inventory, unsafe,
  refactor-budget and naming gates pass. No wheel was supplied to the stub gate.
- `native-source-comparison.log`: shared CIF/mmCIF source, residue parser and
  residue binding paths, plus the two record-to-dict helpers, match the verified
  old native fixture source. No current-main native rebuild or whole-extension
  equivalence is claimed.
- `provenance.json`: source/evidence/runtime hashes and environment.

Inputs are capped at 8 MiB UTF-8, 50,000 parsed atoms, 4,096 model/chain pairs
and 8,192 bytes per native line; callers may only lower these limits. Missing
label chain IDs, invalid limits and native ambiguous/malformed inputs fail
without a partial result. Full source text preserves unknown values/categories,
cell information, additional loops/blocks, comments and numeric spellings;
only the native first atom-site loop is projected into typed records.

No bonded Mol supplier, biological assembly expansion, atom-site rewrite,
streaming parser, Rust/WASM interface, real-structure corpus qualification,
performance or published-wheel validation is part of this change.

## Reproduce

Place a compatible native extension beside this revision's Python package and
add its parent directory to PYTHONPATH. Add an independent Gemmi 0.7.3 install
to PYTHONPATH only for the oracle; unit tests do not import Gemmi/Open Babel.

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3.13 -m pytest \
  crates/chematic-py/tests/test_mmcif_chains.py \
  crates/chematic-py/tests/test_mmcif_residues.py \
  crates/chematic-py/tests/test_rdkit_compat.py -q -p no:cacheprovider
PYTHONDONTWRITEBYTECODE=1 python3.13 validation/results/mmcif-chain-records-2026-10-18/oracle.py
python3.13 -m mkdocs build --strict --site-dir /private/tmp/chematic-chains-site
```

Oracle control requires `/opt/homebrew/bin/obabel`, the installed path on this
host. Strict mypy uses `follow_imports=silent`, `follow_imports_for_stubs=True`,
`ignore_missing_imports=True`, and MYPYPATH pointing to the candidate package.
Check the helper with `typing_valid.py`; check `typing_invalid.py` separately
and expect exactly four errors. The typed consumer examples are not input files
to execute, and do not require an actual `input.cif` for type checking.
