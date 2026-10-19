# MOL2 parent/H row ordering — 2026-10-19

Target: [Open Babel #2733](https://github.com/openbabel/openbabel/issues/2733),
verified OPEN on 2026-10-19. Adds an explicit parent/H ordering option rather
than changing the default atom order or generating hydrogens. The native parser
and existing graph writers remain unchanged.

## Checks

- 134 tests passed: 33 new ordering tests plus 101 existing MOL2 Kekule,
  stream/RMSD, section-order and source-preservation tests.
- 14 molecule inputs × 12 atom orders = 168 independent cross-reads. Each
  source MOL2 is authored from an embedded explicit-H RDKit graph; references
  are renumbered together when creating each source permutation.
- RDKit 2026.09.1 validates source identity against the input molecule and
  verifies output element/bond correspondence, explicit-H parent grouping,
  exact coordinates and Standard InChI. Includes opposite tetrahedral and E/Z
  examples; no source or failed comparison was discarded.
- Open Babel 3.2.1 rewrites all 168 grouped outputs. RDKit then verifies the
  mapped graph and Standard InChI against the grouped output. Maximum coordinate
  error is `4.987237599995975e-05` Angstrom, within the expected four-decimal
  rounding. The direct RDKit ordering comparison has zero coordinate error.
- Existing Open Babel `-ismi -:CCO -h -omol2` produces heavy-first rows. Its
  original/grouped text and the mapping are retained in `evidence.json`.
- Mypy on the changed module (`--follow-imports=skip`), Black on the three
  new Python files, Python 3.9 grammar parsing, `git diff --check` and the
  MkDocs site build passed. The committed implementation hash matches the
  module used by the oracle.

The oracle does not compare partial-charge reassignment by external writers;
the source-edit tests and mapped native readback verify that the transformation
retains original partial-charge tokens and separate formal-charge attributes.
This is not general atom-type/valence correctness or NMR naming validation.

## Local runtime boundary

The candidate Python package modules were staged with the published chematic
1.1.8 CPython 3.13 macOS arm64 native reader, SHA-256
`ed5a04679987ee9b0a0a99daadca9f7a907b7af00e9e0c4180eaf49b9dbfcc55`.
The temporary initializer exported native symbols and the new public API. It
also included the unchanged current-source `convert_format`, `_FORMAT_ALIASES`
and `_CONVERT_MAX_INPUT_BYTES` extracted with Python AST for the existing
conversion tests. The full current initializer/native build is a CI integration
check, not a locally verified or published wheel.

An initial broad test attempt had two failures because the minimal temporary
initializer omitted `convert_format`. After including its unchanged source,
all 134 tests passed. The initial 31-test run also had a pytest-cache disk-space
warning; the final run disabled that cache and had no warnings. Earlier
regenerable document build directories created by this task were cleaned up.

## Reproduce

With a built candidate package installed (RDKit is only needed for the oracle):

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -p no:cacheprovider -q \
  crates/chematic-py/tests/test_mol2_hydrogen_order.py \
  crates/chematic-py/tests/test_mol2_kekule.py \
  crates/chematic-py/tests/test_mol2_stream_rmsd.py \
  crates/chematic-py/tests/test_mol2_section_order.py \
  crates/chematic-py/tests/test_mol2_preservation.py
python3 scripts/validate_mol2_hydrogen_order.py \
  --obabel /opt/homebrew/bin/obabel --output /tmp/mol2-hydrogen-order.json
```

The script fixes the embedding and permutation seeds. Every record includes
source/output SHA-256 and output-to-source row mapping. It does not reproduce
the issue's aromatic kekulization diagnostic, which is independent of row order.
