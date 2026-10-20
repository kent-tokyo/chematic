# RDKit-compatible SMILES stability guard — 2026-10-19

Target: [RDKit #9391](https://github.com/rdkit/rdkit/issues/9391), verified OPEN
on 2026-10-19. This change provides an opt-in Python refusal API for observable
legacy-writer non-idempotence. It does not fix stereo perception or establish
chemical identity preservation.

## Results

- 47 tests passed: 27 new guard tests, 13 existing RDKit SMILES tests, and 7
  existing RDKit writer tests.
- The reported ammonium example, explicitly tagged spiro-junction variant,
  and carbon analogue all alternate between two strings. The guard refuses
  all three and retains both candidates in the exception.
- The fourth, fully substituted carbon control is stable and accepted.
- All four cases' first and rewritten outputs match RDKit 2026.09.1 with
  legacy stereo perception enabled.
- Existing canonical-SMILES fixture: 1,131 inputs, 1,125 accepted as stable,
  1 input parser refusal and 5 initial writer refusals. No additional stability
  refusals on this corpus. All parsed inputs retained identical MolJSON before
  and after the check. These results are not a stereochemical identity audit.
- Mypy on the new module (`--follow-imports=skip`), Black on the three new
  Python files, and `git diff --check` passed.
- MkDocs site build completed successfully with the new navigation entry.

`evidence.json` records the complete corpus, rejection reasons, issue cases,
reference outputs and hashes. The corpus is the committed
`crates/chematic-smiles/tests/data/rdkit_canonical_smiles.tsv`; its historical
RDKit expected strings are not used as correctness assertions in this audit.

## Runtime boundary and reproduction

The new Python module was tested using the native extension from the published
CPython 3.13 macOS arm64 chematic 1.1.8 wheel. Native SHA-256:
`ed5a04679987ee9b0a0a99daadca9f7a907b7af00e9e0c4180eaf49b9dbfcc55`.
This is **not** a current-main native rebuild. No Rust code changes are included.

For this isolated test runtime, a temporary package initializer exported
`chematic.chematic.*`, `__version__` and the two new public names. It loaded the
unchanged native extension and the new `_rdkit_smiles_checked.py` module.
The current full package initializer was not exercised locally: unrelated
newer Python helpers require native APIs absent from the published wheel.
Fresh package build/integration remains a CI check. Initial exploratory results
using an older local extension were discarded; committed evidence uses the
1.1.8 extension hash above.

With a built candidate package installed, run:

```sh
PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 python3 -m pytest -q \
  crates/chematic-py/tests/test_rdkit_smiles_checked.py \
  crates/chematic-py/tests/test_rdkit_smiles.py \
  crates/chematic-py/tests/test_rdkit_writers.py
python3 scripts/validation/validate_rdkit_smiles_stability.py --reference \
  --output /tmp/rdkit-smiles-stability.json
```

The validation script's `--reference` option requires RDKit. Ordinary use of
the new API and the tests do not require RDKit. The script temporarily enables
legacy stereo perception for its reference measurements and restores the
previous option afterward.
