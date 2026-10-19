# Position-preserving Standard InChIKey batches

Base: `f3e05eb458bac201dffe1d31eb9e5c1aa2ef0f78`. Checked 2026-10-19.
Open Babel [#2231](https://github.com/openbabel/openbabel/issues/2231) remains
open: omitted failed conversions destroy input/output correspondence.
Chematic lacked a native-backed Python InChIKey batch report API.

`inchi_keys_checked` stages plain input strings before native work, then returns
one immutable indexed success/failure record per accepted input. Failure rows
retain source, status, code and bounded message. A successful read must regenerate
the EXACT input Standard InChI before key generation. This guards the existing
reader's native-to-Rust fallback and never uses approximate `Mol.inchikey`.
It is a conservative conversion profile, not a replacement InChI parser.

## Results

- **191 tests passed**, including **39 new tests** plus existing input-error,
  module-function and molecule tests.
- The **219-row** mixed packet retains all input positions and duplicates.
  Of **144 valid Standard InChI rows** (72 identifiers repeated), **142** produce
  keys exactly equal to RDKit **2026.09.1** `InchiToInchiKey(original_input)`.
  Two rows, duplicates of one valid source, decline with `identity_changed`.
  The full input and reason stay in `oracle-summary.json`; they are not erased
  from coverage or counted as successful key matches.
- Overall statuses: 142 ok, 72 invalid, four unsupported (including those two
  valid-source rows), one resource_limit. Successful-key mismatches: **zero**.
- Open Babel **3.2.1** outputs two keys for three InChI lines with a malformed
  middle record, even using `-e` to continue. The control is an original water /
  malformed / ethanol sequence; stderr and output lines are preserved.
- Scoped strict mypy checks the helper and valid root consumer; the invalid
  consumer must produce exactly three argument errors. Imported legacy stubs
  are silent, not a package-wide strict typing claim.

Tests exercise structural preflight before native work, iterator/count/aggregate
caps, per-record byte/numeric/formula caps including leading zeros, exact source
and duplicate retention, immutable results, typed native categories, no key on
changed identity, native failures/unavailable backend, invalid key shape, bounded
messages, propagation of MemoryError/unexpected exceptions and later valid rows
after a malformed input. A structurally invalid/oversized batch returns no partial
tuple; caller iterator effects are not rolled back.

## Artifact boundary

Local candidate Python is staged in `/private/tmp/chematic-inchi-batch-python`
with the v1.1.8 native extension from `/private/tmp/chematic-public-v118`.
`provenance.json` verifies the wheel SHA256 against live PyPI metadata and
verifies that the extension bytes equal that wheel member. Native parsing,
InChI generation and binding sources are unchanged between tag v1.1.8 and base.
No Rust source changed. This is a published native backend with candidate Python,
not an unmodified-wheel test, fresh native build or release of this feature.

An initial wider run using the old v1.1.7 fixture reported three unrelated missing
`Mol.cleanup` methods; that artifact mismatch is recorded in `fixture-v117.log`.
The final v1.1.8 fixture supports those APIs and the entire selected suite passes.

## Reproduction

```sh
PYTHONDONTWRITEBYTECODE=1 PYTEST_DISABLE_PLUGIN_AUTOLOAD=1 \
PYTHONPATH=/private/tmp/chematic-inchi-batch-python \
python3.13 -m pytest crates/chematic-py/tests/test_inchi_batch.py \
  crates/chematic-py/tests/test_input_errors.py \
  crates/chematic-py/tests/test_module_functions.py \
  crates/chematic-py/tests/test_mol.py -q -p no:cacheprovider

PYTHONDONTWRITEBYTECODE=1 \
PYTHONPATH=/private/tmp/chematic-inchi-batch-python:/private/tmp/chematic-rdkit-2026091/lib/python3.13/site-packages \
python3.13 validation/results/positioned-inchi-keys-2026-10-19/oracle.py

MYPYPATH=/private/tmp/chematic-inchi-batch-python \
python3.13 -m mypy \
  --config-file validation/results/positioned-inchi-keys-2026-10-19/mypy.ini \
  crates/chematic-py/python/chematic/_inchi_batch.py \
  validation/results/positioned-inchi-keys-2026-10-19/typing_valid.py
```
