# CCD molecular templates

The new reader addresses [Open Babel #1516](https://github.com/openbabel/openbabel/issues/1516):
Chemical Component Dictionary (CCD) files describe a molecular template using
`_chem_comp_atom` and `_chem_comp_bond`, rather than the `_atom_site` loop used
by the existing macromolecular mmCIF adapter. This is a separate API.

```python
import json
from pathlib import Path
import chematic

text = Path("GLY.cif").read_text(encoding="utf-8")
record = json.loads(chematic.parse_ccd_json(text))
mol = chematic.from_ccd(text)
assert record["component_id"] == "GLY"
assert mol.formula == "C2H5NO2"
coords = [atom["ideal_coordinates"] for atom in record["atoms"]]
```

## Record and graph contracts

`parse_ccd_component` / `parse_ccd_component_with_limits` in Rust return a
`CcdComponent`. Python and WASM expose the same record through `parse_ccd_json`.
They preserve source atom and bond row order, component/atom IDs, element,
formal charge, aromatic flag, stereo annotation and both optional Cartesian
coordinate families (Å). Coordinates are not inferred or silently substituted
between model and ideal families. All-absent/placeholder triples become `None`
(JSON `null`); a partial finite triple or a nonfinite value rejects.

Unmodeled atom/bond fields remain on `extra_fields`, with lower-case full tag
names. Tags in other categories are listed on `unhandled_tags`; their values
are outside this adapter. This is not a full mmCIF document round trip.

`CcdComponent.to_molecule(allow_stereo_loss)`, Python `from_ccd(text, *,
allow_stereo_loss=False)` and WASM `mol_from_ccd(text, allow_stereo_loss)` project
the explicit template graph. Every source H row remains an atom; no additional
implicit hydrogens, bonds, coordinates or polymer-linking edits are inferred.
Leaving-atom flags remain metadata; this API does not delete leaving atoms to
construct a protein residue. IDs and coordinates remain in the companion CCD
record, rather than being attached to the returned molecule.

R/S and E/Z annotations remain on the record but are **not assigned to graph
parity**. The graph APIs refuse known stereo annotations by default. Explicit
loss consent permits a graph with unassigned stereo; it does not assert stereo
preservation. Reading a chiral component's JSON record requires no such consent.

Supported bond orders are SING, DOUB, TRIP, QUAD and AROM. Aromatic flags or an
AROM order project to aromatic graph bonds, while the original supported order
remains in the CCD bond record. Element-symbol case is normalized (for example CCD `NA` becomes `Na`);
atom names are never used to guess elements. D/T spellings project to isotope 2/3 H.
Absent/placeholder charges remain `None` (JSON `null`); graph projection
refuses unknown charges rather than imputing neutrality. An invalid integer charge
rejects during reading. Absent/placeholder aromatic flags default to false. This is structural
validation, not valence sanitization or stereochemical validation.

## Bounds and refusals

The default native/Python limits are 8 MiB input, 500,000 tokens, 10,000 atoms
and 40,000 bonds. Expanded atom/bond loop key/value bytes are also capped
at the smaller of four times the configured input budget and 32 MiB, before
rows are materialized; long repeated column names cannot amplify that storage
without bound. Rust exposes all limits; Python `parse_ccd_json` accepts
`max_input_bytes`, `max_atoms` and `max_bonds` keyword overrides. WASM applies
its shared 1,000,000-byte input and 16 MiB serialized-output caps as well.

The reader accepts one data block and requires `_chem_comp.id` with matching
per-row component IDs. A second data block rejects, preventing a bulk CCD
file from being silently treated as its first component. Bulk dictionary
streaming, CCD writing, CIF2 collections and STAR save/global/stop controls
are outside this implementation.

Duplicate atom IDs, missing bond endpoints, duplicate/reversed/self bonds,
unknown elements, unsupported bond orders, malformed stereo/flag values,
duplicate columns and incomplete loop rows reject without a partial result.
Singleton category fields are supported for one-atom components. Quoted values,
comments, semicolon text fields and CRLF boundaries are handled independently
of the older mmCIF tokenizer. Python errors carry `format="ccd"` and typed
`category`/`code` attributes, including `ccd_stereo_loss` and `ccd_unknown_charge`.

## Reproduction and evidence

The implementation was independently authored from the wwPDB
[atom](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp_atom.html)
and [bond](https://mmcif.wwpdb.org/dictionaries/mmcif_pdbx_v50.dic/Categories/chem_comp_bond.html)
category definitions. The committed water template and malformed variants are
original synthetic data. No Open Babel or Gemmi implementation code was copied.
Gemmi is an isolated validation dependency, not a product dependency.

```sh
cargo test -p chematic-mol --lib --test ccd_contract --test hin_contract --test public_api_compat_v1_1_3 --offline --quiet
python3 -m pytest -q crates/chematic-py/tests/test_ccd.py crates/chematic-py/tests/test_hin.py crates/chematic-py/tests/test_new_formats.py
python3 scripts/validation/validate_ccd.py --download
node crates/chematic-wasm/tests/ccd.test.mjs
```

The oracle runner uses a source-built extension, pinned Gemmi 0.7.5 and eight
official RCSB CCD downloads: GLY, HOH, ALA, NA, CL, BEN, ATP and HEM. It checks
atom IDs/elements/charges, both coordinate families and bond endpoints/orders
against Gemmi, tests stereo-loss consent, and records both Open Babel CIF and
mmCIF read attempts. Download URLs and exact input/source/runtime hashes are
retained in `validation/results/openbabel-ccd-2026-10-16.json`.

Node uses `pkg-node` by default or `CHEMATIC_WASM_MODULE`. The original fixture
is compared to the Python record in CI. Setting `CHEMATIC_CCD_SOURCES` to the
oracle download directory also verifies all eight official records against
the committed Python/Gemmi evidence in a real WASM runtime. Source candidates,
package publication and general Open Babel parity are separate claims.

Local validation on main `2a0c8b5d` plus this change passed **793 Rust tests**
(12 new CCD tests, including the integrated HIN contracts), **79 Python tests** (13 new CCD tests, including the integrated HIN tests), and the real Node/WASM
contracts. Gemmi 0.7.5 matched **4,493 field checks** over eight official
components (167 atoms, 171 bonds), with zero failed comparisons. All eight
records also matched exactly across Python and WASM. HEM's unknown Fe charge
remains null in its record and refuses graph projection; known stereo labels
in ALA, BEN and ATP require explicit graph-loss consent. The local Open Babel
3.2.1 CIF and mmCIF attempts each returned zero molecules for all eight files.
These are bounded source/runtime observations, not a whole-dictionary claim.
