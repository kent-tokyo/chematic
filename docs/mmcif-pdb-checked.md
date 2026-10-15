# Checked mmCIF to PDB coordinate projection

Open Babel [#3027](https://github.com/openbabel/openbabel/issues/3027)
reports a five-character component identifier being truncated during PDB
export. PDB format 3.3 provides only three columns for a residue name.
`A1AJ1` therefore cannot be represented faithfully in that field.

The new Rust `write_mmcif_pdb_checked` and Python `mmcif_to_pdb_checked`
reject names that do not fit. Metadata-loss consent never permits truncating
names, inventing missing identifiers, renumbering atoms/residues, or merging
identities. Use mmCIF when these identifiers must be retained in the output.

```python
import chematic

result = chematic.mmcif_to_pdb_checked(
    mmcif_text,
    namespace="auth",
    allow_metadata_loss=True,
)
pdb_coordinate_block = result["pdb"]
original_mmcif = result["source_text"]
source_atom_site = result["source"]
```

`namespace` must be `auth` or `label`; there is no fallback to the other
namespace. It selects the chain, residue number, component name, and atom
name together. Both namespaces remain available in `source`, using exactly
the schema of `parse_mmcif_residues`. `source_text` contains the complete
original Python input, including metadata beyond the reader's parsed view.
Retain that text when saving the coordinate block.

In Rust, pass a borrowed `MmcifResidueDocument` and `PdbNamespace::Author`
or `PdbNamespace::Label`. The writer does not mutate it. The original mmCIF
text is not part of that Rust document: callers must retain their input text
separately if they need unparsed categories or original numeric spelling.

## Representability checks

The writer follows the [wwPDB coordinate field definitions](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html).

- Model 1 only, at most 99,999 atoms, unique decimal serials from 1 to 99,999.
- Nonempty ASCII names: residue at most three characters, atom at most four.
- One alphanumeric chain character; one alphanumeric alternate-location or
  insertion-code character when present. Absent alternate/insertion codes
  produce blanks.
- Decimal residue numbers from -999 through 9999. Author numbers must have
  canonical integer spelling: `01`, `+1`, and `-0` are rejected.
- Distinct source residue groups must not collapse to the same output
  `(chain, residue number, insertion code)`. Distinct source atoms must not
  collapse to the same residue/name/alternate-location key.
- Finite coordinates, occupancies, and isotropic B values. Occupancy is in
  [0, 1]. Decimal rounding must still fit the fixed-width output field.
- Formal charge must fit the two-character PDB field, from -9 through +9.
  Element symbols are right-aligned and uppercase.

Output coordinate records have 80 columns. Coordinates use three decimal
places; occupancy and B values use two. The report gives the maximum absolute
rounding error for each of these three categories. Numeric precision remains
in the source sidecar and original text.

PDB cannot carry both identity namespaces or every mmCIF field, so the
function defaults to rejecting metadata loss. Set `allow_metadata_loss=True`
explicitly to obtain a representable coordinate block and a sidecar.
`losses` lists known omissions in the parsed view: `identity_namespaces`,
`entity_identifiers`, `cell_and_symmetry`, `unhandled_atom_site_columns`,
`numeric_rounding`, and `optional_formal_charge_presence` as applicable.
It is not an exhaustive category-by-category comparison of the original text.

## Scope and validation

This is a coordinate block consisting of ATOM/HETATM records and END.
It does not produce TER, CONECT, CRYST1, ANISOU, SEQRES or other deposition
records, infer bonds, reorder polymer chains, validate CCD chemistry, or
support multiple models/hybrid-36 identifiers. It is not a complete PDB
exporter. It inherits the existing mmCIF reader's tokenizer, default input
limits, and integer-only author sequence support.

The independent fixture
`validation/openbabel_issues/mmcif/long-component-name.cif` contains two
synthetic C/O coordinate records, not a chemical reconstruction of A1AJ1 or a
copy of an upstream structure. It is original project data under
MIT OR Apache-2.0. The validation manifest records provenance, minimization,
fixture/source/module hashes, the wwPDB 3.3 oracle, and Open Babel 3.2.1.

With a numeric label-sequence control, Open Babel 3.2.1 still truncates
`A1AJ1` to `A1A`; changing only the name to `UNL` preserves it. The checked
adapter rejects the long name in both namespaces regardless of loss consent.
A six-row NAG/NAG/BMA fixture passes all 720 input-row permutations with zero
fixed-column or source-sidecar mismatches. Open Babel also reads the native
PDB block back with all three residue identities and coordinates retained.
This checks coordinate interoperability, not inferred connectivity.

Reproduce with the built Python extension on `PYTHONPATH`:

```sh
python scripts/validate_mmcif_pdb_checked.py --obabel /opt/homebrew/bin/obabel
cargo test -p chematic-mol --lib --test mmcif_pdb_checked \
  --test mmcif_residue_contract --test mmcif_pqr_integration --offline --quiet
python -m pytest crates/chematic-py/tests/test_mmcif_pdb_checked.py \
  crates/chematic-py/tests/test_mmcif_residues.py \
  crates/chematic-py/tests/test_new_formats.py \
  crates/chematic-py/tests/test_input_errors.py -q
```

Local verification passed 784 Rust tests (12 new checked-writer tests) and
119 Python tests (23 new API tests). Clippy with warnings denied, formatting,
refactor budget, product naming, binding inventory, and strict MkDocs passed.

Local candidate only; this does not imply upstream resolution or a published
chematic package. The patch builds on the previous residue-reader candidate.
