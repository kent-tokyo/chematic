# Distinct PDB hydrogen names

`name_pdb_hydrogens_checked` assigns H1, H2, ... names to existing explicit
hydrogen rows without regenerating the molecule. It addresses the naming request
in [Open Babel #2945](https://github.com/openbabel/openbabel/issues/2945), verified
open on 2026-10-18. The installed Open Babel 3.2.1 also writes four identical
`H` atom names when converting methanol SMILES to PDB with `-h`.

```python
import chematic as ch

source = (
    "HETATM    1  O   HOH A   1       0.000   0.000   0.000  1.00  0.00           O  \n"
    "ATOM      2  H   HOH A   1       0.960   0.000   0.000  1.00  0.00           H  \n"
    "ATOM      3  H   HOH A   1      -0.240   0.930   0.000  1.00  0.00           H  \n"
    "END\n"
)
result = ch.name_pdb_hydrogens_checked(source)
assert result.renamed_count == 2
assert [item.new_name.strip() for item in result.assignments] == ["H1", "H2"]
assert result.pdb.splitlines()[0] == source.splitlines()[0]
```

The API returns text in memory. `output_path.write_bytes(result.pdb.encode("ascii"))`
can save it without changing its line endings. It does not add hydrogens, generate
their positions, infer elements/parents, apply CCD residue-template names, or
change ATOM versus HETATM classification. The issue's record-classification
question is outside this addition; existing record kinds are preserved exactly.

## Assignment and preservation

`hydrogen_indices=None` selects every row whose explicit element column is `H`.
Pass zero-based ATOM/HETATM row indices to select only newly added H atoms, e.g.
`name_pdb_hydrogens_checked(source, [1, 2])`. These are source row indices,
independent of PDB serial numbers. Input selection order does not reorder atoms
or name assignment. An empty sequence returns the original text after checking.
Element columns must be nonempty on all atom rows. Names are never used to infer
an element; D/T element records stay unchanged and cannot be selected as H.

Names restart within each residue, identified by TER segment, chain, numeric
residue ID, insertion code and residue name. Interleaved rows with the same
identity share a name space. All names of unselected atoms in that residue are
reserved first, including names on heavy atoms, isotopes and alternate locations.
Selected H atoms receive the next available H1..H999 in source order. Names with
two or three characters start in column 14; four-character names occupy columns
13..16. Exhausting this finite name space is an error. Existing selected names
may change even when they were already unique; use explicit indices to retain them.

Only columns 13..16 of selected ATOM/HETATM rows and their associated ANISOU
rows can change. All other source bytes remain exact, including serials,
coordinates, occupancy/temperature, isotope/element columns, charge text,
residue identifiers, MODEL/TER/CONECT, blank lines, LF/CRLF/mixed line endings
and the presence or absence of a final newline. No atom is inserted, removed
or reordered. Repeating the same selection is idempotent.

`PdbHydrogenNaming` is frozen and contains `source_text`, `pdb`, `atom_count`,
`assignments`, `renamed_count` and `updated_anisou_lines` (one-based source line
numbers). Each frozen `PdbHydrogenName` records `source_index`, `serial`,
`old_name` and `new_name`; names are raw four-column fields. Assignments include
selected rows whose name was already correct. `renamed_count` counts changed
atom rows; ANISOU updates are reported separately.

## Validation and limits

The function reuses [checked PDB selection](pdb-atom-selection.md) to validate
the complete source and output. Supported input is one implicit coordinate model
or a complete single MODEL/ENDMDL pair, with valid ATOM/HETATM, ANISOU, TER,
CONECT, END and blank records. At least one atom is required. Unsupported
metadata, including LINK/SSBOND/SIGATM and HEADER/TITLE, is refused because it
can carry references or require a separate preservation contract. If exclusion
is intended, call `select_pdb_atoms(..., allow_metadata_loss=True)` separately
and keep its loss report before naming the resulting coordinate block.

Selected H atoms with alternate locations refuse; the API does not guess how
conformers should share labels. Unselected alternative rows remain unchanged.
Malformed fields, bad references, duplicate serials, invalid ANISOU associations,
multiple models and unsupported hybrid-36 serials are errors. No partial output
is returned. Input types use `TypeError`; invalid values/profile use `ValueError`.

Limits: 16 MiB ASCII input/output, 99,999 atoms, 1,000,000 physical records,
1,024 bytes per line and H1..H999 per residue. Integer indices accept the index
protocol and exclude bool; the selection must be an indexed sequence. The API
does not claim complete PDB-entry validation or chemistry/geometry correctness.
Rust/WASM/CLI naming entry points and package publication are outside this change.

The byte-preservation and naming recipe is independently implemented from the
[wwPDB 3.3 coordinate layout](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html).
Reproduction and local validation are in
`validation/results/pdb-hydrogen-names-2026-10-18/`.
