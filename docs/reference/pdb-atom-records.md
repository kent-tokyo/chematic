# Typed PDB atom attributes

`chematic.read_pdb_atom_records(text)` provides immutable, source-ordered
ATOM/HETATM records with coordinates, residue identity, altloc, occupancy and
B factor. It addresses the PDB attribute-access portion of open
[Open Babel #2708](https://github.com/openbabel/openbabel/issues/2708).
The existing general selection API exposes atom identity but not occupancy or
B factor; the separate amino-acid protein projection already exposes some
attributes. This addition covers general coordinate records, including ligands
and water, without turning them into a protein hierarchy or a bonded graph.

```python
from chematic import read_pdb_atom_records

text = (
    "ATOM      1  CA AVAL A 101       1.125  -2.500   3.750  0.50 23.45           C  \n"
    "HETATM    2  O1 BLIG A 101       1.125  -2.500   3.750  0.25 45.67           O  \n"
    "END\n"
)
result = read_pdb_atom_records(text)
assert result.atoms[0].alt_loc == "A"
assert result.atoms[0].occupancy == 0.5
assert result.atoms[0].b_factor == 23.45
assert result.atoms[1].record_type == "HETATM"
assert result.atoms[1].source_index == 1
assert result.atoms[1].line_number == 2
assert result.source_text == text
```

The fixed columns follow the
[wwPDB 3.3 coordinate specification](https://www.wwpdb.org/documentation/file-format-content/format33/sect9.html).
Each `PdbAtomRecord` contains zero-based source atom index, one-based physical
line number, decimal serial, record type, raw four-column atom name, altloc,
trimmed residue name, chain character, residue number, insertion character,
Cartesian angstrom coordinates and optional occupancy/B factor. Character
identifiers retain a blank as `" "`, rather than guessing a value. Atom order
and record type are preserved, including both alternate-location rows.

Blank or absent occupancy/B-factor fields become `None`; explicit zero stays
`0.0`. Malformed, truncated nonblank or nonfinite optional fields are refused
instead of becoming defaults. This is numeric field extraction, not a physical
occupancy or displacement-model assessment. The exact original record is
available in `record_text` without its line ending. `element_field` and
`charge_field` contain original column slices, including spaces/absence; they
are not inferred chemical assignments.

`PdbAtomRead` holds the original text, optional explicit MODEL number, atom
tuple and `unhandled_records`: `(source_line, record_name)` pairs reported by
the native coordinate validator. This is an advisory list of metadata that the
validator does not interpret, not a complete list of all source-only semantics.
ANISOU, CONECT, headers, remarks and all other records remain in `source_text`.
No normalized/filtered PDB text is returned and no files are written. LF/CRLF,
blank fields and source numeric spellings stay intact in the original text.
All returned records, coordinates and collections are immutable.

## Scope and limits

The existing checked native coordinate-selection reader validates framing,
identity, finite numbers, serial uniqueness and references. An empty internal
selection avoids exporting rewritten coordinates; the Python view then extracts
all validated original atom rows. This retains that reader's bounded,
single-model profile and strict ANISOU/TER/CONECT checks. Multiple MODEL blocks,
hybrid-36/extended serials, non-ASCII text, malformed fields and unsupported
framing raise errors without partial results. Source metadata is retained, so
reading requires no metadata-loss consent.

Hard/default caps are 8 MiB input, 50,000 ATOM/HETATM records, 250,000 physical
lines and 1,024 bytes per native line. Callers may lower them with
`max_input_bytes`, `max_atoms`, `max_records` and `max_line_len`. Limits accept
integer-index values and reject Python bool. Input must be text contents,
not a path or bytes; empty input returns an empty atom tuple.

No bonds, chemical elements/charges, residue templates, conformers or altloc
selection are inferred. This does not add a full MDAnalysis converter,
trajectory ownership API, PDB writer or complete deposition validation. The
existing graph/protein/selection APIs retain their contracts and defaults.
There is no Gemmi, Open Babel or MDAnalysis runtime dependency.

## Validation

New tests include all 720 permutations of six atom rows, occupancy/B-factor
retention and blanks versus zero, altloc and insertion codes, ATOM/HETATM,
model numbers, CRLF/LF identity, immutable results, source-only metadata,
ANISOU/TER/CONECT, invalid optional numbers and input/resource failures.

An independent Gemmi 0.7.3 PDB parser checks 96 seeded cases / 2,592 atoms,
including shuffled rows, alternate locations, blank chains, insertion codes,
model numbers and mixed ATOM/HETATM. Coordinates agree exactly. Maximum
occupancy/B-factor differences are `2.87e-8` / `3.67e-6`, consistent with
Gemmi's float32 storage; the new API retains source-derived Python floats.

For the original two-atom control, local Open Babel 3.2.1 preserves occupancy
but rewrites altloc A/B to blanks and B factors 23.45/45.67 to zeros. This is
one observed read/write path, not a claim about every upstream version or API.
The chematic read retains both attributes and complete original text.

Evidence, reproducible oracle and provenance are under
`validation/results/pdb-atom-records-2026-10-18/`. Tests reuse a verified existing
native extension; relevant selection/strict-PDB source and binding files match
its source revision. No fresh current-main native build, published-wheel,
large real-structure corpus or MDAnalysis integration validation is claimed.

## Python reference

::: chematic._pdb_atom_records
