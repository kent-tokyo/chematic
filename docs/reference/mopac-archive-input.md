# MOPAC archive input

`parse_mopac_archive` reads one MOPAC `.arc` final-geometry restart deck in
Rust and Python. This addresses the archive-reading request in
[Open Babel #2253](https://github.com/openbabel/openbabel/issues/2253), verified
open on 2026-10-17. The `.arc` extension also names unrelated formats; this API
identifies a MOPAC archive by its final-geometry marker, not its extension.

```python
import chematic

archive = chematic.parse_mopac_archive(text, max_atoms=4096)
for atom in archive["atoms"]:
    if atom["references"] is None:
        xyz_angstrom = atom["values"]
    else:
        distance_angstrom, angle_degrees, dihedral_degrees = atom["values"]
        one_based_references = atom["references"]
```

The public Rust APIs are `parse_mopac_archive` and
`parse_mopac_archive_with_limits`, returning `MopacArchive` with typed
`MopacArchiveAtom` rows. They reuse `MopacParseLimits` and `MopacError`.
Python returns a dictionary with the same field names, using canonical element
symbols and optimization flag strings.

The reader accepts exactly one `FINAL GEOMETRY OBTAINED` line, optionally
followed by `CHARGE`, then the keyword/title/comment lines and contiguous atom
rows. A blank line ends the atom block. Each atom contains three finite values
and three optimization flags. An optional group of three integer reference
numbers identifies MOPAC internal-coordinate fields. Reference numbers remain
one-based, with zero retained as an unused reference. Nonzero references must
identify other atoms in the parsed block, without repetition within a row;
forward references are retained. No reference-order normalization is performed.
The early internal-coordinate rows contain unused angle/dihedral fields;
these values are preserved without interpreting them as physical angles.

Rows without references retain Cartesian coordinates in Å. The reader does
not convert internal coordinates to Cartesian coordinates. A `CHARGE` marker
requires exactly one additional finite charge value per row; without that
marker, extra fields are rejected. Canonical element symbols, atomic numbers,
parenthesized labels, Fortran D exponents, and LF/CRLF are supported.

`summary_text` contains the exact text before the marker. `restart_text` is
the complete input deck after the marker. `trailing_text` retains the blank
separator and any ancillary sections, including `SYMMETRY` instructions.
`source_text` retains the whole archive. These raw fields are not interpreted
or executed. A parsed archive does not establish convergence or an energy,
and does not infer topology, bond orders, formal charge or spin. Ancillary
sections are preserved without syntax or chemical validation.

Atomic-unit input (`A0`), keyword continuations, periodic translation vectors,
dummy/sparkle encodings, symbolic coordinates, unflagged rows, unannounced
extra columns and multiple final geometries are refused. This is a bounded
archive reader, not a general MOPAC output reader or restart executor.

Input-byte and line-length limits apply to the entire archive, including raw
sections; atom limits apply before retaining another row. Errors in atom rows
refer to original archive line numbers. No partial result is returned.

The format interpretation follows the official
[MOPAC input/archive guide](https://openmopac.net/guides/basics/),
[Cartesian coordinate definition](https://openmopac.net/Manual/cart_coord.html)
and [geometry constraint documentation](https://openmopac.net/Manual/Applying_geometry_constraints.html).
Validation evidence is recorded in
`validation/results/mopac-archive-2026-10-17/` with original synthetic fixtures
and a separately fetched upstream attachment diagnostic. No MOPAC calculation
was run as part of this validation.

Local validation passed 805 existing library tests, 7 new Rust tests and 277
Python tests (9 new). The 512-record preservation corpus and 16 independent
Open Babel Cartesian cross-reads passed; the separately fetched issue
attachment retained all four atoms and its ancillary symmetry section.
