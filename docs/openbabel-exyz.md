# ORCA-AICCM EXYZ geometry

This candidate adds an independently implemented EXYZ reader and writer to
Rust and Python. It follows the [Open Babel EXYZ format description](https://openbabel.org/docs/FileFormats/Extended_XYZ_cartesian_coordinates_format.html)
(retrieved 2026-10-15). No Open Babel implementation code is incorporated.

The previous candidate already supported ASE ExtXYZ in `xyz/extxyz.rs`.
ORCA-AICCM EXYZ uses a different grammar: `%PBC`, `%VIRTUAL`, and a separate
`Vector1`/`Vector2`/`Vector3`/`Offset` block. The existing ASE reader cannot
represent this complete document. This implementation adds a separate API.

## Python

```python
import chematic

text = chematic.write_exyz(
    ["Si", "Cl"],
    [[0.12345678901234567, -2, 3], [-1, 0, 2]],
    header="original synthetic example",
    virtual_atoms=[False, True],
    lattice_angstrom=[[4, 1, 0], [-1, 5, 2], [1, 0, -6]],
    offset_angstrom=[0.5, -1, 2],
)
geometry = chematic.parse_exyz(text)
assert geometry["virtual_atoms"] == [False, True]
assert geometry["offset_angstrom"] == [0.5, -1.0, 2.0]
```

`parse_exyz` returns `header`, `elements`, `positions_angstrom`,
`virtual_atoms`, `lattice_angstrom`, `offset_angstrom`, and `source_text`.
Nonperiodic documents have `None` for both cell fields. The original second
line is retained, including whitespace and directive order. Original source
text is a sidecar; writing uses current structured values.

The writer adds `%PBC` for an explicit lattice and `%VIRTUAL` when any atom
is marked virtual. A supplied `%VIRTUAL` remains present even when all atom
flags are false. A lattice without an explicit offset gets `[0, 0, 0]`.
An offset without a lattice is an error. Parallel arrays must match exactly.
Canonical element symbols are required; isotope aliases and dummy symbols
are refused. Errors are Python `ValueError` exceptions.

Use `from_extxyz` / `to_extxyz` for the existing Python ASE dialect.
The new APIs do not add EXYZ to the common format conversion bridge.

## Rust

The crate root exports `ExyzGeometry`, `ExyzAtom`, `ExyzCell`,
`ExyzParseLimits`, `ExyzError`, `parse_exyz`, `parse_exyz_with_limits`, and
`write_exyz`. Each atom stores an `Element`, a Cartesian position and its
virtual flag. The optional cell stores three Cartesian row vectors and an
offset. All numeric quantities use angstroms.

```rust
use chematic_mol::{ExyzError, parse_exyz, write_exyz};

fn edit_geometry() -> Result<String, ExyzError> {
    let text = "1\n%VIRTUAL\nC 1 2 3 VIRTUAL\n";
    let mut geometry = parse_exyz(text)?;
    geometry.atoms[0].position_angstrom[0] = 0.125;
    write_exyz(&geometry)
}
```

Public mutable Rust fields are validated before writing, including header
consistency. `ExyzError` has a detail and one-based source line when a field
is rejected during parsing. Whole-document consistency errors and errors in
constructed geometry use line zero.

## Supported boundary

- Exactly one document with a positive, exact atom count. Plain XYZ without
  EXYZ directives is accepted as a nonperiodic subset. Trailing nonblank
  data, including another frame, is refused.
- Optional atom `VIRTUAL` markers require `%VIRTUAL`. Unknown or duplicate
  percent directives and extra atom columns fail explicitly. ASE `Lattice=`
  or `Properties=` headers are refused with a dialect-specific error.
- `%PBC` requires a blank separator followed by exactly four ordered rows:
  `Vector1`, `Vector2`, `Vector3`, `Offset`. The offset is retained and is
  never added to atom positions. Rows are not rotated or wrapped.
- Numeric fields must be finite. Lattice row lengths must be finite and
  positive, and the absolute determinant of the normalized row matrix must
  exceed `1e-6`. Both right- and left-handed cells are accepted.
- Defaults: 64 MiB input, 65,536 bytes per line and 1,000,000 atoms. Python
  accepts smaller limits through keyword arguments. The atom cap remains
  a hard maximum. A tiny input with an enormous count fails before any
  count-dependent allocation. Canonical output must fit 64 MiB.
- Output uses scientific notation with 17 decimal places to retain finite
  double precision, including subnormals and negative zero. Source line
  endings may be LF or CRLF; canonical writing uses LF.
- This is explicit geometry interchange. Bond perception, isotope metadata,
  charge, spin, topology, implicit hydrogens, symmetry expansion, trajectories,
  and WASM/CLI bindings are outside this API.

## Validation

Local candidate validation on 2026-10-15 used Open Babel **3.2.1** and Python
3.13.6. `scripts/exyz_openbabel_gate.py` pins that oracle version and generates
192 original documents from 12 small geometry recipes, periodic/nonperiodic
settings, four virtual-flag styles, and LF/CRLF variants. Every document is
checked against its original arrays and an exact structured round trip.
The 96 periodic documents additionally pass through Open Babel CJSON and
EXYZ output with bonding disabled (`-ab`).

| Check | Result |
|---|---|
| Original structured round trips | 192 passed |
| Periodic Open Babel comparisons | 96 passed |
| CJSON coordinate maximum error | 0 Å |
| CJSON lattice maximum error | 8.882e-16 Å |
| Open Babel EXYZ output maximum error | 4.996e-6 Å |
| Unexpected mismatches | 0 |
| Rust library and new integration tests | 753 + 9 passed |
| Python EXYZ/TURBOMOLE/GZMAT tests | 124 passed |

Open Babel EXYZ output uses five decimal places. The comparison permits
`5.00001e-6` Å for its rounded output and `2e-12` Å for CJSON lattice values.
Those tolerances do not weaken chematic's exact internal round-trip check.

The pinned oracle has observed differences, recorded separately from passing
common geometry checks:

- Its EXYZ writer omits per-atom `VIRTUAL` markers: 132 marked atoms across
  the generated periodic cases lose their marker. chematic retains them.
- Its EXYZ reader rejects the tested plain nonperiodic EXYZ document, asking
  for unit-cell data even though the header has no `%PBC`.
- Converting that plain XYZ document to EXYZ adds `%PBC` and an identity cell.
  chematic's writer retains an absent cell.

These measurements establish the stated periodic geometry contract, with
explicit differences for virtual flags and nonperiodic behavior. They do not
establish general Open Babel feature parity.

All input geometries are original synthetic project data under
`MIT OR Apache-2.0`; no external molecule corpus is redistributed. The report
records per-case provenance, minimization, oracle settings, expected results
and hashes of normalized LF UTF-8 payloads.

Regenerate with a built candidate Python extension on `PYTHONPATH`:

```sh
python scripts/exyz_openbabel_gate.py --obabel /path/to/obabel
```

The standalone candidate and portable patch retain the earlier feature work.
No primary-checkout source edits or commit/push are part of this delivery.
