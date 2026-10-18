# Checked Maestro output

`write_mae_checked` writes a single ASCII Maestro (`.mae`) structure from a
concrete native graph and caller-supplied 3D coordinates. `v2000_to_mae_checked`
provides a bounded single-record MOL/SD adapter. This implements the missing
chematic capability associated with [Open Babel #2774](https://github.com/openbabel/openbabel/issues/2774).
The issue concerns format availability on Windows; installed Open Babel 3.2.1
on this host also has no registered `mae` format. This does not establish that
all Open Babel builds lack it.

```python
import chematic as ch

water = ch.from_smiles("[H]O[H]")
result = ch.write_mae_checked(
    water,
    [[0.96, 0.0, 0.0], [0.0, 0.0, 0.0], [-0.24, 0.93, 0.0]],
    title="water",
    allow_metadata_loss=True,
)
assert result.atom_count == 3 and result.bond_count == 2
assert "m_atom[3]" in result.text
assert "m_bond[2]" in result.text
```

Coordinates are in angstroms, in original atom order, including every explicit
H. No coordinates, atoms, bonds or force-field types are generated. The API
returns text in memory; `path.write_text(result.text, encoding="ascii")` can
save it. It requires explicit metadata-loss consent because source records,
SD fields, residue identity, atom/bond labels and force-field metadata are
outside this graph/geometry export. The input graph/properties stay unchanged.

## Supported profile

The writer reuses native `write_cdjson_graph` for checked indexed atom/bond
projection. It additionally requires all H to be explicit, no isotopes, and a
small ordinary closed-shell valence profile: H, C, N, O, F, P, S, Cl, Br, I and
selected Li/Na/K/Be/Mg/Ca ions. Neutral C has valence 4; N uses 3/4/2 for charge
0/+1/-1; O uses 2/3/1; neutral P allows 3/5 and S allows 2/4/6. Neutral halogens
have valence 1; halide ions and the listed metal ions have valence 0. H+ and H-
are supported isolated ions. Other combinations refuse.

Bond orders are plain single/double/triple. Aromatic/query/dummy/mapped atoms,
nonstandard bonds, declared stereo, enhanced groups, isotopes and insufficient
explicit valence refuse. No kekulization or hydrogen completion is performed.
The reader may perceive stereochemistry from coordinates; this API does not
promise stereo-label parity or physical geometry validity.

`title=None` uses `mol.GetName()`; an explicit title overrides it. Titles may be
empty or contain at most 256 printable ASCII characters. Quotes/backslashes are
escaped; controls/newlines and non-ASCII titles refuse. The frozen `MaeExport`
contains `text`, `title`, `atom_count`, `bond_count`, optional `source_text`, and
an immutable `omitted_features` tuple.

Positions serialize with 17 significant digits. Every coordinate must be finite
and within +/-1e8 angstrom; nonzero magnitudes below 1e-100 refuse. Input is an
indexed sequence of exact three-value rows. Booleans are excluded from numeric
values. Limits: 1..4096 atoms, <=8192 bonds, <=2 MiB ASCII output. Duplicate
edges, incorrect dimensions and exhausted limits return no partial result.

## V2000 MOL/SD adapter

```python
from pathlib import Path

source = Path("single.sdf").read_text(encoding="ascii")
report = ch.v2000_to_mae_checked(source, allow_metadata_loss=True)
assert report.source_text == source
Path("single.mae").write_text(report.text, encoding="ascii")
```

The adapter keeps the original source text and reads all three fixed-column
coordinates, including Z. It accepts one V2000 block with <=999 atoms/bonds,
optional M CHG records, and an optional opaque SD data tail and final `$$$$`.
It refuses multiple records, V3000, raw stereo/query/map flags, isotope/radical
flags, unsupported ancillary CTAB records, and valence overrides other than
0/15 (default/no implicit H). Native projection counts must match source counts.
Metadata is explicitly omitted, not silently converted into guessed chemistry.

Input is <=16 MiB ASCII, <=20000 lines and <=4096 bytes per line. SD values
remain opaque metadata; this adapter is not a complete SDF schema validator.
It never skips a malformed molecular record. A batch reader, Maestro reader,
`.maegz`, Rust/WASM/CLI API, GUI testing and package publication are outside this
addition. Wrong types raise `TypeError`; unsupported/invalid values raise
`ValueError`.

The independently implemented encoder follows the public block/property names
in [Schrodinger maeparser](https://github.com/schrodinger/maeparser) and its
[constants](https://github.com/schrodinger/maeparser/blob/master/MaeConstants.hpp).
Independent RDKit 2026.09.1 Maestro reading checks atom order/element/charge/H
counts, bonds/orders, title and coordinates for 1404 V2000 and 1404 direct-graph
cases, plus five coordinate-edge and two hypervalent graph controls. All 18959 atom and 16216 bond checks
pass, with zero observed coordinate difference. These are 39 authored
nonaromatic ordinary-valence molecules and selected ions, not all chemical space.
Reproduction: `validation/results/maestro-output-2026-10-18/`.
