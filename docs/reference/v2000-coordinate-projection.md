# V2000 source coordinate projection

MOL/SDF atom records carry three coordinate fields. The topology-only conversion
bridge previously discarded them for V2000 input: converting a 3D MOL to CJSON
omitted `atoms.coords.3d`, and Python MOL→XYZ required a separate coordinate
argument even though the source already contained geometry.

The V2000→CJSON bridge retains the finite source-order x/y/z values in both
Python and WASM when every parsed graph property is representable by the local
CJSON profile. It applies to one MOL block or a one-record SDF document; the
existing SDF framing boundary refuses extra records. It does not generate or
recenter coordinates. The numerical projection is independent of the header's
`2D`/`3D` declaration and the tolerance used to select a perception conformer.
If 3D perception assigns atom chirality, the checked CJSON bridge refuses the
conversion instead of dropping that assignment; the low-level coordinate reader
still returns the original coordinates.

```python
import json
import chematic

mol, name, xyz = chematic.from_mol_block_with_coords_3d(mol_text)
cjson = json.loads(chematic.convert_format(mol_text, "mol", "cjson"))
assert cjson["atoms"]["coords"]["3d"] == [v for point in xyz for v in point]
```

## API and boundaries

Rust exports `parse_mol_with_coords_3d(&str)`, returning
`Result<(Molecule, MolMetadata, Vec<(f64, f64, f64)>), MolParseError>`.
Python `from_mol_block_with_coords_3d(str_or_utf8_bytes)` returns
`(Mol, name, [[x, y, z], ...])`. These low-level readers accept a V2000 block;
use the SDF document/iterator APIs for multi-record input. They do not audit or
retain later SD fields or opaque CTAB metadata. Existing graph/2D readers keep
their historical behavior.

The new reader:

- Requires finite numbers in all three 10-byte coordinate fields. Empty fields,
  invalid tokens, NaN and infinities fail before geometry perception; they are
  never replaced with zero.
- Keeps every parsed binary64 coordinate, including negative zero, tiny z,
  flat/all-zero blocks and translated planes, in atom-block order.
- Limits source text to 16 MiB and reuses native atom/bond/count validation.
  Errors identify the atom line and failed axis. WASM's bridge additionally
  retains its 1 MB input and 10000-atom limits.
- Delegates graph and stereo interpretation to the existing V2000 reader.
  Retaining source coordinates does not prove that a dimension declaration or
  stereochemical assignment is correct.

Python also supplies these coordinates to PDB, XYZ and PDBQT output. Each target
keeps its existing precision/metadata boundary; these routes do not promise an
exact molecular round trip. Explicit `coords=` retains precedence. For CJSON,
`coords=[]` explicitly omits geometry. Source fields must still be valid when
an override is supplied. MOL/SDF→MOL/SDF source retention and other topology-only
routes are unchanged.

Use `xyz` with `write_xyz_checked(..., notation="scientific", decimal_places=16)`
when binary64 coordinate round trips through XYZ are required. See
[XYZ precision](xyz-output-precision.md); not every external XYZ reader accepts
scientific notation. CJSON output retains coordinate values directly, alongside
the locally supported indexed elements, formal charges, isotopes and bond orders.
CJSON projection does not retain SD properties, source layout, all CTAB query
semantics or stereochemical declarations. Full chemical equivalence is a
separate contract.

## Bounded comparison

The authored packet `validation/v2000_coordinate_projection/cases.json` contains
10 profiles, each exercised as MOL/SDF with LF, CRLF and no final newline:
60 routes. Fifty-four graph/coordinate projections succeed; the six tetrahedral
routes return a typed atom-chirality refusal. It includes ordinary four-decimal
coordinates and deliberately extended precision/scientific spellings inside the
10-byte fields. The latter are accepted-input boundary tests, not a claim that
every V2000 reader supports those spellings.

A separate fixed-column decoder checks the source coordinates, and the shared
CJSON observer checks output coordinate bits and indexed graph fields without
using chematic readback as its only oracle. Python and real Node/WASM agree on
all 54 representable projections and all six typed refusals. The old local
package emits no coordinate array on all 60.

Open Babel 3.2.1 is compared on the **same source bytes and CJSON target**, using
both default output and `-xs` (single geometry). It retains coordinate bits on
54/60 projections. The six tiny-z/2D-header variants return no CJSON, with an
explicit no-3D error rather than silently fabricated geometry. Its isotope
array omits the source labels on the six ionic-isotope variants. It matches all
declared fields on 48/60 projections. Its inferred partial charges and other
metadata are outside this comparison. No parser-speed or overall-toolkit rank
is inferred from this small synthetic packet.
The isotope check uses chematic's existing local `atoms.isotopes` profile; it
does not assert that every ChemicalJSON implementation promises that field.

Exact commands, outputs, binary/source hashes and every comparison row are in
`validation/v2000-coordinate-projection-2026-10-24/`. The
[Open Babel CJSON writer](https://github.com/openbabel/openbabel/blob/openbabel-3-2-1/src/formats/json/chemicaljsonformat.cpp)
documents the `-xs` option (also shown by `obabel -L cjson`); its
[MOL profile](https://openbabel.org/docs/FileFormats/MDL_MOL_format.html)
documents V2000/V3000 support and coordinate-dependent stereo interpretation.
