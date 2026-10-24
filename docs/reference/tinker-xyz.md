# Tinker XYZ: typed coordinates and explicit connectivity

Tinker XYZ (`txyz` in Open Babel) now has a dedicated editable read/write
record in Rust, Python, WASM/Node and the CLI. This closes a format gap found
by checking source APIs against Open Babel's installed format inventory.
Ordinary XYZ and ASE extended XYZ are separate formats; sharing `.xyz` does
not make their readers interchangeable.

## Usage

```python
import chematic

text = "3 water\n1 O 0 0 0 41 2 3\n2 H .9572 0 0 42 1\n3 H -.2399872 .927297 0 42 1\n"
record = chematic.parse_tinker_xyz(text)  # UTF-8 bytes also accepted
record["atoms"][0]["name"] = "OW"
record["atoms"][0]["atom_type"] = 90001
record["atoms"][0]["position_angstrom"][0] = 0.12345678901234567
output = chematic.write_tinker_xyz(record)
assert chematic.parse_tinker_xyz(output) == record
```

Rust exports `parse_tinker_xyz`, `parse_tinker_xyz_with_limits`,
`write_tinker_xyz`, `tinker_xyz_from_json`, `TinkerXyz`, `TinkerAtom`,
`TinkerCell`, `TinkerParseLimits` and `TinkerError` from `chematic_mol`.
The public structs are mutable and implement serde. The writer validates
current fields every time; changing a record changes its output.

WASM exports `tinker_xyz_record_json(text)` and
`tinker_xyz_from_record_json(json)`. Parse the first result with `JSON.parse`,
edit it, and pass `JSON.stringify(record)` to the second function.

```sh
chematic tinker-xyz -i water.xyz -o normalized.xyz
chematic tinker-xyz -i water.xyz --inspect > record.json
chematic tinker-xyz --from-json -i record.json -o edited.xyz
```

The JSON shape is identical across bindings:

```json
{
  "title": "one site",
  "cell": {"lengths_angstrom": [10, 11, 12], "angles_degrees": [70, 80, 90]},
  "atoms": [
    {"index": 1, "name": "OW", "position_angstrom": [0, 0, 0], "atom_type": 41, "neighbors": []}
  ]
}
```

## Preservation contract and supported profile

- One Cartesian frame, with an optional six-parameter cell line between the
  count/title and atom table. Cell lengths and positions are in angstroms;
  cell angles are in degrees. A missing cell stays `null`.
- IDs must be sequential, starting at 1. Names contain 1–3 printable ASCII
  characters and start with a letter. Names remain opaque, including names
  such as `OW` and `HW` that are not element symbols.
- Nonnegative force-field type numbers fit `u32`; zero denotes an unassigned
  type. Types are retained as identifiers. This does not validate their
  existence in an external `.prm`/`.key` force field.
- Every connection must be reciprocal and reference an existing other atom.
  Self connections and duplicate neighbors are refused. Neighbor-list order
  is preserved, but has no chemical significance in the comparison.
- Finite coordinates and cell parameters use round-trip f64 output. Invalid,
  singular or near-singular cell metrics (normalized determinant <= 1e-12)
  are refused. Atom coordinates are never wrapped or transformed by the cell.
- Title content after the count separator is retained, including Unicode.
  Leading separator whitespace is normalized. Edited titles with leading
  whitespace or control characters are refused. Other layout is canonicalized.
- Unknown JSON keys at every record level are refused. Non-finite numbers,
  invalid topology, truncated atom tables and extra nonblank content are errors.

Default Rust/Python/CLI limits are 64 MiB input/output, 65,536 bytes per line,
100,000 atoms and 64 neighbors per atom. Rust provides custom parser limits;
the writer uses the defaults. WASM limits text input to 1,000,000 bytes,
atoms to 10,000, and record JSON input/output and text output to 16 MiB. A
writer result exceeding the text-input cap needs the Rust/Python reader.
Errors include a one-based source line
where available; line zero indicates structured validation or a global limit.

ARC archives containing concatenated frames are refused. Tinker `.key`/`.prm`,
fractional `.frac`, internal-coordinate `.int` and velocity files are outside
this profile. Bond orders, isotopes, charges and element assignments are not
declared by this record; converting through a generic `Mol` or SMILES would
require a separate checked projection. This API is not part of the common
graph conversion bridge.

## Bounded Open Babel comparison

The original MIT-licensed packet has 12 valid records and 16 refusal cases.
Rust, Python, Node/WASM and CLI exercise the same packet, including edited
write/read cycles. Rust also checks resource limits and extreme finite f64
values. The comparator uses an independent table observer, not chemical
perception. Refusal cases are not fed to Open Babel or scored as format wins.

On the local 2026-10-24 run with Open Babel **3.2.1**, `-itxyz -otxyz -as -xc`
explicitly selected single-bond reading and custom atom-class output. This
preserves more source types than its default MM2 writer. Of 12 valid inputs,
10 produced observable records; both cell-bearing inputs produced no molecule
(`0 molecules converted`).

| Declared field | chematic, all 12 | Open Babel, 10 emitted records |
|---|---:|---:|
| Atom IDs | 12/12 | 10/10 |
| Atom names | 12/12 | 9/10 |
| Force-field type numbers | 12/12 | 9/10 |
| Exact f64 coordinate values | 12/12 | 1/10 |
| Explicit connectivity, ignoring neighbor order | 12/12 | 10/10 |
| Title | 12/12 | 0/10 |

Both cell-bearing records retain their six parameters in chematic (2/2).
Open Babel adds a `Custom parameters` suffix to titles, canonicalizes the
opaque names in one case, rewrites type `4294967295` to `0`, and rounds
coordinates to six decimal places. No-output cases are kept separate from
observed field changes. These results cover this packet and configuration;
they do not establish overall superiority or physical simulation accuracy.

Reproduce after building `chematic-cli`:

```sh
python3 scripts/validation/check_tinker_xyz_preservation.py \
  --cli target/debug/chematic --obabel /path/to/obabel --output /tmp/tinker-comparison.json
```

The checked-in packet and receipts are under
[`validation/tinker_xyz`](https://github.com/kent-tokyo/chematic/tree/main/validation/tinker_xyz).
The receipts include exact commands, input/source/binary hashes, stdout,
stderr and per-case field observations. They record a local source candidate,
not a released package or hosted-CI result.

## Sources and remaining format gaps

The implementation follows the public
[Tinker file-type description](https://dasher.wustl.edu/tinker/distribution/doc/sphinx/tinker/_build/html/text/file-types.html)
and [Tinker-GPU file-type description](https://tinker9-manual.readthedocs.io/en/latest/m/tutorial/filetypes.html).
Open Babel documents its reader/writer options on the
[Tinker XYZ page](https://openbabel.org/docs/FileFormats/Tinker_XYZ_format.html).
No Open Babel implementation code was copied.

The inventory also identified Molden (`molden`), MSI BGF (`bgf`) and DL-POLY
CONFIG/HISTORY as remaining gaps in the inspected Rust/Python/WASM source.
They are listed in the
[Open Babel format inventory](https://openbabel.org/docs/FileFormats/Overview.html).
POSCAR/CONTCAR already exists in `chematic-crystal`, and Python already has
an XSF writer: absence from the short matrix alone is not proof of a gap.
This change adds one format with explicit preservation boundaries rather than
treating aliases or read-only/output-only entries as separate complete formats.
