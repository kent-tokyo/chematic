# DL-POLY CONFIG: retain simulation vectors and cell data

DL-POLY CONFIG is a new read/write format in chematic. REVCON uses the same
record layout. Rust, Python, WASM/Node and the CLI expose an editable record
that retains positions, optional velocities and forces, Cartesian cell vectors,
opaque names, optional IDs and the count declaration. It does not project the
record through a chemical graph.

This closes a gap found by comparing Open Babel's installed format inventory
with actual source APIs. Tinker XYZ is a separate pending change; Molden, MSI
BGF and DL-POLY HISTORY remain outside this implementation.

## Usage

```python
import chematic

source = "original\n2 0 1\nDPD-A 1\n0 0 0\n1 2 3\n4 5 6\n"
record = chematic.parse_dlpoly_config(source)  # also accepts UTF-8 bytes
record["atoms"][0]["velocity"] = [1.2345678901234567, 2, 3]
record["atoms"][0]["force"] = [4, -123.45678901234567, 6]
output = chematic.write_dlpoly_config(record)
assert chematic.parse_dlpoly_config(output) == record
```

```sh
chematic dlpoly-config -i CONFIG -o normalized.CONFIG
chematic dlpoly-config --inspect -i CONFIG > record.json
chematic dlpoly-config --from-json -i record.json -o edited.CONFIG
```

Rust exports `DlpolyConfig`, `DlpolyAtom`, `DlpolyConfigLimits`, `DlpolyError`,
`parse_dlpoly_config`, `parse_dlpoly_config_with_limits`,
`dlpoly_config_from_json` and `write_dlpoly_config` from `chematic_mol`.

WASM provides `dlpoly_config_record_json(text)` and
`dlpoly_config_from_record_json(json)`. Use `JSON.parse` on the reader result,
edit it and pass `JSON.stringify(record)` to the writer. All bindings use this
same record schema:

```json
{
  "title": "original",
  "levcfg": 2,
  "imcon": 0,
  "declared_atom_count": 1,
  "cell_vectors": null,
  "atoms": [{
    "name": "DPD-A", "index": 1, "header_extra": [],
    "position": [0, 0, 0], "velocity": [1, 2, 3], "force": [4, 5, 6]
  }]
}
```

## Units and external context

`position`, `velocity`, `force` and `cell_vectors` retain source numbers.
CONFIG does not encode whether those numbers use DL-POLY internal or DPD units.
The API therefore adds no unit label and performs no conversion. FIELD/CONTROL
context must travel separately with the caller's workflow.

For ordinary internal units, lengths are angstroms, velocities angstroms/ps
and forces angstrom daltons/ps². DPD uses another context. Names such as `OW`
are force-field labels, not element declarations. No masses, charges, bonds,
isotopes or atomic numbers are inferred from them or from suffix words.
For example, Open Babel adds a third atom-header word; `header_extra` retains
such words without assigning a physical meaning.

## Supported profile

The parser supports `levcfg` 0 (positions), 1 (positions and velocities) and
2 (positions, velocities and forces). `imcon` 0 has no cell. Modes 1/2 use
positive axis-aligned cubic/orthorhombic cells in this profile; mode 3 retains
general nonsingular vectors, including basis orientation and handedness.
No cell normalization, wrapping or atom reordering is performed.

Three-field control records include an atom count. The two-field Classic
variant keeps `declared_atom_count: null`. Explicit IDs must uniquely cover
1..N, but may appear in any order. Name-only Classic records keep `index: null`;
mixed indexed/name-only headers are refused. Atom-header suffix words require
an explicit index and are retained in source order.

Format bounds are 72 ASCII bytes per record, 40 bytes per non-title word and
8 bytes per atom name. An empty title is permitted; blank or comment body
records are refused. CRLF and Fortran `D`/`d` exponents are accepted. Whitespace
and number spelling are normalized on output. Titles and suffix token strings
are retained exactly.

The writer emits current fields, validates level/vector consistency, and uses
the shorter of ordinary and scientific round-trip f64 text. It refuses a
constructed vector whose exact output cannot fit the 72-byte format limit;
it never resolves that overflow by rounding. If atoms are added/removed, update
a present `declared_atom_count` as well; stale declarations are errors.

Rust/Python/CLI default limits are 64 MiB and 100,000 atoms. Rust offers custom
parser limits; output uses the defaults. WASM has a 1,000,000-byte text-input
limit, 10,000 atoms, and 16 MiB for record JSON input/output and text output.
Output above the WASM text-input cap needs the Rust/Python reader. Errors carry
a one-based line where available, or zero for record validation/global limits.

Boundary modes 4/5/6/7, extended parallel-I/O control headers, HISTORY, binary
restart files, FIELD and CONTROL parsing are unsupported. Unknown JSON fields,
non-finite vectors, duplicate/missing IDs, inconsistent levels, extra content
and truncated blocks are refused. A name-only legacy input without an atom
count has no record delimiter; every syntactically valid atom block belongs to
that one record. This is not the common graph conversion bridge.

## Local comparison and verification

The original MIT-licensed packet contains 17 valid cases and 20 refusals.
Rust, Python, WASM/Node and CLI test the same expectations and edited records.
The packet covers all supported levels/modes, reordered IDs, optional IDs/counts,
opaque labels/suffixes, precise vector values and Fortran exponents. Rust tests
also cover finite f64 extremes, line limits and inconsistent constructed data.
WASM tests a 10,000-atom record whose JSON exceeds the text-input size cap.

The local Open Babel **3.2.1** comparison uses `-iCONFIG -oCONFIG -l1` on exactly
the same valid inputs. An independent table observer scores raw fields, without
chemical perception. `-l1` selects the single-record workflow. Each process has
a 60-second limit; an initial 15-second bound was insufficient for the periodic
inputs. The checked-in receipts
include commands, input/source/binary hashes, stdout/stderr, terminal status
and per-case observations. Refusal cases are not fed to the comparator.

Both engines emit observable records for **17/17** inputs with the final bound.
chematic preserves all declared fields. Open Babel emits `levcfg: 0` and
`imcon: 0` in every case, omitting the input velocity/force blocks and cells.

| Data present in the input | chematic preserved | Open Babel preserved |
|---|---:|---:|
| Velocities | 14/14 | 0/14 |
| Forces | 13/13 | 0/13 |
| Cell vectors and boundary mode | 5/5 | 0/5 |
| Exact position values | 17/17 | 16/17 |
| Atom names | 17/17 | 16/17 |
| Explicit/absent IDs in source order | 17/17 | 15/17 |

Open Babel also replaces atom-header suffix words and omits a present count
declaration. The count remains inferable from its atom table; the report tracks
that declaration separately from vector/identity preservation. An absent input
vector is not counted as evidence of support for a present vector. These are
bounded packet observations, not general format/physics correctness or overall
superiority claims. Timeouts/execution failures have their own status and are
not interpreted as observed numeric loss.

Reproduce after building `chematic-cli`:

```sh
python3 scripts/validation/check_dlpoly_config_preservation.py \
  --cli target/debug/chematic --obabel /path/to/obabel --output /tmp/config-comparison.json
```

The packet and receipts are under
[`validation/dlpoly_config`](https://github.com/kent-tokyo/chematic/tree/main/validation/dlpoly_config).
They describe a local source candidate, not a published release or hosted CI.

The implementation follows the public
[DL_POLY CONFIG description](https://ccp5.gitlab.io/dl-poly/UserManual/data_files/input_files.html#the-config-file)
and [DL_POLY unit definitions](https://ccp5.gitlab.io/dl-poly/UserManual/introduction.html#units).
Open Babel lists the format on its
[CONFIG page](https://openbabel.org/docs/FileFormats/DL-POLY_CONFIG.html).
No Open Babel implementation source was copied.
