# Molden geometry and vibration records

Molden (`mold`, `molden`, `molf` in Open Babel) now has dedicated editable
geometry/vibration read/write APIs in Rust, Python, WASM/Node and the CLI.
This is a bounded profile: typical electronic-structure files containing GTO,
STO or MO sections are refused. They are never read as geometry while silently
discarding orbital information.

## Usage

```python
import chematic

text = """[Molden Format]
[FREQ]
-123.45678901234567
[FR-COORD]
H 0 0 0
H 0 0 1
[FR-NORM-COORD]
vibration 1
0 0 0.1
0 0 -0.1
[INT]
12.5
"""
record = chematic.parse_molden(text)  # UTF-8 bytes also accepted
record["intensities"][0]["raman"] = 0.0  # distinct from absent Raman (None)
record["normal_modes"][0]["displacements_bohr"][0][2] = 0.12345678901234567
written = chematic.write_molden(record)
assert chematic.parse_molden(written) == record
```

```sh
chematic molden -i vibrations.molden -o normalized.molden
chematic molden --inspect -i vibrations.molden > record.json
chematic molden --from-json -i record.json -o edited.molden
```

Rust exports `parse_molden`, `parse_molden_with_limits`, `write_molden`,
`molden_record_from_json`, `MoldenRecord`, `MoldenAtoms`, `MoldenAtom`,
`MoldenVibrationAtom`, `MoldenNormalMode`, `MoldenIntensity`,
`MoldenLengthUnit`, `MoldenLimits` and `MoldenError` from `chematic_mol`.
The public structs are editable and implement serde.

WASM exports `molden_record_json(text)` and `molden_from_record_json(json)`.
Use `JSON.parse` on the former result, edit it and send `JSON.stringify(record)`
to the latter. All bindings share this shape:

```json
{
  "atoms": null,
  "frequencies": [-123.45678901234567],
  "fr_coords": [
    {"symbol": "H", "position_bohr": [0, 0, 0]},
    {"symbol": "H", "position_bohr": [0, 0, 1]}
  ],
  "normal_modes": [{"index": 1, "displacements_bohr": [[0, 0, 0.1], [0, 0, -0.1]]}],
  "intensities": [{"ir": 12.5, "raman": null}]
}
```

An optional `atoms` object contains `unit: "angstrom"` or `"bohr"`, plus
`atoms: [{symbol, index, atomic_number, position}, ...]`.

## Preservation contract

The supported sections are `[Atoms] Angs|AU`, `[FREQ]`, `[FR-COORD]`,
`[FR-NORM-COORD]` and `[INT]`, following `[Molden Format]`. Atom coordinates
keep their explicit source unit; vibration coordinates and displacements keep
their specified Bohr unit. Values are never converted, normalized or mass-scaled.
Frequency and intensity values retain the source scale without an added unit
label; negative/zero frequencies stay negative/zero.

Mode IDs address the one-based frequency list independently of mode-block order.
Present modes must cover every frequency exactly once and contain one triplet
per vibration atom. Present intensity rows must match the frequency count.
Missing Raman columns stay `null`, while an explicit zero stays zero. Complete
vibration data needs its frequency/geometry/vector relationship; the API also
retains partial inspection records without mode vectors, or FR coordinates
alone, without claiming they form a complete vibration calculation.

`[Atoms]` requires canonical element symbols agreeing with atomic numbers 1–118.
IDs must uniquely cover 1..N. Both atom channels, when present, must agree on
the source-order element sequence. Their coordinates and units remain separate;
the API does not collapse or infer equality between coordinate channels.
Source atom and mode order is preserved. Layout, numeric spelling and section
order are canonicalized. Output uses round-trip f64 precision.

The writer always serializes current structured fields. Unknown JSON keys,
duplicate/empty sections, omitted/unknown atom units, non-finite numbers,
inconsistent IDs/dimensions and truncated vectors are refused. Electronic
sections (`GTO`, `STO`, `MO`, angular-basis markers), convergence sections,
`GEOMETRIES` trajectories and unknown sections are refused even when they follow
otherwise valid geometry. This is not part of the common molecule/SMILES bridge;
no bonds, charge, spin, isotopes or orbital data are inferred.

Default Rust/Python/CLI limits are 64 MiB input/output, 65,536 bytes per line,
100,000 atoms, 100,000 frequencies/modes and 1,000,000 total displacement
triplets across all modes. Rust offers custom parser limits. WASM uses 1,000,000
bytes of text input, 10,000 atoms/frequencies/modes, 100,000 displacement triplets
and 16 MiB for JSON input/output and text output. Text output above the WASM
text-input cap needs the Rust/Python reader. Aggregate record collection is
also bounded before section-row allocations can grow without limit.

Errors include one-based source lines where available, or zero for record
validation/global limits. CRLF, blank separator lines, case-insensitive section
names and Fortran `D`/`d` exponents are accepted. Dummy/ghost atoms, arbitrary
atom labels and trailing annotations are outside this profile.

## Bounded Open Babel comparison

The original MIT packet contains 17 accepted profile cases and 29 refusals.
There are 12 complete vibration cases, three atom-geometry cases and two partial
inspection cases. Rust/Python/Node/CLI test the same expectations, plus edited
data. Additional tests cover malformed matrices, total displacement limits,
unknown JSON metadata, finite f64 extremes and a 10,000-atom JSON round trip.

The local 2026-10-24 comparison uses Open Babel **3.2.1** with
`-imolden -omolden -ab -l1`: no inferred bonds, one input record, 30 seconds per
process. An independent output-table observer compares fields present in the
source. Derived atom geometry added to a vibration-only file is allowed.
Unit changes are tracked as encoding changes, not automatically called physical
errors. Source-absent optional fields are excluded from preservation denominators.

Both engines emitted observable records for all **12 complete vibration cases**:

| Present source field, complete vibration cases only | chematic exact | Open Babel exact |
|---|---:|---:|
| Bohr vibration coordinates | 12/12 | 8/12 |
| Frequency values | 12/12 | 8/12 |
| Indexed mode blocks and Bohr displacements | 12/12 | 8/12 |
| IR intensities | 11/11 | 0/11 |
| Present Raman intensities | 5/5 | 0/5 |

Open Babel omits `[INT]` in these outputs and rounds high-precision frequency
and displacement values. Other differences, including zero-frequency handling,
mode IDs, duplicate coordinate channels and Fortran exponent input, remain
visible in per-case receipts. Exact value equality is a source-preservation
contract, not a claim that every changed value is chemically significant.

Across all 17 profile inputs, chematic preserves every source field in 17/17.
Open Babel emits 16/17 observable records and preserves every measured source
field in 1/17. Its rejection of the partial frequency/geometry-without-modes
case is recorded separately; it is not evidence of nonconformance for a complete
Molden vibration file. This packet does not establish overall superiority or
physical correctness, and says nothing about orbital support.

```sh
python3 scripts/validation/check_molden_preservation.py \
  --cli target/debug/chematic --obabel /path/to/obabel --output /tmp/molden-comparison.json
```

Inputs, expectations and receipts are under
[`validation/molden_vibrations`](https://github.com/kent-tokyo/chematic/tree/main/validation/molden_vibrations).
Receipts contain commands, source/input/binary hashes, stdout/stderr and terminal
states. They describe a local candidate, not a released artifact or hosted CI.
The implementation follows the public
[Molden format description](https://www.theochem.ru.nl/molden/molden_format.html);
the comparator options are documented on Open Babel's
[Molden page](https://openbabel.org/docs/FileFormats/Molden_format.html).
No Open Babel or Molden implementation code was copied.
