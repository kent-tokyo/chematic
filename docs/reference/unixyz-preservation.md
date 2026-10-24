# UniChem XYZ (`unixyz`)

The checked UniChem XYZ reader/writer adds a format previously absent from
chematic. Rust, Python, WASM/Node and CLI expose a single-record editable model
that retains the first-line title/comment, source atom order, atomic numbers
and Cartesian coordinate decimal strings in Angstrom.

[Open Babel lists UniChem XYZ](https://openbabel.org/docs/FileFormats/UniChem_XYZ_format.html).
Its [XYZ reference](https://openbabel.org/docs/FileFormats/XYZ_cartesian_coordinates_format.html)
identifies UniChem's atomic-number column. The UniChem title/count order below
is verified against installed Open Babel's `xyz → unixyz` output; the ordinary
XYZ page's count/title order describes ordinary XYZ. No complete UniChem
specification is assumed. This implementation and its fixtures were independently
authored without copying an external implementation or fixture corpus.

## Python example

```python
import chematic

text = """retained title
2
6 0.12345678901234567 -0.0 .25
8 .5 .75 1
"""
record = chematic.parse_unixyz(text)  # UTF-8 bytes also accepted
assert record["title"] == "retained title"
assert record["atoms"][0]["position_angstrom"][0] == "0.12345678901234567"
record["title"] = "edited title"
record["atoms"][0]["position_angstrom"][0] = "-.23456789012345678"
assert chematic.parse_unixyz(chematic.write_unixyz(record)) == record
```

Coordinates are strings so editing and JSON transport retain the original
decimal precision and signed zero. `atomic_number` is an integer. The record
contains exactly `title` and `atoms`; each atom has `atomic_number` and
`position_angstrom`. The output count is regenerated from current atoms.

## Supported profile

| Field | Contract |
|---|---|
| First line | UTF-8 title/comment; empty, numeric, leading/trailing spaces and tabs retained |
| Second line | Positive decimal atom count, with optional surrounding spaces/tabs |
| Atom rows | Exactly the declared count; atomic number plus three space/tab-separated coordinates |
| Atomic numbers | Integers 1..118; element symbols and pseudoatom zero refuse |
| Coordinates | Raw Cartesian Angstrom; no unit conversion |
| Numeric grammar | Optional sign, integer/decimal mantissa, optional `e`/`E` signed integer exponent |
| Output | Current fields, canonical integer spelling, single atom-row spaces, LF and final newline |

The title is not replaced by a filename. Each coordinate token stays verbatim,
including long decimal precision, exponent spelling and signed zero. Rust
`unixyz_nominal_value(token)` returns its checked, rounded binary64 nominal
value while the decimal string stays authoritative.

Nominal coordinates must be finite. Nonzero mantissas rounding to zero refuse;
representable subnormals and signed zero are supported. Fortran `D` exponents,
uncertainty suffixes, hexadecimal values and numeric suffix garbage refuse.
The format profile does not check whether a coordinate is chemically plausible.

One complete record is supported. Missing rows, extra atoms, second records,
blank atom rows and every extra trailing line, including blank lines, refuse.
A normal final newline is allowed, as is no final newline. Atom rows are counted
physically and are never skipped to repair an incorrect count.

Connectivity, charge, isotope, spin, stereochemistry, cell/units extensions and
extra columns are outside this profile. Unknown JSON fields and recursively
duplicate keys, including escaped aliases, refuse. JSON coordinate numbers,
booleans and nulls refuse; callers must supply validated decimal strings.
No bonds, formal charges, implicit H counts or other chemical state are inferred.
This record has a separate API from the common `Molecule` conversion bridge.

## APIs and limits

- **Rust:** `parse_unixyz`, `parse_unixyz_with_limits`, `parse_unixyz_json`,
  `write_unixyz`, `write_unixyz_with_limits`, `unixyz_nominal_value`;
  `UnixyzRecord`, `UnixyzAtom`, `UnixyzParseLimits`, `UnixyzError`.
- **Python:** `parse_unixyz(text, *, max_input_bytes=None, max_atoms=None,
  max_line_len=None)` and `write_unixyz(record, ...)` with the same limits.
- **WASM / Node:** `parse_unixyz_json(text)` returns editable record JSON;
  `write_unixyz_from_json(json)` returns checked text.
- **CLI:** `chematic unixyz -i input.unixyz -o output.unixyz`; `--inspect`
  emits JSON, `--from-json` reads it and emits UniChem text. These flags conflict.

Native/Python defaults are 64 MiB input/text output, 4,096 bytes per physical
line, one million atoms and 256 bytes per coordinate token. Writer byte limits
use the `max_input_bytes` setting too. Strict JSON input uses the same byte
limit, and its encoding may be larger than equivalent text. CLI exposes
`--max-input-bytes`, native atom/line defaults and its shared 64 MiB output cap.
WASM fixes input/text output at 1,000,000 bytes, atoms at 10,000 and JSON output
at 16 MiB; line/token limits match native defaults.

Native errors have a line/detail; zero denotes record-wide validation. Python
raises `ValueError`, WASM throws a format error and CLI exits 2 without partial
stdout. Python invalid UTF-8/type errors remain `UnicodeDecodeError`/`TypeError`.

## Local validation and Open Babel comparison

The [original packet](https://github.com/kent-tokyo/chematic/blob/main/validation/unixyz_preservation/cases.json)
has **15 positive and 41 refusal** cases, covering title variants, atomic order,
heavy elements, scientific/long decimals, signed zero, subnormals, finite
extremes, CRLF/tabs and missing final newline. Separate tests cover edited fields,
strict JSON, both-direction limits and the distinct ordinary-XYZ header layout.
These are small authored cases, not a production corpus.

The [raw observation](https://github.com/kent-tokyo/chematic/blob/main/validation/unixyz_preservation/observed.json)
records source/runtime hashes, versions, commands and outputs. Each candidate
lane checks initial fields and all three successive read/write cycles:

| Workflow | Entire token record retained | Explicit refusals |
|---|---:|---:|
| Python native extension | 15/15 | 41/41 |
| Actual WASM under Node | 15/15 | 41/41 |
| Native CLI | 15/15 | 41/41 |

Installed **Open Babel 3.2.1** is observed through fresh-process
**`unixyz → unixyz`** rewriting with `-ab` to disable bond perception. No external
title is injected. An independent field splitter observes the output.
Fourteen of the 15 original inputs produce an observable atom record. On those
14, atomic numbers/order are exact in **14/14**, coordinate nominal binary64
bits in **10/14**, and first-line title in **0/14**: this workflow writes the
input path instead of the supplied first line. Precision, scientific,
subnormal and long-decimal coordinate cases change. Numeric token spelling
alone is a separate lexical metric and is not scored as numeric meaning loss.

The finite-extreme stress case completes with `1 molecule converted`, but its
output atom row has joined fields that cannot be split into atomic number and
three coordinates. This is recorded as **unobservable external output**, not a
safe refusal or numeric loss score. Its coordinates are deliberately extreme
format stress data, not a plausible molecular geometry. It remains visible in
the packet and all raw receipts.

Regenerated candidate text also yields 14/15 observable same-format Open Babel
outputs, with the same finite-extreme output limitation. A separate
`unixyz → CJSON` probe observes all **15/15** regenerated inputs with matching
atomic number order and finite coordinate arrays of the correct shape. That
probe avoids relying on Open Babel's UniChem output columns for input acceptance;
it does not establish isolated reader accuracy or coordinate precision.

The layout probe observes `title`, then count, then atomic-number rows in
Open Babel's `xyz → unixyz` output. A unit probe observes `(.125, .25, .5)`
unchanged in its XYZ output, confirming the Angstrom convention on that input.
The observer gate requires the candidate contract, layout/unit probes and CJSON
read-shape check; external same-format unobservable outputs are retained and
excluded from preservation denominators.

This is a local source candidate and a small same-format workflow comparison.
It does not establish every UniChem dialect, multi-record support, general
parser/writer superiority, speed, released-package behavior or broad toolkit
superiority.

Reproduce with freshly built artifacts:

```sh
python3 scripts/validation/check_unixyz_preservation.py \
  --python-path /path/to/staged-python \
  --wasm-module /path/to/node/chematic_wasm.js \
  --cli /path/to/chematic \
  --openbabel /path/to/obabel \
  --output /tmp/unixyz-observed.json
```

Eight observer controls check changed coordinates/titles/atomic numbers, signed
zero, harmless numeric spelling changes, empty/truncated/joined/extra output
and crash/refusal classification. Panic exit 101, signals, wrong error prefixes
and partial stdout cannot count as typed refusals.
