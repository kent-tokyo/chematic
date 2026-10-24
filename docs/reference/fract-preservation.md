# Free Form Fractional (`fract`)

The checked `fract` reader/writer adds a single-cell crystallographic format
previously absent from chematic. Rust, Python, WASM/Node and CLI expose the same
editable record. Decimal strings retain precision, signed zero and trailing
parenthesized uncertainty digits through serialization.

The format description comes from the
[Open Babel fract documentation](https://openbabel.org/docs/FileFormats/Free_Form_Fractional_format.html).
The implementation and validation cases were authored independently; no Open
Babel implementation or external fixture corpus was copied.

## Python

```python
import chematic

text = """example
4.000(3), 5, 6, 90, 90, 90
C 0.12345678901234567(32) -0.0 1.25
"""
record = chematic.parse_fract(text)  # UTF-8 bytes also accepted
assert record["cell"][0] == "4.000(3)"
assert record["atoms"][0]["fractional"][0] == "0.12345678901234567(32)"
record["atoms"][0]["fractional"][2] = "2.75(2)"
assert chematic.parse_fract(chematic.write_fract(record)) == record
```

The JSON-compatible record contains exactly `title`, six ordered `cell` tokens
and ordered `atoms`. Each atom has a canonical `element` symbol and exactly
three `fractional` tokens. Numeric fields are strings, allowing the decimal and
its uncertainty suffix to stay together when edited. Python numbers, booleans,
nulls, unknown fields and unsupported chemical annotations refuse.

## Supported profile

| Field | Contract |
|---|---|
| Title | First physical line; UTF-8, including leading/trailing spaces, tabs or empty title |
| Cell | Second line: a, b, c, alpha, beta, gamma; whitespace and/or comma delimiters |
| Cell units | Lengths in Angstrom; angles in degrees |
| Atom rows | Canonical element plus three whitespace-separated fractional coordinates |
| Coordinate units | Dimensionless fractions of the cell basis, including values outside [0, 1] |
| Decimal grammar | Optional sign; integer or decimal mantissa; optional `e`/`E` signed integer exponent |
| Uncertainty | Optional trailing `(ASCII digits)`; retained verbatim without interpreting an absolute standard deviation |
| Output | Current edited fields, space delimiters, LF and final newline; every number token retained verbatim |

For example `1.23e-2(04)` is supported; `1.23(04)e-2` and Fortran `D` exponents
are outside this profile. Arbitrary spacing and blank atom-section lines are
accepted and normalized on write. Original source bytes are not retained.

Nominal numbers must be finite when read as binary64. Nonzero mantissas that
round to zero refuse; representable subnormals and signed zero are supported.
`FractNumber::value()` gives the rounded binary64 nominal value, while
`token()` remains authoritative for the original decimal precision.
`uncertainty_digits()` exposes the uninterpreted suffix digits.

Cell lengths must be positive; angles must be strictly between 0 and 180
degrees. The shared CIF `UnitCell` validator requires a dimensionless Gram
determinant above `1e-12` and finite positive volume. This rejects inconsistent,
singular, numerically near-singular and overflowing cells. No default cell is
invented. Atom count must be positive.

Connectivity, formal/partial charges, isotopes, site labels, occupancy, symmetry,
multiple records and extra columns are outside this profile and refuse when
present. No wrapping, basis rotation, coordinate conversion, symmetry expansion
or bond inference occurs. The common molecular conversion bridge does not
accept this record.

## APIs and limits

- **Rust:** `parse_fract`, `parse_fract_with_limits`, `write_fract`,
  `write_fract_with_limits`, `parse_fract_json`; `FractRecord`, `FractAtom`,
  `FractNumber`, `FractParseLimits`, `FractError` in `chematic_mol`.
- **Python:** `parse_fract(text, *, max_input_bytes=None, max_atoms=None,
  max_line_len=None)` and `write_fract(record, *, max_input_bytes=None,
  max_atoms=None, max_line_len=None)`.
- **WASM / Node:** `parse_fract_json(text)` returns record JSON;
  `write_fract_from_json(json)` emits validated text.
- **CLI:** `chematic fract -i input.fract -o output.fract`; `--inspect` emits
  editable JSON, `--from-json` reads it. These flags conflict.

Rust/Python defaults are 64 MiB serialized input/output, 4,096 bytes per physical
line, one million atoms and 256 bytes per numeric token. Limit overrides apply
in both read/write directions. CLI has the same defaults and
`--max-input-bytes`; its shared output cap is 64 MiB. JSON input uses the byte
limit too, so its larger encoding may need a larger limit than equivalent text.
WASM fixes input/text-output at 1,000,000 bytes, atoms at 10,000 and JSON output
at 16 MiB; its line/token limits match native defaults.

Strict JSON input rejects duplicate keys at every depth, including escaped
aliases, as well as unknown fields. Native errors carry a line and detail;
Python raises `ValueError`, WASM throws a format error and CLI exits 2 without
partial stdout. Record-wide errors use line zero. A JSON or writer field diagnostic uses
logical record rows; text reader diagnostics use physical source rows.

## Local validation and comparison

The [original packet](https://github.com/kent-tokyo/chematic/blob/main/validation/fract_preservation/cases.json)
contains 17 accepted and 43 refused cases. It covers orthogonal/triclinic cells,
comma/tab layouts, Unicode/empty titles, source order, unwrapped coordinates,
long decimal tokens, scientific notation, subnormals, signed zero and uncertainty
digits. Separate tests cover edited output, unknown/duplicate JSON keys and
limits. These are small authored fixtures, not a production corpus.

The [raw observation](https://github.com/kent-tokyo/chematic/blob/main/validation/fract_preservation/observed.json)
records local runtime/source hashes, commands and outputs. Each candidate lane
checks the initial record and all three successive parse/write cycles:

| Workflow | Full token record preserved | Explicit refusals |
|---|---:|---:|
| Python native extension | 17/17 | 43/43 |
| Actual WASM under Node | 17/17 | 43/43 |
| Native CLI | 17/17 | 43/43 |

Installed Open Babel **3.2.1** is observed through fresh-process `fract → fract`
conversion with `-ab` to disable bond perception. An independent field splitter,
not the candidate reader, observes its output. Of 17 inputs, 16 produce an
observable record. The blank-atom-line case produces no record and reports
`0 molecules converted`; this is recorded as unobservable, not safe refusal or
a silent-loss score.

On those 16 observable cases, Open Babel retains nominal binary64 number bits
in 10, title/element order in 13, and uncertainty tokens in **0 of 3** applicable
cases. Seven retain all declared fields together. The full original number-token
spelling is a separate lexical metric; changes in spelling alone are not counted
as numeric meaning loss. Open Babel documents five-decimal output and ignored
uncertainty suffixes. The observation includes an orthogonal-cell unit probe:
cell `(4, 5, 6)` and fractional `(.25, .5, .75)` become Cartesian Angstrom
`(1, 2.5, 4.5)` in its XYZ output.

These results concern a same-format read/write workflow. They do not isolate
Open Babel's parser from its writer, establish speed, demonstrate all `fract`
dialects, or prove broad crystallographic/toolkit superiority. This feature is a
local source candidate until its PR is merged and released.

Reproduce after building fresh Python, Node WASM and CLI artifacts:

```sh
python3 scripts/validation/check_fract_preservation.py \
  --python-path /path/to/staged-python \
  --wasm-module /path/to/node/chematic_wasm.js \
  --cli /path/to/chematic \
  --openbabel /path/to/obabel \
  --output /tmp/fract-observed.json
```

The observer distinguishes format errors from crashes, wrong error prefixes,
partial stdout and missing external output. It also feeds regenerated candidate
text to Open Babel as a separate acceptance check: all 17 produce an observable
record. This does not claim isolated
reader accuracy or override the numeric rewrite results above.
