# Ghemical GPR records

The checked GPR reader/writer adds a format previously absent from chematic.
Rust, Python, WASM/Node and CLI support the **`!Header gpr 100` / `!Info 1`**
profile: atom IDs/numbers, literal bond codes, one Cartesian coordinate set in nm
and optional partial charges. Decimal tokens stay strings, preserving source
precision and signed zero when edited and written.

[Open Babel documents Ghemical support](https://openbabel.org/docs/FileFormats/Ghemical_format.html).
The upstream [Ghemical example](https://sources.debian.org/src/ghemical/3.0.0-5/examples/bromobenzene.gpr/)
shows the version-100 section layout. This implementation and its fixtures were
authored independently; no external implementation or fixture corpus is copied.
The public examples establish a layout, not a complete specification of every
Ghemical dialect. Acceptance is limited to the profile below.

## Python example

```python
import chematic

text = """!Header gpr 100
!Info 1
!Atoms 2
0 6
1 8
!Bonds 1
0 1 S
!Coord
0 0.12345678901234567 0 -0.0
1 .25 0 0
!Charges
0 .125
1 -.125
!End
"""
record = chematic.parse_gpr(text)  # UTF-8 bytes also accepted
assert record["coordinates"][0]["position_nm"][0] == "0.12345678901234567"
assert record["partial_charges"][1]["value"] == "-.125"
record["coordinates"][0]["position_nm"][0] = "0.23456789012345678"
assert chematic.parse_gpr(chematic.write_gpr(record)) == record

# Preserve an absent charge declaration; do not calculate or invent zero values.
record["partial_charges"] = None
assert "!Charges" not in chematic.write_gpr(record)
```

## Profile and meaning

| Section / field | Contract |
|---|---|
| Header | Exactly `!Header gpr 100` and `!Info 1` |
| `!Atoms N` | Positive count, followed by `id atomic_number`; numbers 1..118 |
| Atom IDs | Consecutive 0..N-1 in definition order; sparse/permuted IDs refuse |
| `!Bonds M` | Nonnegative count; `begin_id end_id code` rows |
| Bond codes | Literal `S`, `D`, `T`, `C`; retain row order and endpoint orientation, without perception or Kekulé conversion |
| Bond validity | Existing distinct endpoints, at most one unordered pair |
| `!Coord` | Exactly N `atom_id x y z` rows, in consecutive atom-ID order |
| Coordinate units | Raw Cartesian nm; no conversion during read/write |
| `!Charges` | Optional; if present, exactly N `atom_id value` rows in consecutive atom-ID order |
| Partial charges | Declared values only; absent stays absent, including when edited |
| `!End` | Required; no further nonblank rows or records |

The JSON-compatible record has `atoms`, `bonds`, `coordinates` and optional
`partial_charges` (null/omitted means absent). Atom IDs and atomic numbers are
JSON integers; bond codes are strings. Coordinates and partial-charge values
are decimal **strings**, not JSON numbers. Tokens may have an optional sign,
integer/decimal mantissa and `e`/`E` integer exponent. Their binary64 nominal
value must be finite; a nonzero mantissa that rounds to zero refuses.
Representable subnormals and signed zero are supported.

The token remains authoritative even when its decimal precision exceeds
binary64. Rust `gpr_nominal_value(token)` returns its checked, rounded nominal
value without changing units. No numeric reformatting occurs during writing.

Partial charges are not formal charges. No formal charge, isotope, spin,
implicit hydrogen count or stereochemistry is inferred. This is an editable
format record, separate from the common `Molecule` conversion bridge. It does
not promise valence checking or chemical validity from literal bond codes.

Other versions, the old `mm1gp` header, multiple coordinate sets, sparse or
permuted atom/coordinate/charge rows, extra columns, comments, force/energy/cell
extensions and unknown sections refuse. Unknown JSON keys, recursively duplicate
keys (including escaped aliases), numeric coordinate/charge values and invalid
IDs also refuse. These restrictions prevent unrepresented fields from being
silently discarded.

Blank physical lines and tabs/CRLF are accepted; writing normalizes layout to
spaces and LF with a final newline. Original whitespace and integer spelling
are not retained. Current edited fields are written, never stale source text.

## APIs and bounds

- **Rust:** `parse_gpr`, `parse_gpr_with_limits`, `parse_gpr_json`, `write_gpr`,
  `write_gpr_with_limits`, `gpr_nominal_value`; `GprRecord`, `GprAtom`, `GprBond`,
  `GprBondOrder`, `GprCoordinate`, `GprPartialCharge`, `GprParseLimits`, `GprError`.
- **Python:** `parse_gpr(text, *, max_input_bytes=None, max_atoms=None,
  max_bonds=None, max_line_len=None)` and `write_gpr(record, ...)` with the same
  keyword limits.
- **WASM / Node:** `parse_gpr_json(text)` returns editable record JSON;
  `write_gpr_from_json(json)` emits checked text.
- **CLI:** `chematic gpr -i input.gpr -o output.gpr`; `--inspect` emits JSON,
  `--from-json` reads JSON and emits GPR. These flags conflict.

Native defaults are 64 MiB input/text output, 4,096 bytes per physical line,
one million atoms, two million bonds and 256 bytes per decimal token. Python
exposes the same overrides. CLI exposes `--max-input-bytes` and uses native
atom/bond/line defaults plus its shared 64 MiB output cap. Strict JSON input
uses the byte limit too, and can require more bytes than equivalent text.
WASM fixes input/text-output at 1,000,000 bytes, atoms at 10,000, bonds at
20,000 and JSON output at 16 MiB; line/token bounds match native defaults.

Native errors carry line/detail; line zero denotes record-wide validation.
Python raises `ValueError`, WASM throws a format error, and CLI exits 2 without
partial stdout. Invalid UTF-8/type boundaries in Python retain their normal
`UnicodeDecodeError`/`TypeError`.

## Local verification and Open Babel comparison

The [original packet](https://github.com/kent-tokyo/chematic/blob/main/validation/gpr_preservation/cases.json)
contains **16 positive and 55 refusal** cases. It covers all four bond codes,
reversed endpoints, source bond order, absent versus explicit-zero charges,
nonintegral charges, long/scientific decimals, signed zero, subnormals, heavy
and repeated elements, tabs/CRLF and blank lines. Separate tests exercise edits,
strict JSON and both-direction bounds.

The [raw observation](https://github.com/kent-tokyo/chematic/blob/main/validation/gpr_preservation/observed.json)
contains source/runtime hashes, tool versions and commands/output. Each candidate
lane checks the initial parsed fields and all three subsequent read/write cycles:

| Local workflow | Entire token record retained | Explicit refusals |
|---|---:|---:|
| Python native extension | 16/16 | 55/55 |
| Actual WASM under Node | 16/16 | 55/55 |
| Native CLI | 16/16 | 55/55 |

Installed **Open Babel 3.2.1** is observed in a fresh-process **GPR → GPR**
workflow. A small independent section/field observer reads its output. Fifteen
of the 16 original inputs produce an observable record; the leading/section
blank-line case produces none (`0 molecules converted`). That case is
unobservable, not counted as safe refusal or a field-loss score.

On the 15 observable cases, the atom/bond graph is exact in **15/15**,
coordinate nominal binary64 bits in **11/15**, and partial-charge declaration
presence plus nominal bits in **11/15**. Ten retain all three together. One
no-charge input receives an explicit zero-charge section. High-precision,
scientific and long-decimal coordinates/charges change; subnormal coordinates
also change. Changes in token spelling alone are reported as a separate lexical
metric, not numeric meaning loss. Bond orientation/order normalization alone
is not counted as graph loss.

All **16 regenerated candidate outputs** also produce observable Open Babel
GPR output. This is an acceptance check, not isolated reader accuracy. A separate
unit probe feeds nm `(.125, .25, .5)` and observes Angstrom `(1.25, 2.5, 5)`
in Open Babel XYZ output. The Ghemical application itself was not executed.

These are small authored fixtures and a local source candidate. They do not
establish full Ghemical project/dialect support, parser-only or writer-only
accuracy, speed, released-package behavior or broad toolkit superiority.

Reproduce with freshly built Python, Node WASM and CLI artifacts:

```sh
python3 scripts/validation/check_gpr_preservation.py \
  --python-path /path/to/staged-python \
  --wasm-module /path/to/node/chematic_wasm.js \
  --cli /path/to/chematic \
  --openbabel /path/to/obabel \
  --output /tmp/gpr-observed.json
```

Eight observer controls check numeric/charge/bond changes, invented charges,
signed zero, harmless endpoint reversal, missing output and crash/refusal
classification. Partial stdout, panic exit 101, signals and wrong error prefixes
cannot count as typed refusals.
