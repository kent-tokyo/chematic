# ViewMol geometry (`vmol`)

The checked ViewMol geometry reader/writer adds a format previously absent from
chematic. Rust, Python, WASM/Node and CLI retain an optional title, raw coordinate
conversion factor, raw coordinate decimals, element order and fixed-axis flags.
An explicit checked API projects positions to Angstrom without changing the
editable source record.

[Open Babel lists ViewMol](https://openbabel.org/docs/FileFormats/ViewMol_format.html).
The [ViewMol input-filter specification](https://viewmol.sourceforge.net/documentation/node32.html)
defines the coordinate factor as a multiplier to Angstrom and the optional
`x`/`y`/`z` combination as directions fixed during optimization. This bounded
implementation and its fixtures were authored independently; no external
implementation or fixture corpus is copied.

## Python

```python
import chematic

text = """$title
retained geometry
$coord 2.0
.125 .25 .5 C xz
$end
"""
record = chematic.parse_viewmol(text)  # UTF-8 bytes also accepted
assert record["factor_to_angstrom"] == "2.0"
assert record["atoms"][0]["fixed_axes"] == "xz"
assert chematic.viewmol_positions_angstrom(record) == [[.25, .5, 1.0]]
record["atoms"][0]["position"][0] = "-.12345678901234567"
assert chematic.parse_viewmol(chematic.write_viewmol(record)) == record
```

The JSON-compatible record contains optional `title`, `factor_to_angstrom` and
ordered `atoms`. Each atom has `element`, three `position` strings in source
units and optional `fixed_axes`. Coordinates and the factor remain decimal
strings, preserving precision, exponent spelling and signed zero through edits.
Unknown JSON fields and recursively duplicate keys, including escaped aliases,
refuse. JSON numbers, booleans or nulls cannot replace numeric strings.

## Supported profile

| Data | Contract |
|---|---|
| `$coord factor` | Exactly one coordinate group; explicit finite positive factor required |
| Atom rows | Three raw coordinate tokens, canonical element symbol, optional fixed-axis token |
| Elements | Standard case, e.g. `C`, `Cl`, `Og`; labels, pseudoatoms and case aliases refuse |
| Fixed axes | A nonempty combination of distinct lowercase `x`, `y`, `z`; token order retained |
| `$title` | Optional; exactly one following physical line, including an empty line |
| Title state | Absent/null differs from an explicitly empty title; no filename is invented |
| Group order | Title may precede or follow the single contiguous coordinate block |
| `$end` | Required final nonblank group; later records/groups/data refuse |
| Layout | Spaces/tabs, blank layout lines and CRLF accepted; writer uses canonical spaces and LF |

Titles retain leading/trailing spaces and UTF-8, but control characters except
tab and a title beginning with a data-group marker refuse. Empty fixed-axis
strings, duplicate flags, uppercase flags and separate axis columns refuse.

The writer emits the title first when present, then unchanged factor/coordinate
strings, optional flags and `$end` with a final newline. It serializes current
fields; it does not replay stale source text. Group order, whitespace and line
endings are normalized. Atom order and represented data stay intact.

Numbers use an optional sign, integer/decimal mantissa and optional `e`/`E`
signed exponent. Their nominal binary64 value must be finite; nonzero mantissas
rounding to zero refuse. Representable subnormals and signed zero are supported.
Fortran `D` exponents, uncertainty suffixes and trailing numeric garbage refuse.

Cell, vibrational, gradient/trajectory, orbital/basis, energy and other data
groups are outside this profile and refuse. There is no bond perception,
chemical-state inference, valence checking, cell expansion or common `Molecule`
conversion. Unknown groups are never stripped to manufacture a geometry-only
success. General ViewMol projects and TURBOMOLE control files are not supported.

## Explicit Angstrom projection

Raw read/write preserves the source factor and coordinates. Projection computes
binary64 `raw * factor` separately. It rejects a nonfinite result or a nonzero
coordinate rounded to zero during multiplication. A raw record can therefore
round-trip successfully while its requested binary64 projection refuses.

For example the finite raw pair `factor=2`, `x=1.7976931348623157e308` exceeds
the binary64 range when multiplied; `factor=.5`, `x=4.9406564584124654e-324`
rounds to zero. Both raw records remain serializable. These are deliberate
numerical stress cases, not plausible molecular geometries. No conversion occurs
implicitly during ordinary writing. Projection returns only coordinate arrays;
use the record for titles, species and constraints.

## APIs and limits

- **Rust:** `parse_viewmol`, `parse_viewmol_with_limits`, `parse_viewmol_json`,
  `write_viewmol`, `write_viewmol_with_limits`; `ViewmolRecord`, `ViewmolAtom`,
  `ViewmolParseLimits`, `ViewmolError`. Records expose `positions_angstrom()` and
  `positions_angstrom_with_limits(...)`.
- **Python:** `parse_viewmol(text, ...)`, `write_viewmol(record, ...)`,
  `viewmol_positions_angstrom(record, ...)`. Each accepts keyword limits
  `max_input_bytes=None`, `max_atoms=None`, `max_line_len=None`.
- **WASM / Node:** `parse_viewmol_json(text)`, `write_viewmol_from_json(json)`,
  `viewmol_positions_angstrom_json(json)`.
- **CLI:** `chematic viewmol -i input.vmol -o output.vmol`; `--inspect` emits
  editable JSON, `--from-json` reads it, and `--angstrom` explicitly emits only
  projected positions as JSON. `--from-json --angstrom` is supported; `--inspect`
  conflicts with both flags.

Native/Python defaults are 64 MiB input/text output, 4,096 bytes per physical
line, one million atoms and 256 bytes per numeric token. Writer byte limits use
`max_input_bytes` too. Strict JSON input uses that byte limit; its encoding may
be larger than equivalent text. CLI exposes `--max-input-bytes`, native atom/line
defaults and its shared 64 MiB output cap. WASM fixes input/text output at
1,000,000 bytes, atoms at 10,000 and JSON output at 16 MiB; line/token bounds
match native defaults.

Native errors carry line/detail; line zero denotes record-wide validation.
Python raises `ValueError`, WASM throws a format error and CLI exits 2 without
partial stdout. Python invalid UTF-8/type boundaries remain
`UnicodeDecodeError`/`TypeError`.

## Local validation and comparison

The [original packet](https://github.com/kent-tokyo/chematic/blob/main/validation/viewmol_preservation/cases.json)
contains **28 accepted records, 53 refused inputs and 2 raw records whose
projection refuses**. It includes explicit unit factors, absent/empty/spaced and
Unicode titles, both title positions, all seven nonempty fixed-axis sets plus
token-order variants, source atom order, heavy elements, high/scientific/long
precision, subnormals, signed zero and layout variations.

The [raw observations](https://github.com/kent-tokyo/chematic/blob/main/validation/viewmol_preservation/observed.json)
include source/runtime hashes, commands and outputs. Python, actual WASM under
Node and native CLI each retain all 28 records through three read/write cycles,
with the initial and three subsequent Angstrom projections matching nominal
binary64 bits. All 53 refusals and both projection refusals are explicit. The
latter two records still retain their raw fields through ordinary read/write.

Installed **Open Babel 3.2.1** is observed through seekable file-based
**`vmol → vmol`** rewriting with `-ab` to disable bond perception. An independent
data-group observer inspects its output. Eighteen of 28 original inputs produce
an observable record. The nine fixed-axis cases and one blank-layout case
produce none (`0 molecules converted`); they are unobservable, not scored as
safe refusals or silent-loss cases. The fixed-axis syntax comes from ViewMol's
own specification; these observations identify a bounded coverage gap in the
installed Open Babel workflow.

Of the 18 observable records, physical Angstrom position bits are exact in
**14/18**, element order in **18/18**, and title state/text in **15/18**. Eleven
retain all declared physical/metadata fields together. High/scientific/long
precision and subnormal position cases change; absent/empty/spaced title states
also change. All observable inputs lack fixed axes, so their matching empty
constraint state does not demonstrate constraint handling by Open Babel.

Raw factor/coordinate spelling is a separate representation metric. Changing
the factor to 1 while scaling coordinates correctly is **not** counted as
physical meaning loss. Neither a different equivalent axis-token order nor
numeric spelling alone counts as physical loss. The observer has ten controls
for normalization, scale loss, fixed-axis loss/order, signed zero, title/species
changes, missing output and crash/refusal accounting.

Regenerated candidate text yields 19/28 observable Open Babel rewrites; the nine
fixed-axis records remain unobservable while canonicalized blank layout is
accepted. A separate unit probe confirms that raw `(.125, .25, .5)` with factor
2 becomes Angstrom `(.25, .5, 1)`. The ViewMol application itself was not run.

This is a small authored local-source packet, not a production corpus or proof
of all ViewMol semantics, isolated parser/writer accuracy, speed, released-package
behavior or broad toolkit superiority. Unobservable external outputs remain in
the raw receipts and are excluded from preservation denominators.

```sh
python3 scripts/validation/check_viewmol_preservation.py \
  --python-path /path/to/staged-python \
  --wasm-module /path/to/node/chematic_wasm.js \
  --cli /path/to/chematic \
  --openbabel /path/to/obabel \
  --output /tmp/viewmol-observed.json
```
