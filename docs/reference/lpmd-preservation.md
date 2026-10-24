# LPMD trajectory preservation

Rust, Python, WASM/Node and the CLI read and write the bounded **LPMD 2.0 L**
text profile. This adds an adapter absent from the checked base `b2d6d3a4`.
All frames, atom order, column order, oriented cell vectors and original decimal
tokens remain editable data. The writer uses current fields; it never replays
stale source text. This is field preservation, not byte-identical source replay.

## Profile and meaning

The first two lines are `LPMD 2.0 L` and `HDR <columns>`. Each subsequent frame
has a positive atom count, nine cell-vector components and exactly that many
atom rows. The document contains at least one frame; atom counts and cells can
change between frames. The header applies to every frame.

Required columns are `SYM`, `X`, `Y`, `Z`, each exactly once. `SYM` uses canonical
element symbols through Og, without label, isotope or pseudoatom inference.
`X/Y/Z` are dimensionless fractional coordinates. Cell vectors are **rows** in
angstroms; they are never reconstructed from cell lengths and angles. Explicit
Cartesian projection uses `x*a + y*b + z*c`, retaining the supplied orientation,
including left-handed cells. Unwrapped fractional coordinates are retained;
neither raw writing nor projection wraps them into a unit cell.

Optional columns are complete `VX/VY/VZ`, `AX/AY/AZ` and `FX/FY/FZ` triplets;
`CHG`, `PHG`, `FHG`; `rgb` or `RGB`; and scalar color `C`. Their source tokens
are preserved independently, including when both acceleration and force
triplets occur. There is **no** force/acceleration/mass conversion, velocity
rescaling or interpretation of charge columns as molecular formal charges.
The charge/force labels include Open Babel extensions; this does not establish
support by the LPMD simulation engine. Colors are `<r,g,b>` or scalar `C`, with
finite components in `[0,1]`. Column order may vary; duplicate names and the
simultaneous `rgb`/`RGB` aliases refuse.

All numeric fields are ASCII decimal strings with optional signs and `e/E`
exponents. This preserves precision, signed zero and exponent spelling without
a binary64 round trip. Each token is at most 256 bytes, its nominal binary64
value must be finite, and a nonzero mantissa rounding to zero refuses. Finite
subnormals are accepted. Cell vectors must be nonzero and their row-normalized
binary64 determinant nonzero; numerically degenerate cells refuse. Projection
additionally refuses nonfinite products/sums and nonzero products underflowing
to zero. A record unsafe to project can still be read and written unchanged.

Compressed `Z`, LPMD 1.0, six-parameter cells, missing or incomplete required
columns, partial vector triplets, unknown columns, malformed rows, comments,
extra headers and incomplete/trailing frames refuse. Tabs, CRLF, a missing final
newline and trailing blank layout lines are accepted. Blank lines inside a
frame or before another frame refuse. No topology, periodic-boundary flags,
timestamps, simulation units beyond the stated coordinate convention, or
molecular chemistry are fabricated.

## APIs

```python
from pathlib import Path
import chematic

record = chematic.parse_lpmd(Path("input.lpmd").read_bytes())
# Numeric fields remain strings. Rows follow record["columns"].
x = record["columns"].index("X")
record["frames"][0]["atoms"][0][x] = "0.12500000000000000001"
Path("output.lpmd").write_text(chematic.write_lpmd(record), encoding="utf-8")
cartesian_frames = chematic.lpmd_positions_angstrom(record)
```

The dictionary/JSON schema is `{"columns": [...], "frames": [{"cell":
[["a_x","a_y","a_z"], ...], "atoms": [["C","x","y","z", ...], ...]}]}`.
The example's atom order assumes the corresponding column order; arbitrary
column arrangements remain supported. Unknown fields, numeric JSON values
instead of strings, duplicate JSON keys (including escaped aliases), and wrong
array shapes refuse. Python input is `str` or UTF-8 `bytes`; invalid UTF-8 and
wrong argument types use normal Python decode/type errors.

Rust exports `LpmdDocument`, `LpmdFrame`, `LpmdParseLimits`, `LpmdError`,
`parse_lpmd`, `parse_lpmd_with_limits`, `parse_lpmd_json`, `write_lpmd` and
`write_lpmd_with_limits`. Document methods validate and explicitly project
positions. WASM exports `parse_lpmd_json`, `write_lpmd_from_json` and
`lpmd_positions_angstrom_json`; JSON projection is grouped by frame, atom and
Cartesian component. Python returns the corresponding nested lists.

```sh
chematic lpmd -i input.lpmd -o output.lpmd
chematic lpmd --inspect -i input.lpmd > record.json
chematic lpmd --from-json -i record.json -o edited.lpmd
chematic lpmd --angstrom -i input.lpmd > positions.json
```

`--from-json --angstrom` is allowed. `--inspect` conflicts with both mode flags.
Projection output intentionally contains positions only; keep the record to
retain cell and optional columns. A format refusal is `LpmdError`, Python
`ValueError`, a thrown WASM error string, or CLI exit 2 with empty stdout.
Format errors begin `LPMD line `; line zero means a document/edited-record error.

Native/Python defaults are 64 MiB text/JSON input and text output, 8,192-byte
physical lines, 1,024 frames and **1,000,000 total atom rows** across all frames.
Python exposes keyword overrides for these four limits. WASM limits input and
text output to 1,000,000 bytes and total rows to 10,000; its JSON output cap is
16 MiB. CLI uses native defaults and a 64 MiB final output cap, with
`--max-input-bytes` controlling input and text output. This is a materialized
trajectory API, not an unbounded streaming reader.

## Local comparison and verification

Original synthetic MIT/Apache-2.0 inputs are in
`validation/lpmd_preservation/cases.json`. The retained observation is
`validation/lpmd_preservation/observed.json`; runtime hashes, implementation
source hashes, base commit and each Open Babel subprocess receipt are recorded.
This identifies a local source candidate, not a clean release artifact.

The packet has **26 positive records**, **53 malformed/unsupported inputs** and
**two raw-valid records whose explicit projection must refuse**. Python, real
Node/WASM and the CLI preserve every field of all 26 records through initial
read plus three write/read cycles. Every projection agrees with the independent
row-vector observer as binary64 bytes. All 53 refusals and both unsafe
projection cases pass on all three surfaces. Coverage includes all 118 element
symbols, reordered columns, every optional column, changing cells and atom counts
across frames, left-handed/rotated bases, raw precision and subnormals. Rust
integration tests exercise the same packet; observer controls test dropped
columns/frames, malformed output, physical equivalence and wrong geometry.

Installed **Open Babel 3.2.1**, with bond perception disabled (`-ab`), can expose
source-to-XYZ geometry for **22/26** inputs. Frame/atom counts and elements match
in those 22; Cartesian binary64 bytes match in **20/22**. XYZ rounding changes
the high-precision and subnormal cases. This route does not observe cells or
velocity/force/color preservation, and is not a full-record preservation score.

Its default LPMD self-output has a specific framing defect: **22 nonempty
outputs declare `HDR X Y Z` but write a leading element symbol in each row**.
An independent inspector rejects the missing `SYM` declaration. Reading those
outputs back with Open Babel produces wrong elements/positions in all 22
observable XYZ outputs. Four other inputs (tabs, trailing blanks and two
multi-frame cases) produce no output with error diagnostics; they are recorded
as unobservable, not safe refusals or silent-loss wins. Exit status is zero in
these runs, so output inspection is necessary. No output options were enabled
for this default-lane observation; alternative options or versions are not
ranked. The original files and the re-read generated files are separate lanes.

These results establish a bounded local preservation improvement for this
adapter. They do not establish toolkit-wide superiority, full LPMD engine
interoperability, performance, a published package or hosted CI completion.

Reproduce with freshly built runtimes:

```sh
python3 scripts/validation/check_lpmd_preservation.py \
  --python-path /path/to/python-package-directory \
  --wasm-module /path/to/node/chematic_wasm.js \
  --cli /path/to/chematic --openbabel /path/to/obabel \
  --output /tmp/lpmd-observed.json
python3 -m unittest scripts.tests.test_lpmd_observer
```

## Primary sources

The format layout follows the authors' [LPMD paper, Figure 10 and section 5](https://repositorio.uchile.cl/bitstream/handle/2250/119202/Davis_Sergio.pdf?sequence=1)
and the [official LPMD format plugin description](https://github.com/lpmd/plugins/blob/master/lpmd.cc).
The [Open Babel adapter documentation](https://openbabel.org/docs/FileFormats/LPMD_format.html)
separates output level, acceleration/force mode, cell and charge options.
External GPL implementation files were consulted as references and were not
copied or translated into the implementation. The official LPMD executable
was not built or run; its historical format-site URL is unavailable.
