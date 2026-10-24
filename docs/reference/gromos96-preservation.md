# GROMOS96 block and unit preservation

Rust, Python, WASM/Node and CLI expose bounded GROMOS96 read/write operations.
This adapter was absent from base `a73204f9`; Open Babel's `gr96` adapter is
[write-only](https://openbabel.org/docs/FileFormats/GROMOS96_format.html).
The new profile retains every supported block, frame, label, comment, numeric
token and fixed-field padding as current editable data. It never guesses atom
elements, chemical topology or length units from labels or filenames.

## Supported profile

A document has one initial `TITLE` block and one or more coordinate frames.
The supported frame order is optional `TIMESTEP`, required `POSITION` or
`POSITIONRED`, optional corresponding `VELOCITY` or `VELOCITYRED`, then optional
`BOX`. Each block ends with `END`. New coordinate/time blocks start subsequent
frames. Frame atom counts may vary; no missing time, cell or velocity values
are propagated from a previous frame.

The authoritative model is an ordered list of blocks with `kind`, `leading`
layout lines and current payload `lines`, plus final `trailing` layout lines.
Column-one `#` comments remain in their original locations. Title text may be
empty, multiline or UTF-8. Coordinates and identity prefixes are ASCII. CRLF
and a missing final newline are accepted; output canonicalizes block markers
and LF termination while retaining all payload/layout line contents. This is
field/payload preservation, not a promise of byte-identical source replay.

Full coordinate and velocity rows use the GROMACS-compatible 24-character
identity prefix: residue number `[0:5]`, separator, residue name `[6:11]`,
separator, atom name `[12:17]`, atom number `[17:24]`. The three numeric fields
are each 15 characters, for a 69-character row. Reduced rows contain only the
three fields, for 45 characters. A checked variant with an extra space before
each numeric field (72/48 characters), including Open Babel output, is also
accepted and retained. Arbitrary free-format rows are outside this profile.

Names are nonempty contiguous ASCII labels up to five characters, **not element
symbols**. Signed residue numbers and zero atom numbers are accepted. Velocities
must have the same full/reduced style, row count and trimmed identity fields as
the corresponding positions; mismatches refuse. Numeric identity spelling is
retained, so differing spellings such as `1` and `00001` are not reconciled.

`TIMESTEP` has one row: a 15-character signed integer and a 15-character decimal.
`BOX` has one row of three or nine 15-character decimals. Its order is
`ax, by, cz[, ay, az, bx, bz, cx, cy]`; projection reconstructs row vectors
without reduction, rotation or wrapping. Zero/singular BOX values remain source
data; this API does not assert that they define a usable periodic simulation.
Numbers support signed decimal and `e/E` notation within their field widths.
Nonfinite nominal binary64 values, nonzero decimals underflowing to zero,
Fortran `D` notation, incomplete fields and overflowed fields refuse.

Unsupported blocks (including `GENBOX` and `POSITIONREF`), compressed input,
repeated titles, out-of-order or duplicate velocity/BOX blocks, missing `END`,
orphan timesteps and trailing non-layout data refuse. Blank lines inside numeric
blocks refuse. Opaque comments are retained but not interpreted as simulation
parameters. This is a GROMACS-compatible block subset plus the stated spaced
variant, not the complete GROMOS ecosystem or arbitrary trajectory syntax.

## Units stay explicit

GROMOS96 text has no mandatory length-unit declaration. Open Babel's default
writer emits angstrom-valued coordinates, while `-xn` emits nanometers. The
same unmarked numbers therefore cannot safely establish their physical unit.

`length_unit` is caller-owned context: `"nm"`, `"angstrom"` or `"unspecified"`.
Reading without a declaration preserves the raw document with `unspecified`.
Raw writing does not rescale data or insert a proprietary unit marker. On a
subsequent text read, supply the same authoritative unit context again; keep it
in the JSON record or an application sidecar. **Text alone cannot round-trip
this external context.** Changing the record's unit field reinterprets the
numbers; it is not a physical conversion operation.

The separate `geometry_angstrom` projection requires a declared unit and
converts both positions and BOX vectors by the same factor (10 for nm, 1 for
angstrom). It refuses nonfinite products and underflow to zero. It returns all
frames, with original label identities or `None` for reduced positions. It does
not convert velocities/time units or return molecular chemistry. The source
record remains available, including velocity and time payloads. Unsafe raw
records can still be read/written even when projection refuses.

## APIs

```python
from pathlib import Path
import chematic

# Use nm only when the producing workflow establishes that convention.
record = chematic.parse_gromos96(Path("input.g96").read_bytes(), length_unit="nm")
geometry = chematic.gromos96_geometry_angstrom(record)
Path("output.g96").write_text(chematic.write_gromos96(record), encoding="utf-8")
# Retain record["length_unit"] alongside the unmarked output file.
```

The JSON/dictionary schema is `{"length_unit": "nm", "blocks": [{"kind":
"TITLE", "leading": [], "lines": ["title"]}, ...], "trailing": []}`. Writer
edits are made to current block payloads and revalidated; stale source is never
replayed. For example, a full atom row's first coordinate is `[24:39]` in the
69-character style. Numeric strings that do not fit a field must be rejected
or explicitly reformatted by the caller; this writer never rounds/truncates
them to make them fit. Unknown JSON fields, duplicate keys (including escaped
aliases), non-string payload entries and invalid shapes refuse.

Rust exposes `Gromos96Document`, `Gromos96Block`, `Gromos96LengthUnit`,
`Gromos96Identity`, `Gromos96Geometry`, `Gromos96ParseLimits`, `Gromos96Error`,
`parse_gromos96`, `parse_gromos96_with_limits`, `parse_gromos96_json`,
`write_gromos96` and `write_gromos96_with_limits`. Document methods validate and
project geometry. Python accepts text or UTF-8 bytes and returns dictionaries
and nested lists. WASM exports `parse_gromos96_json(text, length_unit)`,
`write_gromos96_from_json` and `gromos96_geometry_angstrom_json`.

```sh
chematic gromos96 -i input.g96 -o output.g96
chematic gromos96 --inspect --length-unit nm -i input.g96 > record.json
chematic gromos96 --from-json -i record.json -o edited.g96
chematic gromos96 --angstrom --length-unit nm -i input.g96 > geometry.json
```

`--from-json --angstrom` uses the JSON record's unit. `--length-unit` conflicts
with `--from-json`, and `--inspect` conflicts with both output-mode flags.
An unspecified unit refuses projection before file output. Format errors use
`GROMOS96 line ` diagnostics, Python `ValueError`, thrown WASM strings or CLI
exit 2 with empty stdout. Line zero denotes a whole-document/payload/edited-record
check; only framing/line-limit diagnostics promise a physical source line.
Wrong Python argument types and invalid UTF-8 use normal type/decode errors.

Native/Python defaults: 64 MiB input/JSON/text output, 8,192-byte physical lines,
1,024 frames and **1,000,000 position rows in total**, with optional velocity
rows checked against those positions. Python exposes overrides for all four
limits. WASM caps input/text output at 1,000,000 bytes and total position rows
at 10,000; JSON output has a 16 MiB cap. CLI uses native defaults and a 64 MiB
final output cap, with `--max-input-bytes` controlling input/text output. The
document API materializes supported frames; it is not an unbounded stream.

## Local evidence and limits

Original synthetic MIT/Apache-2.0 fixtures and receipts are retained in
`validation/gromos96_preservation/cases.json` and `observed.json`. Implementation
and runtime hashes identify a local source candidate, not a clean package
release. All three actual Python, Node/WASM and CLI runtimes preserve all fields
of **26 positive records** through initial read and three write/read cycles,
with caller context resupplied on every text read. Projected identities,
positions and BOX values match the independent nominal-binary64 observer.
**45 malformed/unsupported inputs** refuse; **three raw-valid records** preserve
their data while refusing projection (unknown unit, position overflow or BOX
overflow). Rust tests exercise the same packet. Observer controls cover unit
mistakes, BOX scaling, dropped frames/comments, unobservable output, label
ambiguity, signed zero and physical equivalence across unit conventions.

**ASE 3.26.0** independently reads positions/BOX from **18/26** original files
and 18/26 regenerated files, with count agreement and geometry error within
`2e-14 * max(1, max(abs(expected geometry components)))`. Reduced rows, some
comments/dense fields, unresolved labels and the multi-frame workflow are
unobservable in this reader lane. They are not classified as data loss or safe
refusal. ASE atom-type guesses are not scored, and this lane does not observe
title, time, velocity, label preservation or a frame sequence.

A separate **one-way writer** lane gives both writers six authored, common
element/Cartesian-coordinate inputs in both angstrom and nm modes: **12 cases**.
The native writer receives validated fixed-field records constructed without
rounding; Open Babel 3.2.1 receives the corresponding XYZ data and writes `gr96`,
with `-xn` only for the nm mode. Independent decimal arithmetic applies the
known unit convention and compares the resulting physical coordinates to the
original decimal inputs: chematic **12/12**, Open Babel **8/12**. Open Babel's
five-decimal output changes the precision and tiny-coordinate cases in both
modes. This compares representable coordinate preservation, not equal parsing
work, molecular metadata, round-trip behavior or speed.

ASE exposes **10/12** native writer outputs and **12/12** Open Babel outputs,
matching each file's declared numerical geometry within the above tolerance.
The two dense high-precision native field cases are unobservable to this
whitespace-based ASE reader. Matching rounded output does not establish source
precision preservation. No GROMACS/GROMOS executable was built or run, and
complete engine compatibility, all external readers, toolkit-wide superiority
and performance are not established.

Reproduce with fresh runtimes and ASE 3.26.0 installed in a separate directory:

```sh
python3 scripts/validation/check_gromos96_preservation.py \
  --python-path /path/to/python-package-directory \
  --wasm-module /path/to/node/chematic_wasm.js --cli /path/to/chematic \
  --ase-path /path/to/ase-directory --openbabel /path/to/obabel \
  --output /tmp/gromos96-observed.json
python3 -m unittest scripts.tests.test_gromos96_observer
```

The profile follows the official [GROMACS G96 format documentation](https://manual.gromacs.org/current/reference-manual/file-formats.html#g96)
and fixed-field declarations in the [GROMACS G96 I/O source](https://github.com/gromacs/gromacs/blob/main/src/gromacs/fileio/g96io.cpp).
The [Open Babel adapter documentation](https://openbabel.org/docs/FileFormats/GROMOS96_format.html)
defines its direction and nm output option. External implementation sources
were consulted as format/behavior references and were not copied or translated
into this implementation.
