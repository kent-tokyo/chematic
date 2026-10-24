# Extended XYZ output integrity

`write_extxyz`, `ExtxyzWriter::write_frame`, Python `to_extxyz` and WASM
`to_extxyz_json` now validate the complete frame before writing it. Previously,
hand-built native properties could panic or emit inconsistent columns;
nonfinite geometry/data could be emitted by Python; some literal metadata was
changed when read by ASE. No new public export or file-format alias is added.

```python
import chematic

mol = chematic.from_smiles("CO")
text = chematic.to_extxyz(
    mol, [[0.125, -0.0, 0.25], [1.5, 0.75, -0.25]],
    properties={"vectors": [[0.1, 0.2, 0.3], [0.4, 0.5, 0.6]]},
    info={"path": "alpha\\beta", "note": "O'Brien"},
)
assert chematic.from_extxyz(text)["info"]["path"] == "alpha\\beta"
```

The writer still uses source atom order and round-tripping finite binary64
spelling for XYZ, cell and real atom data. Metadata containing backslashes,
apostrophes, brackets/braces or `=` is now quoted and escaped, so those literal
characters survive the measured independent ASE consumer. Metadata quote style
can change. Canonical reserialization does not preserve source bytes.

## Preflight contract

- All coordinates, lattice entries and `R` property components must be finite.
- Every property has exactly one row per atom, exactly its declared positive
  component count, and values matching its `R`/`I`/`S`/`L` declaration. Aggregate
  column-count overflow is refused. Native string tokens must be nonempty and
  contain no whitespace, controls, quotes, backslashes, `=`, brackets or braces.
- Info keys and property names must be unique and writable without quoting.
  Empty names, whitespace, controls and grammar delimiters are refused;
  property names additionally exclude `:`. `Lattice` and `Properties` info keys
  remain reserved. Newline/control injection in info values is refused.
- Extra atom properties cannot be named `species`, `pos`, `symbols`, `positions`,
  `numbers` or `Z`: those names collide with element/position fields in the local
  or ASE projection. `charge` and `charges` cannot coexist because ASE maps both
  to the same atom array. A single such charge column remains allowed.
- With structured metadata, the header is rebuilt from lattice/properties/info;
  `frame.comment` is the original lexical header and is superseded, as before.
  With none of those fields present, the free comment is retained but must not
  inject lines or declare metadata/columns missing from the structured frame.
  A species/pos-only `Properties=` header remains writable.
- Bare empty info strings are refused. In ASE 3.26.0, `empty="" note="text"`
  can be read as one combined `empty` value. If an empty string is intentional,
  explicitly pass the raw typed tag `'_JSON ""'`. ASE decodes it to an empty
  string; the native reader retains that raw tag, matching its existing raw-info
  contract. Numeric/boolean/PBC/JSON-tagged info may be typed by ASE; native info
  remains a list/dict of raw string values. Arbitrary reader-wide string typing
  equivalence is not promised.

Native streaming validation returns `io::ErrorKind::InvalidData` before the count
line, leaving an existing sink prefix untouched. `write_extxyz` returns the
existing `XyzError::UnwritableMetadata` variant for invalid/unrepresentable frame
data, including geometry and property errors. Sink I/O errors still propagate
and can leave part of a valid frame written. This is frame preflight, not an
atomic-file or whole-trajectory transaction.

Python/WASM still expose real-valued atom properties only. Rust also supports
integer, logical and the checked single-token string profile. WASM options must
be an object containing only `lattice`, `properties` and `info`; malformed types,
unknown fields and duplicate decoded keys now fail instead of being ignored.
`null` optional fields still mean absence. Existing public signatures are unchanged.

These are deliberate producer-side refusals of previously accepted unsafe data.
Readers are unchanged and can accept a wider profile than the writer. Molecular
bonds, formal charges, isotope identity and stereo are not transferred by these
geometry writers. Supply an atom-data column explicitly when it is needed;
plain XYZ output remains a separate API. Output materialization is still governed
by the existing APIs; no new configurable output-size limit is introduced.

## Independent consumer evidence

The original synthetic packet `validation/extxyz_output_integrity/cases.json`
contains nine positive frames and 19 refusal inputs, without duplicated newline
variants or an unseen-holdout claim. It includes a plain-geometry control,
binary64 precision/signed zero, an asymmetric cell, vector/scalar atom data,
backslash/apostrophe/bracket metadata, Unicode/quotes, typed external metadata
and explicitly tagged empty strings. Native tests separately check integer,
logical and string columns, numeric extremes, shape/type failures, duplicate
names, stale headers, sink failures and refusal before any bytes are written.

The binding writers receive structured fixtures. Open Babel receives independently
authored XYZ text containing the same `Lattice`/`Properties`/info data. ASE 3.26.0
reads that authored source and each output; source geometry, cell convention,
property bits/shapes and expected metadata are also asserted against the fixture.
An independent raw four-column observer measures geometry even when the emitted
header/row combination is not readable as ASE ExtXYZ.

| Path | Geometry bits | ASE-readable | Complete ASE consumer state |
|---|---:|---:|---:|
| Candidate Python / WASM, each | 9/9 | 9/9 | 9/9 |
| Earlier local Python / WASM, each | 9/9 | 9/9 | 6/9 |
| Open Babel 3.2.1 XYZ default / `-ab`, each | 8/9 | 8/9 | 7/9 |

The earlier local writer changes three metadata strings: backslash, apostrophe
and bracket cases. The candidate refuses all 19 outside-profile inputs, compared
with seven explicit refusals in the earlier local bindings. Refusals are not
counted as successful preservation. Separate Python tests cover NaN and both
infinities; separate real-WASM tests cover malformed/duplicate JSON options.

Open Babel rounds one precision case and drops the atom-data columns from the
vector frame while retaining its `Properties=` comment. That output is not
readable by ASE as the original frame. The other seven controls preserve the
measured consumer state. This is a **workflow comparison through Open Babel's
plain XYZ route**, whose [documented scope](https://openbabel.org/docs/FileFormats/XYZ_cartesian_coordinates_format.html)
is element/XYZ rows with five decimal places. It is not evidence that Open Babel
claims support for ASE ExtXYZ, and its separate `exyz` format is not an alias for
this schema. No toolkit-wide, chemical-graph, speed or general ExtXYZ coverage
ranking follows from these nine frames.

Records retain complete fixtures, authored text, all outputs/errors, command
arguments, a 10-second/1-MiB combined-pipe process bound, source/runtime hashes,
Open Babel version and the installed ASE reader hash/version. The baseline is a
preserved local `fda4dfac` build from PR #1187, not a release or main build. Its
ExtXYZ writer/binding source was verified identical to main `255edf73`. Candidate
records identify working-tree source hashes on that parent, before commit.

- `validation/results/extxyz-output-integrity-baseline-2026-10-24.json.gz`
- `validation/results/extxyz-output-integrity-candidate-2026-10-24.json.gz`

Reproduce using ASE 3.26.0, a freshly built Python extension on `PYTHONPATH` and a
fresh Node-target WASM module:

```sh
python scripts/validation/check_extxyz_output_integrity.py \
  --wasm-module /path/to/chematic_wasm.js \
  --output /tmp/extxyz-output.json.gz
python -m unittest scripts.tests.test_extxyz_output_integrity
```

ASE's [pinned reader source](https://gitlab.com/ase/ase/-/blob/3.26.0/ase/io/extxyz.py)
provides the independent metadata and atom-array interpretation. It is a test
dependency; its implementation is not copied into chematic.
