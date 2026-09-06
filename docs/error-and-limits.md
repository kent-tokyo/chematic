# Error Handling & Parse Limits

What actually limits input size, how typed errors surface in each language,
and which operations are explicitly named as lossy. This describes real,
specific, checkable behavior — bounded input validation on the formats
that have it, not a claim that chematic is broadly "secure." Formats
without a documented limit type are named as such below, not silently
folded into a "handled" bucket.

See also: [`format-capabilities.md`](format-capabilities.md) for the
per-format matrix these limits are drawn from, and
[`language-bindings.md`](language-bindings.md) for the `ValueError`/
JS-error mapping in more general terms.

---

## Parse limits, by format

The formats below have dedicated `*ParseLimits` or reader-options types. The
rest do not — this is not
inconsistent oversight to be quietly worked around; it reflects that most of
these formats have no natural analogue of "grid points" or a comparably
unbounded substructure, or simply have not had a limits type added yet.

| Format | Limits type | Fields |
|---|---|---|
| mmCIF | `MmcifParseLimits` | `max_input_bytes`, `max_atoms`, `max_line_len` |
| PQR | `PqrParseLimits` | `max_input_bytes`, `max_atoms`, `max_line_len` |
| Gaussian Cube | `CubeParseLimits` | `max_input_bytes`, `max_atoms`, `max_grid_points` |
| OpenDX | `OpenDxParseLimits` | `max_input_bytes`, `max_grid_points` (no `max_atoms` — the format has no atom section) |
| SMILES | none | — |
| SMARTS | none | — |
| MOL/SDF | none | — |
| PDB | none | — |
| CIF (plain) | none | — |
| XYZ / Extended XYZ | `chematic_3d::XyzParseLimits`, `chematic_mol::XyzParseLimits` | 3D XYZ: `max_input_bytes`, `max_atoms`, `max_line_bytes`; extended XYZ: input/atom/frame/line/property limits |
| KET | `KetParseLimits` | `max_input_bytes`, `max_atoms`, `max_bonds` |
| SMILES table | `SmiFileParseLimits` | `max_input_bytes`, `max_line_bytes`, `max_records` |
| Streaming SMILES table | `SmilesReaderOptions` | `max_line_bytes`, `max_records`, `max_fields` |
| Streaming TDT | `TdtReaderOptions` | `max_line_bytes`, `max_records`, `max_tags_per_record` |
| InChI (pure Rust) | `InchiParseLimits` | `max_input_bytes`, `max_atoms` |
| MDL RXN | `RxnFileParseLimits` | `max_input_bytes`, `max_reactants`, `max_products`, `max_molecules` |
| QCSchema | none | — (JSON-size limits, if any, are whatever the caller or `serde_json` impose) |
| ORCA input | none | — |
| ORCA output | none | — |
| LAMMPS data | none | — |
| LAMMPS dump/trajectory | none | — |

In Python, each `*ParseLimits` field is exposed as an optional keyword
argument on the corresponding `parse_*_with_limits`-style function, using
the Rust `Default` values when omitted (e.g. `parse_mmcif(text,
max_input_bytes=None, max_atoms=None, max_line_len=None)`).

WASM convenience renderers also apply binding-level bounds before building
intermediate molecule vectors: `depict_svg_grid`,
`depict_svg_grid_highlighted`, and `batch_report_html` accept at most 1 MiB of
input, 1,024 non-empty records, and 10,000 atoms per parsed molecule. Since
these legacy APIs return strings rather than `Result`, a violation returns an
explicit SVG/HTML error document and never a partial rendering.

The high-level WASM workflow APIs apply the same 1 MiB and 1,024-record
boundaries, plus a 10,000-atom per-molecule bound, before comparison,
screening, reporting, or 3D generation. APIs returning `Result` surface the
violation as a JavaScript error; `screen_smiles_json` returns an explicit JSON
error object.

WASM JSON writers for mmCIF, PQR, and LAMMPS trajectories also cap every
record/frame array at 1,024 items before converting values into Rust vectors.
The same cap applies to nested record arrays such as ORCA coordinate atoms.
LAMMPS trajectory readers apply the same cap while retaining parsed frames;
the 1,025th frame is rejected rather than returned as a partial success.
WASM ORCA output parsing likewise applies the binding input, line, geometry
frame, and geometry atom limits through `OrcaOutputParseLimits` before
serializing the trajectory.
ORCA input and QCSchema WASM parsers also override the larger Rust defaults
with binding-sized line/block/atom and JSON depth/array/string limits before
decoding structured documents.
Cube/OpenDX, mmCIF, PQR, LAMMPS data, and single-frame dump WASM parsers
likewise override their larger Rust defaults for input bytes, atom/record
counts, line size, section counts, columns, and volumetric grid points.
Cube/OpenDX typed-array accessors use those same limits rather than delegating
to the larger parser defaults, so JSON and typed-array entry points reject the
same oversized grids.
The WASM `convert_common_format` bridge also rejects parsed molecules over the
10,000-atom binding limit before serializing a target format.
The MCP PubChem lookup caps the upstream response body at 1 MiB and rejects
oversized or invalid UTF-8 responses before JSON parsing.
MCP SMILES tool arguments are limited to 100,000 bytes and 10,000 parsed
atoms; SMARTS matching applies the 100,000-byte query limit before compiling
the query.
MCP MolJSON input uses explicit limits of 100,000 bytes, depth 64, 1,024 array
items, 10,000 atoms, and 20,000 bonds before molecule construction.
MCP SMARTS matching is bounded to 10,000 embeddings and 1,000,000 VF2 visits;
exhaustion is reported as an explicit domain error rather than a partial list.
MCP stdio JSON-RPC responses are capped at 1 MiB after envelope serialization;
an oversized response is replaced by an explicit protocol error.
MCP stdio request framing also reads at most the 1 MiB request boundary before
JSON parsing; an oversized frame terminates the loop rather than being split
into a later request.
CLI output is capped at 64 MiB for both stdout and file destinations; an
oversized result is rejected before any bytes are written.
Format-oriented WASM JSON helpers cap serialized results at 16 MiB and return
an explicit JS error when the boundary is exceeded.
The Extended XYZ frame JSON helper uses the same shared output boundary.
WASM tautomer enumeration also caps the result array at 1,024 entries and
returns an explicit error object when the cap is exceeded.
WASM `run_reactants` bounds reaction input to 1 MiB, reactants to 1,024
entries, each reactant to 10,000 atoms, and each product set to 1,024 entries.
WASM MCS configuration JSON is also limited to the 1 MiB binding input
boundary before deserialization.
WASM R-group decomposition applies the same 1 MiB input limit to its
`core_smarts` query before compilation.
WASM MMP pair output is capped at 1,024 pairs before JSON materialization;
larger result sets return an explicit JS error.
WASM R-group queries are limited to 10,000 SMARTS query atoms after
compilation, in addition to the 1 MiB `core_smarts` source limit.
WASM reaction and library products are also limited to 10,000 atoms before
canonicalization and JSON materialization.
WASM library templates and fragment inputs are limited to 1 MiB, and each
scaffold/building-block molecule is limited to 10,000 atoms before enumeration.

---

## Checked arithmetic and NaN/Infinity rejection

- `VolumetricGrid::checked_index` (shared by Cube and OpenDX) computes the
  flat-array index from `(i, j, k)` with explicit bounds checking — it
  returns `None` rather than panicking or wrapping on out-of-range input,
  including index-arithmetic overflow.
- QCSchema rejects every non-finite numeric leaf (`NaN`/`Infinity`) rather
  than letting it pass through — `serde_json::Number` itself cannot
  represent a non-finite value, and the parser treats an attempt to smuggle
  one in as a typed error rather than silently coercing it to `null` or `0`.
- LAMMPS dump's box-bounds/triclinic conversion (`box_bounds_to_true`/
  `true_to_box_bounds`) is typed to reject a non-finite coordinate, not just
  a malformed one.

---

## Unsupported-format-subset behavior: reject, don't guess

Every one of these is a **typed rejection**, not a silent best-effort
fallback or a silent truncation:

| Situation | Behavior |
|---|---|
| LAMMPS data with an `atom_style` outside `atomic`/`charge`/`molecular`/`full` | `LammpsDataError::UnsupportedAtomStyle` — `charge` and `molecular` rows are both genuinely ambiguous 6-field rows by column count alone, so there is no safe guess to fall back to. |
| LAMMPS "Type Labels" sections (`Atom Type Labels`, `Bond Type Labels`, etc.) | rejected — any section name ending in `"Type Labels"` fails the parse rather than being silently misread as a different section shape. |
| Gaussian Cube file with more than one dataset | `CubeError::MultiDatasetUnsupported` — rejected, not silently truncated to the first dataset. |
| OpenDX write of a `Bohr`-tagged grid via `write_opendx` | `OpenDxError::NonAngstromUnits` — see the fail-closed section below. |
| OpenDX write of a grid with any atoms (either writer) | `OpenDxError::AtomsNotSupported` — there is no lossy-atom-dropping path; the format simply has no atom section to write to. |

---

## Typed-error taxonomy → `ValueError` / JS-error mapping

Every format module defines its own error enum (`LammpsDataError`,
`OpenDxError`, `CubeError`, `MmcifError`, `PqrError`, `OrcaInputError`,
`OrcaOutputError`, `QcSchemaError`, `CifError`, `XyzError`,
`MolParseError`, `CdxmlError`, `CmlError`, ...). The mapping across
language boundaries is uniform but lossy in one specific sense:

| Language | Representation |
|---|---|
| Rust | the original typed enum variant, with its own fields (e.g. `CubeError::MultiDatasetUnsupported { natoms_field, nval }`) |
| Python | `ValueError`, message = the Rust error's `Display` text |
| WASM | thrown JS error / `JsValue::from_str(...)`, message = the Rust error's `Display` text |

Both bindings currently flatten the structured Rust variant down to a
string message — neither language exposes the original enum's fields
programmatically across the boundary today. If you need to distinguish
error *kinds* (not just read a message) in Python or JS, you currently have
to do it by matching on message text, which is a real limitation, not an
oversight this page is hiding.

---

## Fail-closed writers

`write_opendx` is the clearest example in the codebase: given a grid whose
`units` field is `GridUnits::Bohr`, it refuses to write rather than
producing an OpenDX file (a format with no in-file unit tag) that a
downstream reader would silently misinterpret as Ångström. The explicit
opt-in for the lossy path is `write_opendx_lossy` — see the naming
convention below. This same fail-closed principle applies to the
unsupported-format-subset rejections in the table above: every one of them
raises a typed error instead of writing/returning a best-effort,
potentially-wrong result.

---

## The `*_lossy` naming convention

Every explicitly lossy write path in this codebase is named with a
`_lossy` suffix and is **opt-in only** — never the default, never triggered
implicitly by a missing flag:

- `write_opendx_lossy` (Rust) / `VolumetricGrid.to_opendx_lossy()` (Python)
  / `write_opendx_lossy_json` (WASM) — the sole named lossy operation
  across all 15 formats in this pass. It rescales `origin`/`axes` from
  Bohr to Ångström; it never rescales `values`, and it still refuses
  (`AtomsNotSupported`) a grid with atoms.

No other format in this documentation pass has a `_lossy`-suffixed
function. Where a format has a *documented* lossy characteristic that isn't
behind an explicit opt-in — e.g. mmCIF's `label_*`/`auth_*` tag-pair
collapse on read, which is unconditional, not a named lossy function — it
is called out in [`format-capabilities.md`](format-capabilities.md)'s
per-format detail instead, precisely because it doesn't fit this
opt-in-only convention.

---

## Formats that never fabricate bonds

mmCIF, PQR, ORCA (input and output), Gaussian Cube, OpenDX, LAMMPS data,
and plain CIF never infer or fabricate a bond table. PDB and
`chematic_3d::parse_xyz` are the two exceptions among the 15 covered here —
both infer bonds from 3D geometry (distance-based), a disclosed,
documented choice specific to those two entry points, not a default
behavior you should assume elsewhere. See
[`format-capabilities.md`](format-capabilities.md#formats-that-never-fabricate-bonds)
for the full per-format connectivity notes.

---

## Streaming vs. full materialization, per format/language

LAMMPS dump/trajectory is the one format with a real Rust-level
streaming-vs-materializing distinction: `LammpsDumpReader<R: BufRead>` is a
true streaming `Iterator`, while both the Python (`parse_lammps_dump_all`)
and WASM (`lammps_trajectory_to_json`) bindings materialize the whole
trajectory instead — a documented scope choice. MOL/SDF's `SdfFileReader<R:
BufRead>` is a true streaming reader too, but this is not currently called
called out separately. Python's `iter_sdf` and `iter_sdf_batched` now expose
the file-backed reader; the batch iterator has bounded batch sizes, explicit
cancellation, and a JSON progress manifest. See
[`language-bindings.md`](language-bindings.md#streaming-vs-materialization-by-language)
for the full per-format table.

---

## Misuse-prevention examples

These are short illustrations of the fail-closed/typed-rejection behavior
above — not new functionality, just what already happens if you try the
tempting-but-wrong thing.

**OpenDX: writing a Bohr-tagged grid without `_lossy`**

```rust
// grid.units == GridUnits::Bohr
let err = write_opendx(&grid).unwrap_err();
// err is OpenDxError::NonAngstromUnits { units: GridUnits::Bohr } —
// the file is not written. Use write_opendx_lossy(&grid) if you
// intend the Bohr->Angstrom conversion (rescales origin/axes only).
```

**Gaussian Cube: a multi-dataset file**

```rust
// input .cube has natoms field encoding "2 datasets"
let err = parse_cube(text).unwrap_err();
// err is CubeError::MultiDatasetUnsupported { natoms_field, nval } —
// not silently parsed as dataset 1 of N.
```

**LAMMPS: guessing `atom_style` instead of stating it**

```rust
// A 6-field Atoms row is ambiguous: could be `charge` (id type q x y z)
// or `molecular` (id mol-id type x y z). Passing the wrong style
// doesn't error at read time -- it silently mis-assigns which
// column is charge vs. molecule-id. Always pass the atom_style the
// simulation actually used; an out-of-set value (e.g. a typo) is
// caught (`UnsupportedAtomStyle`), but a *wrong-but-valid* style is not.
let data = parse_lammps_data(text, LammpsAtomStyle::Full)?; // be explicit
```

**QCSchema: assuming `connectivity` is always present**

```rust
let qc: QcMolecule = parse_qcschema_molecule(text)?;
// qc.connectivity: Option<Vec<(usize, usize, f64)>> -- QCSchema does
// not require a bond list. Unwrapping without checking panics on any
// spec-valid molecule that omits it (most do).
if let Some(bonds) = &qc.connectivity {
    // only reachable if the source document actually included one
} else {
    // no bonds available -- do not fabricate them here
}
```
