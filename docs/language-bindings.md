# Language Bindings Cross-Reference

chematic ships three language surfaces over the same Rust core: the Rust
crates directly, Python (PyO3), and WASM (wasm-bindgen). They are not
identical APIs — this page documents where they agree, where they diverge,
and why, using real function names (not illustrative examples) so the
divergences are checkable against source.

See also: [`format-capabilities.md`](format-capabilities.md) for the
per-format read/write/streaming/limits matrix this page's examples are
drawn from.

## RDKit-compatible Python import surface

The Python package includes an explicit, lightweight compatibility namespace at
`chematic.rdkit_compat`. The common Morgan entry points below resolve to the
same chematic-backed fingerprint implementation, so migration does not depend
on one particular RDKit import style:

```python
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import AllChem, rdFingerprintGenerator

mol = Chem.MolFromSmiles("CCO")
legacy = AllChem.GetMorganFingerprintAsBitVect(mol, 2)
stable = rdFingerprintGenerator.GetMorganGenerator(
    radius=2, fpSize=2048
).GetFingerprint(mol)
assert legacy.GetOnBits() == stable.GetOnBits()
```

`Chem.CanonSmiles` is also available. Unsupported RDKit options such as
count-simulation Morgan fingerprints raise `NotImplementedError`; they are not
silently substituted with a different algorithm. This is a Python API
compatibility surface, not a claim of bit parity with RDKit's C++ implementation.

---

## The three surfaces, in one sentence each

- **Rust** (`chematic-*` crates): the source of truth. Typed errors
  (`thiserror`-style enums), zero-copy where the type system allows it,
  `BufRead`-based streaming where a reader type exists.
- **Python** (`chematic-py`, PyO3): free functions mirroring the Rust free
  functions where practical (`parse_mmcif`, `parse_pqr`, ...), plus `Mol`
  and a handful of dedicated pyclasses (`VolumetricGrid`, `LammpsDumpFrame`,
  `PipelineV2Config`). Typed Rust errors map to Python `ValueError`.
  NumPy arrays are used for large numeric payloads (fingerprints, grid
  values) — **always fresh copies, never views**.
- **WASM** (`chematic-wasm`, wasm-bindgen): mostly free functions returning
  JSON strings (`Result<String, JsValue>`), with selected functions also
  returning `js_sys::Float64Array`/`Uint32Array` for
  large numeric grid/row payloads — **also always fresh copies, never
  views.** Typed Rust errors map to a thrown JS error / `JsValue`.

**No zero-copy exists anywhere in this codebase today** — not in NumPy
arrays, not in WASM typed arrays. Every array crossing a language boundary
is a fresh allocation and copy of the underlying Rust data. This is a
current-state fact, not a promise about the future.

**No format auto-detection/dispatch exists anywhere.** Every parse function
name states the format it parses; there is no `parse_any(text)` that
sniffs format from content.

---

## Copy-vs-view semantics

| Boundary | Owns a view of Rust memory? | Notes |
|---|---|---|
| Python NumPy arrays (`ecfp4_numpy()`, `VolumetricGrid.values`, `VolumetricGrid.values_3d`, ...) | No — fresh copy | `values_3d` reshapes the flat `values` copy to `(nx, ny, nz)`, third axis (`k`) fastest, matching `chematic_mol::VolumetricGrid::checked_index`. |
| WASM `js_sys::Float64Array`/`Uint32Array` (`cube_values_f64`, `opendx_values_f64`, `lammps_dump_rows_f64`, `lammps_dump_cartesian_positions_f64`, `cube_shape_u32`, `opendx_shape_u32`) | No — fresh copy | Additive alongside the JSON-string functions; neither replaces the other. |
| WASM JSON strings (`cube_grid_json`, `opendx_grid_json`, `mmcif_to_json`, ...) | No — serialized copy | The oldest/default binding shape in `chematic-wasm`; large numeric arrays here round-trip through a JSON number array, a disclosed perf trade-off versus the typed-array functions above. |

---

## Unit-conversion ownership: who converts what

| Concept | Owned by | Where |
|---|---|---|
| Gaussian Cube Bohr/Ångström tag | Rust core | `GridUnits::{Bohr, Angstrom}` — read and preserved as-is, never silently converted, in Rust, Python, and WASM alike. |
| OpenDX Ångström assumption | Rust core | No in-file unit tag exists; every parsed grid is tagged `Angstrom` by convention in `chematic_mol::opendx`. This is a core-crate assumption, not a binding-layer one — Python and WASM inherit it unchanged. |
| QCSchema geometry Bohr → Ångström | **Binding layer**, not core | `chematic_mol::QcMolecule.geometry` stays in Bohr (a0), matching the QCSchema spec. Python and WASM bindings convert to Ångström for convenience when exposing coordinates — this conversion does **not** happen in `chematic-mol` itself. |
| mmCIF/PQR/ORCA coordinate units | Rust core | Always Ångström as documented per-format; no binding-layer conversion involved since there is nothing to convert. |

---

## `None` / `null` / `Err` divergence points

These are real, disclosed differences in how "no value" is represented
across the three languages for the *same* underlying Rust computation —
not accidents.

### LAMMPS dump Cartesian positions — the clearest case

`chematic_mol::LammpsDumpFrame::cartesian_positions()` returns
`Option<Vec<[f64; 3]>>` — `None` when the frame has no `x y z` or `xs ys zs`
columns (it deliberately does not fall back to `xu yu zu`; see
[`format-capabilities.md`](format-capabilities.md#lammps-dumptrajectory)).

| Language | Function | "no value" representation |
|---|---|---|
| Rust | `LammpsDumpFrame::cartesian_positions()` | `None` |
| Python | `LammpsDumpFrame.cartesian_positions()` | (Python-native `None`, same `Option` mapping) |
| WASM (JSON) | `lammps_dump_cartesian_positions_json` | JSON `null` |
| WASM (typed array) | `lammps_dump_cartesian_positions_f64` | **`Err`**, not `null` — a `Float64Array` has no `null` representation, so the unresolvable case becomes a thrown error with a message naming the columns it looked for. |

This is the one place in the codebase where the *same* WASM concept has two
different "no value" behaviors depending on which of its two sibling
functions you call — both are doc-commented at their definitions
(`crates/chematic-wasm/src/format_io.rs`), and this is deliberate: it is
not safe to assume every `*_json` / `*_f64` sibling pair behaves the same
way just because one does.

### Typed-error mapping

| Rust | Python | WASM |
|---|---|---|
| A typed error enum (e.g. `LammpsDataError::UnsupportedAtomStyle`, `OpenDxError::NonAngstromUnits`, `MmcifError::...`) | `ValueError` with the Rust error's `Display` text as the message | Thrown JS error / `JsValue::from_str(...)` with the same `Display` text |

Python and WASM represent the **same underlying Rust error type**
differently at the language boundary — `ValueError` vs. a JS exception —
even when the Rust-side error and its message text are identical. Neither
language currently exposes the original Rust error's structured
variant/fields across the boundary; both flatten to a string message.

### QCSchema unknown fields

`chematic_mol::qcschema`'s `unknown_fields: JsonObject` bag (kept distinct
from the spec's own open extensibility bags — `extras`/`keywords`/
`protocols`/`native_files`/`properties`) round-trips losslessly in **all
three** language bindings — Rust, Python (routed through Python's own
`json` module as the dict↔text boundary, not a hand-written field mapper),
and WASM. This is one of the few places where all three surfaces behave
identically rather than diverging; noted here for contrast with the
divergent cases above.

---

## Worked examples: same concept, three languages

### Gaussian Cube parse + grid values

| Language | Parse | Access grid values |
|---|---|---|
| Rust | `chematic_mol::parse_cube(text)` / `parse_cube_with_limits(text, limits)` → `VolumetricGrid` | `grid.values: Vec<f64>` (flat, row-major, third-axis-fastest — see `checked_index`) |
| Python | `chematic.VolumetricGrid.from_cube(text)` (staticmethod) | `grid.values` (flat NumPy copy) or `grid.values_3d` (reshaped `(nx, ny, nz)` NumPy copy, `k` fastest) |
| WASM | `cube_grid_json(text)` → JSON string | `cube_values_f64(text)` → `Float64Array` (flat, same ordering) or parse the JSON string's `values` field |

### OpenDX strict vs. lossy write

| Language | Strict (fails closed on non-Ångström units) | Explicit lossy (Bohr→Ångström, opt-in) |
|---|---|---|
| Rust | `write_opendx(&grid)` → `Err(OpenDxError::NonAngstromUnits)` for a `Bohr`-tagged grid | `write_opendx_lossy(&grid)` — rescales `origin`/`axes` only, never `values` |
| Python | `grid.to_opendx()` → raises `ValueError` for `Bohr` | `grid.to_opendx_lossy()` |
| WASM | `write_opendx_json(grid_json)` → `Err`/thrown error | `write_opendx_lossy_json(grid_json)` |

No language collapses these into a single lossy-by-default function — the
strict/lossy split from the Rust core is preserved identically in every
binding.

### LAMMPS Cartesian positions

| Language | Function | Notes |
|---|---|---|
| Rust | `LammpsDumpFrame::cartesian_positions()` | `Option<Vec<[f64; 3]>>` |
| Python | `LammpsDumpFrame.cartesian_positions()` | pyclass method, same resolution rules |
| WASM (JSON) | `lammps_dump_cartesian_positions_json(frame_json)` | `null` on unresolvable |
| WASM (typed array) | `lammps_dump_cartesian_positions_f64(frame_json)` | `Err` on unresolvable — see divergence table above |

---

## Streaming vs. materialization, by language

| Format | Rust | Python | WASM |
|---|---|---|---|
| MOL/SDF | `SdfFileReader<R: BufRead>` — true streaming `Iterator` | `iter_sdf` and `iter_sdf_batched` are file-backed streaming iterators; batches are lazy, cancellable, and expose a progress manifest | materializes |
| LAMMPS dump/trajectory | `LammpsDumpReader<R: BufRead>` — true streaming `Iterator` | `parse_lammps_dump_all` materializes the whole trajectory as a list (disclosed scope choice, not a silently dropped capability) | `lammps_trajectory_to_json` materializes (same disclosed choice) |
| Gaussian Cube | `CubeFileReader<R: BufRead>` streams the *input reading* only — the returned `VolumetricGrid.values` is still one fully-materialized `Vec<f64>` (single-dataset format, nothing to iterate across) | via `VolumetricGrid.from_cube()`, materializes | materializes |
| Other documented formats | no `BufRead`-backed streaming reader type exists | materializes | materializes |

LAMMPS dump is the one format where a real Rust-level streaming/
materializing distinction exists **and** both bindings deliberately choose
materialization. This is a documented compatibility choice.

---

## Fail-closed behavior (binding-independent)

These fail-closed checks live in the Rust core and are inherited unchanged
by both bindings — no binding loosens or bypasses them:

- `write_opendx` on a `Bohr`-tagged grid → typed error in all 3 languages (see OpenDX table above).
- Both OpenDX writers on a grid with any atoms → `OpenDxError::AtomsNotSupported` in all 3 languages (no lossy-atom-dropping path exists at all).
- Gaussian Cube multi-dataset input → `CubeError::MultiDatasetUnsupported` in all 3 languages (typed rejection, not silent truncation).
- LAMMPS data with an unrecognized `atom_style` → `LammpsDataError::UnsupportedAtomStyle` in all 3 languages (no best-effort guess).
- LAMMPS "Type Labels" sections → rejected in all 3 languages.

See [`error-and-limits.md`](error-and-limits.md) for the full typed-error
taxonomy and parse-limit reference.
