# Format Capability Matrix

What each of chematic's supported file formats actually does, in Rust,
Python, and WASM — read/write coverage, streaming, coordinate
units, connectivity handling, round-trip fidelity, lossy operations, parse
limits, and known limitations.

See also: [`rdkit-migration.md`](rdkit-migration.md) for a feature-by-feature
Supported/Partial/Not-supported breakdown against RDKit, and
[`language-bindings.md`](language-bindings.md) for cross-language API
mapping (function names, return shapes, `None`/`null`/`Err` divergences).

This page describes what exists today. It does not claim "full support" for
any format — each row states exactly what is implemented and names the gaps.

### Common format conversion bridge

The Python package provides `convert_format(text, input_format, output_format)`
and the WASM package provides `convert_common_format(text, input_format,
output_format)` for explicit conversion between the most interoperable molecular formats:
SMILES, MOL/SDF (V2000), MOL V3000, MOL2, CML, ChemicalJSON, MolJSON, CDXML,
PDB, XYZ, PDBQT, and Gaussian input. Extensions are accepted as aliases (for
example `.mol2`, `smi`, and `com`). PDB/XYZ/PDBQT output requires 3D
coordinates; graph-only conversions do not claim to preserve format-specific
metadata. Unsupported formats fail before conversion with a `ValueError`.

The WASM bridge intentionally exposes the topology-only subset (SMILES, MOL,
MOL V3000, MOL2, CML, ChemicalJSON, MolJSON, and CDXML) so that JavaScript
and Python share a predictable core contract. Use the existing coordinate APIs
for PDB/XYZ/PDBQT in browser workflows.

This is intentionally a bounded interoperability layer, not an Open Babel
replacement. The format-specific APIs remain available when coordinates,
metadata, multiple records, or domain-specific options must be preserved.

---

## Support matrix

Read (R) / Write (W) per language. A dash means that operation is not
exposed in that language. "—" in the WASM columns for CIF means CIF has
**zero** WASM exposure at all, not an oversight in this table.

| Format | Rust R | Rust W | Python R | Python W | WASM R | WASM W |
|---|---|---|---|---|---|---|
| SMILES | Y | Y | Y | Y (`Mol.smiles` property) | Y | Y |
| SMARTS | Y | — (matching, not writing) | validation only | — | Y (match) | — |
| MOL/SDF | Y | Y | Y | Y (`Mol.to_mol_block()`/`.to_mol_block_2d()`) | Y | Y |
| PDB | Y | Y | Y | Y (`Mol.to_pdb(coords)`) | Y | Y (`ConformerHandle.get_conformer_pdb(idx)`) |
| mmCIF | Y | Y | Y | Y | Y | Y |
| CIF (small-molecule) | Y | Y | Y | — | — | — |
| PQR | Y | Y | Y | Y | Y | Y |
| XYZ | Y | Y | Y | — | Y | Y |
| Extended XYZ | Y | Y | Y | Y | Y | Y |
| QCSchema | Y | Y | Y | Y | Y | — (validation only, no writer) |
| ORCA input | Y | Y | Y | Y | Y | Y |
| ORCA output | Y | — (no writer exists in Rust) | Y | — | Y | — |
| Gaussian Cube | Y | Y | Y (via `VolumetricGrid`) | Y (via `VolumetricGrid`) | Y | Y |
| OpenDX | Y | Y (strict + `_lossy`) | Y (via `VolumetricGrid`) | Y (via `VolumetricGrid`, strict + `_lossy`) | Y | Y (strict + `_lossy`) |
| LAMMPS data | Y | Y | Y | Y | Y | Y |
| LAMMPS dump/trajectory | Y | Y | Y | Y | Y | Y |

Notes on the cells above that need qualification:

- **SMILES Python write / MOL-SDF Python write**: exposed as `Mol` object
  methods/properties (`Mol.smiles`, `Mol.to_mol_block()`/
  `.to_mol_block_2d()`), not as free `formats.rs` functions in the same
  shape as the read side (`from_smiles`, `from_mol_block`, etc.). Both
  directions exist; the binding shape just differs from the read side —
  see [`language-bindings.md`](language-bindings.md) for the exact names.
- **QCSchema WASM write**: `qcschema_validate_atomic_input` /
  `qcschema_validate_atomic_result` validate a JSON document; the only WASM
  writer present is `to_qcschema_molecule_json` (for `QcMolecule`, not
  `AtomicInput`/`AtomicResult`).

---

## Per-format detail

### SMILES

- **Rust**: `chematic_smiles::{parse, write, canonical_smiles, parse_cxsmiles, write_cxsmiles, parse_smi_file, write_smi_file, random_smiles}`.
- **Python**: `from_smiles`, `from_cxsmiles`, `is_valid_smiles`, `parse_smi_file`, `write_smi_file`.
- **WASM**: `parse_smiles`/`write_smiles`, `mol_block_from_smiles`, `smiles_array_to_sdf`, `parse_cxsmiles_json`, `normalize_cxsmiles`, `is_valid_smiles`.
- **Streaming**: `parse_smi_file`/`write_smi_file` materialize fully. Separately,
  `chematic-mol`'s `smiles_table` module (`SmilesRecordReader`/`SmilesRecordWriter`,
  a **different sub-format** — delimited SMILES-table files with name/property
  columns, not plain `.smi`) implements a true `BufRead`-backed streaming
  `Iterator`. It is not exposed to Python/WASM and is out of scope for the
  "15 formats" table this page otherwise tracks; noted here only because it
  bears directly on any "no format in this codebase streams" claim.
- **Coordinate units**: N/A — SMILES encodes graph topology and stereo parity,
  not coordinates. CXSMILES can carry optional 2D coordinates as an extension.
- **Connectivity**: native — bonds are the primary content of the format.
- **Round-trip**: graph-canonical, not byte-identical (whitespace/atom order
  are not preserved). `canonical_smiles()` is **not** guaranteed to be a safe
  dedup/cache key today — see the Known Limitations section below.
- **Lossy operations**: none inherent to the format itself.
- **Parse limits**: `SmilesParseLimits` controls input bytes, atom count, and
  bond count; `parse` applies finite safe defaults and
  `parse_with_limits` accepts a stricter policy.
- **Known limitations**: `canonical_smiles()` has a documented residual —
  isolated/simple E/Z double bonds can still produce two different, both-valid
  canonical strings for the same molecule in ~1 in 18 stereo-bearing
  molecules (measured on a 5,000-mol ChEMBL corpus; see README's "Known
  Limitations" section for the exact figures). Do not use it as a dedup key
  without accounting for this.

### SMARTS

- **Rust**: `chematic_smarts::{parse_smarts, find_matches, find_matches_with_rings, ...}`.
- **Python**: `is_valid_smarts` — validation only. Full parse/substructure
  matching is exposed through `Mol.find_matches()`/`Mol.has_substructure()`
  (on the molecule object, not as a `formats.rs`-style free function) — do
  not read this as "SMARTS has no Python support," it means the shape of the
  binding differs from the other formats in this table.
- **WASM**: `smarts_match_atoms`, `smarts_match_atoms_with_chirality`, `match_smarts_smiles`, `parse_cxsmarts_json`.
- SMARTS is a query language, not a molecule storage format — streaming,
  coordinate units, connectivity, round-trip, and lossy-operation columns
  don't apply in the same sense as the other 14 formats.
- **Parse limits**: `PdbParseLimits` bounds input bytes, physical line length,
  ATOM/HETATM records, and MODEL records; use
  `parse_pdb_atoms_with_limits` for a typed resource-limit error contract.

### MOL/SDF

- **Rust read**: `parse_mol`, `parse_mol_with_coords`, `read_mol_with_diagnostics`, `parse_mol_v3000*`, `SdfReader`, `SdfFileReader`, `SdfRecordReader`.
- **Rust write**: `write_mol`, `write_mol_with_conformer[_checked]`, `write_sdf*`, `write_mol_v3000*`.
- **Python**: `from_mol_block`, `from_mol_block_with_coords`, `from_mol_block_with_diagnostics`, `parse_sdf_with_coords`, `from_mol_v3000[_with_coords|_with_diagnostics]`.
- **WASM**: `mol_from_v3000_block`, `mol_from_sdf_block`, `to_mol_block`, `to_mol_v3000_block`, `sdf_to_smiles_json`, `sdf_to_records_json`, `sdf_from_records_json`, `mol_block_stereo_diagnostics_json`, `mol_v3000_stereo_diagnostics_json`, `mol_block_coords_json`.
- **Streaming**: `SdfFileReader<R: BufRead>` is a true I/O-streaming `Iterator`
  (does not require the whole file in memory up front). `SdfReader`/
  `SdfRecordReader` are lazy iterators over an already-loaded `&str` (do not
  eagerly collect every record into a `Vec`, but do require the full text in
  memory). Python and WASM bindings materialize (no streaming reader is
  bound in either language).
- **Coordinate units**: Ångström, per the Ctab standard.
- **Connectivity**: native Ctab bond table.
- **Round-trip**: semantic, not byte-identical — `write_mol` regenerates a
  Ctab from the `Molecule`'s own perceived bonds/coordinates, not a copy of
  the original text.
- **Lossy operations**: none named beyond the stereo-write scope below.
- **Square-planar stereo (`@SP1`/`@SP2`/`@SP3`)**: read automatically on
  every MOL/SDF parse. Write is **opt-in only**, via exactly 3 functions —
  `write_mol_with_conformer_checked`, `write_mol_v3000_with_conformer_checked`,
  `write_sdf_record_with_conformer_checked`. Do not describe MOL/SDF writing
  in general as "square-planar supported" — the plain `write_mol`/`write_sdf`
  path does not perceive or emit it.
- **Parse limits**: `SdfParseLimits` bounds streaming input bytes, individual
  record bytes, and yielded record count; use `SdfFileReader::with_limits` for
  an explicit policy. The in-memory readers remain intentionally separate.

### PDB

- **Rust**: `chematic_3d::{PdbAtom, parse_pdb_atoms, pdb_to_molecule, write_pdb}` —
  this format lives in `chematic-3d`, not `chematic-mol`; the one exception
  among these 15.
- **Python**: `from_pdb` (delegates to `chematic_3d`) for reading;
  `Mol.to_pdb(coords)` for writing.
- **WASM**: `mol_from_pdb`, `pdb_coords_json` for reading. Writing is
  exposed as a method, not a free function:
  `ConformerHandle.get_conformer_pdb(idx)` returns conformer `idx` as a PDB
  string (or `null` if `idx` is out of range), delegating to
  `chematic_3d::write_pdb` internally — this is real write capability, just
  shaped around an existing conformer handle rather than a
  `write_pdb_json(mol, coords)`-style free function.
- **Coordinate units**: Ångström (PDB standard).
- **Connectivity**: `pdb_to_molecule` perceives bonds from geometry
  (distance-based), unlike most of the other formats in this table, which
  either carry an explicit bond table or never infer one.
- **Round-trip**: not characterized in this pass — no documented round-trip
  guarantee found; treat as best-effort.
- **Lossy operations**: none named.
- **Parse limits**: no `*ParseLimits` type exists.

### mmCIF

- **Rust**: `chematic_mol::{parse_mmcif, parse_mmcif_with_limits, write_mmcif, MmcifAtomRecord, MmcifError, MmcifParseLimits, MmcifResult}`.
- **Python**: `parse_mmcif`, `write_mmcif`.
- **WASM**: `mol_from_mmcif`, `mmcif_coords_json`, `mmcif_to_json`, `write_mmcif_json`.
- **Streaming**: no streaming reader type exists for mmCIF.
- **Coordinate units**: always orthogonal Ångström (unlike small-molecule CIF's fractional coordinates).
- **Connectivity**: **no bond table, ever** — this module never infers bonds.
- **Round-trip**: lossy on read for the `label_*`/`auth_*` tag-pair collapse —
  `auth_*` is preferred when both are present, both are treated as one
  logical field, except `label_seq_id` (kept as its own field since it is
  legitimately `.` for non-polymer atoms even when `auth_seq_id` is normal).
  `write_mmcif` writes the single stored value into both columns, so a
  read→write round trip loses the distinction if the source file's `label_*`
  and `auth_*` genuinely disagreed. Unrecognized `_atom_site` columns are
  collected into `unhandled_columns` rather than silently dropped, but are
  not necessarily re-emitted verbatim on write.
- **Lossy operations**: the `label_*`/`auth_*` collapse above.
- **Parse limits**: `MmcifParseLimits { max_input_bytes, max_atoms, max_line_len }`.

### CIF (plain, small-molecule)

- **Rust**: `chematic_mol::{parse_cif, write_cif, CifError, CifResult, UnitCell}`, plus `#[cfg(feature = "crystal")]`-gated periodic-structure variants.
- **Python**: `parse_cif` — **read only**, no `write_cif` binding exists in `chematic-py`.
- **WASM**: **none.** Zero WASM exposure — confirmed by grep, not by omission from this table.
- **Coordinate units**: fractional cell coordinates converted to orthogonal Ångström on read (IUCr convention).
- **Connectivity**: **no bond table at all** — CIF's `_atom_site` loop carries only positions; this module returns atoms with no bonds.
- **Round-trip**: symmetry expansion is **not** performed on the plain (non-`crystal`-feature) path — only atoms literally listed in `_atom_site_*` are returned (effectively P1 treatment). See `docs/crystal_scope.md` for the `crystal`-feature periodic-structure variants.
- **Lossy operations**: symmetry-operation information is not applied/round-tripped on the plain path.
- **Parse limits**: `CifParseLimits { max_input_bytes, max_line_bytes,
  max_tokens, max_atoms }`; use `parse_cif_with_limits` for an explicit
  policy. The existing `parse_cif` parser uses finite defaults.

### PQR

- **Rust**: `chematic_mol::{parse_pqr, parse_pqr_with_limits, write_pqr, infer_element, PqrAtomRecord, PqrError, PqrParseLimits, PqrResult}`.
- **Python**: `parse_pqr`, `write_pqr`, `infer_element`.
- **WASM**: `mol_from_pqr`, `pqr_coords_json`, `pqr_to_json`, `write_pqr_json`, `pqr_infer_element`.
- **Coordinate units**: Ångström.
- **Connectivity**: **no bond table, none inferred** (unlike PDB, which does infer from geometry).
- **Element**: PQR has no element column in the format; `Element` is
  inferred from the atom name via `infer_element`, a documented heuristic
  matching classic PDB-parser convention for blank element columns.
- **Round-trip**: not characterized beyond the element-inference note above.
- **Lossy operations**: none named beyond element inference (which is a
  read-side reconstruction, not a write-side loss).
- **Parse limits**: `PqrParseLimits { max_input_bytes, max_atoms, max_line_len }`.

### XYZ / Extended XYZ

Two distinct functions share the name `parse_xyz` in different crates —
disambiguate by crate, not by name alone:

- `chematic_3d::parse_xyz -> Result<(Molecule, Coords3D), XyzError>` —
  infers bonds by distance, returns a full `Molecule`.
- `chematic_mol::parse_xyz -> Result<XyzFrame, XyzError>` — no `Molecule`,
  no bonds; the Extended-XYZ-superset document type. Also has
  `parse_extxyz`/`parse_extxyz_all`/`write_extxyz`/`write_xyz`,
  `XyzReader`/`ExtxyzReader`/`ExtxyzWriter`.
- **Python**: `from_xyz` uses `chematic_3d`'s bond-inferring version.
  `from_extxyz`/`from_extxyz_all`/`to_extxyz` use `chematic_mol`'s
  non-bond version.
- **WASM**: `mol_from_xyz`/`to_xyz` use `chematic_3d`'s version.
  `mol_from_extxyz`/`extxyz_frame_json`/`to_extxyz_json` use
  `chematic_mol`'s version.
- **Streaming**: `XyzReader`/`ExtxyzReader` are lazy iterators over an
  already-loaded `&str` (same category as `SdfReader`, not a `BufRead`-based
  reader). No `BufRead`-backed streaming type exists for XYZ. Python/WASM
  materialize.
- **Coordinate units**: Ångström (standard XYZ/extended-XYZ convention).
- **Connectivity**: `chematic_3d::parse_xyz` infers bonds by distance;
  `chematic_mol::parse_xyz`/`parse_extxyz` never do (no `Molecule` is even
  produced on that path).
- **Round-trip**: not characterized beyond the two-crate split above.
- **Lossy operations**: none named.
- **Parse limits**: `XyzParseLimits` bounds input bytes, atoms per frame, frame
  count, and physical line length. `parse_xyz_with_limits`,
  `parse_xyz_all_with_limits`, and `parse_extxyz_with_limits` accept an
  explicit policy; default single/all-frame parsers use finite defaults.

### QCSchema

- **Rust**: `chematic_mol::{QcMolecule, AtomicInput, AtomicResult, parse_qcschema_molecule, write_qcschema_molecule, parse_atomic_input, write_atomic_input, parse_atomic_result, write_atomic_result, chematic_to_qc_molecule, qc_molecule_to_chematic}`.
- **Parse limits**: `QcSchemaParseLimits` bounds input bytes, JSON nesting depth,
  array entries, and string bytes. `parse_qcschema_molecule_with_limits`,
  `parse_atomic_input_with_limits`, and `parse_atomic_result_with_limits` accept
  an explicit policy; the existing parsers use finite defaults.
- **Python**: same names, routed through Python's stdlib `json` module rather than a hand-mapped dict.
- **WASM**: `mol_from_qcschema_molecule`, `qcschema_molecule_coords_json`, `to_qcschema_molecule_json`, `qcschema_validate_atomic_input`, `qcschema_validate_atomic_result`.
- **Coordinate units**: `QcMolecule.geometry` is explicitly **Bohr (a0)** in Rust; Python and WASM bindings convert to Ångström for convenience.
- **Connectivity**: `connectivity: Option<Vec<(usize, usize, f64)>>` is the
  **only** one of these 15 formats where bonds are ever an optional,
  spec-native part of the document — used if present, never fabricated if
  absent.
- **Round-trip**: open extensibility bags (`extras`/`keywords`/`protocols`/
  `native_files`/`properties`) round-trip losslessly via `BTreeMap`, kept
  distinct from an `unknown_fields` bag for genuinely undocumented
  top-level keys.
- **Lossy operations**: none named — every numeric leaf is rejected if
  non-finite rather than silently coerced.
- **Parse limits**: no `*ParseLimits` type exists (size limits, if any, are
  whatever the caller or `serde_json` impose).

### ORCA input

- **Rust**: `chematic_mol::{parse_orca_input, write_orca_input, OrcaInput, OrcaInputError, OrcaBlock, OrcaCoords, OrcaAtom}`.
- **Python**: `parse_orca_input`, `write_orca_input`.
- **WASM**: `mol_from_orca_input`, `orca_input_coords_json`, `orca_input_to_json`, `write_orca_input_json`.
- **Coordinate units**: Ångström.
- **Connectivity**: no bond perception anywhere in this module.
- **Round-trip**: lossless preservation of unknown `%block ... end` blocks, including nested sub-blocks.
- **Lossy operations**: none named.
- **Parse limits**: no `OrcaParseLimits` type exists (unlike mmCIF/PQR/Cube/OpenDX).

### ORCA output

- **Rust**: `chematic_mol::{parse_orca_output, OrcaOutput, OrcaOutputError, OrcaTermination, OrcaOptConvergence}` — **read-only, no writer, in all 3 languages.**
- **Python**: `parse_orca_output`.
- **WASM**: `orca_output_to_json`.
- **Coordinate units**: Ångström.
- **Connectivity**: no bond perception.
- **Semantics**: `termination` and `optimization_convergence` are two
  independently-reported fields — `ORCA TERMINATED NORMALLY` alone does
  **not** imply a requested geometry optimization converged; check both.
- **Parse limits**: no `*ParseLimits` type exists.

### Gaussian Cube

- **Rust**: `chematic_mol::{parse_cube, parse_cube_with_limits, write_cube, CubeError, CubeFileReader, CubeParseLimits}`, plus the shared `VolumetricGrid` type.
- **Python**: no free functions — only via the `VolumetricGrid` pyclass: `VolumetricGrid.from_cube()` (staticmethod), `.to_cube()`, `.values_3d`.
- **WASM**: `mol_from_cube`, `cube_grid_json`, `write_cube_json`, plus typed-array `cube_values_f64`, `cube_shape_u32`.
- **Streaming**: `CubeFileReader<R: BufRead>` streams the *input reading*
  only (reduces peak parse-time memory) — the returned
  `VolumetricGrid.values` is always one fully-materialized `Vec<f64>`. This
  is a different sense of "streaming" than `SdfFileReader`/`LammpsDumpReader`
  (Cube is a single-dataset format, so there is nothing to iterate across).
- **Coordinate units**: `GridUnits::{Bohr, Angstrom}` explicit tag; values
  are never silently normalized between them. The sign of the first-axis
  voxel count in the header disambiguates Bohr vs. Ångström (documented
  against a real spec source).
- **Connectivity**: no bonds, ever. `GridAtom.charge` is effective nuclear
  charge, not partial/formal charge.
- **Multi-dataset files**: **typed-rejected**, not silently truncated —
  `CubeError::MultiDatasetUnsupported`. Cube is single-dataset-only in this
  codebase today.
- **Round-trip**: not byte-identical (no claim of verbatim footer/comment
  preservation is made in the module docs).
- **Lossy operations**: none named on the read/write path itself (the
  multi-dataset case is rejected, not silently truncated).
- **Parse limits**: `CubeParseLimits { max_input_bytes, max_atoms, max_grid_points }`.

### OpenDX

- **Rust**: `chematic_mol::{parse_opendx, parse_opendx_with_limits, write_opendx, write_opendx_lossy, OpenDxError, OpenDxParseLimits}`.
- **Python**: no free functions — only via the `VolumetricGrid` pyclass: `.from_opendx()`, `.to_opendx()`, `.to_opendx_lossy()`. No `mol_from_opendx` — OpenDX has no atom section at all.
- **WASM**: `opendx_grid_json`, `write_opendx_json`, `write_opendx_lossy_json`, plus typed-array `opendx_values_f64`, `opendx_shape_u32`.
- **Coordinate units**: OpenDX has **no in-file unit tag** — every parsed
  grid is tagged `Angstrom` by convention (an assumption about the dominant
  real-world producer, not information recovered from the file itself).
- **Connectivity**: no bonds — OpenDX carries no atom section.
- **Round-trip**: footer boilerplate is regenerated on write, not preserved byte-for-byte.
- **Lossy operations, explicitly named**:
  - `write_opendx` **fails closed** (`OpenDxError::NonAngstromUnits`) on a `Bohr`-tagged grid — it will not silently write a wrong-unit file.
  - `write_opendx_lossy` is the explicit opt-in Bohr→Ångström conversion. It rescales `origin`/`axes` only — it **never** rescales `values`.
  - Both writers refuse (`OpenDxError::AtomsNotSupported`) any grid with non-empty atoms — there is no lossy-atom-dropping path.
- **Parse limits**: `OpenDxParseLimits { max_input_bytes, max_grid_points }` — no `max_atoms` (the format has none).

### LAMMPS data

- **Rust**: `chematic_mol::{parse_lammps_data, write_lammps_data, LammpsData, LammpsDataError, LammpsAtom, LammpsAtomStyle, LammpsBond, LammpsBox, LammpsMass, LammpsVelocity}`.
- **Python**: `parse_lammps_data`, `write_lammps_data` (plain dict).
- **WASM**: `lammps_data_to_json`, `write_lammps_data_json`.
- **`atom_style` handling**: `atom_style` **cannot be inferred from the
  file** and must be passed explicitly. `atomic`/`charge`/`molecular`/`full`
  are supported; anything else maps to `LammpsAtomStyle::Other` and is
  rejected with `LammpsDataError::UnsupportedAtomStyle` — there is no
  best-effort guess, because `charge` and `molecular` rows are both
  genuinely ambiguous 6-field rows by column count alone.
- **Connectivity**: raw index topology (atom-type/bond indices), not
  chemically perceived — this is a standalone document type, not integrated
  with `Molecule`.
- **Section handling**: 4 sections are typed (`Masses`/`Atoms`/`Velocities`/
  `Bonds`); everything else is preserved byte-for-byte verbatim in
  `unparsed_sections` (an ordered `Vec<(String, String)>`, not a
  `HashMap`). The writer emits typed sections first, then opaque sections in
  their original relative order — interleaving between typed and opaque
  sections is **not** preserved, though each opaque section's own content is.
- **LAMMPS "Type Labels" extension**: explicitly **unsupported**, fails
  closed — any section whose name ends in `"Type Labels"` is rejected
  rather than silently mishandled.
- **Coordinate units**: none stated in-file — LAMMPS units come from the
  simulation's own `units` command, not this file type; chematic does not
  infer or convert them.
- **Round-trip**: opaque sections round-trip verbatim; typed sections
  round-trip through the typed representation (not guaranteed
  byte-identical, e.g. whitespace).
- **Parse limits**: no `LammpsParseLimits` type exists.

### LAMMPS dump/trajectory

- **Rust**: `chematic_mol::{parse_lammps_dump_frame, write_lammps_dump_frame, write_lammps_trajectory, box_bounds_to_true, true_to_box_bounds, LammpsDumpFrame, LammpsDumpReader, LammpsDumpError}`.
- **Python**: `parse_lammps_dump_frame`, `parse_lammps_dump_all` (**materializes the whole trajectory** — does not expose `LammpsDumpReader`'s streaming; a disclosed scope choice, already stated in CHANGELOG), `write_lammps_dump_frame`, `write_lammps_trajectory`, `box_bounds_to_true`, `true_to_box_bounds`; pyclass `LammpsDumpFrame` (`.column()`, `.cartesian_positions()`).
- **WASM**: `lammps_dump_frame_to_json_str`, `lammps_trajectory_to_json` (**also materializes, does not stream** — same disclosed choice), `write_lammps_dump_frame_json`, `write_lammps_trajectory_json`, `lammps_dump_cartesian_positions_json`, plus typed-array `lammps_dump_rows_f64`, `lammps_dump_cartesian_positions_f64`.
- **Streaming**: `LammpsDumpReader<R: BufRead>` is a **true streaming
  `Iterator`** over any `BufRead` — reads and yields one frame at a time
  without holding the whole trajectory in memory. Together with
  `SdfFileReader`, this is one of two `BufRead`-backed streaming readers
  among these 15 Rust formats, and the only one whose Python/WASM bindings
  both deliberately materialize instead of streaming (a disclosed
  trade-off, not an oversight — see CHANGELOG `[0.18.0]`/`[0.17.0]`).
- **Box bounds**: triclinic box-bounds↔true-box conversion
  (`box_bounds_to_true`/`true_to_box_bounds`) independently verified
  against LAMMPS's own documentation. `LammpsBox` is reused from the
  `lammps_data` module; only this module performs the bound↔true conversion.
- **`cartesian_positions()`**: resolves `x y z` (passthrough) or
  `xs ys zs` (triclinic-aware scaled transform), but **deliberately does
  not** resolve `xu yu zu` ("unwrapped" coordinates) — returns `None` if
  neither of the first two triples is present, even if unwrapped
  coordinates are, since conflating "current position" with "unwrapped
  position" would mislead. Callers use `.column("xu")` directly for that case.
- **Divergence between JSON and typed-array WASM bindings**: when a frame
  has no resolvable coordinate columns,
  `lammps_dump_cartesian_positions_json` returns JSON `null`, but
  `lammps_dump_cartesian_positions_f64` returns `Err` — a `Float64Array`
  has no `null` representation. Both are doc-commented; this is disclosed,
  not accidental. See [`language-bindings.md`](language-bindings.md).
- **Coordinate units**: not stated in-file (same as LAMMPS data).
- **Parse limits**: no `*ParseLimits` type exists.

---

## Formats that never fabricate bonds

mmCIF, PQR, ORCA (input and output), Gaussian Cube, OpenDX, LAMMPS data, and
plain CIF never infer or fabricate a bond table — they either have no bond
concept in the format at all, or (PQR specifically) have an element-inference
step that is documented separately from connectivity. PDB and
`chematic_3d::parse_xyz` are the exceptions among these 15: both infer bonds
from 3D geometry (distance-based), which is a documented, disclosed choice,
not a default you should assume applies elsewhere.

## Cross-reference

For the exact per-language function signatures, return-value shapes, and
`None`/`null`/`Err` divergence points summarized above, see
[`language-bindings.md`](language-bindings.md). For how these facts map to
RDKit's equivalent APIs, see [`rdkit-migration.md`](rdkit-migration.md).
