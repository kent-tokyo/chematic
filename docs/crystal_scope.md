# `chematic-crystal` scope

`chematic-crystal` is a structural/geometric foundation for periodic
(crystalline) materials, meant to be consumed by future higher-level
projects (`mikiwame` for explainable materials diagnostics, `gugen` for
synthesis/process planning) without those projects, or their concerns,
living in this crate.

## In scope (v0.1)

- `Lattice`: 3x3 matrix + inverse, cell-parameter constructors, volume,
  reciprocal vectors, fractional<->Cartesian conversion, validation
  (rejects `NaN`/`Infinity`/non-positive lengths/degenerate angles/singular
  and near-singular matrices).
- `FractionalCoord` / `CartesianCoord`: distinct newtypes so the two spaces
  can't be silently mixed up at a call site.
- `PeriodicSite` / `SiteSpecies` / `Occupancy`: one or more (element,
  occupancy) pairs per site, occupancy-sum validation, room for future
  disorder without a breaking change.
- `PeriodicStructure`: lattice + sites, `validate()`, coordinate wrapping,
  periodic-neighbor enumeration, diagonal supercell generation.
- Exact PBC displacement/minimum-image (bounded exhaustive search derived
  from reciprocal vectors -- see `docs/rfcs/chematic_crystal_foundation.md`),
  verified against a from-scratch brute-force oracle, including triclinic
  cells.
- Diagonal supercells (`[nx, ny, nz]`, each `>= 1`).
- Optional `serde` feature (JSON round-trip through validated constructors,
  not raw field derives).
- `wasm32-unknown-unknown` build target.
- POSCAR/CONTCAR (VASP structure file format) read/write (`src/poscar.rs`,
  added after v0.1's initial scope) -- VASP 5 only (explicit species-name
  line required, VASP 4's implicit POTCAR-derived ordering rejected with a
  typed error), single- and 3-component scale factors, Direct/Cartesian
  coordinate modes, selective dynamics, ion velocities, and CONTCAR's
  verbatim-preserved predictor-corrector MD-restart data. No VASP
  INCAR/KPOINTS/POTCAR parsing. See `crates/chematic-crystal/src/poscar.rs`'s
  module docs for the full list of format-fidelity decisions.

## Out of scope (v0.1) -- do not assume these exist

- Symmetry: no space-group determination, no symmetry-operation *search*
  (deriving operations from a name/International Tables number), no
  Wyckoff positions, no primitive/conventional cell conversion, no Niggli
  reduction, no spglib (or any) FFI. This still holds inside
  `chematic-crystal` itself. **Not** the same thing as *applying*
  operations a CIF already states literally in its own text: that
  capability was added in `chematic-mol` (`crates/chematic-mol/src/
  cif_symmetry.rs`, consumed by `parse_cif_periodic_structure_with_options`)
  and consumes only this crate's existing public API
  (`minimum_image`/`PeriodicSite::new`/`PeriodicStructure::new`, all
  unchanged) — `chematic-crystal` gained no new symmetry code and remains
  as scoped above.
- No CIF parser rewrite in `chematic-crystal` itself.
  `chematic_mol::cif::{parse_cif, write_cif, UnitCell}` (the original,
  non-periodic small-molecule reader) are untouched. The
  `parse_cif_structure() -> PeriodicStructure` adapter this section used
  to describe as a future sketch has since been implemented, in two
  stages, entirely inside `chematic-mol`: `parse_cif_periodic_structure`
  (see `docs/rfcs/chematic_crystal_foundation.md`'s "CIF migration"
  section for the original sketch) and, later,
  `parse_cif_periodic_structure_with_options`'s explicit
  symmetry-operation expansion described above.
- No XRD simulation, no crystal fingerprinting.
- No VASP INCAR/KPOINTS/POTCAR parsing (POSCAR/CONTCAR structure I/O is in
  scope, see above).
- No oxidation-state inference, no coordination-geometry classification.
- No DFT, no formation-energy or phase-diagram computation, no band
  structure, no phonons.
- No defect generation, no surface/slab construction.
- No materials-ML models.
- No prediction of stability, synthesizability, or any other materials
  property -- this crate stores and computes pure geometry, nothing more.
- No `mikiwame` diagnostics, no `gugen` synthesis planning -- those are
  separate crates/projects layered on top, not part of this repository's
  scope in this PR.
- No Python (`chematic-py`), WASM (`chematic-wasm`), or MCP
  (`chematic-mcp`) bindings.
- No arbitrary 3x3 integer supercell transform, only diagonal
  `[nx, ny, nz]`.
- No spatial-partitioning / cell-list neighbor-search optimization --
  correctness first, the baseline neighbor search is an exact bounded
  enumeration, not necessarily the fastest one possible.

If a task or PR touching `chematic-crystal` starts to need any of the above,
that is a signal to stop and open a new RFC rather than grow this crate's
scope inline.
