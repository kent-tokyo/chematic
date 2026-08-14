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

## Out of scope (v0.1) -- do not assume these exist

- Symmetry: no space-group determination, no symmetry-operation search, no
  Wyckoff positions, no primitive/conventional cell conversion, no Niggli
  reduction, no spglib (or any) FFI.
- No CIF parser rewrite. `chematic_mol::cif::{parse_cif, write_cif,
  UnitCell}` are untouched; a future `parse_cif_structure() ->
  PeriodicStructure` adapter is sketched (not implemented) in
  `docs/rfcs/chematic_crystal_foundation.md`.
- No POSCAR/VASP I/O, no XRD simulation, no crystal fingerprinting.
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
