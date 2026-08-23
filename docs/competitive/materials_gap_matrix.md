# Materials-Science Gap Matrix (pymatgen / ASE / spglib)

Method: every row is grounded in `crates/chematic-crystal/`, `crates/chematic-mol/`'s
CIF/POSCAR/LAMMPS readers, and `crates/chematic-ewald/` as they exist on
`main` at 2026-08-23 (v0.19.0). No row is based on README/marketing text.

**Scope note, stated once rather than per-row**: pymatgen, ASE, and spglib
are not one thing to be "ahead of" or "behind" on a single axis. pymatgen
is a broad, Python-native materials-analysis ecosystem (structure
manipulation, thermodynamics, electronic-structure-code interop, phase
diagrams) built over ~15 years by a large community. ASE is a
simulation-orchestration framework (the `Calculator` abstraction, MD/
optimization drivers, dozens of DFT/force-field code adapters) rather than
a structure-representation library per se. spglib is a single-purpose,
highly-optimized C library for symmetry/space-group determination, wrapped
by both of the above. chematic-crystal is a young, single-Rust-crate,
single-maintainer effort (test count: 115, vs. these projects' orders-of-
magnitude larger surfaces) that currently exists to give chematic's
Rust/Python/WASM molecular runtime a periodic-structure foundation, not to
replace any of the three. A row classified "Missing" here means "does not
exist in chematic today," not "chematic is failing at its stated goals" —
most of these capabilities were never in chematic-crystal's original scope
(`docs/crystal_scope.md`'s own non-goal list already discloses most of
this honestly). Read the matrix as a map of what a *materials-science
expansion* would need to build, not a scorecard chematic is currently
losing.

Classification: **Ahead** / **Competitive** / **Partial** / **Missing** /
**Not planned**.

| # | Capability | Status | Evidence |
|---|---|---|---|
| 1 | Periodic structure representation | Competitive | `Lattice` (`lattice.rs:69`, 17 tests, triclinic-capable via `from_parameters`), `FractionalCoord`/`CartesianCoord` (`site.rs:17`/`:59`, distinct newtypes), `SiteSpecies`/`Occupancy`/`PeriodicSite` (`site.rs:96-157`, multi-species disorder support with `Occupancy::SUM_TOLERANCE = 1e-6` validation), `PeriodicStructure` (`structure.rs:17`). 115 tests across the crate incl. `tests/serde_roundtrip.rs` (14) and `tests/periodicity.rs` (9). |
| 2 | Minimum image (periodic distance) | Ahead (exactness, not speed) | `periodic::minimum_image` (`periodic.rs:96`) is an **exact**, reciprocal-vector-bounded exhaustive search over candidate periodic images — not the common `frac -= frac.round()` approximation, which is known to fail for skewed triclinic cells. 6 tests including a real tie-break regression fixture (`cubic_half_cell_tie_uses_lexicographically_smallest_image`) found via triclinic testing. This exactness-first design is a genuine, deliberate strength — matched by pymatgen/ASE in correctness but the point is chematic didn't take the common shortcut either. |
| 3 | Neighbor search | Partial | `neighbors_within` (`neighbor.rs:57`) is confirmed **O(n²) all-pairs** (nested nested loop over every site pair, `neighbor.rs:94-95`), not a cell-list/Verlet-list. The crate's own scope doc states this directly: "correctness first... not necessarily the fastest one possible" (`docs/crystal_scope.md:80-82`). A safety cap (`MAX_NEIGHBOR_IMAGE_CANDIDATES = 5,000,000`) guards against pathological cutoff/cell mismatches, not performance. 14 tests (10 in `neighbor.rs` + 4 in `tests/neighbor.rs`). pymatgen/ASE both use cell-list-based neighbor search for large structures. |
| 4 | Arbitrary supercell construction | Missing (diagonal-only exists) | `make_supercell` (`supercell.rs:27`) takes only `[nx, ny, nz]` — the module's own doc states outright "Arbitrary 3x3 integer supercell transforms are out of scope" for this version, confirmed again in `docs/crystal_scope.md:78-79`. 10 tests, all diagonal-scaling scenarios. pymatgen's `make_supercell` accepts arbitrary integer transformation matrices. |
| 5 | Cell reduction (Niggli/Minkowski) | Missing | Zero hits workspace-wide for "Niggli"/"Minkowski"/"reduced cell"/"lattice reduction" outside the crate's own non-goal disclaimer (`lib.rs:15`). No implementation anywhere. |
| 6 | Symmetry operations (apply declared operations) | Competitive | `expand_sites` (`chematic-mol/src/cif_symmetry.rs:783`) expands a CIF's asymmetric unit into a full cell using operations **literally parsed from the CIF text** (`_space_group_symop_operation_xyz`/`_symmetry_equiv_pos_as_xyz`, both tag-spelling conventions). Public entry: `parse_cif_periodic_structure_with_options` (`cif.rs:966`). Landed v0.17.0. 34+39 tests across the two files. Golden fixtures: `iucr_p21c_four_operation_golden_fixture` (P2₁/c No. 14, 4 IUCr-sourced operations, self-verified for group closure) and `C2C_CIF` (synthetic C2/c No. 15, 8-operation loop driving actual end-to-end expand+special-position-dedup, 8→4 sites). |
| 7 | Space-group determination (arbitrary structure, spglib-equivalent) | **Missing, and deliberately so** | Confirmed absent and actively *tested against*: `declared_space_group_with_no_operation_list_at_all_stays_unexpanded` (`cif.rs:1758`) proves a CIF that names a space group but doesn't list operations stays `UnexpandedSymmetry` by design — chematic will not look up operations from a name/number table. `docs/crystal_scope.md:41-52` states this explicitly: "no space-group determination... no spglib (or any) FFI." `ExpandedExplicitOperations`'s own doc disclaims completeness against any space-group database. This is the single largest, hardest gap in this entire matrix (see RFC's B3 discussion) — do not conflate row 6 (apply-what's-declared) with this row (derive-from-nothing); they are different problems of very different difficulty. |
| 8 | Wyckoff positions | Missing | Depends on row 7 (space-group determination); not attempted. |
| 9 | Primitive/conventional cell | Missing | Zero hits workspace-wide for "primitive_cell"/"conventional_cell"; explicitly disclaimed (`docs/crystal_scope.md:41-44`). |
| 10 | Structure matching / duplicate detection (periodic) | Missing | No `StructureMatcher`-equivalent or lattice-aware-tolerance structure comparison exists for `PeriodicStructure`. The only "duplicate detection" hits in the workspace are molecular-graph/conformer concepts (`chematic-3d/conformer.rs`, `chematic-chem/hash.rs`) — unrelated to periodic-lattice comparison. |
| 11 | Crystal fingerprint | Missing | No periodic-structure fingerprint exists. (Extensive molecular fingerprint infrastructure exists in `chematic-fp` — irrelevant to this row.) |
| 12 | XRD / neutron diffraction pattern calculation | Missing | Zero hits for "XRD"/"diffraction"/"structure_factor" as a diffraction concept. The one `compute_structure_factor` hit (`chematic-ewald/src/pme.rs:305`) computes `\|S(k)\|²` as an internal term of Particle Mesh Ewald electrostatics summation — unrelated to X-ray form factors, Bragg angles, or powder-pattern output. Explicitly disclaimed (`docs/crystal_scope.md:63`). |
| 13 | Coordination environment analysis | Missing | Zero hits for periodic coordination-geometry classification. (Molecular "coordination number" hits elsewhere in the workspace — e.g. MMFF94 atom typing, a Pt coordination-number benchmark — are unrelated node-degree concepts, not periodic coordination-environment analysis.) |
| 14 | Oxidation-state inference | Missing | Zero hits for real oxidation-state chemistry; the only "oxidation state" grep hit is a CIF atom-label-parsing comment about stripping trailing digits, unrelated to chemistry. Explicitly disclaimed (`docs/crystal_scope.md:66`). |
| 15 | Surfaces / slabs | Missing | Zero functional hits for "slab"/"surface" construction in `chematic-crystal`. Explicitly disclaimed (`docs/crystal_scope.md:69`). |
| 16 | Defects (vacancies, substitutions) | Missing | Zero functional hits for "defect"/"vacancy" construction logic (the only "vacancy" hits are a doc comment describing sub-1.0 `Occupancy` as vacancy *modeling*, and an unrelated POSCAR test-fixture comment string — neither is defect-*construction* logic). Explicitly disclaimed (`docs/crystal_scope.md:69`). |
| 17 | Phase diagrams / convex hull | Missing | Zero hits workspace-wide for "convex_hull"/"phase_diagram," including `chematic-ff` — no reusable general convex-hull utility exists anywhere to build on. Explicitly disclaimed (`docs/crystal_scope.md:67`). |
| 18 | Calculator abstraction (pluggable energy/force backend) | Missing | Zero hits for a `Calculator`/`EnergyModel`/`ForceField` trait abstraction anywhere in the workspace. ASE's `Calculator` interface has no chematic analog. |
| 19 | Geometry optimization (building block, not periodic-aware) | Partial | Real minimizers exist and are well-tested: `minimize_uff` (steepest descent, adaptive step), `minimize_mmff94_full` (steepest descent), `minimize_mmff94_lbfgs` (real L-BFGS, quasi-Newton, Armijo line search) — but all in `chematic-ff`, which has **zero dependency on `chematic-crystal`** (confirmed via `Cargo.toml` + source grep) — none can accept a `Lattice`/`PeriodicStructure` or compute PBC-aware forces/energy. These are the nearest reusable building blocks for a future periodic optimizer, not existing periodic capability. |
| 20 | Molecular dynamics (building block, not periodic-aware) | Partial | Real velocity-Verlet MD exists (`chematic-3d/src/md.rs::run_md`, genuine half-step velocity update → position update → force recompute → velocity completion, NVE + NVT Berendsen thermostat, 3 tests) — but operates on `chematic_core::Molecule` with a hand-rolled bond+angle+DREIDING-vdW+Gasteiger-Coulomb energy, **zero coupling to `chematic-crystal`**, no periodic boundary conditions. |
| 21 | NEB (nudged elastic band) | Missing | Zero hits for "NEB"/"nudged elastic" anywhere in the workspace. |
| 22 | Phonon calculation | Missing | Zero hits for "phonon" anywhere in the workspace. |
| 23 | VASP interoperability | Competitive | POSCAR/CONTCAR read/write (`poscar.rs`, 1587 lines, 27 tests): `parse_poscar`/`parse_contcar`/`write_poscar`, VASP 5 explicit-species format only (VASP 4 implicit-ordering rejected with a typed error). Landed v0.16.0. No INCAR/KPOINTS/POTCAR parsing (deliberately out of scope — those are calculation-control files, not structure files). |
| 24 | Quantum ESPRESSO interoperability | Missing | Zero hits for "quantum espresso"/"pwscf"/`.pwi`/`.pwo` anywhere in the workspace. |
| 25 | LAMMPS interoperability | Competitive | Data-file (`chematic-mol/src/lammps_data.rs`, `parse_lammps_data`/`write_lammps_data`, 17 tests) and dump/trajectory (`lammps_dump.rs`, streaming reader + writer, 18 tests) both exist, landed v0.17.0, including a real triclinic box-bounds↔true-box conversion (`box_bounds_to_true`/`true_to_box_bounds`) with a cross-language (Rust/Python/WASM) parity fixture that caught a real bug (PR #343). Both modules are standalone (no `chematic_core::Molecule`/bond-perception integration — raw atom-index topology only). |
| 26 | Python/WASM bindings for periodic-structure types | Partial | PyO3 bindings for `Lattice`/`PeriodicStructure` exist (`crates/chematic-py/src/crystal.rs`, added v0.17.0) — this corrects a stale line in `docs/crystal_scope.md` claiming no Python bindings exist. **WASM and MCP bindings for any crystal type are confirmed still absent** (zero references in `chematic-wasm`/`chematic-mcp` source) — a real, currently-open gap for any future browser-facing materials feature. |

---

## Summary

| Area | Overall | Priority for Track B |
|---|---|---|
| Representation + exact minimum-image (rows 1-2) | Ahead/Competitive | Foundation — reused, not redone |
| Neighbor search performance (row 3) | Partial | **B1** |
| Arbitrary supercell (row 4) | Missing | **B1** |
| Cell reduction, primitive/conventional, structure matching (rows 5, 9-10) | Missing | **B2** |
| Symmetry determination + Wyckoff (rows 7-8) | Missing, hardest item in the matrix | **B3** |
| Fingerprint/coordination/oxidation-state/XRD (rows 11-14) | Missing | **B4** |
| Slabs/defects (rows 15-16) | Missing | **B5** |
| Phase diagrams, Calculator trait, MD/NEB/phonons (rows 17-22) | Missing/Partial | **B6** |
| VASP/LAMMPS I/O (rows 23, 25) | Competitive | Reused as-is by B6, not rebuilt |
| QE I/O (row 24) | Missing | **B6** |
| WASM/MCP crystal bindings (row 26) | Partial | Not in B1-B6; a smaller, separate binding-layer item |
