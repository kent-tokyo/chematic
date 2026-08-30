# chematic-crystal

Periodic crystal structure representation and geometry for the
[chematic](https://github.com/kent-tokyo/chematic) ecosystem: lattices,
fractional/Cartesian coordinates, periodic sites, periodic boundary
conditions, exact minimum-image distance, periodic neighbor enumeration,
diagonal supercells, and POSCAR/CONTCAR (VASP structure file format)
read/write. Pure Rust, `#![forbid(unsafe_code)]`, WASM-compatible, one
required dependency (`chematic-core`, for `Element`).

`v0.1` is a structural/geometric foundation, not a symmetry or materials-
property library -- see "Out of scope" below and
`docs/crystal_scope.md` in the workspace root for the public scope and
limitations.

## What this crate is not: `Molecule`

`chematic_core::Molecule` is a bond graph: atoms, typed bonds, aromaticity,
implicit hydrogens, stereo perception. `PeriodicStructure` is a lattice +
sites (element/occupancy + fractional position) with **no bond graph at
all**. These are deliberately separate first-class types, not variants of
each other:

- A crystal's "neighbor" relationship is a geometric/distance fact under
  periodic boundary conditions, not chemical bonding. This crate never
  auto-derives a `Bond` from periodic proximity.
- `chematic-crystal` does not depend on `chematic-core::Molecule`, and
  nothing in `chematic-core::Molecule` gains periodic-structure fields from
  this crate.

## Quick start

```rust
use chematic_core::Element;
use chematic_crystal::{
    FractionalCoord, Lattice, PeriodicSite, PeriodicStructure, SiteSpecies,
};

let lattice = Lattice::cubic(5.64).unwrap();
let sites = vec![
    PeriodicSite::new(
        vec![SiteSpecies::full(Element::NA)],
        FractionalCoord::new([0.0, 0.0, 0.0]),
        Some("Na1".to_string()),
    )
    .unwrap(),
    PeriodicSite::new(
        vec![SiteSpecies::full(Element::CL)],
        FractionalCoord::new([0.5, 0.5, 0.5]),
        Some("Cl1".to_string()),
    )
    .unwrap(),
];
let structure = PeriodicStructure::new(lattice, sites).unwrap();

let neighbors = structure.neighbors_within(5.0).unwrap();
let supercell = structure.make_supercell([2, 2, 2]).unwrap();
```

See `examples/basic_structure.rs` (`cargo run -p chematic-crystal --example
basic_structure`) for a runnable version with printed output.

## Conventions

### Lattice matrix

`Lattice::matrix()` is `[[f64; 3]; 3]` with **rows = lattice vectors a, b,
c** (`matrix[0]` = a, `matrix[1]` = b, `matrix[2]` = c). Cartesian is
obtained from fractional by row-vector x matrix:

```text
cartesian = fractional . matrix
cartesian_k = sum_j fractional_j * matrix[j][k]
```

This matches `chematic_ewald::BoxVectors`'s existing row convention
(`chematic-crystal` does not depend on `chematic-ewald`; the convention is
just kept compatible for a possible future bridge). Reciprocal vectors use
the crystallographic (no `2*pi`) convention `a_i . b_j = delta_ij`; under
the row convention, reciprocal row `b_j` is **column `j`** of
`inverse_matrix()`, not row `j`.

### Coordinates and units

- `FractionalCoord([f64; 3])`: dimensionless, lattice-relative. Not
  range-restricted by the type -- use `.wrapped()` to reduce into `[0, 1)`
  (via `rem_euclid(1.0)`; `1.0` maps to `0.0`, negative values wrap up).
- `CartesianCoord([f64; 3])`: Angstrom -- the same unit convention as the
  rest of chematic (`chematic_core::coords3d::Point3`, MOL/SDF, PDB).
- Public constructors (`Lattice::from_matrix`/`from_parameters`/`cubic`/
  `orthorhombic`, `PeriodicSite::new`, `PeriodicStructure::new`) reject
  `NaN`/`Infinity`. The coordinate newtypes' fields are `pub` (matching the
  crate's immutable-by-default design), so a value built via direct struct
  literal isn't automatically checked -- call `.is_finite()` or run it
  through a structure's `validate()` if it didn't come from a validated
  constructor.

### PBC distance / minimum image

`periodic::minimum_image(lattice, from, to)` returns the **exact** nearest
periodic image -- not an approximation. A naive `delta_frac -=
delta_frac.round()` reduction is only guaranteed correct for cells close to
orthogonal; for sufficiently skewed triclinic cells it can pick the wrong
image (verified concretely by a pinned regression fixture in
`tests/periodicity.rs`: a real triclinic cell where naive rounding reports
8.92 Angstrom for a pair whose true nearest-image distance is 1.73
Angstrom). Instead, `minimum_image` derives a finite, provably sufficient
search box from the lattice's own reciprocal vectors and checks every
candidate inside it -- exact for any lattice this crate's validation
accepts. See `docs/crystal_scope.md` for the public scope and limitations.

`PeriodicDisplacement` reports which image was chosen (`image: [i32; 3]`),
the Cartesian and fractional displacement, and the distance. Ties (multiple
images at the same minimal distance) are broken deterministically by
ascending iteration order.

### Neighbor images

`PeriodicStructure::neighbors_within(cutoff)` (Angstrom, **inclusive**:
`distance <= cutoff`) returns a **full** neighbor list:

- `(same site, image [0,0,0])` is always excluded (the trivial self-pair).
- `(same site, nonzero image)` is **not** auto-excluded -- physically valid
  when the cell is smaller than the cutoff.
- For distinct sites `i != j`, both `(center=i, neighbor=j, image=n)` and
  `(center=j, neighbor=i, image=-n)` can appear as independent entries --
  each describes a different center's neighbor shell, not a duplicate.
- Distance `0` (or near-`0`) for distinct sites/images is kept, not
  filtered -- e.g. two disorder-split sites at the same position.
- Output is sorted by `(center_index, neighbor_index, image[0], image[1],
  image[2])` for determinism.
- Uses the same reciprocal-vector search-box machinery as `minimum_image`,
  with the cutoff as the bound; a cutoff far larger than the cell (a likely
  unit-conversion mistake) is rejected with `CrystalError::
  NeighborSearchTooLarge` rather than silently examining an enormous number
  of candidate images.

### Occupancy

`SiteSpecies { element, occupancy }` -- a `PeriodicSite` holds `Vec
<SiteSpecies>` so multi-species (disordered) sites are representable.
`Occupancy` is validated finite and `>= 0`; the *sum* of occupancies at one
site must not exceed `1.0 + Occupancy::SUM_TOLERANCE` (summing to less than
`1.0` is legal -- a vacancy). An empty species list is rejected
(`CrystalError::EmptySpeciesList`), not silently treated as a vacant site.
`v0.1` stores this data faithfully but does not implement disorder-aware
chemistry (no averaged scattering, no partial-occupancy energetics) on top
of it.

### POSCAR/CONTCAR

`poscar::{parse_poscar, parse_contcar, write_poscar}` and the
`PoscarDocument`/`PoscarError`/`PredictorCorrector` types read and write
VASP's plain-text structure format (`parse_contcar` is a thin alias --
CONTCAR is the same format). VASP 5 only (an explicit species-name line is
required; a VASP 4-style file with implicit POTCAR-derived ordering is
rejected with `PoscarError::Vasp4NotSupported` rather than silently
mis-parsed). Both scale-factor forms from the VASP wiki are supported on
read (single value, including the negative "target cell volume" form, and
the 3-component per-axis form); `write_poscar` always emits `1.0` with
pre-scaled vectors and `Direct` (fractional) coordinates -- the simplest
form that's always exactly correct, since `PeriodicStructure` stores
fractional coordinates canonically. Selective dynamics and ion velocities
round-trip; CONTCAR's predictor-corrector MD-restart section is preserved
verbatim (its numeric layout is not documented by VASP itself, "cannot be
entered by hand" per the wiki, so this reader doesn't attempt to interpret
it) and refuses to be written back if doing so would require silently
reordering that opaque data. POSCAR has no disorder/partial-occupancy
concept, so every parsed site is a single, fully-occupied species, and
`write_poscar` rejects a multi-species or partially-occupied site rather
than dropping data. No per-atom labels (VASP has no field for one beyond
the species symbol). Full format-fidelity decision list:
`crates/chematic-crystal/src/poscar.rs`'s module docs.

## Out of scope (v0.1)

No space-group determination, symmetry-operation search, Wyckoff
positions, primitive/conventional cell conversion, Niggli reduction, or
spglib (or any) FFI. No CIF parser rewrite (`chematic_mol::cif` is
untouched; see `docs/crystal_scope.md` for the supported migration path).
No VASP INCAR/KPOINTS/POTCAR parsing (POSCAR/CONTCAR structure I/O
is in scope, see above). No XRD, crystal fingerprinting, oxidation-state
inference, coordination-geometry classification, DFT, formation-energy or
phase-diagram computation, band structure, phonons, defect generation, or
surface/slab construction. No materials-ML models and **no prediction of
stability or synthesizability** -- this crate stores and computes pure
geometry. No WASM/MCP bindings (Python bindings for `Lattice`/
`PeriodicStructure`/`Site` shipped in v0.17.0 via `chematic-py`, see
`CHANGELOG.md`). No arbitrary 3x3 integer supercell
transform (diagonal `[nx, ny, nz]` only). No spatial-partitioning neighbor-
search optimization -- the baseline is an exact bounded enumeration, not
necessarily the fastest possible one.

`mikiwame` (explainable materials diagnostics) and `gugen` (materials
synthesis/process planning) are separate downstream projects layered on top
of this crate; neither is implemented here. Full list:
`docs/crystal_scope.md`.

## Benchmarks

Baseline measurements only -- `v0.1` prioritizes correctness over
performance (see "no spatial-partitioning optimization" above). Method:
`criterion` (`cargo bench -p chematic-crystal`), except the ~10000-site
neighbor search, which is a single wall-clock timing (see
`benches/crystal_bench.rs`'s module doc for why: the naive O(n^2 *
search-box) baseline takes long enough per call at that size that
criterion's minimum sampling would make the benchmark run impractically
slow, for a number this PR isn't claiming is fast).

CPU: Apple M4. rustc 1.97.0. Commit: `fa63822` (chematic-crystal Phase 4).

| Operation | Cell shape | Sites | Cutoff | Neighbors found | Method | Time |
|---|---|---|---|---|---|---|
| `frac_to_cart` (1000 points, 1 lattice) | triclinic | -- | -- | -- | criterion | 2.14 us total (~2.1 ns/point) |
| `cart_to_frac` (1000 points, 1 lattice) | triclinic | -- | -- | -- | criterion | 4.81 us total (~4.8 ns/point) |
| `neighbors_within` | cubic | 125 (~100) | 6.5 A | 4000 | criterion | 10.09 ms |
| `neighbors_within` | triclinic | 125 (~100) | 6.5 A | 7500 | criterion | 9.02 ms |
| `neighbors_within` | cubic | 1000 | 6.5 A | 32000 | criterion | 499.4 ms |
| `neighbors_within` | triclinic | 1000 | 6.5 A | 60000 | criterion | 585.3 ms |
| `neighbors_within` | cubic | 10648 (~10000) | 6.5 A | 340736 | single-run wall clock | ~54.3 s |
| `neighbors_within` | triclinic | 10648 (~10000) | 6.5 A | 638880 | single-run wall clock | ~47.2 s |
| `make_supercell([2,2,2])` | cubic | 125 -> 1000 | -- | -- | criterion | 101.4 us |

The O(n^2) scaling is expected and visible above (125 -> 1000 sites, an 8x
increase, is roughly a 50-65x slowdown -- cubic 49.5x, triclinic 64.9x --
bracketing the 64x pair-count increase (8^2) expected from O(n^2) scaling
at this fixed cutoff/density ratio); a spatial-partitioning neighbor list
(cell lists / Verlet lists) is
the natural next step if this becomes a bottleneck for a real workload, but
is explicitly deferred past v0.1.

## Testing

`cargo test -p chematic-crystal --all-features` runs unit tests (lattice
construction/validation, coordinate wrapping, occupancy, structure
validation, POSCAR/CONTCAR read/write -- including a triclinic
parse-write-parse round trip, both scale-factor forms, selective dynamics,
velocities, and predictor-corrector verbatim round-tripping), integration
tests (`tests/periodicity.rs`, `tests/neighbor.rs`: brute-force-oracle
comparisons on cubic/orthorhombic/triclinic cells, two pinned regression
fixtures for skewed/near-singular triclinic minimum-image, a randomized
fixed-seed property test, site-order-permutation invariance), and (with
`--features serde`) `tests/serde_roundtrip.rs` (JSON round-trip, field-name
stability, invalid-value rejection). `cargo test -p chematic-crystal --doc`
runs the doc examples above.

## Optional `serde` support

Enable the `serde` feature for `Serialize`/`Deserialize` on every public
type. Every implementation is hand-written, not derived: this crate's
types enforce invariants at construction (finite coordinates, non-negative
occupancy, occupancy-sum tolerance, matrix non-singularity), and a plain
`#[derive(Deserialize)]` would silently skip all of them by assigning
fields directly -- every `Deserialize` impl here routes through the type's
own validated constructor instead. `SiteSpecies` round-trips its `Element`
through `symbol()`/`from_symbol()` as a string (`chematic-core` has no
serde support of its own). `Lattice` persists only `matrix`, not the
derived `inverse` -- deserializing re-derives and re-validates it via
`from_matrix`.

If you serialize/deserialize through `serde_json` specifically and need
exact float round-tripping, enable serde_json's `float_roundtrip` feature
-- its default float parser is not guaranteed correctly-rounded (observed
directly while writing this crate's own round-trip tests: a
`from_parameters`-derived matrix entry serialized as
`3.453503973595536e-16` parsed back as `3.4535039735955363e-16`, 1 ULP
off, without that feature).
