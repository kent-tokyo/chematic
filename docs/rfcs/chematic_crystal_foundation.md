# chematic-crystal foundation

Branch: `feat/chematic-crystal-foundation`. Scope: a new crate providing pure
structural/geometric representation of periodic (crystal) materials --
lattice, fractional/Cartesian coordinates, periodic sites, PBC displacement,
periodic neighbor enumeration, diagonal supercells. No symmetry, no CIF
rewrite, no `Molecule` changes. This document is the Phase 0 audit +
design record the implementation phases (1-5) build against.

## Why a new crate, not new `Molecule` fields

`chematic_core::Molecule` is a bond graph: atoms + typed edges, with
aromaticity, implicit-hydrogen inference, and stereo perception all defined
in terms of that graph. A crystal's "neighbors" are a function of geometry
and a periodic lattice, not chemical bonding -- two atoms 3.2 Angstrom apart
across a cell boundary are a coordination-sphere fact, not evidence of a
`Bond`. Encoding periodicity onto `Molecule` would force every consumer of
`Molecule` (SMILES, aromaticity, CIP, fingerprints, depiction -- the entire
existing crate graph) to reason about a periodicity dimension that is
meaningless to them, and would tempt future code into auto-deriving `Bond`
edges from periodic proximity, silently smuggling molecular semantics into
what should stay a geometry fact.

`PeriodicStructure` is therefore a distinct first-class type: lattice +
sites (element/occupancy + fractional position), no bond graph at all. It
depends on `chematic-core` only for `Element` (and, per audit below, cannot
depend on it for anything else without adding an unwanted dependency).
`Molecule` is untouched by this work -- no new fields, no new variants.

## Type ownership

| Type | Crate | Concept |
|---|---|---|
| `Molecule` | `chematic-core` | bond graph, chemistry semantics |
| `Element` | `chematic-core` | shared, reused as-is |
| `Coords3D`/`Point3` | `chematic-core` | atom-indexed 3D storage for `Molecule` conformers -- **not** reused here (see below) |
| `Lattice`, `PeriodicSite`, `PeriodicStructure`, `FractionalCoord`, `CartesianCoord`, `Occupancy` | `chematic-crystal` (new) | periodic structure |
| `UnitCell` (existing) | `chematic-mol::cif` | legacy CIF-parser-local cell parameters; unchanged (see CIF section) |
| `BoxVectors` | `chematic-ewald` | MD-simulation periodic box; unchanged, no dependency edge to/from `chematic-crystal` |

**`chematic_core::coords3d::Coords3D` was audited and deliberately not
reused.** It stores one `Point3` per atom indexed by `AtomIdx` (insertion
order into a `Molecule`) -- its own doc comment states it exists so
`chematic-mol` doesn't have to pull in `chematic-3d`'s heavier dependency
chain for a plain per-atom Cartesian point. `PeriodicStructure` has no
`Molecule`/`AtomIdx` in the picture at all, and needs **fractional**
coordinates as the primary representation (Cartesian is derived), plus a
type-level distinction between the two coordinate spaces that `Point3` does
not carry. Reusing it would mean either (a) storing fractional coordinates
in a type named/shaped for Cartesian atom storage, confusing readers, or (b)
depending on `chematic-core::coords3d` for a type whose invariants don't
match this crate's needs. Writing `FractionalCoord`/`CartesianCoord` as
small local newtypes costs a few lines and keeps the distinction the spec
requires (never take a fractional value where Cartesian is expected, or vice
versa) enforceable by the type checker.

## Dependency direction

```
chematic-core  <---  chematic-crystal  <---  chematic (facade, optional feature)
chematic-core  <---  chematic-mol      (existing CIF parser, unchanged)
```

`chematic-crystal` depends on `chematic-core` only (plus optional `serde`).
It does **not** depend on `chematic-mol`, `chematic-ewald`, or any other
domain crate. `chematic-mol` does not depend on `chematic-crystal` in this
PR either -- see "CIF migration" below for why the natural
`chematic-mol -> chematic-crystal` adapter direction is deferred to RFC-only
this round.

## Matrix convention (binding for all of `chematic-crystal`)

`Lattice::matrix()` is `[[f64; 3]; 3]` with **rows = lattice vectors a, b,
c** (`matrix[0]` = a, `matrix[1]` = b, `matrix[2]` = c). A Cartesian point is
obtained from fractional coordinates by **row-vector x matrix**:

```
cartesian = fractional . matrix          (fractional treated as a row vector)
cartesian_k = sum_j fractional_j * matrix[j][k]
```

This is the same row-vector convention `chematic_ewald::BoxVectors` already
uses (`BoxVectors(pub [[f64; 3]; 3])`, rows = box vectors, audited in
`crates/chematic-ewald/src/lib.rs`) -- chosen deliberately so a future
`Lattice <-> BoxVectors` conversion helper (not implemented in this PR, out
of scope) is a straight field copy, no transpose. `chematic-crystal` does
**not** depend on `chematic-ewald` to get this; the convention is just kept
compatible in case a bridge is written later, by either side, without a
breaking transpose.

It's also the convention `chematic_mol::cif::UnitCell::frac_to_cart` already
implements algebraically (`X = a.fx + b.cos(gamma).fy + c.cos(beta).fz`,
i.e. `fx * a_vec + fy * b_vec + fz * c_vec`) -- `Lattice::from_parameters`
reproduces the same IUCr formula, so a CIF-parsed cell and a
`Lattice::from_parameters` cell built from the same six numbers agree
bit-for-bit-modulo-fp on Cartesian output.

**Reciprocal vectors** use the crystallographic (no `2*pi`) convention:
`a_i . b_j = delta_ij`. Under the row convention above, `matrix . inverse =
I`, so `b_j` (row `j` of the reciprocal matrix) is **column `j` of
`inverse`**, not row `j` -- this is a real transposition, not a naming
choice, and is asserted directly in `lattice.rs`'s test suite
(`a_i . b_j == delta_ij` for all 9 pairs) rather than trusted by inspection.

## Coordinate convention

- `FractionalCoord([f64; 3])`: dimensionless, lattice-relative. Not
  range-restricted by the type itself (a raw fractional value can legally
  fall outside `[0, 1)` -- e.g. mid-calculation before wrapping); `.wrapped()`
  reduces each component into `[0, 1)` via `rem_euclid(1.0)`.
- `CartesianCoord([f64; 3])`: Angstrom, same units as the rest of chematic
  (`chematic_core::coords3d::Point3`, MOL/SDF, PDB all use Angstrom).
- Constructors and structure-level `validate()`/`new()` reject non-finite
  (`NaN`/`Infinity`) components; the newtypes' fields stay `pub` per the
  spec sketch; unwrapped direct field construction is possible but every
  *validated* entry point (`Lattice::from_matrix`, `PeriodicSite::new`,
  `PeriodicStructure::new`) runs a finiteness check, and `validate()` is
  callable on-demand for structures assembled by hand.

## PBC displacement / minimum-image algorithm

**Decision: exact bounded enumeration, not `round()`-only minimum image.**
`delta_frac -= delta_frac.round()` is only exact for cells close to
orthogonal; for a sufficiently skewed triclinic cell the true nearest
periodic image is not the one obtained by rounding each fractional
component independently, because "nearest in fractional space per axis"
and "nearest in Cartesian space" diverge once lattice vectors are not close
to mutually perpendicular. The chosen algorithm derives a *provably
sufficient* finite search box from the lattice's own reciprocal vectors, so
it is exact for every valid (non-singular, `condition_indicator` above
threshold) lattice, not an approximation:

1. Reduce the raw fractional delta component-wise into `[-0.5, 0.5]`:
   `base = delta_frac - round(delta_frac)`. This is one legitimate periodic
   image (the "naive" one) with Cartesian distance `D0`.
2. For a displacement `r = (base + n) . M` (row convention) and reciprocal
   rows `b_j` (`a_i . b_j = delta_ij`), the identity `(base + n)_j = r . b_j`
   plus Cauchy-Schwarz gives `|base_j + n_j| <= |r| * |b_j|`. Any image
   at least as good as the naive one (`|r| <= D0`) must therefore have
   `n_j` inside `[-D0*|b_j| - base_j, D0*|b_j| - base_j]` for every axis `j`
   -- a finite, explicitly computed integer box.
3. Enumerate every integer point in that box (small for any well-conditioned
   cell -- a handful of candidates per axis), evaluate the true Cartesian
   distance for each, keep the minimum. Ties broken by lexicographically
   smallest image vector for determinism.

This is the spec's "recommended" option (derive a safe search range
mathematically, verify against a brute-force oracle) rather than the
"reduced-cell-only, fallback elsewhere" alternative -- the bound holds for
any non-singular lattice, so there is no fallback branch to maintain.
`neighbors_within` reuses the same bound with `D = cutoff` instead of `D0`
(a cutoff-radius neighbor search needs every image within the cutoff, not
just the closest one; the box-then-filter shape is identical).

Verification: `tests/periodicity.rs` and `tests/neighbor.rs` compare this
algorithm against a **from-scratch** brute-force oracle (hardcoded
`-4..=4` per-axis triple loop, sharing no code with the production path) on
cubic, orthorhombic, and randomized triclinic cells (fixed-seed PRNG, no new
dependency -- see "Property-based tests" in the crate's test files), plus a
pinned regression fixture: a specific skewed triclinic cell where naive
`round()`-only minimum image disagrees with the true nearest image, found by
randomized search and committed so the exact-vs-naive distinction is
regression-tested, not just asserted in prose.

## Occupancy model

`Occupancy(f64)`, finite, `>= 0`, no fixed upper bound on a single value
(placeholder/vacancy modeling can have a single low-occupancy species). Site
validity is a **sum** constraint: `SiteSpecies` occupancies at one
`PeriodicSite` must sum to `<= 1.0 + OCCUPANCY_SUM_TOLERANCE`. Summing to
less than 1.0 is legal (vacancy). An empty `species: Vec<SiteSpecies>` is
rejected as a construction error (`CrystalError::EmptySpeciesList`) rather
than silently treated as a vacant/ghost site, per the spec's requirement to
either forbid or explicitly define the empty case. `v0.1` does not implement
any disorder-aware chemistry on top of multi-species sites (no averaged
scattering, no partial-occupancy energetics) -- the type only guarantees the
data model won't need a breaking change when that lands.

## serde design note: `chematic-core::Element` has no serde support

`chematic-core` is an intentionally zero-external-dependency crate (its own
`lib.rs` doc: "Zero external dependencies: compiles to wasm32-unknown-unknown
without modification") and `Element` does not derive `Serialize`/
`Deserialize`. Adding a `serde` feature to `chematic-core` to support this
one new downstream crate would be a real change to a foundational,
widely-depended-on crate -- out of the declared scope ("Moleculeを拡張し
ない" extends in spirit to "don't touch chematic-core's dependency profile
either", and the task brief lists only `chematic-crystal`'s own files as
in-scope). Instead, `SiteSpecies` (the only type holding an `Element`
directly) implements `Serialize`/`Deserialize` **by hand**, round-tripping
through the element's existing public `symbol()`/`from_symbol()` API as a
string (`"Na"`, `"Cl"`, ...) rather than deriving through `Element`. Every
other `chematic-crystal` type with private/derived-state fields (`Lattice`
stores a redundant `inverse` alongside `matrix`; `PeriodicStructure`'s
constructor is where validation lives) also gets a hand-written
`Deserialize` that round-trips through the crate's own validated
constructors (`Lattice::from_matrix`, `PeriodicStructure::new`) rather than
a field-literal derive, so a deserialized value can never skip validation
that a normally-constructed one goes through. `Occupancy`,
`FractionalCoord`, `CartesianCoord` do the same for the same reason
(reject non-finite / out-of-range values on the way in, not just on the way
out).

## Neighbor semantics

- `(same site, image [0,0,0])` excluded (trivial self-pair).
- `(same site, nonzero image)` **not** auto-excluded -- physically valid for
  small cells (an atom interacting with its own periodic image).
- Distance `0` or near-`0` for *distinct* sites/images is kept, not filtered
  -- e.g. two disorder-split sites at the same fractional position.
- Cutoff comparison is `distance <= cutoff` (inclusive boundary).
- Output order is sorted by `(center_index, neighbor_index, image[0],
  image[1], image[2])` -- integers only, no float comparison, so the order
  is deterministic independent of platform/optimization-level float
  rounding differences.
- The neighbor list is a **full** (not half/bonded) list: for two distinct
  sites `i != j` within cutoff, both `(center=i, neighbor=j, image=n)` and
  `(center=j, neighbor=i, image=-n)` are legitimate, independent entries
  (each describes a different center's neighbor shell), not duplicates of
  each other.

## Facade feature: `crystal` added, and included in `full` (decided)

`crates/chematic/Cargo.toml` gains an optional `crystal =
["dep:chematic-crystal"]` feature and a `pub use chematic_crystal as
crystal` re-export gated on it, following the exact shape every other
facade feature already uses (`smiles`, `mol`, `chem`, ...).

Whether `crystal` also joins the `full` aggregate feature was initially
left open in this PR pending explicit human confirmation -- the project
owner has since decided **yes**: `full` is operated as the aggregate of
*every* optional facade capability, not scoped to `Molecule`-only
features, and a user who writes `features = ["full"]` expects everything,
`crystal` included. `chematic-crystal`'s architectural independence from
`Molecule` (see above) is why it's its own crate, not a reason to exclude
it from `full`. This does **not** change the facade's `default` feature
(still `[]`) -- only explicit `--features full` (or `--all-features`,
which already built `crystal` regardless of `full` membership either way)
users are affected. Convention going forward: a new optional facade
feature should default to joining `full` too, unless there's a specific
reason not to.

## WASM constraint

`chematic-crystal`'s only dependency is `chematic-core` (zero-dependency,
WASM-clean) plus optional `serde` (also WASM-clean, used elsewhere in this
workspace's WASM build already via `chematic-wasm`). `cargo check -p
chematic-crystal --target wasm32-unknown-unknown` is part of this PR's
quality gates (Phase 5). `chematic-wasm`/`chematic-py`/`chematic-mcp` are
not touched -- no bindings are added in this PR.

## CIF migration (RFC only -- no code change to `chematic-mol` in this PR)

`chematic_mol::cif` today:

- Parses `_atom_site_*` loops into a flat `Molecule` (atoms only, no bonds)
  plus a `Vec<(f64, f64, f64)>` of orthogonal Cartesian coordinates and an
  optional `UnitCell { a, b, c, alpha, beta, gamma }`.
- Performs no symmetry expansion -- effectively P1 (every atom in the file
  is taken literally; symmetry-equivalent positions implied by a space group
  are not generated).
- `UnitCell` carries only the six cell parameters (no matrix, no inverse, no
  reciprocal vectors) and exposes `volume()`/`frac_to_cart()` computed
  ad hoc from those six numbers each call.
- Does not preserve occupancy (`_atom_site_occupancy` is not read) or model
  disorder.
- Returns a `Molecule`, not a periodic-aware first-class type -- a
  CIF-parsed structure today has no distinction between "this is a periodic
  crystal" and "this is a normal isolated molecule with 3D coordinates
  attached", other than the caller happening to also have an `Option
  <UnitCell>` in hand.

**Semantic difference from the new `Lattice`/`PeriodicStructure`:**
`UnitCell` is parameters-only (recomputes derived quantities every call,
no reciprocal-vector API, no validation beyond what `parse_cif` itself
enforces at parse time); `Lattice` is a validated matrix + cached inverse
with reciprocal vectors, condition-number rejection, and the exact
minimum-image machinery above. They are not drop-in replacements for each
other today.

**Future API shape (not implemented, sketch only):**

```rust
// crates/chematic-mol/src/cif.rs -- future, not in this PR
pub fn parse_cif_structure(input: &str) -> Result<PeriodicStructureResult, CifError>;

pub struct PeriodicStructureResult {
    pub structure: chematic_crystal::PeriodicStructure,
    pub labels: Vec<String>,           // _atom_site_label, if present
    // occupancy would come from _atom_site_occupancy, read but currently discarded
}
```

This would require `chematic-mol` to depend on `chematic-crystal` (the
correct direction per this RFC's dependency-direction section -- never the
reverse). `chematic-mol` currently has **no `[features]` table at all** and
is a direct, unconditional dependency of `chematic-wasm`,
`chematic-fp`/`chematic-chem`/`chematic-3d` (via their own dep chains), and
`chematic-py` -- adding a new unconditional dependency (or introducing
`chematic-mol`'s first-ever optional feature, purely to gate this one
addition) is a structural change to a widely-depended-on crate's dependency
graph and build surface, which the task brief calls out explicitly as a
case to leave as RFC-only rather than implement
("この追加で既存feature graphやWASM buildへ影響が出る場合は、実装せずRFC
だけに留めて報告"). **Not implemented in this PR.** `parse_cif`/`write_cif`/
`UnitCell` are unchanged, not deprecated, not touched.

## Future symmetry boundary (non-goal, noted for later RFCs)

Space-group determination, symmetry-operation search, Wyckoff positions,
primitive/conventional cell reduction, Niggli reduction, and any spglib FFI
are explicitly out of scope for `chematic-crystal` v0.1 and are not stubbed
or scaffolded here. A future symmetry layer would most naturally sit
*above* `PeriodicStructure` (consuming it, producing symmetry metadata
alongside it or a reduced structure) rather than inside this crate, to keep
`chematic-crystal` a stable, dependency-light base that `mikiwame`/`gugen`
and a future symmetry crate can all build on independently. No interface for
that layer is designed in this PR.

## Non-goals (this PR)

space group / symmetry operations / Wyckoff positions / primitive-conventional
cell conversion / Niggli reduction / spglib FFI / CIF parser rewrite / POSCAR
parser / XRD / crystal fingerprinting / oxidation-state inference /
coordination-geometry classification / DFT / formation energy / phase
diagrams / band structure / phonons / defect generation / surfaces & slabs /
materials ML / `mikiwame` diagnostics / `gugen` synthesis planning /
Python/WASM/MCP bindings / arbitrary 3x3 integer supercell transforms
(diagonal only) / spatial-partitioning neighbor-search optimization
(baseline is a correct, not necessarily fast, exact search).
