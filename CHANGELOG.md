# Changelog

All notable changes to chematic will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

### Added

#### chematic-chem — CIP stereochemistry (Phase 3 completion)
- `assign_cip(mol: &Molecule) -> CipAssignment` — assigns R/S (tetrahedral) and E/Z (double bond) CIP codes:
  - BFS sphere expansion with phantom atoms for double bonds and ring revisits.
  - Tetrahedral R/S via OpenSMILES @/@@ parity with correct bracket-H insertion rule.
  - E/Z from Up/Down stereo bonds on double-bond endpoints.
- `CipAssignment::get(idx: AtomIdx) -> Option<CipCode>` accessor.
- `CipCode` enum (R, S, E, Z) added to `chematic-core`; re-exported from both crates.
- 19 new tests; chematic-chem total: 67.

#### chematic-smarts — MCS (Phase 4)
- `find_mcs(mols: &[&Molecule]) -> QueryMolecule` — McGregor connected-growth MCS.
- `find_mcs_with_config(mols, config) -> QueryMolecule` with `McsConfig { match_bonds, min_atoms, timeout_ms }`.
- Branch-and-bound pruning via element-count upper bound; `std::time::Instant` timeout.
- `QueryMolecule::atom_count()` accessor added.
- 12 new tests; chematic-smarts total: 46.

#### chematic-chem — tautomer normalization (Phase 4)
- `canonical_tautomer(mol: &Molecule) -> Molecule` — fixed-point rule-based canonical form.
- `enumerate_tautomers(mol: &Molecule) -> Vec<Molecule>` — BFS enumeration, max 32.
- 5 rules: keto-enol, amide-iminol, imine-enamine, 1,3-H-shift N→O, 1,3-H-shift N→N.
- 10 new tests.

#### chematic-mol — MOL V2000 stereo bond parsing
- Bond block stereo field (columns 9-11) now parsed: stereo=1/4 → `BondOrder::Up`, stereo=6 → `BondOrder::Down`.
- Backward compatible: lines shorter than 12 chars default to stereo=0.
- 2 new tests; chematic-mol total: 36.

#### chematic-fp — MACCS and topological path fingerprints (Phase 4)
- `maccs(mol) -> BitVec2048` — MACCS 166-bit structural keys fingerprint (`maccs.rs`):
  - All 166 SMARTS patterns evaluated via the existing `chematic-smarts` VF2 engine.
  - Bit `i` set when MACCS key `i+1` matches the molecule (at least one occurrence).
  - Key 164 corrected to `[!#6;!#1]` (standard MDL heteroatom detector); fixes zero
    fingerprint for simple alcohols like ethanol.
  - Silent fallback on unparseable patterns (rare; none currently fail).
  - `chematic-smarts` promoted from dev-dep to production dep in `chematic-fp/Cargo.toml`.
- `topo_path(mol, &TopoPathConfig) -> BitVec2048` — topological path fingerprint (`topo_path.rs`):
  - Enumerates all simple paths of 2–`max_len` atoms via DFS (default `max_len = 7`).
  - Path encoded as interleaved `[atomic_num, bond_order, atomic_num, ...]` bytes.
  - Canonicalized by taking the lexicographically smaller of forward and reverse encodings.
  - Hashed with FNV-1a 64-bit, folded into `BitVec2048` via `hash % nbits`.
- `TopoPathConfig { max_len: usize, nbits: usize }` — configurable path length and output size.
- Both modules now exported from `chematic-fp/src/lib.rs` as `pub mod maccs`, `pub mod topo_path`
  with `pub use` re-exports (`maccs`, `topo_path`, `TopoPathConfig`).
- 13 new tests across `maccs` (7) and `topo_path` (6) modules; total test count: 250 → 263.

#### chematic-mol (extended)
- `parse_mol_v3000(input) -> Result<(Molecule, MolMetadata), MolParseError>` in `mol3000.rs`:
  - Two-phase parser: pre-pass collects and joins `M  V30 ` continuation lines (trailing `-`).
  - State machine: `BeforeCtab` → `InCtab` → `InAtomBlock` → `AfterAtomBlock` → `InBondBlock` → `Done`.
  - Supports `CHG=`, `MASS=`, `HCOUNT=`, and `aamap` key-value fields.
  - Errors on missing `END ATOM` or `END BOND`.
- `V3000ParseError { line: usize, msg: String }` variant added to `MolParseError`.
- `#![forbid(unsafe_code)]` added crate-wide.

#### chematic-depict (new crate)
- `compute_layout(mol) -> Layout` — rule-based 2D coordinate generation:
  - Ring placement: regular polygon with radius `BOND_LEN / (2 sin(PI/n))`.
  - Fused ring placement: centroid-based outward direction, signed-angle CW/CCW selection.
  - Zigzag chain placement: ±30° alternating DFS traversal, `BOND_LEN = 40.0` px.
  - Fragment offset: components separated by 2×BOND_LEN gap.
- `render_svg(mol, layout) -> String` — SVG serializer:
  - Single bonds: `<line stroke-width="1.5">`.
  - Double/triple bonds: parallel offset lines (±2 px / ±3 px).
  - Aromatic bonds: solid + dashed parallel lines.
  - Wedge (Up): filled `<polygon>` triangle.
  - Dash (Down): series of short transverse bars.
  - Atom labels: element symbol + H count for non-C atoms; white background rect.
  - Rendering order: bonds → background rects → labels.
- `depict_svg(mol) -> String` — convenience wrapper: calls `compute_layout` then `render_svg`.

#### chematic-chem (new crate)
- `molecular_weight(mol) -> f64` — average isotopic mass including implicit H.
- `exact_mass(mol) -> f64` — monoisotopic mass; respects `atom.isotope`.
- `heavy_atom_count(mol) -> usize`.
- `hbd_count(mol) -> usize` — N/O atoms with H count > 0.
- `hba_count(mol) -> usize` — all N and O atoms.
- `rotatable_bond_count(mol) -> usize` — non-ring single bonds between non-terminal atoms; amide C–N excluded.
- `tpsa(mol) -> f64` — Ertl (2000) atom-type lookup table.
- `logp_crippen(mol) -> f64` — simplified Crippen-Wildman atom contributions.
- `lipinski_passes(mol) -> bool` — MW ≤ 500, HBD ≤ 5, HBA ≤ 10, LogP ≤ 5.
- Key design: kekulize before H-count-sensitive descriptors (aromatic bonds `order_int=1` overcounts).

#### chematic-fp (new crate)
- `BitVec2048` — 2048-bit bitvector (`[u64; 32]`) with `set`, `get`, `popcount`, `and`, `or`, `fold`, `tanimoto`, `dice`.
- `EcfpConfig { radius: u32, nbits: usize }` — configurable radius and bit count.
- `ecfp(mol, config) -> BitVec2048` — FNV-1a 64-bit Morgan iteration:
  - Initial invariants: `atomic_number`, `degree`, `h_count`, `charge+8`, `is_in_ring`, `is_aromatic`.
  - Double-buffered ID arrays to avoid intra-pass contamination.
  - Canonical neighbor ordering: sorted `(bond_type_int, neighbor_id)` pairs.
- `ecfp4(mol) -> BitVec2048` — radius=2, 2048 bits.
- `ecfp6(mol) -> BitVec2048` — radius=3, 2048 bits.
- `tanimoto_ecfp4(a, b) -> f64` — convenience similarity function.

#### chematic-smarts (new crate)
- `QueryMolecule` — query graph with `AtomQuery`/`BondQuery` logical trees.
- `AtomPrimitive` variants: `AtomicNum`, `Symbol`, `Aromatic`, `Charge`, `HCount`, `Degree`, `RingMembership`, `RingSize`, `Wildcard`.
- `BondPrimitive` variants: `Single`, `Double`, `Triple`, `Aromatic`, `Any`, `Ring`.
- `parse_smarts(s) -> Result<QueryMolecule, SmartsError>` — recursive-descent parser:
  - Organic-subset shorthands: `C` → `And(Symbol("C"), Aromatic(false))`, `c` → aromatic.
  - Bracket atoms with full precedence: `!` > juxtaposition/`&` > `,` > `;`.
  - Ring closures, branches, and explicit bond tokens.
- `find_matches(query, mol) -> Vec<HashMap<usize, AtomIdx>>` — VF2 subgraph isomorphism:
  - `EvalCtx` caches `find_sssr` once per call.
  - Injective mapping; bond compatibility checked against already-mapped neighbors.

#### chematic-3d (new crate)
- `Point3 { x, y, z }` — 3D vector with full linear-algebra ops (add, sub, scale, dot, cross, norm, normalize).
- `Coords3D` — indexed by `AtomIdx`; wraps `Vec<Point3>`.
- `generate_coords(mol) -> Coords3D` — rule-based DFS 3D coordinate builder:
  - Ideal bond lengths by element-pair + bond order.
  - Rodrigues rotation formula for bond-angle placement (sp3=109.5°, sp2=120°, sp=180°).
  - Ring templates placed as regular polygons in XY plane (aromatic C–C = 1.40 Å).
  - Disconnected components offset +5 Å along X.
- `parse_pdb_atoms(s) -> Vec<PdbAtom>` — parses ATOM/HETATM fixed-column records.
- `pdb_to_molecule(atoms) -> (Molecule, Coords3D)` — distance-based bond inference (1.3× sum of covalent radii).
- `write_pdb(mol, coords) -> String` — HETATM records, fixed-column PDB format.
- `parse_xyz(s) -> Result<(Molecule, Coords3D), XyzError>` — XYZ format parser.
- `write_xyz(mol, coords, comment) -> String` — XYZ format writer.

### Planned
- Phase 5 remaining: UFF force field minimization
- Phase 6 remaining: WASM package (npm: chematic), ChEMBL-scale validation

---

## [0.1.0] — 2026-05-26

Initial release covering Phase 1 (foundation) and Phase 2 (molecular perception + file I/O).

### Added

#### chematic-core 0.1.0
- `Element` newtype (`Element(u8)`) covering all 118 elements of the periodic table.
  - `from_symbol(s)` case-sensitive lookup; `symbol()` returns canonical symbol string.
  - `atomic_number()`, `is_organic_subset()`, `normal_valences()` accessors.
  - Organic subset: B, C, N, O, F, P, S, Cl, Br, I.
- `Atom` struct with fields: `element`, `isotope`, `charge` (i8), `hydrogen_count` (Option<u8>),
  `aromatic` (bool), `chirality` (Option<Chirality>), `wildcard` (bool), `atom_map` (u16).
  - Constructors: `Atom::new()`, `Atom::organic()`, `Atom::aromatic()`, `Atom::bracket()`, `Atom::wildcard()`.
- `BondOrder` enum: `Single`, `Double`, `Triple`, `Quadruple`, `Aromatic`, `Up`, `Down`.
  - `order_int()` method mapping aromatic/single=1, double=2, triple=3.
- `Bond` and `BondEntry { atom1: AtomIdx, atom2: AtomIdx, order: BondOrder }`.
- `Molecule` with adjacency-list graph (no petgraph); `AtomIdx(u32)` and `BondIdx(u32)` newtypes.
  - `atom()`, `bond()`, `neighbors()`, `atom_count()`, `bond_count()`, `formula()` (Hill order).
- `MoleculeBuilder` with `add_atom()`, `add_bond()`, `build()`, `atom_at()`, `atom_neighbors()`.
- `implicit_hcount(mol, idx) -> u8` in `valence` module.
  - Bracket atoms: returns stored explicit H count.
  - Organic-subset atoms: computes from normal valence table with formal charge adjustment.
  - Wildcard atoms and non-organic-subset atoms: returns 0.
- `kekulize(mol) -> Result<KekuleResult, KekuleError>` in `kekulization` module.
  - Augmenting-path maximum matching on the aromatic subgraph.
  - Lone-pair donors (O, S, Se, pyrrole-type N) are optional in the matching.
  - `apply_kekule(mol, kekule) -> Molecule` rebuilds molecule with double/single bonds assigned.
- 30 unit tests covering element lookups, valence calculations, and kekulization of
  benzene, pyridine, furan, pyrrole, and naphthalene.

#### chematic-smiles 0.1.0
- OpenSMILES parser (`parse(s) -> Result<Molecule, SmilesError>`):
  - Organic subset atoms (B, C, N, O, P, S, F, Cl, Br, I) with implicit aromaticity inference.
  - Aromatic atoms (c, n, o, p, s) with automatic aromatic bond inference between adjacent aromatics.
  - Bracket atoms `[isotope?symbol±charge:hcount@chirality:map]` with full field parsing.
  - Wildcard atom `[*]` via `Atom::wildcard()`.
  - Ring closures: single-digit (`C1...C1`) and two-digit (`C%10...C%10`).
  - Branch notation (`C(CC)CC`).
  - Disconnected fragments (`.` separator).
  - Tetrahedral stereo (`@`, `@@`) parsed and stored on Atom.
  - Bond types: `-`, `=`, `#`, `$`, `:`, `/`, `\`.
- SMILES writer (`write(mol) -> String`):
  - Depth-first traversal with correct ring-closure numbering.
  - Branches wrapped in parentheses; canonical child ordering.
  - Bond order symbols elided for single bonds (except explicit hydrogen notation).
- Canonical SMILES (`canonical_smiles(mol) -> String`):
  - Morgan rank algorithm: FNV-1a hash propagation over atomic invariants.
  - Initial invariants: atomic number, degree, formal charge, isotope, H count, aromaticity.
  - Tie-breaking: atomic number, isotope, charge, aromaticity, degree (no atom-index dependence).
  - Stable across roundtrips for aspirin, caffeine, glucose, naphthalene, disconnected molecules.
- 50 tests: roundtrip parsing for aspirin, caffeine, glucose, NaCl; canonical SMILES stability;
  wildcard atoms; stereo; multi-ring systems.

#### chematic-perception 0.1.0
- `find_sssr(mol) -> RingSet` — Smallest Set of Smallest Rings:
  - BFS spanning forest to find r = edges - atoms + components fundamental cycles.
  - LCA-based path reconstruction to get cycle bond sets.
  - GF(2) Gaussian elimination (XOR on sorted `Vec<BondIdx>`) selects r linearly independent rings.
  - `RingSet` API: `rings()`, `ring_count()`, `contains_atom()`, `atoms_in_ring_count()`.
- `assign_aromaticity(mol) -> AromaticityModel` — Hückel 4n+2 aromaticity:
  - Calls `find_sssr` internally; checks sp2 compatibility of each ring atom.
  - Pi electron contribution: C(double bond neighbor)=1, pyridine-N=1, pyrrole-N(H)=2, O=2, S=2.
  - Hückel criterion: `pi_count >= 2 && (pi_count - 2) % 4 == 0`.
  - Supports: benzene, pyridine, pyrrole, furan, thiophene, imidazole, naphthalene, indole, quinoline.
  - `AromaticityModel { aromatic_atoms: HashSet<AtomIdx>, aromatic_bonds: HashSet<BondIdx> }`.
- 14 tests covering benzene, pyridine, pyrrole, furan, cyclopentadiene, cyclohexane,
  naphthalene, indole, and non-aromatic ring systems.

#### chematic-mol 0.1.0
- MOL V2000 (CTfile) parser (`parse_mol(s) -> Result<(Molecule, MolMetadata), MolParseError>`):
  - Header block: molecule name, program/timestamp, comment lines.
  - Counts line: atom count, bond count, chiral flag.
  - Atom block: fixed-column x/y/z coordinates, element symbol, mass difference, charge code.
  - Bond block: atom indices (1-based), bond type (1-4), stereo flag.
  - Charge codes: 0=0, 1=+3, 2=+2, 3=+1, 5=-1, 6=-2, 7=-3.
  - Bond types: 1=Single, 2=Double, 3=Triple, 4=Aromatic.
  - `M  END` terminator.
- MOL V2000 writer (`write_mol(mol, metadata) -> String`):
  - Outputs valid CTfile with zero 2D/3D coordinates.
  - Charge code back-conversion from formal charge.
  - Correct 1-based atom indexing in bond block.
- SDF multi-molecule reader:
  - `SdfReader<'a>` iterator splitting on `$$$$` delimiter.
  - `parse_sdf(s) -> Result<Vec<(Molecule, MolMetadata)>, MolParseError>` for bulk loading.
- `MolMetadata { name, comment, extra_lines }` carrying header information.
- 19 tests: MOL parsing, charge handling, aromatic bonds, multi-molecule SDF iteration,
  writer roundtrip, error cases.

### Technical decisions
- Zero C/C++ FFI: entire codebase is pure Rust.
- WASM-compatible: no `std::fs`, no threads in core or perception crates.
- No petgraph: custom adjacency-list graph with chemical semantics embedded in types.
- `AtomIdx(u32)` / `BondIdx(u32)` newtypes prevent index-confusion bugs at compile time.
- `#![forbid(unsafe_code)]` on all crates.
- FNV-1a hashing for reproducible, deterministic canonical SMILES across platforms.

[Unreleased]: https://github.com/chematic-rs/chematic/compare/v0.1.0...HEAD
[0.1.0]: https://github.com/chematic-rs/chematic/releases/tag/v0.1.0
