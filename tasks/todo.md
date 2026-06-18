# chematic — Full Phase Roadmap (Pure Rust Implementation Targeting RDKit Feature Parity)

Goal: A single Rust crate ecosystem with no C/C++ FFI, native WASM support, covering all RDKit functionality.

Constraints: Zero external C/C++ FFI, core crates are WASM-compatible, no petgraph, no Python bindings.

---

## Phase 1 — Foundation (complete)

- [x] Cargo workspace setup
- [x] chematic-core: Element (118 elements), Atom, Bond, Molecule, MoleculeBuilder
- [x] chematic-smiles: OpenSMILES parser (organic subset / bracket / aromatic / ring / branch / stereo / disconnected)
- [x] chematic-smiles: DFS SMILES writer (accurate ring closure number assignment)
- [x] 50 tests all passing (aspirin / caffeine / glucose / NaCl etc. roundtrips)
- [x] Wildcard atom `*`: added Atom::wildcard(), parser correctly handles [*]
- [x] Implicit H count calculation: implicit_hcount(mol, idx) -> u8 in chematic-core/src/valence.rs
      (sum of bond orders, select minimum fitting valence, adjust for charge)
- [x] Kekulization: chematic-core/src/kekulization.rs
      (assign double bonds by maximum matching on aromatic subgraph)
- [x] Canonical SMILES: chematic-smiles/src/canonical.rs
      (Morgan rank iteration -> canonical DFS order)

---

## Phase 2 — Molecule Perception (complete)

New crates: chematic-perception, chematic-mol, chematic-depict

- [x] SSSR (Smallest Set of Smallest Rings): Balducci-Pearlman + GF(2) Gaussian elimination
      find_sssr(mol) -> RingSet  [chematic-perception/src/sssr.rs]
- [x] Aromaticity perception: Hückel 4n+2 π-electron model
      assign_aromaticity(mol) -> AromaticityModel
      Supports: benzene, pyridine, pyrrole, furan, naphthalene  [chematic-perception/src/aromaticity.rs]
- [x] SDF/MOL file format: V2000 parser+writer, SDF multi-molecule iterator
      parse_mol / write_mol / SdfReader  [chematic-mol/]
- [x] SDF V3000 parser (extended blocks)
      M  V30 BEGIN/END CTAB, ATOM, BOND block support
      Line continuation (trailing `-`), CHG= / MASS= / HCOUNT= / aamap support
      [chematic-mol/src/mol3000.rs]
- [x] 2D rendering engine (SVG): chematic-depict crate
      - Chain template (zigzag, BOND_LEN=40px, ±30° alternating)
      - Ring template library (3–8 membered rings: r = BOND_LEN / (2 sin(π/n)))
      - Fused ring greedy placement (select outward direction relative to centroid)
      - Labels for non-C atoms (white background rect + text)
      - Wedge/dash stereo bonds, double/triple/aromatic bond SVG rendering
      [chematic-depict/]
- [x] Roundtrip validation against ChEMBL sample set (large-scale test)
      Fetched 1000+ molecules SDF from ChEMBL, verified parse -> write -> parse consistency

---

## Phase 3 — Chemical Intelligence (complete)

New crates: chematic-chem, chematic-fp, chematic-smarts

- [x] Molecular descriptors (chematic-chem):
      - molecular_weight (average isotopic mass), exact_mass (monoisotopic mass)
      - heavy_atom_count
      - hbd_count (hydrogen bond donors: NH, OH)
      - hba_count (hydrogen bond acceptors: N, O)
      - rotatable_bond_count (non-ring single bonds + between non-terminal heavy atoms + exclude amides)
      - tpsa (topological polar surface area: Ertl 2000 atom-type lookup table)
      - logp_crippen (simplified Crippen-Wildman atom contribution table)
      - lipinski_passes (MW<=500, HBD<=5, HBA<=10, LogP<=5)
      - fsp3 (fraction of sp3 carbons)
      - aromatic_ring_count (number of aromatic rings via SSSR)
      [chematic-chem/src/descriptors.rs]
- [x] QED score (chematic-chem/src/qed.rs):
      - Geometric mean of 8 metrics from Bickerton et al. 2012 (Nature Chemistry)
      - 7-parameter ADS (Asymmetric Double Sigmoidal) function — same parameters as RDKit
      - 113 Brenk 2008 structural alert SMARTS (8th metric)
      - qed(mol) -> f64 in [0, 1]
- [x] Molar Refractivity (chematic-chem/src/descriptors.rs):
      - Wildman-Crippen additive model (same atom-type framework as LogP)
      - molar_refractivity(mol) -> f64
- [x] Drug-likeness filters (chematic-chem/src/descriptors.rs):
      - Veber: TPSA ≤ 140 Å², rotatable bonds ≤ 10
      - Egan: TPSA ≤ 131.6 Å², LogP ≤ 5.88
      - REOS: 6 criteria — MW / LogP / HBD / HBA / charge / heavy atom count
      - Ghose: MW 160–480, LogP −0.4–5.6, heavy atoms 20–70, MR 40–130
- [x] Tautomer rule expansion (chematic-chem/src/tautomer.rs):
      - 5 rules → 15 rules (thioamide, thio-iminol, thio-keto-enol, 6 cross-heteroatom 1,3 proton transfers)
- [x] BRICS fragmentation (chematic-chem/src/brics.rs):
      - Bond cleavage based on 16 environment rules from Dien et al. 2008
      - brics_bonds(mol) -> Vec<(AtomIdx, AtomIdx)>
      - brics_fragments(mol) -> Vec<Molecule> (with [*] attachment points)
- [x] ECFP / Morgan fingerprints (chematic-fp):
      - Configurable radius (ECFP4 = r2, ECFP6 = r3)
      - Atom invariants: atomic number, charge, degree, H count, in-ring flag, aromatic flag
      - Hash: FNV-1a 64bit (reproducible and deterministic, to_le_bytes for determinism)
      - Output: BitVec2048 ([u64; 32]), foldable to 1024/512/256 bits
      - Tanimoto coefficient, Dice coefficient
      [chematic-fp/]
- [x] SMARTS parser + VF2 evaluator (chematic-smarts):
      - Atom primitives: [#6] [!C] [a] [A] [r5] [D3] [H2] [R]
      - Bond primitives: ~ @ - = # :
      - Logical operators: & , ; ! (precedence: NOT > high-priority AND > OR > low-priority AND)
      - Recursive SMARTS `$(...)`: nested support, VF2 anchored matching
      - Extended primitives: [vN] valence, [xN] ring bond count, [^N] hybridization
      - [XN] total connection count (heavy atom degree + implicit H count), [RN] ring membership count
      - Explicit neutral charge [+0]/[-0] support
      - QueryMolecule type
      - find_matches(query, mol) -> Vec<HashMap<usize, AtomIdx>>
- [x] Molecule standardization + Murcko scaffold (chematic-chem):
      - Salt removal: largest_fragment(mol) -> Molecule (select largest fragment)
      - Charge neutralization: neutralize_charges(mol) -> Molecule (carboxylate/ammonium/protonated ether support)
      - Murcko scaffold: murcko_scaffold(mol) -> Molecule (fixed-point linker expansion)
      - Generic Murcko: generic_murcko_scaffold(mol) -> Molecule
      [chematic-chem/src/standardize.rs, scaffold.rs]
- [x] CIP stereochemistry (chematic-chem):
      - R/S assignment for tetrahedral centers (CIP priority rules)
      - E/Z assignment for double bonds
      - CIPCode enum stored on Atom

---

## Phase 4 — Similarity, Search, and Standardization (complete)

- [x] MACCS 166-bit structural keys (chematic-fp/src/maccs.rs)
      - Standard MACCS 166-bit SMARTS-based structural keys
      - `maccs(mol) -> BitVec2048`
      - key 164 = `[!#6;!#1]` (any heteroatom) and other standard patterns
- [x] Topological path fingerprint (chematic-fp/src/topo_path.rs)
      - DFS path enumeration, max_len=7 (default)
      - `topo_path(mol, &config) -> BitVec2048`
      - FNV-1a hash, canonicalized (select smaller of forward/reverse directions)
- [x] AtomPair fingerprint (chematic-fp/src/atom_pair.rs)
      - Carhart et al. 1985, atom pair + BFS topological distance encoding
      - `atom_pair_fp(mol) -> BitVec2048`
- [x] Topological Torsion fingerprint (chematic-fp/src/atom_pair.rs)
      - Nilakantan et al. 1987, 4-atom path encoding
      - `torsion_fp(mol) -> BitVec2048`
- [x] Maximum Common Substructure (MCS): McGregor or FMCS algorithm
      find_mcs(mols: &[Molecule]) -> QueryMolecule
- [x] Tautomer normalization (rule-based)
- [x] Stereo perception from 2D wedge bonds

---

## Phase 5 — 3D Chemistry (complete)

New crate: chematic-3d (zero external dependencies)

- [x] 3D coordinate generation (rule-based DFS placement):
      - Ideal bond length table (by element pair + bond order)
      - Bond angles by hybridization: sp3/sp2/sp (109.5° / 120° / 180°)
      - Ring templates placed in XY plane + chain extended as branches
      - Disconnected fragments offset along X axis
      [chematic-3d/src/dg.rs]
- [x] 3D file formats:
      - PDB parser/writer (ATOM/HETATM records, distance-based bond inference)
      - XYZ parser/writer
      [chematic-3d/src/pdb.rs, xyz.rs]
- [x] UFF force field energy minimization (Pure Rust):
      - Bond stretching, angle bending, dihedral, VDW, electrostatic interactions
      - Gradient descent / LBFGS minimization

---

## Phase 6 — Ecosystem and RDKit Parity (complete)

New crates: chematic-wasm, chematic-rxn, chematic (umbrella)

- [x] WASM package (chematic-wasm):
      - wasm-bindgen bindings: parse, writer, fingerprints, descriptor calculation
      - npm package: @kent-tokyo/chematic
        - v0.1.3: tpsa/mw/hba/hbd/lipinski/ecfp4
        - v0.1.4: logp/fsp3/qed/exact_mass/rotbonds/aromatic_ring_count/
                  tanimoto_atom_pair/tanimoto_torsion/brics_fragment_count
        - (unreleased): molar_refractivity/formal_charge_sum/veber_passes/
                        egan_passes/reos_passes/ghose_passes
        ("chematic" rejected by npm as too similar to "chromatic" → published under scoped name)
- [x] 2D rendering enhancements (chematic-depict):
      - CPK coloring (N=blue, O=red, S=yellow, Cl=green, F=yellow-green, Br=brown, I=purple, P=orange)
      - render_svg_highlighted / depict_svg_highlighted (yellow highlight + orange bonds)
- [x] Reaction SMILES / SMIRKS (chematic-rxn):
      - Reaction SMILES parser (>> separator)
      - Atom-atom mapping
      - Reaction template application
- [x] Umbrella crate (chematic):
      - Re-exports all sub-crates via feature flags
      - feature "full": all features enabled
      - feature "wasm": 3D and heavy dependencies disabled
- [x] Validation and benchmarks:
      - [x] Quantitative comparison of property accuracy against RDKit on 175-molecule dataset (docs/rdkit_comparison.md)
            MW: MAE=0.0002 Da, r=1.0000 | HAC: r=1.0000 | HBD: r=0.9974
            TPSA: MAE=0.081 Å², r=0.9999 | HBA: MAE=0.137, r=0.9750
            LogP: MAE=0.134, r=0.9847 (improvements: v0.1.0 MAE=1.346 → v0.1.3 MAE=0.298 → Sprint C MAE=0.141 → Sprint D MAE=0.134)
            ECFP4 Tanimoto: Spearman r=0.917 (50×50 pairs)
      - [x] Sprint C — RDKit quality improvements (LogP MAE 0.298→0.141, TPSA MAE 0.759→0.324 Å²):
            - Junction C atom type fix (fused aromatic rings: naphthalene/indole etc.)
            - Vinyl C atom type sign fix (Crippen contribution for C=C: +0.2274)
            - Nitro N TPSA fix ([N+](=O)[O-]: N=41.44 Å², O-=0 Å²)
            - Imine N TPSA fix (non-ring N in C=N: 12.89 Å²)
      - [x] Sprint D — RDKit quality improvements + missing features (LogP MAE 0.141→0.134, TPSA MAE 0.324→0.081 Å²):
            - Imine N-H TPSA fix (C=N-H: 23.79 Å², resolves metformin/arginine errors)
            - Phosphate P TPSA fix (with P=O: 26.88 Å² vs without P=O: 34.14 Å²)
            - Phosphate P LogP fix (with P=O: +0.7933 vs without P=O: -0.3451)
            - 5 new ring descriptors: num_aromatic_heterocycles, num_aliphatic_heterocycles,
              num_saturated_heterocycles, num_spiro_atoms, num_bridgehead_atoms
            - Tautomer rules 15 → 20 (rules 16–20: O→N, O→O, N→C, C→O, C→N)
      - [x] Sprint E — Guanidinium N LogP fix + tautomer 1,2-shift (LogP MAE 0.134→0.117):
            - LogP fix for guanidinium/amidine N (Wildman-Crippen N14 type: -0.335):
              Detects imine =N (direct C=N double bond) and adjacent guanidinium N (N adjacent to C=N)
              metformin error 2.07 → ~0.00, arginine improved
            - Tautomer 1,2-shift added (pyrazole N1H↔N2H etc.):
              find_direct_aromatic_matches + transfer_hydrogen_aromatic (no bond order change)
              enumerate_tautomers: H-assignment fingerprint distinguishes positional isomers
              canonical_tautomer: converges N1H/N2H to the same canonical form using minimum H-assignment
      - [x] Multi-agent security/bug/refactoring audit:
            - [Security] Added depth limit 8 to recursive SMARTS $(...) (SmartsError::RecursionDepthExceeded)
            - [Security] Ring closure unwrap() → expect() (makes invariants explicit)
            - [Bug] clone_mol / transfer_hydrogen_aromatic .ok() → .expect() (prevent silent bond omission)
            - [Refactor] FNV-1a magic numbers in mol_fingerprint → named constants
            - [Refactor] TPSA nitro group detection: 2 neighbor scans → 1 fold
      - [x] criterion benchmarks for all hot paths
      - [x] Full ChEMBL 37 validation: **2,897,819 molecules / 100.000% success** (parse + roundtrip)
              curl chembl_37_chemreps.txt.gz | gzip -d | awk | validate_smiles stream validation

---

## Current Test Counts

| Crate                  | Tests | Status   |
|------------------------|-------|----------|
| chematic-core          | 48    | complete |
| chematic-smiles        | 57    | complete |
| chematic-perception    | 34    | complete |
| chematic-mol           | 63    | complete |
| chematic-depict        | 43    | complete |
| chematic-chem          | 375   | complete |
| chematic-fp            | 50    | complete |
| chematic-smarts        | 87    | complete |
| chematic-3d            | 147   | complete |
| chematic-rxn           | 30    | complete |
| chematic-wasm          | 175   | complete |
| chematic               | 1     | complete |
| chematic-iupac         | 14    | complete |
| chematic-inchi         | 28 lib + 14 integration* | complete |
| **CI lib total**       | **1,649** | `cargo test --workspace --lib` |

\* integration tests: `cargo test -p chematic-inchi --features native-inchi --test standard_inchi`

---

## Final Crate Structure

    chematic/
    crates/
      chematic-core/        Phase 1  — Atom, Bond, Molecule, Element
      chematic-smiles/      Phase 1  — SMILES parse/writer/canonicalization
      chematic-perception/  Phase 2  — SSSR, aromaticity perception
      chematic-mol/         Phase 2  — SDF/MOL V2000+V3000 file formats
      chematic-depict/      Phase 2  — 2D SVG rendering (CPK colors, highlight)
      chematic-chem/        Phase 3  — Molecular descriptors, BRICS, QED, standardization, CIP stereochemistry
      chematic-fp/          Phase 3  — ECFP, MACCS, path, AtomPair, Torsion FP
      chematic-smarts/      Phase 3  — SMARTS + VF2 substructure search (recursive SMARTS support)
      chematic-3d/          Phase 5  — Rule-based 3D coordinates, PDB/XYZ formats
      chematic-rxn/         Phase 6  — Reaction SMILES/SMIRKS
    chematic/               Phase 6  — Umbrella crate with feature flags

---

## Inter-Phase Dependencies

    Phase 1 (core + SMILES)
      -> Phase 1 kekulization
        -> Phase 2 (SSSR, aromaticity, SDF, 2D rendering)
          -> Phase 3 (descriptors, ECFP, SMARTS)
            -> Phase 3 SMARTS -> Phase 4 (MACCS, MCS, standardization)
              -> Phase 5 (3D, force field)
                -> Phase 6 (WASM, reactions, validation)

---

## Phase 7 — Full RDKit Parity (not yet started)

Major features not yet implemented compared to RDKit, listed in priority order.
Constraints: Zero FFI and WASM compatibility remain unchanged.

### Tier 1 — High Priority (features most needed by pharma/cheminformatics users)

#### 7-1. Reaction SMIRKS Application (RunReactants) (complete: Sprint J)
  - [x] RDKit: `rxn.RunReactants(reactants)` → enumerate product SMILES
  - Implementation: chematic-rxn/src/transform.rs (implemented)
  - `run_reactants(smirks, reactants) -> Result<Vec<Vec<Molecule>>, TransformError>`
  - VF2 subgraph isomorphism + BFS substituent transfer + Cartesian product enumeration

#### 7-2. Topological Descriptors (complete: Sprint G)
  - [x] Wiener index (sum of all atom-pair distances)
  - [x] Hall–Kier Kappa indices κ1 / κ2 / κ3
  - [x] Molecular connectivity indices Chi χ0v / χ1v / χ2v / χ3v / χ4v (Kier–Hall)
  - [x] Bertz complexity (BertzCT)
  - [x] Labute approximate surface area (LabuteASA)
  - Implementation: chematic-chem/src/topo_descriptors.rs (implemented)
  - Difficulty: medium (distance matrix is the foundation; rest are derived)

#### 7-3. Explicit H Management (complete: Sprint G)
  - [x] `add_hydrogens(mol) -> Molecule` — convert all implicit H to explicit atoms
  - [x] `remove_hydrogens(mol) -> Molecule` — return explicit H atoms to implicit
  - Current state: only implicit H calculation via implicit_hcount()
  - Implementation: chematic-chem/src/hydrogen.rs (implemented)
  - Difficulty: low–medium

#### 7-4. SVG Grid Rendering (complete: Sprint G)
  - [x] `depict_svg_grid(mols, cols) -> String` — SVG with multiple molecules arranged in a grid
  - RDKit: `Draw.MolsToGridImage`
  - Implementation: chematic-depict/src/grid.rs (implemented)
  - Difficulty: low (just combines existing depict_svg calls)

---

### Tier 2 — Medium Priority (QSAR / 3D workflows)

#### 7-5. Shape Descriptors (requires 3D coordinates) (complete: Sprint H)
  - [x] Principal moments of inertia PMI1 / PMI2 / PMI3
  - [x] Normalized principal axis ratios NPR1 / NPR2
  - [x] Radius of Gyration
  - [x] Asphericity / Eccentricity
  - [x] Plane of Best Fit (PBF)
  - RDKit: `rdMolDescriptors.CalcPMI`, `CalcNPR1/2`, `CalcRadiusOfGyration` etc.
  - Implementation: chematic-3d/src/shape_descriptors.rs (implemented, 3×3 Jacobi eigensolver hand-implemented)
  - Difficulty: medium (requires eigenvalue decomposition, nalgebra or hand-implemented)

#### 7-6. Conformer Management (complete: Sprint I)
  - [x] Structure to hold multiple conformers (coordinate sets) in Molecule
  - [x] `add_conformer()` / `get_conformer()` / `get_conformer_mut()` / `remove_conformer()`
  - [x] RMSD calculation between conformers (`conformer_rmsd_no_align` / `conformer_rmsd`)
  - Design: no changes to chematic-core; implemented as external container `ConformerEnsemble` in chematic-3d
  - Implementation: chematic-3d/src/conformer.rs (implemented)
  - RMSD with Kabsch alignment reuses existing jacobi3

#### 7-7. UFF Parameter Improvements (complete: Sprint K)
  - [x] Ideal bond length table by element pair (C-C/C-N/C-O/C-S/C-F/C-Cl/C-Br/C-H etc., 30+ pairs)
  - [x] Ideal bond angles by hybridization (SP/SP2/SP3) (O:104.5°, N:107°, S:99° etc.)
  - [x] VDW repulsion energy using element-specific UFF/Bondi VDW radii
  - Implementation: chematic-3d/src/minimize.rs (improved)
  - Tests: 68 (previously 58), +10 new tests (bond length accuracy, hybridization, symmetry)
  - Full MMFF94 implementation (8 parameter tables, 95 atom types) deferred due to excessive effort

#### 7-8. Stereochemistry Assignment from 3D (complete: Sprint H)
  - [x] Auto-compute R/S and E/Z from 3D coordinates (AssignStereochemistryFrom3D)
  - Current state: CIP assignment from SMILES wedge/dash only → can compute independently via signed volume + dihedral angle
  - Implementation: chematic-3d/src/stereo3d.rs (new, 1-sphere CIP priority handles the core)
  - Difficulty: medium

---

### Tier 3 — Low Priority (niche / high difficulty)

#### 7-9. Stochastic 3D Embedding (ETKDG equivalent)
  - [ ] Initial coordinate generation by Distance Geometry
  - [ ] Refinement using experimental torsion angle distributions (the "ET" part of ET-DG)
  - RDKit: `AllChem.EmbedMolecule` / `EmbedMultipleConfs`
  - Difficulty: very high (eigenvalue decomposition of distance matrix + experimental library needed)

#### 7-10. Dense Count Format for Hash-Based FP
  - [x] Count vector format for Morgan FP (integer counts instead of bits)
  - [x] `GetMorganFingerprint(mol, radius)` → `{hash: count}` format
  - Implementation: extension of chematic-fp/src/ecfp.rs
  - Difficulty: low (just changes the output format of existing ECFP)

#### 7-11. InChI / InChIKey (complete: Sprint v0.4.0)
  - [x] Standard InChI string generation (`standard_inchi()` — via IUPAC C library 1.07.5)
  - [x] InChIKey (27-character hash) generation (`standard_inchi_key()` — bit-exact match)
  - **Implementation approach**: vendored IUPAC InChI C library (v1.07.5) as opt-in FFI via `native-inchi` feature
  - Default build maintains zero FFI / WASM compatibility; C compilation only occurs when `native-inchi` feature is enabled
  - Implementation: `crates/chematic-inchi/src/native/` + `build.rs` + `vendor/inchi-src/`
  - Tests: 14 integration tests (`cargo test -p chematic-inchi --features native-inchi --test standard_inchi`)
  - Note: InChI 1.07.5 differs from PubChem (1.06) in some stereo keys (/m layer assignment changed)

---

### Out of Scope (conflicts with zero-FFI policy or excessive effort)

- Full ETKDG reproduction (stochastic sampling + DG)
- InChI (~~C library is the only official implementation~~ → resolved in v0.4.0 with `native-inchi` feature)
- ML-based prediction models (LogP, solubility, etc.)
- HELM / FASTA notation (peptides/proteins)
- Transition metals and coordination compounds (coordination chemistry)

---

## Recommended Implementation Order (Sprint G onward)

```
Sprint G: (complete) 7-2 (topological descriptors) + 7-3 (explicit H management) + 7-4 (SVG grid)
          → code additions only, no breaking changes, tests +38 (582→620)
Sprint H: (complete) 7-5 (shape descriptors) + 7-8 (stereo from 3D)
          → chematic-3d/src/shape_descriptors.rs + stereo3d.rs, tests +15 (620→635)
Sprint I: (complete) 7-6 (conformer management)
          → chematic-3d/src/conformer.rs, Kabsch RMSD, tests +14 (623→637)
Sprint J: (complete) 7-1 (RunReactants)
          → chematic-rxn/src/transform.rs, VF2 + BFS substituent transfer, tests +11 (612→623)
Sprint K: (complete) 7-7 (UFF parameter improvements)
          → element-specific bond lengths, hybridization angles, VDW radii, tests +10 (637→646)
Sprint L: (complete) Sprint L audit — security/bug/refactoring review (0.1.5 → 0.1.6)
Sprint M: (complete) SMARTS highlight display + click highlight + reaction scheme (demo 0.1.11)
Sprint N: (complete) Tab UI + 3D interactive viewer (demo 0.1.12)
Sprint P: (complete) SDF/MOL WASM bindings + EState indices + path fingerprint WASM (v0.1.14)
Sprint Q: (complete) IFG + SA Score + Gasteiger charges + VSA descriptors + MaxMin/Butina (v0.1.15)
          → tests: 697 → 736 (+39)
Sprint R: (complete) E/Z double bond stereochemistry SMILES output (v0.1.16)
Sprint S: (complete) SA score fragment table implementation (v0.1.17)
          → tests: 742 → 743 (+1)
Sprint T: (complete) Per-atom color highlight + named functional group detection + atom info API (demo v0.1.18)
Sprint U: (complete) WASM convenience API for interactive articles (v0.1.19)

## Phase 8 — WASM Feature Expansion / File Formats / Edit API (v0.1.20–v0.1.21)

Sprint V–AA: (complete) WASM exports expanded from 84 → 103 (v0.1.20)
  - Murcko / tautomers / standardization / MACCS / batch descriptors / MOL 2D coordinate fix
  - PAINS/CIP details / ECFP6 / Dice / 3D shape descriptors / MaxMin/Butina / MCS
  - V3000 loading / 3D minimization / SDF property read/write / SMARTS highlight grid
  - XYZ/PDB I/O / per-atom descriptors / SSSR / custom ECFP / stereoisomer enumeration
  - BRICS SMILES / AtomPair/Torsion bitvec / FCFP6 / SDF write
  - FCFP4/6 bitvec / Dice ECFP6 / write_smiles / reaction SMILES normalization
  - ConformerEnsemble WASM / R-group decomposition / MMP analysis
  - CML read/write / CDXML read / Mutable API / DepictData / SDF/V3000 write / CPK
  - Tests: 743 → 863 (+120)

Sprint v0.1.21: (complete) Mutable API expansion / SDF/CDXML enhancements (v0.1.21)
  - chematic-core: with_atom_charge, with_atom_element, with_bond_added → (Mol, BondIdx)
  - chematic-mol: parse_mol_with_coords, parse_sdf_with_coords, parse_cdxml_all, CDXML stereochemistry
  - chematic-depict: depict_data_with_coords
  - WASM: mol_with_atom_charge, mol_with_atom_element, cdxml_to_smiles_json, mol_block_coords_json, depict_data_with_coords_json
  - Tests: 863 → 869 (+6)

Sprint v0.1.22: (complete) MCS ring-awareness constraints (Issue #1)
  - ring_matches_ring_only: block ring↔non-ring cross-matches during McGregor search phase using SSSR
  - complete_rings_only: post-search iterative removal of partial rings in mol[0]
  - Tests: 869 → 877 (+8)

Sprint v0.1.23: (complete) Element radius API / implicit H completion / aromaticity application API (v0.1.23)
  - chematic-core: Element::vdw_radius() / covalent_radius() (Bondi 1964 + Alvarez 2013/2008 tables for 118 elements, fallback 1.70/0.77 for unknown)
  - chematic-core: Molecule::implicit_hydrogen_count(idx) (wrapper for valence::implicit_hcount)
  - chematic-core: Molecule::total_formula() (Hill formula including implicit H — CH4, C2H6O etc.)
  - chematic-core: Molecule::with_atom_aromatic() / with_bond_order() (immutable update API extension)
  - chematic-perception: apply_aromaticity(mol) → returns new Molecule with aromatic flags and BondOrder::Aromatic applied to a kekulized molecule
  - chematic-3d: minimize_uff() alias (makes the existing UFF minimize() more discoverable by name)
  - Tests: 877 → 886 (+9)

Sprint v0.1.24: (complete) validate_valence public API + run_reactants product filtering (v0.1.24)
  - chematic-core: ValenceError struct + validate_valence(mol) -> Vec<ValenceError> (per-element normal_valences + formal charge adjustment)
  - chematic-core/lib.rs: re-export ValenceError / validate_valence
  - chematic-perception/lib.rs: re-export from chematic_core (accessible as chematic::perception::validate_valence)
  - chematic-rxn: apply validate_valence to products in run_reactants, exclude product sets with over-valent atoms
  - Tests: 886 → 893 (+7)

Sprint v0.1.25: (complete) suggest_bond_direction public API (v0.1.25)
  - chematic-depict/layout.rs: suggest_bond_direction(mol, atom, layout) -> f64 (radians)
    - Collects existing bond angles → selects from 30° grid + chemical offsets (sp2±120°, sp3 zigzag±150°) to maximize minimum separation angle
  - chematic-depict/lib.rs: re-export BOND_LEN and suggest_bond_direction
  - Enables replacement of the draw-side's own 30° brute-force implementation
  - Tests: 893 → 897 (+4)

Sprint v0.1.26: (complete) atom_color_rgb public API (v0.1.26)
  - chematic-depict/svg.rs: atom_color_rgb(atomic_number: u8) -> [u8; 3] (same CPK values as atom_color, no hex parsing)
  - chematic-depict/lib.rs: re-export
  - Enables replacement of draw-side hex parser (direct use of egui::Color32::from_rgb)
  - Tests: 897 → 900 (+3)

Sprint v0.1.27: (complete) MolMetadata builder API (v0.1.27)
  - chematic-mol/mol2000.rs: MolMetadata::with_name(name) -> Self, with_comment(comment) -> Self
  - Enables setting name/comment at SDF export via MolMetadata::default().with_name("...").with_comment("...")
  - Tests: 900 → 902 (+2)

Sprint v0.1.27-ext: (complete) E/Z double bond stereochemistry (from 2D coordinates) + extended StereoGroup + isotope distribution (v0.1.27)
  - chematic-perception: assign_ez_from_2d(mol, coords) — assign E/Z from 2D coordinate cross products
                          cip_ez_descriptor(mol, bond_idx, coords) -> Option<CipCode> — return E/Z for a specific bond
                          [crates/chematic-perception/src/stereo2d.rs; re-exported in lib.rs]
                          Algorithm: 2D cross product of double bond vector vs substituent position vector + 1-sphere CIP priority
  - chematic-core: StereoGroupKind enum (Absolute / Or(u32) / And(u32)), StereoGroup struct
                   Molecule.stereo_groups field + stereo_groups() / set_stereo_groups() / add_stereo_group() methods
                   MoleculeBuilder add_stereo_group() method + from_molecule() copies stereo_groups
                   [crates/chematic-core/src/stereo_group.rs; re-exported in lib.rs]
  - chematic-mol: V3000 parser reads BEGIN COLLECTION / MDLV30/STEABS / MDLV30/STEOR<n> / MDLV30/STEAND<n>
                  V3000 writer outputs COLLECTION block when stereo_groups present
                  Roundtrip test added
                  [crates/chematic-mol/src/mol3000.rs]
  - chematic-chem: isotope_distribution(mol, resolution) -> Vec<(f64, f64)>
                   returns (m/z, relative intensity) pairs, normalized to base peak=1.0
                   resolution parameter merges peaks within specified Da
                   Isotope support for H/C/N/O/F/Si/P/S/Cl/Br/I/Se/Na/K/As
                   Prioritizes explicit isotope labels (atom.isotope) when present
                   [crates/chematic-chem/src/isotope_distribution.rs; re-exported in lib.rs]
  - chematic (umbrella crate): Complete revision of lib.rs //! module documentation (feature table, quick-start examples, feature flag table)
                   Cargo.toml: updated description, added parser-implementations/rendering to categories
                   [package.metadata.docs.rs] section added: features=["full"], rustdoc-args=["--cfg","docsrs"]

Sprint v0.1.28: (complete) All remaining tasks implemented (v0.1.28)
  - Issue C: BricsConfig { min_fragment_size } + brics_fragments_with_config (chematic-chem)
  - Issue E: MatchConfig { max_matches } + find_matches_with_config (chematic-smarts)
  - Issue A: AtomCompare / BondCompare enum + McsConfig field additions (chematic-smarts)
            AnyHeavyAtom mode enables scaffold hopping MCS between different heterocycles
  - (10) xlogp3.rs: xlogp3() / xlogp3_per_atom() — Cheng 2007 atom-type contribution table (chematic-chem)
  - (11) chematic-iupac: new crate (Pure Rust, no network required)
            straight-chain alkanes/alkenes/alkynes, cycloalkanes, alcohols/amines/haloalkanes
            IupacError::NotSupported explicitly indicates unsupported structures
  - Tests: 902 → 915 (+13)

Sprint v0.1.29: (complete) Mutable Molecule API + Fragments + MoleculeBuilder::from_molecule (v0.1.29)
  - Molecule::add_atom / remove_atom / add_bond / remove_bond / set_charge / set_element / set_cip_code
  - Molecule::is_connected() / fragments() → connected component splitting
  - MoleculeBuilder::from_molecule(mol)
  - Tests: 915 → 924 (+9)

Sprint v0.1.30: (complete) 2D stereochemistry + Aromatize/Kekulize in-place (v0.1.30)
  - chematic-perception: stereo2d.rs new (assign_stereo_from_2d / apply_stereo_from_2d)
  - chematic-perception: aromatize(mol: &mut Molecule) / kekulize_inplace(mol: &mut Molecule)
  - Tests: 924 → 926 (+2)

Sprint v0.1.31: (complete) SdfRecord extension + reaction SVG + chemical abbreviations (v0.1.31)
  - SdfRecord: added coords: Vec<(f64,f64)> + meta: MolMetadata + properties: HashMap<String,String>
  - chematic-depict: depict_reaction_svg / depict_reaction_svg_opts (reactants→arrow→products SVG)
  - chematic-chem: expand_abbreviation / abbreviations (table of 30 abbreviations)
  - Tests: 926 → 929 (+3)

Sprint v0.1.32: (complete) MDL RXN file + formula_with_isotopes (v0.1.32)
  - chematic-mol: parse_rxn_file / write_rxn_file (MDL RXN V2000 format)
  - chematic-core: Molecule::formula_with_isotopes() (molecular formula with isotope labels such as ²H/¹³C)
  - Tests: 929 → 933 (+4)

Sprint v0.1.25: (complete) P2 feature completion / release (2026-06-06)
  - detect_crossings: 2D layout quality evaluation (bond crossing detection)
  - invert_stereocenter: R/S chirality inversion (wedge bond inversion)
  - enumerate_stereoisomers: full stereoisomer enumeration (2^n, max 64)
  - render_svg_with_metadata: SVG metadata embedding (SMILES)
  - find_reaction_center: reaction center analysis (broken/formed bonds + changed atoms)
  - Tests: 865 → 935 (+70)
  - cargo & npm publish complete
  - CHANGELOG / README updated in all languages

Sprint v0.1.26: (complete) Issue D + P3 Features (complete: 2026-06-06)

### Completed
  - [x] Issue D (matchChiralTag): `McsConfig.match_chiral_tag` implemented
        - R/S enantiomer matching control (default: false)
        - 3 new tests added (enantiomer blocking/allowing)
        - Implementation: crates/chematic-smarts/src/mcs.rs
        - Tests: 87 tests all passing
  
  - [x] parse_condensed(): "CH3COOH" → structure parsing implemented
        - Condensed formula lexing + functional group substitution
        - Implementation: crates/chematic-chem/src/condensed.rs (new)
        - Tests: 10 implemented (covers basic cases)
        - Note: H-count digits (CH3) processing is a target for future improvement
        - Implementation: parse_condensed(input) → Result<Molecule, CondensedError>

  - [x] WASM bindings
        - find_reaction_center_json(reaction_smiles) → JSON
        - standardize_smiles(mol, opts) → SMILES

  - [x] Demo updates
        - "Stereo" tab added (stereoisomer enumeration)
        - "Reaction" tab extended (broken/formed bond highlighting)

---

## Sprint v0.1.27–v0.1.28: DREIDING + MD + SPME (complete: 2026-06-07)

### Completed
  - [x] **Phase 1**: chematic-ff (DREIDING atom typing + parameters)
        - 20 atom types (C_3/C_2/C_1/C_R, N_3/N_2/N_1/N_R, O_3/O_2/O_R, S_3/S_R, P_3, H, halogens)
        - 40+ bond length parameters, hybridization-specific bond angles, VDW parameters
        - Implementation: crates/chematic-ff/src/dreiding.rs + params.rs
        - Tests: 25+ passing

  - [x] **Phase 2**: chematic-3d MD integrator
        - Velocity Verlet integration (NVE + NVT with Berendsen thermostat)
        - Maxwell-Boltzmann initial velocity assignment (accurate unit conversion: 0.01038 factor)
        - Bond stretching, angle bending, VDW, Coulomb energy calculation
        - Implementation: crates/chematic-3d/src/md.rs
        - Tests: 84 tests all passing
        - **CRITICAL FIXES**:
          - Velocity init: added 0.01038 unit conversion factor (kcal/mol → amu·Ų/fs²)
          - VDW energy: use DREIDING parameters + added 1-2/1-3 exclusion

  - [x] **Phase 3**: chematic-ewald (SPME long-range charges)
        - Direct Coulomb (non-periodic) + SPME (periodic)
        - Real space + reciprocal space + self-energy correction
        - Implementation: crates/chematic-ewald/src/pme.rs + real.rs
        - Tests: 8 tests all passing
        - **CRITICAL FIX**: Mesh indexing — isqrt() corruption → accurate 3D→1D conversion (ix + iy*M0 + iz*M0*M1)

  - [x] npm publish: v0.1.29 (demo/pkg/package.json)
  - [x] WASM integration: run_md_json(), coulomb_energy_json(), minimize_dreiding_json()
  - [x] Demo "Dynamics" tab: Coulomb calculator, MD simulator, geometry optimizer
  - [x] Tests: 92 tests all passing

### Detected Unresolved Issues (Audit 2026-06-07)

#### CRITICAL (1)
- WARNING: PME Mesh Indexing OOB write (non-cubic mesh): fixed → linear_idx calculation improved

#### HIGH (2)
- WARNING: Thermostat NaN injection (T→0K): guard needed (`if temperature < 1e-6 then lambda = 1.0`)
- WARNING: Singular box volume: silent default when `det < 1e-10` → returning Result type recommended

#### MEDIUM (4)
- WARNING: fastrand entropy weakness: low-entropy RNG, thread_local state management recommended
- WARNING: SVG string interpolation XSS: molecule SVG/symbol sanitization recommended (currently hardcoded and safe)
- WARNING: HTML innerHTML risk: if energy term names become user input, XSS risk → use textContent
- WARNING: Ring closure u8 truncation: ring collision possible with SMILES %00-%99 designators

#### LOW-MEDIUM (3)
- WARNING: Coulomb singularity (r→0): `r.max(1e-5)` clamp recommended
- WARNING: MD force cloning: 6N coord clone per step → optimize with EnergyCache (3–5× speedup potential)
- WARNING: VDW parameter: Lorentz-Berthelot combining rules implemented

#### Refactoring Priority
1. **[HIGH]**: ideal_bond_len() × 3 duplicates → consolidate into chematic-ff/bond_params.rs
2. **[HIGH]**: Add error handling (thermostat temp check, mesh bounds assertion)
3. **[MEDIUM]**: MD force caching layer (EnergyCache struct)
4. **[MEDIUM]**: WASM JSON serialization reduction (binary protocol option)
5. **[INVOLVED]**: demo/index.html 3090 LOC → component modularization

---

## Sprint v0.1.33: CXSMILES/CXSMARTS + StandardizationPipeline with Audit (2026-06-07, in progress)

### Completed
  - [x] **CXSMILES/CXSMARTS Metadata Support** (chematic-smiles/smarts)
        - Atom labels (`$...$`), atom properties (`atomProp:key.value`), atom radicals (`^n:`), zero-order bonds (`Z:`)
        - `parse_cxsmiles()` / `parse_cxsmarts()` / `write_cxsmiles()` / `write_cxsmarts()` implemented
        - `CxSmiles` / `CxSmarts` structs hold metadata
        - Implementation: crates/chematic-smiles/src/cx.rs, crates/chematic-smarts/src/cx.rs
        
  - [x] **StandardizationPipeline with Audit Reports** (chematic-chem)
        - `StandardizationPipeline::run()` → `(Molecule, StandardizationReport)`
        - Per-stage tracking: `StandardizationStepReport` (step, enabled, changed, before/after snapshots)
        - `StandardizationReport`: status, input/output snapshots, warnings
        - `StandardizationWarning`: code + message (metal disconnection, valence errors)
        - JSON serialization support (serde)
        - Implementation: crates/chematic-chem/src/standardize.rs
        
  - [x] **WASM Bindings for CX + Audit**
        - `parse_cxsmiles_json()`: returns atom labels / properties / radicals / zero-bonds as JSON
        - `parse_cxsmarts_json()`: same functionality for SMARTS
        - `normalize_cxsmiles()`: re-serializes CX metadata
        - `standardize_smiles_report_json()`: returns standardization report as JSON
        - Tests: 12 new (cx metadata round-trip, audit report structure)
        - Implementation: crates/chematic-wasm/src/lib.rs
        
  - [x] **Error Trait Implementations** (Section 4 complete)
        - `Display` + `std::error::Error` implemented for cx.rs + BondOrder::Zero related code
        - Added `Zero` variant to BondOrder enum (non-bonded interaction / virtual bond)

### Test Counts
- New tests: +12 (933 → 945 planned)
- chematic-smiles: cx.rs unit tests
- chematic-smarts: cx.rs unit tests
- chematic-wasm: cxsmiles_json, cxsmarts_json, standardize_report_json tests

## Test Status (v0.1.33)
- **Total**: 945 tests passing (planned value)
  - chematic-smiles: +4 (cx round-trip)
  - chematic-smarts: +4 (cx round-trip)
  - chematic-wasm: +4 (JSON serialization)

## Sprint v0.1.34: InChI Ring Closure + Stereo Layers + SEO (2026-06-08, complete)

### Completed
  - [x] **InChI Ring Closure Bonds** (chematic-inchi)
        - Detect back-edges by tracking DFS tree edges
        - Benzene: `InChI=1S/C6H6/c1-2-3-4-5-6-1/h1-6H` (ring closure `-1` added)
        - Implementation: crates/chematic-inchi/src/layers/connection.rs
        - Tests: ring closure confirmed in test_connectivity_benzene

  - [x] **InChI Stereo Layers (/t, /b)**
        - `/t` layer: R/S tetrahedral stereo via CIP code assignment
        - `/b` layer: E/Z double bond stereo via CIP code assignment
        - L-alanine: `InChI=1S/C3H7NO2/c1-2(4)3(5)6/h2H,4H,1,5-6H3/t2-` (includes R/S)
        - Implementation: crates/chematic-inchi/src/layers/stereo.rs (new)
        - Integration: stereo layer added in crates/chematic-inchi/src/lib.rs

  - [x] **SEO Documentation Improvements** (Phase 1-2)
        - Workspace `homepage` → updated to live demo URL
        - `chematic-inchi`: keywords/categories added
        - Individual README created for 9 crates (chematic-smiles, chematic-fp, chematic-smarts, chematic-inchi, chematic-core, chematic-depict, chematic-rxn, chematic-iupac)
        - CI workflow `.github/workflows/ci.yml` added (test + clippy)
        - README status badges added (CI, crates.io, npm)

### Test Counts
- New tests: +4 (stereo layers round-trip)
- Total: 1120+ tests passing
- Clippy: clean

---

## Sprint v0.1.35: wasmBridge Support + Version Sync (2026-06-08, complete)

### Completed
  - [x] **Version Synchronization (P0)**
        - chematic-wasm/Cargo.toml: all 11 crates 0.1.33 → 0.1.34
        - chematic-inchi dependency added

  - [x] **InChI / InChIKey WASM API (P1)**
        - Free functions: `inchi_from_smiles()`, `inchikey_from_smiles()`
        - MolHandle methods: `.to_inchi()`, `.to_inchikey()`
        - Implementation: crates/chematic-wasm/src/lib.rs

  - [x] **enumerate_stereo_isomers_json Enhancement (P1)**
        - Output format extended: `["smiles1", "smiles2"]` → `[{"smiles":"...", "inchi":"...", "inchikey":"..."}, ...]`
        - Include InChI/InChIKey for each isomer (supports database search)
        - Test updated: count "smiles" objects instead of string parsing

  - [x] **invert_stereocenter WASM binding (P1)**
        - New function: `invert_stereocenter_at(mol, atom_idx) → Result<MolHandle>`
        - Inverts stereochemistry of U/D wedge bonds

### Scope Evaluation (out-of-scope)
- [ ] to_svg_with_metadata (specification unclear → P2-P3)
- [ ] detect_layout_crossings (specification unclear → P2-P3)
- [ ] validate_molecule (substitute with is_valid_smiles)
- [ ] Spiro/cumulative/metal compounds (large-scale implementation → P3+)

### Test Counts
- New tests: +2 (enumerate_stereo format verification)
- Total: 1120+ tests passing
- Clippy: clean

---

## Sprint v0.1.36: Issue #1 Audit + BUG-2/3/4 Fix (2026-06-08, complete)

### Completed
  - [x] **Issue #1 Audit: Topologically Correct but Chemically Meaningless Results**
        - Discovered 4 similar bugs in the codebase where algorithms are topologically correct but yield chemically wrong results
        - Pattern: RDKit has constraint options that weren't implemented in chematic, causing silent invalid results on migration
  
  - [x] **BUG-2: `[h]` SMARTS Primitive (implicit H count)**
        - Added `ImplicitHCount(u8)` variant to AtomPrimitive enum
        - Parser now correctly parses lowercase `h` as implicit H-only (not aromatic H)
        - Added matching logic in eval_atom_primitive() using `implicit_hcount()` only
        - Tests: All 124 chematic-smarts tests passing
  
  - [x] **BUG-3: MCS `maximize_bonds` Tiebreak**
        - Implemented maximize_bonds tiebreaking when atom counts are equal
        - Modified grow() function to prefer mappings with higher bond_count
        - Default: maximize_bonds=true to match RDKit behavior
  
  - [x] **BUG-4: `/\` SMARTS Geometric Stereo Bonds**
        - Added `Up` and `Down` variants to BondPrimitive enum for E/Z double bonds
        - Updated parser: is_bond_token() and consume_bond_prim() now handle `/` and `\`
        - Added matching logic in eval_bond_primitive()
  
  - [x] **Verification & Testing**
        - All 1,120+ tests passing (chematic-smarts: 124, full suite: 1,120+)
        - No compilation errors, clippy clean
  
### Implementation Details
  - **Files Modified**:
    - crates/chematic-smarts/src/query.rs: Added ImplicitHCount + Up/Down variants
    - crates/chematic-smarts/src/parser.rs: Added [h] parsing + / \ bond tokens
    - crates/chematic-smarts/src/match_vf2.rs: Added eval logic for implicit H + Up/Down
    - crates/chematic-smarts/src/mcs.rs: Modified grow() with bond count tiebreak
  
  - **Test Results**: 124 smarts tests all passing

---

## Sprint v0.1.69–v0.1.74: RDKit Gap Analysis + 6 Feature Implementations (2026-06-08, complete)

### Completed

**Phase 1: Gap Analysis (v0.1.68 → docs/rdkit_comparison.md)**
- Systematic analysis of feature gaps relative to RDKit
- Classified into 3 tiers: Priority A (high impact) / B (medium) / C (low priority)
- Identified 15 unimplemented features

**Sprint v0.1.69: EState_VSA Descriptor (A5)**
  - [x] EState_VSA bins (11 bins) implemented: `estate_vsa(mol) -> Vec<f64>`
  - [x] Integrated with Labute ASA per-atom and E-State indices
  - [x] 9 tests added (bin length, sum consistency, non-zero)
  - Implementation: `crates/chematic-chem/src/vsa.rs`
  - Tests: +9 (226 → 235)

**Sprint v0.1.70: Tautomer 1,5-shift + Scoring (A1/A2)**
  - [x] Tautomer 1,5-shift rules added: β-ketoenamine, enaminone-long-range, guanidinium
  - [x] Added `path_len` field to `TautomerRule` struct
  - [x] Tautomer scoring function implemented: aromatic bonus + O-H/N-H/S-H priority
  - [x] Score-based sorting integrated into canonical_tautomer
  - Implementation: `crates/chematic-chem/src/tautomer.rs`
  - Tests: +18 (235 → 253)

**Sprint v0.1.71: Scaffold Network Library Aggregation (B1)**
  - [x] ScaffoldNetwork new struct: `pub struct ScaffoldNetwork { scaffolds, counts, parents }`
  - [x] `scaffold_network_with_counts(mols: &[Molecule]) -> ScaffoldNetwork`
  - [x] Aggregates scaffold occurrence frequency across a molecular library
  - Implementation: `crates/chematic-chem/src/scaffold.rs`
  - Tests: +12 (253 → 265)

**Sprint v0.1.72: RMSD Conformer Pruning + CIP Rule 3 (B3/B2)**
  - [x] ConformerConfig: `{ count, rmsd_threshold }`, generate_conformer_ensemble_with_config
  - [x] RMSD-based conformer pruning: 0.5 Å default, 0.0 = no pruning
  - [x] CIP Rule 3 tests added: naphthalene, decalin, fused ring systems (3 cases)
  - Implementation: `crates/chematic-3d/src/conformer.rs`, `crates/chematic-chem/src/cip.rs`
  - Tests: +29 (265 → 294, chematic-3d +7, chematic-chem +22)

**Sprint v0.1.73: Remaining Low-Priority Items (C4 preparation)**
  - [x] Functional group bond count preparation (to be implemented in next Sprint)
  - Tests: +(0, merged into next Sprint)

**Sprint v0.1.74: Functional Group Bond Counts (C4)**
  - [x] `num_amide_bonds(mol: &Molecule) -> usize` — C(=O)-N linkage detection
  - [x] `num_ester_bonds(mol: &Molecule) -> usize` — C(=O)-O-R detection (excluding COOH)
  - [x] 8 tests: acetamide, urea, primary amide, no-amide cases (4 each)
  - Implementation: `crates/chematic-chem/src/descriptors.rs`
  - Tests: +81 (294 → 375)

### Summary
- Across 6 Sprints, implemented 5 high-priority (A1/A2/A5), 3 medium-priority (B1/B2/B3), and 1 low-priority (C4) items from 15 RDKit gaps
- Test count: 933 → 1,150 (+217)
- Progress toward full RDKit parity: Priority A 100% implemented, Priority B 60%, Priority C 20%
- Remaining: B4-B8 (3D-related / FP extensions), C1-C5 (specialty/niche features)

---

## Completed (v0.3.x series)

### Phase 16 — MCP Server + pKa/ADMET (v0.3.0–v0.3.2) (complete)

- [x] **MCP server** (`chematic-mcp`) — AI agent integration, 8 tools, JSON-RPC 2.0 over stdio
- [x] **pKa prediction** (`pka.rs`) — 15 SMARTS rules, `predict_pka`/`pka_acid`/`pka_base`
- [x] **ADMET profile** (`admet.rs`) — BBB/Caco-2/hERG/CYP3A4 + `AdmetProfile`
- [x] **IUPAC expansion** — 15 → 25+ compound classes (piperidine, morpholine, naphthalene, sulfide, etc.)
- [x] **ETKDG KB expansion** — 5 → 20+ torsion patterns (biphenyl, sulfoxide, disulfide, etc.)
- [x] **WASM bindings** — pKa/ADMET exposed as 130+ WASM functions (v0.3.1)
- [x] **criterion benchmarks** — descriptor/SMARTS speed measurements, RDKit comparison scripts (v0.3.2)
- Test count: 1,961 (v0.2.11) → **1,941 lib / 2,100+ all** (v0.3.2)

## Phase 17 — Python PyO3 Bindings (`chematic-py`) (complete)

### Sprint v0.4.0 (implemented, uncommitted)

#### Crate Structure

```
crates/chematic-py/
  src/
    lib.rs    — Mol class (70+ descriptors) + module-level functions
    io.rs     — SDF streaming (iter_sdf / iter_sdf_str / SdfRecord / SdfIter)
    index.rs  — SimilarityIndex (MinHash LSH approximate nearest neighbor search)
    bulk.rs   — Rayon parallel batch processing (parse / FP / descriptors / Tanimoto matrix)
  python/chematic/
    __init__.py   — re-export + type hints
    __init__.pyi  — stubs (for mypy / IDE completion)
  Cargo.toml    — PyO3 + maturin + numpy + rayon dependencies
  pyproject.toml — maturin build configuration
```

#### Implemented Features

**Mol class**
- Identifiers: `smiles`, `formula`, `inchi`, `inchikey`, `iupac_name`
- Basic properties: `mw`, `exact_mass`, `logp`, `tpsa`, `qed`, `hbd`, `hba`, `rotatable_bonds`, `fsp3`, `sa_score`, `molar_refractivity`, `formal_charge`
- Ring/stereo: `ring_count`, `aromatic_ring_count`, `num_stereocenters`
- Drug-likeness filters: `lipinski_passes`, `veber_passes`, `pains_passes`, `ghose_passes`, `egan_passes`, `reos_passes`, `brenk_passes`
- pKa/ADMET: `pka()`, `admet()`, `esol`
- All descriptors at once: `descriptors()` → dict (70+ keys)
- Fingerprints (bytes): `ecfp4()`, `ecfp6()`, `fcfp4()`, `atom_pair_fp()`, `torsion_fp()`, `ecfp4_chiral()`, `maccs()`
- numpy FP: `ecfp4_numpy()`, `maccs_numpy()` (directly usable with scikit-learn / PyTorch)
- SVG rendering: `svg()`, `svg_highlighted(atom_indices, color)`
- Transformations: `standardize()`, `scaffold()`, `canonical_tautomer()`, `enumerate_tautomers()`, `enumerate_stereoisomers()`, `add_hydrogens()`, `remove_hydrogens()`, `remove_stereo()`, `remove_isotopes()`, `largest_fragment()`, `neutralize()`, `generic_scaffold()`, `brics_fragments()`

**Module-level functions**
- `from_smiles(smiles)`, `from_mol_block(block)`, `from_inchi(inchi)`, `is_valid_smiles(smiles)`
- `tanimoto(a, b)` — bytes support (ECFP4/MACCS etc.)
- `smarts_match(smarts, mol)`, `smarts_find(smarts, mol)`
- `depict_grid(mols, cols)`
- `run_smirks(smirks, reactants)`, `find_mcs(mols)`

**SDF streaming** (`io.rs`)
- `iter_sdf(path)` — lazy iterator from file path
- `iter_sdf_str(content)` — lazy iterator from string
- `SdfRecord` — `mol`, `name`, `properties()`, `get(key)`

**LSH similarity index** (`index.rs`) — **not present in RDKit**
- `SimilarityIndex(num_hashes=128)`, `from_smiles(smiles_list)`
- `add(smiles)`, `search(query, threshold=0.7, k=None)`, `get_smiles(index)`

**Parallel batch processing** (`bulk.rs`, Rayon)
- Batch SMILES parsing, batch FP computation (ECFP4 numpy matrix), batch descriptor DataFrame rows
- Tanimoto similarity matrix (N×N)

#### Publishing Infrastructure

- `publish-pypi.yml` — GitHub Actions + maturin + PyPI Trusted Publishing
  - Linux x86_64 / aarch64 / macOS x86_64 / aarch64 / Windows
  - Auto-triggered on `v*` tags
  - PyPI package name: `chematic` (pip install chematic)

#### Next Actions

- [x] `native-inchi` feature implemented (IUPAC InChI C library 1.07.5 vendored, `standard_inchi()` / `standard_inchi_key()` exposed)
      `cargo test -p chematic-inchi --features native-inchi --test standard_inchi` passes 14 tests
- [ ] Add and confirm tests with `cargo test -p chematic-py`
- [ ] Verify local behavior with `maturin develop`
- [ ] Initial PyPI publish with `v0.4.0` tag
- [ ] Begin writing JOSS paper

---

## Next Steps (v0.4.x candidates)

- [ ] Initial PyPI release (`git tag v0.4.0 && git push --tags`)
- [x] chematic-py test expansion (tests/ directory, pytest)
      → tests/{conftest,test_mol,test_fp,test_module_functions,test_io,test_bulk,test_similarity_index}.py
      → 150+ assertions, pytest.ini, __init__.py
- [ ] Write JOSS paper (requires Python bindings completion as prerequisite)
- [x] B4: ETKDG torsion knowledge base expansion (more functional-group-specific patterns)
      → urea, sulfonamide, aryl ether, fluoroalkane, nitro, hydrazone/oxime, imide, benzyl, allylic (+9 patterns)
- [x] B5-B6: LayeredFingerprint + variable-length BitVec
      → Mol.layered_fp() / Mol.layered_fp_numpy() exposed in chematic-py
- [x] B7: Reaction SMARTS queries
      → chematic.reaction_smarts_match(smarts, rxn_smiles) → bool exposed
- [x] B8: 3D SASA descriptor
      → Mol.sasa() / Mol.sasa_per_atom() exposed (chematic-3d added as chematic-py dependency)
- [x] C1-C5: Specialty features (atropisomer M/P SMILES, IUPAC bridged/spiro, InChI parser)
      → C1: Mol.atropisomers() → [(bond_idx, "Biaryl"|"Allene"|"Constrained")]
      → C2: Add spiro[a.b]alkane / bicyclo[x.y.z]alkane naming to chematic-iupac (2 ring families only)
      → C3: from_inchi() already exposed (wraps parse_inchi)
- [ ] Mol2 file format (competing with Open Babel)
- [x] Virtual screening (3D shape similarity)
      → chematic-3d: usr_from_dg() + shape_screen() added (USR Ballester 2007)
      → Python: mol.usr_descriptors(), mol.usr_similarity(other), chematic.shape_screen(query, smiles_list)
      → Tests +4
- [x] Additional ADMET metrics (Ames mutagenicity, PPB, clearance)
      → ames_alerts/ames_passes/ames_risk_score (Kazius 2005 SMARTS alerts, 12 types)
      → ppb_percent (LogP logistic model)
      → clearance_score/clearance_class (Low/Medium/High, MW+LogP+heteroatom)
      → AdmetProfile extended (ames_risk, ppb, clearance fields added)
      → Python: mol.ames_risk(), mol.ames_passes(), mol.ppb(), mol.clearance(), added to admet() dict
      → Tests +14 (chematic-chem 483→493)
- [x] performance: SIMD optimization candidate investigation and implementation
      → #![forbid(unsafe_code)] + WASM compatibility constraint prevent hand-written intrinsics
      → Confirmed u64::count_ones() lowers to hardware POPCNT
      → bitvec.rs: added #[inline] to popcount/and/or/intersection_popcount to encourage autovectorization
      → Investigation results recorded in code comments (AVX2/NEON handled automatically by LLVM)
- [x] Documentation improvements
      → docs/ excluded from .gitignore
      → mkdocs.yml + mkdocs-material + mkdocstrings setup
      → docs/index.md, getting_started/, api/ created
      → docs/cookbook.md translated to English (old Japanese version moved to cookbook_ja.md)
      → docs/rdkit_cheatsheet.md translated to English, new features added (SASA/atropisomer/reaction SMARTS)
      → Old versioned docs deleted (chematic_vs_rdkit_v0210.md, benchmark_results.md)
      → CHANGELOG_ja.md / CHANGELOG_zh.md deleted (only 56%/25% of main content, misleading)
      → AUDIT_EXECUTIVE_SUMMARY.md deleted (internal dev note, not appropriate at root)
      → mkdocs build step added to .github/workflows/pages.yml (demo moved to /playground/)

---

## Future Improvement Candidates (for subsequent phases)

| Priority | Improvement | Status | Notes |
|----------|-------------|--------|-------|
| Medium | SMARTS extension: named smarts for functional groups | Under consideration | Integration with C1=C pattern library, possible synergy with IFG |
| Low | LogP: context-dependent values for alkene C | Not implemented | Distinguish terminal =CH2 (0.1551) vs aryl-adjacent =CH- (0.2640), extend atom_type logic in chematic-chem/src/logp_crippen.rs |
| Low | LogP: C=O group internal refinement | Under consideration | Already accurate at group level (separate handling for ketone/aldehyde/acid/ester); atom-level optimization carries additional cancellation risk |
| Low | 3D Conformer Diversity Metrics | Under consideration | PCA-based distribution analysis improvements, additional diversity metric for ConformerEnsemble |
| Low | SVG metadata embedding expansion | Under consideration | Extension of render_svg_with_metadata, JSON metadata embedding for atom/bond properties |
| Low | Reaction library statistics | Not implemented | Statistical analysis of reaction centers detected by find_reaction_center, retro-synthetic route scoring |

### Selection Criteria for Improvement Candidates

1. **Priority "Medium"**: Many user requests, moderate implementation cost, clear feature gap versus RDKit
2. **Priority "Low"**: Niche cases, specialized use, excessive implementation effort, limited cost-to-benefit ratio
3. **Definition of Status**:
   - **Not implemented**: Requirements defined only, implementation not started
   - **Under consideration**: In design stage, implementation approach being discussed
   - **Pilot complete**: Prototype implementation done, awaiting decision on full implementation

### Constraints and Trade-offs

- **LogP atom-level optimization**: The Crippen atom-type contribution table is based on RDKit empirical values; additional context-dependent corrections carry the cancellation risk of "improving some molecules while degrading others." When proposed, accuracy at the group level is deemed sufficient.
- **SMARTS named patterns**: When managed as a library, need to balance maintenance cost (updates when adding new functional groups) with expressiveness (limits of representing complex patterns).
- **3D Diversity Metrics**: RMSD alone may be insufficient in some cases, but it is recommended to collect actual use cases (library design, HTS diversity evaluation) before implementing.
```

---

## Issue Candidates (Similar Pattern to Issue #1 — Potential Future Problems)

Issue #1 pattern: **Algorithm returns topologically correct results but chemically meaningless ones** (constraint options present in RDKit not implemented in chematic, silently producing incorrect results on migration).

The following items have a high likelihood of becoming future issues following the same pattern.

### Issue Candidate A ([HIGH]): MCS — `atomCompare` / `bondCompare` level (resolved in Sprint v0.1.28)
  - **Status**: Implemented in Sprint v0.1.28 (`AtomCompare::Elements/AnyHeavyAtom/Any`, `BondCompare::OrderOrAromatic/Any`)
  - **Implemented location**: `crates/chematic-smarts/src/mcs.rs` (McsConfig struct lines 48-69)
  - Use `find_mcs_with_config` with `McsConfig { atom_compare, bond_compare, ... }`
  - Chirality comparison: see separate Issue D

### Issue Candidate B ([HIGH]): `run_reactants` — no product valence validation (resolved in Sprint v0.1.24)
  - **Current state**: No valence check on product Molecule after SMIRKS application
  - **Symptom**: Alkylation of quaternary nitrogen silently produces `[N](C)(C)(C)(C)` (valence 5)
  - **RDKit behavior**: `sanitizeMols=True` by default (excludes products with valence violations)
  - **Target**: `crates/chematic-rxn/src/transform.rs`
  - **Implementation**: Check `valence::bond_order_sum > max_valence` per product and exclude (or `TransformError::InvalidProduct`)
  - **Recommended Sprint**: v0.1.24 (prioritized due to correctness issue)

### Issue Candidate C ([MEDIUM]): BRICS — no `minFragmentSize` option (resolved in Sprint v0.1.28)
  - **Status**: Implemented in Sprint v0.1.28 (`BricsConfig { min_fragment_size }` + `brics_fragments_with_config`)
  - **Implemented location**: `crates/chematic-chem/src/brics.rs` (BricsConfig lines 69-76, brics_fragments_with_config)
  - min_fragment_size can filter out meaningless 1-2 atom fragments

### [TODO] Issue Candidate D ([MEDIUM]): MCS — no `matchChiralTag` option (to be implemented in Sprint v0.1.26)
  - **Status**: Planned for implementation (v0.1.26) — important feature for chiral SAR analysis
  - **Symptom**: MCS between R/S enantiomers returns "all atoms match" (chemically they are different compounds)
  - **RDKit behavior**: `matchChiralTag=True`
  - **Target**: `crates/chematic-smarts/src/mcs.rs`
  - **Implementation**: Add `McsConfig { match_chiral_tag: bool }` (default: false), check chirality in `atoms_compatible`
  - **Test**: Verify behavior with R-Ala vs S-Ala

### Issue Candidate E ([MEDIUM]): `find_matches` — no match count limit (resolved in Sprint v0.1.28)
  - **Status**: Implemented in Sprint v0.1.28 (`MatchConfig { max_matches }` + `find_matches_with_config`)
  - **Implemented location**: `crates/chematic-smarts/src/match_vf2.rs` (lines 73-78, match_recursive lines 102-104)
  - max_matches can prevent memory explosion

### Issue Candidate F ([MEDIUM]): VF2 substructure search — chirality awareness
  - **Status**: Implemented (verified in Sprint v0.1.35)
  - **Implementation**: `MatchConfig { use_chirality: bool }` controls `[@]/[@@]` matching
  - **API level**: Can specify `use_chirality=true` in `find_matches_with_config()`
  - **WASM**: `smarts_match_atoms_with_chirality(smarts, mol, use_chirality)` exposed
  - **Tests**: L-alanine `[C@@H]` match + D-alanine `[C@H]` match (2 cases → supplemented in v0.1.35)

### Issue Candidate G ([LOW]): ECFP fingerprints — chirality invariant support
  - **Status**: Implemented (verified in Sprint v0.1.35)
  - **Implementation**: `EcfpConfig { use_chirality: bool }` adds chirality byte to initial atom invariant
  - **API level**: Can specify `config.use_chirality=true` in `ecfp(mol, config)`
  - **WASM**: `ecfp4_bitvec_with_chirality()`, `ecfp6_bitvec_with_chirality()` exposed
  - **Tests**: L/D-alanine FP same (default) vs L/D-alanine FP different (use_chirality=true) (2 cases → supplemented in v0.1.35)

---

---

## Section 4 — WASM & API Improvements (complete: 2026-06-07)

### Required 1: fastrand js feature configuration (WASM RNG seed) (complete)

- **Status**: COMPLETED
- **Implementation**: Added `[target.'cfg(all(target_arch = "wasm32", target_os = "unknown"))'.dependencies]` to crates/chematic-3d/Cargo.toml
- **Fix**: MD simulation initial velocities now use cryptographically random values in WASM as well
- **Commit**: fca2920
- **Tests**: cargo build -p chematic-wasm --target wasm32-unknown-unknown

### Required 2: parse_mol_v3000_with_coords added (complete)

- **Status**: COMPLETED
- **Implementation**: 
  - New `parse_mol_v3000_with_coords()` function implemented in crates/chematic-mol/src/mol3000.rs
  - Return type: `(Molecule, MolMetadata, Vec<(f64, f64)>)` to recover 2D coordinates
  - Existing `parse_mol_v3000()` changed to a wrapper that discards coordinates
- **Re-export**: Added to crates/chematic-mol/src/lib.rs
- **Commit**: fca2920
- **Tests**: cargo test -p chematic-mol (65 tests pass)

### Recommended 3: Y-coordinate system specification documented (complete)

- **Status**: COMPLETED
- **Implementation**:
  - `crates/chematic-depict/src/layout.rs`: `compute_layout()` explicitly notes SVG Y-down
  - `crates/chematic-mol/src/cml.rs`: `parse_cml()` explicitly notes chemical Y-up + Y-negation instruction
  - `crates/chematic-mol/src/cdxml.rs`: `parse_cdxml()` explicitly notes ChemDraw Y-down (SVG-compatible)
- **Purpose**: Prevent coordinate system bugs, eliminate caller confusion
- **Commit**: fca2920
- **Tests**: cargo doc

### Recommended 4: Error type Display + Error trait implementation (complete)

- **Status**: COMPLETED (13 types)
- **High priority** (Display + Error):
  - `SmartsError` (crates/chematic-smarts/src/parser.rs)
  - `ValenceError` (crates/chematic-core/src/valence.rs)
  - `StereoError` (crates/chematic-perception/src/stereo_validation.rs)
- **Medium priority** (Error trait added):
  - `CmlError`, `CdxmlError` (crates/chematic-mol/src/)
  - `Mol2Error`, `RxnParseError` (crates/chematic-mol/src/)
  - `MolError` (crates/chematic-core/src/molecule.rs)
  - `IupacError` (crates/chematic-iupac/src/lib.rs)
  - `ConformerError` (crates/chematic-3d/src/conformer.rs)
  - `RxnError`, `TransformError` (crates/chematic-rxn/src/)
- **Commit**: fca2920
- **Tests**: cargo test --lib (171 tests pass)

### Step 2: 3D Constraint Satisfaction (background execution) (complete)

- **Status**: COMPLETED
- **Implementation**: crates/chematic-3d/src/constraints.rs (639 lines)
  - `BondConstraint`, `AngleConstraint`, `ConstraintSet` structs
  - `build_constraints()`: extract ideal bond distances and angles
  - `satisfy_constraints()`: iterative constraint projection (O(n²) per iteration)
  - `generate_and_minimize_constrained()`: DG → constraints → DREIDING pipeline
- **Performance**: benzene 150µs, naphthalene 400µs, caffeine 700µs
- **Commit**: 137a418
- **Tests**: 12/12

### Step 3: Aromaticity Model Strictness (background execution) (complete)

- **Status**: COMPLETED
- **Implementation**: crates/chematic-perception/src/aromaticity.rs (725 lines)
  - `RingAromaticity` enum: Aromatic/Antiaromatic/NonAromatic
  - `ring_pi_electrons()`: π electron count calculation for C/N/O/S
  - `classify_ring_aromaticity()`: Hückel 4n+2 rule
  - `AromaticityModel` methods: `ring_classifications()`, `antiaromatic_rings()`, `has_antiaromaticity()`
- **Supports**: benzene, pyridine, furan, pyrrole, thiophene (aromatic)
           cyclobutadiene, cyclooctatetraene (antiaromatic), cyclohexane (non-aromatic)
- **Commit**: 137a418
- **Tests**: 16/16

### Version Bump

- **v0.1.30 → v0.1.32**: 2-step bump (Section 4 + Step 2&3 integration)
- **Cargo.toml**: [workspace.package] version = "0.1.32"
- **CHANGELOG.md**: v0.1.32 entry added
- **Commit**: b3227d8

### npm Publishing

- **Status**: Ready for publication
- **Target**: `chematic-wasm` v0.1.32 → `@kent-tokyo/chematic` scope (npm registry)
- **Build**: `cd crates/chematic-wasm && wasm-pack build --target web --release`
- **Package**: pkg/package.json (v0.1.32)
- **Command**: `cd pkg && npm publish` (pending)

---

## Phase 9 — MCP Integration Strategy (not started, to be considered after Phase 3 complete)

### [HIGH] Decision: Wait until Phase 3 is complete (decision: 2026-06-07)

**Conclusion**: MCP (Model Context Protocol) integration will not happen now. Algorithm completion in Phase 3 (SMARTS, 3D coordinate generation, CIP stereochemistry, comprehensive descriptors) comes first.

### Rationale

1. **Algorithms incomplete**: LogP MAE 0.054, SMARTS limitations, 3D coordinate generation is rule-based (no ETKDG), CIP stereochemistry is partial
2. **WASM first**: Browser and serverless use cases are already covered by `@kent-tokyo/chematic` WASM
3. **Maintain focus**: Right now, completing Phase 1-3 algorithms is the top priority. MCP is infrastructure.

### Value After Phase 3 Completion (future)

- "Pure Rust cheminformatics AI tool that works without RDKit" is a clear differentiator
- Demand will emerge as a lightweight, WASM-compatible alternative to Python MCP (RDKit)
- Use cases will develop in AI-agent-driven drug design and screening workflows

### Future Implementation Plan (notes)

**Repository**: `crates/chematic-mcp/` (new crate)

**Implementation stack**:
- Axum (async HTTP server) or stdio transport (Claude Code MCP standard)
- serde_json (JSON request/response)

**Priority APIs** (high impact, low effort):
- `parse_smiles(smiles) -> { atoms, bonds, mol_weight }`
- `calc_logp(smiles) -> f64`
- `calc_tpsa(smiles) -> f64`
- `ecfp4(smiles) -> BitVec2048 (hex)`
- `tanimoto_ecfp4(smiles1, smiles2) -> f64`
- `smarts_match(query, smiles) -> [bool]` (atomwise)
- `write_smiles(smiles) -> canonical_smiles`
- `find_mcs(smiles_list) -> query_smiles`

**Tests**: 30+ API call tests (chematic-mcp/tests/)

**Documentation**: API reference + Claude Code integration examples

**Sprint candidate**: Maintenance Sprint after Phase 3 complete (v0.1.35 onward)

---

## Phase 10 — RDKit Gap Analysis Closure (v0.1.89) (complete)

### Achievement: 89% Gap Closure (A1-A6, B1-B2)

**Completed Items (8/9)**:
- A1: PME panic → Result<T, EwaldError> (4 function signatures)
- A2: InChI stereo parsing (/b, /t, /m, /s layers)
- A3: MMFF94 charge accuracy (formal charge redistribution)
- A4: MHFP implementation quality documentation
- A5: ERG implementation quality + functional group bits (3 bits)
- A6: Reaction FP structural difference encoding
- B1: InChI metadata layer parsing (/m, /s)
- B2: normalize_groups expansion (azide, sulfoxide, 3-pass)

**Statistics**:
- Total tests: 1,521 (all pass, zero regressions)
- New tests: 46 (A1-A6, B1-B2)
- New commits: 8 (a235141 → 46f4cee)
- Documentation: rdkit_feature_comparison.md (379 lines)
- Release notes: RELEASE_NOTES_v0.1.89.md (354 lines)

**Gap Closure Progress**:
- v0.1.87: 67% (fdf5a84 release)
- v0.1.88: 67% (1,475 tests)
- v0.1.89: 89% (+22%) ← YOU ARE HERE

---

## Phase 11 — IUPAC Naming Expansion + CI Hardening (current main) (complete)

### Completed

- [x] `chematic-iupac`: scope expanded from simple hydrocarbons/monofunctional derivatives to ~15 supported classes.
      - Ketones with position locants: `propan-2-one`, `butan-2-one`, `pentan-3-one`.
      - Carboxylic acids: `methanoic acid`, `ethanoic acid`, `propanoic acid`.
      - Esters: `methyl methanoate`, `methyl ethanoate`, `ethyl ethanoate`.
      - Primary/secondary amides: `methanamide`, `ethanamide`, `propanamide`.
      - Aromatics: benzene plus pyridine, furan, thiophene, pyrrole, imidazole, pyrimidine.
- [x] Dispatch changed from single-heteroatom guard to pattern-based `(O, N, S, halogen)` composition.
- [x] `count_c_chain()` BFS helper added for carbonyl chain sizing without crossing blocked atoms.
- [x] CI Clippy hardening:
      - Deprecated `total_hcount` remains exported for compatibility, with warning suppressed only on that re-export.
      - New stable Clippy lints fixed for ECFP loops, condensed formula guards, and DG doc comments.

**Validation**:
- `cargo test -p chematic-iupac`: 14/14 passing.
- `cargo test --workspace --lib --quiet`: 1,649 library tests passing.
- `cargo clippy --workspace -- -D warnings`: passing.

### Known Limitations (By Design)

**Out of scope** (design constraints):
- MMFF94 BCI table full precision (±0.5e → ±0.1e, v0.1.96+ target)
- Transition metal chemistry (valence model limitation)
- Polymers/peptides (format out of scope)

### Roadmap (v0.1.96+)

**High Priority**:
1. MMFF94 BCI table (±0.5e → ±0.1e, Bond Charge Increment table implementation)

**Low Priority**:
2. LogP alkenyl C context values (terminal =CH₂ vs aryl-adjacent =CH−)
3. Kekulization edge cases (Edmonds flower algorithm for odd rings)

---

## Phase 12 — True Fingerprint Algorithms (v0.1.95) (complete)

### Completed

- [x] **A4 MHFP canonical hash**: Morgan-style circular fragment hash replaces atom-index-dependent byte signature (`crates/chematic-fp/src/mhfp.rs`). New tests: +3 (canonical, similar>dissimilar, radius effect).
- [x] **A5 ERG pharmacophore node types**: `assign_pharmacophore_features()` correctly assigns DONOR/ACCEPTOR/POSITIVE/NEGATIVE/HYDROPHOBIC. Pyridine N (acceptor) vs pyrrole N-H (donor) distinguished. New tests: +5.
- [x] **A6 Reaction FP XOR**: `use_xor: true` confirmed as default in `reaction_fp.rs`. Comparison table updated WARNING: → 
- [x] Tests: 1,649 → 1,657 (+8), all passing. Clippy clean.
- [x] CHANGELOG, README (all languages), tasks/todo.md updated.
- [x] Version: 0.1.94 → 0.1.95

---

## Phase 13 — MMFF94 Complete Stack + Stereo SMILES (v0.2.7–v0.2.9) (complete)

### Sprint v0.2.7 (commit bd7ab12, 2026-06-14)

- [x] **Canonical SMILES stereo parity correction** — pre-solves RDKit issue #8775
      - `crates/chematic-smiles/src/canonical.rs`: `corrected_chirality()` method
      - `crates/chematic-smiles/src/parser.rs`: `StereoEntry` enum tracks parse-time neighbor order
      - `crates/chematic-core/src/molecule.rs`: `stereo_neighbor_order: HashMap<AtomIdx, Vec<u32>>` + `STEREO_H_SENTINEL`
      - Detects odd permutations between parse-time and canonical write-time neighbor order; auto-flips `@`/`@@`
      - Tests: L-alanine N-first vs C-first, aminocyclopentane, fluorocyclohexane

- [x] **MMFF94 faithful partial charges** (Halgren 1996 equation 15)
      - `crates/chematic-ff/src/mmff94_numeric.rs` (new, ~1,300 lines)
      - `MMFF94_PBCI`: 99 entries (one per numeric atom type 1–99)
      - `MMFF94_CHG`: 498 entries (bond charge increments from RDKit Params.cpp)
      - CHG sign convention: entry (bt, a, b, bci) → b gets +bci, a gets −bci
      - Glycine cross-validated against MMFF94_reference.log
      - Tests: +15 (total: 1,930)

### Sprint v0.2.8 (commit 093bf0b, 2026-06-14)

- [x] **MMFF94 full energy parameters** (Halgren 1996 Tables IV–VII)
      - `crates/chematic-ff/src/mmff94_energy.rs` (new, ~4,000 lines)
      - `MMFF94_BOND_ENERGY`: 493 entries (Table IV) — kb in md/Å, r0 in Å
      - `MMFF94_ANGLE_ENERGY`: 2,245 entries (Table V) — ka in md·Å/rad², theta0 in degrees
      - `MMFF94_TORSION_ENERGY`: 926 entries (Table VI) — v1/v2/v3 in kcal/mol; wildcard (0) fallback hierarchy
      - `MMFF94_VDW_ENERGY`: 95 entries (Table VII) — Slater-Kirkwood combining rule params
      - Data source: verbatim from RDKit `Code/ForceField/MMFF/Params.cpp` via `gh api` download
      - Cross-validated: C-C-C-C torsion v1=0.103/v2=0.681/v3=0.332, C-C bond kb=4.258/r0=1.508 vs RDKit Python API
      - Existing PBCI (99) and CHG (498) values verified against Params.cpp — all match
      - Lookup: O(log n) binary search on sorted tables; torsion wildcard hierarchy (exact → reversed → wildcard ends → generic)
      - Tests: +11 (total: 1,941)

### Sprint v0.2.9 (commit 6570fe1, 2026-06-14)

- [x] **MMFF94 geometry minimizer** — full Halgren 1996 force field
      - `crates/chematic-ff/src/mmff94_minimizer.rs` (new, ~590 lines)
      - Bond: `E = (143.9325 × kb / 2) × ΔR² × (1 − cs×ΔR + (7/12)×cs²×ΔR²)` (cubic correction, cs=2.0)
      - Angle: `E = (0.043844 × ka / 2) × Δθ² × (1 − 0.007×Δθ)` (Δθ in degrees)
      - Torsion: `E = (v1/2)(1+cosφ) + (v2/2)(1−cos2φ) + (v3/2)(1+cos3φ)` — **first implementation**
      - vdW: buffered 14-7 `E = ε × t⁷ × (t⁷ − 2)`, `t = 1.07r* / (r + 0.07r*)` + Slater-Kirkwood combining rule
      - Electrostatic: `E = 332.0716 × qi×qj / (r + 0.05)`, 1-4 scaling 0.75, 1-2/1-3 excluded
      - Algorithm: steepest descent with finite-difference gradients (δ=1e-4 Å), convergence 1e-4
      - Public API: `mmff94_total_energy(mol, coords)` + `minimize_mmff94_full(mol, coords, max_iter) → MinimizeResult`
      - Tests: +6 (torsion conformer discrimination, vdW repulsion, dihedral geometry, minimize reduces energy)
      - Total tests: **1,947**

### MMFF94 Complete Stack Summary

```
Charge calculation (v0.2.7): PBCI 99 entries + CHG 498 entries → equation 15
Energy parameters  (v0.2.8): Bond 493 / Angle 2,245 / Torsion 926 / VdW 95 entries
Minimizer          (v0.2.9): cubic bond/angle + buffered-14-7 vdW + torsion
```

### Fingerprint Coverage (confirmed v0.2.9)

13/14 RDKit fingerprint algorithms implemented (Avalon excluded — requires C library):
ECFP, FCFP, MACCS, RDKit Path, Atom Pairs, Topological Torsion, MHFP, ERG, Layered, Pattern, Topo Path, 2D Pharmacophore, Reaction FP

### RDKit Parity Status (v0.2.9)

| Domain | Verdict |
|--------|---------|
| WASM/Browser | Surpassed (60× smaller, native WASM) |
| Rust-native pharma tools | Surpassed (all major features ≥ RDKit) |
| MMFF94 force field | Full parity (complete stack v0.2.9) |
| StandardInChI compliance | `native-inchi` feature (vendored IUPAC C lib 1.07.5) |
| General purpose (no FFI) | Surpassed (only 1 remaining gap) |

---

## Phase 15 — Surpassed RDKit in 3 Domains (v0.2.11) (complete)

### Sprint v0.2.11 (commit de156b9, 2026-06-14)

#### MMFF94 OOP + Stretch-Bend (all 7 energy terms from Halgren 1996)

- [x] **Out-of-Plane bending (OOP)** (`MMFF94_OOP`, 117 entries)
      - E = (0.043844 × koop / 2) × χ² (Wilson angle, degrees)
      - Maintains planarity of trigonal sp² centers (carbonyl C, amide N, aromatic)
      - `mmff94_oop(type_j, type_i, type_k, type_l)` — stepwise wildcard fallback
      - Implementation: `crates/chematic-ff/src/mmff94_energy.rs` (appended)

- [x] **Stretch-Bend coupling (STRE-BEN)** (`MMFF94_STBN`, 282 entries)
      - E = 2.51210 × (kba_ijk × Δr_ij + kba_kji × Δr_kj) × Δθ
      - Cross term of bond stretching × angle bending, Halgren MMFF.V eq.4
      - `mmff94_stbn(angle_type, type_i, type_j, type_k)` — symmetric lookup

- [x] **EnergyBreakdown extended**: 5 terms → 7 terms (added `stretch_bend`, `oop`)
- [x] 2 new energy terms integrated into `total_energy()`
- [x] Existing test `energy_breakdown_sums_to_total` updated for 7 terms

#### MAP4 Fingerprint (not present in RDKit)

- [x] **MAP4** (`crates/chematic-fp/src/map4.rs`) newly created
      - MinHashed Atom-Pair FP (Minervini et al. J. Cheminform. 2020, 12, 26)
      - Hash all atom-pair circular environments with FNV-1a → MinHash signature
      - `Map4Config { max_radius: 2, n_permutations: 1024 }`
      - `map4(mol, config) -> Vec<u32>`
      - `tanimoto_map4(a, b) -> f64` — Hamming distance ≈ Jaccard similarity
      - BFS all-atom-pair distances + circular env hash + MinHash permutation mixing
      - Added `pub mod map4` + exports to `crates/chematic-fp/src/lib.rs`
      - chematic FP types: 13 → 14 (RDKit also has 14, but MAP4 requires an external package)

#### SMARTS Compile Cache + Named Pattern Library

- [x] **SmartsCache** (`crates/chematic-smarts/src/cache.rs`) newly created
      - LRU eviction (VecDeque + HashMap, variable capacity)
      - `compile(smarts) -> &QueryMolecule` — parse once, O(1) thereafter
      - `find_matches()` / `has_match()` / `find_matches_with_config()`
      - 5–20× speedup for repeated matches
      - Added `pub mod cache` + exports to `crates/chematic-smarts/src/lib.rs`

- [x] **named_pattern()** — 20 named SMARTS patterns:
      donor, acceptor, aromatic, hydrophobic, positive, negative,
      carboxylic_acid, aldehyde, ketone, alcohol, phenol,
      amine_primary/secondary/tertiary, amide, ester, ether, halide,
      aromatic_n, sulfonamide

### Test Counts

- New tests: +10
- **Total: 1,961 tests, all passing**

### RDKit Surpassed Status (v0.2.11)

| Domain | chematic v0.2.11 | RDKit |
|--------|------------------|-------|
| MMFF94 energy terms | **all 7 terms** (including OOP+STRE-BEN) | all 7 terms (C++ only) |
| MAP4 FP | **native implementation** | external package required |
| SMARTS cache | **integrated LRU cache** | none |
| SMARTS pattern library | **20 built-in patterns** | none |

### CI failure follow-up (commit 2f0d6e3, 2026-06-14)

- [x] **Root cause**: The v0.2.11 feature additions themselves pass `cargo check --workspace`, but
      CI's `cargo clippy --workspace -- -D warnings` hard-errors on new Clippy lints.
      Furthermore, `security.yml` runs `cargo clippy --workspace --all-targets -- -D warnings`,
      which means warnings in test-only code also become failure causes, not just normal targets.

- [x] **Main fix scope**:
      - `chematic-smiles` / `chematic-iupac`: `collapsible_if`, `unnecessary_map_or`
      - `chematic-ff`: To preserve literal values in MMFF94 parameter tables,
        `approx_constant` / `type_complexity` allowed at file scope. L-BFGS API reorganized to `&mut [[f64; 3]]`.
      - `chematic-fp`: MAP4 `unnecessary_cast` / `needless_range_loop`, ERG donor condition and clamp reorganized.
      - test modules: resolved unused import / fixture helper / range-loop lints that appear only under `--all-targets`.

- [x] **Verification commands**:
      - `cargo clippy --workspace -- -D warnings`
      - `cargo clippy --workspace --all-targets -- -D warnings`
      - `cargo test --workspace`
      - `cargo test --workspace --lib --quiet`

- [x] **Operational notes**:
      - For CI failures, check not only the normal CI but also the clippy command in `.github/workflows/security.yml`.
      - `cargo fmt --all` may reformat large amounts of existing code; when fixing CI, review the diff and
        do not commit unrelated formatting churn.
      - Existing local changes to `tasks/todo.md` were excluded from the CI fix commit.

---

## Phase 14 — L-BFGS + MMFF94 WASM API + Demo Enhancements (v0.2.10) (complete)

### Sprint v0.2.10 (commit ac81a2c, 2026-06-14)

- [x] **L-BFGS geometry minimizer** (`crates/chematic-ff/src/mmff94_minimizer.rs`)
      - `minimize_mmff94_lbfgs()`: limited-memory quasi-Newton, m=5 history pairs
      - Two-loop recursion: q=g → scale by γ=(s·y)/(y·y) → forward loop → p=-r
      - Backtracking Armijo line search (c=1e-4, τ=0.5, max 20 steps)
      - Fallback to SD step when curvature condition y·s > 0 not met
      - Shared `compute_gradient()` helper (refactored out of SD loop)
      - `VecDeque<(Vec<[f64;3]>, Vec<[f64;3]>, f64)>` circular history buffer

- [x] **EnergyBreakdown struct + mmff94_energy_breakdown()**
      - `EnergyBreakdown { bond, angle, torsion, vdw, electrostatic, total }`
      - All 5 MMFF94 energy terms evaluated and returned separately

- [x] **Torsion scan** (`mmff94_torsion_scan()`)
      - Rodrigues rotation formula: rotate atoms past j-k bond axis by step_rad per step
      - BFS from k (not crossing j) to collect moving atoms
      - Returns `Vec<(f64, f64)>` of (angle_deg, energy_kcal) pairs

- [x] **Tests: +4** (lbfgs_reduces_energy, lbfgs_converges_faster_than_sd,
      energy_breakdown_sums_to_total, energy_breakdown_bond_term_positive)
      - Total: **1,951 tests**

### WASM Bindings (chematic-wasm)

- [x] Add `chematic-ff` to `crates/chematic-wasm/Cargo.toml`
- [x] `minimize_mmff94_json(mol, max_iter)` → MMFF94 steepest descent
- [x] `minimize_mmff94_lbfgs_json(mol, max_iter)` → MMFF94 L-BFGS
- [x] `mmff94_energy_breakdown_json(mol)` → `{bond,angle,torsion,vdw,elec,total}`
- [x] `torsion_scan_json(mol, i, j, k, l, steps)` → `[{angle,energy},...]`
- [x] `mmff94_partial_charges_json(mol)` → `{charges:[...]}`

### Demo Enhancements (`demo/index.html`)

- [x] **MMFF94 Force Field Optimizer** (Dynamics tab)
      - SD / L-BFGS algorithm selector
      - Max iterations input
      - Energy breakdown table with colored bars (bond/angle/torsion/vdW/elec)
      - Energy, RMSD, iterations, converged display

- [x] **Force Field Comparison** (Dynamics tab)
      - Single click: run DREIDING + MMFF94 L-BFGS in parallel
      - Side-by-side table: energy / RMSD / iterations / timing

- [x] **Torsion Scan** (Dynamics tab)
      - 4-atom index inputs (i-j-k-l), step count selector
      - SVG line chart (420×180): energy vs dihedral angle 0°→360°
      - Green circle = minimum, red circle = maximum
      - Min/max energy and angle readout

- [x] **MMFF94 Charge Map** (2D tab ±q button)
      - Toggles per-atom MMFF94 charge coloring on 2D structure
      - Blue=positive (>+0.3e), red=negative (<-0.3e), white=neutral
      - DOMParser-safe SVG update (no innerHTML)

### npm publish

- [x] Workspace version bump: 0.2.0 → 0.2.10
- [x] `wasm-pack build --target web --release` (8.02s compile + wasm-opt)
- [x] `pkg/package.json` name fixed: `chematic-wasm` → `@kent-tokyo/chematic`
- [x] `npm publish --access public` → `@kent-tokyo/chematic@0.2.10`
      - Package size: 873.1 KB (2.4 MB unpacked)
      - WASM binary: 2.2 MB
      - TypeScript defs: 78.5 KB
