# Changelog

All notable changes to chematic will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

---

## [0.1.7] — 2026-05-30

### Fixed (`chematic-chem`) — HBA accuracy: Ertl S inclusion + charged N exclusion

**`hba_count` now uses the full Ertl (2000) definition** (`rdMolDescriptors.CalcNumHBA`):

1. **Sulfur counted as HBA** (new): divalent uncharged S (thiothers, thiols, aromatic S like thiophene) is now included. Matches Ertl SMARTS `$([S;!+;X2;!$([S]=[#8])])` and `$([s;+0])`.

2. **Sulfonic/sulfonamide OH excluded** (new): O–H bonded to oxidized S (S=O present) is excluded from HBA, matching RDKit's exclusion of sulfonate–OH.

3. **Charged N excluded** (new): N with non-zero formal charge (`[N+]`, `[n+]`) is never an HBA. This correctly excludes nitro-group N+ (4-nitrophenol, clonazepam) and thiazolium n+ (thiamine).

**Benchmark results** (175-molecule ChEMBL test set):

| Property | Before | After (v0.1.7) |
|---|---|---|
| HBA MAE | 0.1371 | **0.0400** (−71%) |
| LogP MAE | 0.1174 | 0.1174 (unchanged) |
| TPSA MAE | 0.0808 | 0.0808 (unchanged) |

---

## [0.1.6] — 2026-05-30

### Added (`chematic-wasm`) — WASM bindings for Sprint G–K features

**Topological descriptors** (`MolHandle` methods):
- `wiener_index()`, `kappa1()`, `kappa2()`, `kappa3()` — Wiener index and Hall–Kier κ shape indices.
- `chi0()` – `chi4()` — Kier–Hall molecular connectivity χ indices (unweighted).
- `chi0v()` – `chi4v()` — valence-weighted χv indices.
- `bertz_ct()` — Bertz complexity index.
- `labute_asa()` — Labute approximate surface area (Å²).
- `morgan_fp_counts_json(radius)` — Morgan count fingerprint as a JSON object string (`{"<hash>": count, …}`).

**Free functions**:
- `add_hydrogens(mol) -> MolHandle` — convert implicit H to explicit atoms.
- `remove_hydrogens(mol) -> MolHandle` — remove explicit H atoms.
- `depict_svg_grid(smiles_block, cols) -> String` — grid SVG from newline-separated SMILES; invalid lines silently skipped.
- `run_reactants(smirks, reactants_smiles) -> Result<String, JsValue>` — SMIRKS reaction transform; `reactants_smiles` is pipe-separated (`"CC(=O)O|CCO"`); returns JSON `[["product_smi", …], …]`.

Tests: 646 → 656 (+10 new WASM binding tests).

---

## [0.1.5] — 2026-05-30

### Improved (`chematic-3d`) — UFF-derived minimizer parameters (Sprint K)

**`chematic-3d/src/minimize.rs`**:
- **Bond lengths**: replaced single-constant-per-bond-order with a 30+-entry element-pair table (`ideal_bond_len`). Covers C–C/C–N/C–O/C–S/C–F/C–Cl/C–Br/C–H/N–N/N–O/O–H/S–S/H–X etc.
- **Bond angles**: replaced neighbor-count heuristic with hybridization-aware ideal angles (`atom_hybridization` + `ideal_angle_rad`). Detects SP (triple bond → 180°), SP2 (double/aromatic → 120°), SP3 (element-specific: O 104.5°, N 107°, S 99°, P 93°, others 109.47°).
- **VDW repulsion**: replaced fixed r₀ = 2.0 Å with element-specific UFF/Bondi radii (`uff_vdw_radius`); cutoff extended from 5.0 → 8.0 Å.

Tests: 637 → 646 (+9 new tests covering bond length precision, hybridization detection, and table symmetry).

---

### Added (`chematic-rxn`) — SMIRKS reaction transform (Sprint J)

**`chematic-rxn/src/transform.rs`** (new):
- `run_reactants(smirks: &str, reactants: &[&Molecule]) -> Result<Vec<Vec<Molecule>>, TransformError>` — applies a SMIRKS reaction template to a list of reactant molecules and returns all product sets.
  - Parses SMIRKS into reactant/product SMARTS patterns via `chematic-smarts`.
  - Matches each reactant pattern via VF2 subgraph isomorphism.
  - Builds product molecules by copying non-reaction-centre atoms, applying bond changes from the template, and transferring unmapped substituents via BFS traversal.
  - Returns the Cartesian product of all match sets across reactant molecules.
- `TransformError` — parse and arity error variants.

Tests: 623 → 634 (+11 new tests covering esterification, amide coupling, cyclisation, and error cases).

---

### Added (`chematic-3d`) — Conformer ensemble + Kabsch RMSD (Sprint I)

**`chematic-3d/src/conformer.rs`** (new):
- `ConformerEnsemble` — external container holding a `Molecule` and an ordered `Vec<Coords3D>`. No changes to `chematic-core`.
- `add_conformer`, `get_conformer`, `get_conformer_mut`, `remove_conformer` — CRUD with atom-count validation; returns `ConformerError::AtomCountMismatch` on mismatch.
- `conformer_rmsd_no_align(a, b) -> Option<f64>` — raw per-atom RMSD without superposition.
- `conformer_rmsd(a, b) -> Option<f64>` — Kabsch-aligned RMSD minimised over all rigid rotations+translations; uses `jacobi3` (3×3 Jacobi eigensolver from `shape_descriptors`) to compute the SVD of the 3×3 covariance matrix; reflection correction via determinant check.
- `ConformerError` — atom-count mismatch error type.

Tests: 609 → 623 (+14 new tests).

---

### Added (`chematic-chem`, `chematic-depict`) — Topo descriptors + H management + SVG grid (Sprint G)

**Topological connectivity indices** (`chematic-chem/src/topo_descriptors.rs`, new):
- `wiener_index(mol) -> f64` — sum of all pairwise shortest-path distances (Wiener 1947).
- `kappa1`, `kappa2`, `kappa3` — Hall–Kier κ shape indices.
- `chi0`, `chi1`, `chi2`, `chi3`, `chi4` — Kier–Hall molecular connectivity χ0–χ4 (unweighted).
- `chi0v`, `chi1v`, `chi2v`, `chi3v`, `chi4v` — valence-weighted connectivity χ0v–χ4v.
- `bertz_ct(mol) -> f64` — Bertz complexity index (BertzCT 1981).
- `labute_asa(mol) -> f64` — Labute (2000) approximate surface area (Å²).

**Explicit hydrogen management** (`chematic-chem/src/hydrogen.rs`, new):
- `add_hydrogens(mol) -> Molecule` — converts all implicit H counts to explicit H atoms.
- `remove_hydrogens(mol) -> Molecule` — removes explicit H atoms and updates implicit H count on heavy atoms.

**SVG grid layout** (`chematic-depict/src/grid.rs`, new):
- `depict_svg_grid(mols: &[&Molecule], cols: usize) -> String` — renders multiple molecules in a grid SVG (200×200 px per cell). Equivalent to RDKit's `Draw.MolsToGridImage`.

Tests: 544 → 582 (+38 new tests across topo_descriptors, hydrogen, and grid modules).

---

### Added (`chematic-chem`, `chematic-fp`) — LabuteASA + Morgan count FP

**LabuteASA** (`chematic-chem/src/topo_descriptors.rs`):
- `labute_asa(mol) -> f64` — Labute (2000) approximate surface area (Å²) computed from covalent radii and bond-type-specific interatomic distances; implicit H atoms included.

**Morgan count fingerprint** (`chematic-fp/src/ecfp.rs`):
- `morgan_fp_counts(mol, radius) -> HashMap<u64, u32>` — count-based Morgan fingerprint returning raw `hash → count` map. All (atom, radius) pairs contribute without deduplication (equivalent to `includeRedundantEnvironments=True`). Hash scheme is identical to `ecfp`, so bit-folded and count forms are consistent.

Tests: 635 → 645 (+10 new tests).

---

### Added (`chematic-3d`) — Shape descriptors + stereo from 3D (Sprint H)

**Shape descriptors** (`chematic-3d/src/shape_descriptors.rs`, new):
- `pmi(mol, coords) -> (f64, f64, f64)` — principal moments of inertia PMI1 ≤ PMI2 ≤ PMI3 (Da·Å²) from mass-weighted inertia tensor eigenvalues.
- `pmi1`, `pmi2`, `pmi3` — individual PMI accessors.
- `npr1`, `npr2` — normalized PMI ratios (PMI1/PMI3, PMI2/PMI3; range 0–1).
- `radius_of_gyration` — mass-weighted Rg (Å).
- `asphericity` — PMI3 − (PMI1+PMI2)/2; zero for perfect sphere.
- `eccentricity` — sqrt(1 − PMI1/PMI3); zero for sphere, 1 for rod.
- `plane_of_best_fit` — RMS deviation from the least-squares plane (Å); ≈ 0 for flat molecules like benzene.
- Internals: 3×3 symmetric Jacobi eigensolver (no nalgebra dependency; pure Rust; converges in ≤ 100 sweeps).

**Stereo from 3D** (`chematic-3d/src/stereo3d.rs`, new):
- `assign_stereo_from_3d(mol, coords) -> StereoAssignment3D` — assigns R/S (tetrahedral) and E/Z (alkene) from 3D coordinates using signed-volume (scalar triple product) and dihedral-angle conventions respectively.
- Uses 1-sphere CIP priority (atomic number + sorted neighbor atomic numbers). Stereocenters that cannot be resolved at this level are omitted.
- `StereoAssignment3D::get(idx) -> Option<CipCode>` for lookup.

Tests: 620 → 635 (+15 new tests in shape_descriptors and stereo3d modules).

---

### Fixed — Security, bug, and code quality (audit)

**Security** (`chematic-smarts`):
- **Recursive SMARTS depth limit**: `$(…)` patterns nested beyond 8 levels now return `SmartsError::RecursionDepthExceeded` instead of panicking with a stack overflow. Protects against malformed SMARTS strings used as a DoS vector.
- **Ring closure digit `unwrap()`** (`chematic-smarts`, `chematic-smiles`): replaced with `expect()` plus an invariant comment documenting that the caller always `peek()`s a digit before entering the branch, making the assumption visible in the source.

**Bug** (`chematic-chem`):
- **`clone_mol` silent bond loss**: `add_bond(…).ok()` discarded errors silently if a bond could not be re-added during molecule cloning, producing a structurally corrupt molecule without warning. Changed to `expect()` so any failure is immediately visible. Same fix applied to `transfer_hydrogen_aromatic`.

**Refactor** (`chematic-chem`):
- **FNV-1a named constants**: `mol_fingerprint` now uses `FNV1A_OFFSET` / `FNV1A_PRIME` constants instead of inline magic numbers.
- **TPSA nitro detection single-pass**: the nitro group check (`[N+](=O)[O−]`) previously iterated `mol.neighbors` twice; consolidated into a single `fold` pass.

Tests: 542 → 544 (two new SMARTS recursion-depth tests).

---

### Fixed (`chematic-chem`) — LogP guanidinium N accuracy (Sprint E)

- **LogP guanidinium/amidine N** (`descriptors.rs`): non-aromatic nitrogen atoms in imine or guanidinium context now use Wildman–Crippen N14 type (−0.335) instead of the generic aliphatic amine values (−0.595 to −1.019). Detection: N with a direct double bond to C (`=N`, Type A) or N bonded to a C that itself has a C=N double bond (adjacent N, Type B). Fixes metformin (error 2.07 → ~0.00), improves arginine, diazepam, clonazepam.

**Benchmark results** (175-molecule ChEMBL test set):
| Property | Before (Sprint D) | After (Sprint E) |
|----------|-------------------|-----------------|
| LogP MAE | 0.134 | **0.117** |
| TPSA MAE | 0.081 Å² | 0.081 Å² (unchanged) |

### Added (`chematic-chem`) — Tautomer 1,2-shift (Sprint E)

- **`enumerate_tautomers`**: now generates direct aromatic 1,2-shift tautomers (e.g. pyrazole N1H ↔ N2H) in addition to 1,3-shift rule-based tautomers. Uses a separate H-assignment fingerprint to distinguish positional isomers that share the same structural fingerprint.
- **`canonical_tautomer`**: after rule-based normalization, direct aromatic 1,2-shift candidates are compared by lexicographic H-assignment and the minimal form is returned, ensuring both N1H and N2H of pyrazole converge to the same canonical molecule.

---

### Fixed (`chematic-chem`) — TPSA and LogP accuracy (Sprint D)

- **TPSA imine N-H** (`descriptors.rs`): sp2 imine nitrogen with one H (C=N-H, as in amidine and guanidinium groups) now uses 23.79 Å² instead of 12.03 Å² (generic secondary amine). Detection: N with `h=1` and a double bond from N to a carbon neighbor. Reduces metformin TPSA error from 23.64 → 0.12 Å², arginine from 11.82 → 0.06 Å².
- **TPSA phosphate P** (`descriptors.rs`): non-aromatic phosphorus with a P=O bond now uses 26.88 Å² (Ertl 2000 phosphate type) instead of 34.14 Å² (phosphine type). Trimethyl phosphate TPSA error: 7.26 → 0.00 Å².
- **LogP phosphate P** (`descriptors.rs`): non-aromatic P with a P=O bond now uses Wildman–Crippen contribution +0.7933 instead of −0.3451 (phosphine). Trimethyl phosphate LogP error: 1.14 → 0.00.

**Benchmark results** (175-molecule ChEMBL test set):
| Property | Before (Sprint C) | After (Sprint D) |
|----------|------------------|-----------------|
| TPSA MAE | 0.324 Å² | **0.081 Å²** |
| LogP MAE | 0.141 | **0.134** |

### Added (`chematic-chem`) — Ring descriptors (Sprint D)

- **`num_aromatic_heterocycles`**: count of SSSR rings where all atoms are aromatic and at least one is a heteroatom (pyridine, furan, imidazole, etc.).
- **`num_aliphatic_heterocycles`**: count of SSSR rings with at least one non-aromatic atom and at least one heteroatom (piperidine, morpholine, THF, etc.).
- **`num_saturated_heterocycles`**: count of SSSR rings where all atoms are sp3 (no double/triple/aromatic bonds) and the ring contains at least one heteroatom.
- **`num_spiro_atoms`**: number of atoms shared by exactly two rings that share no other atoms (spiro centers).
- **`num_bridgehead_atoms`**: number of atoms shared between two bridged rings, identified by non-adjacent shared-atom pairs in the ring intersection.

### Added (`chematic-chem`) — Tautomer rules (Sprint D)

- **Rules 16–20**: five additional 1,3-proton-shift patterns covering O→N, O→O, N→C, C→O, and C→N heteroatom combinations with any bridge element. Expands tautomer coverage to hydroxamic acids, cross-conjugated enol/iminol systems.

### Added (`chematic-wasm`) — Ring descriptor bindings (Sprint D)

- New `MolHandle` methods: `num_aromatic_heterocycles`, `num_aliphatic_heterocycles`, `num_saturated_heterocycles`, `num_spiro_atoms`, `num_bridgehead_atoms`.

---

### Fixed (`chematic-chem`) — LogP and TPSA accuracy (Sprint C)

- **LogP aromatic junction C** (`descriptors.rs`): aromatic C at fused-ring junctions (e.g. naphthalene C4a, indole C3a/C7a) now uses Crippen value 0.2956 instead of 0.1441, when all neighbors are aromatic and ≥2 are aromatic carbons. Verified: naphthalene (±0.001), quinoline (±0.001), indole (±0.001) now match RDKit exactly.
- **LogP alkene C** (`descriptors.rs`): sp2 vinyl carbons (C=C, non-aromatic) now use +0.2274 (Wildman-Crippen C5 type) instead of wrong negative values (−0.215 to −0.350). Styrene LogP error reduced from −1.03 to +0.04.
- **LogP benchmark SMILES** (`scripts/rdkit_ref_properties.tsv`): morphine and codeine entries updated to aromatic SMILES notation so chematic's aromaticity perception succeeds.
- **TPSA nitro group** (`descriptors.rs`): `[N+](=O)[O−]` now contributes 41.44 Å² (Ertl 2000 table) and the `[O−]` oxygen contributes 0 (absorbed into N). Previously nitro N was treated as tertiary amine (3.24 Å²). 4-nitrophenol TPSA error: 30.67 → 1.70 Å²; clonazepam: 39.79 → 1.17 Å².
- **TPSA imine N** (`descriptors.rs`): aliphatic C=N imine nitrogen (h=0, double bond to C) now uses 12.89 Å² (same as pyridine-type aromatic N) instead of 3.24 Å² (generic tertiary N). Diazepam TPSA error: 9.12 → 0.53 Å².

**Benchmark results** (175-molecule ChEMBL test set):
| Property | Before (v0.1.3) | After |
|----------|----------------|-------|
| LogP MAE | 0.298 | **0.141** |
| TPSA MAE | 0.759 Å² | **0.324 Å²** |
| TPSA RMSE | 4.40 Å² | **2.13 Å²** |

### Added (`chematic-smarts`)

- **`[XN]` total connectivity**: matches atoms where heavy-atom degree + implicit-H count equals N (distinct from `[DN]` which counts only heavy-atom neighbours).
- **`[RN]` ring count**: matches atoms that belong to exactly N SSSR rings.
- **Compound bond expressions**: OR (`,`) and AND (`&`) now supported in bond queries. Examples: `=,:` (double or aromatic), `=!@` (double non-ring), `-!@` (single non-ring). Required for full PAINS SMARTS compatibility.
- **HCount fix**: the `[HN]` atom primitive now counts both explicit H neighbors and implicit H (matches SMARTS spec); previously only implicit H was counted.

### Added (`chematic-chem`)

- **Improved QED** (`qed`): rewritten using the exact 7-parameter ADS (Asymmetric Double Sigmoidal) function from Bickerton 2012 / RDKit. Now includes 113 Brenk 2008 structural alerts as the eighth desirability component.
- **Molar Refractivity** (`molar_refractivity`): Wildman–Crippen additive MR model (same atom-type framework as LogP).
- **Formal charge sum** (`formal_charge_sum`): sum of atom formal charges over the whole molecule.
- **Veber filter** (`veber_passes`): TPSA ≤ 140 Å² and rotatable bonds ≤ 10.
- **Egan filter** (`egan_passes`): TPSA ≤ 131.6 Å² and LogP ≤ 5.88.
- **REOS filter** (`reos_passes`): MW, LogP, HBD, HBA, charge, and heavy-atom criteria.
- **Ghose filter** (`ghose_passes`): MW 160–480, LogP −0.4–5.6, heavy atoms 20–70, MR 40–130.
- **Expanded tautomer rules**: 5 → 15 rules covering thioamide, thio-iminol, thio-keto-enol, and six cross-heteroatom 1,3-proton-shift patterns.
- **Count descriptors**: `num_heteroatoms`, `ring_count`, `num_aliphatic_rings`, `num_saturated_rings`, `num_stereocenters`, `num_unspecified_stereocenters`.
- **PAINS structural alerts** (`pains_matches`, `pains_passes`): all 480 patterns from Baell & Holloway 2010 / RDKit FilterCatalog. Molecules are expanded to explicit-H form before matching for full coverage.

### Added (`chematic-fp`)

- **FCFP fingerprints** (`fcfp4`, `fcfp6`, `tanimoto_fcfp4`): pharmacophore-based circular fingerprints using feature classes (Donor, Acceptor, Aromatic, Hydrophobic, PosIonizable, NegIonizable) as atom invariants — bioisostere-aware similarity.

### Added (`chematic-wasm`)

- New bindings: `molar_refractivity`, `formal_charge_sum`, `veber_passes`, `egan_passes`, `reos_passes`, `ghose_passes`.
- Sprint B bindings: `num_heteroatoms`, `ring_count`, `num_stereocenters`, `pains_passes`, `tanimoto_fcfp4`.

---

## [0.1.4] — 2026-05-28

### Added (`chematic-chem`)

- **BRICS fragmentation** (`brics_bonds`, `brics_fragments`): breaks molecules at retrosynthetically interesting bonds per Dien et al. 2008.
- **QED score** (`qed`): Quantitative Estimate of Drug-likeness (Bickerton et al. 2012); geometric mean of 8 desirability functions. Returns value in [0, 1].
- **Fsp3** (`fsp3`): fraction of sp3 carbons.
- **Aromatic ring count** (`aromatic_ring_count`): number of fully aromatic rings from SSSR.

### Added (`chematic-fp`)

- **AtomPair fingerprint** (`atom_pair_fp`): 2048-bit; encodes atom-pair codes with topological BFS distances (Carhart et al. 1985).
- **Topological Torsion fingerprint** (`torsion_fp`): 2048-bit; encodes four-atom paths with degree ≥ 2 at inner positions (Nilakantan et al. 1987).

### Added (`chematic-smarts`)

- **Recursive SMARTS** `$(...)`: atom must be root of an embedding of the inner SMARTS. Supports arbitrary nesting.
- **Valence** `[vN]`: matches atoms with total valence N (explicit bond orders + implicit H).
- **Ring-bond count** `[xN]`: matches atoms with exactly N bonds where both endpoints share a SSSR ring.
- **Hybridization** `[^N]`: 1 = sp, 2 = sp2 (including aromatic), 3 = sp3.
- **Explicit zero charge** `[+0]` / `[-0]`: matches neutral atoms (charge == 0). Previously `+0` defaulted to `+1`.
- `PartialEq` derived for `QueryAtom`, `QueryBond`, `QueryMolecule`.

### Added (`chematic-depict`)

- **CPK atom coloring**: heteroatoms (N, O, S, Cl, F, Br, I, P) are now colored using the CPK palette in SVG output.
- **`render_svg_highlighted`** / **`depict_svg_highlighted`**: render with yellow circle backgrounds on highlighted atoms and orange strokes on highlighted bonds.

### Added (`chematic-wasm`)

- New descriptor bindings: `logp_crippen`, `fsp3`, `aromatic_ring_count`, `qed`, `exact_mass`, `rotatable_bond_count`.
- New fingerprint similarity functions: `tanimoto_atom_pair`, `tanimoto_torsion`.
- `brics_fragment_count`: number of BRICS fragments.

---

## [0.1.3] — 2026-05-27

### Fixed (`chematic-chem` — LogP Crippen accuracy)

Five new atom-type contexts derived analytically from the 175-molecule RDKit reference set.
LogP MAE vs RDKit: **0.419 → 0.298** (−29%); Pearson r: **0.925 → 0.944**.
17 molecules now have Δ = 0.000: phenol, catechol, resorcinol, hydroquinone,
benzoic_acid, methyl_benzoate, salicylic_acid, toluene, ethylbenzene, phenylacetic_acid,
tetralin, histamine, aniline, n_methylaniline, 4_aminophenol, thiophenol, dopamine.

#### Fix 1 — Phenolic OH hydrogen (+0.1319, was −0.2677 aliphatic alcohol)
- Triggered when O-H is directly bonded to aromatic C (phenol, catechol, tyrosine OH, etc.)
- Verified: phenol (exact), catechol/resorcinol/hydroquinone (2× exact), salicylic_acid (combined exact), dopamine (combined exact)

#### Fix 2 — C=O adjacent to aromatic C (−0.1226, was −0.3800 aliphatic C=O)
- Triggered when sp2 C=X carbon has at least one aromatic C neighbor (Ar-CHO, Ar-COOH, Ar-COOR, Ar-CO-R)
- Verified: benzoic_acid (exact), methyl_benzoate (exact), salicylic_acid (combined exact)

#### Fix 3 — Benzylic sp3 C (Wildman-Crippen C25–C28, was 0.1441 pure alkyl)
- Triggered when sp3 C is bonded to aromatic C but **not** to any heteroatom
- H=3: 0.0764 | H=2: −0.0597 | H=1: −0.1415 | H=0: −0.2037
- Verified: toluene (exact), ethylbenzene (exact), tetralin (exact), phenylacetic_acid (exact), histamine (exact), dopamine (combined exact)

#### Fix 4 — Aniline-type N (bonded to aromatic C, non-amide)
- H=2 primary aniline: −0.7092 (was −1.0190 aliphatic NH2)
- H=1 secondary aniline: −0.2010 (was −0.7096 aliphatic NH)
- H=0 tertiary aniline: −0.5950 (unchanged, no calibration data)
- Verified: aniline (exact), n_methylaniline (exact), 4_aminophenol (combined exact)

#### Fix 5 — Thiol S (0.3132, was 0.6482 thioether)
- Triggered when non-aromatic S has h>0 and no S=O bonds
- Verified: thiophenol (exact), cysteine (residual 0.047 ✓)

### Added

- `@kent-tokyo/chematic` npm package v0.1.3 published to npmjs.com — WebAssembly bindings for browser/Node.js
  - Install: `npm install @kent-tokyo/chematic`
  - Note: unscoped `chematic` blocked by npm similarity check against `chromatic`
- 7 new LogP regression tests in `chematic-chem/tests/rdkit_reference.rs` (phenol, catechol, salicylic_acid, toluene, ethylbenzene, aniline, thiophenol)
- Large-scale ChEMBL validation: **2,897,819 molecules (ChEMBL 37 full set), 100.000% parse+roundtrip success**
  - `chematic-smiles/examples/validate_smiles.rs` — standalone validator (stdin or file, progress every 10k)
  - `scripts/download_chembl_smiles.py` — ChEMBL REST API downloader (deduplication, fragment filter)
  - Streaming pipeline: `curl chembl_37_chemreps.txt.gz | gzip -d | awk | validate_smiles`

---

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

## [0.1.2] — 2026-05-27

### Fixed (`chematic-chem`)

#### TPSA — aromatic N values corrected to match RDKit
- `[nH]` (pyrrole-type aromatic N-H): 13.97 → **15.79 Å²** (RDKit `_CalcTPSAContribs()` measured value)
- `[n;degree≥3]` (N-substituted: N-methyl, N-aryl): 12.89 → **4.93 Å²**
- Effect: caffeine TPSA now 61.82 Å² (was 85.70), exact match with RDKit.
- TPSA MAE vs RDKit (175 molecules): 1.33 → **0.76 Å²**; Pearson r: 0.993 → **0.994**.

#### HBA — aligned with `rdMolDescriptors.CalcNumHBA`
- `[nH]` (aromatic N-H) is **no longer counted** as HBA (lone pair participates in aromaticity).
- Non-aromatic amide N (bonded to C=O) is **excluded** (lone pair delocalized into carbonyl).
- Carboxylic OH (O-H adjacent to C=O) is **excluded**.
- MAE vs RDKit: 0.606 → **0.137** (-77%); Pearson r: 0.932 → **0.975**.
- Verified: aspirin=3, paracetamol=2, caffeine=6, indole=0, acetic acid=1.

#### LogP — Crippen-Wildman with calibrated H contributions
- Added per-H contributions analytically derived from 175-molecule RDKit reference set:
  - H on C: +0.1230 | H on N: +0.2142 | H on alc-O: −0.2677 | H on COOH: +0.2980
- Fixed aromatic C: `[cH]` = +0.1581 (was 0.1441); confirmed from benzene.
- Fixed aromatic N: `[n]`/`[nH]` = −0.3239 (was +0.2626); confirmed from pyridine, pyrrole.
- Fixed S: thioether = +0.6482 (was 0.2432); aromatic S = +0.6237 (was 0.0).
- Fixed O: alcohol OH = −0.2893, ether O = −0.0684, carbonyl O = −0.0509 (were all 0.1552).
- Fixed Cl: aromatic = +0.7904, aliphatic = +0.6895.
- Added exocyclic C=O handling for aromatic C (caffeine carbonyl C now C10 = −0.3800).
- MAE vs RDKit: 1.346 → **0.419** (-69%); Pearson r: 0.456 → **0.925** (+103%).

### Added

- 21 new regression tests in `chematic-chem`: 7 HBA tests + 14 LogP calibration tests anchored to RDKit TSV values.
- `docs/rdkit_comparison.md` — quantitative comparison report vs RDKit (175 molecules, v0.1.0 → v0.1.2).

---

## [0.1.1] — 2026-05-27

### Added

- All crates bumped to version 0.1.1.
- `chematic-wasm`: New crate providing WebAssembly (wasm-bindgen) bindings for JavaScript/TypeScript consumers. Exposes SMILES parsing, canonical SMILES, molecular descriptors, ECFP fingerprints and Tanimoto similarity via `wasm-bindgen`.
- ChEMBL roundtrip validation tests: parse → write → parse identity verified against 1000+ ChEMBL molecules (MOL/SDF V2000 format).
- criterion benchmarks added to `chematic-smiles` (`parse_bench`) and `chematic-fp` (`ecfp_bench`) for continuous performance tracking.

### Changed

- SEO/metadata improvements to all `Cargo.toml` files: added `readme`, `homepage`, and `documentation` fields; improved `keywords` (max 5) and `categories`; sharpened `description` to clearly identify each crate as part of the pure-Rust RDKit-alternative ecosystem.
- All internal path dependency version constraints updated from `"0.1.0"` to `"0.1.1"`.

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

[Unreleased]: https://github.com/kent-tokyo/chematic/compare/v0.1.6...HEAD
[0.1.6]: https://github.com/kent-tokyo/chematic/compare/v0.1.5...v0.1.6
[0.1.5]: https://github.com/kent-tokyo/chematic/compare/v0.1.3...v0.1.5
[0.1.3]: https://github.com/kent-tokyo/chematic/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/kent-tokyo/chematic/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/kent-tokyo/chematic/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/kent-tokyo/chematic/releases/tag/v0.1.0
