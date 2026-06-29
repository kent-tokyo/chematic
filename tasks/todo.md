# chematic — Status & Roadmap

Current version: **v0.4.24** (2026-06-29)

---

## Crate map

| Crate | Role | Tests |
|-------|------|-------|
| `chematic-core` | Atom, Bond, Molecule, Element, kekulization (4-pass incl. Edmonds' blossom + boron fix) | 69 |
| `chematic-smiles` | OpenSMILES parser/writer, canonical SMILES | 48 |
| `chematic-perception` | SSSR, 2-pass Hückel aromaticity, CIP stereo, `count_aromatic_rings` | 34 |
| `chematic-smarts` | SMARTS, VF2 subgraph, MCS (McGregor), LRU SMARTS cache; **atom map `:N`** in SMARTS (`[O;D1;H0:3]`); **`find_matches_with_rings`** — shared SSSR across multi-pattern batches | 142 |
| `chematic-chem` | **190+ descriptors**, ADMET, BOILED-Egg, QED, SA Score, PAINS/Brenk, pKa, ESOL; Schultz/Gutman MTI, VABC, Gravitational index; **HBD 100% / TPSA ±0.1 Å² / LogP ±0.3** RDKit parity; **`logp_and_mr`** + `chi_all` + `pains_passes_and_matches` perf APIs | 639 |
| `chematic-fp` | ECFP/FCFP, MACCS, MAP4, AtomPair, Torsion, MHFP, ERG, Tanimoto | 87 |
| `chematic-ff` | MMFF94 full stack (7 terms), DREIDING, L-BFGS minimizer | 51 |
| `chematic-3d` | ETKDG (80 torsion rules, chair/envelope ring conf), WHIM (22-dim), GETAWAY (19-dim), RDF (20-dim), SASA, USR shape screen | 45 |
| `chematic-depict` | 2D SVG, grid rendering, **PDF (`pdf` feature)**, **EPS (pure Rust)** | 34 |
| `chematic-rxn` | Reaction SMILES/SMIRKS, `run_reactants`/`run_reactants_strict`, RECAP/BRICS; **`retro_disconnect()` — 60 retro-SMIRKS** (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ranking; **parity-aware `@`/`@@` SMIRKS stereo filtering** | 25 |
| `chematic-inchi` | InChI/InChIKey: pure-Rust approx (inline SHA-256, no sha2 dep) + IUPAC-exact (`native-inchi` feature, v1.07.5) | 28+16* |
| `chematic-iupac` | IUPAC name generation, 25+ compound classes | 45 |
| `chematic-mcp` | MCP server, **15 tools** (JSON-RPC 2.0 over stdio, `name_to_smiles` via PubChem) | 28 |
| `chematic-mol` | SDF/MOL V2000/V3000, CML, CDXML, **ChemicalJSON (.cjson)** | 40 |
| `chematic-wasm` | 160 WASM exports, npm `@kent-tokyo/chematic` (**504 KB gzip**, -38.5%) | 209 |
| `chematic-py` | PyO3 Python bindings (`pip install chematic`); Sprint 18–26+: 300+ API endpoints | 300+ |
| `chematic-ewald` | PME Ewald summation, B-spline interpolation | 12 |

`cargo test --workspace --lib --quiet` → **2319 tests** (lib only), all passing

---

## MCP tools (chematic-mcp, 15 total)

`parse_smiles` · `canonical_smiles` · `calc_properties` · `lipinski_check` · `sa_score` · `pains_check` · `brenk_check` · `admet_profile` · `boiled_egg` · `ecfp4` · `tanimoto` · `smarts_match` · `find_mcs` · `generate_3d` · `name_to_smiles`

---

## Known limitations

| Item | Status |
|------|--------|
| Kekulization failures | **1/5000** — only pure H₂ `[H][H]` (no heavy atoms; IUPAC InChI library constraint, not a kekulization issue) |
| Aromatic ring count vs RDKit | **~100%** (222/222 bench5k failures fixed in v0.4.11 — `augmented_ring_set` XOR guard `min`→`max`) |
| Bridgehead atoms vs RDKit | **100%** (4999/4999 — bond-intersection algorithm, 2026-06-28) |
| HBA agreement vs RDKit | 99.98% (4999/5000) |
| HBD agreement vs RDKit | **100%** (175/175 TSV bulk test; S-H thiol fix in v0.4.12+) |
| TPSA accuracy vs RDKit | **±1.0 Å²** (175-mol bulk test; nitro-N, oxide bridge, Kekulé-N fixes) |
| LogP accuracy vs RDKit | **±0.3** (175-mol bulk test; oxide bridge O fix) |
| InChI E/Z `/b` layer | done: implemented (v0.4.5) |
| Kekulization blossom | done: implemented (v0.4.5, 128→2 failures) |
| Boron aromatic kekulization | done: fixed (v0.4.7, 2→1 failure) |
| BOILED-Egg | done: implemented (v0.4.5 Rust, v0.4.6 Python, v0.4.7 WASM) |

---

## Out of scope (constraints)

- Kekulization: all corpus failures resolved except pure H₂ (IUPAC InChI library requires at least one heavy atom)
- Full ETKDG stochastic sampling (requires ML distance geometry)
- Transition metals / coordination compounds
- ML-based property prediction
- HELM / FASTA (peptides/proteins)
- IBM RXN4Chemistry / commercial APIs

---

## WASM バンドルサイズ削減メモ (2026-06-23計測)

**最終結果 (2026-06-23): 2156 KB raw / 819 KB gzip → 1309 KB raw / 510 KB gzip (−38%)**

| 対策 | 効果 | 状況 |
|------|------|------|
| `tiny_skia` を optional `png` feature に移動、WASM で無効化 | −220 KB raw / −80 KB gzip | ✅ done |
| `sha2` クレートをインライン SHA-256 (60行) に差し替え | ~−15 KB gzip | ✅ done |
| `[profile.release] opt-level="z" lto=true codegen-units=1` | −541 KB raw / **−172 KB gzip** (最大効果) | ✅ done |
| `wasm-opt -O3` を CI (pages.yml) に統合 | ~−5 KB gzip | ✅ done |
| `run_md_json` / `coulomb_energy_json` / `torsion_scan_json` / `determine_bonds_from_xyz_json` 除去 + `chematic-ewald` 依存削除 | −5 KB gzip | ✅ done |
| PAINS/Brenk SMARTS 圧縮 | **逆効果** (raw −26 KB だが gzip +4 KB) — 差し戻し | ❌ not worth it |

**現状**: 819 KB → **504 KB gzip** (−38.5%)

**次フェーズ候補**: 504 → ~450 KB gzip
- コード削減: 使われていない WASM export の体系的な分析
- IUPAC 命名モジュール削除オプション（ブラウザ用 lite ビルド）

---

## Next candidates

| Priority | Item |
|----------|------|
| done: | Name→SMILES via PubChem REST proxy — `name_to_smiles` tool added (v0.4.8) |
| done: | AutoDock PDBQT format (parse_pdbqt / write_pdbqt / autodock_atom_type) — chematic-mol v0.4.9 |
| done: | UFF force field (assign_uff_types / uff_total_energy / minimize_uff) — chematic-ff v0.4.9 |
| done: | SDF partial charge writing (write_sdf_with_charges) — chematic-mol v0.4.9 |
| done: | Sprint 18–26: 50+ new Python API endpoints (PyPI gap analysis) — tanimoto_matrix, ring_families, stereo_from_coords/2d_coords, from_cxsmiles, from_rxn_file/to_rxn_file, parse_sdf_with_coords, etc. |
| done: | MAP4 Python binding (mol.map4, bulk.map4, tanimoto_map4) — Python gap vs chemfp/mordred |
| done: | LogD(pH) Python binding (mol.logd, mol.logd_profile) — ADMET key descriptor |
| done: | MQN descriptors Python binding (mol.mqn) — 42-element Ertl 2009 |
| done: | Butina clustering + MaxMin diversity Python binding — closes gap vs chemfp |
| done: | generate_3d / conformer_ensemble / WHIM / GETAWAY / PDB-XYZ I/O Python bindings |
| done: | mmff94_energy_breakdown + mmff94_torsion_scan Python bindings |
| done: | functional_groups + scaffold_network Python bindings |
| Medium | **JOSS paper** (Journal of Open Source Software) — 提出目標: **2026-11-20**。必要: `paper.md`, `CITATION.cff`, `LICENSE-MIT`/`LICENSE-APACHE` テキストファイル, `CODE_OF_CONDUCT.md`, Zenodo DOI |
| done: | Aromatic ring count improvement — `augmented_ring_set` XOR guard `min`→`max` fix (v0.4.11); 95.6% → ~100% RDKit agreement (222件全修正) |
| done: | MMFF94 BCI precision — already at ±0.05e (better than ±0.1e target); no action needed |
| done: | LogP alkenyl C — already implemented (terminal 0.1551, aryl-adjacent 0.2640); no action needed |
| done: | `retro_disconnect()` — 60 retro-SMIRKS templates (6 reaction classes) + SA Score ranking; Python `mol.retro_disconnect(reaction_class=...)` |
| done: | TPSA/LogP/HBD descriptor accuracy: nitro-N fix, oxide bridge fix, Kekulé-N fix, S-H HBD fix; 175-mol bulk regression tests (±1.0 Å² / ±0.3 / exact) |
| done: | bench5k.py extended with TPSA, LogP, HBD comparison vs RDKit |
| done: | OSS credibility sprint (2026-06-27): README 3-row badges, SECURITY.md v→v0.4.23, pyproject.toml 190+ descriptors + Python 3.13, security.yml hardened (no continue-on-error, dtolnay, cargo-deny job), deny.toml, codeql.yml (deleted — conflicts with Default setup), scripts/check.sh, scripts/bump_version.py, CLAUDE.md release flow |
| done: | Bridgehead atoms 98.5% → **100%** (2026-06-28): bond-intersection algorithm — an atom is bridgehead iff some ring pair shares ≥ 2 bonds and the atom is incident to exactly 1; fixes both over-count (peroxide cages) and under-count (spiro+bridgehead co-occurrence) |
| done: | `chematic.screen(smiles, profile="druglike")` — bundles lipinski/veber/pains/brenk/qed/sa_score into one call; profiles: druglike/fragment/leadlike; returns list[dict] with per-filter _pass flags and overall_pass |
| done: | `scripts/analyze_logp_mismatches.py` — compare chematic vs RDKit Crippen LogP, output buckets + TSV; confirmed 0 mismatches (100% at ±0.01 on 5k corpus) |
| done: | `chematic.rdkit_compat` Sprint 2 (v0.4.24): `ExplicitBitVect` (GetNumBits/GetBit/SetBit/GetOnBits/ToBitString/IndexError); `DataStructs.BulkTanimotoSimilarity`; `GetMorganFingerprintAsBitVect` returns `ExplicitBitVect`; unsupported opts fail loudly (useFeatures/bitInfo/nBits≠2048/unknown kwargs) |
| done: | `chematic.rdkit_compat` Sprint 3 (v0.4.24): Rust `atom_table`+`bond_table` getters; `BondType` constants; `Atom` wrapper (GetIdx/Symbol/AtomicNum/FormalCharge/IsAromatic/Degree/TotalDegree/TotalNumHs/IsInRing/IsInRingSize); `Bond` wrapper (GetIdx/BeginAtomIdx/EndAtomIdx/BeginAtom/EndAtom/BondType/BondTypeAsDouble/IsAromatic/OtherAtomIdx); `Mol` GetNumBonds/GetAtoms/GetBonds/GetAtomWithIdx/GetBondWithIdx |
| done: | `chematic.rdkit_compat` Sprint 4 (v0.4.24): Rust `sssr_atom_rings` getter; `RingInfo` (NumRings/AtomRings/BondRings/NumAtomRings/NumBondRings); `Mol.GetRingInfo()` with lazy cache; `Bond.IsInRing()` implemented; `Atom.IsInRingSize(n)` uses ring info cache |
| done: | SDF true streaming: `SdfFileReader<R: BufRead>` + Python `iter_sdf(path)` → lazy streaming, `iter_sdf_batched(path, batch_size=1000)` new; `iter_sdf_str()` unchanged; `MolParseError::Io` added |
| done: | `bulk.descriptors_array(smiles, columns)` — columnar numpy output, ~25% faster than `descriptors()+DataFrame()`; float64/bool/NaN for optional |
| done: | Stereocenters 99.8% → 99.98% (legacy) / 98.7% (new CIP) — CIP Rule 5 tie-breaking via provisional R/S + equality-based signature comparison; cage false-positives correctly avoided |

---

## Version history (recent)

| Version | Date | Highlights |
|---------|------|-----------|
| v0.4.24 | 2026-06-29 | `chematic.rdkit_compat` Sprints 2-4: ExplicitBitVect, DataStructs, BulkTanimoto, Mol/Atom/Bond surface, RingInfo; RDKit-compatible SDF I/O (SDMolSupplier/SDWriter + property roundtrip) |
| v0.4.23 | 2026-06-28 | `screen()`+LogP analyzer+SDF streaming+`descriptors_array()`; bridgehead/spiro/rotatable 100% |
| v0.4.23 | 2026-06-26 | LogP 96.5% → 99.7% (symmetric triple bond VF2 dedup fix); OSS sprint: badges, CI hardening, cargo-deny, check.sh, bump_version.py |
| v0.4.19 | 2026-06-23 | PDF/EPS 出力、ChemicalJSON、Schultz/Gutman MTI・VABC・Gravitational index、bulk.substructure_match、bulk.generate_3d、bulk.tanimoto_matrix、bulk.standardize、inline SHA-256 (sha2 dep 除去)、WASM 819→504 KB gzip (-38.5%) |
| v0.4.18 | 2026-06-23 | Jupyter `_repr_svg_`、`from_smiles_list`、`descriptors_df`、`chi_all`、`pains_passes_and_matches`、CNS MPO perf |
| v0.4.17 | 2026-06-23 | PAINS/Brenk dedup perf sprint |
| v0.4.16 | 2026-06-22 | **Perf**: shared SSSR in SMARTS matching (117→1 per Crippen, ~480→1 per PAINS, ~300→1 per BRENK); `logp_and_mr()` combined Crippen pass; `logd_from_logp()` helper; `cns_mpo_score` logP dedup; `eccentric_connectivity_index` reuses `graph_eccentricities`; `heavy_degrees()` pre-comp in randic/zagreb. New public API: `find_matches_with_rings`, `find_matches_with_rings_and_config`, `logp_and_mr`, `logd_from_logp`. CI: setup-python v6, upload-artifact v7 |
| v0.4.15 | 2026-06-21 | Tautomer tetrazole 1H/2H normalization — BFS 1,2-shift + canonical SMILES tiebreaker; CDXML Order=1.5→Aromatic |
| v0.4.13 | 2026-06-21 | HBD S-H fix; TPSA nitro-N / oxide bridge / Kekulé-N fixes; LogP oxide bridge fix; `retro_disconnect()` 60 retro-SMIRKS; ETKDG 40 torsion patterns; bulk TPSA ±1.0/LogP ±0.3/HBD 100% |
| v0.4.12 | 2026-06-21 | SMARTS atom map `:N` support (`[O;D1;H0:3]`); fix aromatic-bond false MapNumberMismatch; fix `[C:]` parse error; propagate atom_map in mol_to_query |
| v0.4.11 | 2026-06-21 | Aromatic ring count ~100% RDKit parity (222 bench5k fixes); CIF/Gaussian parser 8 safety fixes; clippy CI fixes |
| v0.4.10 | 2026-06-20 | Sprint 18–26: 50+ new Python bindings (PyPI gap analysis p.12–p.20): tanimoto_matrix, ring_families, stereo_from_coords, CXSMILES, RXN file I/O, 2D stereo, SDF batch coords, DREIDING minimize, etc. |
| v0.4.9 | 2026-06-19 | AutoDock PDBQT format, UFF force field, SDF partial charge writing; Python+WASM bindings |
| v0.4.8 | 2026-06-19 | Iterative `augmented_ring_set`, `name_to_smiles` MCP tool (PubChem proxy) |
| v0.4.7 | 2026-06-19 | Boron aromatic kekulization fix (2→1 failure), WASM `admet_profile_json` + BOILED-Egg |
| v0.4.6 | 2026-06-19 | `boiled_egg()` Python binding, `admet()` extended, WASM `boiled_egg_json()`, stubs updated |
| v0.4.5 | 2026-06-19 | Kekulization blossom (128→2), InChI E/Z `/b` layer, 6 new MCP tools, BOILED-Egg |
| v0.4.4 | 2026-06-18 | (skipped tags v0.4.2–v0.4.3) |
| v0.4.1 | 2026-06-18 | `aromatic_ring_count` fix + HBA rewrite — closes issue #12 |
| v0.4.0 | 2026-06-17 | native-inchi (IUPAC C lib 1.07.5), Python PyO3 bindings |
| v0.3.2 | 2026-06-15 | criterion benchmarks |
| v0.3.1 | 2026-06-15 | WASM pKa/ADMET bindings |
| v0.3.0 | 2026-06-15 | MCP server (8 tools), pKa, ADMET |
| v0.2.11 | 2026-06-14 | MMFF94 OOP+STRE-BEN, MAP4, SMARTS cache |
| v0.2.10 | 2026-06-14 | L-BFGS, MMFF94 energy breakdown, torsion scan |
| v0.2.9 | 2026-06-14 | MMFF94 full minimizer (bond/angle/torsion/vdW/elec) |
| v0.2.7-8 | 2026-06-14 | MMFF94 charges + energy parameters complete |

For detailed sprint history see `git log --oneline`.
