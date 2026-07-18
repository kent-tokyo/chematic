# chematic — Status & Roadmap

Current version: **v0.4.30** (2026-07-17)

---

## Crate map

| Crate | Role | Tests |
|-------|------|-------|
| `chematic-core` | Atom, Bond, Molecule, Element, kekulization (4-pass incl. Edmonds' blossom + boron fix) | 71 |
| `chematic-smiles` | OpenSMILES parser/writer, canonical SMILES | 109 |
| `chematic-perception` | SSSR, 2-pass Hückel aromaticity, CIP stereo, `count_aromatic_rings` | 101 |
| `chematic-smarts` | SMARTS, VF2 subgraph, MCS (McGregor), LRU SMARTS cache; **atom map `:N`** in SMARTS (`[O;D1;H0:3]`); **`find_matches_with_rings`** — shared SSSR across multi-pattern batches | 142 |
| `chematic-chem` | **190+ descriptor values (71 functions)**, ADMET, BOILED-Egg, QED, SA Score, PAINS/Brenk, pKa, ESOL; Schultz/Gutman MTI, VABC, Gravitational index; **HBD 100% / TPSA ±0.1 Å² / LogP ±0.3** RDKit parity; **`logp_and_mr`** + `chi_all` + `pains_passes_and_matches` perf APIs | 662 |
| `chematic-fp` | ECFP/FCFP, MACCS, MAP4, AtomPair, Torsion, MHFP, ERG, Tanimoto | 185 |
| `chematic-ff` | MMFF94 full stack (7 terms), DREIDING, L-BFGS minimizer | 98 |
| `chematic-3d` | ETKDG (80 torsion rules, chair/envelope ring conf), WHIM (22-dim), GETAWAY (19-dim), RDF (20-dim), SASA, USR shape screen | 265 |
| `chematic-depict` | 2D SVG, grid rendering, **PDF (`pdf` feature)**, **EPS (pure Rust)** | 64 |
| `chematic-rxn` | Reaction SMILES/SMIRKS, `run_reactants`/`run_reactants_strict`, RECAP/BRICS; **`retro_disconnect()` — 60 retro-SMIRKS** (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ranking; **parity-aware `@`/`@@` SMIRKS stereo filtering** | 137 |
| `chematic-inchi` | InChI/InChIKey: pure-Rust approx (inline SHA-256, no sha2 dep) + IUPAC-exact (`native-inchi` feature, v1.07.5) | 96+16* |
| `chematic-iupac` | IUPAC name generation, 25+ compound classes | 47 |
| `chematic-mcp` | MCP server, **20 tools** (JSON-RPC 2.0 over stdio, `name_to_smiles` via PubChem) | 31 |
| `chematic-mol` | SDF/MOL V2000/V3000, CML, CDXML, **ChemicalJSON (.cjson)** | 130 |
| `chematic-wasm` | 160 WASM exports, npm `@kent-tokyo/chematic` (**719 KB gzip** as of v0.4.29, 2026-07-17 — up from 504 KB, see `benchmarks/2026-07-17.md`) | 211 |
| `chematic-py` | PyO3 Python bindings (`pip install chematic`); Sprint 18–26+: 300+ API endpoints | 300+ |
| `chematic-ewald` | PME Ewald summation, B-spline interpolation | 16 |

`cargo test --workspace --lib --quiet` → **2,366 tests** (lib only), all passing

---

## MCP tools (chematic-mcp, 20 total)

`parse_smiles` · `canonical_smiles` · `calc_properties` · `lipinski_check` · `sa_score` · `pains_check` · `brenk_check` · `admet_profile` · `boiled_egg` · `ecfp4` · `tanimoto` · `smarts_match` · `find_mcs` · `generate_3d` · `name_to_smiles` · `retrosynthesis` · `smiles_to_moljson` · `moljson_to_smiles` · `representation_router` · `molecule_context_pack`

---

## Known limitations

| Item | Status |
|------|--------|
| PyPI/npm publish lag behind crates.io | **Known gap**: git tags `v0.4.23`–`v0.4.26` were never created/pushed to origin (only the release commits exist), so `publish-pypi.yml` / `publish-npm.yml` / `release.yml` never fired for them. crates.io stays current via a separate manual `cargo publish`. Fix: push the missing tags (`git tag vX.Y.Z && git push origin vX.Y.Z`) when ready to catch PyPI/npm/GitHub Releases up — this is a real external-publish action, done deliberately by a maintainer, not automatically. |
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
| done: | `chematic-cip` new crate, Milestone 1 → 3B closeout: hierarchical-digraph CIP engine (`assign_cip_accurate_experimental`), sphere-by-sphere Rules 1a/1b/2 comparator, RDKit-compatible MANCUDE fractional atomic numbers; full-corpus accuracy on this experimental engine 96.68% → 99.14% (4047/4186 → 4150/4186), 0 regressions; not yet wired into `chematic_chem::assign_cip()` — see `docs/cip_accurate_rfc.md` |
| done: | `chematic-cip` perf fix: `ring_bond_set` was calling full `find_sssr` for a boolean ring-bond check (30ms vs 27us on a 182-atom molecule, ~1000x); replaced with an O(V+E) bridge-edge DFS, byte-identical output verified |
| done: | CI: Criterion regression gate bootstrap fix (PR #69) — a new/removed bench target no longer aborts the whole gate job |
| found, not yet fixed: | Criterion regression gate pseudo-replication — ~100 samples/side from one process aren't independent trials, single-run environment differences amplify into false-looking uniform failures across unrelated benchmarks; gate is non-required and its fail verdicts aren't currently trusted; redesign tracked in issue #70 |
| done: | v0.4.30 release: chematic-cip published to crates.io (was `publish = false`, blocked the whole 18-crate publish graph), `scripts/check_publish_graph.py` added (publish-graph + dev-dependency-deadlock checks, wired into `scripts/check.sh`); full crates.io/PyPI/npm/GitHub Release cycle completed, fresh-install smoke-tested on all 3 channels |
| done: | CIP-Perf-A0 (PR #106, diagnosis only): measured whether chematic has RDKit PR #9171's bug (auxiliary/Rule-4b-equivalent descriptors computed for every stereocenter instead of only tied ones) — it doesn't: 99.3% of stereocenters resolve at Rules 1a/1b/2 alone, only 0.7% ever reach Rule 4b/5 (already gated on `SkipReason::Tied`). CIP-Perf-S1 closed as unnecessary. Surfaced a separate real cost tail in `rank_children` on symmetric/duplicate-heavy digraphs (median 24 comparisons, p95 954, max 89,250) — tracked as issue #107 (CIP-Perf-A1), not started |
| done: | MANCUDE-Decision-A0 (PR #106, diagnosis only): `CompareContext::fractional_decisions` found nonzero (248) on the same frozen corpus Milestone 3B-1b's closeout claimed would be 0 — that claim was only ever checked against one curated molecule, never at full-corpus scale. Classified by isolating Kekulé-respelling structure from the fractional atomic-number values (first attempt conflated both and wrongly got 3 rows as label-changing; corrected): all 36 affected stereocenters are D (fraction locally load-bearing, final-label-inert), 0 are E. Milestone 3B-2 closed by diagnosis, no corrective implementation needed; `compare.rs`/`docs/cip_accurate_rfc.md` updated |
| done: | EZ-S1 — fixed the EZ-A0 gap directly (no separate diagnosis PR): `substituent_is_up` (`crates/chematic-chem/src/cip.rs`, legacy engine) now reads `mol.bond_direction(bond_idx)` before falling back to `bond.order`, same pattern `canonical.rs` already used. Reading the stash exposed two further, independent, pre-existing bugs in `highest_stereo_sub`: (a) it returned the highest-CIP-priority substituent *among those carrying an explicit `/`/`\` marker*, not the true highest-priority substituent overall, silently using a lower-priority marked substituent's raw side instead of the true-highest one's geometric complement (affects all trisubstituted alkenes, aromatic or not, e.g. non-aromatic `Cl/C(F)=C(\Br)I`); (b) no tie guard — a genuine CIP priority tie between an alkene end's only two substituents (e.g. the two ring branches of an *unsubstituted*, symmetric ring's ipso carbon; confirmed via `compare_branches(...) == Ordering::Equal` in both argument orders) was silently assigned an arbitrary, adjacency-order-dependent E/Z instead of correctly reporting no stereogenic bond, flipping under `canonical_smiles` renumbering. All three fixed together after user sign-off on expanding scope past the original "direction-read only" ask, then a follow-up review caught two measurement gaps in the first verification pass (see PR #108 review): the original 678/678 figure only covered bonds chematic assigned *something* to, silently excluding any bond chematic missed entirely, and the "0 allene changes" claim used an `"=C=" in smiles` string heuristic instead of structural classification. `scripts/ez_stash_gap_report.py` rewritten oracle-first (scans every corpus SMILES via RDKit directly, not just chematic's own output) and `crates/chematic-chem/examples/ez_stash_gap_snapshot.rs` now tags each row's `kind` (`tetra`/`ez`/`allene`/`parse_fail`) structurally, mirroring `is_allene_central`'s shape check rather than string-matching. Full 5,000-molecule corpus: `rdkit_ez_total`=678, `candidate_correct`=678, `candidate_wrong`=`candidate_missing`=`candidate_extra`=0 (baseline: 492/169/17/0) — full completeness proven, not just correctness-where-attempted. `ez` newly=17/lost=0/flipped=169 (all RDKit-agree), `allene` newly=lost=flipped=0 (structural classification), `tetra` (R/S) changed=0. Round-trip audit (`canonical_stereo_d0_roundtrip_audit.rs`): 5,000/5,000 stable both before and after |
| done: | EZ-S1 sibling gap fixed — `crates/chematic-inchi/src/native/convert.rs`'s `find_stereo_sub` closure (native-inchi feature, `standard_inchi()`) mirrored the pre-fix `substituent_is_up` exactly: read only `nb_bond.order`, never `mol.bond_direction`, for InChI `Stereo0D` double-bond descriptors. Only the stash-read fix applies here (`mol.bond_direction(bond_idx).unwrap_or(nb_bond.order)`, same pattern) — the CIP-priority-fallback/tie-guard bugs found alongside `substituent_is_up` don't have an analogue: InChI's Stereo0D format doesn't need a CIP-priority-highest substituent, any one determinate substituent per end is sufficient (parity is defined relative to whichever neighbor is fed in). `apply_kekule` preserves the stash verbatim across the mandatory pre-InChI kekulization, so the gap applied post-kekulization too. Full 5,000-molecule corpus: 12/5,000 changed, all pure `/b`-layer gains (0 lost, 0 changed, 0 other side effects), all 12 match RDKit `Chem.MolToInchi` byte-for-byte. New `crates/chematic-inchi/tests/ez_stash_stereo0d.rs` (3 tests, RDKit-pinned specific values, incl. a mobile-H-tautomer negative control confirming InChI legitimately drops stereo for the exact repro-class ring when a tautomeric N-H is present) |
| done: | ECFP RDKit atom-invariant parity, PR 1 (chematic-fp core only): new `EcfpInvariantMode` enum (`Chematic` default / `RdkitMorgan`) threaded through a new `ecfp_with_invariant_mode`/`ecfp_with_bitinfo_and_mode` shared core, with `ecfp`/`ecfp_with_bitinfo`/`morgan_fp_counts` unchanged thin wrappers — `EcfpConfig` itself untouched (public-struct-literal breakage risk avoided) — plus `ecfp4_rdkit_invariants`/`ecfp6_rdkit_invariants`/`atom_invariants` convenience API. `RdkitMorgan` mode invariant (heavy-atom degree, total H count, formal charge, ring membership, isotope delta — no aromaticity byte) empirically reverse-engineered against RDKit 2026.03.3's `GetConnectivityInvariants` via `SmilesParserParams(removeHs=False)`-preserved explicit-H probes, not assumed from ECFP literature. Gates: (1) non-regression — full 5,000-mol corpus, `ecfp4`/`ecfp6`/`ecfp4`-chiral/`ecfp_with_bitinfo` (fp+origins)/`morgan_fp_counts` MD5-identical before/after: 0 bytes changed. (2) atom-invariant partition agreement (hash-value-independent, same methodology as the existing ECFP4-vs-RDKit measurement) — **100.0000%** on the 5,000-mol corpus and the new 20-case `scripts/ecfp_rdkit_edge_fixtures.csv` (isotopes, explicit H, aromatic/Kekulé pairs, charged aromatics, fused/bridged rings). One classified (not unclassified) residual found and kept, not chased further: chematic's `isotope - round(atomic_mass)` delta formula doesn't reproduce every one of RDKit's undocumented isotope-delta rounding cases (e.g. carbon-13 vs carbon-12— real corpora carry 0 isotope labels), pinned as `rdkit_isotope_delta_known_gap_vs_rdkit`. Reference (non-gating): Tanimoto-correlation vs RDKit r=0.9428 (≈ the existing Chematic-mode r=0.94 baseline — expected, since aromaticity rarely flips overall pairwise similarity ranking even though it changes exact bit agreement). Named `RdkitMorgan` not `RdkitCompatible` (atom-invariant parity only, not RDKit hash/folding/environment-dedup compatibility — that name is reserved for a possible future full-compat milestone). PR 2 (Python/WASM bindings) not started |

---

## Version history (recent)

| Version | Date | Highlights |
|---------|------|-----------|
| v0.4.30 | 2026-07-17 | chematic-cip published to crates.io, all 18 crates live; CIP-Perf-A0/MANCUDE-Decision-A0 diagnosis closeout (see rows above); README badge cleanup. Versions v0.4.25–v0.4.29 not individually detailed here — see `CHANGELOG.md`/`git log` for the full PR history in between |
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
