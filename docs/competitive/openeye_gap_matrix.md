# OpenEye Gap Matrix

Method: every row is grounded in chematic's actual source, tests, and
measured results as of `main` at 2026-08-23 (v0.19.0). No row is based on
README/marketing text alone. OpenEye's own behavior is described only from
publicly documented capability (product pages, published papers) — this
project has not obtained, run, or reverse-engineered any OpenEye binary,
and does not compare against OpenEye's actual output. Comparisons are
therefore capability-shape comparisons, not head-to-head benchmark numbers,
except where chematic's own gap (e.g. "no Gaussian-volume shape code
exists at all") is unambiguous regardless of what OpenEye's numbers are.

Classification: **Ahead** / **Competitive** / **Partial** / **Missing** /
**Not planned**.

Full technical detail behind every row: `docs/rfcs/openeye_materials_advantage_rfc.md`.

## OEChem (molecule model, SMILES/SMARTS, formats, standardization)

| Capability | Status | Evidence |
|---|---|---|
| Core molecule model, SMILES/SMARTS parse+write | Competitive | `chematic-smiles`/`chematic-smarts`; RDKit-differential-validated elsewhere in this project (not re-audited here — out of this RFC's scope, already covered by the 100-point ladder). |
| File format breadth | Partial | ~20 formats total (per `docs/rdkit-comparison.md` line 143's own self-reported figure) vs. OpenEye/OEChem's broad multi-format support and OpenBabel's 146 (`ROADMAP.md` line 44). Materials/simulation formats (CIF, POSCAR/CONTCAR, LAMMPS, QCSchema, ORCA, Cube, OpenDX) are a deliberate, differentiated strength (`docs/format-capabilities.md`, 15-format materials matrix) — breadth is intentionally not chased for its own sake (`ROADMAP.md` line 45-47). |
| Standardization / Parent identity | Competitive | `crates/chematic-chem/src/standardize.rs`: `FragmentPolicy`+`select_fragment` (line 285/470), `fragment_parent`/`charge_parent`/`isotope_parent`/`stereo_parent` (638/655/668/689), `tautomer_parent`/`super_parent` (`tautomer.rs:1753`, `parent.rs:118`). Explainable, typed audit trail (`ParentResult`/`ParentComputationStatus`/`ParentAudit`) — this is a genuine differentiator, not just parity. |

## Quacpac (protonation, tautomer, microstates, partial charges)

| Capability | Status | Evidence |
|---|---|---|
| Tautomer canonicalization | Partial | Aromatic lactam/lactim (2-pyridone/4-pyridone/uracil) fixed (round 2C, PR #365). Ring-N-H non-uniqueness (cytosine/guanine/hypoxanthine, §1.7) and nitroso/oxime (§1.6) confirmed open — `docs/rfcs/tautomer_parent_identity_phase2_rfc.md`. `TautomerScoringConfig` (customizable scoring) is scoped, not implemented (0 occurrences in `crates/*/src`). |
| pKa prediction | Partial | `crates/chematic-chem/src/pka.rs`: 23-entry static SMARTS rule table (6 acid + 17 base), returns the single most-acidic/most-basic site only (`reduce(f64::min/max)`) — a scalar prediction, not distribution/enumeration. Doc comment states accuracy "±1-2 pKa units," suitable for triage not lead optimization. Note: `docs/rdkit-comparison.md` says "15 SMARTS rules" — stale vs. the code's actual 23 (a documentation-freshness finding, not a capability gap; worth a separate small fix outside this RFC's scope). |
| Protonation/ionization microstate enumeration | **Missing** | Workspace-wide search for "protomer"/"microstate"/"ionization state" returns zero enumeration code. `reionize`/`uncharge` (`standardize.rs:1235`/`1306`) each apply one fixed heuristic and return exactly one molecule — never a candidate set. This is the single cleanest, most complete Quacpac gap found in this audit: OpenEye's whole microstate-enumeration-plus-ranking concept has no chematic analog at all. |
| Tautomer/protomer dominant-state ranking | Partial | Tautomer scoring exists (`tautomer_score`, position-blind by construction per round 2C-N's own finding) and ranks *tautomers*. No equivalent ranking exists across *ionization* states (there is nothing to rank, since no enumeration exists). |
| Partial charges | Competitive | Gasteiger-Marsili PEOE (`gasteiger.rs`, RDKit-parameter-derived), MMFF94 BCI charges (`mmff94_bci.rs` in both `chematic-chem` and `chematic-ff` — noted duplication, unreconciled, a code-hygiene item outside this RFC's scope). No AM1-BCC/RESP/EEM equivalent (a real but lower-priority gap — AM1-BCC requires semiempirical QM infrastructure chematic doesn't have and isn't planned to build). |

## Omega (conformer ensemble generation, torsion drive, macrocycle)

| Capability | Status | Evidence |
|---|---|---|
| Single-conformer 3D generation | Competitive | `embed_pipeline_v2` (`pipeline_v2.rs:600`): distance geometry + torsion-knowledge (82 tests across `etkdg_knowledge/*.rs`) + stereo verify/repair + policy-gated force field, fully typed 12-stage failure enum, explicit deterministic `random_seed: u64` (tested: `same_seed_reproducible`/`different_seed_gives_different_output`). Last measured `pipeline_v2_mmff94_strict` = 241/265 on the project's own 265-molecule corpus (v0.17.0 State 3) — not yet re-validated after v0.18.0's atom-typing fix, but that fix's own blast-radius accounting makes a regression unlikely. |
| Multi-conformer ensemble generation | **Missing** (as a *sound* capability) | A `conformer_ensemble()` API exists (`crates/chematic-py`, wrapping `generate_conformer_ensemble_with_config`) but sits on the older, less rigorous `etkdg.rs` path, not `embed_pipeline_v2`. **Confirmed live defect, corrected 2026-08-23 during A1's own re-verification**: the originally-cited decane repro (fixed 3-iteration constraint satisfaction; ring-unaware torsion rotation) no longer reproduces on current `main` (20/20 fresh attempts stayed sound ~1.5-1.6 Å) — likely an incidental side effect of an unrelated ring-placement fix (PR #253), not established to generalize beyond that one example. The third mechanism is independently re-confirmed still live: MMFF94 silently contributes zero energy/gradient for atom-type pairs its tables don't cover — on a halomethane case (`[C@H](F)(Cl)Br`), central-C-to-halogen distances land at 8.9-11.3 Å under `conformer_ensemble(1, 0.0, 'mmff94', 0.0)`, vs. a real ~1.4-1.9 Å. See `validation/openeye_advantage_fixtures.jsonl`'s `oe-01` row. No seed control (process-global atomic counter, not caller-reproducible); no energy ranking across attempts (RMSD-uniqueness only); minimizer failures are silent (no `Result`). `embed_pipeline_v2` (the sound engine) has never been wired to an ensemble loop — `prune_rms_threshold` exists on `EmbedParameters` but is documented as unconsumed, reserved for exactly this purpose. |
| Torsion parameter coverage | Competitive | Was the dominant MMFF94-bonded-term gap through v0.16.0 (257 missing-torsion instances, 71% of `complete_bonded_term_gated` failures) — root-caused (v0.17.0) to a bond-order-classification bug, not table gaps, and fixed: 257→0 missing-torsion instances on the 265-corpus. No regression through v0.19.0. |
| Macrocycle conformer generation | Partial | 9-rule macrocycle torsion-preference table (`etkdg_knowledge/rules_macrocycle.rs`, adapted from RDKit's own macrocycle preferences) feeds `pipeline_v2`'s scoring stage, but macrocycle/ring torsions can only be *scored*, never mechanically *rotated* (`FailClosed` by default; ring bonds have no well-defined single rigid-side split for the existing bridge-bond rotation mechanism). `chematic-ff`'s force fields (MMFF94/UFF/DREIDING) have zero macrocycle-specific terms — same potential regardless of ring size. |
| RMSD/TFD/stereo-retention/runtime reported separately | Ahead of chematic's own current public claims | Already true in `validation/results/mmff94_bci_gap_227_phase2_report.md` (RMSD mean 1.685 Å, TFD mean 0.2228, stereo satisfied/violated counted separately, per-stage runtime) — but `docs/benchmark.md`/`docs/rdkit-migration.md`/even the v0.19.0 CHANGELOG still say "no RMSD/TFD comparison exists yet." This is a documentation gap, not a capability gap — flagged for a small follow-up fix outside this RFC's own scope. |
| Best-of-N conformer selection | **Missing** | The project's own benchmark script gives RDKit a `BEST_OF_N=10` UFF-optimized-then-lowest-energy arm; chematic has zero equivalent arm — every chematic measurement in that benchmark is single-embed. |

## Shape / FastROCS (Gaussian volume, shape/color overlay, query, search)

| Capability | Status | Evidence |
|---|---|---|
| Gaussian-volume shape/color Tanimoto overlay | **Missing** | Confirmed zero hits workspace-wide for "gaussian volume"/"shape tanimoto"/"ROCS"/"volume overlay". Chematic's own `ROADMAP.md` (lines 1038-1044, B-tier, explicitly gated "don't start until S/A land") already names the exact unbuilt API surface: `shape_tanimoto()`, `color_tanimoto()`, `align_shape()`, `screen_shapes()`. |
| Lighter-weight 3D shape/similarity methods | Partial | `usr.rs` (Ultrafast Shape Recognition — 12 statistical moments over 4 reference points, Soergel-distance similarity, `shape_screen()` batch search, 7 tests incl. a fixed order-dependence bug), `spectrophores.rs` (Silicos-it-style 48-element property-encoded fingerprint), `pharmacophore_fp_3d.rs`, `o3a.rs` (Open3DAlign-style MMFF-type correspondence search + Kabsch rigid alignment). All real and tested, none competitive with Gaussian-volume overlap in accuracy or the rigid-body shape-optimization step ROCS performs. |
| 2D fingerprint similarity/search infrastructure (reusable substrate) | Ahead | `crates/chematic-fp/`: ECFP/FCFP, MACCS, atom-pair, topological-torsion, layered, Avalon-style, pharmacophore-2D, MAP4, MHFP (true Lowe & Sayle circular-SMILES MinHash) — 20+ algorithms. Real Tanimoto (`bulk.rs`: `tanimoto_slice`/`tanimoto_matrix`/`top_k_similar`) and a MinHash LSH index (`lsh.rs`). Note for future Shape work: the LSH index itself is single-threaded in the reusable Rust core; Rayon parallelism currently lives only in the Python-binding layer (`chematic-py/src/bulk.rs`), not the core crate — a future 3D shape index would need its own parallelism story, not inherit this one for free. |
| Substructure/similarity query serialization for shape search | **Missing** | No shape-specific query format exists (expected — depends on the shape engine above existing first). |

## OEFF / Szybki (force field, minimization, energy ranking)

| Capability | Status | Evidence |
|---|---|---|
| MMFF94 force field | Competitive | Full bond+angle+stretch-bend+OOP+torsion+VdW+electrostatic (`crates/chematic-ff/src/mmff94*.rs`), 372 tests crate-wide. Bonded-term coverage on the 265-corpus: `complete_bonded_term_gate` 249/265 (torsion coverage fixed v0.17.0; residual is 1 bond + 46 angle type-only gaps). |
| MMFF94 minimization (steepest descent + quasi-Newton) | Competitive | `minimize_mmff94_full` (steepest descent) and `minimize_mmff94_lbfgs` (real L-BFGS, history-5, Armijo line search, "typically 2-5x fewer iterations" per its own doc, falls back to steepest descent on line-search failure) both exist, both tested (43+59 tests). Only finite-difference gradients exist (no analytic gradient anywhere in `chematic-ff` — a real, disclosed cost: ~6n energy evaluations per gradient). |
| UFF force field | Partial | Bond+angle+VdW+electrostatic only — **no torsion or out-of-plane term at all**. Confirmed permanent, non-budget-fixable defect: fused polycyclic aromatics (anthracene) converge (`converged: true`) at a structurally invalid 3.39 Å bond, stable through 30,000+ iterations, because the potential itself has no term to resist this — not a step-count/timeout problem. `UffMinimizeResult.sound: bool` (v0.19.0) discloses this rather than hiding it; the root-cause fix (add real torsion/OOP terms) is not started. |
| Robustness / catastrophic-failure handling | Competitive (disclosed, not silent) | A documented rescue path (`run_uff_bridge`) improved a 58-molecule robustness corpus from 41/58 sound (17 blown-up) to 53/58 sound (5 blown-up); the remaining 5 are declared-stereo molecules where the rescue geometry doesn't reliably preserve chirality — now gated by a `verify_stereo(...).is_fully_satisfied()` check before accepting a rescue, falling through to a typed failure otherwise rather than silently returning a wrong-stereo structure. |
| Energy-ranked conformer selection | **Missing** | Depends on the Omega-equivalent ensemble gap above — there is no multi-conformer set to rank by energy in the first place on the sound (`embed_pipeline_v2`) path. |
| Performance vs. reference implementations | Partial | Documented, not yet addressed: UFF/DREIDING/MMFF minimizer baseline ~2.2x slower end-to-end than RDKit, up to ~11x slower on FF-only arms (`ROADMAP.md` item 4, not started). |

## Docking / POSIT (receptor, pose prediction, scoring)

| Capability | Status | Evidence |
|---|---|---|
| Receptor preparation | **Missing** | See Spruce row below — no protein-specific structure handling exists at all, a prerequisite for any docking work. |
| Pose generation/prediction | **Missing** | Confirmed zero hits for "dock"/"pose"/"receptor" as functional code anywhere in the workspace (the only hits are doc-comments on `chematic-mol::pdbqt`, which *writes* ligand PDBQT files for use with **external** docking software — chematic implements no docking engine of its own). |
| Pose scoring | **Missing** | No scoring function exists (depends on pose generation, receptor prep, and — per this RFC's Track A design — sound conformer/shape/microstate foundations first). |
| **Not planned** this round or the next (A5, explicitly gated) | — | Per this RFC's Track A design: docking is deliberately sequenced last, gated on A1-A4's acceptance criteria. Building a docking engine before conformer/shape/microstate quality is solid would just relabel today's 3D-quality gap as a docking-quality gap. |

## OEDepict (2D depiction)

| Capability | Status | Evidence |
|---|---|---|
| 2D layout algorithm | Competitive | `crates/chematic-depict/src/layout.rs`: rule-based (SSSR via Balducci-Pearlman → ring-as-polygon placement, reflecting fused rings over shared edges → DFS zigzag chains → fragment spacing → bond-crossing detection), no physics simulation. Two real layout bugs (long-chain rotation drift/self-wrap; exocyclic-bond-angle misplacement at ring junctions) found and fixed in v0.19.0 (issue #347) — explicitly disclosed as a breaking coordinate-output change (no fixture pins the old buggy output). No further open layout issue after that fix. |
| Rendering formats | Competitive | SVG (with highlighting), PDF (via `svg2pdf`), EPS (pure-Rust), PNG (tiny-skia, optional feature, auto-disabled on WASM). |
| Multi-molecule grids, reaction depiction, SimilarityMap | Competitive | `grid.rs` (`depict_svg_grid`/`_with_opts`/`_highlighted`), `reaction_svg.rs` (`depict_reaction_svg`), `similarity_map.rs` (atom-colored SVG by scalar weight, e.g. LogP/TPSA/fingerprint-contribution visualization). |
| Overall assessment | **Competitive, not a priority gap.** | This module does not need Track A/B investment; it is already at parity for the common depiction use cases OpenEye's OEDepict covers. |

## Spruce (protein/structure preparation)

| Capability | Status | Evidence |
|---|---|---|
| Protein/ligand/water separation | **Missing** | `crates/chematic-3d/src/pdb.rs`'s `pdb_to_molecule` treats every ATOM and HETATM record identically, building one flat `Molecule` via pure interatomic-distance bond inference. `chain_id`/`res_name` are captured on the intermediate `PdbAtom` struct but never used to split chains or separate ligand/protein/water at the molecule-building step. `mmcif.rs` preserves the equivalent fields for round-trip format fidelity, but no function anywhere (`is_hetatm`/`separate_chains`/`extract_ligand`/`protein_only` — zero hits outside test assertions) acts on them. |
| Missing-atom / missing-residue detection | **Missing** | No such logic exists anywhere in the workspace. |
| Protein-residue protonation | **Missing** | `pka.rs` is small-molecule SMARTS-pattern-only, with no residue-context awareness. |
| Binding-site detection | **Missing** | No such logic exists. |
| Even ligand-extraction from PDB/mmCIF (a prerequisite, far short of full Spruce-style prep) | **Missing, explicitly future work** | `ROADMAP.md` line 1055 lists "PDB + PDBx/mmCIF ligand-subset I/O" as unstarted C-tier work — i.e. chematic's own roadmap already acknowledges it hasn't even reached the prerequisite step. |

---

## Summary

| Module | Overall | Priority for Track A |
|---|---|---|
| OEChem | Competitive/Partial (format breadth intentionally deprioritized) | Not a Track A target — covered by the existing 100-point ladder |
| Quacpac | Partial→Missing (microstate enumeration) | **A4** |
| Omega | Partial (sound single-conformer; broken/unsound ensemble) | **A1, A2** |
| Shape/FastROCS | Missing | **A3**, gated behind A1/A2 |
| OEFF/Szybki | Competitive (MMFF94), Partial (UFF) | Supporting work for A2; UFF torsion/OOP fix is a separate, smaller item not in this RFC's Track A/B lists |
| Docking/POSIT | Missing | **A5**, explicitly gated last, not planned this round or next |
| OEDepict | Competitive | No investment needed |
| Spruce | Missing | Not in Track A/B — noted as a real gap but out of scope for this round's recommendation |
