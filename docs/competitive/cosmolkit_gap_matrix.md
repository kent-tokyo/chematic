# COSMolKit Gap Matrix

Full directive: `docs/rfcs/cosmolkit_advantage_rfc.md`. Method: every chematic
row is grounded in this repo's actual source, tests, CI config, and measured
results as read on 2026-08-26 (see snapshot below). COSMolKit rows are drawn
from their own public README, `VALIDATION.md`, `Cargo.toml`, and CHANGELOG —
this project has not built, installed, or run COSMolKit itself, so COSMolKit's
own quantitative claims (validation counts, parity percentages) are reported
as **self-reported**, not independently re-verified, exactly as this
directive's own discipline requires (§15: "測定していない項目を「対応」と書
いてはいけません" — never write "supported" for something unmeasured).

```
COSMolKit comparison snapshot:
version: 0.2.13 (crates.io + PyPI, confirmed matching)
commit: 68ac32968abed60ed5f21e73667c5d91406e9859 (main, 2026-08-26T00:34:56Z)
checked_at: 2026-08-26
source: github.com/cosmol-studio/COSMolKit README.md, VALIDATION.md, Cargo.toml
        (fetched via `gh api`/`gh repo view`; not cloned/built/run locally)
```

Classification: **Ahead** / **Competitive** / **Partial** / **Missing** /
**Not planned**. A row can be "Partial" on both sides at once (e.g. neither
project has shipped ML/array export yet) — this is noted explicitly rather
than forced into a winner.

## Headline correction to this track's original premise

The directive that started this track assumed WASM/browser deployment was
COSMolKit's clearly unfinished area. **That premise is now out of date**:
COSMolKit's current README documents "Browser-native deployment... through
WebAssembly" and a live public site, COSMolKit Web Tools
(tools.cosmol.org/tools); their `Cargo.lock` confirms a real `wasm-bindgen`
dependency (not just marketing text). chematic remains ahead on *measured
bundle size* and *exposed API surface breadth* (see Deployment section), but
"has a live browser tool at all" is now a **parity** point between the two
projects, not a unique chematic advantage. This correction should propagate
to any future messaging about this track.

## Self-reported COSMolKit validation scale (not independently verified)

COSMolKit's `VALIDATION.md` claims a complete ChEMBL 37 profile (2,897,819
source records) checked against pinned RDKit `2026.03.1` across chemistry,
descriptors, fingerprints, 3D/UFF/MMFF, InChI, and I/O -- **2,960,559,232
matching checks, self-reported zero blocking mismatches**, plus a
5,000-molecule "maintained strict" corpus and a 152-molecule daily corpus.
This is a self-reported claim from their own repository, not something this
project has reproduced or audited. It is nonetheless evidence of a
seriously-resourced validation methodology (pinned reference version,
git-tracked phase definitions, documented tolerances: `1e-8` for matrix
entries, `1e-6` for coordinates/energies/gradients) that materially outweighs
chematic's own largest single-round validation efforts (5,000-molecule dev
corpora, a 4,999-molecule NCI holdout, a 265-molecule 3D benchmark corpus) by
roughly three orders of magnitude in corpus size. This is the single most
important fact for Axis 1 (化学的正確性・identity, the highest-weighted axis)
and should not be understated in any future round's scoring.

## Core Chemistry

| Capability | chematic status | COSMolKit status | Evidence | Known limitation |
|---|---|---|---|---|
| Molecule/atom/bond graph model | Competitive | value-style, explicit `_`-suffixed in-place API | `chematic-core/src/molecule.rs:46` -- `Molecule` exposes `with_atom_added`/`with_atom_charge`/`with_bond_order`/etc. (value semantics, `&self -> Molecule`); separate `MoleculeBuilder` (mutation-based construction) and Python-only `RWMol` (explicit in-place edits, RDKit-compat) | Split across 3 types rather than COSMolKit's single type + `_`-suffix convention; same "value ops never mutate input" guarantee holds, different ergonomics |
| SMILES parsing | Competitive | shipped | `chematic-smiles`; #395/#399 fixed this session (see issue history) | #402/#403/#407 open, small, traced residuals |
| Canonical SMILES writing | Competitive | shipped, RDKit-style writer options | same | same |
| Random SMILES generation | Ahead | not documented | `chematic-smiles/src/random_smiles.rs:14,30` -- real implementation, not a stub | -- |
| CXSMILES | Ahead (tentative) | not documented in README | `chematic-smiles/src/cx.rs`, `chematic-smarts/src/cx.rs`, plus Python/WASM bindings | Not deeply stress-tested against a CXSMILES-specific corpus this round |
| SMARTS parsing + recursive SMARTS | Competitive | shipped, Python SMARTS parse metadata | `chematic-smarts/src/parser.rs:83-153`, `MAX_RECURSION_DEPTH`, tested | No query-serialization API found (zero grep hits) -- COSMolKit documents "Python SMARTS parse metadata" as a shipped feature; likely a gap |
| Substructure matching, uniquify, chirality-aware matching | Unconfirmed | shipped | Matching itself confirmed real; chirality-awareness and a dedicated "uniquify" pass not directly verified this round | Needs a follow-up audit before claiming parity either way |
| Ring perception / SSSR | Competitive | not explicitly claimed as symmetrized | `chematic-perception/src/sssr.rs:88`, explicit symmetry-equivalent-ring handling (line 843) | -- |
| Valence/aromaticity/Kekulization/sanitization | Partial | shipped as a single named "Sanitization" feature | Operations exist (standardize.rs, kekulize, aromaticity perception) but scattered -- zero hits for any dedicated `sanitize` module/entrypoint | API-ergonomics gap: no single named "sanitize" pipeline entrypoint the way COSMolKit documents one |
| Tetrahedral + double-bond stereo | Competitive | shipped (inspection only, no square-planar/generalized-geometry claim) | Core stereo model | -- |
| Non-tetrahedral / square-planar stereo, generalized rotation-group geometry | Ahead | not claimed | `chematic-core/src/stereo_geometry.rs` (PR #313/#326: square-planar `@SP1/@SP2/@SP3`, A4/D4 rotation groups) | Genuine, distinct chematic differentiator -- COSMolKit's README only claims tetrahedral-shape stereo inspection |
| Isotope, atom maps, substance groups (stereo groups) | Ahead | not itemized | `Atom.atom_map: Option<u16>` (`atom.rs:108`), `StereoGroup` (`stereo_group.rs:39`) both first-class | -- |

## Notation and Search

(See Core Chemistry table above -- SMILES/SMARTS/CXSMILES rows cover this
section; COSMolKit's README doesn't split these into a separate category.)

## Fingerprints and Descriptors

| Capability | chematic status | COSMolKit status | Evidence | Known limitation |
|---|---|---|---|---|
| Morgan/ECFP | Competitive | shipped, source-backed exact-parity claimed | `chematic-fp/src/ecfp.rs`, `morgan_environment.rs`, `rdkit_morgan_ecfp4.rs` -- two modes: default FNV-hashed (invariant-parity only) + `RdkitMorgan` (claimed hash/environment/dedup parity). Session memory cites 5,046/5,046 environment parity | Not re-verified this pass -- treat prior figure as historical, not reconfirmed |
| MACCS | Competitive | shipped | `maccs.rs`, 166 keys per official RDKit `MACCSkeys.py` definitions | -- |
| RDKit topological / RDKFingerprint | Partial | shipped, exact-parity claimed at ChEMBL-37 scale | `path.rs` (RDKit-compatible defaults) + separate `topo_path.rs` (independent, no parity claim) -- two overlapping implementations, no bit-exact test found | Real gap on the bit-exact-parity axis specifically |
| Avalon | Partial (self-flagged) | shipped, claimed bit-exact at ChEMBL-37 scale | `avalon.rs`'s own doc comment: "not a bit-exact reimplementation... the bar is deterministic, isomorphism-invariant, similarity-preserving hashing" | Honest, material gap vs. COSMolKit's claimed exact parity |
| Pattern fingerprint (tautomer-aware) | Partial | shipped, tautomer-aware, claimed exact parity | `pattern.rs`: from-scratch atom-centric hash, no RDKit-parity claim, no tautomer-awareness found | Same shape as Avalon |
| AtomPair (2D/3D) / Topological Torsion | Partial | shipped, 2D+3D distances, exact parity claimed | `atom_pair.rs`: 2D-topological-distance only found; no 3D-distance variant located this pass | 3D AtomPair variant may exist in `chematic-3d`, not confirmed |
| Layered fingerprint | Partial (both sides caveat) | experimental / upstream-classified on COSMolKit's own side | `layered.rs`: 7-layer bit-packing, no RDKit-parity claim | Both sides have caveats here, different ones |
| MHFP / MAP4 / pharmacophore (2D+3D) / ERG | Ahead | not documented | `mhfp.rs` (literature-faithful Lowe/Sayle), `map4.rs` (Minervini 2020), `pharmacophore_fp.rs` + `chematic-3d`'s 3D variant, `erg.rs` (Sheridan/Rarey-Dixon + Ertl 2017 FG detection) | Genuine breadth advantage, not confirmed absent on COSMolKit's side (only confirmed absent from their public README) |
| Reaction fingerprint | Partial (self-flagged) | not documented | `reaction_fp.rs`'s own doc comment: placeholder, "OR of reactant/product ECFP4... not the actual transformation" | Known, self-admitted gap vs. RDKit's real algorithm |
| MW/formula/HBD/HBA/TPSA/Crippen logP+MR/Fsp3/QED/rotatable bonds | Competitive-to-Ahead | shipped, documented parameter space | `chematic-chem/src/descriptors.rs`, `qed.rs` -- session memory: TPSA +/-0.1, HBA/HBD 99.98% on 4999 molecules | Not re-verified this pass |
| Topological descriptors (Wiener, kappa1-3, chi0-4/chi0v-4v, Bertz, LabuteASA), VSA bins, xlogp3, drug-likeness filter suite (Lipinski/Veber/Egan/REOS/Ghose/Rule-of-3/Pfizer), activity cliffs | Ahead | not documented | `topo_descriptors.rs`, `vsa.rs`, `xlogp3.rs`, `drug_score.rs`, `activity_cliff.rs` | Genuine breadth advantage |

## 3D and Force Fields

| Capability | chematic status | COSMolKit status | Evidence | Known limitation |
|---|---|---|---|---|
| ETKDG-equivalent embedding | Partial | shipped, DG/KDG/ETDG/ETKDG presets, `etkdg_v3()` | Production default: `embed_pipeline_v2` (`chematic-3d/src/pipeline_v2.rs`), 12-stage typed pipeline. A better engine, `generate_coords_connectivity_ordered` (PR #391, fixes ring/chain-placement defects, measured 10/33 to 33/33 on a targeted corpus), exists but is not wired into any production caller | Real routing gap: the better engine is available but not the default |
| Multi-conformer / ensemble | Competitive | shipped, `with_3d_conformers(n)` | `embed_ensemble_v2` (PR #371) -- energy-ranked, automorphism-aware pruning; bound to Python | Not yet bound to WASM (verified: zero references in `chematic-wasm/src/*.rs`) |
| Deterministic seeding | Competitive | shipped, `random_seed` | `PipelineV2Config.random_seed: u64`, tested | -- |
| Typed failure/attempt provenance | Competitive | shipped, `track_failures`/`params.failures` | 12-stage typed failure enum, per-attempt provenance | -- |
| UFF | Partial | shipped, claimed ChEMBL-37-scale parity | `chematic-ff/src/uff.rs`; v0.20.0 partially fixed issue #210 (chirality-ignoring rescue path) -- CHANGELOG explicitly states naproxen_S/ibuprofen_S/testosterone/cholesterol remain unfixed, #210 stays open | `CatastrophicBondBlowup` failure mode still live on at least one corpus molecule |
| MMFF94 | Competitive | shipped, claimed ChEMBL-37-scale parity, 524,288 records | `chematic-ff/src/mmff94*.rs` (6 files) + separate DREIDING FF (not claimed by COSMolKit at all); measured `pipeline_v2_mmff94_strict_complete_bonded_term_gated` = 241/265 on chematic's own 265-molecule corpus | Chematic's own validation corpus is ~2000x smaller than COSMolKit's self-reported ChEMBL-37 scale |
| Stereo-safe 3D | Ahead | not claimed | `PipelineV2Config::stereo_safe()`: 144/145 (99.3%) correct_and_ok, 0 silently-wrong, 0 loud-failure-stereo on a 29-molecule x 5-seed corpus; bound in both Python and WASM | Precisely scoped, real, shipped v0.20.0 headline feature |
| Chiral-tag assignment FROM a 3D conformer | Competitive | shipped, 77 fixed full-state oracle records, claimed exact pinned-RDKit parity | `assign_stereo_from_3d` (`chematic-3d/src/stereo3d.rs:228`) -- direct equivalent exists | No equivalently-sized fixed-oracle validation set confirmed for this specific function |
| DG bounds matrix as standalone public API | Unconfirmed | shipped | `bounds_conformance()` (`distance_geometry_v2.rs:656`) checks coordinates against bounds; a standalone bounds-matrix-returning function not confirmed to exist as public API this pass | Needs follow-up grep before claiming parity or gap |
| RMSD / symmetry-aware RMSD | Competitive | shipped | `rmsd_no_align` (`align.rs:34`), `rmsd_symmetric` (`conformer.rs:343`, genuine symmetry-aware) | -- |
| TFD (torsion fingerprint deviation) | Missing | shipped (implied via general FF/conformer parity claims) | Benchmark docs' "median TFD 0.344" figure is computed by calling RDKit's own `TorsionFingerprints.GetTFDBetweenMolecules` from a Python benchmarking script -- chematic has no native TFD implementation | Precise, real gap -- not a native chematic capability at all |
| Shape similarity / shape screening | Unconfirmed | not directly itemized | `shape_screen`/`o3a_align` are real Python-exposed functions (`chematic-py/src/similarity.rs:46,460`) but their algorithmic backing (Gaussian-volume overlap vs. something simpler) was not traced this round | Do not claim parity/superiority without checking the actual algorithm behind these names |

## File I/O

| Capability | chematic status | COSMolKit status | Evidence | Known limitation |
|---|---|---|---|---|
| MOL/SDF read+write (V2000/V3000) | Ahead/Competitive | shipped | `chematic-mol/src/mol2000.rs`, `mol3000.rs`, `sdf.rs` -- `SdfRecord` carries rich per-record diagnostics (stereo/E-Z/stereo3d/geometry-rank) beyond a bare read | Exact V3000 branch coverage vs. RDKit not independently re-verified this round |
| MOL2 read+write | Ahead | read only | `chematic-mol/src/mol2_tripos.rs`: both `parse_mol2` and `write_mol2` exist; COSMolKit documents read only | -- |
| XYZ read+write, extended XYZ | Ahead | read only (block read) | `chematic-mol/src/xyz.rs`: parse+write, plus extended XYZ (materials-focused) -- COSMolKit doesn't mention extended XYZ at all | -- |
| SDF indexed/streaming dataset | Partial | shipped, `SdfDataset`: byte-range index, random access, `n_jobs` batches | `SdfFileIter`/`SdfBatchIter` (`chematic-py/src/io.rs`) genuinely stream (verified: sequential `BufReader`-backed, not load-then-iterate) but are forward-only -- no byte-range index, no random access (`dataset[100]`-equivalent), no built-in parallel-worker parameter | Real, working streaming but strictly narrower than COSMolKit's indexed+random-access+parallel design |
| PDB/mmCIF -> structured protein model | Missing | shipped, `Protein`/`BioStructure`: full model/chain/residue/entity/hetero/water/ligand/nucleic-acid graph | `chematic-3d/src/pdb.rs`'s `pdb_to_molecule()` builds a single flat `Molecule` via distance-based bond inference; chain/residue/model identity is read then discarded. `chematic-mol/src/mmcif.rs`/`pqr.rs` are the same flat shape. Repo-wide grep for `struct Chain`/`struct Residue`/`struct Protein`/`struct BioStructure`: zero matches | Complete structural gap, not a maturity gap -- no such layer exists at all |
| Binary serialization round-trip for Molecule | Missing | shipped, claimed ChEMBL-37-scale parity | Repo-wide grep for Python object-serialization hooks in `chematic-py/src/`: zero matches | -- |

## Batch and ML Readiness

| Capability | chematic status | COSMolKit status | Evidence | Known limitation |
|---|---|---|---|---|
| Ordered batch with per-record errors | Partial | shipped, `MoleculeBatch`, composable chained transforms | `screen_smiles`/`screen_smiles_with_options` (`chematic-chem/src/workflow.rs`): genuinely ordered (`input_index` preserved), genuine per-record error retention | Single-purpose (always produces a full `MoleculeReport`, not a composable chained-transform pipeline); sequential only, no `n_jobs`; no progress-bar hook; no `errors="raise"/"skip"/"keep"` mode selection |
| ML/array export (atom/bond feature matrices, graph tensors) | Missing (both sides) | planned, not shipped (COSMolKit's own README marks this not yet public) | Repo-wide grep for feature-matrix/graph-tensor helpers: zero matches | Neither project has shipped this -- a race not yet run, not a loss |

## Protein and Structural Biology

| Capability | chematic status | COSMolKit status | Evidence | Known limitation |
|---|---|---|---|---|
| High-level structure model (chain/residue/atom hierarchy) | Missing | shipped, `BioStructure` full structure + `Protein` amino-acid projection | See "PDB/mmCIF -> structured protein model" row above | Complete gap |
| Selection utilities (chain/residue/neighborhood queries) | Missing | planned, not shipped either (COSMolKit's own README marks this not yet public) | Follows from the above -- nothing to select over | Neither project has this yet, but COSMolKit at least has the base hierarchy to eventually build selection on top of |
| Ligand extraction from a structure | Missing | planned (COSMolKit's own README marks ligand/mixed-structure ergonomic APIs not yet public) | `pdb_to_molecule()` already discards HETATM/ATOM distinction | Both sides incomplete, chematic further behind (no hierarchy at all vs. COSMolKit having the hierarchy but not yet the ergonomic API on top) |

## Deployment (Rust/WASM/Browser)

| Capability | chematic status | COSMolKit status | Evidence | Known limitation |
|---|---|---|---|---|
| WASM crate existence + API surface breadth | Ahead | shipped, surface breadth not itemized in their README | `chematic-wasm` wraps 11 of ~19 Rust crates (core/smiles/perception/chem/fp/depict/rxn/smarts/3d/mol/inchi/iupac/ff) | `chematic-crystal` and `chematic-cip` are NOT exposed to WASM |
| Live public browser tool | Parity (corrected from "chematic-only") | shipped, tools.cosmol.org, confirmed `wasm-bindgen` dependency | `kent-tokyo.github.io/chematic/` confirmed live (HTTP 200), auto-deployed via `.github/workflows/pages.yml` on every relevant push to `main` | Both projects have a real, live, deployed site -- no longer a unique chematic differentiator |
| npm package | Partial (naming issue) | shipped, `pip install cosmolkit` | Actual current publish path is `@kent-tokyo/chematic` (confirmed live at 0.20.0); the unscoped `chematic-wasm` name in `package.json` is stale/abandoned at 0.1.33 on the registry | Real, fixable doc/consistency gap, not functional |
| Bundle size (measured) | Ahead | not published | `demo/pkg/chematic_wasm_bg.wasm` ~= 2.94 MB; README documents "~2.3x smaller than RDKit.js's `RDKit_minimal.wasm` (6.91 MB)"; no equivalent COSMolKit figure found | -- |
| Web Worker support | Missing | not documented either | Only a textual, likely-aspirational mention in one demo's README; no actual Worker message-passing code found | Do not claim as shipped |
| JS/TS test coverage, dual Node+browser CI target | Ahead | not itemized | 9 `.test.mjs` files including a dedicated `pipeline_v2_web_target.test.mjs`; CI builds and tests both `--target nodejs` and `--target web` | -- |

## Summary read (qualified, not a verdict -- Phase A snapshot only)

Per this directive's own section 3 win conditions, **no "chematic exceeds
COSMolKit" claim is warranted from this round** -- Phase A is a
capability-shape audit, not a run head-to-head benchmark (that is Phase C's
job, not yet started). Qualified reads only:

- **Axis 1 (chemical correctness/identity, weight 30)**: COSMolKit's
  self-reported ChEMBL-37-scale validation is the single largest fact in this
  whole matrix and should not be dismissed; chematic's own validation, while
  real and methodologically sound (independent oracles, holdouts, permutation
  invariance), operates at roughly 1/1000th the corpus scale. chematic has
  distinct strengths within this axis (non-tetrahedral stereo, explainable
  standardization with typed provenance) COSMolKit doesn't claim at all.
- **Axis 2 (3D/FF/conformer, weight 20)**: roughly competitive; each side has
  narrow, real advantages (chematic: stereo-safe 3D; COSMolKit: claimed
  validation scale + native TFD).
- **Axis 3 (Python/batch/ML, weight 20)**: COSMolKit ahead -- real
  `MoleculeBatch` + `SdfDataset` vs. chematic's narrower, single-purpose
  equivalents. ML/array export is unshipped on both sides.
- **Axis 4 (Rust/WASM/deployment, weight 15)**: chematic ahead on breadth and
  measured bundle size; "has WASM/a live site" itself is now parity.
- **Axis 5 (feature breadth/docs/operational quality, weight 15)**: mixed --
  chematic ahead on materials/crystal/reaction chemistry (COSMolKit doesn't
  attempt this at all) and fingerprint/descriptor breadth; COSMolKit ahead on
  protein/structural-biology (chematic has zero hierarchy layer) and
  validation-documentation rigor.

## Not investigated this round

Query serialization for SMARTS; uniquify/chirality-aware substructure
matching depth; DG bounds-matrix standalone public API; shape-similarity
algorithmic backing; exact V3000 branch coverage; re-verification of several
session-memory-cited percentages (Morgan 5,046/5,046, TPSA +/-0.1, HBA/HBD
99.98%) against current `main`.
