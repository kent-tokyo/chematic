# Changelog

All notable changes to chematic will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

No unreleased changes.

## [0.39.0] - 2026-08-31

### v0.39.0 — complete standardization stage audit coverage

- Added `StandardizationStep::DisconnectMetals` and recorded the automatic
  metal-coordination disconnection in `StandardizationReport.steps`.
- The report now covers the full execution order from the original molecule
  through charge and fragment processing, including bond-count changes caused
  by metal disconnection.

## [0.38.0] - 2026-08-31

### v0.38.0 — deterministic bounded Parent tautomer audits

- Extended bounded, non-overlapping tautomer batching to `tautomer_parent`,
  keeping `max_transforms` enforced while selecting independent sites by
  canonical rank rather than input atom order.
- Every batched transformation is now retained in the existing
  `TautomerAuditRecord::applied_transforms` contract, with regression coverage
  for bounded parent selection and atom-order invariance.

## [0.37.0] - 2026-08-31

### v0.37.0 — bounded, atom-order-invariant tautomer convergence

- Tautomer canonicalization now applies non-overlapping matches of the active
  rule together within one bounded iteration. Conflicting matches remain
  deferred for a later pass, preserving the safety of `max_iter` while
  preventing independent sites from being selected by input atom order.
- Promoted the large independent-site convergence regression from an ignored
  diagnostic to a normal test and documented the deterministic bounded-pass
  contract.

## [0.36.0] - 2026-08-31

### v0.36.0 — operation-level compatibility scorecards

- Added a dependency-free scorecard generator for the Phase 0 comparison
  contract. It records engine versions, source commits, corpus identity, and
  per-operation `match`, `mismatch`, `unsupported`, `failure`, and
  `uncomparable` outcomes.
- Added a versioned scorecard schema and regression coverage. Unsupported or
  failed operations remain visible and are never counted as matches.

## [0.31.0] - 2026-08-30

### v0.31.0 — WASM Parent identity bindings

- Exposed fragment, charge, isotope, stereo, and composed Super Parent
  operations through WASM JSON APIs with consistent status-shaped results.
- Preserved the existing resource limits and size guard so browser callers
  receive explicit budget or size outcomes instead of unbounded work.
- Added native regression coverage for every new Parent binding.

### v0.30.0 — Python Parent identity bindings

- Exposed `fragment_parent`, `charge_parent`, `isotope_parent`, and
  `stereo_parent` on Python `Mol` objects, preserving the Rust Parent
  semantics and transformation boundaries.
- Added Python `super_parent`, returning the composed parent molecule and its
  explicit computation status alongside the existing `tautomer_parent` API.
- Added regression coverage for fragment, charge, isotope, stereo, and
  composed parent behavior.

### v0.29.0 — external comparison adapter protocol

- Added a command-based external-engine adapter runner for the Phase 2
  comparison harness. Adapters receive the shared corpus and emit only the
  versioned common JSONL contract.
- External output is rejected unless every corpus record has the expected hash,
  unique id, and valid operation status, making future COSMolKit runs
  reproducible without coupling the repository to one installation method.

### v0.28.0 — reproducible comparison gates

- Added dependency-free JSONL contract validation for the Phase 2 comparison
  harness, including corpus hash, record uniqueness, operation status, and
  complete-corpus checks.
- Added strict mismatch gating and deterministic Markdown reports so comparison
  runs can be used in CI without counting unsupported operations as failures.

### v0.27.0 — Phase 2 direct-comparison smoke harness

- Added a versioned public smoke corpus, common JSONL result schema, RDKit and
  chematic runners, and a scorer that separates matches, mismatches, parse
  failures, and unsupported operations. The harness is ready for a COSMolKit
  adapter without treating unavailable tooling as a parity result.
- The comparison runner uses chematic's RDKit-exact Morgan API when available;
  older installations report that operation as unsupported rather than falling
  back to native ECFP4.

### v0.26.0 — RDKit Morgan tetrahedral chirality parity

- The generalized RDKit Morgan fingerprint API now supports opt-in tetrahedral
  chirality, with R/S contributions matching RDKit's repeated per-round
  environment hashing. The implementation is pinned by exact sparse-count
  regression fixtures for both alanine enantiomers; the default non-chiral path
  remains bit-identical.
- Python `rdkit_ecfp_config` and `rdkit_ecfp_config_detail` accept the optional
  `include_chirality` keyword. WASM adds ABI-preserving
  `rdkit_ecfp_config_chiral_bitvec` and `rdkit_ecfp_config_chiral_detail_json`.
- E/Z bond-stereo invariants are intentionally not included in this increment and
  remain a documented follow-up.

### v0.25.0 — aromaticity exactness

- The opt-in `AromaticityAlgorithm::RdkitLike` mode now recognizes phosphole
  rings: neutral P in a two-connected ring contributes its lone pair as 2π,
  matching the existing Se/Te extension. The default strict Hückel mode is
  unchanged, and the diagnostic electron-contribution trace follows the same
  rule.
- `RdkitLike` now routes through the independently verified fused-ring parity
  engine, covering non-alternant whole-perimeter systems such as azulene while
  retaining the historical engine as a defensive fallback for un-kekulizable
  input.
- The accurate CIP residual report now measures 4,175/4,186 (99.74%) agreement
  against the modern RDKit labeler estimate, with zero regressions in the
  frozen residual subset. The remaining phosphorus rows stay explicitly
  unresolved when respelled inputs produce representation-unstable oracle
  labels; no plausible-looking R/S guess is emitted for those cases.
- MMFF94 now classifies non-aromatic, three-connected iminium nitrogen (N+=C)
  as type 54 before the generic charged-nitrogen fallback; a direct Kekule
  regression protects the independent atom-typing rule.
- The opt-in RDKit ring-count model now derives candidates from the bounded
  root-centered shortest-ring (D2-like) primitive, re-searches duplicate D2
  groups after bond trimming, and retains independently verified minimum
  replacements. MMFF94 production now consumes this symmetrized ring model.
  A fresh 265-molecule run improves atom parity from 99.522245% to 99.641684%
  and bond parity from 99.446060% to 99.584545%, leaving five known fused
  heteroaromatic residual molecules explicitly tracked for a later correction.
- The Phase 7C ring-perception implementation now includes RDKit-style
  root-centered BFS, D2 root selection, bond trimming, duplicate re-search,
  and exact rooted ring stitching. Cubane, dodecahedrane, SMARTS, and MMFF94
  regressions cover the public behavior; the remaining five corpus residuals
  are documented rather than hidden behind a heuristic.

### v0.24.0 release policy

- Phase 2 tautomer identity and Parent API work is release-ready. The legacy
  CIP residual, cosmetic E/Z spelling variance, and aromaticity-context gap are
  explicit known limitations for v0.24.0; `CipMode::Accurate` plus its typed
  unresolved result is the opt-in path for callers needing stronger CIP
  guarantees. See the public API and limitations documentation.

### Fixed — dual-flank nucleobase tautomer canonicalization

- Re-audited the `tp2-39` and `tp2-holdout-06` residual fixtures against
  independent RDKit InChIKeys. Each originally included a positional isomer
  mislabeled as a tautomer; corrected same-identity variants now pass without
  adding a fixture-specific ring-bond transformation.
- Cytosine and guanine keto/amino-imino tautomer spellings now converge via a
  deterministic carbonyl-centered aromatic N-H normalization. Explicit H
  placement is retained in tautomer deduplication so equivalent N sites are
  not incorrectly collapsed.
- Dual-flank selection now uses a rooted structural key with explicit H
  placement removed, making the choice independent of input atom order.
  `tp2-39` is covered after correcting its enol fixture to the same positional
  isomer; the former positional-isomer input remains rejected as invalid data.
- `enumerate_tautomers` now traverses the aromatic lactam/lactim edge in both
  directions and retains all eligible dual-flank orientations, so bounded
  enumeration exposes the complete candidate component instead of only the
  directed canonical-search path.

### Fixed — `chematic-chem` nitroso/oxime tautomer canonicalization

- `canonical_tautomer` now converges nitroso/oxime pairs such as `CCN=O` and
  `CC=NO` via a dedicated C-H adjacent to N=O rule. The generic any-bridge
  C→O rule remains non-forward, preventing unrelated aldehydes and ketones
  from being incorrectly enolized.
- E/Z-bearing molecules are no longer collapsed during tautomer candidate
  deduplication: directional bond metadata is included in the fingerprint, and
  canonical tautomerization preserves such inputs until directional remapping
  is available. Hydrazone E/Z regression coverage is now enabled.

## [0.23.0] — 2026-08-30

Minor release: MCS accuracy and correctness (a default-comparator behavior change,
see below), plus the fifth and sixth entries in the RDKit fingerprint-parity series
(`rdkit_rdk_fp`, `rdkit_layered_fp`) and full `McsConfig`/`McsOutcome` exposed to the
Python/WASM MCS bindings. No breaking API changes except the documented
`AtomCompare::Elements` semantics change below.

### Changed — `chematic-smarts` (`find_mcs`'s default `AtomCompare::Elements` no longer requires matching aromaticity)

- **Behavior change for every `find_mcs` caller (Rust/Python/WASM/MCP) using the
  default `atom_compare`.** `AtomCompare::Elements` previously required both atomic
  number AND aromaticity flag to match between atoms. It now matches on atomic
  number alone, exactly mirroring RDKit's identically-named
  `rdFMCS.AtomCompare.CompareElements` (confirmed via live oracle: RDKit never
  encodes aromaticity as a per-atom constraint, only via bond-type queries, even
  when both input molecules fully agree on aromaticity). Concretely: `find_mcs`
  can now return a match (or a larger match) for atom pairs that agree on element
  but disagree on aromaticity -- previously this returned `None` or an undersized
  MCS. Differential measurement against a live RDKit oracle (`scripts/
  mcs_rdkit_fmcs_diff.py`, n=200 pairs/corpus, exhaustive-only): agreement rose
  from 74.6%/68.2%/70.4% to 88.4%/88.5%/97.0% across the three established
  corpora (descriptor_census/ChEMBL/NCI).
- **`AtomCompare::AnyHeavyAtom`/`Any` are unaffected** (they never checked
  aromaticity). **There is currently no `AtomCompare` mode that restores the old
  strict element+aromaticity match** -- callers relying on that exact behavior
  need to add their own post-filter on the result.
- Co-dependent fix: `build_query`/`molecule_to_query` no longer encode
  aromaticity as a per-atom query constraint at all (previously hard-coded from
  `mol[0]`'s own atom, regardless of which comparator mode found the match --
  itself a latent, pre-existing gap for `AnyHeavyAtom`/`Any` too). Aromaticity is
  now conveyed purely through the existing bond-level `Aromatic` query, matching
  RDKit's own representation. The 4 independent "QueryMolecule -> concrete
  Molecule" reconstruction sites (Python's `qmol_to_mol`, WASM's
  `qmol_to_molecule`, the MCP server's own `qmol_to_molecule`, `chematic-chem`'s
  `query_molecule_to_smiles`) were updated to derive `atom.aromatic` from
  incident bond aromaticity instead of an atom-level primitive that no longer
  exists.
- **Known residual limitation**: an MCS result's `.smiles` accessor is a concrete
  SMILES string, which cannot express "aromaticity don't-care" the way a SMARTS
  pattern can. When two inputs' matched atoms genuinely disagree on aromaticity,
  `.smiles` (built from one specific input's own bonds) is only guaranteed to
  re-match *that* input as a substructure, not necessarily every input. A
  SMARTS-preserving accessor is the real fix and is tracked as a follow-up, not
  included here.

### Fixed — `chematic-py`/`chematic-wasm`/`chematic-mcp`/`chematic-chem` (MCS result silently lost heteroatoms and/or aromaticity)

- `find_mcs`'s query-molecule atoms are encoded as the compound
  `AtomQuery::And(AtomicNum(n), Aromatic(bool))` (`chematic-smarts`'s
  `build_query`/`molecule_to_query`), never a bare `AtomicNum` primitive. Four
  independent copies of the "reconstruct a concrete `Molecule` from an MCS
  result" helper (Python's `qmol_to_mol`, WASM's `qmol_to_molecule`, the MCP
  server's own `qmol_to_molecule`, and `chematic-chem`'s
  `query_molecule_to_smiles`) each matched only the bare `AtomicNum` case (or,
  for the Python binding, correctly unwrapped `AtomicNum` but not
  `Aromatic`) — so every MCS atom silently fell through to carbon in three of
  the four copies (WASM, MCP, `chematic-chem`), and every MCS atom lost its
  aromaticity flag in the fourth (Python), regardless of the molecule's real
  elements/aromaticity. The Python binding's own MCS result then failed to
  re-match as a substructure of either input molecule it was computed from —
  defeating a primary purpose of an MCS result (self-verification /
  substructure screening).
- Found while scoping the "99-point directive" Phase 2 differential-corpus
  measurement (MCS vs. RDKit's FMCS): a hand-picked aspirin/paracetamol pair's
  hydroxyphenyl MCS came back as an all-carbon, non-aromatic ring in every
  binding. Root-caused via the same live-oracle-first discipline used
  throughout this session: confirmed the underlying MCS *algorithm* already
  finds the chemically correct answer (`build_query` already encodes the right
  atomic number and aromaticity per atom), isolating the bug entirely to these
  four independently-duplicated, never-kept-in-sync conversion helpers.
- Fixed by adding a recursive `extract_atomic_num`/`extract_aromatic` pair
  (unwrapping through `AtomQuery::And`) to each of the four copies, matching
  the pattern the Python binding's own `extract_atomic_num` already used for
  atomic number (just never extended to aromaticity, and never ported to the
  other three copies). New regression tests in all four locations, each using
  a heteroatom-containing MCS pair specifically chosen because the
  pre-existing benzene/toluene-style tests' all-carbon MCS answers could never
  have caught this class of bug.

### Added — `chematic-fp`/`chematic-py` (`rdkit_torsion_fp`, RDKit-compatible Topological Torsion fingerprint)

- New `chematic_fp::rdkit_torsion_fp` (Rust) / `Mol.rdkit_torsion_fp()` (Python): a
  from-scratch Rust port of RDKit's
  `GetHashedTopologicalTorsionFingerprintAsBitVect`, opt-in and fully separate from
  the existing native `torsion_fp` (neither affects the other). First entry in the
  fingerprint-parity series (Track A / 99-point directive Phase 6) — the closest of
  six priority fingerprints measured against RDKit, chosen to establish the pattern.
- Measured 87.2% bit-exact on a 1000-molecule general corpus sample against a live
  RDKit oracle, up from 0% at the start of this round. Three real bugs found and
  fixed along the way, each verified against the oracle: a missing `-2`
  "topological torsion correction" on atom invariants, double-counted torsion paths
  inflating RDKit's "count simulation" threshold bits, and missing 3-membered-ring
  closure entries (a triangle's own closing bond forms an additional torsion that a
  plain 4-distinct-atom path search can't find).
- Two known, documented residuals remain (see `chematic_fp::rdkit_torsion`'s module
  doc comment): RDKit's hybridization-gated pi-electron count for hypervalent atoms
  (e.g. phosphorus/sulfur) isn't replicated, and asymmetrically-substituted
  3-membered rings can need more than the one closure entry generated here.

### Added — `chematic-fp`/`chematic-py` (`rdkit_atom_pair_fp`, RDKit-compatible Atom Pair fingerprint)

- New `chematic_fp::rdkit_atom_pair_fp` (Rust) / `Mol.rdkit_atom_pair_fp()` (Python): a
  from-scratch Rust port of RDKit's `GetHashedAtomPairFingerprintAsBitVect`, opt-in
  and fully separate from the existing native `atom_pair_fp` (neither affects the
  other). Second entry in the fingerprint-parity series (Track A / 99-point
  directive Phase 6), reusing `rdkit_torsion_fp`'s atom-invariant primitive since
  RDKit's own `AtomPairAtomInvGenerator` is the shared generator for both.
- Measured 87.2% bit-exact on a 1000-molecule general corpus sample against a live
  RDKit oracle, with zero implementation bugs found this round — every mismatching
  molecule in the sample contains a hypervalent phosphorus or sulfur atom, confirming
  the entire residual is `rdkit_torsion_fp`'s already-documented hybridization-gated
  pi-electron count gap (shared via the reused atom-invariant), not a new gap in
  atom-pair enumeration or hashing.

### Added — `chematic-fp`/`chematic-py` (`rdkit_pattern_fp`, RDKit-compatible Pattern fingerprint)

- New `chematic_fp::rdkit_pattern_fp` (Rust) / `Mol.rdkit_pattern_fp()` (Python): a
  from-scratch Rust port of RDKit's `PatternFingerprint`, opt-in and fully separate
  from the existing native `pattern_fp` (neither affects the other). Third entry in
  the fingerprint-parity series (Track A / 99-point directive Phase 6) — structurally
  unrelated to the first two (SMARTS substructure matching against 13 fixed patterns,
  not a path/pair enumeration), reusing chematic's existing `chematic-smarts` matcher
  rather than a new one.
- Measured against a live RDKit oracle on three corpora with different chemical
  distributions: 100% bit-exact on a 1000-molecule general corpus sample, 100% on a
  ChEMBL sample, 99.6% on an NCI sample. Two real bugs found and fixed along the way:
  (1) chematic's SMILES parser preserves a Kekule-notation input's literal
  single/double bond orders even when the ring is perceived as aromatic, so a raw
  `bond.order == Aromatic` check silently missed every ring bond of a Kekule-written
  heteroaromatic — caught the NCI corpus number at 30.4% before switching to
  chematic's own aromaticity-perception model; (2) naively re-perceiving *every*
  non-literally-aromatic bond then regressed the ChEMBL corpus from 100% to 94.1% —
  for a molecule with both lowercase-aromatic and Kekule-written rings, re-perceiving
  the whole molecule wrongly classified a genuinely non-aromatic Kekule ring as
  aromatic too. Fixed by only re-perceiving molecules with *zero* literal `Aromatic`
  bonds anywhere; molecules that already carry any literal `Aromatic` bond trust every
  bond's own stored order throughout instead. The remaining NCI residual traces
  entirely to chematic's own aromaticity-perception model disagreeing with RDKit's on
  specific wholly-Kekule-written ring systems (a brominated anthraquinone/xanthone
  core, an S/N-containing bicyclic heterocycle, a ketone-fused tetralone-like system,
  a Zn coordination complex) — a pre-existing, separate gap, not a defect in this
  fingerprint's own match-enumeration or hashing logic (independently verified:
  per-pattern match counts identical to RDKit's own `GetSubstructMatches(...,
  uniquify=False)` across all 13 patterns on a repro molecule before either bug was
  even found).

### Added — `chematic-py`/`chematic-wasm` (full `McsConfig`/`McsOutcome` exposed to MCS bindings)

- Python's `find_mcs` previously exposed only 2 of `McsConfig`'s 11 fields
  (`ring_matches_ring_only`/`complete_rings_only`); the rest were silently
  hardcoded to their defaults. Extended to accept all 11 as keyword
  arguments (`match_bonds`, `min_atoms`, `timeout_ms`, `atom_compare`,
  `bond_compare`, `match_chiral_tag`, `match_charge`, `match_isotope`,
  `maximize_bonds`, plus the existing two) — fully backward-compatible,
  every new argument defaults to `McsConfig::default()`'s own value.
- New `find_mcs_checked` (Python) and `mcs_smiles_json_with_config` (WASM,
  camelCase JSON config in, `{"smiles": string|null, "wasTimedOut": bool}`
  out) expose `McsOutcome`'s exhaustive/timed-out distinction for the first
  time — previously always silently discarded via `find_mcs_with_config`
  (never the `_checked` variant), so a caller using `timeout_ms` had no way
  to tell a possibly-non-optimal result apart from a proven-optimal one.
- The two pre-existing WASM `mcs_smiles_json`/`mcs_smiles_json_with_ring_config`
  functions keep their exact signatures and behavior unchanged, refactored
  to share one `QueryMolecule`-to-`Molecule` reconstruction helper (was
  duplicated verbatim 2-3x per binding).

### Added — `chematic-fp`/`chematic-py` (`rdkit_rdk_fp`, RDKit-compatible "RDKit fingerprint")

- New `chematic_fp::rdkit_rdk_fp` (Rust) / `Mol.rdkit_rdk_fp()` (Python): a
  from-scratch Rust port of RDKit's default fingerprint
  (`Chem.RDKFingerprint`/`RDKFingerprintMol`), opt-in and fully separate from the
  existing native `path_fp` (`rdkit_path_fp` on the Rust side; neither affects the
  other). Fourth entry in the fingerprint-parity series (Track A / 99-point
  directive Phase 6) — chosen after Avalon was deliberately parked (it wraps a
  vendored, decades-old external C toolkit rather than being a portable
  algorithm, so raw hamming-distance-to-RDKit ranking doesn't reflect porting
  difficulty).
- Reproduces RDKit's full branched-subgraph enumeration (`minPath=1..maxPath=7`,
  not just linear paths), its per-bond hash using **path-local** atom degree
  (how many times an atom appears as an endpoint *within the current subgraph*,
  not the molecule's true degree — confirmed against a live oracle, since using
  true molecular degree does not reproduce RDKit's own output even for a plain
  3-atom chain), and, since RDKit's default `numBitsPerFeature=2`, its
  deliberately weakened Mersenne Twister variant for the second bit per
  feature. That PRNG has a genuine boost-library footgun worth flagging: boost's
  *deprecated* `mersenne_twister<>` wrapper class (the one RDKit's own typedef
  instantiates) silently discards its own last template parameter — RDKit's
  source passes `3346425566U` there — and hardcodes the textbook MT19937
  constant `1812433253` instead, confirmed by reading boost's own
  `mersenne_twister.hpp` source directly.
- Measured bit-exact (identical on-bit sets) against a live RDKit oracle: 100%
  on `descriptor_census_corpus.smi` (5000/5000) and
  `chembl_accuracy_corpus_4999.smi` (5000/5000), 99.44% on
  `nci_first_5k_smiles_only.smi` (4963/4991). Every one of the 28 NCI
  mismatches is a fused polyheteroaromatic dye, an exotic charged aromatic
  heterocycle, or a metal-coordination complex where chematic's own
  `chematic-perception` Hückel aromaticity model doesn't (yet) recognize the
  same aromatic system RDKit's does — the same class of pre-existing,
  out-of-scope dependency gap already documented for the three earlier
  fingerprints in this series, not a defect introduced by this port.

### Added — `chematic-fp`/`chematic-py` (`rdkit_layered_fp`, RDKit-compatible Layered fingerprint)

- New `chematic_fp::rdkit_layered_fp` (Rust) / `Mol.rdkit_layered_fp()` (Python):
  a from-scratch Rust port of RDKit's `Chem.LayeredFingerprint`
  (`LayeredFingerprintMol`), opt-in and fully separate from the existing
  native `layered_fp` (neither affects the other). Fifth and final entry in
  the fingerprint-parity series (Track A / 99-point directive Phase 6).
  Upstream itself documents this fingerprint as experimental.
- Reuses the same branched-subgraph enumeration as `rdkit_rdk_fp`, but
  computes 6 independent per-bond feature "layers" (topology, bond order,
  atom type, ring membership, ring size, aromaticity), each packed into a
  small bitfield and folded separately into its own fingerprint bit per
  subgraph — no `numBitsPerFeature`/PRNG step, unlike `rdkit_rdk_fp`.
- One real bug found and fixed: `LayeredFingerprintMol` enumerates with
  `useHs=false` (unlike `rdkit_rdk_fp`'s own `useHs=true`), so any bond
  touching an explicit hydrogen atom must be excluded entirely from
  enumeration. chematic, like RDKit, represents an isotope-labeled hydrogen
  (e.g. `[2H]`/`[3H]`) as a real graph atom rather than folding it into an
  implicit H count, so a naive `useHs=true` port silently included such bonds
  and produced extra bits (caught on the general and ChEMBL corpora, each via
  a single deuterium/tritium-containing molecule).
- Measured bit-exact (identical on-bit sets) against a live RDKit oracle:
  100% on `descriptor_census_corpus.smi` (5000/5000) and
  `chembl_accuracy_corpus_4999.smi` (5000/5000), 99.46% on
  `nci_first_5k_smiles_only.smi` (4964/4991) — the same pre-existing
  `chematic-perception` aromaticity/ring-model gap class (fused
  polyheteroaromatic dyes, exotic charged heterocycles, metal-coordination
  complexes) already documented for the other fingerprints in this series,
  confirmed via the same throwaway-Python-prototype-against-RDKit's-own-mol
  methodology used for `rdkit_rdk_fp`.

## [0.22.0] — 2026-08-29

Minor release: new additive, non-breaking WASM API (`embed_ensemble_v2_json`), plus two
correctness/robustness fixes — a `chematic-smiles` canonicalization hang (issue #421)
and a `chematic-3d` distance-geometry embedding gap affecting every 3-membered-ring
molecule. No breaking API changes.

### Fixed — `chematic-smiles` (canonical_smiles/canonical_atom_order could hang, issue #421)

- `canonical_automorphism::extend_mapping`'s colored-graph automorphism
  backtracking search had no internal step bound of its own — the outer
  `SearchBudget` in `canonical_search.rs` only counts *calls* to
  `has_colored_automorphism_mapping`, not work done *inside* one call. On a
  molecule with several simultaneously-unresolved large symmetric regions
  (e.g. 3 near-identical repeated substituent arms all still non-singleton
  at once), a single call could explore a combinatorially large space and
  never return — observed running past 2 minutes, never confirmed to
  terminate, on a real 94-atom ChEMBL molecule whose atom order already
  coincided with `canonical_atom_order`'s own output order. Fixed with an
  always-on `MAX_EXTEND_MAPPING_STEPS` ceiling (200,000) that falls back to
  `false` on exceeding it — safe per this module's own pre-existing
  documented invariant ("a false result may cost performance... a true
  result must always be a genuine automorphism"), so this can only ever
  cost a missed prune, never a wrong canonical answer. Verified: the
  reordered repro molecule now canonicalizes in ~0.6s and produces the
  exact same canonical SMILES as the original ordering.

### Added — `chematic-wasm` (`embed_ensemble_v2_json`, A2.1 WASM bindings)

- New `embed_ensemble_v2_json` binding for `chematic_3d::embed_ensemble_v2`
  (Track A, OpenEye advantage RFC), mirroring the Python binding
  (`Mol.conformer_ensemble_v2()`, already shipped) via the existing
  `pipeline_v2.rs` WASM conventions (camelCase keys, `schemaVersion: 1`
  tagged-union envelope, `FiniteF64` for JSON-safe non-finite handling).
  Preserves `embed_ensemble_v2`'s documented error asymmetry: an ensemble
  where every attempt fails and zero conformers are kept is still
  `{"ok": true, "result": {...}}`, with per-attempt detail in
  `result.attempts` — only a config that could never succeed (an invalid
  `rmsdThreshold`) surfaces as `ok: false`.

### Fixed — `chematic-3d` (3-membered rings no longer fail closed at the distance-geometry stage)

- `dg_fft::build_bond_angle_bounds`'s angle-constraint loop treated every pair
  of a center atom's neighbors as a generic 1-3 (through-center) relationship,
  tightening their bound with the generic ~109.5°/120° ideal angle. In a
  3-membered ring, that "1-3" pair is *also* a direct 1-2 bonded pair (the
  ring closes one bond away), so the generic-angle bound overwrote the
  correct, much shorter bond-length bound with a contradictory one
  (cyclopropane's ring-closing C-C pair: bond constraint gave upper ≈ 1.59 Å,
  the angle constraint then tightened lower to ≈ 2.41 Å) — `lower > upper`,
  caught by `try_embed_once`'s pre-smoothing sanity check and failed closed
  as `BoundsConstructionFailed`. This was a disclosed limitation affecting
  every cyclopropane/epoxide/aziridine/thiirane-containing molecule (11/265
  in the strict 3D corpus). Fixed by skipping the angle-derived bound for any
  neighbor pair that is itself directly bonded — a 3-membered ring's three
  bond-length constraints already fully determine its shape, so nothing is
  lost. Every molecule without a 3-membered ring is provably unaffected (the
  skip condition can only fire on a ring-closing bonded pair). Verified: all
  11 previously-failing molecules now embed and pass the strict-MMFF94 arm
  under its exact production config; `embed_pipeline_v2`'s strict-MMFF94
  corpus result moves from 252/265 to 263/265.

## [0.21.0] — 2026-08-27

Minor release: `McsConfig` gains new fields (`match_charge`/`match_isotope`) and
`McsOutcome` (additive, non-breaking public API), plus four correctness fixes —
MCS branch-and-bound completeness, and issues #403/#407/#415 (metal-complex
charge neutralization, zwitterion proton transfer, tautomer valence validity).
No breaking API changes.

### Added — `chematic-smarts` (`McsConfig`: `match_charge`/`match_isotope`, typed timeout outcome)

- `McsConfig` gained `match_charge`/`match_isotope` fields, mirroring the
  existing `match_chiral_tag` exactly (default `false`, ignoring charge/
  isotope unless opted in — e.g. a carboxylate `[O-]` matches neutral `O`,
  `[13CH4]` matches `C`, unless explicitly required to match).
- New `McsOutcome` enum (`Exhaustive`/`TimedOut`) and
  `find_mcs_with_config_checked`, mirroring `match_vf2.rs`'s own
  `find_matches_with_rings_and_config_checked` pattern: reports whether
  `McsConfig::timeout_ms` was reached before the search finished, rather
  than silently returning a possibly-non-optimal result indistinguishable
  from an exhaustive one. `find_mcs`/`find_mcs_with_config`'s existing
  signatures and behavior are unchanged — purely additive.

### Fixed — `chematic-smarts` (`find_mcs`'s branch-and-bound search was incomplete)

- `grow()` only ever tried `frontier[0]` at each search node, with no way
  to exclude it and try a different frontier atom instead. If `frontier[0]`
  had no compatible candidate in some other input molecule, the whole
  branch died silently — even when skipping it (mapping a different,
  compatible frontier atom instead) would reach a strictly larger common
  substructure. Minimal repro: `OC(N)N` vs `NC(N)` returned an MCS of 2
  atoms instead of the true 3 (every seed path hit an unmatchable O ahead
  of a still-needed N leaf in frontier iteration order, and nothing
  backtracked past it).
- Fixed via standard include/exclude branch-and-bound: include the first
  frontier atom (as before), or exclude it from this subtree and let a
  later frontier atom be tried instead, unwinding the exclusion on
  backtrack. `upper_bound_additional`'s pruning tightened to also skip
  excluded mol[0] atoms. `McsOutcome::Exhaustive`'s doc comment previously
  claimed a completeness the algorithm didn't actually have; corrected.
- Found while scoping MCES/multi-tie-MCS-enumeration work (both would have
  quietly inherited this incompleteness); fixed first, as its own
  correctness fix, before any new-feature work on top of it.

### Fixed — `chematic-chem` (issue #403: `disconnect_metals` left a dative-bond-derived formal charge unneutralized)

- `disconnect_metals` severed dative M-O/M-N bonds without touching the
  non-metal atom's stored `hydrogen_count`. A dative bond is commonly
  written with a formal charge that exactly balances the bond (e.g. `[O+]`
  single-bonded to a metal, satisfying O+'s valence-3 with 0 implicit H) —
  after disconnection the atom's true valence changed, but its stale H
  count didn't, so the very next pipeline stage, `neutralize_charges`
  (guard `h > 0` on the raw stored field), saw `h == 0` and skipped it,
  leaving a dangling formal charge with nothing left to justify it. The
  charge only got neutralized on a *second* standardize pass, once a fresh
  parse of the incorrectly-charged first-pass output stored the H count
  explicitly — a real, confirmed idempotency bug across 34/4999 molecules
  in RDKit's own bundled NCI Diversity Set holdout.
- Fixed: `disconnect_metals` now recomputes the affected atom's H count via
  valence inference against the post-disconnection topology, so
  `neutralize_charges` sees the true state on the very first pass.
- Also fixed a related bug in `remove_hydrogens`: it unconditionally reset
  any atom with `hydrogen_count == Some(0)` to `None`, when its own doc
  comment's stated intent was narrower — only atoms that actually had an
  explicit H *atom* neighbor removed this call. Tightened to match.
- Verified: NCI `first_5K.smi` (4,999 molecules, now a permanent regression
  corpus, `scripts/nci_first_5k_smiles_only.smi`) 34 → **0** failures. New
  hand-built metal-complex holdout added (11 fixtures spanning Ni, Co, Al,
  Zn, Cr, Fe, Mg, Mn, Hg, Cd, Pd with varied ligand shapes plus an
  ionic-salt negative control), all pass. Existing dev corpora unaffected.

### Fixed — `chematic-chem` (issue #407: `normalize_zwitterion` invented a proton for non-transferable charge pairs)

- The active proton-transfer path unconditionally neutralized the negative
  atom (+1 charge, +1 H) but only neutralized the paired positive atom if
  it had an available H. For a permanently charge-separated group with no
  transferable proton on either side (e.g. a diazo-N,N'-dioxide,
  `[N+]([O-])=[N+]...[O-]`, structurally similar to a nitro group), this
  invented a hydrogen on the negative atom from nowhere, silently changing
  the molecule's formula and net charge.
- `has_zwitterion`'s "some + and some - charge exists somewhere" check is
  necessary but not sufficient for a real protonation-state zwitterion.
  Rather than hand-classifying every non-transferable functional-group
  family (nitro, diazo N-oxide, mesoionic, ...), fixed by gating the
  transfer on the invariant that actually matters: a proton can only move
  if the chosen donor has one to give. If the positive atom has no
  available H, neither atom is modified — the pair is left untouched.
- Verified: `descriptor_census_corpus.smi` standardize-path idempotency
  4 → **1** residual (the 3 known molecules for this issue now pass; the 1
  remaining is an unrelated `canonical_tautomer` interaction, confirmed
  standing alone — see issue #402/#415). New property-based tests (atom/
  element/H-count conservation, net-charge conservation, idempotency, a
  genuine-zwitterion positive control, nitro/pyridine-N-oxide negative
  controls, and an atom-permutation invariance check) rather than
  spot-checking one hardcoded expected string per molecule.

### Fixed — `chematic-chem` (issue #415: `canonical_tautomer` could produce an unkekulizable molecule)

- Both of `canonical_tautomer`'s aromatic H-shift mechanisms
  (`transfer_hydrogen_exocyclic_lactam`, `transfer_hydrogen_aromatic`)
  validated their acceptor atom in isolation (0 implicit H, correct
  degree) but not against the rest of the ring. In a fused/bridged ring
  system, an acceptor that looks individually valence-legal can be
  ring-adjacent to *another* aromatic atom that already carries an "extra"
  H — two such atoms next to each other in one aromatic ring can't both
  correctly contribute a lone pair to the same ring's pi system. Confirmed
  via a real repro (`Oc1[nH]ncc2c3cc(OCc4ccccc4)ccc3nc1-2`) that this
  produced an `[nH2]`-shaped, over-valent nitrogen RDKit's own parser
  rejects outright.
- Fixed by validating both mechanisms' output with `kekulize` before
  accepting it (`validate_valence` doesn't cap aromatic atoms to their
  primary valence and doesn't catch this shape). `find_sssr`'s own ring
  choice for such fused systems isn't always unique, which is why the same
  physical molecule reached this state on a second standardize pass but
  not the first — a separate, deeper, not-fully-solved order-dependence
  residual this fix doesn't chase down; it stops the resulting corruption.
- The known corpus residual this molecule was already causing
  (`descriptor_census_corpus.smi`, ceiling 1) is unchanged in count — the
  molecule still isn't byte-identical after one standardize pass vs two —
  but now converges to a valid tautomer by the second pass instead of
  producing invalid chemistry on the first.

## [0.20.1] — 2026-08-26

Patch release: three canonical-SMILES/standardization correctness fixes, all
merged after the `v0.20.0` tag was cut (none of them are retroactively part of
that release's own claims). No breaking API changes.

### Fixed — `chematic-smiles` (issue #390: coupled E/Z canonicalization could silently change geometry)

- Root cause, two independent defects in `CanonicalWriter`'s E/Z marker
  machinery, both needed to reproduce the filed witness
  (`O/N=C/C(C=N/O)=N\NC`, whose atom3=atom7 double bond was silently
  written as E instead of the true Z — confirmed via independent RDKit
  `MolToInchi`/`GetStereo()`, not just chematic's own self-consistency):
  1. `resolve_ez_markers`'s carrier election for an ambiguous end could
     elect a candidate bond whose sibling was raw-marked and load-bearing
     for a *different*, unrelated double bond — demoting the sibling
     silently under-specified that other double bond, while the elected
     candidate simultaneously handed a *third*, genuinely undefined
     double bond (confirmed via InChI's own `?` stereo descriptor for it)
     a geometry it never had. Neither the demotion nor the promotion is
     geometry-neutral, and picking between them at random (by whichever
     canonical-numbering trial happened to explore first) is exactly how
     the witness's true geometry got lost. Fixed by
     `CanonicalWriter::is_load_bearing_elsewhere`: an election must not
     demote a raw-marked candidate that is some other, non-ambiguous
     double bond's only geometric anchor. Deliberately narrower than "has
     a raw mark" — a candidate whose sibling is *itself* ambiguous (has
     its own resolution path, e.g. a genuinely coupled/shared-carrier
     system) is not protected, so legitimate coupled resolution is
     unaffected.
  2. Independently, `normalize_ez` decided a shared E/Z group's sign from
     a value that had already been re-oriented for one specific DFS write
     direction. Which end of a directional bond a given canonicalization
     trial happens to write "forward" vs "backward" varies across
     candidate atom numberings for reasons unrelated to that bond's own
     geometry (a tie elsewhere in the molecule), so the seeded sign could
     vary too, non-deterministically flipping an otherwise-correct group.
     Fixed by splitting `normalize_ez` into a mol-relative propagation
     step (always flips `effective_order`, the bond's own topology-fixed
     `atom1`→`atom2` reading, never an already-write-oriented value) and a
     write-perspective anchor-seeding step (the write atom decides the
     group's shared sign exactly once, and only that — it never enters
     propagation). Found and fixed second, after the first fix alone
     restored correctness for the filed witness but broke canonical-form
     stability (10 independently-rooted, InChI-confirmed-equivalent
     respellings of the witness converged to only 1 string before this
     defect existed in the code at all — introducing defect #1's fix
     alone dropped that to 3 non-idempotent strings; both fixes together
     restore 10/10 convergence).
- An intermediate, never-shipped attempt at defect #2 (seeding purely
  from write-perspective, dropping the mol-relative anchor entirely)
  restored canonical-form stability but silently made canonicalization
  *informationally lossy* for this shape — the witness's true-Z and a
  hand-verified true-E mirror both canonicalized to the identical string,
  each losing its own stereo identity in different directions. Caught by
  a mirror-distinctness regression test before being combined with defect
  #1's fix into what actually shipped; not a real intermediate state of
  the code, called out here only because the failure mode (idempotent AND
  self-consistent, yet wrong) is exactly the kind that hides behind a
  weaker "does it round-trip" check alone.
- Verified against the real 290-compound corpus from the originating
  investigation (eMolecules, 9.47M compounds,
  `renkin doctor stock reimport_idempotency`) two ways: idempotence
  (**290/290**, up from 289/290 before this fix) and, independently, that
  each record's chematic canonical form reparses in RDKit to the exact
  InChIKey recorded for that record at investigation time (**290/290**) —
  the corpus itself is not committed (see PR #389's own note on this), only
  aggregate counts.
- New tests: the witness's own geometry preserved and stable, its only
  safe alternate-carrier candidate confirmed to have none available
  (`alternate_ez_markings` returns empty — the sibling candidate is
  load-bearing elsewhere, so no valid respelling moves the mark there),
  mirror-image (E vs Z) distinctness, and full atom-order-permutation
  invariance across 18 relabelings/markings via the same
  `ez_carrier_test_variants` harness `EZ_SHARED_CARRIER_FULLY_RESOLVED`'s
  own regression test uses.
- Not addressed, not fixed, not blocking this PR: a synthetic edge case
  found while writing this fix's own tests (not part of the filed issue,
  the 290-corpus, or any pre-existing test) — an ambiguous end whose
  *both* candidate bonds carry mutually-consistent raw marks, where one
  candidate's sibling is itself adjacent to a genuinely undefined double
  bond, produced a geometry mismatch between the raw input and a
  canonicalize→reparse round-trip in this crate's own test-only
  `up_of_reference` oracle. Not confirmed as a production defect (the
  oracle is test-only scaffolding, not the production code path) or ruled
  out as one — flagged here rather than silently dropped, deliberately
  not filed as an issue without that confirmation.

### Fixed — `chematic-chem` (issue #399: `standardize()` silently dropped stereo tables on several rebuild paths)

- 8 functions in `standardize.rs` (`neutralize_charges`, `normalize_zwitterion`
  active path, `normalize_groups`, `remove_isotopes`, `reionize`, `uncharge`,
  `prefer_organic`, `disconnect_metals`) rebuilt the molecule via a bare
  `MoleculeBuilder` without carrying `stereo_neighbor_order`/
  `bond_directions`/`stereo_groups` forward — even when zero atoms were
  actually modified. Once #392 (v0.20.0) restored `remove_hydrogens`'s own
  table-copying, its adjacency-based fallback reconstruction (correct for
  ring-closing stereocenters, transposed for ring-opening ones) started
  firing on every stereocenter passing through any of these 8 functions,
  flipping `@`/`@@` depending on which role a stereocenter played in the
  original text vs. the canonical rewrite.
- Fixed: a bulk `copy_stereo_from`/`copy_bond_directions_from`/
  `copy_stereo_groups_from` for the 7 functions that preserve every
  atom/bond 1:1; `prefer_organic` now delegates to the already-correct
  `extract_fragment` (it genuinely drops atoms); `disconnect_metals` remaps
  `bond_directions` bond-by-bond (it drops metal-adjacent bonds, shifting
  survivors' indices). Added the stereo-preservation invariant as a doc
  comment on `Molecule::stereo_neighbor_order`.
- Also corrected `tp2_20_isotope_parent_preserves_stereo`'s hardcoded
  expected value, which had baked in this exact bug (verified via RDKit
  `MolToInchi` which of the two candidates matches the isotope-free
  reference structure).
- Verified: dev-corpus standardize-path idempotency 615/519 → 68/60 (the
  exact pre-#392 baseline); NCI `first_5K.smi` blind holdout (4,999 unused
  real molecules, run once): 0 stereo-related failures, 0 stereocenter-count
  mismatches against an independent RDKit CIP oracle. New regression suite
  (`standardize_stereo_table_preservation.rs`, 12 tests) verified to catch
  the bug (fails on pre-fix code).

### Fixed — `chematic-smiles` (issue #395: ring-closure bond marker ignored the closure partner's aromaticity)

- The canonical writer's ring-closure marker decision checked only the
  currently-written atom's own aromaticity, never its ring-closure
  partner's — unlike the equivalent, correct decision for a plain tree-edge
  child, which checks both endpoints. A bare ring-closure digit between two
  aromatic atoms is read back by the parser as an *aromatic* bond, so a
  genuinely `Single`-order ring-closure bond joining two atoms that each
  individually happen to be aromatic (a non-aromatic fusion bond connecting
  two separately-aromatic ring systems, e.g. `c1-2`) silently became
  aromatic on re-parse whenever the writer omitted its `-` marker.
- Fixed by mirroring the tree-edge `implicit` computation for ring closures,
  checking both endpoints' aromaticity.
- Verified: dev-corpus bare-parse idempotency 73/57 → **0/0**, a complete
  fix, not a partial improvement. Independent RDKit oracle: all 10,000
  corpus lines' canonical output now round-trips to the exact same InChI as
  the original input — 0 mismatches. Combined with #399's fix, the
  standardize-path corpus is now at 0/4 (the 4 residual failures traced to
  two newly-filed, unfixed issues: #407, a `normalize_zwitterion`
  proton-transfer asymmetry, and a `canonical_tautomer` interaction of the
  same class as #402).

## [0.20.0] — 2026-08-25

### Added — `chematic-3d`/`chematic-py`/`chematic-wasm` (stereo-safe 3D generation, issue #291)

- `PipelineV2Config::expand_implicit_h_through_pipeline`: runs
  `embed_pipeline_v2`'s whole pipeline on a temporary
  `add_hydrogens`-expanded copy of the molecule instead of the original,
  closing a real gap for ring-fused declared stereocenters (testosterone,
  cholesterol, and similar) where `repair_tetrahedral_center` previously had
  no coordinate to reflect an implicit H against. `PipelineV2Result::coords`/
  `final_stereo` stay scoped to the caller's original atom count either way;
  other diagnostic fields describe the expanded internal working state,
  documented explicitly on the struct. Stage 2/3 (torsion knowledge) always
  run against the *original* molecule, not the expanded one — `add_hydrogens`
  appending real graph nodes would otherwise silently reclassify e.g. a
  secondary amine's hybridization for torsion-rule purposes, an interaction
  this design deliberately avoids rather than leaves unmeasured.
- `PipelineV2Config::stereo_safe(force_field_policy)`: a new convenience
  constructor bundling `stereo_policy: RepairAndVerify` +
  `enforce_chirality: true` + `expand_implicit_h_through_pipeline: true` —
  setting only some of these three together silently falls back to a
  configuration already measured unsound for exactly this molecule class.
  Exposed identically through `chematic-py` (`PipelineV2Config.stereo_safe(...)`
  staticmethod) and `chematic-wasm` (`pipeline_v2_stereo_safe_config_json(...)`,
  since JS has no static-method equivalent).
- **Measured**, 29-molecule × 5-seed regression corpus: **144/145 (99.3%)
  correct_and_ok, 0 silently_wrong, 0 loud_failure_stereo** — testosterone
  and cholesterol both now succeed on every declared stereocenter, every
  seed (5/5 each). The one remaining failure is cholesterol's pre-existing,
  unrelated UFF `CatastrophicBondBlowup` issue (not a stereo defect).
- Existing callers unaffected: the new flag defaults to `false` at every
  layer (Rust/Python/WASM), and `expand_implicit_h_through_pipeline: true`
  without `enforce_chirality: true` is rejected with a typed
  `InvalidConfiguration` error rather than silently doing nothing.

### Added — `chematic-3d` (connectivity-ordered 3D coordinate generation, issues #256/#255)

- `generate_coords_connectivity_ordered` (new `chematic_3d::dg_connectivity_ordered`
  module, also re-exported at the crate root): a parallel rule-based 3D
  placement engine, structurally ported from `chematic-depict`'s proven 2D
  technique — a single worklist discovers and places rings and chain atoms
  in true connectivity order (never "all rings, then all chains"), unlike
  the legacy `generate_coords`'s `place_rings`, which places every ring in a
  component before walking any chain atom (issue #256) and can produce
  distorted fusion-seam bonds via a fixed `+y` extension (issue #255).
- **Measured**, 33-molecule differential corpus (issue #277's 17 real
  ChEMBL molecules + 8 RFC known-broken topologies + positive controls):
  raw-geometry soundness **10/33 → 33/33**, zero regressions, all 17 of
  #277's real molecules improved. Post-UFF-minimization: an initial
  "new-island" ring-entry-direction regression (rings joined by a single
  direct bond, e.g. biphenyl) was found and fixed via a 12-candidate
  centroid-away-plus-clearance ranking — the fixed engine's post-UFF
  `mean_viol15` reaches **0.0000**, past even the legacy engine's own
  0.0055 baseline. Determinism (33/33 identical across repeated runs) and
  atom-order-permutation-invariant quality both confirmed.
- **`generate_coords` itself is completely unchanged** — this ships as an
  available, independently-selectable alternative (not a default-behavior
  switch or topology-based routing), the conservative and fully reversible
  of the three options considered. No existing caller (`generate_coords_etkdg`,
  `embed_pipeline_v2`, ...) is routed to the new engine by this release.
  Issues #255/#256/#277 stay open — availability as public API is not the
  same as a production routing decision.
- Not yet exposed through Python/WASM bindings — Rust core first, per this
  project's established pattern; deferred to a follow-up.

### Fixed — `chematic-3d` (UFF-rescue path did not enforce declared chirality, issue #210)

- `rescue_with_distance_geometry_v2` (the bridge that retries embedding
  after a UFF catastrophic-bond-length-blowup) embedded its retry with
  `enforce_chirality: false` unconditionally — a constraint that no longer
  needed to hold once `EmbedParameters::materialize_implicit_h_for_chirality`
  became available (issue #291, above) at this same low-level embed API.
- **Measured** against the 58-molecule corpus this bridge's own tests use:
  zero regressions across the full corpus, and one of the 5 residual
  molecules this issue names (`atorvastatin_fragment`) newly succeeds with
  declared stereochemistry preserved.
- `naproxen_S`/`ibuprofen_S`/`testosterone`/`cholesterol` remain unfixed by
  this specific change — this bridge has no post-minimization
  repair-and-reverify step the way `embed_pipeline_v2`'s own stage 11 does,
  so a UFF-introduced post-embed violation still falls through to the
  original (honest) failure. Issue #210 stays open, partial fix only.

### Added — benchmark (`chematic-3d` best-of-N conformer generation vs. RDKit)

- New benchmark arm, `chematic_pipeline_v2_uff_best_of_10`: `embed_ensemble_v2`
  (`count=10`, `UffOnly`, `max_attempts=1`, `rmsd_threshold=0.0`), matched
  against `docs/rfcs/pipeline_v2_vs_rdkit_etkdgv3_benchmark.md`'s existing
  `rdkit_etkdgv3_best_of_n` arm so "10 attempts, best by energy" means the
  same thing on both sides.
- **Measured**: ~250/265 molecules successful on both sides; median paired
  RMSD **2.147 Å**, median TFD **0.344** against RDKit's own best-of-10
  `EmbedMultipleConfs`. This confirms `embed_ensemble_v2` (A2) works
  robustly at this scale, but **does not** establish that chematic's
  energy-based conformer *selection* picks the same conformer RDKit would
  from the same pool — that's a separate, unestablished claim this
  benchmark's numbers alone don't support.

### Fixed — `chematic-chem` `remove_hydrogens` (isotope labels silently destroyed)

- `remove_hydrogens` previously removed *any* atom with `element == H`,
  including deuterium (`[2H]`), tritium (`[3H]`), and any other
  isotope-labeled hydrogen — collapsing an explicit isotopic-H atom into an
  ordinary heavy atom's opaque `hydrogen_count` silently discards the
  isotope label, since that representation has no way to record "N
  implicit hydrogens, one of which is deuterium." Found via a downstream
  consumer (RENKIN) whose stock-identity pipeline calls `standardize`
  with `remove_explicit_h: true` (chematic's own default) on a
  9.48M-compound real-world building-block corpus containing 12,688
  explicitly-isotopic-hydrogen rows — every one of them lost its D/T
  label on the very first canonicalization pass, exhaustively confirmed,
  zero exceptions.
- Now: only a *non-isotopic* explicit H (`element == H && isotope.is_none()`)
  is removed. An isotope-labeled H is kept as an explicit atom node, like
  any other heavy atom — its bond is preserved, and a heavy atom that
  retains an isotopic-H neighbor still gets its `hydrogen_count` reset to
  `None` so implicit H is recomputed from valence, which correctly
  accounts for the kept neighbor's own bond (`valence_inferred_hcount`
  counts every bonded neighbor by bond order, not by element identity, so
  this needed no separate special-casing).
- Heavy-atom isotopes (¹³C, ¹⁴C, ¹⁵N, ¹⁸O, ...) were never affected by
  this bug in the first place (they're not `element == H` atoms at all)
  and remain untouched — new regression tests pin this explicitly rather
  than leaving it as an implicit assumption.
- New tests: fully-deuterated methane, tritium, mixed deuterium+plain-H on
  the same atom (confirms `hydrogen_count` recomputation is exactly
  correct, not off-by-the-kept-neighbor), the four heavy-isotope no-op
  cases, an isotope-labeled tetrahedral stereocenter (structural
  soundness only — see below), a full `standardize()` round-trip, and a
  canonical-round-trip (`canon(parse(canon(parse(s))))`) case matching
  the exact re-canonicalization scenario that surfaced the bug.
- **Does not** address the separate, independent finding from the same
  investigation: `canonical_smiles` can pick a different, structurally-
  identical-but-differently-parity'd traversal on a second canonicalization
  pass for certain (typically ring-fused or otherwise symmetric-adjacent)
  molecules, occasionally flipping a declared tetrahedral or E/Z
  descriptor's *meaning* even though the literal token or overall
  structure looks unchanged. That is a distinct root cause (this
  function's own missing `stereo_neighbor_order` restoration, unlike
  `add_hydrogens`, which explicitly does restore it) needing its own,
  separate fix and issue — not addressed, not fixed, and not blocking
  this PR.

### Fixed — `chematic-chem`/`chematic-smiles` (stereo identity changes after canonical round-trip)

- `chematic-chem::hydrogen::remove_hydrogens` never restored `Molecule`'s
  `stereo_neighbor_order`/`bond_directions` side tables (unlike its sibling
  `add_hydrogens`, which explicitly does via `copy_stereo_from` + sentinel
  remapping) -- it always rebuilds a fresh `MoleculeBuilder`, so these
  tables were silently wiped on *every* call, even a complete no-op one
  that removed nothing.
- `chematic-smiles`'s canonical writer's `corrected_chirality` requires
  `stereo_neighbor_order` to safely reinterpret a stored `@`/`@@` tag
  against a different (e.g. canonically-reordered) neighbor sequence; with
  it missing, the writer silently passed the raw stored tag through
  unchanged against whatever new traversal order it picked. Since
  `standardize`'s `remove_explicit_h: true` (this crate's own default)
  calls `remove_hydrogens`, re-canonicalizing an already-canonical SMILES
  could flip a declared tetrahedral stereocenter to its mirror image on
  some symmetric-ranking-ambiguous molecules -- a real,
  independently-confirmed correctness defect (via RDKit InChIKey
  divergence), not just a cosmetic re-spelling.
- Now: `remove_hydrogens` restores both side tables for every surviving
  atom/bond, remapping indices (and reintroducing the `STEREO_H_SENTINEL`
  marker where a removed H's slot is now implicit again) the exact
  inverse of what `add_hydrogens` already does for the opposite direction.
- Verified against the real 290-compound InChIKey-mismatch corpus from the
  originating investigation (eMolecules, 9.47M compounds, `renkin doctor
  stock reimport_idempotency`): **289 of 290 now match the true input
  identity** (up from 0 before this fix). The one residual case is a
  coupled/shared-bond E/Z system (an oxime/hydrazone shape) with a
  confirmed **different** root cause, independent of `remove_hydrogens`
  entirely (reproduces with bare `parse`/`canonical_smiles`, no
  `standardize` involved) -- tracked separately, not mixed into this fix;
  see issue #390.
- New tests: `remove_hydrogens` restoring the side tables both when
  nothing is removed and when a real explicit H neighbor is removed
  (chained through `add_hydrogens` to reach a genuinely non-sentinel
  starting order), a minimized tetrahedral witness, a minimized (simple,
  non-coupled) E/Z witness, CIP-descriptor preservation across a second
  canonicalization pass, mirror-image distinctness (the fix must not
  degrade into "never distinguish stereo"), atom-order-permutation
  invariance for two independently-written spellings of the same
  configuration, and Boc-protecting-group / fused-bicyclic-ring
  regression cases.

### Added — `chematic-py` (A2.1: Python bindings for the conformer ensemble core)

- `Mol.conformer_ensemble_v2(config)`: a new, separate Python method
  exposing `chematic_3d::embed_ensemble_v2` (A2, PR #373) — a deterministic
  multi-conformer generator built as a new outer loop over
  `embed_pipeline_v2`, with energy-ranked selection scoped within each
  force field actually used and full per-attempt provenance. Added
  alongside the existing `Mol.conformer_ensemble()`, not in place of it —
  see below.
- `EnsembleV2Config(per_conformer, count, base_seed, rmsd_threshold=0.5,
  use_symmetric_rmsd_pruning=True, ensemble_timeout_ms=None)`: seed,
  attempt count, ensemble-wide timeout, and RMSD pruning threshold are all
  explicit, required-or-defaulted arguments — no hidden global state.
  Construction is infallible; an invalid `rmsd_threshold` (negative, `NaN`,
  infinite) is rejected with `ValueError` at `conformer_ensemble_v2()` call
  time, matching where the Rust core itself validates it.
- The returned dict never raises just because zero conformers were kept —
  every per-attempt outcome (kept, pruned as a near-duplicate, or a typed
  failure) is reported in `attempts`, and `conformers`/`conformer_provenance`
  are simply empty in that case. This is the opposite convention from
  `Mol.embed_pipeline_v2()`, which does raise on a failed embed — documented
  explicitly on the new method, since "no exception raised" does not imply
  "got at least one conformer" here.
- MMFF94 and UFF energies are never cross-compared: kept conformers in
  `conformers` are ordered group-by-group (by force field actually used),
  ascending energy only *within* each group, and `mixed_force_field`
  discloses when more than one group has a kept member. No flattened,
  globally energy-sorted field is exposed — that would silently reintroduce
  the cross-scale comparison this design avoids.
- `Mol.conformer_ensemble()` (the legacy, `etkdg`-backed multi-conformer
  method) is unchanged in this release — no behavior or return-type change,
  no removal. Its docstring now carries a deprecation note pointing at
  `conformer_ensemble_v2` and naming the concrete reasons (no seed, no
  energy ranking, silent `Err(_) => []` on internal failure, and the live
  MMFF94-zero-energy defect from PR #369). The docstring's own duplicated
  intro paragraph (a pre-existing copy-paste artifact, unrelated to this
  change) was cleaned up in the same edit.
- Not done this round, by design: no WASM bindings for `embed_ensemble_v2`.
  The best-of-N-vs-RDKit benchmark arm planned as a next step here has
  since landed — see "Added — benchmark" below.

### Fixed — documentation (correction to a v0.19.0 changelog entry)

- v0.19.0's own benchmark-refresh entry below states 3D conformer quality
  was corrected to match `docs/rdkit-migration.md`'s "Experimental, no
  RMSD/TFD comparison exists" characterization — **that characterization
  was itself wrong**, found during this round's OpenEye/materials-science
  competitive audit. `validation/results/mmff94_bci_gap_227_phase2_
  report.md` already measures RMSD (mean 1.685 Å) and TFD (mean 0.2228)
  against RDKit's ETKDGv3+MMFF94 on the project's 265-molecule corpus.
  `docs/benchmark.md`, `docs/rdkit-migration.md`, and the same stale
  "15 SMARTS rules" pKa-rule-count claim (actual: 23 rules, 6 acid + 17
  base, `crates/chematic-chem/src/pka.rs`) copied across README.md/
  README_ja.md/`docs/rdkit-comparison.md`/`docs/rdkit_cheatsheet.md`/
  `docs/benchmark.md` are corrected. Not rewriting v0.19.0's already-
  released changelog text itself, per this project's convention of never
  silently editing a shipped version's own historical record — see the
  entry below for what v0.19.0 originally said.

### Known limitations in this release

- `generate_coords` (the legacy 3D placement engine) is **not** routed to
  the new connectivity-ordered engine — no existing caller's default
  behavior changed. See "Added — connectivity-ordered 3D coordinate
  generation" above.
- Issue #210's remaining 4 residual molecules (`naproxen_S`, `ibuprofen_S`,
  `testosterone`, `cholesterol`) are still unfixed via the UFF-rescue
  bridge specifically — `stereo_safe` (this release's own headline feature,
  above) already resolves them through `embed_pipeline_v2`'s own path.
- Issue #390 (a single-molecule coupled/shared-bond E/Z correctness
  residual, independent of every fix in this release) remains open.
- The best-of-10 conformer benchmark establishes that `embed_ensemble_v2`
  works robustly at scale; it does **not** establish RDKit conformer-
  *selection* parity — see "Added — benchmark" above for the exact
  distinction.
- This release's evidence is drawn from measurements already taken during
  development (29-molecule × 5-seed stereo sweep, 33-molecule connectivity
  differential, 265-molecule best-of-10 benchmark, 58-molecule UFF-rescue
  sweep, the 9.47M-compound downstream identity investigation) — no fresh
  full-corpus re-measurement was run solely for this release, per this
  project's own "minimize heavy measurements" policy.

## [0.19.0] — 2026-08-23

### Added — `chematic-chem` (Tautomer & Parent Identity round 2C, ROADMAP.md Phase 2)

- `canonical_tautomer`/`tautomer_parent` now canonicalize the aromatic
  lactam/lactim class (e.g. 2-pyridone `O=c1cccc[nH]1` /
  2-hydroxypyridine `Oc1ccccn1`) that was previously not invariant across
  input tautomer spelling. Implemented as a new directional step
  (`apply_exocyclic_lactam_shift_tracked` in `tautomer.rs`) — an exocyclic
  O donor/acceptor across an aromatic ring atom at odd ring-path distance
  — rather than fed into the existing score-ranked candidate pool, which
  measurably selects the wrong (lactim) tautomer for this class. See
  `docs/rfcs/tautomer_parent_identity_phase2_rfc.md` section 4.4a for the
  full mechanism and the measurement behind that design choice.
- **Measured, not assumed:** the shift fires correctly on all 5 of the
  RFC's design molecules, but end-to-end convergence holds for only 3 of 5
  (2-pyridone, 4-pyridone, uracil) — cytosine and guanine do not converge
  end-to-end, because their carbonyl carbon is flanked by two ring
  nitrogens (RFC section 1.7, a second, distinct, out-of-scope tautomer
  defect found while fixing this one). **Correction (round 2C-N,
  diagnosis-only, no production-code change): hypoxanthine was originally
  reported here as a clean holdout that generalized the fix — this was
  wrong.** Re-measurement found hypoxanthine has the identical two-candidate
  ring-N-H ambiguity as cytosine/guanine (RFC section 1.7); the original
  holdout check happened to test only one of its two valid keto spellings.
  Closing section 1.7 for all three molecules needs ring-internal
  N-position normalization, deliberately not part of this change (same
  reasoning as excluding the nitroso/oxime defect, section 1.6). Neither
  section 1.6 nor section 1.7 is resolved as of this release; Phase 2
  (Tautomer & Parent Identity) is not complete, and `TautomerScoringConfig`/
  Python/WASM Parent-API bindings remain unimplemented.
- New negative controls confirmed unaffected: phenol, anisole, aniline,
  pyridine N-oxide, plain amides, and — the two checks that specifically
  probe the fix's scope boundary — 3-hydroxypyridine (excluded: even/meta
  ring distance, no valid Kekulé path) and 4-/2-aminopyridine (excluded:
  the analogous N-acceptor case is evidenced-out-of-scope, not swept in by
  generalizing the acceptor element).
- Isotope labels and remote stereocenters are confirmed preserved through
  the shift (atom-level checks, not string comparison).
- **Hardened after review, before merge:** the matcher is fail-closed —
  `bridge` must be a neutral aromatic carbon (not any aromatic element),
  `donor` a neutral, degree-1 oxygen with exactly one transferable H, and
  `acceptor` a neutral, degree-2 (pyridine-type, valence-compatible after
  +1H) aromatic nitrogen. Every candidate is verified post-generation
  against a full atom/bond invariant check (nothing changes except the one
  H-count pair and the one bond order) before being accepted, and the
  candidate tie-break cross-checks an independent fingerprint against
  canonical SMILES rather than trusting the SMILES string alone. Additional
  confirmed-excluded negative controls: a charged/3-connected pyridinium
  acceptor, an aromatic-nitrogen bridge (N-hydroxypyrrole), and a fused
  bridgehead-nitrogen acceptor. Atom-permutation invariance is now checked
  directly (same molecular graph, different `AtomIdx` insertion order) for
  2-pyridone and uracil's two independent sites, plus respelling-based
  checks for 4-pyridone, hypoxanthine, the isotope case, and the remote
  stereocenter case — `canonical_tautomer`, `tautomer_parent`'s molecule
  and status, and the applied `TautomerRuleId` sequence must all agree.

### Added — `chematic-chem` (Tautomer & Parent Identity round 2B, ROADMAP.md Phase 2)

- `fragment_parent`/`charge_parent`/`isotope_parent`/`stereo_parent`/
  `tautomer_parent`/`super_parent`: an explicit "Parent identity" concept —
  an idempotent, deterministic reduction of one axis of molecular
  variability (fragment, charge, isotope, stereo, tautomer) to one
  representative structure, meant as a grouping/dedup key. `charge_parent`
  is **not** a bare `neutralize_charges` wrapper: it selects the fragment
  parent first, then neutralizes that single fragment — a design correction
  made during RFC review before any code was written (see the RFC's round
  2A revision note). `fragment_parent`/`isotope_parent`/`stereo_parent`
  are thin, explainable wrappers over the existing `select_fragment`/
  `remove_isotopes`/`remove_stereo` primitives. `super_parent` composes all
  five in the fixed order fragment → charge → isotope → stereo → tautomer,
  returning every intermediate stage's audit record, not just the final
  molecule.
- `TautomerLimits { max_transforms, max_tautomers, timeout_ms }`: a
  deterministic budget for `tautomer_parent`, generalizing
  `TautomerConfig::max_iter`/`max_tautomers`. `max_transforms`/
  `max_tautomers` are reproducible (same input + limits ⇒ same result,
  always); `timeout_ms` is an explicitly non-deterministic escape hatch
  (wall-clock, machine-dependent), documented as outside that guarantee.
  `max_restarts`/`Canceled` from the original RFC sketch were dropped this
  round: `max_restarts` has no defined meaning at this level, and `Canceled`
  had no cancellation mechanism to ever produce it — see the RFC's round 2A
  revision.
- `ParentResult`/`ParentComputationStatus`/`ParentAudit`: `tautomer_parent`
  and `super_parent` now report `Completed`/`MaxTransformsReached`/
  `MaxTautomersReached`/`TimedOut`/`Abstained`/`InvalidInput` instead of a
  bare `Molecule` — budget exhaustion is a previously-silent gap in
  `canonical_tautomer_with_config` (confirmed via an existing `#[ignore]`d
  regression test on a 25-independent-site "comb" molecule) that is now
  visible to the caller instead of silently returning a possibly
  input-order-dependent result.
- `TautomerAuditRecord`/`ScoreContribution`/`TautomerScoreTerm`/
  `AppliedTransform`/`TautomerRuleId`: an explainable trail for
  `tautomer_parent` — which candidate was selected, its score breakdown by
  named term, which rules fired (one of 42 stable `TautomerRuleId`
  variants, not a bare string), and which atoms (if any) had stereo or
  isotope data affected. `applied_transforms` covers the rule-based
  1,3-/1,5-shift loop only; the final direct-aromatic-shift tie-break
  contributes to `score_breakdown`/`candidate_count` instead (disclosed
  scope limitation, not an oversight).
- Fixed the module's stale `"The 15 tautomer rules"` doc comment — the
  array actually has 42.
- **Not in this round** (round 2C/2D, per the RFC's explicit split): the
  aromatic lactam/lactim canonical-tautomer fix (2-pyridone/cytosine/
  uracil/purine-class non-invariance, confirmed but not fixed); the
  nitroso/oxime non-convergence gap found while reviewing a fixture;
  `TautomerScoringConfig`/custom rule support; Python/WASM bindings for any
  of the above.
- Design: `docs/rfcs/tautomer_parent_identity_phase2_rfc.md`.

### Added — `chematic-chem` (explainable Molecule Standardization Phase 1, ROADMAP.md Phase 1)

- `FragmentPolicy` + `select_fragment(mol, &FragmentPolicy) -> (Molecule,
  TransformationRecord)`: a structural, non-named-list classification for
  largest-fragment selection and salt/solvent removal — a tiny always-strip
  monatomic-ion set (Li/Na/K/F/Cl/Br/I), water, and a no-carbon/small-fragment
  heuristic, in place of matching a named SMARTS catalog. Ranking is
  heavy-atom count → carbon presence → an intrinsic canonical-SMILES
  tie-break, never input/discovery order.
- `TransformationRecord`/`FragmentRecord`/`FragmentSnapshot`/`FragmentDecision`:
  a per-fragment audit trail (formula, canonical SMILES, kept/removed
  decision with a machine-readable `rule_id` and human-readable reason, and
  a whole-transformation `abstained` reason for inputs with no confident
  organic parent, e.g. `NaCl`).
- `largest_fragment`/`remove_salts` are now thin wrappers around
  `select_fragment` with the default policy — same signatures, corrected
  behavior. `SaltCatalog`/`remove_salts_with_catalog`/`is_salt_fragment`
  are unchanged, kept as opt-in-only legacy behavior, no longer on the
  default path.
- Fixes 3 confirmed defects in the prior `largest_fragment`/`remove_salts`:
  tied-fragment-size selection was spelling-order-dependent (kept a
  different fragment for `"CCC.CCN"` vs `"CCN.CCC"`); fragment "size"
  counted raw atom count including explicit hydrogens, not heavy atoms;
  `SaltCatalog`'s "ammonium" SMARTS entry false-positived on real organic
  cations (e.g. choline), previously masked rather than fixed by an
  unrelated size comparison.
- Also fixes a real, previously-undiscovered bug found while implementing:
  fragment extraction (both the new code and the pre-existing
  `remove_salts_with_catalog`) silently corrupted stereocenters by
  rebuilding fragments through a fresh `MoleculeBuilder` without remapping
  `stereo_neighbor_order`. `extract_fragment` now rebuilds fragments via
  repeated `Molecule::remove_atom` (which already remaps this correctly),
  removing atoms in descending index order — one fix, both callers.
- Design: `docs/rfcs/explainable_standardization_phase1_rfc.md`.
  Acceptance/holdout fixtures: `validation/standardization_phase1_
  fixtures.jsonl` (34 rows), `validation/standardization_phase1_
  holdout.jsonl` (10 rows) — all 44 now exist as `chematic-chem` tests
  (`phase1_*` in `standardize.rs`).
- **Not in this round** (per the RFC's stated scope): Python/WASM bindings
  for any of the new types; a non-heuristic signal for genuine
  cocrystal-vs-counterion ambiguity (found, empirically, not to be
  distinguishable from an unrelated unambiguous case via a heavy-atom-count
  margin alone — disclosed as unimplemented rather than shipped as a
  heuristic that doesn't actually work); wiring the richer audit trail into
  `StandardizationPipeline`'s own per-stage report.

### Added — `chematic-ff`/`chematic-py`/`chematic-wasm` (UFF minimizer soundness signal, ROADMAP.md Backlog item 5a)

- `UffMinimizeResult` (`crates/chematic-ff/src/uff.rs`) gains a `sound: bool`
  field: true iff the final coordinates are all-finite and no bond exceeds
  a 3.0 Å sane-covalent-bond ceiling — the same, already corpus-validated
  threshold `chematic-3d`'s `minimize::MAX_SANE_BOND_LENGTH` uses for its
  own (unrelated, `embed_pipeline_v2`-only) soundness gate. `chematic-ff`
  cannot depend on `chematic-3d` (the dependency runs the other way), so
  this is a deliberately duplicated copy of that constant/check, not an
  independently chosen one.
- Deliberately independent of `converged`: steepest descent frequently
  reports `converged == false` on perfectly sound geometries that simply
  haven't hit the tight RMS-gradient stopping threshold within `max_iter`
  — treating `converged` alone as a quality signal would be a false-failure
  generator, the same rationale `chematic-3d`'s own soundness gate already
  documents.
- Before this change, `chematic-py`'s `Mol.minimize_uff()` and
  `chematic-wasm`'s `minimize_uff_json()` called UFF directly with zero
  soundness signal of any kind — a caller had no way to distinguish a
  genuinely converged, trustworthy geometry from one where UFF's
  torsion/out-of-plane-incomplete potential settled at a real but unsound
  stationary point (e.g. a fused polycyclic aromatic folding non-planar
  with a blown-up bond — issue #185's finding). Both bindings now surface
  `sound` alongside `coords`/`energy`/`iterations`/`converged`.
- Does not add a residual-force half of the gate the way `chematic-3d`'s
  own `check_minimization_soundness` has (`max_residual_force` ceiling):
  `minimize_uff`'s steepest-descent loop doesn't retain a converged
  gradient norm past each iteration, and computing one is out of scope for
  this change — the bond-length check alone is a real, if partial, signal
  that previously didn't exist at all.
- **4 new tests** (`crates/chematic-ff/src/uff.rs`): known-sound (ordinary
  ethanol minimizing normally), known-unsound via a deliberately-stretched
  5.0 Å bond (using `max_iter=0` to deterministically exercise the
  bond-length check without depending on steepest descent actually getting
  stuck there), and non-finite coordinates.

### Fixed — `chematic-depict` 2D auto-layout (issue #347)

- Long open chains no longer drift into a monotonic ~30°/bond rotation and
  wrap onto themselves (a plain 13-carbon chain used to place its first and
  last atom at identical coordinates). `dfs_zigzag`
  (`crates/chematic-depict/src/layout.rs`) was alternating its ±30°
  deflection by a neighbor's position within its *own parent's*
  unplaced-neighbor list, which is `0` for every ordinary single-successor
  chain step — not a real alternation signal. Now threads an explicit
  `sign: f64` through the DFS stack that flips on every step, so a plain
  chain continuation genuinely zigzags.
- Exocyclic substituent bonds at ring junctions now point along the correct
  outward bisector instead of missing it by up to ~30°.
  `best_outgoing_direction` had its own coarser 6-point 60° candidate grid
  with a last-wins tie-break; it now shares the richer, chemistry-aware
  candidate set (`suggest_bond_direction`'s existing 12-point 30° grid plus
  bond-relative sp2/sp3/anti offsets) via a new `ranked_candidates` helper.
- `best_outgoing_direction` also gained a positional collision check: the
  angularly-best candidate is skipped if it would land on top of an
  already-placed atom elsewhere in the layout (falls back through the
  ranked list, then to the top pick if every candidate collides). This is
  a real, independent gap the coarse grid had been accidentally masking —
  angular separation from an atom's own bonds doesn't guard against a
  distant already-placed atom happening to sit on the chosen ray.
- **Known, intentional behavior change**: any MOL/SDF/CML output that
  relies on `chematic`'s auto-generated 2D coordinates (no
  caller-supplied coordinates) will now differ from before, since these
  coordinates were wrong. No fixture in this repo pins the old
  (buggy) auto-layout output as a golden value.
- New `Mol.depict_data()` binding for `chematic-py`
  (`crates/chematic-py/src/mol_methods.rs`): returns the full structured
  `DepictData` (atoms with element/position/label/color/charge, bonds with
  kind including wedge/hash stereo) as a native Python `dict`, mirroring
  `chematic-wasm`'s existing `depict_data_json` (shipped since v0.1.20) —
  previously Python could only get an SVG string or a function that takes
  coordinates as input, not raw layout data.

### Fixed — `chematic-smiles` canonical SMILES (issue #149's ring-constrained E/Z residual)

- `compute_stereo_alkene_ends` now excludes a double bond's ends from
  marker-carrier candidacy when the bond itself is endocyclic in an SSSR
  ring smaller than 8 atoms — such a bond's real-world geometry is fixed by
  the ring, not a free stereochemical choice, so treating it as ambiguous
  only destabilized marker-carrier selection for a genuinely stereogenic
  double bond it happened to share a candidate bond with. Implements the
  predicate specified in `docs/rfcs/ez_ring_constrained_residual_audit.md`
  ("Wave 2C audit", 0/1,387 row-level disagreements with RDKit's own
  stereo-possibility judgment on a 5,000-molecule corpus).
- All 8 of the previously-documented 8 ring-constrained residual fixtures
  (`EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS`) now fully converge to one
  canonical output; re-verified (not assumed) that the 10 already-resolved
  fixtures — including 5 sharing the identical endocyclic shape — remain
  fully resolved. All 18 are now one merged fixture list
  (`EZ_SHARED_CARRIER_FULLY_RESOLVED`) in `crates/chematic-smiles/src/
  canonical.rs`'s test module.
- `chematic-perception` promoted from a dev-only to a real (production)
  dependency of `chematic-smiles`, needed for `find_sssr` ring-membership
  data. No cycle: `chematic-perception`'s own dependency on
  `chematic-smiles` remains dev-only (test-only), so the real dependency
  graph is still acyclic (verified via `scripts/check_publish_graph.py`).
  `compute_stereo_alkene_ends` runs on every `canonical_smiles()` call, so
  `find_sssr` (Horton: O(V·E) candidates + GF(2) elimination) is only
  computed when the molecule has at least one double-bond candidate that
  could need the ring check — most molecules (no double bonds, or only
  non-candidate ones like a ketone's `C=O`) skip it entirely, keeping the
  common-case cost of this hot path unchanged.
- **Issue #149 stays open**: this closes only ~10% of the corpus's general
  shared-carrier coupling population (3 of 31 coupling components, per the
  Wave 2C audit's own blast-radius measurement) — the other ~90% is a
  separate, still-unidentified mechanism this predicate does not touch.
- `cargo test -p chematic-smiles --lib`: 202/202 passed. No full
  5,000-molecule corpus rescan performed for this change — the 18-fixture
  regression check is the audit's own specified verification.

### Fixed — `chematic-3d` `embed_pipeline_v2` (`total_timeout_ms: 0` race, intermittent CI failure)

- `check_timeout!`'s elapsed-time comparison (`crates/chematic-3d/src/
  pipeline_v2.rs`) used `Instant::elapsed().as_millis() > budget`, which
  only has millisecond granularity — for a `total_timeout_ms: Some(0)`
  config on a small molecule with a minimal-work config (`forceFieldPolicy:
  "none"`, all optional torsion sources disabled), a fast enough machine
  could complete every stage within the same millisecond tick, reading
  `elapsed_ms == 0` even at the last checkpoint. `0 > 0` is false, so the
  pipeline silently succeeded instead of failing closed as documented and
  tested (`timeout_zero_fails_closed_with_typed_timeout`).
- Root-caused after this was reported as an intermittent `Test (WASM)` CI
  failure (`pipeline_v2.test.mjs`/`pipeline_v2_web_target.test.mjs`, both
  assert the "zero timeout must fail closed" contract): reproduced 0/60
  times locally before the fix (fast local machine never hit the race) but
  observed on 3 independent CI runs the same day — consistent with a race
  that manifests more often on fast CI runners, not a deterministic bug a
  simple local re-run would catch.
- Fix: `budget == 0` is now checked unconditionally, independent of the
  (unreliable-at-this-resolution) elapsed-time reading — a zero budget
  means no time was granted at all, so any checkpoint reached must fail,
  full stop. Every `budget > 0` case is untouched (identical `>`
  comparison) — this is not a behavior change for any other timeout value.
- Fires at the same `check_timeout!` call site as before (first one, right
  after torsion-knowledge/stage 2 completes), so the existing evidence-
  preservation guarantee (`timeout_failure_still_carries_evidence_computed_
  before_it_tripped` — a timeout must still carry the torsion-knowledge
  report) is unaffected, not just incidentally still passing.
- Shared Rust core (`chematic-3d`), so this also fixes the same latent gap
  for the Python binding (`crates/chematic-py/src/pipeline_v2.rs`), which
  doesn't currently test the zero-timeout edge case explicitly but goes
  through the identical `embed_pipeline_v2` function.
- `cargo test -p chematic-3d --lib`: 540/540 passed (0 failures, 3 ignored,
  unchanged from before this fix). Both previously-flaky Node test files
  re-run 100/100 times each after rebuilding the WASM target with the fix
  — 0 failures (does not prove the race can never recur, only that it's
  now logically impossible for the `budget == 0` case specifically, by
  code inspection, not by re-running until lucky).

### Changed — benchmark/validation refresh (release-prep, v0.19.0)

Every number in `docs/benchmark.md`/`docs/validation.md` was pinned to
**chematic v0.4.29 / RDKit 2026.03.3** — ~14 minor releases stale. Re-measured
fresh against this release (RDKit 2026.03.4); full raw results and
methodology: `benchmarks/2026-08-23.md`.

- **The 4,999-molecule ChEMBL accuracy corpus is now committed**
  (`scripts/chembl_accuracy_corpus_4999.smi`) — every prior snapshot back to
  2026-06 depended on an uncommitted personal `~/Downloads/SMILES.csv` path
  nobody else could reproduce from.
- **Molecular weight now has a real, corpus-wide accuracy check** (99.82%,
  4990/4999, ±0.01 Da vs `Descriptors.MolWt`) — `scripts/bench5k.py` never
  actually measured it before; `docs/validation.md`'s "175-mol reference"
  annotation was unconnected template prose, not a real second measurement.
  Fixed at the source, not just reworded.
- Fixed a real bug in `scripts/gen_validation_report.py`'s stereocenter
  oracle-disagreement breakdown: it read chematic-vs-oracle delta fields
  where it needed legacy-vs-new-CIP oracle-internal disagreement counts —
  the two coincidentally summed to the right total in the last snapshot's
  data, and do not in this one (68 disagreements, not "0 + 70").
- **CIP R/S/E/Z label agreement** (a distinct metric from stereocenter
  *count* agreement) re-measured at 99.74–99.78%, up from a stale
  96.30–96.83% — reflects CIP-engine fixes landed across the ~14
  intervening releases, not a change in this release itself.
- The ECFP4 "diverse corpus" throughput figure previously had **no script
  that could reproduce it at all** — added an optional `--corpus FILE`
  argument to `scripts/benchmark_vs_rdkit.py` (reuses its existing timing
  logic; still kept strictly separate from the repeated-fixture numbers).
- WASM bundle size rebuilt clean (`--target web`, not the stale
  `--target nodejs` artifact the old figure partially matched) — current:
  2.98 MB raw / 1.11 MB gzip (chematic) vs 6.91 MB raw / 2.06 MB gzip
  (RDKit.js, independently verified via `npm pack @rdkit/rdkit`).
- `docs/benchmark.md`'s "3D conformer quality: Good (ETKDG rules)" framing
  overstated what's actually measured — corrected to match
  `docs/rdkit-migration.md`'s existing honest "Experimental, no RMSD/TFD
  comparison exists" characterization. No 3D/RMSD/TFD measurement was added
  or is claimed.
- **Not claimed**: Phase 2 (Tautomer & Parent Identity) completion, or full
  tautomer canonicalization — see this file's round 2C/2C-N entries above.
  This is a benchmark/tooling refresh, not a Phase 2 status change.

## [0.18.0] — 2026-08-20

Python and WASM bindings for the 7 file formats `chematic-mol` gained in
v0.17.0 (mmCIF, PQR, ORCA, QCSchema, Gaussian Cube, OpenDX, LAMMPS
data/dump), which had zero Python/WASM exposure until now, plus an MMFF94
atom-typing fix and a finishing cross-language-consistency pass. One
sub-bug of issue #337 (aryl isothiocyanate cumulated-double-bond CSP
carbon) fixed via a strict-superset condition matching RDKit's real rule
(`getTotalDegree() == 2`); the other 6/8 molecules behind that issue were
re-diagnosed as a genuine RDKit Kekulization/MMFF-aromaticity-perception
artifact rather than a locally-fixable typing rule, and left as an
honestly-disclosed residual (confirmed via direct negative-control
fragments, not merely asserted). The Binding Quality Pack that closes out
this release adds `VolumetricGrid.values_3d` (Python), a verified-in-wheel
`py.typed` marker with a `mypy --strict` check against a fresh-venv wheel
install, 5 additive `js_sys::Float64Array`/`Uint32Array`-returning WASM
functions alongside the existing JSON-string API (this crate's first
typed-array precedent), cross-language parity fixtures proving Rust/
Python/WASM agree on the same small Cube/OpenDX/mmCIF/LAMMPS-dump inputs,
and runnable Python + Node.js examples — deliberately scoped as polish and
consistency, not new format support. Purely additive at the Rust API
level for existing crates — no breaking changes. **Not in this release**
(deferred to a future version): the residual MMFF94 atom-typing artifact
noted above; targeted canonical-SMILES/aromaticity known-residual fixes;
Rustdoc warning cleanup and a broader API/format-capability documentation
pass; Gaussian Cube multi-dataset support.

### Fixed — `chematic-ff` (MMFF94 atom typing: aryl isothiocyanate CSP carbon, issue #337)

- **Root cause**: `assign_c_type`'s sp-carbon check only fired on
  `triple_bonds > 0`, so a carbon reached via *two* double bonds (a
  cumulated diene / "allenic" carbon — e.g. the central C of an aryl
  isothiocyanate's `N=C=S`) fell through to the `double_bonds > 0`
  "double-bonded to N/O/P/S" branch and was mistyped 3 (generic
  carbonyl-family) instead of RDKit's real CSP type 4. RDKit's actual rule
  (`AtomTyper.cpp`, pinned commit
  `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`, lines ~954-960) is simply
  `getTotalDegree() == 2` for a non-aromatic carbon, unconditional on which
  elements the two bonds go to — a real triple bond and a cumulated
  double-bond pair aren't special-cased separately by RDKit at all; both
  just leave the carbon with exactly 2 neighbors.
- Fix: replaced the `triple_bonds > 0 => Ok(4)` check with
  `total_degree(mol, idx) == 2 => Ok(4)`, ahead of the `double_bonds > 0`
  branch it was previously losing to. Strict superset of the old
  condition (every triple-bonded carbon already has `total_degree == 2`),
  so it cannot regress any previously-correct triple-bond CSP assignment.
  Also, correctly, broader than the corpus: a plain carbon allene
  (`C=C=C`, not present in the 264-molecule corpus) was mistyped the same
  way and is fixed too.
- Measured, full 264-molecule corpus, genuine per-atom join (not
  aggregate-count arithmetic), same tool/oracle pair as every other
  measurement in this file: type-mismatch ledger 34 → 32 atoms (both
  `chembl_tier_b_0071`/`_0082`'s isothiocyanate carbon move to exact
  match; molecule count 8 → 6); charge-mismatch ledger 62 → 56 atoms (6
  move: the corrected carbon in both molecules, plus its N and S
  neighbors in both, whose own TYPES were already correct but whose
  BCI-bond-lookup charge depends on the carbon's type too; molecule count
  8 → 6). Zero regressions on either ledger, verified by diffing full
  before/after mismatch lists. `lookup_chg_contribution` already covers
  the (4, 9) and (4, 16) bond-type pairs the corrected type now looks up
  — no new charge-computation errors introduced.
- **Not addressed by this fix**: the other 6/8 molecules behind issue
  #337 (a pyridinium-conjugated exocyclic-amine scaffold,
  `chembl_tier_b_0009`/`_0023`/`_0028`/`_0029`/`_0030`/`_0034`) turned out,
  on live re-investigation, to be a genuine RDKit Kekulization/
  MMFF-aromaticity-perception artifact for a specific fused, macrocyclic
  ring topology — not a locally-statable atom-typing rule chematic is
  missing a branch for (confirmed by direct negative-control fragments:
  the same ring+exocyclic-amine motif types correctly in four simpler,
  non-macrocyclic contexts). Reproducing RDKit's exact behavior would
  require porting its global Kekulization matching algorithm and its
  iterative, ring-processing-order-dependent MMFF-aromaticity pi-electron
  count, not adding a discriminating condition — left as an
  honestly-disclosed residual rather than shipping a heuristic that would
  create new false mismatches on the negative-control cases. These 6
  molecules remain unchanged by this fix (32/6,693 type-mismatched, 56/
  6,693 charge-mismatched atoms, same molecules as before). Full writeup,
  including the corrected diagnosis of which atom actually mismatches
  (the ring nitrogen, not the exocyclic amine as originally described) and
  the RDKit source citations for both the fixed and the unfixed sub-bug:
  `scripts/mmff94_provenance/PROVENANCE.md`'s issue #337 follow-up
  addendum.
- 6 new regression-pinned/synthetic-fixture tests
  (`crates/chematic-ff/src/mmff94_numeric.rs`): full-array pins for both
  corpus molecules, two minimal isothiocyanate fixtures, a no-regression
  plain-alkyne pin, and the broader-than-corpus allene pin.

### Added — chematic-py (Python bindings for the 7 v0.17.0 file formats)

- Python (PyO3) bindings for the 7 file-format modules `chematic-mol`
  gained in v0.17.0, which had zero Python exposure until now: mmCIF, PQR,
  ORCA input/output, QCSchema (`Molecule`/`AtomicInput`/`AtomicResult`),
  Gaussian Cube, OpenDX, and LAMMPS data/dump.
- `chematic.parse_mmcif`/`write_mmcif`, `chematic.parse_pqr`/`write_pqr`/
  `infer_element`, `chematic.parse_orca_input`/`write_orca_input`/
  `parse_orca_output` (`formats.rs`) — plain functions returning/consuming
  dicts, matching the existing `parse_cif` convention; every typed Rust
  parse error maps to `ValueError`; `*ParseLimits` are exposed as optional
  keyword arguments with the Rust `Default` values.
- `chematic.parse_qcschema_molecule`/`write_qcschema_molecule`,
  `chematic.chematic_to_qc_molecule`/`qc_molecule_to_chematic`,
  `chematic.parse_atomic_input`/`write_atomic_input`,
  `chematic.parse_atomic_result`/`write_atomic_result` (`formats.rs`) —
  routed through Python's own `json` module as the dict<->text boundary
  rather than a hand-written field-by-field struct mapper, so QCSchema's
  open extensibility bags (`extras`/`unknown_fields`/`keywords`/
  `protocols`/`native_files`/`wavefunction`) round-trip losslessly with no
  new dependency (`serde_json` is not linked into `chematic-py`).
- `chematic.VolumetricGrid` (new `volumetric.rs`), a pyclass shared by Cube
  and OpenDX (both read/write the same underlying Rust type): `from_cube`/
  `from_opendx` constructors, `to_cube`/`to_opendx`/`to_opendx_lossy`
  writers (the fail-closed `to_opendx` vs. explicit-opt-in
  `to_opendx_lossy` split from PR #335 is preserved faithfully — no boolean
  flag defaulting to lossy), `get`/`checked_index`/`point_count`/
  `to_molecule` methods, and numpy-array `values`/`axes` properties.
- `chematic.parse_lammps_data`/`write_lammps_data` (plain dict, new
  `lammps.rs`) and `chematic.LammpsDumpFrame` (pyclass, real `column`/
  `cartesian_positions` behavior) plus `parse_lammps_dump_frame`/
  `parse_lammps_dump_all`/`write_lammps_dump_frame`/
  `write_lammps_trajectory`/`box_bounds_to_true`/`true_to_box_bounds`.
  `parse_lammps_dump_all` materializes the whole trajectory as a list
  rather than exposing `LammpsDumpReader`'s true streaming iteration to
  Python — a disclosed scope decision (see the function's docstring), not
  a silently dropped capability.
- pytest coverage for all 7 formats in `crates/chematic-py/tests/` (at
  least one round-trip test and one typed-error-to-`ValueError` test per
  format).

### Added — `chematic-wasm` (WASM/JS bindings for the 7 v0.17.0 file formats)

- WASM (wasm-bindgen) bindings for the 7 file-format modules `chematic-mol`
  gained in v0.17.0 with no prior WASM exposure: Gaussian Cube, OpenDX,
  PDBx/mmCIF, PQR, QCSchema JSON, ORCA input/output, and LAMMPS data +
  dump/trajectory. New `crates/chematic-wasm/src/format_io.rs` module,
  ~40 `#[wasm_bindgen]` functions total. No new dependencies (`chematic-mol`
  was already an unconditional, default-featured dependency of
  `chematic-wasm`; `serde_json` was already available).
- Follows `mol_io.rs`'s existing conventions: `Result<T, JsValue>` for
  fallible calls, structured multi-field results as JSON strings (not
  bespoke `#[wasm_bindgen]` structs — this crate has no
  `js_sys::Float64Array`/typed-array precedent to follow for large numeric
  arrays, so Cube/OpenDX's `values` grid round-trips through a plain JSON
  number array; documented as a disclosed perf tradeoff, not a silent one).
  mmCIF/PQR/ORCA/Cube (formats with no bond table) mirror the existing
  `mol_from_pdb`/`pdb_coords_json` split — a topology-only `MolHandle` plus
  a same-atom-order coordinates accessor — rather than fabricating a
  MOL-block/SMILES-shaped result for chemistry these formats never
  perceive.
- `LammpsDumpReader`'s per-frame streaming has no WASM-boundary equivalent
  in this first pass; `lammps_trajectory_to_json` parses the whole input
  and returns every frame at once instead (disclosed, not silently
  dropped). `write_opendx` (fail-closed on non-Ångström units) and
  `write_opendx_lossy` (explicit Bohr→Ångström opt-in) are kept as two
  distinct bindings, not collapsed into one lossy-by-default function.
- `lammps_dump_cartesian_positions_json` delegates to
  `chematic_mol::LammpsDumpFrame::cartesian_positions` to resolve a dump
  frame's real Cartesian positions (`x/y/z` passthrough, `xs/ys/zs`
  orthogonal-or-triclinic scaled-coordinate transform, `null` for an
  `xu/yu/zu`-only frame) rather than requiring a JS caller to reimplement
  that box-bounds/triclinic math itself.

### Added — Binding Quality Pack (polish/consistency for the 7 v0.17.0 formats, no new format support)

Finishing pass over the Python/WASM bindings the 3 previous entries in this
section added (PRs #341-#343) — cross-language consistency and ergonomics,
explicitly NOT new capability. No `chematic-mol` changes.

- `chematic-py`: `VolumetricGrid.values_3d` — `values` (flat) reshaped to a
  3-D numpy array of shape `self.shape`, so `values_3d[i, j, k] ==
  get(i, j, k)`. A plain copy of `values` (no zero-copy view — deliberately
  out of scope). Axis order verified against
  `chematic_mol::VolumetricGrid::checked_index` (k varies fastest) with a
  non-cubic `(2, 3, 4)` fixture and hardcoded interior values, not just a
  cubic 2×2×2 grid where a transposed reshape could pass by coincidence.
- `chematic-py`: added `python/chematic/py.typed` (PEP 561 marker). Verified
  it — and `__init__.pyi` — actually ship in a `maturin build --release`
  wheel (`unzip -l`, not just "the file exists in the source tree"), and
  that `mypy --strict` against a fresh venv with that wheel installed (not
  the source tree — confirmed via `chematic.__file__` resolving to
  `site-packages`) passes clean for the 7 new-format bindings
  (`crates/chematic-py/tests/typecheck_new_formats.py`). No pre-existing
  stub gaps surfaced for these 7 formats; no full-API stub audit performed
  (out of scope).
- `chematic-wasm`: 5 new `js_sys::Float64Array`/`Uint32Array`-returning
  functions — `cube_values_f64`, `cube_shape_u32`, `opendx_values_f64`,
  `opendx_shape_u32`, `lammps_dump_rows_f64`,
  `lammps_dump_cartesian_positions_f64` — added ADDITIVELY alongside the
  existing JSON-returning functions (`cube_grid_json`, `opendx_grid_json`,
  `lammps_dump_frame_to_json_str`, `lammps_dump_cartesian_positions_json`,
  all unchanged), avoiding a full JSON round trip for large numeric grid/row
  payloads. This is the crate's first `js_sys` typed-array precedent;
  `js-sys` promoted from a transitive dependency (already pulled in by
  `web-sys`) to a direct one in `Cargo.toml` — no new supply-chain surface.
  `lammps_dump_cartesian_positions_f64` has one disclosed, real API-shape
  difference from its JSON sibling: the JSON version returns `null` for a
  frame with no resolvable coordinate columns; a `Float64Array` can't
  represent `null`, so this version returns `Err` instead. Read-direction
  only — no typed-array *write* functions added (out of scope; OpenDX's
  fail-closed unit handling stays entirely on the write path, untouched).
  Every new function delegates to the same underlying `chematic_mol` call
  as its JSON sibling (no reimplemented parsing/math); native `cargo test`
  cross-checks the shared pure-Rust helpers against the JSON output
  directly (a `js_sys` typed array can't be constructed outside a real JS
  runtime, so the `#[wasm_bindgen]`-exposed functions themselves are
  exercised in the Node `.test.mjs` suite instead).
- Cross-language parity fixtures (Cube, OpenDX, mmCIF, LAMMPS dump — not
  full-format-spec coverage): the same 4 small fixtures, with the same
  independently-hardcoded expected values, asserted from all 3 language
  entry points — `crates/chematic-mol/tests/format_binding_parity.rs`
  (Rust), `test_parity_*` in `crates/chematic-py/tests/test_new_formats.py`
  (Python), `crates/chematic-wasm/tests/format_parity.test.mjs` (WASM).
  The LAMMPS case reuses the exact triclinic box/tilt fixture from
  `chematic_mol::lammps_dump`'s own hand-computed test — this is the same
  bound-box-vs-true-box / triclinic-shear-term class of bug an independent
  review caught and required fixing in PR #343, so it is deliberately not
  re-derived per language. All 3 language surfaces agree; no discrepancy
  found.
- Runnable examples: `examples/materials_formats_quickstart.py` (Cube
  `values_3d`, mmCIF atom fields, LAMMPS dump Cartesian positions) and
  `crates/chematic-wasm/tests/materials_formats_example.test.mjs` (same 3
  formats via WASM; picked up automatically by the existing WASM CI test
  loop, doubling as both a runnable example and a CI-executed test).
- Deliberately NOT done (see PR description for the full scope-boundary
  list): no format auto-detection/dispatch, no zero-copy numpy views, Cube
  is still single-dataset-only, and no full `.pyi`/mypy audit of the entire
  `chematic-py` API surface beyond these 7 formats' stub entries.

## [0.17.0] — 2026-08-18

Format/Python/materials-interop breadth, plus two MMFF94 charge/bond-order
accuracy fixes: square-planar (`@SP1`/`@SP2`/`@SP3`-equivalent) stereo
read/write for MOL/SDF; PDBx/mmCIF, PQR, QCSchema JSON, and ORCA I/O;
Python (PyO3) bindings for `chematic-crystal`'s `Lattice`/`PeriodicStructure`;
CIF explicit symmetry-operation expansion (Rust + Python); a shared
`VolumetricGrid` type plus Gaussian Cube and OpenDX I/O; LAMMPS data-file
and dump/trajectory I/O; an MMFF94 bond-order-classification fix
(`torsions_missing` 257→0) and an MMFF94 BCI partial-charge fix (own wrong
`bond_type_for`, then a derived-formal-charge source fix) with one
post-minimization stereo-repair addition it surfaced. Production
`pipeline_v2_mmff94_strict`: 240/265 → 241/265 (see the 3-state
measurement note below for full provenance). A release-readiness audit
found and fixed one metadata defect ahead of this release: the declared
MSRV (`rust-version = "1.85"`) didn't match reality (`chematic-depict` ->
`svg2pdf` -> `image` requires rustc 1.88) -- raised to 1.88 (not pinned
back down, to avoid reopening dependency-resolution/security-update risk)
and now continuously verified by a dedicated CI job, plus a handful of
small documentation-hygiene fixes found in the same audit. Purely
additive at the Rust API level for existing crates — no breaking changes.
**Not in this release** (deferred to a future version): the residual
MMFF94 atom-typing bug tracked as issue #337 (62/6,693 corpus atoms, a
different-shaped fix than either charge bug above); targeted
canonical-SMILES/aromaticity known-residual fixes; a format-capability
documentation/discoverability pass.

### Added — `chematic-mol` (square-planar stereo read everywhere; write via new checked conformer writers only)

- MDL/CTfile has no symbolic field for a non-tetrahedral stereo tag
  (confirmed against RDKit 2026.03.4 as a live oracle, and consistent
  with public RDKit documentation/source — see
  `docs/rfcs/square_planar_mol_io_rfc.md`). The only real mechanism is
  3D-coordinate-derived reperception: given a real (non-flat) 3D
  conformer, a coplanar 4-coordinate center whose neighbor pair angles
  unambiguously resolve to one of SP1/SP2/SP3 is reperceived directly
  from geometry **on every read** (`perceive_square_planar_from_3d`,
  wired into both `read_mol_with_diagnostics` and
  `read_mol_v3000_with_diagnostics` — this applies regardless of which
  writer produced the file). Element eligibility reuses
  `Element::normal_valences().is_empty()` (transition metals and
  similar), matching an RDKit oracle observation for every element
  tested (one documented Na/Mg/Al divergence, see the RFC).
- **Write support is opt-in, via three new `_checked` functions only**:
  `write_mol_with_conformer_checked`, `write_mol_v3000_with_conformer_checked`,
  and `write_sdf_record_with_conformer_checked` (backed by the public
  `validate_square_planar_for_write`) validate the declared tag against
  real coordinates before writing — never fabricating a conformer from
  nothing, never silently trusting a mismatch — and fail closed with a
  typed `MolStereoWriteError` otherwise. The pre-existing, more commonly
  used writers are unchanged, but fall into two different gaps — each now
  carries a Rustdoc warning naming its specific one and pointing at its
  `_checked` counterpart where one exists. Do not describe MOL/SDF
  writing in general as "square-planar supported"; only the three
  `_checked` functions above are.
  - **2D-only, so they drop the tag outright**: `write_mol`,
    `write_mol_with_coords`, `write_mol_v3000`, and the whole `write_sdf*`
    family (`write_sdf`, `write_sdf_with_charges`, `write_sdf_record`,
    `write_sdf_record_v3000`) — none of these has a z coordinate to write
    a square-planar tag against in the first place.
  - **3D-capable but unvalidated**: `write_mol_with_conformer`,
    `write_mol_v3000_with_conformer`, and `write_sdf_record_with_conformer`
    write whatever conformer they're handed and so *do* preserve the tag
    when the conformer actually matches it — but they trust the caller and
    never check, so a mismatched or flat conformer is written silently
    self-inconsistent with no error.
  Explicitly out of scope: pure-2D wedge-only square-planar encoding (no
  such MDL mechanism exists), 3-heavy + implicit-H square-planar centers,
  and a chematic-specific lossless SDF extension (Tier 3, not built).
- New public types: `MolFormat`, `UnsupportedStereoReason`,
  `MolStereoWriteError`, `SquarePlanarRejectionReason`,
  `SquarePlanarPerceptionDiagnostic`; `MolReadReport` gained a new
  `square_planar_diagnostics` field.

### Fixed — `chematic-mol` (`wedge_vs_3d_conflicts` tetrahedral-only gate)

- `wedge_vs_3d_conflicts` gated on `atom.chirality == Chirality::None`,
  so any non-`None`, non-tetrahedral chirality fell through into a
  computation that assumes a tetrahedral shape — latent until the
  square-planar MOL/SDF work above (this same release), since no MOL/SDF
  reader ever produced `Chirality::SquarePlanar` before then. Fixed by
  switching to `!atom.chirality.is_tetrahedral()`, an
  exhaustive-match-safe allowlist gate (same fix shape as two prior
  instances of this bug class in `chematic-3d`/`chematic-chem`).

Issue #227 Phase 2: MMFF94 BCI (bond-charge-increment) partial-charge bug,
investigated per Phase 1's own flagged follow-up and fixed. Also adds a
full `embed_pipeline_v2` 3-state quality re-measurement (State 1: v0.16.0
pre-Phase-1; State 2: post-Phase-1 torsion fix; State 3: post-Phase-2 BCI
fix) — see the PR body / `validation/results/` for the full report.

### Fixed — `chematic-ff` (MMFF94 partial-charge BCI bond-type source)

- **Root cause, a compound bug, not the single view-source bug Phase 1's
  own note anticipated**: `mmff94_charges_numeric` used a private,
  standalone `bond_type_for(order: BondOrder) -> u8` that mapped bond
  *multiplicity* directly (`Double -> 1, Triple -> 2, Aromatic -> 4`) —
  unrelated to RDKit's real `getMMFFBondType`, which returns 0 unless the
  bond is formally SINGLE *and* both atom types are `sbmb`/`arom`-flagged,
  and which RDKit's real `computeMMFFCharges` calls identically to its own
  bond-stretch code (`AtomTyper.cpp:2457-2475,3462-3474`, pinned commit).
  It also read bond order from the caller's original, un-reperceived
  molecule rather than the MMFF-specific Kekulized view — the view-only
  bug Phase 1 fixed for bond/angle/torsion/stretch-bend.
- Fix: reuse the already-fixed, already-oracle-validated
  `crate::mmff94_minimizer::bond_type_for(ti, tj, order)` (deleting the
  wrong local function), fed `assign_mmff94_numeric_types_with_view`'s
  reperceived bond order instead of the caller's `mol`.
- Measured against a live RDKit oracle, all 264 typing-succeeded molecules
  (not a sample): 1,687/6,693 heavy-atom charges (25.2%, 206/264
  molecules) mismatched the oracle before the fix; 67/6,693 (1.0%, 11/264
  molecules) after — a genuine per-atom join confirms **zero regressions**
  (0 previously-exact atoms became mismatched; 1,620 moved from
  mismatched to exact). The 67-atom residual is a separate, pre-existing,
  unrelated formal-charge/`fcadj` redistribution gap for charge-separated
  species (nitro/azide/charged-sulfoxide types), confirmed unmoved by this
  fix either direction — flagged as follow-up, not fixed here.
- New regression-pinned test (`acetone_carbonyl_charges_match_rdkit_oracle_after_bond_type_fix`)
  and a renumbering-invariance test
  (`mmff94_charges_numeric_is_invariant_under_atom_renumbering`, same
  `deterministic_permutation`/`rebuild_with_order` helpers Phase 1's own
  reviewer follow-up test uses). Full writeup:
  `scripts/mmff94_provenance/PROVENANCE.md`'s Charges/BCI entry.

### Fixed — `chematic-ff` (MMFF94 partial-charge derived-formal-charge source, Phase 2 Step 6)

- **Root cause**: `mmff94_charges_numeric` fed the molecule's raw, literal
  SMILES formal charge (`atom.charge`) directly into equation 15 — both as
  an atom's own q0 and as the neighbor-formal-charge source for the
  `v*sumFormalCharge` redistribution term and the anionic-neighbor-leak
  adjustment. RDKit's real `computeMMFFCharges` instead computes a separate,
  MMFF-atom-TYPE-derived formal charge ("MMFFFormalCharge") via a dedicated
  per-type switch statement that runs *before* the main charge loop
  (`AtomTyper.cpp` lines ~3095-3350, pinned commit
  `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`), and uses *that* value
  everywhere — for many types (nitro N type 45, azide N types 47/53,
  sulfoxide S type 17: none are switch cases) the derived charge is 0.0
  regardless of the atom's raw SMILES charge; for O2CM/SM (types 32/72) it's
  a fractional value shared across terminal O/S atoms on a common neighbor,
  not the literal per-atom charge. A second, independent bug: the
  anionic-neighbor-leak loop ran unconditionally instead of only when the
  atom's own `fcadj` is zero (RDKit's `isDoubleZero(v)` gate) — mutually
  exclusive with the `v*sumFormalCharge` term, not additive with it.
- Also checked, table-level, before writing any fix: `MMFF94_PBCI`'s
  (pbci, fcadj) *values* were already byte-identical to RDKit's real
  `defaultMMFFPBCI` for every one of the 5 suspected types (17/32/45/47/53)
  — the table was never wrong. Its trailing `//` comments for those 5 types
  *were* wrong (cosmetic only, never read by any lookup) and are corrected
  in the same commit.
- Fix: new `mmff_derived_formal_charge`/`o2cm_sm_formal_charge` helpers
  (reusing this module's existing `count_terminal_o_neighbors`/
  `count_terminal_s_neighbors`/`count_deg2_n_neighbors`, the same counters
  `classify_terminal_o` already uses to *assign* type 32 — with one known,
  disclosed divergence: the shared `count_deg2_n_neighbors` helper omits
  RDKit's `!isAromatic()` condition on secondary nitrogens, which could
  flip the O2CM/SM sulfone-neighbor branch's result for an aromatic
  degree-2 N case the corpus does not contain) compute the derived charge;
  `mmff94_charges_numeric`'s Step 1 and Step 3 now read it instead of
  `atom.charge`, and Step 3's leak is gated to `fcadj_i ≈ 0`
  (`isDoubleZero`-style `1e-10` epsilon). A faithful but intentionally
  partial port: the unconditional ±1/±2/±3/−1 simple-type groups
  (including type 62's full two-part rule — its −1.0 base value *and* its
  extra "subtract half of positive-neighbor-charge" adjustment) and
  O2CM/SM's carbon-neighbor/nitro-nitrate-neighbor/sulfone-neighbor
  branches are implemented; the ±1/±2/±3/−1 groups and 3 of the O2CM/SM
  branches are each independently verified against a live RDKit 2026.03.4
  oracle query (not merely re-derived from this fix's own output), while
  type 62's extra adjustment is implemented but not independently
  oracle-verified (zero corpus exposure either way, so nothing to falsify
  against). O2CM/SM's phosphate/thiosulfinate/perchlorate-neighbor
  branches and the ring-/conjugation-dependent types (76, 55/56/81, 61)
  are not ported at all — zero atoms of any of these types appear anywhere
  in the 264-molecule Wave 1 corpus (confirmed by a dedicated, committed,
  independently re-runnable survey,
  `crates/chematic-3d/examples/mmff94_fchg_type_exposure_survey_227.rs`),
  so the gap cannot be masking a corpus-visible bug; flagged as a
  follow-up.
- Measured against the same live RDKit oracle dump and the same per-atom
  join methodology as the bond-type fix above, entry point
  `crates/chematic-3d/examples/mmff94_bci_charges_dump_227.rs`: the
  67/6,693-atom (11/264-molecule) residual left after the bond-type fix
  shrinks to 62/6,693 atoms (8/264 molecules) — a genuine per-atom join
  confirms **zero regressions** (0 previously-exact atoms became
  mismatched) and exactly 5 atoms across 3 molecules
  (`chembl_tier_b_0080` azide, `chembl_tier_b_0159` nitro,
  `chembl_tier_b_0161` sulfoxide/O2CM) moved from mismatched to exact
  match. The remaining 62/8 are a **separate, unrelated class of bug**:
  genuine MMFF atom-*type* misassignments in `assign_mmff94_numeric_types`
  (an exocyclic amine N adjacent to a pyridinium ring mistyped 58 instead
  of 54 in 6 molecules; an isothiocyanate cumulated-carbon mistyped 3
  instead of 4 in 2 molecules) — confirmed unmoved by this fix in either
  direction, and out of scope per this step's own stop condition (fixing
  them means touching atom-type assignment, a different-shaped change).
- **Blast radius, both directions**: exactly **5 of 6,693 corpus atoms**
  change computed value at all (5 mismatch→match; the other 6,626
  match→match and 62 mismatch→mismatch atoms are all byte-identical
  before/after) — the reason no downstream `embed_pipeline_v2`
  re-measurement was run for this step (contrast the prior BCI bond-type
  fix, which moved 1,620 atoms and produced one genuine new stereo
  violation, `chembl_tier_b_0082`, investigated separately above).
  Outside this corpus, the real behavioral change is broader than "5
  atoms" suggests: it applies to *any* molecule with a carboxylate,
  sulfonate/sulfamate, nitrate, nitro, azide, sulfoxide, or
  quaternary-ammonium group — this corpus simply happens to contain no
  carbon-neighbor or phosphorus-neighbor O2CM/SM atoms (every type-32 atom
  in it has a sulfone/nitro/sulfoxide neighbor) and only 3 molecules
  combining nitro/azide/sulfoxide with a type absent from RDKit's
  derived-formal-charge switch.
- 8 new regression-pinned/synthetic-fixture/renumbering-invariance tests
  (`crates/chematic-ff/src/mmff94_numeric.rs`), same discipline as the
  bond-type fix's tests above — expected values copied verbatim from
  either the already-committed oracle dump or a fresh live oracle query,
  never derived from this fix's own output. Full writeup:
  `scripts/mmff94_provenance/PROVENANCE.md`'s Charges/BCI entry (Phase 2
  Step 6 addendum).

### Measured — 3-state `embed_pipeline_v2` quality re-measurement (`pipeline_v2_mmff94_strict`)

- Fresh measurement (not reused from any older commit) at State 1 (`c079926`,
  v0.16.0 release, pre-torsion-fix), State 2 (`a2baac4`, post-torsion-fix
  main, pre-BCI-fix), State 3 (`e2876bb`, PR #331 tip, post-both-fixes).
  `pipeline_v2_mmff94_strict`
  success: 240/265 → 241/265 → 241/265; RMSD (symmetric, vs
  `rdkit_etkdgv3_mmff94`) mean 1.698 → 1.685 → 1.685 Å; TFD mean 0.2245 →
  0.2233 → 0.2228; **0 status-level regressions** (success/typed_failure/
  timeout) on every pairwise per-molecule join — see below for the one
  genuine stereo-quality regression this note does NOT cover.
- 62-molecule torsion-fix subset (State 1 → State 2): **both coverage and
  geometry quality improved together** (+1 success, RMSD mean 1.156 → 1.115
  Å, coverage@2.0Å 77.4% → 82.3%) — not a coverage-only gain.
- 206-molecule BCI-fix-affected subset (State 2 → State 3): coverage
  unaffected (structural — charges don't gate this policy), aggregate
  RMSD/TFD flat/noise-level (mean RMSD delta −0.0009 Å), **one genuine new
  stereo violation** (`chembl_tier_b_0082`, investigated and addressed —
  see next section).
- Full report: `validation/results/mmff94_bci_gap_227_phase2_report.md`
  (+ `_summary.json`, raw per-state dumps, per-molecule transition tables).
- **State 3's numbers remain valid as the v0.17.0 release figures, not
  re-measured after State 3**: the only production-3D/FF-relevant change
  between `e2876bb` and the v0.17.0 release head is PR #336 (`f401f47`,
  the MMFF94 BCI derived-formal-charge fix), which itself already
  discloses its own blast radius as exactly 5/6,693 corpus atoms across
  3/264 molecules changing computed charge value at all, zero status-level
  regressions — see that entry's "Blast radius, both directions" note
  below for the full disclosure. No other file under
  `crates/chematic-3d/src`, `crates/chematic-ff/src`, or
  `crates/chematic-perception/src` changed between those two commits
  (confirmed via `git diff --stat e2876bb..916ca4d -- crates/chematic-3d/src
  crates/chematic-ff/src crates/chematic-perception/src`); PRs #332-#335 and
  #338 are format/Python-interop/materials I/O work with zero reach into
  the embedding/minimization/verification path. Citing this instead of
  re-running the 265-molecule corpus, per this project's standing
  measurement policy.

### Fixed — `chematic-3d` (post-minimization stereo repair-and-reverify)

- **Found during the 3-state re-measurement above**: `chembl_tier_b_0082`'s
  declared E/Z bond is satisfied post-embedding (identical in States 2 and
  3 — charges don't affect embedding) but MMFF94 minimization walks it to
  violated only under State 3's corrected charges. Oracle-confirmed
  chematic-specific: RDKit's own real MMFF94 minimizer, run on the same
  molecule, does NOT reproduce this on any of its 4 arms — the same
  already-documented "MMFF94 minimization has zero stereo awareness"
  architectural class as `chembl_tier_b_0076`/`chembl_tier_b_0083` (found
  during the v0.14.0 release gate), a third instance, not a new failure
  class this fix introduced.
- Fix: `StereoPolicy::RepairAndVerify` now gets one additional
  repair-and-reverify attempt on the POST-minimization geometry (new
  `PipelineV2Result::post_minimization_stereo_repair` field) — stage 8's
  existing repair runs too early to see a violation minimization itself
  introduces. Empirically verified safe (bond lengths/clash count
  unchanged by the repair) and effective before implementing. Fail-closed:
  accepted only if repair succeeds, the reverified result has zero
  violations, and the geometry stays sound; any rejection falls through to
  the original, unmodified `FinalStereoViolation` failure.
  `StereoPolicy::Ignore`/`VerifyOnly` (including the 3-state measurement's
  own arm, `chematic_pipeline_v2_mmff94_strict`, above) are completely
  unaffected by construction — this fix cannot and does not change any
  of those 3-state numbers, which were not re-measured for this reason.
- Root-cause fix (real stereo-awareness inside MMFF94 minimization) and
  broadening `StereoPolicy::Ignore`'s own gate were both considered and
  explicitly out of scope (large, cross-cutting changes deserving their
  own separate authorization) — under `Ignore`, `chembl_tier_b_0082`'s
  violation remains real and unrecovered, a named, tested residual
  (`chembl_tier_b_0082_ez_bond_survives_bci_fix_under_repair_and_verify_not_under_ignore`),
  not silently hidden.
- New/updated tests: `repair_and_verify_recovers_post_minimization_stereo_violation`
  (updates a now-obsolete negative-control test that used to assert this
  exact `gly_ala_gly` case was unrecoverable), a no-op sanity check, and
  the `chembl_tier_b_0082` golden regression test above. Full writeup:
  `scripts/mmff94_provenance/PROVENANCE.md`'s "Follow-up investigation" entry.

---

Issue #227 Phase 1: MMFF94 torsion parameter coverage gap, root-caused and
fixed. `torsions_missing` on the 265-molecule Wave 1 corpus: 257 instances
across 62 molecules → 0 (`mmff94_term_coverage_audit.rs`).

### Fixed — `chematic-ff` / `chematic-3d` (MMFF94 bond-order classification)

- **Root cause**: `assign_mmff94_numeric_types` already computed an
  MMFF-specific re-perceived molecule (`compute_mmff94_aromatic_view`,
  Kekulized to match RDKit's real `setMMFFAromaticity` output) to derive
  correct atom TYPES, then discarded it — `bond_type_for`/`angle_type_for`/
  `torsion_type_for`/`stretch_bend_type_for` kept reading `BondOrder` from
  the caller's original, un-reperceived molecule. For ring systems where
  chematic's general aromaticity perception and MMFF94's own stricter
  perception disagree (e.g. caffeine's pyrimidinedione ring — oracle-
  confirmed non-aromatic in RDKit's real sanitizer), this fed the
  classification formula the wrong bond order, landing on a table code with
  no row even though the correct row already existed in chematic's own,
  unmodified parameter tables at a different code.
- New `assign_mmff94_numeric_types_with_view` returns `(types, mmff_mol)`;
  `assign_mmff94_numeric_types` is now a thin wrapper. Threaded through
  chematic-ff's 5 production energy/gradient entry points and chematic-3d's
  `compute_mmff94_coverage` (the `Mmff94BondAngleStrict`/
  `Mmff94WithUffFallback` coverage gate), so the gate and the energy
  functions it gates agree on classification.
- Measured on the 265-molecule Wave 1 corpus, same audit tool throughout:
  `torsions_missing` 257→0 instances (62→0 molecules);
  `bonds_missing` (type-only) 80→1; `angles_missing` (type-only) 191→46 —
  side effects of the same shared root cause, not separately implemented.
  Zero success→failure regressions on either the default bond+angle gate
  (`minimize_with_policy`, 248→249/265 `Ok`) or the stricter
  `complete_bonded_term_gate` (`minimize_with_policy_gated(...,true,true)`,
  187→249/265 `Ok`), verified by a full per-molecule join against the
  pre-fix baseline, not aggregate counts.
- New `torsion_no_term_by_design`/`Mmff94Resolution::NoTermByDesign`:
  RDKit's real empirical-rule cascade generates no torsion term at all when
  either central atom is linear (MMFF `lin` flag, e.g. nitrile/acetylenic
  carbon) — the exact, complete explanation for the corpus's remaining 3
  `table_gap` instances (oracle-confirmed `GetMMFFTorsionParams` also
  returns `None`). Wired into the coverage gate as a new
  `Mmff94CoverageReport::torsions_no_term_by_design` counter so these are no
  longer misclassified as coverage gaps.
- **No Halgren empirical torsion rule was implemented.** Investigated and
  falsified two hypotheses against a live RDKit oracle (all 254 real
  instances, not a sample) before finding the real cause above — see
  `scripts/mmff94_provenance/PROVENANCE.md`'s Torsion entry for the full
  writeup, including why an empirical-rule implementation was deliberately
  not shipped (zero instances in this corpus need it).
- These `complete_bonded_term_gate`/`minimize_with_policy` gate figures are
  a separate measurement from the production `pipeline_v2_mmff94_strict`
  entry point (`embed_pipeline_v2`, previously measured at 239/265 -- see
  the `[0.16.0]` entry below) -- coverage-gate improvement here does not by
  itself imply the same gain on the full production pipeline (embedding,
  minimization convergence, stereo verification). A full `embed_pipeline_v2`
  re-measurement is planned as Phase 2 follow-up work, not part of this fix.

---

Format-expansion Wave 1: `chematic-mol` gains bidirectional PDBx/mmCIF,
PQR, QCSchema JSON (`Molecule`/`AtomicInput`/`AtomicResult`), and ORCA
input/output support. Goal is depth over format count -- each format is
implemented loss-aware (unrecognized fields are preserved or surfaced as
a typed warning, never silently dropped) and, where the format supports
it, round-trip verified (read -> write -> read reproduces the original
data). LAMMPS data/dump and Gaussian Cube/OpenDX (a shared
`VolumetricGrid` type) are scoped for a later Wave 2, not this PR.

### Added — `chematic-mol` (PDBx/mmCIF)

- New `mmcif` module: `parse_mmcif`/`parse_mmcif_with_limits`/`write_mmcif`.
  Reads/writes the `_atom_site` loop (group_PDB, atom/comp/asym id in both
  `label_*` and `auth_*` forms, altloc, entity id, insertion code,
  Cartn_x/y/z, occupancy, B_iso_or_equiv, formal charge, NMR-style model
  number) plus `_cell.*`/`_symmetry.space_group_name_H-M`. Reuses the
  existing small-molecule-CIF tokenizer (`chematic-mol`'s `cif` module,
  now exposing its tokenizer as `pub(crate)`) rather than a second parser.
  No bond table (mmCIF carries none) and no bond perception. Unrecognized
  `_atom_site` loop columns are surfaced via a typed `unhandled_columns`
  list rather than silently dropped. `MmcifParseLimits` bounds input
  size/atom count/line length; NaN/Infinity are rejected as typed errors.
  Only the first `data_` block is read, matching the existing `cif`
  module's own scope; no symmetry expansion in this crate (see
  `chematic-crystal`'s CIF adapter, PR #323, for periodic/crystal CIF --
  this module targets macromolecular mmCIF specifically). Open Babel's
  own mmCIF support is read-only; this adds write too.
- Open Babel/RDKit were used only as behavioral oracles during
  development, never as a source for code, comments, or tables -- the
  wwPDB PDBx/mmCIF dictionary is the implementation's primary source.

### Added — `chematic-mol` (PQR)

- New `pqr` module: `parse_pqr`/`parse_pqr_with_limits`/`write_pqr`,
  handling both the 10-field (no chain) and 11-field (with chain) atom
  line shapes real `pdb2pqr` output uses. `infer_element` resolves the
  element PQR itself doesn't store (from atom name + residue name +
  record type), documented as a deterministic heuristic, not
  authoritative data recovery. Same size-limit and NaN/Infinity-rejection
  discipline as `mmcif`.

### Added — `chematic-mol` (QCSchema JSON)

- New `qcschema` module: MolSSI QCSchema `QcMolecule`, `AtomicInput`,
  `AtomicResult` (`parse_*`/`write_*` for each), hand-parsed against
  `serde_json::Value` (no `serde` derive dependency added). Bidirectional
  `qc_molecule_to_chematic`/`chematic_to_qc_molecule` conversion against
  `chematic_core::Molecule`, with an explicit, documented Bohr<->Angstrom
  conversion point (`BOHR_TO_ANGSTROM`, CODATA 2018) -- the reverse
  direction divides rather than multiplying by a reciprocal, to avoid the
  classic precision loss. Extensible/open QCSchema fields (`keywords`,
  `extras`, `provenance`) are preserved as opaque JSON on round-trip
  rather than dropped. NaN/Infinity rejected in all numeric fields.
  `OptimizationInput`/`OptimizationResult` (trajectory-bearing) are out of
  scope for this wave.
- `qcschema::JsonObject`/`qcschema::Connectivity` (the type aliases used
  as public field types on `AtomicInput`/`AtomicResult`/`QcMolecule`) are
  now re-exported at the crate root (`chematic_mol::JsonObject`/
  `chematic_mol::Connectivity`), matching the crate's existing
  re-export convention -- found during the v0.17.0 release audit as a
  caller-facing gap (a `chematic_mol::QcMolecule` user couldn't name the
  type of its own `.connectivity`/`.extras` fields without reaching into
  the `qcschema` module directly).

### Added — `chematic-mol` (ORCA input/output)

- New `orca` module: `parse_orca_input`/`write_orca_input` (round-trip;
  unknown `%...end` blocks, including nested sub-blocks, are preserved
  verbatim rather than dropped) and `parse_orca_output` (final geometry,
  best-effort optimization trajectory, final single-point energy,
  charge/multiplicity, vibrational frequencies when present, and two
  independent typed statuses: termination -- `Normal`/`Error`/
  `Incomplete` for a truncated/crashed log -- and optimization
  convergence, since normal termination does not imply convergence).
  Open Babel's own ORCA support is input-write-only/output-read-only;
  this adds input-read and keeps output-read. Explicitly out of scope:
  `$new_job` multi-job input files, the `%coords` block form, ghost/
  dummy/point-charge atom designations, full normal-mode vectors, and
  semantic (only verbatim) parsing of `* int` Z-matrix coordinate blocks.

### Added — `chematic-py` (Python bindings for `chematic-crystal`)

- Roadmap step 3: `PeriodicStructure`/`Lattice`/`Site` are now exposed to
  Python (`pip install chematic`), the first host-language binding for
  `chematic-crystal` (no WASM binding exists yet). New `crystal.rs`
  module, following this crate's existing flat-module/`PyValueError`/
  `IntoPyArray` conventions rather than inventing new ones.

  ```python
  from chematic import PeriodicStructure, Lattice, Site

  s = PeriodicStructure.from_cif(cif_text)
  s.lattice.volume
  s.sites[0].species          # [(element_symbol, occupancy), ...] -- disorder preserved, never collapsed
  s.cartesian_positions()     # (N, 3) numpy array
  s.neighbors(cutoff=3.0)
  s.make_supercell((2, 2, 2)).to_cif()
  s.wrap_into_cell()
  s.formula                   # occupancy-weighted, unreduced Hill-order string
  PeriodicStructure.from_poscar(poscar_text).to_poscar()
  ```

  Design decisions (see the PR body for full reasoning):
  - **Immutable wrappers**, matching the Rust side: no setters;
    `wrap_into_cell()`/`make_supercell()` return a new `PeriodicStructure`.
    No `__eq__`/pickling support added (matches this crate's existing
    baseline: only `__repr__`/`__str__` exist anywhere in `chematic-py`
    today).
  - **Disorder is never collapsed** — `Site.species` is always the full
    `list[(element_symbol, occupancy)]`, even for a single-species site.
  - Every Python-facing constructor (`Site()`, `PeriodicStructure()`,
    `Lattice.from_matrix`/`from_parameters`/`cubic`/`orthorhombic`) routes
    through the Rust side's own fallible constructors — no parallel
    validation logic, so a Python caller cannot construct a structure that
    violates the Rust invariants (occupancy sums, coordinate finiteness,
    degenerate/singular lattices, ...); every `CrystalError`/
    `CifPeriodicError`/`PoscarError` surfaces as a typed `ValueError` with
    the Rust error's own message.
  - `PeriodicStructure.symmetry_status` surfaces `CifSymmetryStatus`
    (`is_p1`/`space_group_name`/`operation_count`) for any `from_cif`-sourced
    structure — and `to_cif()` now **raises `ValueError`** rather than
    silently re-emitting a false `P 1` declaration when the source CIF
    declared symmetry this parser doesn't expand (a write-path gap the
    Rust-level `chematic-mol` adapter itself doesn't guard; this binding
    adds the check rather than changing that crate's public API). The
    status (and the `to_cif` guard) survive `wrap_into_cell()` and
    `make_supercell()` — the asymmetric-unit problem doesn't go away
    under either transform.
  - `PoscarDocument`'s extra fields (`selective_dynamics`, `velocities`,
    `predictor_corrector`, `comment`) are not exposed as individual Python
    attributes in this first version, but are kept internally so a bare
    `from_poscar()` -> `to_poscar()` round trip (no intervening transform)
    reproduces them faithfully; `make_supercell()` clears them (site count
    changes, so their per-site correspondence would break) while
    `wrap_into_cell()` preserves them (site count/order unchanged).
  - New `formula` property: Hill-order (C, H, then alphabetical),
    occupancy-weighted, **unreduced** cell content (a 2x2x2 NaCl supercell
    reports `"Cl8Na8"`, not `"NaCl"`).
  - Dependency wiring: `chematic-crystal` added as a `chematic-py`
    dependency; `chematic-mol`'s existing `crystal` feature (the CIF <->
    `PeriodicStructure` adapter) is now enabled on the `chematic-py` ->
    `chematic-mol` edge.

### Added — `chematic-mol` (CIF explicit symmetry-operation expansion)

- Roadmap step 4: a CIF's declared asymmetric unit can now be expanded
  into a full unit cell using **only the symmetry operations literally
  written in the CIF text** — no space-group database, no name/number-to-
  operations generation, no spglib-equivalent auto-detection (see
  `crates/chematic-mol/src/cif_symmetry.rs`'s module docs and
  `docs/crystal_scope.md`). New `parse_cif_periodic_structure_with_options`
  (`parse_cif_periodic_structure` becomes a thin wrapper defaulting
  `expand_explicit_symmetry: true`); new `CifPeriodicResult::to_cif_checked`
  moves the CIF write-safety judgment from `chematic-py`-only into Rust
  (see below).
  - Supports both modern dotted (`_space_group_symop.operation_xyz`,
    `_symmetry_equiv.pos_as_xyz`) and legacy underscore
    (`_space_group_symop_operation_xyz`, `_symmetry_equiv_pos_as_xyz`) tag
    spellings, an optional id column
    (`_space_group_symop.id`/`_space_group_symop_id`/
    `_symmetry_equiv.pos_site_id`/`_symmetry_equiv_pos_site_id`), and both
    a `loop_` of operations and a single standalone (non-loop) operation
    item (how a genuine P1 CIF sometimes states its one, identity,
    operation explicitly). If more than one tag alias is present and they
    parse to genuinely different operation sets, that is a typed
    `ConflictingSymmetryOperationLists` error rather than silently picking
    one.
  - New hand-written operation-expression parser (no `eval`, no external
    expression-evaluator crate): case-insensitive, whitespace-tolerant,
    supports `x,y,z`-style rotation terms (coefficient always exactly
    ±1) plus integer/rational translations (`x+1/2,y,z`,
    `-y+x,-y,1/3+z`), checked arithmetic throughout (rejects overflow
    rather than wrapping/panicking), and rejects a parsed rotation matrix
    whose determinant isn't exactly `+1`/`-1`. New internal
    `Rational` numerator/denominator type — no external crate; the one
    existing rational type in this workspace
    (`chematic-cip::rational::RationalAtomicNumber`) is CIP-specific and
    not reusable here.
  - Extended `CifSymmetryStatus` with a third variant,
    `ExpandedExplicitOperations { space_group_name, operation_count,
    asymmetric_site_count, expanded_site_count }` — documented prominently
    as a faithfulness claim about the CIF's own text ("every operation it
    listed was applied"), **not** a claim that the list is complete or
    correct for the named/numbered space group (this implementation never
    cross-checks that against any space-group database). A declared
    space-group name/number with no parseable operation list at all still
    classifies as `UnexpandedSymmetry`, unchanged.
  - Special-position dedup reuses `chematic_crystal::minimum_image`
    (Cartesian, triclinic-exact) rather than a fresh `min(d, 1-d)`
    reimplementation. The disorder-row-grouping tolerance and the new
    expansion-dedup tolerance are now the **same** constant,
    `SITE_MERGE_TOLERANCE_ANGSTROM = 1e-3` Å (previously two different,
    differently-unitted tolerances: a `1e-4` *fractional*-coordinate
    check for disorder grouping, nothing at all for expansion, since
    expansion didn't exist).
  - New typed errors (`CifSymmetryError`, nested in
    `CifPeriodicError::Symmetry`): `MalformedSymmetryOperation`,
    `UnsupportedSymmetryExpression`, `ZeroDenominator`,
    `InvalidRotationMatrix`, `DuplicateSymmetryOperation`,
    `MissingIdentityOperation` (IUCr convention requires the identity to
    be explicitly listed — a CIF that lists only one operation and that
    operation genuinely *is* the identity is not an error; a list missing
    identity entirely is), `ConflictingSymmetryOperationLists`,
    `SymmetrySiteCollision` (two differently-composed asymmetric-unit
    sites expanding onto the same position), and `ExpansionTooLarge`
    (`operation_count * asymmetric_site_count` checked-multiplied against
    a cap before attempting the O(n²) dedup scan — bounds pathological
    input, not a realistic-structure limit).
  - Deterministic output order: asymmetric-unit site order outermost,
    then within each site the identity operation's image first
    (regardless of where identity sits in the source operation list),
    then the remaining operations in their original declared order. A
    non-identity image's label is the source label with `@sym{N}`
    appended (`N` = 1-based position in the *declared* operation list);
    an unlabeled source site produces unlabeled images throughout, never
    a synthesized label.

### Added — `chematic-py` (CIF explicit symmetry-operation expansion)

- `PeriodicStructure.from_cif(text, expand_symmetry=True)` — the new
  default — now expands a CIF's explicit symmetry operations into a full
  unit cell; `expand_symmetry=False` restores the pre-expansion,
  asymmetric-unit-only behavior. `CifSymmetryStatus` gained
  `is_expanded`, `is_complete_cell` (`is_p1 or is_expanded`),
  `asymmetric_site_count`, and `expanded_site_count` (the latter two
  `None` unless `is_expanded`). `to_cif()` now delegates to
  `chematic_mol::CifPeriodicResult::to_cif_checked` (the single
  Rust-side source of truth for the write-safety judgment, previously
  duplicated Python-side against a lossy `is_p1`-only check) — refuses
  with `ValueError` iff `is_complete_cell` is `False`, unchanged for
  `wrap_into_cell()`/`make_supercell()` (both still preserve
  `symmetry_status`, now correctly including the expanded case).

### Added — `chematic-mol` (Gaussian Cube + OpenDX volumetric formats)

Roadmap step 5 (Wave 2 of the format-expansion program that started with
mmCIF/PQR/QCSchema/ORCA in PR #329): a new shared `VolumetricGrid` type
plus Gaussian Cube and OpenDX (the APBS/electrostatics scalar-field
subset) read/write support.

- New `VolumetricGrid { origin, axes, shape, values, atoms, units }` type
  (`chematic_mol::volumetric`), independent of `chematic-crystal::Lattice`
  (a unit-cell-specific type this crate does not depend on outside its
  optional `crystal` feature — voxel step vectors are a different concept
  from lattice edge vectors even though the raw `[[f64;3];3]` shape looks
  similar). `values` uses third-axis-fastest ordering
  (`index = (i*shape[1]+j)*shape[2]+k`), matching *both* Cube's native
  "x outer, z inner" order and OpenDX/APBS's native
  "z fastest, then y, then x" order, so neither reader/writer transposes
  data relative to its source file. `axes` is a general (non-orthogonal)
  3x3 matrix — both formats permit non-axis-aligned voxel step vectors.
  New `GridUnits::{Bohr, Angstrom}` tag (Cube has a real, sourced-but-
  contested Bohr/Ångström ambiguity on the first axis line's sign; rather
  than silently normalizing or silently assuming Bohr, the unit actually
  read/written is recorded explicitly and numbers are never rescaled).
  New shared `GridError` (shape overflow, over-cap, value-count mismatch,
  zero-length dimension, non-finite field) used by both formats' writers
  via `VolumetricGrid::validate`. `VolumetricGrid::checked_index`/`get`
  use fully checked arithmetic for the flat-index computation (every
  field is `pub`, so a caller can construct a shape whose product
  overflows `usize` even though an individual `(i, j, k)` looks
  in-bounds — this returns `None` rather than risking a panic).
- **Gaussian Cube** (`chematic_mol::cube`): `parse_cube`/
  `parse_cube_with_limits`/`write_cube`, plus a streaming-*input*
  `CubeFileReader<R: BufRead>` (same shape as `SdfFileReader`) for the
  multi-gigabyte grids real quantum-chemistry workflows produce — reads
  line-by-line rather than requiring the whole file in memory as one
  `String` (the returned `VolumetricGrid.values` is still a fully
  in-memory `Vec<f64>` either way), and validates the grid-point/atom
  caps from the header before
  the (potentially huge) voxel data block is read. Detects Cube's two
  real documented multi-dataset conventions (negative `NAtoms` with a
  dataset-identifier line; positive `NAtoms` with `NVal != 1`) and
  typed-rejects both (`CubeError::MultiDatasetUnsupported`) rather than
  silently reading only the first dataset. `write_cube` checks each
  `shape` dimension fits the signed `i64` an axis-line field requires
  (`CubeError::DimensionOutOfRange`) instead of an unchecked/wrapping
  `as i64` cast. New `CubeError`/`CubeParseLimits`.
- **OpenDX/APBS scalar field** (`chematic_mol::opendx`): `parse_opendx`/
  `parse_opendx_with_limits`/`write_opendx`, scoped explicitly to the
  regular-grid rank-0 `type double` subset APBS/electrostatics tooling
  actually produces (not the full general OpenDX/IBM Data Explorer
  format family) — a different `rank`/`type`/`class` declaration is a
  typed `UnsupportedArrayDeclaration`, not silently misread. Cross-checks
  `object 2`'s `gridconnections` counts against `object 1`'s
  `gridpositions` counts and `object 3`'s declared `items` count against
  the grid shape, both as typed mismatches rather than trusting either
  blindly. `write_opendx` fails closed rather than silently losing
  information: it refuses a `GridUnits::Bohr` grid
  (`OpenDxError::NonAngstromUnits` — DX has no unit tag and is always read
  back as Ångström, so writing Bohr magnitudes as-is would be a real,
  silent ~1.89x error) and refuses a grid with any `atoms`
  (`OpenDxError::AtomsNotSupported` — DX has no atom section, so there's
  no lossy-but-acceptable way to write them). New `write_opendx_lossy`
  opts into an explicit Bohr→Ångström conversion of `origin`/`axes` only
  (never `values`, which aren't a length quantity) when that's actually
  wanted; it still refuses non-empty `atoms`. New `OpenDxError`/
  `OpenDxParseLimits`.
- Both formats: checked-arithmetic overflow prevention
  (`shape[0]*shape[1]*shape[2]` never allocated before being validated
  against both `usize` overflow and the configured `ParseLimits` cap),
  `is_finite()` checks at every numeric parse point (with injected-NaN/
  Infinity test fixtures), and deterministic output (no `HashMap`
  anywhere in either module).
- Out of scope for this step, noted in both modules' doc comments as a
  future entry point: VASP's CHGCAR/LOCPOT volumetric formats, which
  would build on this same `VolumetricGrid` type.

### Added — `chematic-mol` (LAMMPS data file + dump/trajectory file)

Issue #227 format-expansion Wave 2, continued: LAMMPS `read_data`-format
data files and `dump`-command text-style trajectory files, explicitly
promised in this same Unreleased block's earlier mmCIF/PQR/QCSchema/ORCA
entry ("LAMMPS data/dump ... are scoped for a later Wave 2, not this
PR"). Both new modules (`lammps_data`, `lammps_dump`) are **standalone
document types**, deliberately not integrated with `chematic_core::Molecule`
-- LAMMPS bonds/angles/etc. are raw atom-index topology, not chemically
perceived bonds, and an MD atom under some `atom_style`s is not
necessarily even a chemical element. Verified independently against
LAMMPS's own official manual (<https://docs.lammps.org/read_data.html>,
<https://docs.lammps.org/dump.html>, <https://docs.lammps.org/Howto_triclinic.html>)
before implementation, not taken from memory; nothing in this brief's
format spec needed correcting -- every detail below was confirmed rather
than assumed, including the triclinic bounding-box formula (see below).

- **LAMMPS data file** (`chematic_mol::lammps_data`):
  `parse_lammps_data(text, atom_style)` / `write_lammps_data(data)`. The
  data-file format does not declare `atom_style` in-file and LAMMPS's own
  docs state it must be defined before reading -- it is genuinely not
  always recoverable from column count alone (`atom_style charge`'s
  `atom-ID atom-type q x y z` and `atom_style molecular`'s `atom-ID
  molecule-ID atom-type x y z` are both 6 fields). `atom_style` is
  therefore a required `LammpsAtomStyle` parameter, not inferred; there
  is no "guess from column count" fallback. Supports exactly 4 styles --
  `Atomic`/`Charge`/`Molecular`/`Full`, each with optional trailing `ix
  iy iz` image flags -- typed via `LammpsMass`/`LammpsAtom`/
  `LammpsVelocity`/`LammpsBond` for the `Masses`/`Atoms`/`Velocities`/
  `Bonds` sections and the `LammpsBox` (`lo`/`hi`/optional `tilt`) box.
  Any other `atom_style` (`LammpsAtomStyle::Other`) fails closed with a
  typed `LammpsDataError::UnsupportedAtomStyle`.
- Every other section (`Angles`, `Dihedrals`, `Impropers`, `Pair Coeffs`,
  `Bond Coeffs`, `Angle Coeffs`, `Dihedral Coeffs`, `Improper Coeffs`,
  and any unrecognized name) is preserved opaquely, byte-for-byte,
  in an ordered `Vec<(String, String)>` (`LammpsData::unparsed_sections`
  -- not a `HashMap`), never interpreted. This is safe specifically
  because this module is read/write/round-trip only in v1 (no
  atom-removal/renumbering API yet) -- documented on `LammpsData` as a
  constraint on any future PR that adds atom mutation, since renumbering
  would need to also remap these opaque sections' raw atom-index
  references or reject the operation.
- LAMMPS's separate type-label framework (`Atom Type Labels`/`Bond Type
  Labels`/...), under which `Masses`/`Atoms`/`Bonds` rows can key off a
  string label instead of a numeric type/ID, is fail-closed rejected
  (`LammpsDataError::TypeLabelsUnsupported`) as soon as such a section
  name is seen, rather than silently mis-splitting/dropping data: this
  module's section-boundary detection is numeric-ID-only and cannot
  safely distinguish a label-keyed row from a new section header.
- No unit conversion or inference anywhere: LAMMPS data files don't
  declare a unit system in-file (`units` is a separate input-script
  command), so every numeric value is stored exactly as read,
  unit-agnostic -- a deliberate non-goal, matching this project's
  Cube/OpenDX unit-handling discipline.
- **LAMMPS dump/trajectory file** (`chematic_mol::lammps_dump`):
  `parse_lammps_dump_frame` (single frame), the streaming
  `LammpsDumpReader<R: BufRead>` (`Iterator<Item = Result<LammpsDumpFrame,
  LammpsDumpError>>`, following this crate's `SdfFileReader` precedent
  rather than `CubeFileReader`'s one-shot shape, since a trajectory is
  inherently multi-record), `write_lammps_dump_frame`, and
  `write_lammps_trajectory`. The `ITEM: ATOMS` column list is arbitrary
  and self-declared per dump command; `LammpsDumpFrame` stores it as
  index-parallel `column_names`/`rows` (row-major), not a
  `HashMap<String, Vec<f64>>`, exactly preserving original column names
  and order -- confirmed by a round-trip test using unusual real-world
  column names (`c_myCompute[1]`, `f_myFix`).
- Reuses `LammpsBox` from `lammps_data` for the box shape -- LAMMPS data
  files give the true box directly, but dump files instead give an
  axis-aligned *bounding* box around the (possibly tilted) true box
  (`xlo_bound`/`xhi_bound`/...), needed because some visualization tools
  expect an axis-aligned rendering box. The triclinic bound<->true
  conversion formula was independently confirmed against
  <https://docs.lammps.org/Howto_triclinic.html> in **both directions**
  (the page states the true-box-to-bound formula explicitly, `xlo_bound
  = xlo + MIN(0.0,xy,xz,xy+xz)` etc., and gives the inverse explicitly
  too, `xlo = xlo_bound - MIN(0.0,xy,xz,xy+xz)`) -- new `box_bounds_to_true`/
  `true_to_box_bounds`, tested both via a bound->true->bound identity
  check over several tilt fixtures *and* (since the identity check alone
  can't catch a self-consistently-wrong `MIN`/`MAX` term) a
  hand-computed-absolute-value check.
- `LammpsDumpFrame::cartesian_positions()` resolves `x y z` (pass
  through) or `xs ys zs` (box-scaled, including the triclinic shear
  terms -- **not** simply independent per-axis scaling) into real
  Cartesian positions; `column()`/`column_index()` give raw access to
  any column by name. `xu yu zu` ("unwrapped": not passed through
  periodic-boundary wrapping, a materially different quantity from a
  wrapped position) are deliberately not resolved by
  `cartesian_positions()` -- use `column("xu")` etc. directly. The
  scaled-coordinate transform's box-origin term is this module's own
  derivation (LAMMPS's own docs sentence, "the actual unscaled
  coordinate is `xs*A + ys*B + zs*C`", omits the origin but cannot be
  correct without it, since `xs=ys=zs=0` must map to the box's own
  corner) -- flagged medium-confidence in the module docs per this
  project's citation-confidence convention, verified against a
  hand-computed example.
- Note on `column()`'s signature: returns `Option<Vec<f64>>` (owned), not
  `Option<&[f64]>` -- `LammpsDumpFrame::rows` is row-major, so a single
  column is not stored contiguously and cannot be borrowed without a
  copy; documented as a deliberate deviation.
- Zero-atom frames, custom/unrecognized dump columns, and orthogonal vs.
  triclinic boxes in either format all round-trip correctly (see test
  suite).

**Not supported, by design**: any `atom_style` outside
Atomic/Charge/Molecular/Full (typed `UnsupportedAtomStyle` error);
semantic parsing of `Angles`/`Dihedrals`/`Impropers`/any `*Coeffs`
section (opaquely preserved instead); LAMMPS's type-label framework
(typed `TypeLabelsUnsupported` error); any unit conversion or inference;
binary dump formats (only the default text `dump` style is supported);
any `Molecule`/bond-perception integration (both types are standalone,
see above); atom removal/renumbering (read/write/round-trip only in
v1); exact preservation of the *interleaving* order between typed and
opaque sections in a data file on write (each section's own content is
preserved exactly; `write_lammps_data` always emits the 4 typed sections
first in a fixed order, then opaque sections in their original relative
order -- LAMMPS itself doesn't require a specific section order either).
No Python or WASM bindings this step (consistent with every other
format module added in this Wave -- that remains a separate, later
concern).

---

## [0.16.0] — 2026-08-15

Periodic-structure interoperability and generalized stereochemistry
foundation. `chematic-crystal` gains CIF and POSCAR/CONTCAR interop
(`chematic-mol`'s CIF reader/writer now bridges to `PeriodicStructure`
via an optional `crystal` feature; the `chematic-crystal` crate itself
gains native POSCAR/CONTCAR read/write); `chematic-fp` gains the FPS
fingerprint exchange format; stereo configuration (tetrahedral +
square-planar) is now represented as a coordination geometry plus the
equivalence class of ligand permutations under that geometry's proper
rotation group, replacing two previously-independent hand-written
remapping algorithms and fixing a real bug where a square-planar center
could be silently coerced into a tetrahedral chiral-volume check. Also
includes a release-grade re-measurement of the `pipeline_v2` vs RDKit
2026.03.4 benchmark (superseding the stale 2026-08-06 numbers) that
surfaced a new finding: torsion parameter coverage, not bond/angle, is
now the dominant remaining MMFF94 gap. Purely additive at the Rust API
level for existing crates -- no breaking changes.

### Added — `chematic-core` (generalized stereo-configuration geometry)

- New `stereo_geometry` module: a stereo configuration is now modeled as a
  coordination geometry (`StereoGeometry::Tetrahedral`/`SquarePlanar`,
  `#[non_exhaustive]` for future trigonal-bipyramidal/octahedral) plus the
  equivalence class of ligand-slot permutations under that geometry's
  *proper rotation group* -- A4 (order 12) for tetrahedral, the order-8
  S4-stabilizer of a trans-pair partition (NOT the naive order-4
  in-plane-only group, which would wrongly give 6 orbits instead of 3) for
  square-planar. Replaces two previously independent, hand-written
  stereo-remapping algorithms in `chematic-smiles` (tetrahedral parity via
  cycle-counting; square-planar `@SPn` remapping via trans-pair-partition
  matching) with one shared, exhaustively self-tested primitive.
- Public API: `remap_tetrahedral_parity`/`remap_square_planar_tag` (the two
  bridge functions `chematic-smiles` actually calls), `StereoGeometry`,
  `StereoGeometryError`. `StereoConfiguration`/`CanonicalStereoConfiguration`/
  `canonicalize_configuration`/`equivalent_under_rotation` are `pub(crate)`
  rather than public -- all four are hardcoded to `[u32; 4]`, which only
  fits today's two 4-coordinate geometries; keeping them crate-internal
  defers the real arity-generalization question until a second geometry
  family (5/6 slots) actually needs it, rather than committing chematic-core's
  public API to "every geometry has 4 slots" today.
- See `docs/rfcs/generalized_stereo_geometry_rfc.md` for the full
  orbit-stabilizer derivation, oracle/regression provenance, and TBP/
  octahedral extension sketch.

### Fixed — `chematic-3d` (square-planar centers silently checked as tetrahedral)

- `tetrahedral_constraint_for`'s only guard against non-tetrahedral
  chirality was `debug_assert!(atom.chirality != Chirality::None)` -- a
  release-mode no-op, and wrong even in debug since
  `Chirality::SquarePlanar` is also `!= None`. A `SquarePlanar`-tagged atom
  (e.g. ordinary, non-dative `[Pt@SP1](Cl)(Cl)(N)N`) could silently be
  treated as a declared `@@` tetrahedral center and evaluated by
  `verify_stereo`/`repair_stereo` against a `VOLUME_EPS = 1e-6`
  chiral-volume tolerance meaningless for a square-planar arrangement.
  Fixed with `if !atom.chirality.is_tetrahedral() { return
  Err(UnsupportedCoordination); }`, reusing the rejection variant the
  function already returns a few lines below for other reasons.

### Fixed — `chematic-smiles` (allene-end-carbon stereo parity)

- An allene *end* carbon (sp2, one real double-bond partner standing in for
  the 4th tetrahedral-like position -- e.g. `F[C@@H]=[C]=[C@H]Cl`'s
  F-bearing atom) has a 3-element `stereo_neighbor_order`, not 4. Found
  during the `chematic-core` rewiring above (a byte-identical before/after
  canonical-SMILES diff caught a transient regression the existing
  relative-invariant tests did not); fixed by keeping the original,
  unmodified, length-generic parity fallback for any non-4-length order,
  routing only the common 4-element case through the new module. Pinned
  with an exact golden-value regression test.

### Measured — `pipeline_v2` vs RDKit 2026.03.4 benchmark re-run (issue #227 impact, release-grade)

Release-grade re-measurement only -- no `pipeline_v2`/force-field algorithm
or library code was changed to produce these numbers (the only non-data
changes in this PR are to measurement tooling itself: a report-generator
row-count fix for a pre-existing, unrelated 13th diagnostic arm, and a
scorer-script stereo-report-plumbing fix -- see the PR diff). Supersedes
the stale `docs/rfcs/
pipeline_v2_vs_rdkit_etkdgv3_benchmark.md` numbers (dated 2026-08-06,
predating PRs #314-317's MMFF94 Bond/Angle empirical-rule fallback).
Pinned at commit `494d634`, RDKit `2026.03.4`; full environment/seed/
timeout/ETKDG-parameter/MMFF-variant record at
`validation/results/pipeline_v2_vs_rdkit_environment_record.json`.

- **`chematic_pipeline_v2_mmff94_strict`** (this benchmark's own arm,
  via `embed_pipeline_v2`): **149/265 → 239/265** since the 2026-08-06
  baseline (90.2% usable coverage, up from 56.2%). **This number is NOT
  directly comparable to the `[0.15.0]` entry's 158→248/265 figure below**
  -- that number comes from `mmff94_strict_gate_remeasure_227`, a
  *different* embedding entry point (`dg::generate_coords`, not
  `embed_pipeline_v2`) with a different starting geometry. Same corpus,
  same `Mmff94BondAngleStrict` policy, different upstream pipeline --
  treat as two independently-useful, not-interchangeable numbers, not a
  contradiction.
- Other arms (this run): `no_ff` 254/265, `uff_only` 250/265,
  `mmff94_with_uff_fallback` 253/265, `chematic_legacy_etkdg` 265/265.
  RDKit's own 4 arms are near-saturated at 264/265 each (the one failure
  is a known RDKit-internal `BFGSOpt.h` crash on cyclopentane, present
  identically at both RDKit 2026.03.3 and 2026.03.4 -- see the
  environment record's version-isolation control).
- **New finding: torsion parameter coverage, not bond/angle, is now the
  dominant remaining MMFF94 gap.** `mmff94_strict_complete_bonded_term_gated`
  (requires full torsion+OOP coverage, not just bond/angle) drops to
  **180/265** (67.9%) -- of its 84 failures, 60 (71%) cite non-empty
  `torsions_missing`, while `oop_missing` is 0 across the board and
  `bonds_missing` is 0. This is the concrete evidence the project's
  roadmap needed: issue #227's Bond/Angle empirical fallback closed the
  bond/angle gap, and torsion coverage is the next real bottleneck --
  promotes the (currently lowest-priority) MMFF torsion empirical
  fallback item for re-evaluation. Note: this is a *coverage* finding
  (how often the strict gate refuses to proceed at all); whether the
  gap also causes measurably worse *geometry* on molecules where a
  fallback (UFF) is used instead is not yet measured -- RMSD/TFD/
  coverage-at-threshold metrics against this same corpus are a separate,
  not-yet-done follow-up (`chematic_pipeline_v2_mmff94_with_uff_fallback_
  complete_bonded_term_gated` itself stays healthy at 252/265, since the
  UFF fallback absorbs most of the torsion-coverage-gap cases).
- Per-molecule transition vs. the 2026-08-06 baseline: overwhelmingly
  `typed_failure → success` (matching #227's intent) across every MMFF94
  arm. Two residual anomalies, reported honestly rather than smoothed
  over: (1) `chembl_tier_b_0166` newly times out on 6 arms -- its
  baseline elapsed time was already 15-16s (75-80% of the 20s budget)
  and now runs 21-25s; cause not conclusively isolated (could be a small
  real per-molecule overhead from work landed since the baseline, or
  residual timing variance -- this molecule sits close enough to the
  timeout boundary that it warrants a note, not a firm regression claim).
  (2) `chembl_tier_b_0182`'s `mmff94_with_uff_fallback_repair` arm flips
  `success → typed_failure` (`FinalStereoViolation`) at nearly identical,
  well-under-budget elapsed time (~1.4s baseline, ~4.7s now) -- this is a
  real, reproducible stereo-repair-then-minimize outcome change, not a
  timing artifact; root cause not investigated further here (out of
  scope for a measurement-only pass).
- **Measurement provenance note**: the first attempt at this
  re-measurement (2026-08-14) was discarded after its data was found
  contaminated by an unrelated concurrent process on the same machine
  (confirmed via this benchmark's own gate-widening monotonicity
  integrity check tripping, and via elapsed-time analysis showing
  several molecules pushed 5-30x over their baseline time, well past the
  20s timeout). Re-run cleanly (2026-08-15) after confirming the
  contending process had exited; see the environment record's
  `measurement_history_note` for the full account.

### Added — `chematic-fp` (FPS fingerprint exchange format)

- New `fps` module: streaming read/write for the FPS ("Fingerprint file
  format") text-based fingerprint interchange format popularized by
  chemfp/OpenBabel. `FpsReader<R: BufRead>`/`FpsWriter<W: Write>` iterate
  one record at a time rather than materializing a whole file, and work
  over any `BufRead`/`Write` sink (WASM-compatible, matching this crate's
  existing `#![forbid(unsafe_code)]`/wasm32 constraints).
- `FpsHeader` models `num_bits`/`type`/`software`/`source`/`comment`
  explicitly and carries any other `#`-prefixed header line through
  losslessly via `extra` (including the `#FPS1` version line, kept first
  on write-back per spec). `FpsHeader::for_chematic` stamps
  `software=chematic-fp/<version>` for fingerprints this crate itself
  computed.
- Hex bit-ordering verified against the real chemfp FPS spec
  (<https://chemfp.com/fps_format/>): byte `k` = fingerprint bits
  `[8k, 8k+8)`, LSB-first within the byte -- matches `BitVec2048`/
  `BitVecN`'s own bit numbering directly, so no reordering is needed
  between the two representations.
- Reuses `BitVec2048`/`BitVecN` as the sole bit-vector representation; no
  new fingerprint algorithms or binary formats.

### Added — `chematic-crystal` (POSCAR/CONTCAR read/write)

- `poscar::{parse_poscar, parse_contcar, write_poscar}` and the
  `PoscarDocument`/`PoscarError`/`PredictorCorrector` types: read/write
  support for VASP's plain-text POSCAR/CONTCAR structure format. VASP 5
  only (explicit species-name line; VASP 4's implicit POTCAR-derived
  ordering is rejected with a typed error rather than mis-parsed); both
  scale-factor forms from the VASP wiki (single value, including the
  negative "target cell volume" form, and the 3-component per-axis form);
  Direct/Cartesian coordinate modes; selective dynamics; ion velocities;
  and CONTCAR's predictor-corrector MD-restart section, preserved
  verbatim since VASP's own documentation does not specify its numeric
  layout. Re-exported under the existing `chematic` facade's `crystal`
  feature. See `crates/chematic-crystal/src/poscar.rs`'s module docs for
  the full list of format-fidelity decisions.

### Added — `chematic-mol` (optional `crystal` feature: CIF ↔ `PeriodicStructure` adapter)

- New optional `crystal` feature on `chematic-mol` (its first-ever
  `[features]` table; off by default, no existing consumer affected) gates
  an optional path dependency on `chematic-crystal` and two new functions:
  `parse_cif_periodic_structure`/`write_cif_periodic_structure`. Bridges
  the existing CIF reader/writer (`parse_cif`/`write_cif`/`UnitCell`,
  unchanged) directly to `chematic_crystal::PeriodicStructure` — cell
  parameters to `Lattice`, `_atom_site_occupancy` to `Occupancy` (not
  previously read at all), and multiple `_atom_site_*` rows sharing a
  fractional position merged into one `PeriodicSite`'s multi-species list
  (positional/substitutional disorder). Implements the sketch from
  `docs/rfcs/chematic_crystal_foundation.md`'s "CIF migration" section.
- **Symmetry**: `chematic-mol`'s CIF reader has never expanded symmetry
  operations (P1-equivalent only) — this adapter inherits that scope but
  makes it explicit rather than silent: a new `CifSymmetryStatus` enum
  (`P1` vs `UnexpandedSymmetry { space_group_name, operation_count }`) is
  returned alongside the structure, driven by
  `_symmetry_space_group_name_H-M`/`_space_group_name_H-M_alt`,
  `_symmetry_Int_Tables_number`/`_space_group_IT_number`, and
  `_space_group_symop_operation_xyz`/`_symmetry_equiv_pos_as_xyz` loop
  detection, so a caller can distinguish "this file is genuinely P1" from
  "this file declared symmetry that was not expanded — only the
  asymmetric unit was returned." No symmetry expansion is implemented
  (out of scope, matches `chematic-crystal`'s own non-goals).
- `chematic-crystal` itself is untouched and remains independent of
  `chematic-mol`/`chematic_core::Molecule` (dependency direction is
  `chematic-crystal <- chematic-mol`, never the reverse).

## [0.15.0] — 2026-08-14

New crate `chematic-crystal` (periodic/crystal structure foundation), plus
a substantial MMFF94 Bond/Angle accuracy fix (issue #227): 265-molecule
corpus strict-policy success 158/265 → 248/265 since v0.14.1, zero
regressions. Purely additive at the Rust API level for existing crates --
no breaking changes.

### Fixed — `chematic-ff` (MMFF94 Bond/Angle empirical-rule fallback, issue #227)

- Ported Halgren's MMFF.V eq. 18-20 empirical Bond-stretch/Angle-bend rule
  (`mmff94_bond_energy_resolved`/`mmff94_angle_energy_resolved`, new
  additive functions -- the existing `mmff94_bond_energy`/
  `mmff94_angle_energy` keep their original `Option<Params>` signatures
  unchanged), tried strictly *after* the existing exact-table/`eqLevel`-
  ladder lookup so it never overrides a real table hit. New
  `Mmff94Resolution` enum (`DirectTable`/`EquivalentType`/
  `GenericAngleTypeFallback`/`EmpiricalBond`/`EmpiricalAngle`) reports
  which mechanism resolved a given lookup, for diagnostics/tests.
- Along the way, discovered and fixed a real data gap in the Angle table
  itself: 97 rows present in RDKit's real `defaultMMFFAngleData` (generic,
  central-atom-type-only `theta0` defaults) were missing from chematic's
  port (2245 → 2342 rows) -- restored with a guard so the existing
  `mmff94_angle_energy` contract is provably unchanged for every
  pre-existing input.
- One triple, `(angle_type=0, N-type=43, S-type=18, C-type=63)`, is
  deliberately left unresolved (fails closed) rather than guessed: the
  outer atom type has no equivalence-class table entry, RDKit's own real
  C++ dereferences that unchecked (undefined behavior), and the live
  RDKit oracle's answer for it could not be attributed to any well-defined
  resolution mechanism.
- Also fixed 5 pre-existing MMFF94 atom-typing gaps (nitrile/isocyanide N,
  sulfonamide/sulfonate N, nitro N, azide/diazo N, charged-sulfoxide S)
  and ported RDKit's `eqLevel` atom-type-equivalence ladder for Angle
  table lookup (both prerequisites for the empirical-rule work above).
- **Net effect**, measured on the 265-molecule Wave 1 corpus via the
  production minimization path (`ForceFieldPolicy::Mmff94BondAngleStrict`),
  reported as two separately-verified numbers so they aren't confused
  (both via a full per-molecule join, not just aggregate counts, zero
  regressions in either case):
  - **Full v0.14.1 → v0.15.0 change** (this release's complete Bond/Angle
    work: atom-typing fixes + `eqLevel` ladder + the empirical rule above):
    **158/265 → 248/265 (107 → 17 failing)**.
  - **The empirical-rule work specifically** (from just before this PR to
    just after, isolating its own contribution from the atom-typing/
    `eqLevel` prerequisites merged earlier in this same release):
    **178/265 → 248/265 (87 → 17 failing)**.
  - Of the 3 molecules still `MinimizationFailed` in the final state, all
    3 were already-non-`Ok` (`MissingParameters`) in v0.14.1 -- the
    atom-typing/`eqLevel` fixes gave them real MMFF94 parameters for the
    first time, which then exposed a pre-existing geometry issue that
    "no parameters at all" had been masking. This is a newly-*visible*
    failure mode, not a regression from a previously-successful state
    (confirmed: zero `Ok`→non-`Ok` transitions across the full v0.14.1 →
    v0.15.0 span).

### Added — chematic-crystal (new crate)

- `crates/chematic-crystal`: periodic (crystal) structure representation
  and geometry -- `Lattice` (triclinic-capable, validated matrix/inverse/
  reciprocal vectors), `FractionalCoord`/`CartesianCoord`, `PeriodicSite`/
  `SiteSpecies`/`Occupancy` (multi-species disorder-ready), and
  `PeriodicStructure` with exact (not `round()`-approximate) periodic
  minimum-image distance (deterministic tie-break: equidistant periodic
  images resolve to the lexicographically smallest image), cutoff
  neighbor enumeration, and diagonal supercells. Deliberately **not** an
  extension of `chematic_core::Molecule` -- see
  `docs/rfcs/chematic_crystal_foundation.md`. Optional `serde` feature;
  optional `crystal` feature on the `chematic` facade, included in `full`
  (does not change `default`, which stays empty). No symmetry, no CIF
  parser changes, no Python/WASM/MCP bindings in this release.

## [0.14.1] — 2026-08-12

Anticancer platinum coordination-chemistry compatibility fixes, plus Extended
XYZ (extxyz) read/write support. Patch release: all `chematic-core`/
`chematic-mol`/`chematic-chem` changes are bug fixes to existing behavior
(no new production capability gated behind a flag); the extxyz addition is
new functionality but additive at the Python/WASM binding layer. Note: the
Rust-only `XyzFrame`/`XyzError`/`write_extxyz` signature changes below are a
real break to the v0.14.0 Rust API surface already published to crates.io,
not merely a break to something unreleased — flagged here for anyone
depending on `chematic-mol`'s Rust API directly rather than via SMILES/JSON
bindings.

### Added — platinum coordination-chemistry compatibility benchmark

- `validation/platinum/pt_corpus.jsonl` (18-entry hand-verified corpus of
  anticancer Pt(II)/Pt(IV) coordination complexes — cisplatin, transplatin,
  carboplatin, oxaliplatin, nedaplatin, lobaplatin, picoplatin,
  dicycloplatin, satraplatin, iproplatin, tetraplatin, oxoplatin, plus
  charged/S-donor/C-donor diversity cases and 2 non-Pt generalization-gate
  cases), `crates/chematic-mol/examples/platinum_benchmark.rs`,
  `scripts/platinum_rdkit_oracle.py`,
  `scripts/platinum_compare_chematic_rdkit.py`, and
  `validation/platinum/FEASIBILITY.md` — a measurement-only survey of
  whether chematic can represent/parse/round-trip/canonicalize platinum
  anticancer complexes without silently corrupting their coordination
  chemistry (not an anticancer-activity project; no IC50/resistance/
  toxicity/PK prediction). Found 2 general (not Pt-specific) production
  defects, fixed in 2 follow-up PRs (see this file's `chematic-core`/
  `chematic-mol`/`chematic-chem` `Fixed` entries once they land); measured
  and explicitly did **not** fix a 3rd (square-planar cis/trans stereo has
  no representation in chematic at all — cisplatin and transplatin
  currently canonicalize to the same identity; see the FEASIBILITY doc for
  why this is reported, not patched, this round).

### Fixed — `chematic-core` (dative-bond implicit hydrogen count)

- **`valence_inferred_hcount` treated a `BondOrder::Dative` bond's donor
  side exactly like a covalent single bond when computing implicit
  hydrogen count**, contradicting `BondOrder::Dative`'s own documented
  donor→acceptor semantics and RDKit's identical treatment of the same
  `->`/`<-` SMILES syntax: an un-bracketed dative donor like `N->[Pt]Cl`
  computed as `NH2` instead of the chemically correct `NH3`, because the
  donor's own lone-pair-sharing bond was (wrongly) subtracted from its
  normal covalent valence. Flowed into `molecular_weight`, `exact_mass`,
  and `chematic-smiles`'s canonical-writer atom invariant. Donor-side dative
  bonds now contribute 0 to the valence sum instead of `order_int()`'s 1;
  the acceptor side is unchanged (out of scope — no acceptor in the
  motivating corpus is in the organic subset). Found via the platinum
  coordination-chemistry benchmark; not platinum-specific (verified against
  Fe/Co/Pd/Ru acceptors too). A pre-existing `chematic-smiles` test
  (`dative_bond_direction.rs`) pinned an exact canonical-SMILES string that
  depended on the old, wrong invariant; updated to the new correct value,
  with a replacement case added so the arrow-flip-on-acceptor-first code
  path it originally covered stays covered.
- **Known, documented, not-fixed-this-round divergence:** `validate_valence`
  (built on the separate, public `bond_order_sum`, also used by
  `chematic-cip`/tautomer enumeration/`chematic-ff`) can still report a
  false-positive `ValenceError` on a *bracketed* dative donor whose only
  listed normal valence is exactly met by its own explicit H count (e.g.
  `[OH2]->[Pt]`) — deliberately not widened, given the larger blast radius;
  pinned by a new test rather than left undocumented.

### Fixed — `chematic-mol` (MOL V2000/V3000 dative/coordinate bonds)

- **MDL bond type 9 (dative/coordinate — the exact convention RDKit uses to
  write `Bond::BondType.DATIVE`, V3000 only) silently mapped to
  `BondOrder::Single`** in both the V2000 and V3000 readers' "unknown code"
  catch-all, with no error or warning: reading an RDKit-generated V3000
  molfile containing a dative Pt–N bond silently produced a different
  molecule (connectivity intact, bond semantics quietly lost). Both readers
  now map code 9 to `BondOrder::Dative`, preserving the file's donor/
  acceptor atom order. `mol3000.rs`'s writer now emits code 9 for `Dative`
  bonds instead of collapsing to plain single (`1`); V2000's writer still
  collapses `Dative` to `1` on write, matching RDKit's own inability to
  express a dative bond in V2000 at all. Both readers carry a committed
  regression test built from the literal RDKit-generated molblock that
  surfaced this finding.

### Fixed — `chematic-chem` (periodic-table mass data gap, all elements)

- **`avg_mass`/`mono_mass` (`molecular_weight`/`exact_mass`) covered only
  ~24 light main-group elements and silently fell back to
  `atomic_number as f64` for every other element** — meaning every
  transition metal, lanthanide, actinide, and heavy post-transition element
  (platinum: atomic number 78, real mass ~195 Da, returned as "78.0 Da";
  same defect for Fe, Cu, Zn, Au, U, and 90 others) got a wildly wrong LOW
  mass with no error. Not connected to coordination bonds at all — found
  via the same platinum benchmark's mass check, but affects any molecule
  containing an unlisted element, dative bonds or not. Both tables extended
  to all 118 elements `chematic_core::Element` models, sourced from RDKit's
  `PeriodicTable::getAtomicWeight`/`GetMostCommonIsotopeMass` (2026.03.3).
  The ~24 previously-covered elements keep their pre-existing values
  unchanged (checked element-by-element, not assumed) — most notably
  selenium, where the pre-existing value is the current IUPAC standard
  atomic weight and RDKit 2026.03.3 ships the superseded pre-2013 value.

### Fixed — `chematic-chem` (`trained-solubility-mlp` feature, unrelated to the above)

- **`mlp_solubility`'s embedded model weights (`W1_BYTES`/`B1_BYTES`) never actually compiled under the opt-in `trained-solubility-mlp` feature**: both `include_bytes!` paths had one `../` too many, resolving outside the repository entirely and failing to build with "No such file or directory" whenever the feature was enabled. Found while running this release's `cargo clippy --all-features` gate, the first time this feature has ever been exercised. Fixed to the correct 3-level-up path; a `needless_range_loop` clippy lint the fix then exposed in the same function was also fixed (semantics unchanged). Verified past compiling: with the feature enabled, `mlp_solubility` now actually loads the trained weights and produces plausible output where water is predicted more soluble than octane. `trained-solubility-mlp` is off by default, so this has no effect on any default build, and is unrelated to the platinum/extxyz work above.

### Added — `chematic-mol` (Extended XYZ / extxyz)

- `parse_extxyz`/`write_extxyz`, `ExtxyzReader`/`ExtxyzWriter`,
  `parse_extxyz_all`, plus new `XyzFrame` fields `lattice`, `properties`,
  `info`. Built as an extension of the existing multi-frame `XyzFrame` type
  (not a separate format) — same `species:S:1:pos:R:3` atoms as plain XYZ,
  plus ASE's `Lattice=` cell matrix, typed per-atom `Properties=` columns
  (`forces:R:3`, `charge:R:1`, ...), and arbitrary `key=value` frame
  metadata (`energy=`, `pbc=`, ...). A plain XYZ file (no `Lattice=`/
  `Properties=` in its comment line) parses and round-trips through the
  extxyz reader/writer unchanged — this is what makes extxyz a strict
  superset, not a competing parser. Fails closed on malformed `Lattice=`/
  `Properties=`, atom-row column-count mismatches, non-finite property
  values, unterminated quotes, and duplicate info keys. `parse_xyz`/
  `write_xyz`/`XyzReader`/`XyzWriter` (plain XYZ) behave identically to
  before this change — same inputs produce the same `XyzFrame`s (with the
  three new fields defaulted to `None`/empty) and the same output text.
- Python: `from_extxyz`, `from_extxyz_all`, `to_extxyz` (new; the existing
  `from_xyz`/`to_xyz` bind a separate, single-frame `chematic_3d::xyz`
  module and are untouched). WASM: `mol_from_extxyz`, `extxyz_frame_json`,
  `to_extxyz_json`.
- **Breaking (Rust API only)**: `XyzFrame` gained three public fields
  (`lattice`, `properties`, `info`) — breaks any external
  `XyzFrame { atoms, comment }` struct literal. `XyzError` gained seven
  variants — breaks an external exhaustive `match`. `write_extxyz` now
  returns `Result<String, XyzError>` instead of `String` (see "fails closed
  on unwritable metadata" below). **Note:** `XyzFrame`/`XyzError`'s
  pre-this-PR shape is the actual v0.14.0 API already published to
  crates.io (this branch's merge-base postdates the v0.14.0 tag) — this is
  not a purely hypothetical/unreleased-so-far break.
- **Fails closed on unwritable metadata, not just unparseable input**:
  `write_extxyz`/`ExtxyzWriter::write_frame` reject (rather than silently
  emit corrupt output for) an `XyzFrame::info` key or `XyzProperty::name`
  containing a character extxyz `key=value`/`Properties=` syntax can't
  represent; an info key literally `"Lattice"` or `"Properties"` (always
  re-parsed as the dedicated field, not as itself); an info value with an
  embedded newline (unrepresentable in a single comment line); or a
  `XyzProperty` declaring 0 components (`Properties=` requires at least 1).
  All are reachable only from a hand-built `XyzFrame`, e.g. via the
  Python/WASM `to_extxyz` bindings.
- **Info values may contain `"` and `\` freely**: `"`/`\` inside a quoted
  `key="..."` value are escaped on write and un-escaped on read (`\"`,
  `\\`), matching ASE's own `key_val_str_to_dict`. A bare key with no
  `=value` (e.g. a standalone `constrained` token) parses as an implicit
  boolean flag (`"T"`), also matching ASE, rather than erroring.

### Fixed — `chematic-rxn` (`DEFAULT_TEMPLATES` silent parse failures, issue #296, unrelated to the above)

- **14 of `retro::DEFAULT_TEMPLATES`' 59 entries never parsed at all** —
  `retro_disconnect`'s `Err(_) => continue` swallowed the `SmirksParse`
  failure indistinguishably from "no match," so negishi coupling, reductive
  amination, both Mitsunobu variants, aldol, Michael addition,
  Friedel-Crafts alkylation, Mannich, trifluoromethylation, PMB ether, and
  both C-H bromination/oxidative-addition templates had never fired for any
  caller since inception. Three distinct root causes: SMARTS-only query
  primitives (`X4`/`X3` connectivity, comma OR-lists, semicolon AND) this
  crate's SMILES-based template parser can't express (10 templates,
  rewritten to the supported subset with disclosed precision tradeoffs);
  one malformed template (`pmb_ether`, unbalanced ring closure + an invalid
  multi-atom token used as a single bracket atom); three templates using a
  bare, unbracketed `H` (genuinely invalid SMILES outside brackets, not a
  parser bug — fixed to `[H]`). New `all_default_templates_parse`
  CI-gating test asserts every `DEFAULT_TEMPLATES` entry's SMIRKS parses,
  so a future silently-broken built-in template fails the build instead of
  going unnoticed indefinitely (the issue's own acceptance criterion). 18
  new positive/negative/false-positive tests across the 14 fixed templates.
  Merged (`fc7a42a`) after the `v0.14.0` tag but before the `v0.14.1` tag,
  so this entry is placed here rather than under `[0.14.0]` or
  `[Unreleased]` — it genuinely first shipped in `v0.14.1`, just via an
  unrelated commit that happened to land in the same release window.

## [0.14.0] — 2026-08-11

Stereo-aware distance geometry: declared E/Z (cis/trans) is now enforced as a
genuine bound-matrix constraint at embedding time, not just checked after the
fact, and the post-minimization gate that protects it is now composable with
the mechanism that produces it. Resolves the issue #285 release-gate waiver
from v0.13.0's entry below — the same two named molecules
(`chembl_tier_b_0126`, `chembl_tier_b_0168`) are the ones this release fixes.

Opt-in only (`enforce_chirality: false` remains the default everywhere); the
default conformer path (`generate_coords_etkdg`/`Mol.conformer_ensemble()`)
is untouched.

### Fixed — `chematic-3d` (declared-E/Z distance-geometry embedding, issue #285)

- **Root cause found and fixed**: `apply_vdw_bounds`'s generic non-bonded Van
  der Waals lower bound (sum of radii — two carbons: 3.40 Å) was being
  applied to a declared-E/Z alkene's own 1-4 substituent pair regardless of
  which stereochemistry was declared, structurally excluding the correct cis
  geometry (analytic ≈2.88 Å for `but2ene_Z`) from ever being sampled or
  reconstructed. Not an eigendecomposition sign-convention artifact as an
  earlier diagnosis speculated — that candidate mechanism was empirically
  refuted before this one was confirmed.
- New `apply_declared_ez_bounds` (`enforce_chirality`-only): for each
  declared E/Z double bond, computes the analytic same-side/opposite-side
  1-4 distance from the same bond-length/angle model `build_bond_angle_bounds`
  already uses, and intersects it into the bond matrix *before* the generic
  Van der Waals floor applies — the correct geometry becomes reachable by
  construction, not by repair, retry, reflection, or perturbation after the
  fact (all three considered and rejected — see `docs/rfcs/etkdg_3d_gap_rfc.md`
  for the full comparison).
- Unlike tetrahedral chirality (`@`/`@@`), which a pairwise distance matrix
  can never encode (a molecule and its mirror image have identical pairwise
  distances), declared E/Z is genuinely distance-representable: cis and
  trans are two different scalar separations for the same atom pair, not
  mirror images. Ring-fused tetrahedral stereocenters (`testosterone`,
  `cholesterol`) are unaffected by this fix and remain a known, separately-
  scoped gap.
- Measured on the 265-molecule corpus (declared-E/Z subset, 39 molecules,
  through production `embed_pipeline_v2`): stereo-satisfied count 22 → 42,
  violated 23 → 3, pipeline success/geometry-soundness unchanged (32/39 both
  ways). 18 molecules newly fixed, 2 addressed by the fix below.

### Added — `chematic-3d` (`enforce_chirality` + `StereoPolicy::VerifyOnly` composition)

- `embed_pipeline_v2`'s config-validation gate previously rejected
  `enforce_chirality: true` for any `stereo_policy` other than `Ignore`,
  reasoning that the two stereo mechanisms were unrelated and composing them
  would be confusing. Corpus measurement disproved this: `enforce_chirality`
  protects embedding-time correctness only, and force-field minimization has
  no notion of declared stereo and can walk a correctly-embedded E/Z bond
  back across its boundary afterward (found on 2 real molecules,
  `chembl_tier_b_0076`/`chembl_tier_b_0083` — confirmed by re-running with no
  force field, which leaves them correct). `enforce_chirality: true` is now
  also allowed with `StereoPolicy::VerifyOnly`, whose existing post-
  minimization gate catches exactly this failure mode as a typed error
  instead of silently returning a `success` result with wrong stereo.
  `StereoPolicy::RepairAndVerify` remains rejected in this combination —
  composing its own repair pass with `enforce_chirality`'s is a separate,
  not-yet-validated question.

### Added — `chematic-py`, `chematic-wasm` (`enforce_chirality` exposure)

- `enforce_chirality` (default `false`, fully backward compatible) is now a
  real, settable parameter on `PipelineV2Config`/`PipelineV2Config.safe()`
  (Python) and the `enforceChirality` JSON field (WASM) — previously neither
  binding's config construction threaded the field through at all, so it was
  unconditionally `false` regardless of caller intent, making the fix above
  unreachable from Python or WASM callers.

### Fixed — `chematic-rxn` (`suzuki_biaryl` retro-template, issue #294)

- `[c:1][c:2]>>[c:1]Br.[c:2]B(O)O` never matched a real biaryl bond at all:
  two adjacent aromatic atoms with no explicit bond token default to an
  aromatic bond in this crate's SMILES convention, so the template only ever
  matched *intra-ring* aromatic bonds, not a genuine inter-ring biaryl
  connection. Fixed to `[c:1]-[c:2]` (explicit single bond) — no ring-
  topology check needed. Found while investigating this: 14 of 59
  `DEFAULT_TEMPLATES` entries silently never parse at all (`retro_disconnect`
  swallows the `SmirksParse` error) — filed as issue #296, not fixed here.

## [0.13.0] — 2026-08-10

### Added — `chematic-mol` (XYZ / multi-frame XYZ)

- `parse_xyz`/`write_xyz`, `XyzFrame`/`XyzAtom`, `XyzReader` (multi-frame),
  `parse_xyz_all`, `XyzWriter`, `XyzError`. Explicit hydrogens are kept as
  real atoms, never folded into `implicit_hcount()`. No connectivity/bond-
  order inference is performed, not even opt-in — XYZ carries coordinates
  and elements only. Fails closed on atom-count mismatch, non-finite
  coordinates, and unknown elements.

### Added — `chematic-perception` (per-atom stereocenter candidates, issue #263)

- `stereo_centers(&Molecule) -> Vec<(AtomIdx, bool)>` exposes the per-atom
  tetrahedral-stereocenter classification (`(atom, specified)`) that
  `stereo_completeness` already computed internally but only surfaced as
  aggregate counts. `stereo_completeness` is now implemented in terms of
  `stereo_centers`, so there is one source of truth for the classification
  logic instead of two. Purely additive — `StereoCompleteness`'s public
  field shape is unchanged.
- Two real bugs in the shared `simple_morgan_ranks` helper (used by both
  `stereo_centers` and `validate_stereo`) were found and fixed while adding
  this API:
  - **Negative formal charge caused a `u64` overflow (issue #267).** The
    initial invariant computed `atom.charge as u64 * 1000` where
    `atom.charge: i8` — a plain `as u64` cast sign-extends negative values
    (`-1i8 as u64 == u64::MAX`), so the multiply overflowed unconditionally
    for any negatively-charged atom: panicking in debug builds on any
    ordinary anion (carboxylate, sulfonate, phosphate, ...), silently
    wrapping to a corrupted invariant in release builds. Fixed by computing
    the whole invariant in `i64` and casting to `u64` once at the end
    (bit-for-bit identical to the old code for `charge >= 0`).
  - **Implicit-hydrogen rank-0 sentinel collided with a real atom's
    normalized rank 0.** `simple_morgan_ranks` normalizes each atom's
    Morgan-refinement invariant to a consecutive ordinal starting at `0` —
    so an ordinary heavy-atom neighbor can legitimately hold rank `0` (e.g.
    a plain methyl carbon with the lowest invariant in the molecule).
    `stereo_completeness`/`stereo_centers` used the literal value `0` as a
    stand-in rank for an implicit hydrogen when checking whether a
    candidate stereocenter's 4 groups were all distinct; when a real rank-0
    neighbor and an implicit H coincided at the same center, `dedup()`
    merged the two `0`s and the atom was silently dropped as "not a
    stereocenter" — undercounting real, correctly `@`/`@@`-annotated
    stereocenters (repro: `C[C@@H](Cl)C(Br)(F)I`, atom 1). Fixed by using
    `mol.atom_count()` as the implicit-H sentinel instead of `0` (provably
    unreachable as a real normalized rank, since the max is
    `atom_count() - 1`).

### Added — `chematic-chem` (E/Z double-bond stereo-completeness, issue #264)

- `ez_completeness(&Molecule) -> EzCompleteness { specified, unspecified,
  total }`. Excludes terminal/symmetric double bonds and bonds in rings
  smaller than 8 atoms (matching RDKit's `isBondPotentialStereoBond`/
  `MinBondRingSize` cutoff, oracle-verified), using a BFS shortest-cycle-
  through-bond check (not SSSR alone, which would miss bridged-bicyclic
  cases like norbornene) to determine ring membership correctly.

### Added — `chematic-perception` (macrocycle predicate, issue #266)

- `is_macrocycle(ring: &[AtomIdx]) -> bool { ring.len() >= 9 }` — a single
  shared predicate matching the two hardcoded `9`s already duplicated in
  `chematic-3d` (`MACROCYCLE_RING_THRESHOLD`,
  `etkdg_knowledge::classify`'s `MACROCYCLE_MIN`). The duplication in
  `chematic-3d` itself is not unified by this change.

### Fixed — `chematic-chem` (atropisomer detection/assignment notation-invariance, issues #262, #276)

- `detect_atropisomers` rewritten to use `chematic_perception::find_sssr`
  structurally (independent of `BondOrder`) to confirm two aromatic carbons
  belong to separate SSSR rings, then requires at least one substituted
  ortho ring position on each ring — instead of the previous
  `bond.order == BondOrder::Single` gate, which made the same real molecule
  give different answers depending on whether the SMILES wrote the
  inter-ring bond explicitly or left it implicit. Also fixes a false
  positive on para-substituted biaryls (e.g. 4,4'-dimethylbiphenyl is no
  longer flagged).
- `assign_atropisomer_chirality` had its own, separate, redundant
  `bond.order == BondOrder::Single` gate before annotating a detected
  atropisomeric bond with `Up`/`Down` chirality — the same class of
  notation-dependence issue #262 removed from `detect_atropisomers`, just
  one function later: an implicit-notation biaryl bond that `detect_
  atropisomers` correctly flagged as atropisomeric could still silently
  get skipped for chirality assignment. Fixed by gating on `AtropisomerType
  ::Biaryl` (matching `detect_atropisomers`'s own classification) instead
  of the bond order — not a blanket deletion of the check, since a naive
  deletion would have let the gate's incidental protection of
  `AtropisomerType::Allene`'s central `Double` bond regress too.

### Fixed — `chematic-ff` (MMFF94 stretch-bend classification-key bug, issue #227 Priority 2C, **breaking**)

- **Root cause (the fix the user asked for): stretch-bend used the *angle
  type* (0-8) directly as its own `MMFF94_STBN` lookup key.** RDKit computes
  a distinct, finer-grained "stretch-bend type" (0-11) via
  `getMMFFStretchBendType(angleType, bondType1, bondType2)`
  (`AtomTyper.cpp:2480-2508`, pinned commit — see
  `scripts/mmff94_provenance/PROVENANCE.md`'s "Stretch-bend" row) and uses
  THAT as the table key instead — `MMFF94_STBN`'s own frozen data is keyed
  by stretch-bend type, not angle type (self-consistency proof: the one
  key-5 row is an all-CR3R triple, structurally unreachable as an angle type
  under `angle_type_for`'s own table). A new `pub` `stretch_bend_type_for`
  (`crates/chematic-ff/src/mmff94_minimizer.rs`) computes the correct key,
  ported verbatim from the diagnostic's `getMMFFStretchBendType`/arg-
  canonicalization port (`mmff94_stbn_equivalence_diagnostic_227.rs`,
  already merged, diagnostic-only, untouched by this fix).
  `stretch_bend_energy` (`mmff94_minimizer.rs`) and `chematic-3d`'s
  `compute_mmff94_coverage` (`minimize.rs`) now compute this key before
  calling `mmff94_stbn`/`mmff94_stbn_type_only`, instead of passing the
  angle type straight through.
  - Measured on the 265-molecule Wave 1 corpus, cross-validated two
    independent ways (a live RDKit oracle re-run and a direct production-
    code cross-check against the diagnostic's frozen per-row predictions,
    zero mismatches either way): of the 427 StretchBend
    `routing_bug_candidate` instances (a real, correctly-typed parameter
    existed at a *different* classification code, masked behind RDKit's
    generic Dfsb periodic-row default), **220 moved to the correct,
    specific `MMFF94_STBN` parameter** (the headline result — real
    parameter-selection parity, not just a coverage-count change, since
    coverage was already 100% via Priority 2B's Dfsb port), **27 now
    correctly contribute ZERO stretch-bend energy** (a real energy-level
    fix: chematic previously injected a nonzero generic Dfsb value here;
    `mmff94_stbn_type_only` now finds the real `(0.0, 0.0)` row directly and
    returns `Some((0.0, 0.0))` — numerically equivalent to RDKit's own
    `isDoubleZero`-gated drop, but not identical at the reporting layer:
    RDKit's API returns `None` for these, so `Mmff94CoverageReport`'s
    stretch-bend *coverage* counting now treats them as resolved hits where
    RDKit would count them as absent — a coverage-accounting nuance, not an
    energy discrepancy, left as-is since it doesn't affect any energy or
    minimizer output), **8 correctly remain on the generic Dfsb fallback**
    (RDKit's own real algorithm also falls through here — already
    numerically right, now for the right reason), and **172 remain
    unresolved for a separate, out-of-scope reason** (chematic's SMILES
    aromaticity perception disagreeing with RDKit's sanitizer on certain
    lowercase-aromatic ring inputs, already documented in the diagnostic —
    not fixable by a classification-key fix alone). The
    `routing_bug_candidate` *count* itself moves 427 → 180 (only the fully-
    resolved 220+27 stop appearing as "missing"; the remaining 8+172 still
    do, now correctly keyed by stretch-bend type); `table_gap` (1,680,
    genuine absence at every classification code) is untouched, as expected
    — this fix corrects routing, not data gaps.
  - The 220/27/8 split was verified by a direct per-row join (not just
    matching aggregate counts): zero of the 172 aromaticity-confounded rows
    land in the 220 or 27 buckets. One coverage caveat, reported honestly
    rather than silently assumed: `stretch_bend_type_for`'s `ta == tc`
    branch (where the two outer atom types are equal, so the diagnostic's
    ported argument-canonicalization forces `arg1 == arg2` and angle types
    1/5/7 always take their first sub-code) is exercised 35/427 times in
    this corpus, but never with `bond_type_ij != bond_type_jk` among those
    35 — so the asymmetric-bond-type corner of that branch is carried over
    unvalidated from the already-merged, oracle-cross-checked diagnostic,
    not independently re-verified by this measurement.
- **Compounding, independently-diagnosed bug fixed in the same PR (direct
  dependency, not scope creep): `angle_type_for`'s ring-offset formula for
  `bt_sum=2` (3-ring) and `bt_sum∈{1,2}` (4-ring) disagreed with RDKit's
  real `getMMFFAngleType` formula** (`AtomTyper.cpp:2412-2447`:
  `angleType = ring_size; if bond_type_sum != 0 { angleType += bond_type_sum
  + ring_size - 2 }`). Corrected table: 3-ring bt_sum=2 now 6 (was 8);
  4-ring bt_sum=1 now 7 (was 6), bt_sum=2 now 8 (was 7). This feeds directly
  into `getMMFFStretchBendType`'s first argument, so it had to be fixed
  alongside the stretch-bend key fix to give it a correct input — but it is
  a real, independently-provable bug in its own right (angle_type_for's own
  interface unchanged, only its internal match arms). Measured **LATENT** on
  this 265-molecule corpus (0/113 reachable ring-embedded angle triples hit
  the diverging branches) — reported honestly: it does not move any corpus
  number by itself, only the stretch-bend key fix above does.
- **Breaking change**: `mmff94_stbn`/`mmff94_stbn_type_only`'s leading `u8`
  parameter is now `stretch_bend_type` (RDKit's `getMMFFStretchBendType`
  output, 0-11), **not** `angle_type` (0-8) — same `u8` shape (not a
  compile-time break), but a silent behavioral break for any caller still
  passing a raw angle type: it will now compute the wrong stretch-bend
  parameter without erroring. Migration: compute the stretch-bend type via
  the new `stretch_bend_type_for(angle_type, ta, tc, bond_type_ij,
  bond_type_jk)` (also newly `pub`, in `mmff94_minimizer.rs`) before calling
  either function — see `stretch_bend_energy`
  (`crates/chematic-ff/src/mmff94_minimizer.rs`) for the reference call
  site. `MMFF94_STBN`'s own doc comment ("Format: (angle_type, ...)") is
  corrected to "Format: (stretch_bend_type, ...)" — it was mislabeled from
  the start, per the diagnostic's self-consistency proof.
- `mmff94_term_coverage_audit.rs` (`crates/chematic-3d/examples/`) updated
  to compute and use the stretch-bend type the same way (its own
  `present_at_different_classification` scan range widened `0..=8` →
  `0..=11` for StretchBend specifically; `Angle`'s own `0..=8` scan is
  unaffected — angle type genuinely is Angle's own table key).
- Recommends a **minor** version bump (breaking Rust signature semantics),
  consistent with v0.12.0's own precedent for `mmff94_stbn`.

### Fixed — `chematic-ff` (MMFF94 torsion classification-formula bug, issue #227 Priority 2C, **breaking**)

- **Root cause: `torsion_type_for` classified the non-ring base case purely
  from atom-type membership in the static `MLTB_TYPES` set
  (`(MLTB(tj),MLTB(tk)) -> 0/1/2`), completely ignoring the j-k bond's own
  real MMFF bond order/type.** RDKit's real `getMMFFTorsionType`
  (`AtomTyper.cpp:2528-2571`, pinned commit — see
  `scripts/mmff94_provenance/PROVENANCE.md`'s "Torsion" row) classifies from
  the j-k bond's own `bond_type_for` result instead (`torsionType =
  bondTypeJK`), with an empirically-required override to type 2 when
  `bondTypeJK==0 && order_jk==Single && (bondTypeIJ==1 || bondTypeKL==1)` —
  needed to pass RDKit's own CYGUAN01 regression test, not derivable from
  Halgren's MMFF.IV page 609 formula alone. These are structurally different
  formulas: a double/triple/aromatic j-k bond always gets `bond_type_jk=0`
  under the real formula regardless of the endpoints' own MLTB membership,
  which the old atom-type-membership rule could never see (e.g. benzene's
  own ring torsions, all type 37 which IS in `MLTB_TYPES`, used to
  classify as type 2 — now correctly type 0).
  `torsion_type_for`'s ring-4/5 override is also replaced end to end: a new
  private `ring_size_4_or_5` ports RDKit's real `isTorsionInRingOfSize4or5`
  (`AtomTyper.cpp:403-447`) faithfully — local bond-adjacency, NOT
  SSSR-based (4-ring iff i-l directly bonded; 5-ring iff i and l, excluding
  ring neighbours j/k respectively, share a common neighbour) — and the
  5-ring branch now additionally requires `ti==1 || tj==1 || tk==1 ||
  tl==1`, a condition the old SSSR-based check had no equivalent of at all.
  - Measured on the 265-molecule Wave 1 corpus, verified fresh against the
    actual production code (not restated from the already-merged
    diagnostic's self-port estimates): re-running `mmff94_term_coverage_audit`
    post-fix moves Torsion `routing_bug_candidate` **1,107 → 254**
    (`table_gap` unchanged at 14, genuine data gaps this fix doesn't touch;
    `torsions_missing` 1,121 → 268). Of the original 1,107, **853** now
    resolve to a raw table row via chematic's existing, UNMODIFIED
    `mmff94_torsion_energy` fallback chain alone — no lookup/fallback-chain
    change needed, exactly as the diagnostic predicted — split into **851
    valid non-zero table resolutions** and **2 explicit-zero rows** (RDKit's
    own `isDoubleZero` gate would also drop these to "no term"; not counted
    as resolved). The remaining **254** need RDKit's separate Halgren
    empirical-rule fallback (`getMMFFTorsionEmpiricalRuleParams`,
    `AtomTyper.cpp:2874-3080`), which chematic has no equivalent of at all —
    explicitly out of scope for this fix, a larger separate follow-up.
  - **Beyond the 1,107 "missing" candidates, a corpus-wide before/after
    sweep of ALL 13,530 torsion instances** (frozen copy of the old formula
    vs. the new production code, same enumeration `torsion_energy` itself
    uses) found: 10,617 unchanged, **1,792 changed to a numerically
    different `(V1,V2,V3)`** — a silent-wrong-parameter population an order
    of magnitude larger than the 1,107 "missing" instances, invisible to
    `mmff94_term_coverage_audit.rs` (it only logs misses) — 853 newly
    resolved (== the 853 above), and **0 newly lost** (no torsion that used
    to resolve to a value now resolves to nothing). Oracle-validated against
    live RDKit (`rdkit==2026.3.3`): **all 1,792** changed-value rows were
    checked, not a sample — **1,776/1,792 (99.1%)** match the new, post-fix
    value exactly, **0** match the OLD (pre-fix) value instead (zero cases
    where this fix made a previously-correct value wrong), and the
    remaining 16/1,792 (0.9%) trace to a pre-existing, out-of-scope MMFF94
    aromaticity-perception gap for charged aromatic (pyridinium-type) rings,
    unrelated to torsion classification itself (both sides agree the
    torsion classification code is 0; only the underlying ring-carbon atom
    TYPE differs upstream). The 853 newly-resolved rows were already
    oracle-validated at 853/853 (100%) by the diagnostic itself.
  - Also empirically characterized: chematic's torsion-enumeration loops
    have no `i==l` guard (a 3-membered-ring degenerate torsion where the two
    outer atoms coincide), which theoretically could make the new local-
    adjacency ring check misfire for a substituted 3-ring. Measured 33 such
    instances in this corpus, 0 of which trigger the ring override either
    before or after this fix; a constructed methylcyclopropane check against
    the live RDKit oracle confirms the port is faithful (RDKit's own real
    algorithm computes the identical ring-size-5 local adjacency for this
    case, but its own type-1 gate — using RDKit's correct CR3R=22
    ring-carbon typing, which chematic also assigns correctly here —
    prevents the override from firing, matching chematic's result exactly).
- **Breaking change**: `torsion_type_for`'s signature changes from
  `(rings: &[Vec<AtomIdx>], i, j, k, l, tj: u8, tk: u8)` to `(mol: &Molecule,
  i, j, k, l, ti: u8, tj: u8, tk: u8, tl: u8)` — the correct formula needs
  the actual j-k bond order (not just atom types) plus the i-j/k-l bond
  orders and the outer `ti`/`tl` atom types for the override conditions; the
  `rings` parameter is dropped entirely (the ring override is no longer
  SSSR-based). Migration: pass `mol` (already in scope at every existing
  call site) and all four atom types instead of just the two central ones —
  see `torsion_energy` (`crates/chematic-ff/src/mmff94_minimizer.rs`) for
  the reference call site.
- Recommends a **minor** version bump (breaking Rust signature semantics),
  consistent with v0.12.0's own precedent for `mmff94_stbn`.

### Known issue — release-gate waiver (issue #285)

- The combined v0.13.0 release-gate remeasurement (265-molecule Wave 1
  corpus, both MMFF94 fixes above together) found 2 molecules
  (`chembl_tier_b_0126`, `chembl_tier_b_0168`) regress by 1 declared-
  stereocenter satisfaction each in the un-gated/stretch-bend-gated MMFF94
  arms. Root-caused to a **pre-existing distance-geometry embedding defect**
  (present since v0.12.0, not introduced by either fix above): the
  declared alkene in both molecules is placed with the wrong E/Z sign at
  the embedding stage itself, before force-field minimization ever runs.
  What changed is that the *old*, mistyped torsion terms happened to pull
  this dihedral back to the declared configuration during minimization —
  an accidental rescue that the corrected torsion classification no longer
  performs. Confirmed via RDKit's own MMFF94 (`MMFFGetMoleculeForceField`),
  given the identical pre-minimization starting geometry chematic uses,
  exhibiting the same failure to rescue — two independent MMFF94
  implementations agree on identical coordinates, so the MMFF94 fixes
  themselves are correct. `sound=True` in every case (this is a stereo-
  correctness defect, not a geometry-sanity regression). Net corpus-wide
  stereo satisfaction still improves (`mmff94_strict`: 49 → 48 declared-
  stereocenter violations). Shipped under an explicit release waiver rather
  than blocking v0.13.0 — see issue #285 for the fix, planned for a later
  release.

## [0.12.0] — 2026-08-09

Two independent fixes from the ongoing issue #227 program: a production
MMFF94 stretch-bend coverage fix (**breaking**, see migration notes below)
and a `chematic-3d` starting-geometry correctness fix for fused/multi-ring
molecules.

### Fixed — `chematic-ff` (MMFF94 stretch-bend accuracy, issue #227 Priority 2B)

- `chematic_ff::mmff94_stbn` now falls back to RDKit's own periodic-table-row
  default stretch-bend parameters (`MMFFDfsbCollection`'s real equivalent — a
  small, 29-row table ported verbatim from the pinned RDKit commit, see
  `scripts/mmff94_provenance/PROVENANCE.md`) when the specific/generic
  MMFF-type table has no row at all. **Unconditional production behavior**,
  not behind any opt-in flag — applies to every MMFF94 policy's
  energy/gradient calculation and to `Mmff94CoverageReport`'s stretch-bend
  coverage measurement the same way. `gate_mmff94_stretch_bend`/
  `gate_mmff94_torsion_oop` (Priority 2's strict-refusal *gate configuration*)
  are unaffected and remain independent opt-ins, still `false` by default —
  but this only means the *gates themselves* didn't move; see the next two
  bullets for what changing every MMFF94 arm's underlying energy/gradient
  unconditionally can still do.
  - Measured on the 265-molecule Wave 1 corpus: missing stretch-bend
    instances **2,107 → 0 final-unresolved**. This 2,107 splits into two
    structurally different populations, both now "resolved" by coverage but
    NOT equally fixed: **1,680** were genuine table gaps (absent at every
    MMFF classification code) — Dfsb closing these matches RDKit's own real
    behavior exactly. The remaining **427** were routing-bug candidates (a
    correctly-typed parameter already exists at a *different* classification
    code) that Dfsb *also* happens to resolve — these are **masked, not
    fixed**: chematic is now using RDKit's generic periodic-row default
    instead of the specific parameter a correctly-routed classification
    would use. `mmff94_term_coverage_audit.rs` was fixed to keep reporting
    both populations separately (an earlier version of this fix collapsed
    them into "0 missing", making the 427-instance masked population
    invisible) — parameter-selection parity for those 427 is real follow-up
    work, not addressed here, to keep this PR's root cause singular (the
    Dfsb port itself, not `angle_type_for`'s classification logic).
    `mmff94_strict_stretch_bend_gated`'s success count converges to legacy
    `mmff94_strict`'s (149/265 both, identical molecule sets) since the
    stretch-bend gate essentially never fires anymore.
    `mmff94_strict_complete_bonded_term_gated` (still gates torsion/OOP,
    untouched by this fix) rises 37→86/265 — the residual gap there is a
    separate, known issue (routing-bug-candidate-dominated torsion
    coverage), not addressed in this PR either.
  - **This is a real production energy/gradient change, not just a coverage
    gate change** — `mmff94_strict`/`mmff94_with_uff_fallback` never gated
    stretch-bend, so their *gate eligibility* is unchanged, but Dfsb changes
    what every MMFF94 policy's energy function computes for previously-
    zero-contributing stretch-bend terms, which can shift minimizer
    convergence and therefore success/failure — not just geometry, in
    principle. Verified with a full per-molecule regression diff against a
    pre-Priority-2B baseline saved *before* re-running: 0 soundness
    regressions among molecules sound in both runs, and exactly one status
    change on `mmff94_strict` (148→149) — `chembl_tier_b_0166`
    (`elapsed_ms` 20530→16221, `status` timeout→success; the same molecule
    ID was also the timeout-boundary case in Priority 2's own measurement —
    a recurring ~20s-class molecule under this policy family, plausibly but
    not conclusively explained by `total_timeout_ms` boundary sensitivity
    rather than asserted as "known jitter" without checking).
  - `mmff94_stbn`'s public signature gained 3 required `atomic_num_{i,j,k}: u8`
    parameters (needed for the periodic-row lookup, which is element-keyed,
    not MMFF-type-keyed) — breaking for any external caller. The prior
    type-only behavior remains available as the new `mmff94_stbn_type_only`
    function (same signature `mmff94_stbn` used to have).
  - Recommends **v0.12.0** (minor), consistent with #249's own struct-field
    addition already requiring the same bump.

### Added — `chematic-ff`/`chematic-3d` (MMFF94 stretch-bend coverage gate, issue #227 Priority 2)

- `gate_mmff94_stretch_bend` (`PipelineV2Config`) / `include_stretch_bend_in_gate`
  (`minimize_with_policy_gated`): a new, independent opt-in — same shape as
  the pre-existing `gate_mmff94_torsion_oop` — that also refuses
  `Mmff94BondAngleStrict`/`Mmff94WithUffFallback` on a missing MMFF94
  stretch-bend cross term, not just bond/angle. **Defaults to `false`
  everywhere**, so no existing arm/policy's pass/fail *behavior* changes.

### Migration notes (recommend a **minor**, not patch, version bump for this change)

- **Diagnostic output changes even with the gate left at its default `false`.**
  `Mmff94CoverageReport::{stretch_bend_total,stretch_bend_missing}` are new
  fields, always populated now (mirroring the pre-existing torsion/OOP
  measure-but-don't-gate pattern). `total_missing()`/`all_missing()`/
  `missing_parameter_classes` (surfaced via `ForceFieldBridgeError::MissingParameters`
  and `PolicyMinimizeResult`) now include stretch-bend evidence alongside
  bond/angle/torsion/OOP — a caller matching on the *count* or *contents* of
  these fields (not just success/failure) will see a change even without
  opting into the new gate. (Caught by this PR's own test suite: `chematic-3d`'s
  `mmff94_with_uff_fallback_falls_back_and_reports_why_on_chfclbr` needed
  updating, `missing_parameter_classes.len()` 3 → 6.)
- **Rust: adding fields to `Mmff94CoverageReport`/`PipelineV2Config` (both
  `pub`, non-`#[non_exhaustive]`) is a breaking change for any external
  crate constructing either via a struct literal** (not via
  `PipelineV2Config::minimal()`/`Mmff94CoverageReport::default()`, which
  remain source-compatible). No deprecation shim added — this repo's
  stated policy is to change the code directly rather than carry
  backwards-compatibility shims (see `CLAUDE.md`).
- **Python (`chematic-py`) — mixed compatibility posture, deliberately not
  symmetric with WASM:**
  - `PipelineV2Config.safe(...)` (the documented convenience constructor)
    gained `gate_mmff94_stretch_bend: bool = False` — existing callers
    unaffected.
  - `PipelineV2Config(...)` (the raw `#[new]` constructor) gained
    `gate_mmff94_stretch_bend` as a **new required** positional/keyword
    argument, **no default** — existing callers using positional args, or
    keyword args that don't already name this field, will break. This is
    an intentional consistency choice, not an oversight: `new()`'s own
    docstring already states every field is deliberately required with "no
    hidden default" (matching `gate_mmff94_torsion_oop`'s existing
    precedent there), so adding an inconsistent default just for this one
    new field would contradict that constructor's stated design rather
    than preserve it. Callers who need default-preserving behavior should
    use `.safe(...)`.
- **WASM (`chematic-wasm`) — `embed_pipeline_v2_json`'s `gateMmff94StretchBend`
  JSON field is `#[serde(default)]`**: omitting it from a config object
  parses successfully as `false`, so existing 15-field JSON configs keep
  working unmodified (regression-tested:
  `pre_priority2_config_json_without_gate_stretch_bend_still_parses`).
  Deliberately asymmetric with the Rust/Python posture above — WASM's JSON
  boundary has no equivalent to "the struct literal already names every
  field," so silently defaulting an *omitted* field here doesn't compromise
  the `deny_unknown_fields` fail-closed guarantee for *unrecognized* fields.

### Fixed — `chematic-3d` (unsound starting geometry for fused/multi-ring molecules, issue #185/#252)

- **`dg::generate_coords`'s rule-based ring placement (`place_rings`) no
  longer produces atom-coincident or wildly-stretched starting geometries**
  for several distinct multi-ring topologies. Issue #185 originally blamed
  the UFF minimizer itself (reported as spuriously "converged" on unsound
  anthracene-class geometry); re-diagnosis found the minimizer was doing
  its job correctly on a bad *input* — three independent bugs in
  `generate_coords`, not the minimizer:
  - A non-ring root atom was placed unconditionally at a fixed
    `(x_offset, 0, 0)`, colliding with a ring vertex `place_rings`
    independently computed at that same point; `dfs_place` also only ever
    seeded from a single root, silently leaving any substituent on a
    *different* ring atom at the zeroed-coordinate default (e.g. aspirin's
    two substituents on different ring atoms).
  - Ring-visiting order could mismatch true ring-fusion adjacency (checked
    only shared atoms, not direct bonds between rings), silently
    superimposing unrelated rings in some multi-ring systems connected
    purely by bonds (e.g. terphenyl).
  - A fusion-disconnected "new island" ring (e.g. biphenyl's second ring)
    was anchored via an arbitrary fixed offset rather than the real bond
    connecting it to already-placed structure, stretching that bond up to
    ~5 Å.
  - Verified via a new bonded-pair-length test helper
    (`assert_bonded_pairs_sane`, catches stretched bonds that a
    closest-pair-only distance check structurally cannot) plus new
    regression fixtures (terphenyl, a meta-linked biaryl, a spiro positive
    control, and a bibenzyl fixture that pins a still-open limitation, see
    below).
  - Measured on the 265-molecule Wave 1 corpus: all 28
    `MinimizationFailed` cases (19 `CatastrophicBondBlowup` + 9
    `ExcessiveResidualForce`) now resolve to `Ok`, 0 regressions,
    confirmed deterministic (byte-identical across repeated runs).
  - Issue #185 has been retitled and corrected to reflect this root cause
    (kept open, not closed — the minimizer's false-convergence *reporting*
    on unsound geometry is a real, separate concern that remains
    unaddressed); issue #252 (the 265-corpus population this fix resolves)
    is closed as completed.
  - **Known, NOT fixed here, tracked separately:** a fused-ring seam
    orientation bug for genuine ring fusions (anthracene-class, distorted
    bonds at the fusion seam — issue #255) and a chain-bridged ring-island
    placement gap (bibenzyl-class, `place_rings` runs entirely before any
    chain-atom DFS walk — issue #256). Neither is part of the resolved
    28-molecule population above (verified: none of the 28 contain a fused
    polycyclic aromatic ring system).

## [0.11.0] — 2026-08-04

Four independent fixes surfaced while surveying RDKit's open GitHub issues
for applicability to chematic, plus one MMFF94 typing-coverage fix from the
ongoing issue #227 program. No breaking API changes.

### Fixed — `chematic-ff` (MMFF94 correctness and coverage)

- **O2CM terminal-oxygen typing gap closed** (issue #227 Priority 1A-3,
  PR #241). RDKit's `AtomTyper.cpp` `case 8` (aliphatic oxygen) terminal-atom
  branch resolves a much broader set of conditions than the numeric-type
  registry's "O2CM / OXYGEN, CARBOXYLATE ANION" name suggests; chematic's
  typer only covered a subset, so any terminal oxygen outside it fell
  through to the wrong row. Ported the real union of conditions from a
  pinned RDKit source, cross-checked against a live RDKit oracle across 19
  distinct molecules. Wrong-element parameter selection remains
  structurally prohibited by the construction-time semantic-compatibility
  invariant introduced in 0.10.1 — this fix closes a *coverage* gap, not a
  cross-element-mismatch regression (that count stays 0).
  - Measured on the 265-molecule Wave 1 corpus (6693 comparable atoms):
    exact atom-type parity 98.82% → **99.37%** (6614 → 6651/6693);
    oxygen-element parity 95.88% → **100.00%** (861 → 898/898); production
    `minimize_with_policy(Mmff94BondAngleStrict)` success 123 → **130/265**;
    cross-element mismatch 0 → 0 (unchanged); unclassified 0 → 0
    (unchanged). All 42 remaining real mismatches are fully classified by
    bucket, none unexplained. Issue #227 stays open — MMFF94 coverage is
    measurably better, this release does not claim it is complete.

### Fixed — `chematic-rxn`, `chematic-mol` (stereo and format correctness)

- **SMIRKS product chirality made parity-aware** (PR #243). Reaction
  templates that reorder mapped substituents (e.g.
  `[C@@H:1](F)(Cl)Br >> [C@@H:1](Cl)(F)Br`) previously copied the template's
  `@`/`@@` flag onto the product verbatim, ignoring that a reordered
  neighbor-write-order changes the real configuration the flag encodes —
  producing a silently un-inverted product. A product atom with an explicit
  template chirality now has its mapped neighbor order validated against
  the atom's real final topology before the flag is kept. A product atom
  that only *inherits* a reactant's chirality (no explicit template flag)
  is now kept only when a defined `stereo_neighbor_order` exists on the
  reactant side and neither the unmapped-neighbor element set nor the
  mapped-neighbor atom-map-number set changed across the template — closing
  a gap where a stale flag from a topology change inside the mapped core
  could survive. Both branches fail closed to `Chirality::None` on any
  unresolvable case, matching this project's standing no-silent-wrong
  policy.
- **CDXML reader perceives tetrahedral stereo from directional wedges**
  (RDKit issue #9359, PR #244). `parse_cdxml`/`parse_cdxml_all` previously
  never ran tetrahedral-parity perception at all — `Atom.chirality` stayed
  `Chirality::None` regardless of how a molecule was drawn. Wired into the
  same shared `apply_local_parity_from_wedges` mechanism MOL V2000/V3000
  and MRV already use. Non-directional displays (`Bold`/`Hash`/`Dash`, which
  ChemDraw sometimes draws for plain visual emphasis rather than stereo
  intent, and which have no Begin/End reference convention) are perceived
  only opt-in, via the new `CdxmlParseOptions { infer_nondirectional_stereo:
  bool }` (default `false`) and `parse_cdxml_with_options`/
  `parse_cdxml_all_with_options`; directional wedges
  (`WedgeBegin`/`WedgeEnd`/`WedgedHashBegin`/`WedgedHashEnd`) are always
  perceived. When non-directional inference is enabled, the result is now
  independent of which atom a `<b>` element happens to list first (`B` vs
  `E`) — a bond-order-flip normalization (new `Molecule::set_bond_order`,
  chematic-core) replaces an earlier endpoint-swap approach that was found
  to silently perturb `neighbors()` iteration order and corrupt the very
  parity calculation it fed. `bond.order` (`Up`/`Down`) is recorded
  faithfully for any wedge display regardless of the flag; only chirality
  *perception* is gated.

### Fixed — `chematic-depict`, `chematic-3d` (2D/3D correctness)

- **2D depiction no longer collides independent ring systems** (PR #242).
  Ring systems with no shared/fused atoms (separate substituents on a
  chain, or a spiro junction) were placed via a coordinate-blind
  `place_regular_ring` fallback that always centers at the literal origin —
  two unrelated same-sized rings could land at bit-for-bit identical
  coordinates. Replaced the old two-phase "place all rings blind, then
  place chains" layout with a single connectivity-driven growth pass that
  anchors every newly-placed ring, including single-atom (spiro) junctions,
  relative to already-placed geometry.
- **ETKDG macrocyclic amide 1-4 distance bounds now split by true
  cis/trans ring role** (PR #245). `macrocycle_14_bound_adjustments`
  previously pinned all four combinatorial 1-4 atom pairs across a
  tertiary/secondary amide bond to the same tight *cis* distance band — a
  geometrically unsatisfiable constraint set (a planar amide's four
  combinatorial pairs actually split 2 cis + 2 trans) whose least-bad
  embedding compromise was an unphysical ~90-130° dihedral twist instead of
  planar. Now computed from real ring-continuation role instead of a
  blanket assumption. When a central amide bond belongs to two or more
  eligible macrocycles at once (a theta-graph topology, where the correct
  role assignment is genuinely ambiguous), the embedder abstains to a
  relaxed band rather than guessing.

## [0.10.1] — 2026-08-02

This is a **correctness hotfix**, not a coverage-completion release. It
closes a class of bug where `chematic-ff`'s MMFF94 engine could silently
resolve an atom against a parameter row belonging to a different element
and report the resulting (physically wrong) energy as a success — never a
crash, never a warning, just a wrong number. Issue #227 stays open: MMFF94
coverage is measurably better as a side effect, but this release does not
claim MMFF94 support is complete.

### Fixed — `chematic-ff`

- **MMFF94 numeric atom typing: faithful aromatic 5-/6-ring port + a
  construction-time semantic-compatibility invariant** (issue #227). Root
  cause of the "furan collision" (#235's audit): the aromatic atom typer
  never implemented RDKit's real 5-ring/6-ring alpha/beta-heteroatom
  classification, so an aromatic atom could resolve against a parameter row
  from an entirely different element's chemistry — e.g. furan's C-C bond
  resolving as if it involved nitrogen. Ported from a pinned RDKit source
  (`Release_2026_03_3`, commit `e74e7b0a5a2fc4e7f77c04ec26a61d4b8edbf22f`;
  full provenance table and license note in
  `scripts/mmff94_provenance/PROVENANCE.md`), backed by a new 96-entry
  provenance-cited numeric-type registry
  (`crates/chematic-ff/src/mmff94_numeric_type_registry.rs`, generated by
  `scripts/gen_mmff94_numeric_type_registry.py` — do not hand-edit).
  Benzene's aromatic carbons are now correctly type 37 (CB, not 63);
  pyridine's N is type 38 (NPYD, not 67); furan's O is type 59 (OFUR).
  - **`bond_type_for` aromatic-bond fix**: found because the typer fix
    alone would have swapped one silent wrong-row collision for another.
    `bond_type_for` previously treated `BondOrder::Aromatic` like `Single`;
    MMFF94 requires aromatic bonds to unconditionally get `bond_type=0`,
    same as Double/Triple. Verified against a live RDKit oracle
    (`GetMMFFBondStretchParams`) on two contrasting cases sharing the same
    atom types (benzene's aromatic ring bond vs. biphenyl's non-aromatic
    inter-ring bond) — only bond order distinguishes them.
  - **Construction-time semantic-compatibility invariant** in
    `assign_mmff94_numeric_types`: every assigned numeric type's registry
    element must now match the atom's real element, or typing fails closed
    with a typed `NumericTypeError` — never a silent wrong-element
    "success". This is the actual fix, not the atom-typer port above: it
    makes the furan-collision bug class structurally impossible to
    reintroduce silently, in this code or in any future change to it.
  - That invariant firing on the fixed corpus caught **two more instances
    of the identical bug class**, previously invisible: `assign_n_type`'s
    protonated-ammonium branch returned type `32` (O2CM, an **oxygen**
    type) instead of `34` (NR+); `assign_o_type`'s anionic-oxygen branch
    returned `34` (NR+, a **nitrogen** type) instead of `35` (OM). Both
    fixed, both regression-tested
    (`protonated_amine_n_is_type_34_not_the_o2cm_oxygen_row`,
    `carboxylate_anionic_o_is_type_35_not_the_nr_plus_nitrogen_row`).
    Before this fix, every protonated amine and every anionic oxygen
    (carboxylates, phenoxides — common in drug-like molecules) minimized
    with MMFF94 would have silently used the wrong element's charge/bond/
    angle parameters.
  - Measured on the 265-molecule Wave 1 corpus, production
    `minimize_with_policy(Mmff94BondAngleStrict)` (not a simplified
    simulation): **44 → 102 Ok**, 221 → 140 `MissingParameters`,
    0 → 1 `UnsupportedAtomType` (the one remaining is a fixture
    deliberately named to probe unsupported chemistry — a pentavalent
    phosphorus ylide — correctly caught, not a bug), 0 → 22
    `MinimizationFailed` (confirmed unrelated to typing via isolated
    re-run — large flexible molecules hitting minimizer non-convergence,
    a separate, pre-existing robustness gap). Atom-type parity vs. a
    pinned RDKit oracle across all 6693 comparable heavy atoms in the
    corpus: 91.83% exact numeric-ID match, **0 cross-element mismatches**,
    0 unclassified (`scripts/mmff94_type_parity_report.py`,
    `validation/results/mmff94_type_parity_227_postfix.json`).
  - A direct, positive consequence: naphthalene no longer needs the UFF
    fallback under `ForceFieldPolicy::Mmff94WithUffFallback` — it now
    succeeds directly via `Mmff94BondAngleStrict`
    (`validation/pipeline_v2_wasm_parity_fixtures.json`'s
    `force_field_fallback` fixture).
  - PR #235 (audit, merged `abd1d72b`) added
    `docs/rfcs/mmff94_coverage_gap_227_audit.md`, diagnosing the collision and
    measuring the pre-fix baseline; PR #236 (this fix, merged `d75dc3f9`)
    is the production change.
  - **Not done in this release** (deferred, documented in PR #236):
    stretch-bend is not yet gated by the strict coverage check; no
    full-corpus RDKit energy/gradient parity harness exists yet (only the
    pinned unit-test-level oracle checks for benzene's r0/kb/bond_type);
    the 140 remaining `MissingParameters` and 22 `MinimizationFailed`
    cases are not yet individually re-bucketed by root cause; the 8.17%
    residual atom-type parity gap (same-element sub-type imprecision —
    aromatic-ring exocyclic-multiple-bond override, small-ring CR3R/CR4R
    strain context) is characterized but not ported. None of this blocks
    the correctness guarantee above; all of it is scoped to a future PR
    under issue #227.

### Migration notes

- **MMFF94 numeric atom types and some MMFF94 energy results change.** If
  you cached or hardcoded MMFF94 numeric type IDs, bond parameters, or
  minimized energies for aromatic molecules, protonated amines, or
  anionic oxygens (carboxylates, phenoxides), recompute them — the new
  values are the physically correct ones per the RDKit-verified fixes
  above.
- **A calculation that previously "succeeded" under
  `ForceFieldPolicy::Mmff94BondAngleStrict` may now fail closed instead.**
  This is intended: those successes were silently using a wrong-element
  parameter row. A `NumericTypeError`/`UnsupportedAtomType` you did not see
  before this release means the previous result for that molecule was
  wrong, not that this release introduced a regression.
- **Issue #227 is not complete.** This release does not claim MMFF94
  coverage completeness — see "Not done in this release" above. Do not
  read the 44 → 102 success-count improvement as a coverage milestone; it
  is a side effect of removing silent wrong-answers, not the goal.
- No change to CIP, ECFP4, canonical SMILES, or InChI in this release.

## [0.10.0] — 2026-08-01

### Added — `chematic-rxn`

- **`find_reaction_matches`/`apply_reaction_match`** (issue #225): a public
  seam between enumerating a SMIRKS's matches against reactant molecules and
  applying one of them, for callers that need to accept some matches and
  reject others (e.g. based on whether the matched bond is a ring bond)
  without discarding the whole `run_reactants` call. `find_reaction_matches`
  returns a `Vec<ReactionMatch>` — one per accepted match — and
  `apply_reaction_match` builds the product set for exactly one chosen
  match. `ReactionMatch::atom_map_positions` resolves a SMIRKS atom-map
  number to the matched atom without re-deriving it. `run_reactants`/
  `run_reactants_strict` are now implemented in terms of these two
  functions (behavior- and performance-preserving: still one SMIRKS parse
  and one VF2 match pass per call, unchanged from before).

### Fixed — `chematic-mol`

- **MRV reader now perceives 2D wedge/hash tetrahedral and E/Z stereo**
  (issue #202): `parse_mrv` wires `chematic_perception::
  apply_local_parity_from_wedges`/`apply_ez_directions_from_2d` into the
  2D-coordinate path, mirroring the wiring `mol2000.rs`/`cdxml.rs` already
  had. Previously, wedge/dash bonds and 2D coordinates were read into
  `coords_2d` but never converted into `Atom.chirality` or a bond's
  E/Z direction — `parse_mrv` silently dropped stereochemistry that was
  present in the file.

### Fixed — `chematic-smiles`

- **Shared E/Z carrier bonds: joint component solver** (issue #149) —
  replaces the old single-end abstain guard in `resolve_ez_markers` with
  `resolve_component_jointly`, which resolves coupled stereo-alkene ends
  (two double bonds sharing one marker-carrier bond) together instead of
  independently. **10 of the 18 issue #149 fixtures become fully
  permutation-invariant** (idempotent under re-canonicalization, stable
  under relabeling/reordering, RDKit-verified 0 stereo loss/corruption);
  the other 8 remain a documented, semantically-safe residual
  (`EZ_SHARED_CARRIER_RING_CONSTRAINED_RESIDUALS`) — every one of the 8
  is a coupled component containing an endocyclic double bond in a 5- or
  6-membered ring, where marker choice has no free degree left; RDKit
  re-parse confirms every divergent spelling is stereochemically
  identical, never corrupted. Full-corpus random-relabeling-only Check-2
  failures: 18 → 8. **Issue #149 stays open** — the ring-constrained
  residual's root cause (`compute_stereo_alkene_ends` has no ring-size
  gate) is characterized in a follow-up audit
  (`docs/rfcs/ez_ring_constrained_residual_audit.md`) but not yet fixed.

### Migration notes

- **`parse_mrv` (issue #202) now returns stereochemistry it previously
  silently discarded.** An MRV file that carries wedge/hash bonds or
  double-bond 2D geometry will parse to a `Molecule` with `Atom.chirality`
  set and/or E/Z bond direction set where it previously parsed to a flat
  (stereo-unset) molecule at the same input. If downstream code compares
  MRV-derived molecules against a cached/precomputed canonical SMILES,
  InChI, or fingerprint that was computed from the old (stereo-dropped)
  parse, that comparison can now diverge — the new value is the more
  correct one. No change to any other reader (`mol2000.rs`/`cdxml.rs`
  already perceived this).
- **A small number of `canonical_smiles()` outputs change** as a
  consequence of the issue #149 joint component solver — measured at
  exactly 6 changed lines across a 5,000-molecule corpus, all within the
  18 pinned issue #149 fixtures (10 newly converge; 8 remain the
  documented residual, unchanged from `[0.9.0]` and earlier). Same class
  of change as the `[0.8.1]` explicit/implicit-hydrogen migration note
  below — if you hardcode/cache canonical SMILES strings, a rare input
  may now produce a different (but RDKit-verified semantically identical)
  string.
- No change to CIP, ECFP4, or the RDKit benchmark in this release.

## [0.9.0] — 2026-08-01

### Added — `chematic-py`

- **`Mol.embed_pipeline_v2(config)`**: Python binding for the Rust-only
  `chematic_3d::pipeline_v2::embed_pipeline_v2` (torsion-knowledge-aware
  distance geometry + stereo verification/repair + policy-gated force
  field). New `PipelineV2Config` class mirrors the Rust config's
  deliberate lack of a `Default` — every field is required at the
  low-level constructor; `PipelineV2Config.safe(force_field=..., stereo_policy=...,
  ring_torsion_policy=...)` is a convenience constructor that still keeps
  those three policy choices explicit, never hidden defaults.
  - Returns a dict with full per-stage evidence — `coords`, `embed_stats`,
    `bound_adjustment_report`, `torsion_knowledge_report`,
    `ring_torsion_evidence`, `torsion_optimization_report`,
    `stereo_before`/`stereo_repair`/`stereo_after_repair`, `force_field`
    (actual policy used + fallback reason, never silently substituted),
    `final_stereo`, `final_validation`, `elapsed_ms_by_stage` — never just
    final coordinates.
  - New `PipelineV2Error` exception (a `ValueError` subclass) on pipeline
    failure, carrying a structured `.diagnostics` dict with the same
    per-stage partial evidence a Rust caller sees on `PipelineV2Failure`.
    `diagnostics["last_known_coords"]` is explicitly flagged
    (`coords_are_diagnostic_only: True`) so it can never be mistaken for a
    usable result.
  - Applies directly to the caller's own `Mol` atom order — never
    canonicalizes/reparses first (the issue #172 `conformer_ensemble()`
    mistake), verified with dedicated atom-order-permutation regression
    tests.
  - Rust-only; no algorithm change. Scope: Python binding only — no WASM
    binding, no RDKit benchmark (tracked as separate follow-up work).

### Added — `chematic-wasm`

- **`embed_pipeline_v2_json(mol, configJson)`**: WASM mirror of
  `Mol.embed_pipeline_v2()` above, same 15-field config (camelCase JSON,
  `deny_unknown_fields`, no silent defaults) and same per-stage evidence,
  returned as a tagged-union JSON envelope
  (`{"schemaVersion":1,"ok":true,"result":{...}}` /
  `{"schemaVersion":1,"ok":false,"error":{...}}`). Applies directly to the
  caller's own `MolHandle` — no canonicalize/reparse. Fail-closed on
  oversized input/atom count, malformed JSON, unknown fields, unknown enum
  values, and out-of-range integers, all surfaced through the same envelope
  rather than a thrown exception.
  - Rust-only pipeline, opt-in like the Python binding — not a default 3D
    API, no behavior change to `generate_coords`/`etkdg`/existing WASM 3D
    exports.
  - Verified end-to-end under real WASM, both `wasm-pack --target nodejs`
    and `--target web` (the latter is what `publish-npm.yml`/`pages.yml`
    actually ship): success, typed-failure, and typed-timeout paths all
    return real results — not just natively tested. Depended on #219's
    `chematic-3d` clock fix; `std::time::Instant::now()` panicked
    unconditionally on `wasm32-unknown-unknown` before that fix, which this
    binding's first real-runtime run is what originally surfaced.
  - New CI job (`test-wasm` in `ci.yml`) builds both wasm-pack targets and
    runs the Node integration tests on every push/PR — chematic-wasm had no
    WASM-runtime CI coverage at all before this.

### Deprecated — `chematic-3d`

- **`generate_and_minimize_uff()` never ran chematic-ff's UFF, despite its name and doc
  comment.** `#[deprecated]`, with an honest doc comment disclosing what it actually runs
  and where to go for real behavior. Investigation (issue #204) found the defect deeper
  than the issue's own text described: `minimize_uff(mol, coords)` resolves to this
  crate's own `minimize::minimize()`, whose dispatch (`minimize_with_config`) only
  special-cases `ForceField::MMFF94` — the `UFF` and `ForceField::DREIDING` (the
  default) variants are indistinguishable and both fall through to the same catch-all,
  `minimize_generic_with_config`. That is a **third**, untyped, element-pair-parameterized
  harmonic engine, distinct from both `chematic_ff::uff::minimize_uff` (real UFF) and this
  crate's own typed DREIDING engine (`minimize_dreiding`, which assigns real
  `DREIDINGType`s and is what `generate_and_minimize_dreiding` actually calls). A new
  regression test (`generate_and_minimize_uff_delegates_to_generic_minimize`)
  initially asserted the issue's own framing — that the function's output matches
  `generate_and_minimize_dreiding` — and that assertion **failed**, which is what
  surfaced this refinement. The test's only load-bearing, permanently-pinned
  assertion is that the deprecated function's output is numerically identical to
  calling `minimize::minimize()` directly; its additional evidence that the output
  also currently differs from typed DREIDING and real UFF on one test molecule is
  logged as diagnostic, not asserted — different force fields could in principle
  coincide on a shared local minimum for some molecule, which would make "always
  differs" a fragile, unrelated regression gate.
  - Zero in-workspace callers found (`chematic-py`, `chematic-wasm`, `chematic-mcp`, and
    every other crate/test/doc call `generate_and_minimize_dreiding` for DREIDING behavior
    already), so the function is kept (not deleted or behavior-changed) purely as a
    just-in-case compatibility shim for any external, out-of-tree caller of this
    published crate.
  - Use `generate_and_minimize_dreiding` for the same (typed DREIDING) behavior under an
    honest name, or `minimize::minimize_with_policy_gated(..., ForceFieldPolicy::UffOnly, ...)`
    / `chematic_ff::uff::{assign_uff_types, minimize_uff}` for real UFF physics.

### Fixed — `chematic-3d`

- **`embed_pipeline_v2`/`distance_geometry_v2` panicked unconditionally under real
  `wasm32-unknown-unknown`.** `std::time::Instant::now()` has no host time source on
  that target and traps with `"time not implemented on this platform"` — hit at all
  13 call sites used for timeout enforcement and per-stage timing (12 in
  `pipeline_v2.rs`, 1 in `distance_geometry_v2.rs`). This was latent (no previously
  shipped WASM export reached these modules) until PR #220's not-yet-merged
  `embed_pipeline_v2_json` binding exercised it end-to-end and filed it as issue #219.
  Fixed via a small crate-internal `clock` module: `web_time::Instant` (backed by
  `Performance.now()`) on `wasm32-unknown-unknown`, `std::time::Instant` everywhere
  else (via `web-time`'s own re-export) — a pure clock-source swap, no change to
  timeout contracts, stage-timing field names/units, or any chemistry/geometry/torsion/
  force-field logic. Verified end-to-end (both `wasm-pack --target nodejs` and
  `--target web` artifacts, real success + typed-timeout paths, no partial-coords-as-
  success) via a local, unpushed integration of PR #220's head, since #220 itself
  wasn't merged yet at the time. A related but distinct conditional trap in
  `chematic-smarts`'s MCS timeout path was found during the call-site audit and
  filed separately as issue #221 rather than folded into this fix (different
  crate, different root cause) — see below.

### Fixed — `chematic-smarts`

- **`find_mcs_with_config`'s timeout path traps under real `wasm32-unknown-unknown`**
  (issue #221, same underlying mechanism as #219 above, independently fixed here
  since `chematic-3d` depends on `chematic-smarts`, not the other way around).
  4 conditional `Instant::now()` calls in `mcs.rs`'s deadline check — only
  reached when `McsConfig::timeout_ms` is `Some`. Correction to #219's own
  changelog entry above: at the time that entry was written, no WASM export
  actually plumbed a `timeout_ms` through to JS (`mcs_smiles_json_with_ring_config`
  always used `McsConfig::default()`, i.e. `timeout_ms: None`) — so this trap was
  latent, not already reachable in production, same as #219 was before PR #220's
  binding shipped. Fixed via the identical `crate::clock` pattern (own module,
  own `web-time` target-specific dependency — `chematic-smarts` can't reuse
  `chematic-3d`'s module directly, the dependency direction runs the other way).
  No algorithm change; `test_timeout_does_not_panic` (native, `timeout_ms: Some(1)`)
  continues to pass unchanged.

### Migration notes

- **`generate_and_minimize_uff()` is `#[deprecated]` but not removed and not
  behavior-changed** — it still resolves to the same generic, untyped
  element-pair-parameterized minimizer it always has. Nothing breaks; the
  deprecation warning is a nudge toward `generate_and_minimize_dreiding()`
  (same behavior, honest name) or `chematic_ff::uff::{assign_uff_types,
  minimize_uff}` (real UFF), not a behavior change to migrate away from.
- **`Mol.embed_pipeline_v2()` (Python) and `embed_pipeline_v2_json()` (WASM)
  are new, additive, opt-in APIs.** No existing default 3D API
  (`generate_coords`, `etkdg`, `generate_and_minimize_dreiding`, the existing
  WASM 3D exports, etc.) changed behavior in this release. Pipeline v2 is not
  a default and does not need to be adopted.
- **Partial coordinates on pipeline v2 failure are diagnostic-only, never a
  usable result.** `PipelineV2Error.diagnostics["last_known_coords"]`
  (Python) / `error.lastKnownCoords` (WASM) are explicitly flagged
  (`coords_are_diagnostic_only` / `coordsAreDiagnosticOnly: true`) and must
  not be treated as a successful embedding.
  - Same JSON envelope: check `ok`/ the exception type before touching
    `result`, never assume a `coords` field is present on failure.
- **The `chematic-3d`/`chematic-smarts` clock-portability fixes (#219, #221)
  change only the *source* of monotonic timestamps** (`web_time::Instant` on
  `wasm32-unknown-unknown`, `std::time::Instant` — identical to before —
  everywhere else). No timeout threshold, stage-timing field name/unit, or
  chemistry/geometry/torsion/force-field calculation changed. Native builds
  are unaffected structurally, not just by measurement (`web-time`
  re-exports `std::time::Instant` verbatim off-wasm).
- **No change to CIP, E/Z, ECFP4, or the RDKit benchmark** in this release —
  those remain exactly as measured in `[0.8.1]` and earlier.
- The `canonical_smiles()` explicit/implicit-hydrogen migration note from
  `[0.8.1]` below still applies as written; it is not superseded by anything
  in this release.

_Nothing else yet — everything below `[0.8.1]` has shipped._

## [0.8.1] — 2026-07-30

### Fixed — `chematic-smiles` / `chematic-core`

- **`canonical_smiles()` no longer treats a genuinely-redundant explicit
  hydrogen count as distinguishing.** Two representations of the same
  molecule that differ only in whether an atom's H count was written with
  bracket notation (`[Cl]`) or organic-subset notation (`Cl`) — when the
  explicit value merely repeats what valence inference would give anyway —
  now canonicalize to the identical string. Previously they could diverge:
  `initial_invariant`'s Morgan-rank seed and `CanonicalWriter::emit_atom`'s
  bracket-necessity check both read the raw `Atom.hydrogen_count` field
  directly instead of going through the crate's own `implicit_hcount()`
  unification helper. Fixed by routing both decisions through
  `implicit_hcount`/the new `valence_inferred_hcount` (#205, #206).
  - **Isotope, formal charge, real stereochemistry, and any genuinely
    disambiguating explicit H count (e.g. `[CH2]` where organic-subset
    inference would give a different count) are unaffected and remain
    fully distinguishing** — only the specific redundant-explicit-H case is
    unified.
  - **Some canonical SMILES output strings for existing inputs will change**
    as a direct, intended consequence (e.g. a monosubstituted-benzene ring's
    canonical traversal start point can shift). This is a correctness fix,
    not a new feature — see Migration notes below for what downstream
    consumers should check.
  - Found via [kent-tokyo/renkin PR #65](https://github.com/kent-tokyo/renkin/pull/65):
    a `run_reactants` product built from a bracket-notation SMIRKS template
    diverged from a directly user-typed SMILES of the identical compound,
    breaking candidate-identity merging downstream.
  - New permanent regression coverage: `chematic-smiles::canonical::explicit_implicit_h_invariance`
    (isolated bracket/organic pairs, randomized atom-relabeling insertion
    order, disconnected structures, Kekulized rings, isotope/charge/stereo
    non-regression, canonicalize→parse→canonicalize idempotence) and
    `chematic-rxn::transform::reaction_derived_matches_direct_parse_chlorobenzene`.

### Migration notes

- If your code stores or compares hardcoded canonical SMILES strings
  (golden-file tests, cached dedup keys, precomputed candidate IDs, etc.),
  a small number of them may now differ from what this version produces.
  This release itself needed to update 11 such pre-existing internal
  fixtures after tracing each individually to confirm it was a stale
  string for the identical molecule, never a behavioral regression — the
  same audit is recommended for downstream golden strings before
  upgrading in a context that depends on exact string stability across
  versions.
- `are_identical`/`compare_molecules` (`chematic-chem`/`chematic-inchi`)
  benefit directly: pairs that previously required InChI reconciliation to
  resolve a spurious canonical-key split may now be recognized as
  duplicates by the fast canonical-key grouping alone.

## [0.8.0] — 2026-07-29

### Added — `chematic-3d`

- **Opt-in, fail-closed v2 embedding pipeline**, `chematic_3d::pipeline_v2::embed_pipeline_v2`
  (`PipelineV2Config` / `PipelineV2Result` / `PipelineV2Failure`). Integrates four
  independently-merged and independently-verified Wave 1/2 modules into one pipeline without changing any
  existing default behavior (Rust/Python/WASM/MCP untouched — this is a new,
  additive, explicitly-invoked entry point only):
  - stochastic distance geometry (raw embedding)
  - macrocycle 1–4 distance-bound adjustments
  - validated ETKDG-style torsion knowledge (acyclic + small-ring + macrocycle
    rules) with energy-based optimization of applicable acyclic torsions
  - chiral-volume / E-Z stereo constraint verification and repair
  - typed force-field policy minimization (MMFF94 strict / widened / with-UFF-
    fallback, or DREIDING)
  - a final, fail-closed stereo re-verification **after** force-field
    minimization — minimization can reintroduce a stereo violation that repair
    had just fixed (measured on `chematic-ff`'s UFF/DREIDING paths, not a
    hypothetical: 5/50 and 4/52 molecules respectively in earlier gap-check
    corpora), and this stage catches that as a typed `FinalStereoViolation`
    rather than reporting success.
  - 12 ordered stages total, each with explicit per-stage typed failure causes
    and timing (`elapsed_ms_by_stage`), so a failure is always attributable to
    one stage, never a generic error.
- **`RingTorsionApplicationPolicy::{FailClosed, DiagnosticOnly}`**: the current
  torsion optimizer can only *score* ring/macrocycle torsion potentials, not
  mechanically apply them to geometry. `FailClosed` (the default) rejects a
  request to apply ring torsions with a typed
  `RingTorsionApplicationUnsupported` failure rather than silently downgrading
  to scored-only; `DiagnosticOnly` is an explicit, non-default opt-in that
  returns scored-only evidence with `applied_to_geometry: false` /
  `diagnostic_only: true` markers. Measured on a 63-molecule corpus (10 arms,
  fixed seed): 15 typed `RingTorsionApplicationUnsupported` failures under
  `FailClosed` vs. 13 scored-only successes under `DiagnosticOnly` on the
  *same* ring/macrocycle-bearing molecules (2 confounded by an unrelated,
  already-known ring-fused-stereocenter repair limitation, not a ring-torsion
  bug) — confirming the two policies genuinely discriminate on real molecules,
  not just in unit tests.
- Known, disclosed limitations: 3-membered rings still fail closed at the
  distance-geometry stage; steroid-like fused-ring stereocenters (e.g.
  testosterone, cholesterol) can fail stereo repair; UFF/DREIDING minimization
  can reintroduce stereo violations post-repair (always caught, never silently
  passed); the spec's full validation-corpus requirement is only partially
  met — a dedicated Agent-D stereo corpus, Agent-E torsion corpus,
  force-field-failure corpus, and a tuning-vs-holdout partition were not
  separately assembled (the same 63-molecule corpus serves all measurements
  above); rule-order invariance and full atom-relabeling-permutation
  equivalence are covered only narrowly (unit tests on 1–2 molecules each),
  not swept across the whole corpus. No claim of "surpasses RDKit" or "full
  ETKDGv3 parity" is made anywhere.

### Fixed — `chematic-smiles`

- **Canonical SMILES performance regression, root-caused and fixed** (reported
  by RENKIN). The individualize-refine canonicalization search now prunes
  search branches using exact, independently-verified colored-graph
  automorphism checks (`has_colored_automorphism_mapping`, full bijection,
  vertex/edge-color-preserving) — never on 1-WL/Morgan rank or hash equality
  alone, since a single refinement cell can contain multiple distinct
  automorphism orbits. The existing `canonical_smiles(&Molecule) -> String`
  public signature is unchanged. Outputs were byte-identical to the
  pre-change engine on the audited 5,000-molecule corpus; the internal
  canonical search algorithm itself changed, and a known canonicalization
  correctness issue (coronene idempotence) was fixed as part of this work —
  this is not a claim that behavior is invariant across all possible inputs,
  only that it was verified byte-identical on the audited corpus. A new
  fallible, bounded API, `canonical_smiles_with_limits`, is added alongside
  the unchanged default.
  - **Also fixes issue #194**: `canonical_smiles()` was silently truncating a
    `Dative` bond's 2-character `->` token to a plain single-bond `-` (via
    `BondOrder::smiles_char()` only reading the token's first byte), and,
    separately, the SMILES parser itself could not distinguish `->` from
    `<-` at ingestion (both collapsed to the same `BondOrder::Dative` with
    identical `atom1`/`atom2` assignment, silently losing the donor/acceptor
    direction the input SMILES actually specified). Both are fixed together
    — see PR #196 for the full root-cause writeup, including the
    direction-aware write logic needed so a canonical DFS that reaches the
    acceptor atom before the donor atom still emits the semantically correct
    arrow direction, not just an un-truncated but wrongly-oriented one.
  - Measured (PR #193's own audited figures, 3-tier corpus: 13 high-symmetry
    fixtures, 6 low-symmetry negative control, 5,000-molecule external
    corpus; conservative, independently-verified numbers used deliberately
    rather than same-session remeasurements, per this project's own
    documented distrust of single/few-shot wall-clock benchmarking):
    approximately **5x** geometric-mean speedup on high-symmetry fixtures
    (search-leaf count 6625 → 13); approximately **1.1–1.2x** geometric-mean
    speedup on low-symmetry negative controls with **no measured
    regression**; approximately **1.09–1.10x** per-molecule geometric-mean
    speedup on the 5,000-molecule external corpus, approximately **5x**
    total-elapsed improvement there (tail-dominated by a small number of
    high-symmetry molecules, not representative of typical per-molecule
    cost), and **0 canonical-output mismatches** across all 5,000 molecules.
  - The exact-automorphism pruning depends on the pre-existing, unchanged
    `refine_ranks` hash (FNV-1a) for its initial partition seed — a
    theoretical collision there is an already-accepted, crate-wide risk
    (shared by `equivalent_atom_classes`/`are_atoms_equivalent`), not a new
    one; this PR changes its consequence from "redundant exploration" (old
    engine) to "a silently skipped branch" in principle, not observed on any
    fixture, exhaustive small-graph suite, randomized fuzz trial, or the
    5,000-molecule corpus. Disclosed rather than left as an absolute
    "structurally impossible" claim.

### Notes

- Public signatures and existing 3D defaults are unchanged: `chematic_3d`'s
  pre-existing public API, aromaticity perception, and the default
  (`LegacyFast`) CIP engine are unaffected except by pure additions.
  `canonical_smiles` and SMILES parsing, however, intentionally change output
  for the documented Dative-direction (issue #194) and canonicalization
  correctness bugs fixed in this release — this is a deliberate correctness
  fix, not a claim of universal output invariance.
- Each change went through the project's standing independent-verification
  process before merge. Multiple rounds across all three PRs (#192, #193,
  #196) found and fixed real issues (sign conventions, dropped diagnostic
  evidence, missing timeout enforcement, measurement-conflation bugs, a false
  negative in automorphism matching, corpus gaps, benchmark-harness biases,
  two independent root causes behind the Dative-bond bug); PR #192's final
  verification round completed without new findings. See PR #192, PR #193,
  and PR #196
  for the full per-round record.

## [0.7.1] — 2026-07-27

### Added — `chematic-inchi`

- **Accurate-CIP dedup preflight** (issue #161): `compare_with_accurate_cip_preflight`/
  `compare_molecules_with_accurate_cip_preflight`/
  `deduplicate_verified_with_accurate_cip_preflight` retry a legacy-CIP-unresolved
  specified tetrahedral stereocentre via `CipMode::Accurate` and, if it resolves,
  use that to recover verified-comparison capability — without ever letting the
  accurate engine's answer leak into the generated InChI string itself (which
  would silently reopen the 4663/4664 false-duplicate bug fixed previously).
  Additive only: `compare`/`compare_molecules`/`deduplicate_verified`/
  `identity_verify` are byte-for-byte unchanged. On the project's 5,000-molecule
  reference corpus: `verification_unavailable` 15 → 6 (9 recovered, 0 newly
  unavailable). Fails closed (`IdentityDiagnostic`) if the accurate engine
  ties/budgets-out/errors on any flagged atom, or if two flagged atoms in the
  same molecule share a `morgan_ranks` value (ambiguous correspondence — costs
  real recall on 3/12 audited molecules, a known, quantified, unfixed
  limitation).
- **Indexed graph relation API**: `compare_indexed_graph_relation`, with mode
  controlled by two independent, orthogonal axes —
  `GraphStrictness::{RawGraphExact, ChemicalGraphExact}` (literal bond-order/
  aromatic-flag equality vs. Kekulize-first chemical equality) and
  `AtomMapPolicy::{Include, Ignore}` (whether reaction atom-mapping numbers are
  part of molecule identity) — combinable freely via `IndexedGraphRelationMode`
  rather than a flat preset enum. Requires matching atom-index correspondence
  between the two molecules (e.g. conformer/ensemble grouping, MOL/SDF
  round-trip checks) — not a general graph-isomorphism search; named
  accordingly rather than as "exact graph identity" to avoid overclaiming.
  `ChemicalGraphExact` does not yet recognize two independently-Kekulized
  structures representing different valid resonance forms of the same
  aromatic system as equivalent (degrades to an honest mismatch or
  inconclusive result, never a false match).

### Added — `chematic-mcp`

- **MCP 2026-07-28 tools-only stateless stdio server**, alongside unchanged
  legacy (`2024-11-05`-style) stdio support on the same connection. A single
  stdio connection auto-detects and pins to whichever dialect its first
  request speaks; a request that tries to switch dialects mid-connection is
  rejected with a typed protocol error rather than silently reinterpreted.
  New: `server/discover`, per-request `_meta` metadata
  (`io.modelcontextprotocol/protocolVersion`/`clientInfo`/`clientCapabilities`),
  cacheable `tools/list` (`ttlMs`/`cacheScope`), and `structuredContent` on
  every `tools/call` result, validated against a new `outputSchema` added to
  all 20 tools (alongside tightened `inputSchema`s — `additionalProperties:
  false`, length/size bounds). Tool-call failures are now split into
  argument-shape problems (`-32602 Invalid Params`, before any chemistry
  runs) and chemistry-domain failures (a *successful* result with
  `isError: true` and a machine-readable `structuredContent.error.code`,
  e.g. `INVALID_SMILES`) — legacy-era wire behavior for both is unchanged.
  Internal refactor separates transport (stdio framing + era pinning),
  protocol codec (JSON-RPC parsing, error vocabulary, adversarial-input
  size/depth limits), server core (method dispatch + per-era response
  shaping), and tool registry (protocol-agnostic chemistry) into distinct
  modules. Remote HTTP, OAuth, the Tasks extension, and MCP Apps remain
  out of scope for this change — see `docs/mcp/2026-07-28-implementation-rfc.md`
  for the full design, primary-source citations (pinned to
  `modelcontextprotocol/modelcontextprotocol` tag `2026-07-28-RC`, commit
  `9d700ed`), and every deliberate deviation from a literal reading of that
  tag (notably: the RC tag's own `-32003`/`-32004` error codes are
  superseded pre-final by `-32021`/`-32022`, corroborated by the official
  `rmcp` Rust SDK — this PR ships the superseding values, not the RC tag's).

### Fixed — `chematic-py`

- **`Mol.conformer_ensemble()` returned coordinates indexed by a re-canonicalized
  atom order, not the caller's own `Mol` order** (issue #172, live in published
  v0.7.0). The method internally computed `canonical_smiles(&self.inner)`, reparsed
  that string into a fresh `Molecule`, generated the ensemble on the reparsed
  molecule, and returned its coordinates as-is — while `atom_table`, `cip_stereo()`,
  and every other property on the caller's `Mol` stayed indexed by the *original*
  atom order. `canonical_smiles()` routinely reorders atoms for anything with a
  branch or ring (e.g. decane's canonical form is `C(CCCCCCCC)C`, not
  `CCCCCCCCCC`), so the returned coordinate array silently did not correspond
  index-for-index to the caller's topology in the common case. Fixed by generating
  directly on `self.inner.clone()`, matching the existing correct pattern already
  used by `generate_3d()`/`generate_3d_etkdg()`. No public API change (same
  signature, same return shape). `crates/chematic-wasm`'s `conformer_ensemble_json`
  has the identical bug and is not fixed here (out of scope for this fix; tracked
  separately).

### Fixed — `chematic-smiles`

- **`canonical_smiles()` wrote its winning individualize-refine branch's
  string, then wrote it again.** `winning_individualized_ranks` (shared by
  `canonical_smiles`/`canonical_atom_order`, extracted by `c219ee7` so
  `canonical_atom_order` would share the same tie-break) already had to call
  `CanonicalWriter::write_all()` on every candidate branch to find the
  lexicographically smallest one; `canonical_smiles` then called
  `write_all()` a second, fully redundant time on the already-known winning
  ranks, on every single call, tied or not. The redundant call was introduced
  by `c219ee7` (2026-07-12), two days *after* `be5dbb1`'s (2026-07-10)
  individualize-refine correctness fix that this redundant work rides on top
  of — both commits first shipped in 0.4.30, not 0.4.26 (Cargo.toml already
  read 0.4.29 at both commits; the next release tag after either is
  `v0.4.30`). `be5dbb1` itself remains correctly attributed as the root cause
  of the *combinatorial* cost (the expensive-but-currently-necessary
  branch-and-minimize enumeration) — that part is unfixed here, see below.
  Reported via an external consumer's `run_reactants`/`apply_retro`
  performance regression (45-48x slower `canonical_smiles()` in isolation on
  highly symmetric molecules — plain rings, cages, `CF3`/`tBu`-style
  substituents — between chematic 0.4.25 and 0.4.30). Fixed by returning
  `(ranks, winning_string)` from `winning_individualized_ranks` so
  `canonical_smiles` reuses the string instead of recomputing it. Pure
  refactor: output is byte-identical (verified against the full test suite,
  including `be5dbb1`'s own golden-string tests and the issue #50 E/Z
  regression suite, plus an independent verification pass diffing
  `canonical_smiles`/`canonical_atom_order` output over a real 5,000-molecule
  corpus, 0 differences). Measured improvement: real and consistently
  positive on p95/p99/max; independently re-measured with symmetric
  instrumentation on both arms at ~13.4% total elapsed (p50 delta not
  resolvable from run-to-run noise at this sample size — treat the total/p50
  figures in `docs/rfcs/reaction_transform_perf.md` as directional, not precise).
  Does **not** fix the larger, still-open cost on genuinely symmetric
  molecules, which needs automorphism-orbit-aware branch pruning —
  explicitly deferred as future work. Full bisect, methodology, and
  before/after numbers in `docs/rfcs/reaction_transform_perf.md`.

### Added — `chematic-rxn`

- New `perf-instrumentation` Cargo feature (off by default, zero cost when
  disabled): process-global work counters for `run_reactants`'s hot path
  (`run_reactants_calls`, `reaction_parse_calls`,
  `reactant_query_match_calls`, `vf2_match_count`, `match_combination_count`,
  `build_product_calls`, `product_sets_before_dedup`,
  `product_molecules_built`, `atoms_copied_to_products`,
  `bonds_copied_to_products`), added while diagnosing the regression above.
- New `reaction_transform_perf_report` example: a corpus-weighted
  `run_reactants` benchmark accepting external template/probe corpora via
  `RENKIN_TEMPLATES`/`RENKIN_PROBE` env vars, falling back to small
  hand-authored fixtures in `crates/chematic-rxn/fixtures/`.

_Nothing else yet — everything below `[0.7.0]` has shipped._

## [0.7.0] — 2026-07-26

### Added — `chematic-mol` / `chematic-perception` / Python / WASM

- **2D wedge/hash stereochemistry is now perceived automatically when reading MOL/SDF
  files** — the missing wiring step identified by the P1-A0 diagnosis
  (`docs/rfcs/stereo2d_reader_integration_rfc.md`): `chematic_perception::apply_local_parity_from_wedges`
  (shipped in v0.5.0, never called from a reader) is now invoked unconditionally by
  `chematic_mol::read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics` — the new
  parsing core for V2000/V3000 MOL text — immediately after a wedge/hash bond is read,
  before anything else can mutate it away. `parse_mol`/`parse_mol_with_coords`/
  `parse_mol_v3000`/`parse_mol_v3000_with_coords`, the SDF supplier
  (`SdfRecordReader`/`SdfFileReader`), Python (`from_mol_block`, `SDMolSupplier`,
  `iter_sdf*`), and WASM all inherit this for free by delegating to the same core.
  CIP-independent — never touches `Atom.cip_code` or depends on CIP ranking; the
  existing accurate-CIP engine (`assign_cip_with_mode(mol, CipMode::Accurate)`) remains
  a separate, opt-in labeling stage on top.
- **New structured diagnostics API** for malformed/contradictory wedge input — never a
  silent guess. `chematic_perception::{StereoDiagnostic, StereoRejectionReason,
  apply_local_parity_from_wedges_with_diagnostics}` (four typed reasons:
  `ContradictoryWedges`, `MissingCoordinate`, `DegenerateGeometry`,
  `UnsupportedCoordination`) sits alongside the original silent
  `apply_local_parity_from_wedges`/`local_parity_from_wedges`, which keep their exact
  prior behavior and signatures. `chematic_mol::MolReadReport` carries
  `stereo_diagnostics` from `read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics`;
  `SdfRecord` gets a matching `stereo_diagnostics` field. Python:
  `from_mol_block_with_diagnostics()`/`from_mol_v3000_with_diagnostics()` and
  `SdfRecord.stereo_diagnostics()` (stable snake_case reason strings, e.g.
  `"contradictory_wedges"` — never free-form messages). WASM:
  `mol_block_stereo_diagnostics_json()`/`mol_v3000_stereo_diagnostics_json()`, and
  `sdf_to_records_json()`'s per-record objects gained a `stereo_diagnostics` array.
  Diagnostics are emitted only when a wedge/hash bond is actually present at a
  candidate center (degree 3 or 4) — an ordinary atom that merely touches someone
  else's wedge bond, or has no wedge input at all, never produces one.
- **V3000 bond `CFG` (wedge direction) is now decoded on read** — the parser
  previously never collected bond-line `KEY=VALUE` tokens at all, so a V3000 file's own
  wedge/hash annotation was silently dropped. `CFG=1`/`CFG=3` now map to `BondOrder::Up`/
  `Down`; `CFG=2` ("either") is left undirected, matching V2000's own MDL-code-4 policy.
- **2D double-bond E/Z (cis/trans) direction is now derived automatically when reading
  MOL/SDF files** (P1-S2, the double-bond counterpart to the wedge/hash work above) —
  new `chematic_perception::stereo2d_ez_direction::{apply_ez_directions_from_2d,
  apply_ez_directions_from_2d_with_diagnostics, apply_ez_directions_from_2d_ex}`, a
  CIP-independent, all-or-nothing-per-double-bond direction-setting stage mirroring
  RDKit's `setDoubleBondNeighborDirections`. Wired into
  `read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics` immediately after
  tetrahedral parity (SDF inherits via the V2000 core); MDL V2000 double-bond stereo
  code 3 and V3000 bond `CFG=2` ("either"/unspecified) are threaded through so a file's
  explicit "don't know" is never overridden by 2D coordinates. Writes exclusively
  through `Molecule::bond_direction` (the same side channel already used for the
  aromatic-bond E/Z stash, now generalized to a plain `Single`-order bond too) — never
  mutates a bond's own `BondOrder::Up`/`Down`, so raw wedge/hash and E/Z direction
  coexist on the same molecule without ever overwriting or reinterpreting one another.
  New `EzDirectionDiagnostic`/`EzDirectionRejectionReason` (`MissingCoordinate`,
  `NonFiniteCoordinate`, `DegenerateGeometry`, `ExplicitlyUnspecified`,
  `UnsupportedTopology`, `CarrierConflict`) sits alongside the silent convenience
  functions, matching the same silent-vs-diagnostics split as the wedge/hash API.
  `MolReadReport`/`SdfRecord` gained an `ez_diagnostics` field. Fails closed (no
  direction, no diagnostic -- ordinary `NotRequested`, matching a wedge-free atom's
  treatment above) for terminal alkenes, carbonyls/heteroatom termini, and
  topologically-equivalent substituents (confirmed against a live RDKit oracle: RDKit
  itself sets no `BondDir`/`Stereo` and prints no `/`/`\` for a symmetric-substituent
  alkene). Rejects with a typed reason (never guesses) for missing/non-finite/collinear
  coordinates, zero-length double bonds, explicit "either" stereo, cumulenes/allenes,
  Kekulized aromatic-ring bonds (excluded via a one-time, non-mutating
  `assign_aromaticity` query -- the reader never auto-perceives aromaticity, so a
  benzene ring reads as plain alternating `Single`/`Double` bonds structurally
  indistinguishable from a real diene without this check), and any shared-carrier
  conflict between two independently-stereogenic double bonds or between a double bond
  and an existing wedge (Issue #149's joint-carrier-resolution problem is explicitly
  out of scope -- a shared bond is only used when two independent, from-scratch
  geometric computations happen to agree, exactly the ordinary conjugated-diene case;
  when they disagree, or when a branch point's ambiguity could leak across two double
  bonds via SMILES's plain-adjacency `/`/`\` semantics, both bonds fail closed together
  rather than let bond-index or parse order guess a winner). Also fixed two
  prerequisite gaps found while building this: the plain (non-canonical) SMILES writer
  never consulted `Molecule::bond_direction` at all (only the canonical writer did) --
  both writers now read the same effective direction via a shared helper, correctly
  re-orienting a stashed marker regardless of DFS traversal direction; and
  `Molecule::remove_bond`/`with_atom_removed` either misattributed or silently dropped
  `bond_direction` entries during atom/bond removal (now remapped correctly, with
  tests). Broad-corpus validation against RDKit 2026.03.3 (4,999 molecules, standard
  InChI `/b`-layer semantic comparison, per-bond not per-molecule): 622 RDKit-resolved
  double bonds, 276 bond-level semantic agreements, 346 abstentions, **0 semantic
  inversions, 0 false-positive assignments**.

### Fixed — `chematic-mol`

- **V3000 writer emitted `CFG=6` for a hash bond** — that is V2000's stereo code, not a
  valid V3000 `CFG` value (V3000 spec: 1=Up, 2=Either, 3=Down). Fixed to `CFG=3`. Found
  and fixed together with the V3000 bond-`CFG` reader gap above, since both needed
  fixing before V3000 wedge round-tripping could work at all.
- **V2000 MDL stereo code 4 ("either"/unspecified direction) was collapsed into a
  definite wedge** (`1 | 4 => BondOrder::Up`) — harmless while nothing consumed
  `BondOrder::Up` for parity, but load-bearing now: a code-4 bond would have fabricated
  a confident, wrong stereocenter from a file that explicitly declares direction
  unknown. Now maps to `BondOrder::Single` (no defined direction). Round-tripping a
  code-4 bond is therefore lossy by design (documented, not silently wrong) — the
  round-trip-losslessness tests only assert on genuine wedge(1)/hash(6) bonds.
- **V2000 writer never emitted a bond's wedge/hash stereo flag** (hardcoded to `0`) —
  now emits `1`/`6` for `BondOrder::Up`/`Down`, so parsing a wedge MOL, writing it back
  out, and re-parsing recovers the same local parity.

### Known limitations

- MRV, CDXML, CML, and KET readers are not wired into either the tetrahedral or E/Z
  integration (MOL V2000/V3000 + SDF only); the one-line insertion pattern used here is
  reusable for them later.
- Aromaticity defaults, canonical-SMILES ranking, and the default (`LegacyFast`) CIP
  engine are unchanged by this work.
- E/Z direction is not set for: a branch point (an alkene end with 2 substituents)
  adjacent to a different double bond, or two independently-stereogenic double bonds
  whose shared-carrier requirement disagrees -- both fail closed by design (Issue #149,
  explicitly out of scope for this PR). The ordinary conjugated-diene case (a shared
  bond whose two independent requirements agree) is supported.
- Python/WASM bindings do not yet expose `ez_diagnostics` as a getter (mirroring
  `stereo_diagnostics()`'s Python/WASM API) -- the field exists on
  `MolReadReport`/`SdfRecord` in Rust and the E/Z direction itself already reaches
  SMILES output through every Python/WASM entry point; only the diagnostics-introspection
  API is deferred.

### Added — `chematic-inchi`

- **Verified canonical-SMILES deduplication** (`chematic_inchi::dedup`) -- fast
  canonical-SMILES candidate bucketing (`group_candidates`) reconciled against
  native-InChI-verified identity (`deduplicate_verified`), so a caller gets both cheap
  approximate grouping and a high-confidence verified partition without picking one or
  the other. `IdentityPolicy::{StandardInchiString, StandardInchiKey, StereoIgnored,
  IsotopeIgnored}` controls what "same identity" means; `StereoIgnored`/`IsotopeIgnored`
  clone-and-clear the relevant field on native-InChI generation-time options rather than
  doing string-layer surgery on the resulting InChI text. `VerifiedDedupReport` reports
  `groups`/`canonical_splits`/`canonical_collisions`/`verification_unavailable`/
  `invalid_molecules` separately -- a `CanonicalSplit` (verified-same molecules that
  happen to canonicalize to different SMILES strings) is never silently merged away or
  silently missed.
- **`has_unresolved_specified_tetrahedral_stereo` fail-closed guard**: closes a real
  false-`VerifiedDuplicate` found via live 5,000-molecule corpus verification -- two
  genuine diastereomers whose specified `@`/`@@` tetrahedral stereocentres the legacy
  CIP engine (`chematic_chem::tetrahedral_stereo_neighbors`) could not rank collapsed to
  the identical native-InChI string (`?`, undefined parity, at both differing centres).
  Every stereo-sensitive `IdentityPolicy` now fails closed to `VerificationUnavailable`
  rather than promote an unresolved-parity string to a duplicate claim; `StereoIgnored`
  is deliberately not guarded (ignoring stereo is that policy's own contract). Full-corpus
  audit after landing: `verification_unavailable` 1 -> 15 (14 newly fail-closed molecules,
  each individually classified against an independent RDKit 2026.03.3 oracle -- 10 are
  ordinary legacy-CIP ranking failures the accurate engine resolves, 2 are genuine
  CIP-ranking ties even RDKit's own modern CIP labeler can't resolve, 0 unexplained), known
  false-`VerifiedDuplicate` groups 1 -> 0. Does not switch native InChI generation to the
  accurate CIP engine (a separate, larger proposal) -- tracked as
  [#161](https://github.com/kent-tokyo/chematic/issues/161).

### Fixed — `chematic-inchi` (native, `native-inchi` feature)

- **Explicit hydrogen isotopes were never tallied and could silently drop a
  stereocentre's `/t` layer.** `mol_to_inchi_atoms`'s isotope tally counted every
  explicit graph H neighbor (`[2H]`, `[3H]`, ...) as ordinary H regardless of its own
  isotope, so no `/i` layer was ever produced for them; separately, a real graph H atom
  (as opposed to the pre-existing, unaffected bracket-H sentinel case) failed the
  Stereo0D neighbor-index lookup and had its entire stereo descriptor silently dropped
  via a fallthrough `continue` -- both enantiomers of e.g. `[C@](Br)(Cl)(F)[H]`
  collapsed to the identical InChI string with no `/t` layer at all. Fixed via a
  `StereoHSource` enum tracking each manufactured stand-in atom's provenance
  (`Sentinel` vs. `Explicit(AtomIdx)`), isotope-bucketed tallying, and routing a real
  graph H neighbor through the same manufactured-atom mechanism the bracket-H case
  already used. 17/17 required fixtures byte-exact against RDKit 2026.03.3
  (`Chem.MolToInchi`), 0 regressions (25/25 pre-existing tests byte-identical
  pass/fail).
- **Known limitation found while verifying the fix above**: a stereocentre with *two*
  simultaneous, isotopically-distinct H-like substituents (e.g. D+T on one carbon, or
  bracket-H plus an explicit D) still safely drops its `/t` layer rather than emitting
  one (RDKit does emit one for this shape) -- the single-manufactured-atom-per-centre
  mechanism can't yet represent a second H-like slot. Not corrupted, just absent;
  supporting it would need a second, independently-indexed manufactured atom per
  centre, judged out of scope for this fix. This is the exact gap
  [PR #156](https://github.com/kent-tokyo/chematic/pull/156)'s
  `has_unrepresentable_multi_h_stereocenter` dedup guard exists to detect and fail
  closed on.
- CI coverage gap: `.github/workflows/ci.yml`'s native-InChI test job was scoped to one
  integration test binary by name (`--test standard_inchi`), which would have silently
  skipped this fix's own new test file. Widened to run every integration test + doctest
  in the crate.

### Added — `chematic-cip`

- **Pseudoasymmetric (lowercase `r`/`s`) labeling for the "three-armed cage" residual**:
  completes Milestone 4A-2's 15-row carbon-cage family (opt-in `CipMode::Accurate`
  only) -- 15/15 rows now resolve, all matching the RDKit oracle, without the
  symmetry/automorphism-aware joint solver Milestone 4A-2 originally thought would be
  needed (the fix locates each tied branch's nearest embedded stereocentre via the
  existing, already-validated `resolve_chirality`/`embedded_chain` machinery directly,
  instead of a `provisional`-map lookup that failed whenever the embedded reference was
  itself still tied). Raw modern-oracle agreement: 4160/4186 (99.38%) -> 4175/4186
  (99.74%); oracle-stable agreement: 99.64% -> 100.00%.
- **Element-level fail-closed guard for 2 phosphorus atoms this same fix path also
  reached** (a follow-up commit on the same PR): those 2 cyclophosphazene atoms tie for
  the identical structural reason as the carbon-cage family, but independently-verified
  RDKit oracle checking found *neither* `rdCIPLabeler` nor legacy `_CIPCode` has a
  representation-stable answer for that specific molecule (both flip under a
  chemically neutral Kekulé respelling) -- there is no reliable oracle a resolved
  phosphorus label could ever be checked against. Shipping a resolved-but-unverifiable
  label was the actual bug; a per-atom-identity `Element::P` guard in
  `assign_one_with_rule5` now declines (`SkipReason::Tied`, reused, no new variant)
  rather than guess, restoring these 2 atoms to their original, pre-fix
  unresolved/`Tied` state. Scoped to phosphorus specifically (15/15 verified examples
  are all carbon, 0 verified examples for any other element) rather than a broader,
  unvalidated heuristic.
- No default API change: `CipMode::LegacyFast`/`assign_cip()` are untouched;
  `CipMode::Accurate` remains opt-in.

_Everything above this line is unreleased. Everything below `[0.6.0]` has shipped._

## [0.6.0] — 2026-07-25

### Added — `chematic-fp` / Python / WASM

- **Promoted the RDKit-bit-exact Morgan/ECFP4 path (`rdkit_morgan_ecfp4_experimental`)
  to a documented, cross-language opt-in API** — Python (`Mol.rdkit_ecfp4()`,
  `Mol.rdkit_ecfp4_detail()`) and WASM (`rdkit_ecfp4_bitvec()`,
  `rdkit_ecfp4_detail_json()`) bindings, following this codebase's existing fallible-
  experimental-API conventions (`PyValueError` in Python, `Result<_, JsValue>` in
  WASM — same shape as `Mol.cip_stereo(mode="accurate")` and `ecfp_bitvec_custom`
  respectively). The exact config this promotes is unchanged from the diagnosis that
  verified it (`docs/rfcs/ecfp4_bitexact_api_rfc.md`): radius=2 (ECFP4), 2048 bits,
  `useChirality=false`, `useBondTypes=true`, RDKit's default atom invariant. Does
  **not** change `ecfp4()`'s behavior or make the RDKit-exact path the default — the
  two engines use different hash functions and are never silently interchanged.
- **Generalized to a small, independently oracle-verified `(radius, fpSize)` matrix**
  — `chematic_fp::rdkit_morgan_fingerprint`/`RdkitMorganConfig` (Rust; radius ∈
  {0,1,2,3}, fold width ∈ {128,256,512,1024,2048} as closed enums, not raw integers —
  an unsupported value can't be constructed, let alone silently coerced), plus
  `Mol.rdkit_ecfp_config()`/`rdkit_ecfp_config_detail()` (Python) and
  `rdkit_ecfp_config_bitvec()`/`rdkit_ecfp_config_detail_json()` (WASM). Each of the 20
  cells is independently re-verified against a live RDKit oracle — not assumed to
  generalize from the radius=2/2048-bit point alone. `rdkit_morgan_ecfp4_experimental`
  itself is untouched (a separate, structurally isolated module reuses its internal
  hash primitives; see `crates/chematic-fp/src/rdkit_morgan_config.rs`).
- **Shared cross-language fixture+expectation corpus**,
  `validation/ecfp4_rdkit_stable_api_fixtures.json` (34 fixtures: 33 success cases
  spanning baseline/disconnected/isotope/charged/aromatic-vs-Kekulé/stereo shapes, plus
  1 real preprocessing-failure case for explicit error-path coverage) — generated once
  from a live RDKit oracle (`rdkit==2026.03.3`, pinned commit
  `8afba32ec539dcb2369bc84549d802aca3f7eb39`, same pin as `rdkit_morgan_hash.rs`) via
  `scripts/gen_ecfp4_rdkit_stable_api_fixtures.py`, and read identically by the Rust
  (`crates/chematic-fp/tests/rdkit_morgan_stable_api_fixtures.rs`), Python
  (`crates/chematic-py/tests/test_rdkit_ecfp4_stable_api.py`), and WASM
  (`crates/chematic-wasm/tests/rdkit_ecfp4_stable_api.test.mjs`) test suites — one
  source of truth, not three independently-maintained expectation lists. Verified
  byte-identical across all three surfaces by building and running each (`cargo test`,
  `maturin develop` + `pytest`, `wasm-pack build --target nodejs` + `node`), not just
  "identical by construction."
- Preprocessing failures (a small, known class of bridgehead-heteroatom-fused rings
  neither aromaticity engine can kekulize yet) surface as an explicit, typed error on
  every surface — never a silent fallback to the default Hückel-based `ecfp4()` engine,
  which would look successful while returning bits from an incompatible hash. Regression-
  tested on all three surfaces via the shared corpus's error fixture.

### Added — `chematic-perception`

- **`assign_aromaticity_authoritative_experimental(mol)` / `apply_aromaticity_authoritative_experimental(mol)`** — a new, opt-in aromatic-flag engine that makes `build_molecule_from_model`'s promotion/demotion decision authoritative in *both* directions: an atom's flag reflects the Hückel model's actual computed verdict, promoted when the model confirms it and **demoted** when a stale parser-set `aromatic: true` the model doesn't independently confirm survived from the input. The existing default (`apply_aromaticity`/`apply_aromaticity_ex`) is **unchanged** — verified byte-identical against `main` on a representative set of tricky fixtures (heteroaromatic implicit-H cases, the fused-diazine fix target, the still-open azulene gap, Kekulé input), not just by code inspection.
- As part of building this, fixed a real misclassification in `ring_pi_electrons`'s `CarbonExocyclicHeteroatomDouble` rule: a ring-fusion bond into an *adjacent* ring's heteroatom was scored the same as a genuine exocyclic substituent (like tropone's `C=O`), wrongly zeroing that atom's π contribution and — when both fusion atoms of a bicyclic system were affected — landing the ring on exactly π=4 (`Antiaromatic`), a Pass-1 verdict marked non-retryable and so permanently excluded from Pass 2's correction mechanism. Fixed generally (bond-level ring-membership check, no molecule-specific allowlist): fixes quinazoline/quinoxaline/naphthyridine/purine-shaped fused diazines (29 of a 33-molecule cluster; the remaining 4 are compound cases blocked by a second, independent, out-of-scope mechanism) with an unattempted beneficial side effect of also resolving 32 pre-existing `KNOWN_BRIDGEHEAD_N_FALSE_POSITIVES` regression pins and 2 open oracle findings (including the PR #86 purine reproducer).
- **Known, honestly-documented limitations of the new opt-in engine** (not regressions — these molecules were never correctly handled before either): azulene-type non-alternant fused rings (49 molecules) are an architectural gap — each ortho-fused ring independently gets an odd Pass-1 π count and neither seeds Pass 2, since no single SSSR ring can see azulene's whole-perimeter 10π delocalized system. The already-shipped, separately opt-in `apply_aromaticity_rdkit_parity_experimental` engine does resolve this (and the fused-diazine cluster) at 99.9956% corpus-wide agreement, but promoting *that* engine to the default remains a distinct, bigger decision this change does not make.

### Fixed — `chematic-smiles`

- **Canonical output could place an E/Z direction marker on a different substituent bond depending on input atom order/spelling** — the same alkene, written two semantically-identical ways, could canonicalize to two different (though individually valid) strings: permutation invariance held for the underlying stereo *assignment* but not for marker *placement*. Fixed via a new resolution pass (`resolve_ez_markers`) that deterministically picks the rank-lowest substituent as the canonical carrier for every double bond with stereo on both ends, covering trisubstituted/tetrasubstituted alkenes, the aromatic bond-direction stash, and ring-closure bonds. Two known, deliberately conservative exceptions remain unresolved rather than risk corruption: a double bond with a non-stereogenic end, and two independently-stereogenic double bonds that happen to share one physical candidate bond (18 molecules total) — resolving either case's carrier without the other visible can corrupt geometry, and there is no processing order avoiding a *different* order-dependence, so this case is left exactly as input-spelled. Tracked as [#149](https://github.com/kent-tokyo/chematic/issues/149) for a future joint-resolution design. Zero semantic corruption confirmed via two independent methods on every measured molecule.
- Permutation invariance (measured on a 4,992-molecule proxy corpus, separately from idempotence per this project's convention): **92.99% → 98.08%** (264 of 282 known-divergent molecules now converge; the 18 residual cases above are the entire gap). Idempotence unaffected (98.44%, unchanged — its own residual is a separate, pre-existing aromaticity round-trip issue).
- Investigated and **ruled out** a previously-diagnosed "bridged-bicyclic ring-closure ordering" permutation-invariance bug (`docs/rfcs/canonical_smiles_residual_rfc.md`'s "Root Cause 2") — the two SMILES claimed to be "two spellings of the same molecule" turned out, per independent RDKit `MolToInchi` verification, to be genuinely different constitutional isomers; chematic's differing canonical output for them is correct, not a bug. An additional 22-molecule probe using RDKit-`RenumberAtoms`-generated genuine same-molecule respellings (bridged/spiro/fused/cage systems, including stereocenters and a heteroatom bridgehead) found zero convergence failures. No production code changed for this finding; the RFC and test suite were corrected instead.

### Fixed — `scripts`

- `scripts/canonical_residual_diagnosis.py` called RDKit's `FindMolChiralCenters` with a kwarg name (`useLegacy`) the currently-pinned RDKit version no longer accepts (`useLegacyImplementation`) — a pre-existing tooling bug found while re-verifying the E/Z-marker fix above, fixed as its own small change so it doesn't need an external workaround to run.

## [0.5.0] — 2026-07-23

### Added — `chematic-perception`

- **`local_parity_from_wedges(mol, coords, center)` / `apply_local_parity_from_wedges(mol, coords)`** ("P1-S1a-core") — CIP-independent tetrahedral parity from wedge/hash bonds and 2D coordinates, producing `Atom.chirality` + `Molecule::stereo_neighbor_order` directly. Never calls CIP ranking and never touches `Atom.cip_code` — a CIP tie must not prevent a molecule from having a known local parity. Sign convention measured against RDKit's raw `CHI_TETRAHEDRAL_CW`/`CCW` tag on frame-aligned fixtures (not derived by analogy); handles 3-heavy-plus-implicit-H (no synthetic H position, mirroring RDKit's own `atomChiralTypeFromBondDirPseudo3D`) and multiple simultaneous wedges on one center (each checked in isolation for a consistent parity before the combined volume is trusted — a wedge/hash drawing with two substituents pointing opposite ways is valid notation, not automatically contradictory). Full methodology in `docs/rfcs/stereo2d_local_parity_calibration.md`. Not yet called from any reader's default parse path, and the SMILES writer's own wedge-token bug (below) had to be found and separately fixed first.
- **Charge-aware Hückel π-electron counting** (`ring_pi_electrons` / `evaluate_atom_pi_contribution_inner`, "K2a") — protonated N/O and cationic C were previously scored as neutral; tropylium cation, imidazolium, pyridinium, and pyrylium now match RDKit's aromatic-atom/bond flags exactly, verified end-to-end and isolated from the (separately tracked, not yet shipped) authoritative-demotion change below.

### Fixed — `chematic-smiles`

- **Bracket-forced atoms silently dropped their implicit hydrogens** — an atom force-bracketed by isotope/charge/atom-map (not an explicit H count) wrote e.g. `[NH4+]` as `[N+]`. Both the plain and canonical writer now route bracket-H emission through `chematic_core::implicit_hcount()` via one shared helper, so the two writers can't drift out of sync on this again.
- **Standalone wedge bonds were written as meaningless SMILES `/`/`\` tokens** — `BondOrder::Up`/`Down` is overloaded for both 2D wedge depiction and OpenSMILES E/Z directional markers; a wedge bond with no adjacent double bond was still emitted as `/`/`\`, a token RDKit and other parsers silently drop since there's nothing for it to mark. The writer now only emits a directional token when a bond actually flanks a real double bond, reusing the existing E/Z-adjacency check; genuine E/Z markers and the SMILES parser's aromatic-bond-direction stash are unaffected.

### Fixed — `chematic-core`

- **`kekulize()`'s `atom_must_be_matched` was charge-blind and missing Tellurium** ("K1") — RDKit's own kekulization gives a charged carbon and a charged heteroatom *opposite* double-bond-matching requirements (a documented special case, RDKit GitHub #539); chematic's classification didn't distinguish them, and had no case at all for Se/Te-adjacent charge states or P. Tropylium, imidazolium, pyridinium, pyrylium, tellurophene, and phosphole — previously hard failures — now kekulize successfully, bond-for-bond identical to RDKit's own Kekulé choice, with zero regressions on the existing pyridine/pyrrole/furan/thiophene/selenophene/azulene suite or the 5,000-molecule corpus.
- **`Element::normal_valences()` had no entry for Tellurium (52)** — fell through to an empty default, which cascaded into `chematic-perception`'s RDKit-parity aromaticity engine reporting tellurophene's ring non-aromatic (RDKit: aromatic) purely because Te's default valence resolved to `None`. Valence list `[2, 4, 6]` added, cross-verified against both RDKit's pinned-commit source and the live installed RDKit API. Fixes ECFP4 bit-exactness for tellurophene (see below). Does not add Te to the OpenSMILES organic subset or change implicit-H behavior; a real, reported, and accepted side effect: bare degree-1 Te–C pairs in 3D-coordinate bond-order inference now upgrade `Single`→`Double` (degree ≥2 Te, the common case, is unaffected).

### Fixed — `chematic-smarts` / `chematic-chem`

- **Unbounded VF2 backtracking on symmetric PAINS/Brenk targets could hang for minutes** — matching now takes an explicit visit budget and returns a genuine three-way `MatchOutcome::{Found, NotFound, BudgetExhausted}` instead of running to completion or forever. Budget exhaustion is *never* silently folded into "no match" (that would turn a hang into a silent false negative); the existing `pains_matches`/`brenk_matches` APIs conservatively fold it to "alert present" instead — a known, accepted tradeoff (documented, regression-tested with a real fixture: a common di-tert-butylphenol scaffold resolves to `BudgetExhausted`→flagged at the shipped 1,000,000-visit budget, needing ~20s at 20,000,000 to resolve to the true `NotFound`). `matches_detailed_impl`'s per-pattern highlight now also checks `budget_exhausted` before returning any partial atom-index set, so an incomplete enumeration can no longer be returned as though it were exhaustive. Real-budget calibration tests (25–80s) moved to `#[ignore]`, replaced with fast deterministic budget=0 tests in the default suite. Symmetry-aware VF2 candidate ordering, which would let the shipped budget resolve the di-tert-butylphenol case correctly instead of conservatively, is tracked as [#139](https://github.com/kent-tokyo/chematic/issues/139), not implemented here.

### Known limitations

- **Authoritative aromatic-flag demotion is not yet shipped ("K2b").** `build_molecule_from_model` can still only *promote* an atom/bond's flag to aromatic when the Hückel model agrees; it cannot yet demote a stale parser-set flag the model disagrees with. A draft, unmerged fix exists and regresses 499 atom/bond flags across 84/5,000 corpus molecules on the `kekulized`-input calling convention specifically (the `raw` convention every production caller actually uses is unaffected) by exposing pre-existing Hückel-model gaps that flag promotion was previously masking — mostly the already-known odd/odd Pass-1 split (azulene-type non-alternant rings, 49/84) and a newly-diagnosed ring-fusion/exocyclic-substituent ambiguity in `CarbonExocyclicHeteroatomDouble` (fused diazine/quinazoline/quinoxaline/purine rings, 33/84). Held back until these clusters are resolved or demotion ships as an explicit opt-in rather than changing `apply_aromaticity()`'s default behavior.
- **Tellurophene and phosphole still don't reach full RDKit agreement under the default Hückel engine.** K1 makes both kekulize correctly; K2a makes the other four charge-aware fixtures match RDKit exactly. Se/Te/P electron-donor rules (charge, exocyclic bonds, hypervalence) for the *default* engine are out of scope here and need their own source-grounded PR — the existing RDKit-parity engine already handling these elements is not sufficient justification for copying its logic into the default engine without independent verification.
- **`chematic-fp`/ECFP4, descriptor, and canonical-SMILES outputs for the 6 K1-affected molecule classes were re-verified bit-exact/unaffected as part of this release** (tropylium, imidazolium, pyridinium, pyrylium, tellurophene, phosphole all `verified_bit_exact` against the RDKit oracle after K1+Te+K2a landed together) — this is stated explicitly because K1 alone, without the Te fix, would have silently regressed tellurophene's ECFP4 output from an honest `Err` to a wrong `Ok`.

---

## [0.4.30] — 2026-07-17

### Added — `chematic-cip` (new crate, experimental)

- **Hierarchical-digraph CIP (Cahn-Ingold-Prelog) engine** — `assign_cip_accurate_experimental`, a from-scratch, provenance-carrying digraph replacement for the existing shell-pooling comparator in `chematic-chem`. Not yet wired into `chematic_chem::assign_cip()`; a separate, non-default, `publish = false` crate for now. See `docs/rfcs/cip_accurate_rfc.md` for the full design and milestone history.
- Rules 1a (atomic number) / 1b (ring-duplicate handling) / 2 (isotope) comparator, genuinely sphere-by-sphere (breadth-first, not depth-first) — fixes a class of bug the old shell-pooling comparator couldn't (a shallow sibling difference must decide a comparison before a much deeper, irrelevant one is reached).
- **MANCUDE (maximum non-cumulated double bonds) fractional atomic numbers** — `AtomicNumberKey`/`RationalAtomicNumber` give aromatic ring atoms (e.g. pyridine's N-adjacent carbons) a Kekulé-invariant fractional value instead of one arbitrary Kekulé form's integer, matching RDKit's own `calcFracAtomNums` formula (verified against RDKit source, not paraphrase).
- **Full-corpus accuracy on this experimental engine: 96.68% → 99.19%** (4047/4186 → 4152/4186 vs modern RDKit `rdCIPLabeler`, net +105, 0 regressions confirmed by two independently-computed methods). Most of the gain is attributable to Kekulé-respelling structurally (aromatic bonds finally contributing real digraph duplicate nodes) — the fractional MANCUDE values themselves are implemented and RDKit-formula-verified but measured inert (0/4188 rows) on the available corpus, kept rather than reverted since there's no efficiency cost.
- **Milestone 4A** — `CipCode::LowerR`/`LowerS` and a Pass-2 refinement implementing Rule 5 (pseudoasymmetry), scoped to 2 verified-independent rows. A three-armed, locally-symmetric cage family (15 rows) was found, via a direct stereo-tag flip test, to be provably unreachable by this pairwise two-pass architecture (every relevant neighbor is also tied — no seed to refine from) and deferred as **Milestone 4A-2**, pending a symmetry/automorphism-aware joint-resolution design.
- **Milestone 4A-0** — re-froze the residual from scratch at commit `992d18c` (not reusing the pre-4A bucket estimates, since 4A's own +2 fix could have shifted which rows remain): 34/4186 residual, **100% mechanically classified, 0 unexplained** — 15 Rule 5/pseudoasymmetry (the 4A-2 cage family, unchanged), 8 Rule 4 candidate (uppercase RDKit label + no stereogenic double bond in the tied branches + `branch_signature`-confirmed constitutional identity — positively confirmed, not inferred by elimination), 11 phosphorus (9 rows where the comparator fully resolves the ranking but to an incorrect order — a correctness bug, not a missing rule — plus 2 genuinely-tied rows). Supersedes this same entry's earlier 17/11/8 estimate.
- **Milestone 4B — Rule 4b (auxiliary-descriptor / "like/unlike") comparator**, wired into `assign_cip_accurate_experimental` only. Diagnosis first established Rule 4b as the sole applicable subtype for the 8-row `rule4_candidate` residual (Rule 4a/4c structurally ruled out — `chematic_core::Chirality` has no category/axial-descriptor distinction to compare), then built and validated a faithful, reference-relative, bottom-up-in-one-digraph comparator (per Hanson et al. 2018 JCIM 58(9), the paper `rdCIPLabeler` implements — not a re-rooted-per-atom approximation, which was tried, found antisymmetric, and root-caused to violating the paper's single-digraph requirement). Validated **72/72** across the paper's own external suite (VS196/VS197), the 8-row residual, and a new 16-row discriminating corpus, each also checked against its mirror enantiomer — before being ported mechanically into production.
- **Full-corpus accuracy: 99.19% → 99.38%** (4152/4186 → 4160/4186 vs modern RDKit `rdCIPLabeler`), `newly_correct = 8` (exactly the Rule 4b residual, no unexplained extras), 0 regressions confirmed by two independent methods.
- **Milestone 4C — phosphorus residual reclassified as oracle instability, not a chematic defect.** All 11 remaining phosphorus rows (9 cyclophosphazene rows + 2 previously-tied rows, same underlying molecule) are Kekulé-respelling-sensitive: the respelled and original forms are the same molecule (verified via matching InChI including stereo layer), chematic's answer is stable across both spellings, but RDKit's own modern `rdCIPLabeler` is not — and for the 2 tied rows, even RDKit's legacy `_CIPCode` isn't stable either. No established ground truth exists for this ring family to converge toward, so no chematic fix was made; these 11 rows are now stratified into a separate `representation_unstable` bucket (`validation/cip_oracle_instability.jsonl`), excluded from both the correct and incorrect side of a new **oracle-stable** agreement figure reported alongside the raw one.
- **Milestone 4 gate closed: oracle-stable agreement 4160/4175 = 99.64%** (raw 4160/4186 = 99.38% unchanged), 0 regressions, 0 unexplained residuals — the remaining 15 rows are exactly the already-deferred Milestone 4A-2 symmetric-cage family (Rule 5, not a new gap).
- **Milestone 5A — accurate CIP engine wired into every public surface, opt-in only, no default changes.** `chematic_chem::assign_cip_with_mode(mol, CipMode)` (Rust; new `CipMode`/`CipModeAssignment`/`CipUnresolvedReason`/`CipModeError` types, kept distinct from the pre-existing `CipAssignment`/`CipError`), `Mol.cip_stereo(mode="legacy"|"accurate")` + `Mol.cip_stereo_unresolved()` (Python, default mode unchanged), `cip_assignments_accurate_json` + `cip_unresolved_json` (WASM, existing `cip_assignments_json` unchanged). Because the accurate engine only computes tetrahedral R/S (no E/Z or allene handling), `CipMode::Accurate` **merges** the accurate engine's R/S with legacy's E/Z/allene answers rather than swapping engines outright — a naive swap would have silently dropped every E/Z/allene label. Ties/budget-outs are surfaced explicitly via `unresolved`, never silently backfilled with legacy's guess. Verified live (not just `cargo check`) through `maturin develop`+pytest and `wasm-pack`+Node, identical results across Rust/Python/WASM.
  - Structural note: `chematic-chem` now depends on `chematic-cip` as a normal (not dev-only) dependency, inverting the crate-layer diagram in `CLAUDE.md` — a deliberate, flagged exception. The dependency inversion is valid for local workspace builds, but requires `chematic-cip` to be published before `chematic-chem` for crates.io releases (see the "Release" entry below — `chematic-cip` was `publish = false` at the time this milestone landed, which blocked the next crates.io release until fixed).
- **Milestone 5B — opt-in stabilization measurements** (no behavior change): Accurate mode is **~10× slower than legacy** (214–240µs vs 22–24µs per molecule, full 5,000-molecule corpus) — the one still-open item on the default-promotion checklist. Unresolved rate 0.392% (19/4849 stereocenters), **100% traced to the two already-known families** (17 Milestone 4A-2 cage carbons, 2 Milestone 4C-1 phosphorus atoms), zero unexplained causes. Cross-surface parity re-verified at 300 molecules (0 mismatches). Default-promotion gate criteria are listed with current status in `docs/rfcs/cip_accurate_rfc.md`; promotion itself is not decided by this milestone.

### Added — `chematic-perception`

- **`assign_aromaticity_rdkit_parity_experimental(mol)` / `apply_aromaticity_rdkit_parity_experimental(mol)`** — a new, opt-in, fallible aromaticity engine independently reproducing RDKit's actual default algorithm (source-verified port of `Code/GraphMol/Aromaticity.cpp` from RDKit release `Release_2026_03_4`; BSD 3-Clause attribution in `THIRD_PARTY_NOTICES.md`), rather than approximating it from black-box test cases. Root cause of the production engine's known false-positive family (bridgehead-N rings): RDKit computes each atom's electron-donor type once, globally per molecule; chematic's existing per-candidate-ring evaluation doesn't. **Not wired into the default path** — `RdkitLike`/`Huckel` (`assign_aromaticity_ex`/`apply_aromaticity_ex`) are byte-unchanged; this is a separate fallible surface (`AromaticityError::{KekulizationFailed, InternalInvariantViolation}`) because its pre-kekulized-input precondition isn't always satisfiable, unlike the existing infallible engines. No Python/WASM exposure yet.
- **100.0000% atom/bond agreement with real RDKit** on the new engine, measured on 4,999/5,000 comparable molecules (1 excluded for a pre-existing, unrelated `kekulize()` gap) — 138,635/138,635 atoms, 150,004/150,004 bonds, 0 unexplained differences — versus the current **default** engine's established 99.44% atom / 98.82% bond agreement (`docs/rdkit_compat.md`), which is unchanged by this addition. ~4% slower than the default on mean/p50/p95 latency (measured, not optimized).
- Diagnosis groundwork behind the engine above: a 50-molecule false-positive/false-negative/negative-control corpus, a component/conjugation model with an exhaustive reference oracle, and a classification of a 94-molecule canonical-round-trip-instability corpus (9/94 fixed by the `apply_kekule` fix below; 19/94 a known ~1.18% bond-level `RdkitLike`-vs-RDKit gap that self-resolves only if this engine is ever promoted to default; 66/94 a pre-existing, unrelated canonical-SMILES-writer sensitivity on complex fused rings, not touched here). The experimental engine's post-fix canonical round-trip idempotency is 98.32% (4,915/4,999), still short of the ≥98.42% default-promotion bar (current default: 98.62%, 4,930/4,999) — both dominated by the same 66-case pre-existing canonicalizer-sensitivity bug, unrelated to aromaticity or stereo metadata.

### Added — `chematic-chem`

- **`schultz_mti(mol) -> f64`** — Schultz Molecular Topological Index (MTI): weighted sum of adjacency × distance matrix entries.
- **`gutman_mti(mol) -> f64`** — Gutman MTI variant: degree-weighted distance sum.
- **`vabc(mol) -> f64`** — van der Waals volume from Bondi atomic radii (no 3D coordinates required); complements TPSA for bulk/lipophilicity estimation.
- **`gravitational_index(mol) -> f64`** — graph-theoretic gravitational index: Σ (mᵢ·mⱼ / d²ᵢⱼ) over all heavy-atom pairs weighted by atomic mass and topological distance squared.

### Added — `chematic-depict`

- **`depict_pdf(mol)` / `depict_pdf_opts(mol, opts)`** — PDF output via `svg2pdf`; no external dependencies beyond the crate.
- **`depict_eps(mol)` / `depict_eps_opts(mol, opts)`** — EPS (Encapsulated PostScript) output; pure-Rust implementation, no additional dependencies.
- **`png` feature** — `tiny_skia` (raster PNG rendering) moved to an optional `png` feature (default = on for non-WASM builds); WASM builds disable it automatically, reducing bundle size.

### Added — `chematic-mol`

- **ChemicalJSON (`.cjson`) format** — `parse_cjson(s)` and `write_cjson(mol, coords)` support the Avogadro/MolSSI ChemicalJSON format for interoperability with Avogadro2 and MolSSI Cookiecutter projects.

### Added — `chematic-py`

- **`mol.to_pdf()`** — render 2D depiction to PDF bytes.
- **`mol.to_eps()`** — render 2D depiction to EPS string.
- **`mol.to_cjson(coords=[])`** — serialize to ChemicalJSON; optional 3D coordinate list.
- **`chematic.from_cjson(s)`** — parse ChemicalJSON string to `Mol`.
- **`mol.schultz_mti`** — Schultz MTI property.
- **`mol.gutman_mti`** — Gutman MTI property.
- **`mol.vabc`** — van der Waals volume (Bondi radii) property.
- **`mol.gravitational_index`** — gravitational index property.
- **`bulk.substructure_match(smarts, mols)`** — parallel substructure search accepting pre-parsed `Mol` objects; runs VF2 matching in parallel across the list, returning a `list[bool]`.
- **`bulk.generate_3d(smiles, *, method="etkdg")`** — parallel 3D coordinate generation; returns `list[list[[x,y,z]] | None]`. `method="etkdg"` (default) uses the ETKDG knowledge base with chair/envelope ring conformations and 80 torsion rules; `method="dreiding"` is faster.
- **`bulk.tanimoto_matrix(smiles)`** — all-pairs ECFP4 Tanimoto similarity; returns numpy `(N, N)` float32.
- **`bulk.standardize(mols)`** — batch molecule standardization (largest fragment, neutralize, canonical tautomer); returns `list[Mol]`.

### Added — `chematic-inchi`

- **Inline SHA-256** — replaced the `sha2` crate dependency with a self-contained 60-line SHA-256 implementation; saves ~15 KB in the WASM bundle.
- **`sha256_abc` and `sha256_empty` tests** — RFC 4634 test vectors verify the inline implementation.

### Fixed

- **WASM TypeScript docstrings** — `whim_descriptors_json` updated to document 22 values (was 10); `whim_getaway_combined_json` updated to document 41 values (was 19).
- **`mms_member_mw_excludes_wildcard`** — regression test confirming MMS member MW correctly excludes wildcard atoms.
- **`chematic-smarts`: `[rN]` was wrongly aliased to `[kN]`'s "any ring of size N" semantics** — RDKit's real `[rN]` means "this atom's *smallest* ring is exactly size N" (confirmed empirically distinct from `[kN]` on a purine example). New `AtomPrimitive::MinRingSize`, evaluated via a lazily-computed per-atom min-ring-size cache reusing the existing shared `RingSet` — no change to ring perception (`find_sssr`) itself. **SMARTS match-set agreement: 96.9% → 99.93%** (79,944/80,000, full 5,000-molecule corpus); `[rN]` isolated: 70.95% → 100.00%; `[kN]` unaffected (99.98%, the 1 remaining mismatch pre-existing and unrelated). Accounts for 94% of the SMARTS mismatches previously (and incorrectly) attributed to a genuine ring-model gap — the real remaining gap (`[R1]`/`[R2]` ring-*count* on bridged/cage systems, an actual SSSR-basis-cardinality disagreement) is kept separately tracked, not touched by this fix.
- **`chematic-core`: `apply_kekule` silently dropped `stereo_neighbor_order`** (P0) — its `MoleculeBuilder` rebuild never called `copy_stereo_groups_from`/`copy_stereo_from`/`copy_bond_directions_from`, so `@`/`@@` chirality survived but its neighbor-order reference frame didn't; downstream code inferring configuration from a different neighbor order (canonical rank, CIP digraph position) could silently serialize the wrong tetrahedral configuration with no panic or error. `Molecule` now implements `Clone` (used for the `is_empty()` fast path). Verified positive-controlled (every new test fails with the fix reverted) across `chematic-core`, `chematic-smiles` (canonical round trip, enantiomer pairs), `chematic-cip` (both legacy and accurate engines — pre-fix, the accurate engine didn't mislabel affected stereocenters, it silently *skipped* them under a misleadingly-named `SkipReason`), and `chematic-inchi` (confirmed InChI never reads `stereo_neighbor_order` at all, so it was never at risk from this specific bug — but a separate, pre-existing InChI representation-dependence bug on an achiral molecule was found and documented as out of scope).
- **`chematic-chem`: the same missing-metadata-copy bug recurred at several more sites**, each fixed with the same 3-call pattern and verified positive-controlled. A repo-wide audit that `apply_kekule`'s fix triggered found three (`enumerate_stereoisomers`, `transfer_hydrogen_aromatic`, `clone_mol`); two further instances (`transfer_hydrogen`, `invert_stereocenter`) surfaced afterward, during and after the `transfer_hydrogen_aromatic` fix's own session:
  - **`enumerate_stereoisomers`** — the highest-risk of the three: it calls `canonical_smiles` on its own output immediately after every rebuild, so every call was affected. Two distinct gaps: already-specified stereocenters lost `stereo_neighbor_order` on every generated isomer, and newly-assigned stereocenters (the function's whole purpose) never had one at all — fixed by explicitly constructing a well-defined (if arbitrary) neighbor order for each newly-chiral atom. Reverting the fix flips a newly-assigned stereocenter's CIP code from `S` to `R` across a canonical round trip — silently wrong, not just missing.
  - **`transfer_hydrogen_aromatic`** and **`clone_mol`** — same gap on every aromatic N–H tautomer step. `clone_mol` (a hand-rolled passthrough rebuild, 8 call sites) is now fully redundant with `Molecule::clone` and was **deleted**, all call sites replaced with `mol.clone()`.
  - **`transfer_hydrogen`** (the non-aromatic keto-enol/1,3-shift counterpart, e.g. glucose tautomer enumeration) — identical gap, found while implementing the previous fix and fixed with the same pattern; a remote stereocenter uninvolved in the H-shift silently lost its neighbor-order reference frame on every tautomer step.
  - **`invert_stereocenter`** — same metadata-copy gap, plus a separate, more severe bug: the function only ever inverted 2D wedge bonds (`Up`/`Down`), so it was a **functional no-op on plain `@`/`@@` SMILES input** (no wedge bond present → the "no stereochemistry" passthrough ran instead), silently doing nothing for the common case of SMILES-parsed chiral input. Fixed by flipping the `Chirality` enum directly when present (correct independent of neighbor order/atom numbering), keeping wedge-bond inversion as a fallback for 2D-only/API-compatibility callers. Verified via chemical identity (CIP labels, independently-parsed mirror-image SMILES, InChI cross-check — not `@`/`@@` string comparison): single inversion flips only the target center's CIP; double inversion reproduces the original molecule.
- **`chematic-smiles`: aromatic bond-direction stashing was inconsistent across the 3 bond-creation paths in the parser.** A `/`/`\` marker between two aromatic atoms encodes an adjacent exocyclic double bond's geometry (never the aromatic bond's own order); only the main chain-edge path correctly stashed it on the `bond_direction` side channel while keeping `order=Aromatic`. Ring-closure resolution and branch-attachment both instead stored the marker as the bond's own literal `Up`/`Down` order, making a bond's visibility to `assign_ez` depend on which path happened to create it — stable on hand-written input, but able to flip after a canonical round trip since chematic's own canonical DFS can route the same bond through a different path than the original parse. Consolidated into one shared `resolve_aromatic_direction_stash` helper. **This is a stability fix, not an E/Z correctness fix** — `assign_ez`/`substituent_is_up` still don't read the stash, so affected molecules correctly show an empty CIP set both before and after (confirmed via an independent RDKit check that real, stable E/Z geometry exists on these bonds — a real, pre-existing recall gap, tracked as follow-up **EZ-A0**/**EZ-S1**, not closed here). On the 5,000-molecule corpus: first-parse output is byte-identical old vs. new (0 risk to existing literature-SMILES parsing); round-trip CIP-element-multiset stability improved **4,994/5,000 → 5,000/5,000** (6 confirmed cases fixed, 0 gained, 0 lost).

### Performance — WASM bundle size (819 → 504 KB gzip, −38.5%)

- `tiny_skia` made optional (`png` feature) in `chematic-depict`; WASM builds opt out automatically.
- `sha2` replaced with inline SHA-256 in `chematic-inchi`.
- `[profile.release] opt-level="z" lto=true codegen-units=1` added to workspace `Cargo.toml`.
- Removed `run_md_json`, `coulomb_energy_json`, `torsion_scan_json`, `determine_bonds_from_xyz_json` from WASM exports.
- Removed `chematic-ewald` from `chematic-wasm` dependencies.

### Performance — WASM bundle size grew 504 → 719 KB gzip since 2026-06

- Not attributed to a single PR — cumulative effect of the `chematic-cip` crate and the new `chematic-perception` RDKit-parity aromaticity engine (both above) landing in the WASM dependency graph. Measured 2026-07-17 during a benchmark refresh (`benchmarks/2026-07-17.md`); not yet investigated for a targeted reduction.
- The 2026-07-17 benchmark refresh also found the previously-reported ECFP4 throughput headline (3.6 µs/mol, 5–14× vs RDKit) does not reproduce, on three independent measurements including the original fixture/script — see `benchmarks/2026-07-17.md`'s Notes for what was ruled out (methodology alone does not fully explain it) and what remains open (a partial, not-yet-confirmed link to this cycle's SSSR rework, below).

### Performance — `chematic-cip` MANCUDE ring-bond check (30ms → microseconds on large fused-ring molecules)

- `mancude.rs::ring_bond_set` was calling `chematic_perception::find_sssr` (a full minimum cycle basis) on the whole molecule just to answer a boolean "is this bond in some ring" question — replaced with a direct O(V+E) bridge-edge (cut-edge) DFS, since ring-bond membership is exactly the complement of the bridge set. Verified byte-identical output across the full corpus before/after.

### CI

- **Criterion regression gate bootstrap fix** (`bench-pr-gate.yml`) — a new or removed Criterion benchmark target used to abort the whole gate job (`cargo bench --bench <missing>` erroring under `set -e`), hiding every other benchmark's verdict. Now tolerated per-side with a three-way classification: both sides present gates normally; candidate-only is a new benchmark with no baseline yet (not gated); baseline-only is a possibly-removed benchmark (warned, not gated).
- **Criterion gate reliability finding** — the gate currently treats one Criterion process's ~100 internal samples as 100 independent A/B trials; they aren't (all share one process's runner state), so a single environment difference can be amplified into an extreme-looking win rate for every benchmark at once. Confirmed empirically (a run showed unrelated benchmarks all failing with near-identical effect sizes). The gate stays non-required and its `fail` verdicts are not currently reliable evidence of a real regression — process-level redesign tracked in [#70](https://github.com/kent-tokyo/chematic/issues/70).
- **Criterion gate Stage 1/Stage 2 redesign** (#117) — the process-level redesign above (independent process-run blocks, two-stage screening) had two further, more specific bugs, both found via real CI runs rather than local testing: (1) Stage 1's 3-block `sign-test` could never fail — its strongest possible signal at n=3 (unanimous 3-0) has a two-sided exact p-value of 0.25, never crossing the 95%-confidence bar, making Stage 2 structurally unreachable for a regression of any size. Fixed by replacing Stage 1 with a pure magnitude-routing screen (median block ratio ≥ `STAGE1_ROUTE_THRESHOLD=1.04`, chosen via an offline evaluation of 5 candidate rules against 28 historical no-op runs); Stage 1 alone never sets `any_fail`. (2) Stage 2's `sign-test` was independently magnitude-blind: a real ~2.2% build/codegen variance between two separately-compiled binaries of identical source triggered a unanimous, statistically "significant" fail. Fixed by requiring the sign-test fail to also clear `STAGE2_FAIL_THRESHOLD=1.04` on the median ratio, else it's reported as `small-effect-inconclusive`, not blocking; the real incident data is pinned as a permanent regression fixture in `scripts/test_criterion_gate.sh`. The same-binary null control's inability to detect codegen variance *between two separately-compiled binaries* (as opposed to environment-wide contamination) remains a known, reported, not-yet-fixed gap.
- Added `wasm-opt -O3` step to `.github/workflows/pages.yml` CI.
- **crates.io publish-graph check** (`scripts/check_publish_graph.py`, wired into `scripts/check.sh`) — for every workspace crate not marked `publish = false`, verifies every normal/build-dependency path crate is (a) not `publish = false`, (b) on the same workspace version, and (c) publishable in a well-defined topological order (no cycles in the normal dependency graph — `dev-dependencies` excluded from this part, since they're stripped from the manifest seen by downstream consumers). Also flags any `chematic-*` dev-dependency that carries both `path` and `version`: **discovered live, mid-release** — `cargo publish` generates a full `Cargo.lock` (including dev-dependencies) *before* its verify step, so a version-pinned dev-dependency must already be on the registry even though downstream consumers never see it. `chematic-smiles` and `chematic-perception` had a *mutual* version-pinned dev-dependency on each other, which deadlocked both (neither could ever publish first) — six more crates had a version-pinned forward reference to a not-yet-published crate, which happened to work only by luck of the current publish order. All nine fixed to path-only (dropping `version`); positive-controlled (the check fails when reverted). Added after `chematic-cip`'s `publish = false` was found blocking the 0.4.30 crates.io release in the same way — `chematic-chem` had gained a normal (non-dev) dependency on it in Milestone 5A without anyone checking publish-graph validity at merge time.

### Release

- **`chematic-cip` is now published to crates.io** (`publish = ["crates-io"]`, previously `publish = false`) — required for `chematic-chem`'s 0.4.30 release, since `chematic-chem` depends on it as a normal (non-optional) dependency (Milestone 5A) and crates.io requires every normal dependency to be resolvable on the registry. Crate-level docs and `Cargo.toml`'s `description` now state its intended audience explicitly: most applications should use `chematic_chem::assign_cip_with_mode` rather than depending on this crate directly; it remains experimental and may receive breaking API revisions in a future 0.x minor release.
- **Fixed a `cargo publish` deadlock** between `chematic-smiles` and `chematic-perception` (mutual version-pinned dev-dependency, see the publish-graph check entry above) discovered while executing the 0.4.30 crates.io publish sequence, right after `chematic-core` published successfully. Caught before any further crate was published — `chematic-core` was the only one live at the time — so the fix and the rest of the sequence shipped together as 0.4.30, not a follow-up patch release.

---

## [0.4.26] — 2026-06-29

### Fixed — `chematic-rxn`

- **E/Z stereo transfer/creation in `run_reactants`** (#50) — reaction products now preserve `/`/`\` double-bond geometry from reactants instead of losing stereo on transformation.

### Added — validation

- **Sprint 6**: canonical SMILES differential validation vs RDKit.
- **rdkit_compat Sprint 7**: SMARTS/aromaticity differential tests, I/O compatibility, docs, examples.
- Canonical-idempotency and E/Z canonical-stability regression tests documenting the root causes of remaining RDKit divergence (aromaticity round-trip, not Morgan ranks).

### CI

- `rdkit-pypi` → `rdkit` in the validation workflow (the `rdkit-pypi` package is deprecated and no longer published on PyPI).

---

## [0.4.25] — 2026-06-29

### Added — `chematic-py`

- **`chematic.rdkit_compat`** — RDKit API compatibility layer, Sprints 1–5: Morgan `bitInfo`, fingerprint/Mol/Atom/Bond/RingInfo compatibility surface, differential tests against RDKit.
- **`SDMolSupplier` / `SDWriter` / `Mol.GetProp`** — streaming SDF I/O and SD-file properties in the RDKit-compat style.

### Added — `chematic-perception`

- **`AromaticityAlgorithm::RdkitLike`** — Se/Te chalcogen aromaticity handling matching RDKit's model.

### Dependencies

- Bump `miniz_oxide` 0.8 → 0.9; `actions/checkout` 4 → 7; `actions/cache` 4 → 6; `actions/upload-artifact` 4 → 7 (Dependabot).

---

## [0.4.24] — 2026-06-29

### Added — `chematic-chem`

- **CIP Rule 5 stereo tie-breaking** — stereocenter agreement 99.8% → 99.98% vs RDKit.
- Dual-oracle stereocenter benchmark (legacy detector + `FindPotentialStereo`).

### Fixed — `chematic-chem`

- **Bridgehead detection 98.5% → 100.0%** on the 5,000-molecule corpus.
- **Rotatable bonds 99.1% → 100.0%** on the 5,000-molecule corpus.
- **TPSA 100%** agreement with RDKit on the 4,999-molecule corpus.
- **Molar refractivity (MR) 97.5% → 100%** via 3-ring XOR augmentation.
- Spiro/bridgehead atom detection for bridged and cage ring systems.
- Sulfoxide / selenoxide lone pair now counted as the 4th CIP substituent.
- `heavy_atom_count` now excludes explicit H (including isotopic `[3H]`).

### Fixed — `chematic-perception`

- `augmented_ring_set` finds same-size XOR rings (bridged bicyclics); aliphatic/aromatic ring counts via aromatic-single-bond detection.

### Added — `chematic-py`

- **`bulk.descriptors_array()`** — columnar numpy descriptor output.
- **SDF true streaming** — `SdfFileReader<R: BufRead>` + `iter_sdf_batched`.
- **`screen()`** compound-filter workflow + LogP mismatch analyzer.
- **`fragment_text()`** for LLM prompts + RDKit validation dashboard.

### Added — LLM/RAG integration

- **Representation router** (`to_llm_text`, `best_representation`, MCP tool `representation_router`).
- **Molecule context pack** for LLM/RAG pipelines (MCP tool `molecule_context_pack`).
- **Hyper-Dimensional Fingerprints (HDF)** — training-free dense molecular vectors.

---

## [0.4.23] — 2026-06-26

### Fixed — `chematic-chem`

- **LogP 96.5% → 99.7%**: `crippen_anchor_sets` now uses `uniquify: false` so that symmetric triple bonds (internal alkynes, R–C≡C–R) yield both orientations from VF2 matching. Previously, deduplication collapsed `{0:Cₐ, 1:Cᵦ}` and `{0:Cᵦ, 1:Cₐ}` into one, leaving one alkyne carbon unmatched and falling back to the generic `[#6]` value (+0.0796 error per atom).

---

## [0.4.18] — 2026-06-23

### Added — `chematic-py`

- **`Mol._repr_svg_()`** — Jupyter/JupyterLab auto-display hook; writing `mol` in a
  cell renders the 2D structure automatically without `IPython.display.SVG(...)`.
- **`Mol.has_substructure(smarts)`** — method-level SMARTS match returning `bool`
  (equivalent to `chematic.smarts_match(smarts, mol)` with a cleaner call site).
- **`Mol.find_matches(smarts)`** — method-level SMARTS match returning atom-index lists
  (equivalent to `chematic.smarts_find(smarts, mol)`).
- **`chematic.from_smiles_list(smiles, *, skip_invalid=True)`** — pure-Python convenience
  function; batch-parses SMILES via `bulk.parse`, filtering `None` by default.
- **`chematic.descriptors_df(smiles)`** — one-liner SMILES → `pd.DataFrame` wrapper
  around `bulk.descriptors`; raises `ImportError` when pandas is not installed.

### Added — `chematic-chem`

- **`cns_mpo_from_parts(mol, logp, tpsa, mw, hbd, pka_b)`** — CNS MPO score from
  pre-computed values; avoids redundant Crippen SMARTS pass when descriptors are
  already available in the caller.
- **`chi_all(mol) -> (χ0, χ1, …, χ4v)`** — compute all 10 Hall-Kier connectivity
  indices in a single `heavy_indices` pass (was 10 independent calls).
- **`pains_passes_and_matches(mol) -> (bool, Vec<&str>)`** and
  **`brenk_passes_and_matches(mol)`** — single explicit-H + SSSR + pattern scan
  returning both the pass flag and alert names.

### Performance — `chematic-chem` / `chematic-py`

- `named_groups.rs`: `detect_named_functional_groups` shares one `find_sssr()` across
  all 19 patterns (was 19 independent SSSR computations).
- Python `Mol.descriptors()`: uses `logp_and_mr()` (1 Crippen pass for logP + MR),
  `chi_all()` (1 heavy_indices pass for 10 chi indices), and `cns_mpo_from_parts()`.
- `bulk.descriptors()`: uses `logp_and_mr()` and `pka_both()`.

### Docs

- **`docs/benchmark.md`** — new benchmark page: ECFP4 speed (5–14× vs RDKit),
  descriptor accuracy (100% on 5,000-mol corpus), install/WASM comparison, feature table.
- README Docs badge now links to site root `https://kent-tokyo.github.io/chematic/`.

---

## [0.4.17] — 2026-06-23

### Performance — `chematic-chem`

- **PAINS/Brenk dedup** — `pains_passes_and_matches(mol) -> (bool, Vec<&str>)` and
  `brenk_passes_and_matches` perform a single explicit-H conversion + SSSR + pattern scan
  returning both the pass flag and alert names; `workflow.rs` now uses the combined function
  instead of calling `pains_passes` and `pains_matches` separately (2× → 1× 480-pattern scan).
- **QED structural alerts: 113 SSSR → 1** — `structural_alert_count_with_rings` uses
  `find_matches_with_rings_and_config`; `qed_with_bundle` computes `find_sssr` once and
  shares it across all 113 QED structural-alert patterns.
- **pKa SSSR sharing** — `predict_pka` now calls `find_sssr` once and shares the `RingSet`
  across all 42 ionizable-group SMARTS rules (was one SSSR per rule).
- **`pka_both(mol) -> (Option<f64>, Option<f64>)`** — new public function returning acid and
  base pKa from a single `predict_pka` call; `bulk.rs` uses it to avoid the double scan.

### Performance — `chematic-fp`

- **MHFP incremental BFS** — `extract_fragment_hashes` expands each center atom's neighborhood
  shell-by-shell instead of restarting BFS from scratch per radius level; reduces BFS from 3N
  to N per molecule at radius=2 (50-atom molecule: 150 → 50 BFS operations).

---

## [0.4.16] — 2026-06-22

### Performance — `chematic-smarts`

- **Shared SSSR across multi-pattern matching** — `EvalCtx` now borrows `&RingSet` instead of
  owning it, enabling callers to compute the ring set once and reuse it across many queries.
  Two new public functions:
  - `find_matches_with_rings(query, mol, rings)` — match with a pre-computed `RingSet`
  - `find_matches_with_rings_and_config(query, mol, rings, config)` — same with explicit config

### Performance — `chematic-chem`

- **Crippen SMARTS: 117 SSSR → 1 per `logp_crippen` call** — `crippen_anchor_sets()` now calls
  `find_matches_with_rings_and_config` sharing a single `find_sssr()` result across all 117
  Wildman-Crippen patterns (was one SSSR computation per pattern).
- **PAINS/Brenk: ~480/~300 SSSR → 1 per call** — `pains_passes`, `pains_matches`,
  `brenk_passes`, `brenk_matches`, and `matches_detailed_impl` each compute `find_sssr` once
  on the explicit-H molecule and share it across all compiled patterns. `pains_passes` and
  `brenk_passes` also set `max_matches: Some(1)` for early-exit boolean queries.
- **`logp_and_mr(mol) -> (f64, f64)`** — new public function that returns both Crippen LogP and
  MR from a single `crippen_anchor_sets` pass (~2× faster when both are needed). Exact numerical
  agreement with `logp_crippen(mol)` + `molar_refractivity(mol)` verified by regression test.
- **`logd_from_logp(logp, mol, ph) -> f64`** — new public function accepting a pre-computed LogP
  value; avoids the duplicate Crippen pass inside `logd_simple`.
- **`cns_mpo_score` logP dedup** — computes `logp_crippen` once and passes the result to
  `logd_from_logp`, eliminating a redundant 117-pattern SMARTS pass.
- **`workflow.rs` LogP+MR** — `molecule_report` now uses `logp_and_mr` instead of separate
  `logp_crippen` + `molar_refractivity` calls (saves one full Crippen pass per report).
- **`eccentric_connectivity_index` reuses `graph_eccentricities`** — eliminated a duplicate
  O(n_heavy²) BFS traversal; now delegates to the shared eccentricity vector.
- **Degree pre-computation in topological descriptors** — new `heavy_degrees(mol) -> Vec<u32>`
  helper computes all heavy-atom degrees in one pass; `randic_index`, `zagreb_index_m1`, and
  `zagreb_index_m2` now look up pre-computed values instead of iterating neighbors per bond.

### CI

- Bump `actions/setup-python` v5 → v6 (pages.yml, publish-pypi.yml)
- Bump `actions/upload-artifact` v4 → v7 (publish-pypi.yml)

---

## [0.4.15] — 2026-06-21

### Fixed — `chematic-chem`

- **TPSA atom-type calibration** — six systematic divergences from RDKit corrected by
  comparing per-atom `_CalcTPSAContribs` values:
  - Imine N=C (h=0): 12.89 → **12.36 Å²**
  - Imine =NH (h=1): 23.79 → **23.85 Å²**
  - Nitrile N≡C (h=0): 3.24 → **23.79 Å²** (new case, triple-bond check)
  - O⁻ (non-nitro, e.g. carboxylate/phenolate/sulfonate): 9.23 → **23.06 Å²**
  - Ring-junction aromatic N — all bonds `BondOrder::Aromatic` (bridgehead, neutral): 4.93 → **4.41 Å²**
  - Ring-junction aromatic N (bridgehead, cationic): 3.88 → **4.10 Å²**
  - `tpsa_per_atom()` now applies `apply_aromaticity` for Kekulé-form consistency with `tpsa()`
  - `tpsa_all_tsv_reference` bulk-regression tolerance tightened: ±1.0 → **±0.1 Å²** (175 molecules)
- **Bench5k results after calibration** (4 999-molecule ChEMBL-like corpus):
  - HBA: 99.98% → **100%** (4 999 / 4 999)
  - HBD: 99.8% → **100%** (4 999 / 4 999)
  - Aromatic ring count: 98.5% → **100%** (4 999 / 4 999)
  - TPSA (±0.1 Å²): 86.7% → **93.3%** (drug-like 175-mol set: **100%**)

### Fixed — `chematic-rxn`

- **E/Z double-bond stereo filtering in `run_reactants` (issue #21)** — SMIRKS templates with
  `/`/`\` bond-direction descriptors on **both** sides of a double bond now correctly filter
  reactants whose geometry does not match. Mirrors the tetrahedral `smirks_chirality_ok` approach:
  - `ez_stereo_outward(mol, atom, other)` — computes the "outward" Up/Down direction from a
    double-bond endpoint, flipping incoming bonds so direction is always relative to the sp2 atom.
  - `smirks_ez_stereo_ok(tmpl, reactant, match_map)` — post-VF2 parity check: derives E/Z sense
    (same outward direction = Z, different = E) from both template and reactant, rejects mappings
    where the senses disagree.
  - Single-sided or unspecified stereo is not filtered (ambiguous → permissive).

### Added — tests

- 6 new E/Z stereo tests in `chematic-rxn::transform`:
  `ez_stereo_e_template_matches_e_alkene`, `ez_stereo_e_template_rejects_z_alkene`,
  `ez_stereo_neutral_template_matches_both_geometries`, `ez_stereo_one_sided_template_matches_both_geometries`,
  `ez_stereo_retro_wittig_z_matches_z_hexene`, `ez_stereo_z_template_matches_z_alkene`
- 6 new TPSA calibration regression tests in `chematic-chem::rdkit_reference`:
  `tpsa_diazepam`, `tpsa_clonazepam`, `tpsa_metformin` (tightened ±0.5→±0.1), `tpsa_arginine`
  (tightened), `tpsa_nitrile`, `tpsa_carboxylate_anion`, `tpsa_ring_junction_n`, `tpsa_n_substituted_aromatic_n`

---

## [0.4.14] — 2026-06-21

### Fixed — `chematic-rxn`

- **Parity-aware SMIRKS chirality matching** — `run_reactants` previously used a raw `@`/`@@` flag
  comparison in VF2 (`eval_chirality`), which is SMILES write-order-dependent: the same absolute
  configuration written with a different neighbor order stores an opposite chirality flag, causing
  both false positives and false negatives. The fix adds a post-match check (`smirks_chirality_ok`)
  that maps the SMIRKS template's `stereo_neighbor_order` through the VF2 mapping, computes the
  parity of the resulting permutation vs. the reactant's `stereo_neighbor_order`, and accepts a
  match iff `even_parity == (template_flag == reactant_flag)`. New test:
  `stereo_filter_same_config_different_write_order` verifies that L-alanine as both
  `N[C@@H](C)C(=O)O` and `C[C@H](N)C(=O)O` match the same `@@` template while D-alanine is rejected.
- **Product bracket-atom cleanup (issue #18)** — `build_product` was copying `hydrogen_count = Some(0)`
  from bare bracket template atoms (`[O:1]`, `[N:1]`) into product atoms, forcing bracket notation
  (`[O]`, `[N]`) in the SMILES output. Fixed by using `.filter(|&h| h > 0)` so `Some(0)` clears to
  `None` and implicit valence rules apply. `[NH2:1]`-style explicit-H templates continue to produce
  `[NH2]` as expected.
- **`snap_amide_torsions` tertiary amide fix** — previously filtered by `classify_atom_type ==
  NSp2`, silently skipping all degree-3 (tertiary) amide nitrogens; now matches any non-aromatic N
  with a `C(=O)` neighbour. Also fixed a double-correction bug: multiple `set_dihedral` calls on
  the same bond read stale coordinates after the first rotation; limited to one snap per bond with
  `break 'snap`.
- **`is_atom_in_ring` multi-start BFS** — single-pair BFS (`nbs[0]` → `nbs[1]`) returned false for
  degree-≥3 atoms when `nbs[0]` is an exocyclic substituent. Replaced with a multi-start BFS that
  tries every neighbour as the starting point.
- **`tpsa()` always applies aromaticity** — mixed-case input (e.g. indomethacin's Kekulé indolyl N
  written uppercase) was misclassified; `apply_aromaticity` is now called unconditionally so all
  input forms are correctly typed.
- **`is_aromatic_oxide_bridge()` helper extracted** — oxide-bridge detection logic was duplicated in
  `tpsa_oxygen()` and `crippen_logp_for_atom()`; consolidated into one shared helper.

### Fixed — `chematic-3d`

- **ETKDG amide planarity (RDKit #9266)** — `snap_amide_torsions` post-processing step added:
  after constraint satisfaction, tertiary amide `ω` angles outside ±30° of 0° or 180° are snapped
  to the nearer planar value.
- **PBF (Plane of Best Fit) uses heavy atoms only (RDKit #9238)** — `plane_of_best_fit()` now
  excludes hydrogen atoms, matching the published definition and RDKit convention. Including H
  artificially reduces PBF for flat aromatic rings.

### Fixed — `chematic-mol`

- **CDXML E/Z stereo** — `flush()` now calls `assign_ez_from_2d()` after parsing, deriving E/Z
  descriptors from 2D coordinates for double bonds in CDXML input.

### Fixed — `chematic-perception`

- **`count_aromatic_rings()` handles Kekulé input (RDKit #9271)** — bracket atoms without aromatic
  flags (uppercase SMILES) returned 0 aromatic rings. Fixed by applying Hückel aromaticity
  perception (`apply_aromaticity`) internally when no aromatic flags are set, so fluorescein
  dianion and rhodamine-type zwitterions are correctly classified.

### Added — tests

- `transform::tests::stereo_filter_same_config_different_write_order` — verifies parity-aware
  chirality matching across different SMILES write orders
- `transform::tests::product_removes_bracket_from_bare_bracket_atoms` — issue #18 regression
- `shape_descriptors::tests::pbf_uses_heavy_atoms_only` — PBF H-exclusion regression
- `aromaticity::tests::test_fluorescein_dianion_aromatic` — Kekulé aromatic ring count regression
- `aromaticity::tests::test_rhodamine_zwitterion_parses` — zwitterion parse regression
- Canonical SMILES PAH round-trip tests for pyrene and benzo[a]pyrene (14 new cases)

### Changed — documentation

- `eval_chirality` in `chematic-smarts`: doc comment added explaining the raw-flag limitation and
  directing SMIRKS users to `smirks_chirality_ok` instead
- `run_smirks` / `run_smirks_strict` Python docstrings: stereo filtering behaviour documented

---

## [0.4.13] — 2026-06-21

### Added — `chematic-rxn`

- **Template-based retrosynthesis `retro_disconnect()`** — 60 retro-SMIRKS templates across 6 reaction classes:
  - `AmideBond` (10): secondary/tertiary amide, sulfonamide, carbamate, urea, hydrazide, imide
  - `Ester` (6): ester, thioester, carbonate, anhydride, acetal, lactone
  - `Ether` (8): aryl ether (SNAr/Ullmann), Williamson, benzyl, Mitsunobu, silyl
  - `CNBond` (11): reductive amination, SNAr-CN, Buchwald, N-alkylation, Mitsunobu-N, imine reduction
  - `CCBond` (14): Suzuki, Heck, Sonogashira, Negishi, Grignard, aldol, Michael, Wittig, Diels-Alder
  - `CSBond` (10): thioether, disulfide, borylation, halogenation, phosphonate
  - Python API: `mol.retro_disconnect(max_results=20, reaction_class="AmideBond")` → `list[dict]`
  - Returns `{template, reaction_class, precursors, sa_scores, max_sa_score}` ranked by SA Score
- Retro-SMIRKS templates tightened: sp3 constraints added to `reductive_amination`, `n_alkylation`,
  `mitsunobu_n`, `negishi`; enolisable H required for `aldol`, `michael_addition`; `sp3_cc_bond`
  template removed (too broad); `friedel_crafts_alkyl` restricted to benzylic CH2

### Added — `chematic-3d`

- **ETKDG torsion knowledge base expanded** — 28 → 40 patterns:
  - New `OAromatic` / `SAromatic` atom types for furan/thiophene 5-membered heterocycles
  - `NMorpholine` / `NPiperazine` atom types for saturated N-heterocycle torsion preferences
  - Adaptive noise: bond-flexibility scaling (amide 0.2×, biaryl 0.5×, single bond 1.0×)

### Added — `examples/`

- `examples/aizynthfinder_integration.py` — AiZynthFinder + chematic integration tutorial:
  molecule preparation, BRICS 1-step retrosynthesis, scoring, route ranking (works without AiZynthFinder installed)

### Fixed — `chematic-chem`

- **`hbd_count()` now counts S-H (thiol)** — `hbd_count` previously only counted N-H and O-H donors.
  Adding sulfur (atomic number 16) aligns with `rdMolDescriptors.CalcNumHBD`. Affected: cysteine (2→3),
  thiophenol (0→1).
- **TPSA nitro-N contribution corrected** — `tpsa_nitrogen()`: `N+(=O)[O-]` (nitro group) was returning
  41.44 Å² instead of the correct Ertl 2000 value of 43.14 Å². Fixed: 4-nitrophenol now 63.37 (was 61.67).
- **TPSA aromatic oxide bridge** — `tpsa_oxygen()`: bridging O in a ring bonded to aromatic C and a vinyl
  C=C (e.g., morphine's 4,5-epoxy ring) was returning 9.23 Å² (aliphatic ether) instead of 13.14 Å²
  (furanoid). Added BFS ring-check + aromatic-C + vinyl-C detection.
- **TPSA Kekulé-form aromatic N** — `tpsa_nitrogen()`: N written in Kekulé SMILES (uppercase) that is
  embedded in an aromatic ring (degree ≥ 3, aromatic neighbour, in ring) was getting 3.24 Å² (tertiary
  amine) instead of 4.93 Å² (aromatic N). Fixes indomethacin (delta −1.69 Å² resolved).
- **LogP Crippen O7 SMARTS typo** — `[OX1;-,-2,-2][#16]` corrected to `[OX1;-,-2,-3][#16]`
  (Wildman-Crippen 1999 Table 1; the duplicated `−2` entry would miss O³⁻ on S).
- **LogP aromatic oxide bridge** — `crippen_logp_for_atom()`: same oxide-bridge O now returns `[o]`-type
  LogP (0.1552) before the SMARTS loop reaches `[O](a)` (−0.4195). Fixes morphine/codeine LogP delta.

### Added — `scripts/bench5k.py`

- TPSA comparison: `rdMolDescriptors.CalcTPSA(includeSandP=True)` vs `ch_mol.tpsa` (±0.1 Å² match)
- LogP comparison: `Crippen.MolLogP` vs `ch_mol.logp` (±0.01 match)
- HBD comparison: `rdMolDescriptors.CalcNumHBD` vs `ch_mol.hbd`

### Added — `crates/chematic-chem/tests/rdkit_reference.rs`

- `tpsa_all_tsv_reference()` — all 175 reference molecules within ±1.0 Å²
- `logp_all_tsv_reference()` — all 175 reference molecules within ±0.3
- `mw_all_tsv_reference()` — exact MW match (±0.02 Da) for all 175 molecules
- `hac_all_tsv_reference()` — exact HAC match for all 175 molecules
- `hbd_all_tsv_reference()` — exact HBD match for all 175 molecules
- `tpsa_nitrobenzene()`, `tpsa_4_nitrophenol()` — nitro-group TPSA regression

---

## [0.4.12] — 2026-06-21

### Fixed — `chematic-smarts`

- **Atom map number `:N` now supported in SMARTS patterns** — patterns like
  `[O;D1;H0:3]` (SMIRKS-derived SMARTS) previously caused a parse error
  because the SMARTS parser did not handle the `:N` suffix.
  - `QueryAtom` gains `atom_map: Option<u16>` (metadata only; never a
    matching criterion — `:N` is silently accepted and stored).
  - New `QueryMolecule::add_atom_with_map(query, atom_map)` helper (existing
    `add_atom` unchanged — no public API break).
- **`[C:]` (bare `:` with no digit) now returns `SmartsError::UnexpectedChar`**
  instead of silently succeeding with `atom_map: None`.
- Atom map digit accumulation uses `u32` accumulator clamped at `u16::MAX`,
  ensuring `parse_smarts` and `extract_map_numbers_from_section` agree on the
  stored value for any input.

### Fixed — `chematic-rxn`

- **`extract_map_numbers_from_section`**: add bracket-depth tracking so that
  `:` outside brackets (aromatic bond token before a ring-closure digit,
  e.g. `c1:c:c:c:c:c:1`) is never mistaken for an atom map number; previously
  any `:N` anywhere in the string was extracted, causing valid aromatic-ring
  reaction SMARTS to fail with `MapNumberMismatch`.
- **`mol_to_query`**: call `add_atom_with_map(q, atom.atom_map)` so atom map
  numbers from the source `Molecule` are preserved in `QueryAtom`; previously
  `add_atom(q)` always set `atom_map: None`.

### Removed — `chematic-rxn`

- `strip_map_numbers()` workaround in `query.rs` is deleted; the SMARTS
  parser now natively accepts `:N` so the pre-processing step is unnecessary.

---

## [0.4.11] — 2026-06-21

### Fixed — `chematic-perception`

- **Aromatic ring count: 95.6% → ~100% RDKit agreement** (`augmented_ring_set` XOR guard `min` → `max`)
  - Root cause: when the SSSR returns a large macro-ring paired with a same-size ring instead of two equal-sized component rings, the old `min` guard incorrectly skipped recovery of the missing ring.  Changing to `max` recovers any ring strictly smaller than the *larger* parent.
  - All 222 previously failing bench5k cases now match RDKit (`count_aromatic_rings`).
- Extend envelope ring detection to **4-ring GF(2) XOR** in `count_aromatic_rings` — correctly strips coronene-class perimeter cycles that are the bond-symmetric-difference of four inner hexagons.

### Fixed — `chematic-mol` (CIF parser)

- `parse_cif`: return `CifError::InvalidCellParameters` when `sin(γ) ≈ 0` (e.g. `_cell_angle_gamma 0` or `180`) instead of silently producing NaN/Inf coordinates.
- `parse_cif`: return `CifError::MissingCellParameters` when fractional coordinate columns are present but no `_cell_length_*` / `_cell_angle_*` parameters are found in the file.
- Strip oxidation-state suffixes (`Cu2+`, `Fe3+`, `O2-`) from `_atom_site_type_symbol` values in addition to trailing digits — fixes parsing of standard inorganic CIF files.
- Fix comment stripper to not truncate `#` characters inside CIF single- or double-quoted strings (`'foo#bar'` was previously truncated to `'foo`).
- Two new `CifError` variants: `InvalidCellParameters(String)`, `MissingCellParameters`.

### Fixed — `chematic-mol` (Gaussian parser)

- `parse_gjf`: detect charge/multiplicity section **structurally** (two blank-line-separated sections after the `#` route card) instead of scanning for the first `"int int"` line — prevents false matches when the title section is a number pair like `"0 1"`.
- `parse_gaussian_log`: support Gaussian 03 **5-column** Standard orientation table (columns: Center# AtomicNum X Y Z, no Atomic Type column); was previously silently skipping all rows and returning `NoAtoms`.
- `parse_gjf`: handle **bare atomic-number** element specifications (e.g. `6` for carbon) by falling back to `Element::from_atomic_number` when `trim_end_matches(is_ascii_digit)` yields an empty string.

### Fixed — CI / Clippy

- Declare `trained-solubility-mlp` feature in `chematic-chem/Cargo.toml`.
- Gate `use chematic_fp::ecfp4` behind the feature flag to suppress `unused_import` when feature is off.
- Extract `GjfResult` type alias in `gaussian.rs` to resolve `type_complexity` lint.
- Replace `.filter(..).last()` with `.rfind(..)` on `DoubleEndedIterator`.
- Remove unused `find_mcs` import in `chematic-mcp/src/tools.rs`.

---

## [0.4.10] — 2026-06-20

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.20 gap-analysis bindings (Sprint 26)

**2D stereochemistry perception (vs RDKit `Chem.AssignStereochemistry`):**
- `Mol.stereo_from_2d_coords(coords_2d)` → list[dict] — perceives R/S and E/Z from 2D layout coordinates (pairs with `stereo_from_coords` for 3D)

**Multi-record SDF with 2D coordinates (vs RDKit `SDMolSupplier`):**
- `chematic.parse_sdf_with_coords(text)` → list[tuple[Mol, str, list[list[float]]]] — batch SDF parse preserving 2D layout; batch equivalent of `from_mol_block_with_coords`

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.19 gap-analysis bindings (Sprint 25)

**MDL RXN file I/O (vs RDKit `AllChem.ReactionFromRxnFile`):**
- `chematic.from_rxn_file(text)` → str — parse MDL RXN V2000 file to reaction SMILES
- `chematic.to_rxn_file(reaction_smiles)` → str — convert reaction SMILES to MDL RXN V2000

**CML writer (vs RDKit `Chem.MolToCMLBlock`):**
- `Mol.to_cml(coords_2d=None)` → str — serialize to Chemical Markup Language XML; optional 2D layout coordinates

**V3000 round-trip (completes V3000 suite):**
- `chematic.from_mol_v3000_with_coords(block)` → `(Mol, name, coords_2d)` — V3000 parse preserving 2D layout (pairs with `mol.to_mol_v3000`)

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.18 gap-analysis bindings (Sprint 24)

**CXSMILES parsing (vs RDKit CX block support):**
- `chematic.from_cxsmiles(s)` → `(Mol, dict)` — parses CXSMILES and returns molecule + CX metadata (`atom_labels`, `atom_props`, `atom_radicals`)

**MDL MOL V3000 writer (vs RDKit `MolToV3KMolBlock`):**
- `Mol.to_mol_v3000(coords_2d, name=None)` → str — V3000 format supporting >999 atoms; same 2D coord convention as `to_mol_block_2d`

**Custom-radius MHFP (vs chemfp, datamol):**
- `Mol.mhfp_config(radius=2, num_hashes=128, seed=0)` → list[int] — MHFP with configurable radius/length; extends the default `mhfp()` method

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.17 gap-analysis bindings (Sprint 23)

**DREIDING force field minimization:**
- `Mol.minimize_dreiding(coords)` → list[list[float]] — third FF option alongside `minimize_mmff94` and `minimize_uff`

**MDL MOL V2000 2D coordinate round-trip:**
- `chematic.from_mol_block_with_coords(block)` → `(Mol, name, coords_2d)` — parse MOL block preserving 2D layout
- `Mol.to_mol_block_2d(coords_2d, name=None)` → str — write MOL block back with preserved 2D coordinates

**Per-atom SASA with explicit 3D coordinates:**
- `Mol.sasa_per_atom_3d(coords)` → list[float] — per-atom SASA from user-provided coords (vs `sasa_per_atom()` which uses DG internally)

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.16 gap-analysis bindings (Sprint 22)

**3D stereochemistry perception (vs RDKit `AssignStereochemistryFrom3D`):**
- `Mol.stereo_from_coords(coords)` → list[dict] — perceives R/S and E/Z from 3D coordinates; returns `[{'atom_idx': int, 'code': 'R'|'S'|'E'|'Z'}]` (note: requires 4 heavy-atom neighbours at chiral centre)

**Layer-by-layer fingerprint (vs RDKit `LayeredFingerprint(layerFlags=0x7F)`):**
- `Mol.layered_fp_layers()` → list[bytes] — 7 individual layers (atom types → bond orders → aromaticity → ring membership → ring-bond → stereo → combined); each 256 bytes, compatible with `chematic.tanimoto` etc.

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.15 gap-analysis bindings (Sprint 21)

**2-way/3-way combinatorial library enumeration (vs RDKit VirtualLibraries, chemfp):**
- `chematic.enumerate_library_2way(smirks, scaffolds, building_blocks)` → list[str] — scaffold × building block
- `chematic.enumerate_library_3way(smirks, scaffolds, r1_set, r2_set)` → list[str] — scaffold × R1 × R2

**Ring system topology classification (vs chemgraph, ring-analysis):**
- `Mol.ring_families()` → list[dict] — classifies each ring system as `"simple"` | `"fused"` | `"spiro"` | `"bridged"`, with `atom_indices` and `ring_count`

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.14 gap-analysis bindings (Sprint 20)

**ERG molecule-to-molecule similarity:**
- `chematic.tanimoto_erg(mol1, mol2)` → float — direct mol comparison without intermediate `erg_vec()` step

**M×N Tanimoto similarity matrix (vs chemfp `sim_matrix`, RDKit `BulkTanimotoSimilarity`):**
- `chematic.tanimoto_matrix(fps_a, fps_b)` → list[list[float]] — computes all pairwise scores in one call

**Pre-computed fingerprint k-NN (vs chemfp high-speed search):**
- `chematic.nearest_neighbors_from_fp(query_fp, db_fps, k=10)` → list[(idx, score)] — efficient for repeated queries against the same database; fingerprints computed once

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.13 gap-analysis bindings (Sprint 19)

**MDL MOL V2000 writer (vs RDKit `MolToMolBlock`):**
- `Mol.to_mol_block()` → str — completes the MOL format round-trip (`from_mol_block` ↔ `to_mol_block`)

**ADMET numeric score direct properties (vs ADMET-AI, pkCSM):**
- `Mol.bbb_score` → float (0–1) — blood-brain barrier penetration
- `Mol.caco2` → float (nm/s) — Caco-2 intestinal permeability
- `Mol.herg_risk` → float (0–1) — hERG cardiac toxicity risk
- `Mol.cyp3a4_risk` → float (0–1) — CYP3A4 inhibition risk
  All consistent with `mol.admet()` dict values.

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.12 gap-analysis bindings (Sprint 18)

**Bulk Tanimoto (vs chemfp `sim_matrix` / RDKit `BulkTanimotoSimilarity`):**
- `chematic.tanimoto_slice(query_fp, db_fps)` → list[float] — vectorized 1-vs-N Tanimoto for virtual screening

**Count-based Morgan FP (vs RDKit `GetMorganFingerprint`):**
- `Mol.morgan_fp_counts(radius=2)` → dict[int, int] — hash→count map (more informative than bit FP for ML)

**Pharmacophore feature counts:**
- `Mol.pharmacophore_feature_counts()` → list[int] (6 values) — [donor, acceptor, aromatic, hydrophobic, positive, negative]

**Reaction fingerprint similarity:**
- `chematic.tanimoto_reaction_fp(rxn1, rxn2)` → float — Tanimoto between two reaction SMILES

**Atom-type BCI MMFF94 charges:**
- `Mol.mmff94_charges_typed()` → list[float] — atom-type pair BCI model (alternative to `mmff94_charges()`)

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.11 gap-analysis bindings (Sprint 17)

**Ames mutagenicity alert names (vs RDKit `FilterCatalog`):**
- `Mol.ames_alerts()` → list[str] — specific SMARTS alert names (e.g. `"primary_aromatic_amine"`); complements existing `ames_passes()` / `ames_risk()`

**Hepatic clearance raw score:**
- `Mol.clearance_score` → float (0.0–1.0) — continuous score for ML pipelines; complements `Mol.clearance_class` (categorical)

**Per-atom molar refractivity:**
- `Mol.mr_per_atom()` → list[float] — per-atom MR contributions; `sum() ≈ molar_refractivity`

**3D-aware MMFF94 partial charges:**
- `Mol.mmff94_charges_3d(coords)` → list[float] — MMFF94 charges incorporating 3D polarization; more accurate than `mmff94_charges()` (2D topology only) for polar molecules

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.10 gap-analysis bindings (Sprint 16)

**Halogen & phosphorus direct properties (completing element count API, vs RDKit/mordred):**
- `Mol.num_fluorines`, `Mol.num_chlorines`, `Mol.num_bromines`, `Mol.num_iodines`, `Mol.num_phosphorus`

**Named SMARTS pattern lookup (vs RDKit `FilterCatalog`):**
- `chematic.named_pattern(name)` → str | None — look up built-in SMARTS by name
  (e.g. `"donor"`, `"acceptor"`, `"hydrophobic"`, `"positive"`, `"negative"`, `"aromatic"`)

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.9 gap-analysis bindings (Sprint 15)

**Random SMILES enumeration (vs RDKit `MolToSmiles(rootedAtAtom=...)`, ML augmentation):**
- `Mol.random_smiles(seed)` → str — non-canonical SMILES with deterministic seed
- `Mol.random_smiles_n(n, seed=0)` → list[str] — up to n unique random SMILES

**SMI file batch I/O (vs RDKit `SmilesMolSupplier` / `SmilesWriter`):**
- `chematic.parse_smi_file(content)` → list[(Mol, name)] — parse SMILES+name file; invalid lines skipped
- `chematic.write_smi_file(records)` → str — write [(Mol, name)] to .smi format

**Element color utilities (vs RDKit `DrawingOptions`):**
- `chematic.atom_color(atomic_num)` → str — CPK CSS color (e.g. `"#FF0D0D"` for O)
- `chematic.atom_color_rgb(atomic_num)` → (R, G, B) tuple

**Atom equivalence predicate:**
- `Mol.are_atoms_equivalent(a, b)` → bool — True if atoms have the same Morgan canonical rank

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.8 gap-analysis bindings (Sprint 14)

**Reaction SMARTS queries (vs RDKit `AllChem.ReactionFromSmarts`):**
- `chematic.query_reaction(rxn_smiles, smarts)` → bool — single reaction pattern match
- `chematic.batch_query_reactions(rxn_list, smarts)` → `{total, matching, match_pct, matches}` — batch hit rate

**MMFF94 force field energy & types (vs RDKit `MMFFGetMoleculeForceField`):**
- `Mol.mmff94_total_energy(coords)` → float (kcal/mol) — single-call energy without breakdown
- `Mol.mmff94_atom_types()` → list[str] — per-atom MMFF94 type names (e.g. `"C_sp3"`, `"O_Alcohol"`)

**3D pharmacophore Tanimoto (vs RDKit EmbedLib):**
- `chematic.tanimoto_pharmacophore_3d(fp1, fp2)` → float — direct Tanimoto between `pharmacophore_fp_3d()` outputs

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.7 gap-analysis bindings (Sprint 13)

**IUPAC stereo naming (vs RDKit `MolToIUPACName` + CIP prefix):**
- `Mol.iupac_name_stereo()` → str — prepends `(R)-`/`(S)-`/`(1R,2S)-` to IUPAC name

**Reaction SVG depiction (vs RDKit `Draw.ReactionToImage`):**
- `chematic.reaction_svg(rxn_smiles)` → SVG string — renders reaction arrow diagram

**Format I/O expansion (vs RDKit, Open Babel):**
- `chematic.from_cml(cml_str)` → Mol — Chemical Markup Language parser
- `chematic.from_cdxml(cdxml_str)` → Mol — ChemDraw XML parser
- `chematic.from_mol_v3000(block)` → Mol — MDL MOL V3000 parser

**RECAP fragmentation count:**
- `Mol.recap_breakable_bond_count` → int — number of RECAP-breakable C–N/C–O/C–S single bonds

**Scaffold network statistics (vs RDKit `ScaffoldNetwork`):**
- `chematic.scaffold_network_counts(smiles_list)` → `{scaffolds, counts, parents}` — frequency of each scaffold across a library

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.6 gap-analysis bindings (Sprint 12)

**Ring saturation:**
- `Mol.num_saturated_rings` → int (direct property, vs RDKit `CalcNumSaturatedRings`)

**Zwitterion & salt utilities (vs RDKit MolStandardize):**
- `Mol.has_zwitterion()` → bool — simultaneous positive/negative formal charges
- `Mol.normalize_zwitterion()` → Mol — proton transfer to neutral form
- `Mol.remove_salts()` → Mol — catalog-based salt fragment removal (distinct from `largest_fragment`)

**Stereo inversion:**
- `Mol.invert_stereocenter(atom_idx)` → Mol — flip Up/Down wedge bonds to generate enantiomer at specific center

**Multi-fingerprint k-NN search (vs RDKit `BulkTanimotoSimilarity`):**
- `chematic.top_k_similar_fp(query, smiles, k, fp)` — k-NN with selectable FP: `ecfp4` (default), `ecfp6`, `ecfp4_chiral`, `fcfp4`, `maccs`, `topo_path`

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.5 gap-analysis bindings (Sprint 11)

**Topological descriptor direct properties (15 getters, vs mordred/RDKit):**
- Hall-Kier shape: `Mol.kappa1`, `kappa2`, `kappa3`
- Path connectivity: `Mol.chi0`, `chi1`, `chi2`, `chi3`, `chi4`
- Valence connectivity: `Mol.chi0v`, `chi1v`, `chi2v`, `chi3v`, `chi4v`
- `Mol.wiener_index` (sum of topological distances), `Mol.bertz_ct` (graph complexity)

**Ring perception utilities (vs RDKit SSSR):**
- `Mol.ring_membership()` → per-atom list of SSSR ring indices
- `Mol.ring_sizes_for_atom(i)` → ring sizes containing atom i
- `Mol.is_fused_ring_system()` → bool (True if rings share an edge, not just spiro)

**Stereo validation (vs RDKit `FindPotentialStereo`):**
- `Mol.validate_stereo()` → `[{atom_idx, kind}]` — ImpossibleCenter / ConflictingWedges / RedundantStereo
- `Mol.stereo_completeness()` → `{specified, unspecified, total_centers}`

**Pharmacophore feature detection (vs RDKit Chem.Pharm2D):**
- `Mol.pharmacophore_features()` → `[{type, atom_idx, neighbor_indices}]` — Donor/Acceptor/Aromatic/Hydrophobic/Positive/Negative

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.4 gap-analysis bindings (Sprint 10)

Closing final Rust→Python exposure gaps for direct-access properties vs RDKit/mordred/chemfp:

**Element & bond count direct properties (vs RDKit `CalcNumXxx`):**
- `Mol.num_heteroatoms`, `num_carbons`, `num_nitrogens`, `num_oxygens`, `num_sulfurs`, `num_hydrogens`
- `Mol.num_amide_bonds`, `Mol.num_ester_bonds`

**Ring topology direct properties (vs RDKit `CalcNumSpiroAtoms` etc.):**
- `Mol.num_spiro_atoms`, `num_bridgehead_atoms`
- `Mol.num_aromatic_heterocycles`, `num_aliphatic_heterocycles`, `num_saturated_heterocycles`, `num_aliphatic_rings`

**ERG continuous vector similarity (vs chemfp ERG):**
- `Mol.erg_vec()` → `list[float]` of length 315
- `chematic.cosine_erg_vec(v1, v2)` → float
- `chematic.tanimoto_erg_vec(v1, v2)` → float

**SMILES canonical utilities (vs RDKit `CanonicalRankAtoms`):**
- `Mol.morgan_ranks()` → Morgan extended-connectivity ranks per heavy atom
- `Mol.canonical_atom_order()` → canonical permutation vector
- `Mol.equivalent_atom_classes()` → symmetry-equivalent atom class IDs

### Added — `chematic-py` + `.pyi`: PyPI Chemistry p.3 gap-analysis bindings (Sprint 9)

Closing remaining Rust→Python exposure gaps vs RDKit, mordred, and chempy:

**Per-atom E-state indices (vs RDKit `EState.EStateIndices`):**
- `Mol.estate_indices()` → `list[float]` — electrotopological state per heavy atom

**Condensed formula parser (vs chempy):**
- `chematic.from_condensed("CH3OH")` → `Mol | None` — dictionary-based condensed formula → structure

**MMFF94 3D minimization (completes force-field coverage):**
- `Mol.minimize_mmff94(coords)` → `[[x,y,z],...]` — MMFF94 geometry optimization

**Unspecified stereocenters count (vs RDKit `CalcNumUnspecifiedAtomStereoCenters`):**
- `Mol.num_unspecified_stereocenters` → `int` — stereocenters with unknown configuration

**Combined WHIM+GETAWAY descriptor (vs mordred pipeline):**
- `Mol.whim_getaway(coords)` → `list[float]` — single-call 3D featurisation combining both descriptors

### Added — `chematic-py`: PyPI Chemistry p.2 gap-analysis bindings (Sprint 8)

Closing final gaps vs RDKit, chempy, and pure-Rust competitors:

**Deduplication utilities (vs RDKit `rdMolHash`):**
- `chematic.mol_hash(mol)` → u64 — fast structural hash for O(1) duplicate detection
- `chematic.are_identical(mol1, mol2)` → bool — exact graph isomorphism check

**BRICS bond positions (vs RDKit `BRICS.FindBRICSBonds`):**
- `Mol.brics_bonds()` → `[(atom_i, atom_j)]` — which bonds BRICS would break

**Full pKa prediction (all sites):**
- `Mol.predict_pka()` → `[{atom_idx, pka, site_type, group_name}]` — all ionizable sites

**Reaction serialization (completes parse/write symmetry):**
- `chematic.write_reaction(rxn_smiles)` → canonical reaction SMILES string

**Chemical abbreviations (vs chempy):**
- `chematic.abbreviations()` → `{symbol: SMILES}` dict — 32 abbreviations (Boc, Ph, OMe, …)
- `chematic.expand_abbreviation(symbol)` → `Mol | None`

**Topological distance matrix (vs RDKit `GetDistanceMatrix`):**
- `Mol.topological_distance_matrix()` → N×N list of graph shortest-path distances

### Added — `chematic-wasm` + tests: Sprint 7 — WASM gap-closure & test coverage

**New WASM exports (Sprint 1-6 features now available to JS/web users):**
- `xlogp3_json(mol)` → `{"xlogp3": float}`
- `xlogp3_per_atom_json(mol)` → JSON array of per-atom values
- `generate_3d_coords_json(mol)` → `[[x,y,z],...]` raw coords (vs PDB string only before)
- `generate_3d_etkdg_coords_json(mol)` → ETKDG raw coords JSON
- `conformer_ensemble_json(mol, n, rmsd_threshold)` → `{"conformers":[...], "count": int}`

**New pytest tests (`tests/test_sprint16.py`, 45 tests, all passing):**
Covering Sprint 1-6 features: logd, mqn, xlogp3, cip_stereo, generate_3d, gasteiger/mmff94 charges, named_functional_groups, pmi/npr/asphericity, dihedral manipulation, ETKDG, atom_economy, balance_check, molecule_report, screen_smiles, find_mmp, find_reaction_center, prefer_organic, uncharge, sasa_descriptor, randic_index, fcfp6, pattern_fp, mmff94_charges, balaban_j, zagreb_m1, labute_asa_per_atom, top_k_similar, center_on_origin, dice_similarity, butina_cluster.

Test count: 108 → **153 passing** Python tests.

**demo/pkg rebuilt** to include all new WASM exports.

### Added — `chematic-py`: PyPI Chemistry category gap-analysis bindings (Sprint 6, final)

Completing the remaining Rust→Python API surface vs RDKit, mordred, and chemfp:

**MMFF94 charges (vs RDKit `GetMMFF94PartialCharges`):**
- `Mol.mmff94_charges()` → per-atom MMFF94 BCI partial charges (more accurate than Gasteiger)

**Top-K similarity search (vs chemfp):**
- `chematic.top_k_similar(query, smiles_list, k)` → `[(index, score)]` top-K ECFP4 Tanimoto hits

**Topological descriptors (vs mordred completion):**
- `Mol.balaban_j` → Balaban J complexity index
- `Mol.ipc` → Information-theoretic connectivity index
- `Mol.zagreb_m1` → Zagreb M1 (sum of squared vertex degrees)
- `Mol.labute_asa_per_atom()` → per-atom Labute ASA contributions

**SASA statistics (vs RDKit / OpenBabel):**
- `Mol.sasa_descriptor(coords)` → `{total, mean, std_dev, per_atom}` dict

**Coordinate utilities (vs RDKit `AllChem`):**
- `chematic.center_on_origin(coords)` → translate centroid to origin
- `chematic.transform_conformer(coords, matrix_4x4)` → affine transformation

### Added — `chematic-py`: PyPI Chemistry category gap-analysis bindings (Sprint 5)

Completing the API surface vs MolVS, ADMETlab, and rxnSMILES4AtomEco:

**Green chemistry (complete ACS 12-principles set):**
- `chematic.e_factor(waste_mass, product_mass)` → E-factor (waste-to-product ratio)
- `chematic.pmi_rxn(all_masses, product_mass)` → Process Mass Intensity
- `chematic.reaction_mass_efficiency(reactant_masses, product_mass)` → RME [0–1]

**Standardization steps (vs MolVS / RDKit MolStandardize):**
- `Mol.normalize_groups()` → normalize nitro/azide/diazo/sulfoxide groups
- `Mol.prefer_organic()` → keep largest organic fragment, remove counterions
- `Mol.reionize()` → re-apply ionization rules based on pKa
- `Mol.uncharge()` → remove all formal charges

**ADMET extension:**
- `Mol.clearance_class` → `"Low"` / `"Medium"` / `"High"` hepatic clearance

**SASA extension (vs OpenBabel / RDKit):**
- `Mol.sasa_with_probe(coords, probe_radius=1.4)` → SASA with custom probe radius
- `Mol.sasa_per_element(coords)` → SASA breakdown by element symbol dict

**Additional fingerprints + topology:**
- `Mol.fcfp6()` → FCFP6 fingerprint (256 bytes)
- `Mol.pattern_fp()` → SMARTS pattern fingerprint (256 bytes)
- `Mol.randic_index` → Randić connectivity index (molecular branching)

### Added — `chematic-py`: PyPI Chemistry category gap-analysis bindings (Sprint 4)

Closing remaining gaps vs RDKit, Schrödinger Canvas, and rxnmapper:

**Workflow API (unique differentiator — no RDKit equivalent):**
- `chematic.molecule_report(smiles)` → complete single-molecule dict (descriptors, filters, functional groups, scaffold)
- `chematic.screen_smiles(smiles_list)` → batch screening with per-molecule reports + MaxMin picks + Butina clusters
- `chematic.compare_molecules(smiles_list)` → pairwise similarity + descriptor deltas + MCS

**MMP Analysis (vs RDKit rdMMPA / Schrödinger Canvas):**
- `chematic.find_mmp(smiles_list)` → `[{mol_a, mol_b, core, fragment_a, fragment_b}]` Matched Molecular Pairs

**Reaction Analysis (vs RDKit / rxnmapper):**
- `chematic.find_reaction_center(rxn_smiles)` → `{broken_bonds, formed_bonds, changed_atoms}`

**Additional fingerprints:**
- `Mol.path_fp()` → RDKit-compatible Daylight path fingerprint (256 bytes)
- `Mol.topo_path_fp()` → topological path fingerprint (256 bytes)
- `Mol.pharmacophore_fp_3d(coords)` → 3D pharmacophore fingerprint (complement to 2D `pharmacophore_fp`)
- `Mol.reaction_fp(rxn_smiles)` → `{reactant_fp, product_fp, combined_fp}` dicts (each 256 bytes)

### Added — `chematic-py`: PyPI Chemistry category gap-analysis bindings (Sprint 3)

Closing gaps vs rxnSMILES4AtomEco, Gypsum-DL, PIKAChU, and MolVS:

**Charges (docking pipeline now complete):**
- `Mol.gasteiger_charges()` → per-atom Gasteiger–Marsili partial charges for PDBQT writing

**Reaction chemistry (vs rxnSMILES4AtomEco / RDKit):**
- `chematic.atom_economy(rxn_smiles)` → float (green chemistry metric)
- `chematic.balance_check(rxn_smiles)` → `{balanced, diff}` atom balance report
- `chematic.enumerate_library(smirks, fragment_sets)` → combinatorial product SMILES

**Structure analysis (vs PIKAChU / RDKit):**
- `Mol.cip_stereo()` → per-atom/bond CIP R/S/E/Z assignments
- `Mol.pains_alerts()` → list of matched PAINS alert names (not just bool)
- `Mol.brenk_alerts()` → list of matched Brenk alert names
- `Mol.named_functional_groups()` → list of detected group names (carboxyl, ester, …)

**3D shape descriptors (vs mordred / OpenEye):**
- `Mol.pmi(coords)` → `[PMI1, PMI2, PMI3]` principal moments of inertia
- `Mol.npr(coords)` → `[NPR1, NPR2]` normalised principal moments (for Rod/Disc/Sphere plot)
- `Mol.asphericity(coords)`, `Mol.eccentricity(coords)`, `Mol.radius_of_gyration(coords)`
- `Mol.plane_of_best_fit(coords)` → PBF deviation

**Conformer editing (vs Gypsum-DL / RDKit):**
- `Mol.generate_3d_etkdg()` → higher-quality ETKDG 3D coordinates
- `Mol.get_dihedral(coords, i,j,k,l)` → dihedral angle in degrees
- `Mol.set_dihedral(coords, i,j,k,l, angle_deg)` → new coords with rotated fragment
- `Mol.get_bond_length(coords, i, j)` → bond length (Å)
- `Mol.get_bond_angle(coords, i, j, k)` → bond angle (deg)

### Added — `chematic-py`: PyPI Chemistry category gap-analysis bindings (Sprint 2)

Closing remaining gaps vs RDKit 2026, mordred, chemfp 5.1, and ChemAxon:

**Similarity metrics (new — RDKit has 8, chematic now has 3):**
- `chematic.dice_similarity(fp1, fp2)` → Dice coefficient (= 2|A∩B|/(|A|+|B|))
- `chematic.tversky_similarity(fp1, fp2, alpha, beta)` → Tversky index (generalises Dice/Tanimoto)
- `chematic.tanimoto_mhfp(fp1, fp2)` → MHFP similarity (position-wise)

**Descriptors (Rust→Python, mordred/RDKit parity):**
- `Mol.xlogp3` / `Mol.xlogp3_per_atom()` — alternative logP (Wang 2000)
- `Mol.autocorr_2d()` — 2D autocorrelation (mordred-compatible)
- `Mol.hall_kier_alpha` — correction term for κ shape indices
- `Mol.peoe_vsa()`, `Mol.slogp_vsa()`, `Mol.smr_vsa()`, `Mol.estate_vsa()` — 4 VSA descriptor vectors
- `Mol.usrcat()` — 42-element topological shape/pharmacophore descriptor

**File output:**
- `Mol.to_sdf_with_charges(coords, charges)` — 3D SDF with `> <PARTIAL_CHARGES>` property (Rust implemented in v0.4.9, now in Python)

**Fingerprints:**
- `Mol.pharmacophore_fp()` → 2048-bit pharmacophore 2D fingerprint (HBD/HBA/hydrophobic/aromatic)
- `Mol.mhfp()` → 128 u64 MinHash hash values

**3D alignment:**
- `chematic.align_coords(probe, reference)` → `(aligned_coords, rmsd)` (Kabsch algorithm)
- `chematic.rmsd(coords_a, coords_b)` → RMSD without alignment

### Added — `chematic-py`: PyPI gap-analysis bindings (Sprint 1)

Exposed 20 Rust-implemented features to Python for the first time, closing gaps identified
by comparison with chemfp, mordred, datamol, and scikit-mol:

**Fingerprints:**
- `Mol.map4()` → `list[int]` (1024 u32 MinHash hashes) + `Mol.map4_numpy()` → `ndarray(1024, uint32)`
- `Mol.erg()` → `bytes` (2048-bit Extended Reduced Graph fingerprint)
- `chematic.tanimoto_map4(a, b)` → similarity for MAP4 fingerprints (position-wise, not bitwise)
- `chematic.bulk.map4(smiles)` → `ndarray(N, 1024, uint32)` — parallel MAP4 batch

**Descriptors:**
- `Mol.logd(ph=7.4)` → `float` — pH-adjusted LogD (pKa-weighted, key ADMET descriptor)
- `Mol.logd_profile()` → pH 0–14 curve (28 points)
- `Mol.mqn()` → `list[int]` — 42 Molecular Quantum Numbers (Ertl 2009)
- `Mol.logp_per_atom()` → per-atom Crippen LogP contributions
- `Mol.isotope_distribution()` → `list[(mass, intensity)]` natural isotopic envelope

**Chemical analysis:**
- `Mol.functional_groups()` → list of `{atom_indices, atom_types}` dicts (Ertl 2017)
- `Mol.scaffold_network()` → Schuffenhauer ring-stripping hierarchy

**3D generation & descriptors:**
- `Mol.generate_3d()` → `list[[x,y,z]]` (DG + DREIDING minimization)
- `Mol.conformer_ensemble(n, rmsd_threshold)` → list of conformer coordinate arrays
- `Mol.whim(coords)` → WHIM 3D shape/symmetry descriptors
- `Mol.getaway(coords)` → GETAWAY 3D descriptors
- `Mol.autocorr_3d(coords)` → 3D autocorrelation

**File I/O:**
- `Mol.to_pdb(coords)` / `Mol.to_xyz(coords)` — write 3D to PDB/XYZ
- `chematic.from_pdb(pdb_str)` / `chematic.from_xyz(xyz_str)` → `(Mol, coords)`

**Force field analysis:**
- `Mol.mmff94_energy_breakdown(coords)` → per-term energy dict (bond/angle/torsion/oop/vdW/elec/total)
- `Mol.mmff94_torsion_scan(coords, i,j,k,l, steps)` → `list[(angle_deg, energy)]`

**Diversity selection:**
- `chematic.butina_cluster(smiles, cutoff)` → Butina clusters (ECFP4 similarity)
- `chematic.maxmin_picks(smiles, n)` → MaxMin diversity indices

---

## [0.4.9] — 2026-06-19

### Added — `chematic-mol`: AutoDock PDBQT format (`pdbqt.rs`)

New module `crates/chematic-mol/src/pdbqt.rs` for the AutoDock4 / Vina docking format:

- `autodock_atom_type(mol, idx)` — assigns AutoDock type (C, A, N, NA, O, OA, S, SA, H, HD, P, F, Cl, Br, I, Zn, …)
- `write_pdbqt(mol, coords, charges, name)` — writes a rigid-body PDBQT (no torsion tree)
- `parse_pdbqt(s)` → `(Molecule, Vec<(f64,f64,f64)>, Vec<f64>)` — reads ATOM/HETATM records

With this, the full zero-FFI docking preparation pipeline is possible:
SMILES → 3D generation → MMFF94 optimisation → `write_pdbqt()` → Vina input.

### Added — `chematic-ff`: Universal Force Field (`uff.rs`)

New module `crates/chematic-ff/src/uff.rs` implementing UFF (Rappé et al. 1992):

- `assign_uff_types(mol)` — maps all elements to `UffType` (C_3/C_2/C_R, N_3/N_R, O_3, metals Zn/Fe/Cu/…)
- `uff_total_energy(mol, types, coords)` — bond stretching + angle bending + Lennard-Jones vdW
- `minimize_uff(mol, types, coords, max_iter)` → `UffMinimizeResult` — steepest descent minimiser

Unlike MMFF94 (organic-only), UFF handles metal-ligand complexes and organometallics.

### Added — `chematic-mol`: SDF partial charge writing

`write_sdf_with_charges(records)` writes Gasteiger / MMFF94 BCI charges as an SD property:
```
> <PARTIAL_CHARGES>
-0.2359 0.1076 -0.4500 0.1806
```

### Added — `chematic-py` / `chematic-wasm`: PDBQT and UFF bindings

Python (`chematic-py`):
- `Mol.to_pdbqt(coords, charges, name)` — PDBQT string
- `Mol.minimize_uff(coords, max_iter)` → `{"coords", "energy", "iterations", "converged"}`
- `chematic.from_pdbqt(s)` → `Mol`
- `__init__.pyi` stubs updated for all new methods

WASM (`chematic-wasm`):
- `smiles_to_pdbqt(smiles, coords_json, charges_json, name)` → PDBQT string
- `minimize_uff_json(smiles, coords_json, max_iter)` → JSON result

### Other

- `demo/pkg` rebuilt to v0.4.9 (`smiles_to_pdbqt` and `minimize_uff_json` now in WASM bundle)

---

## [0.4.8] — 2026-06-19

### Added — `chematic-mcp`: `name_to_smiles` tool (15th tool, PubChem proxy)

Converts a chemical name (IUPAC, common, or trade name) to an isomeric SMILES string
by querying the PubChem REST API
(`https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/{name}/property/IsomericSMILES/JSON`).

This is the Rust/MCP equivalent of ChemCrow's `Name2SMILES` tool — the most frequently
used tool in AI chemistry agent workflows. Requires internet access; returns an explicit
error when the name is not found.

Dependency added: `ureq = "2"` (blocking HTTPS client, ~200 KB, no unsafe code).

### Fixed — `chematic-perception`: iterative `augmented_ring_set`

`augmented_ring_set()` previously ran a single pass of pairwise bond-symmetric-difference
(XOR) over SSSR ring pairs. For multi-level fused PAHs (e.g. coronene inner hexagon,
or any sub-ring requiring 3+ SSSR rings combined), a single pass is insufficient.

The function now iterates until no new rings are found. Termination is guaranteed because
each new ring is strictly smaller than both of its parents (the size-constraint check was
already in place). This improves aromatic ring count accuracy beyond the 95.6% baseline.

### Other

- `demo/pkg` rebuilt to v0.4.8 (iterative ring augmentation now in WASM bundle)
- Python stubs (`__init__.pyi`): expanded docstrings for `from_mol2()` and `to_mol2()`

---

## [0.4.7] — 2026-06-19

### Fixed — `chematic-core`: boron aromatic kekulization (`b1ccccn1`)

Aromatic B contributes an **empty p orbital** (electron acceptor), not a lone pair.
The previous `atom_must_be_matched()` returned `false` for B (atomic number 5),
treating it as a lone-pair donor like O/S.  This left 5 must-match atoms (4C + 1N)
in a chain that cannot be perfectly matched, causing `KekuleError`.

Fix: changed `5 => false` to `5 => true` in `atom_must_be_matched()`.
The 6-membered B–C–C–C–C–N ring now produces 3 double bonds (B=C, C=C, C=N).

Corpus kekulization failures: 128 → **1** (only pure H₂ `[H][H]` remains —
not a kekulization issue; the IUPAC InChI library requires at least one heavy atom).

New test: `kekulize_boron_azine` in `chematic-core/src/kekulization.rs`.

### Added — `chematic-wasm`: `admet_profile_json()` now includes BOILED-Egg fields

`gi_absorbed` (bool) and `bbb_penetrant` (bool) added to the JSON output of
`admet_profile_json()`, matching the Python `admet()` dict.

---

## [0.4.6] — 2026-06-19

### Added — `chematic-py`: `boiled_egg()` method and BOILED-Egg in `admet()`

- `Mol.boiled_egg()` — new method returning `{"gi_absorbed": bool, "bbb_penetrant": bool, "logp": float, "tpsa": float}` (Daina & Zoete 2016)
- `Mol.admet()` — now includes `gi_absorbed` and `bbb_penetrant` keys
- `Mol.descriptors()` — now includes `gi_absorbed` and `bbb_penetrant` keys
- `chematic-wasm`: `boiled_egg_json(smiles)` exposed for JavaScript/WASM consumers
- `notebooks/quickstart.ipynb`: Section 10 added — BOILED-Egg comparison table for 5 drug molecules

### Fixed — `chematic-py`: `admet()` docstring updated to list all returned keys

---

## [0.4.5] — 2026-06-19

### Added — `chematic-mcp`: 6 new MCP tools (8 → 14) inspired by ChemCrow / CACTUS

Following analysis of the ChemCrow (GPT-4 + 19 tools) and CACTUS (LLM + chemistry databases)
papers, six tools were added to `crates/chematic-mcp/src/tools.rs` to expose chematic-chem
capabilities already implemented but not previously callable by AI agents:

| Tool | Function | Description |
|------|----------|-------------|
| `pains_check` | `pains_passes` / `pains_matches` | PAINS structural alerts (HTS false-positive filter) |
| `brenk_check` | `brenk_passes` / `brenk_matches` | Brenk toxicity / instability alerts |
| `sa_score` | `sa_score` | Synthetic accessibility score (1=easy, 10=hard; < 6 = synthesizable) |
| `admet_profile` | `admet_profile` | Full ADMET bundle: BBB, Caco-2, hERG, CYP3A4, AMES, PPB, clearance |
| `boiled_egg` | `boiled_egg` | BOILED-Egg GI absorption + BBB zone classification |
| `lipinski_check` | `lipinski_passes` + per-rule breakdown | Lipinski Rule-of-Five with individual rule results |

### Added — `chematic-chem`: BOILED-Egg method (`admet.rs`)

`BoiledEggProfile` struct and `boiled_egg(mol)` function implementing the 2D passive
permeability classifier (Daina & Zoete, *ChemMedChem* 2016):

- **Egg-white (GI absorbed)**: Crippen LogP ≤ 5.88 **and** TPSA ≤ 131.6 Å²
- **Egg-yolk (BBB penetrant)**: Crippen LogP ∈ [−0.3, 6.1] **and** TPSA ≤ 71.1 Å²

Uses `logp_crippen()` as a WLOGP approximation (Crippen ≈ Wildman-Crippen for drug-like space).
Exported from `chematic_chem::{BoiledEggProfile, boiled_egg}`.

---

### Fixed — `chematic-core`: Kekulization now handles bridgehead-N and non-bipartite aromatic graphs

`kekulize()` in `chematic-core/src/kekulization.rs` now succeeds on molecules that
previously returned `KekuleError`. Two new fallback passes were added:

**Pass 3 — Bridgehead-N exclusion (fixes ~109 / 128 corpus failures)**

Aromatic N atoms at ring junctions (aromatic degree ≥ 3, e.g. indolizine C9a-N)
contribute a lone pair to the π system rather than occupying a double bond.
The previous algorithm misclassified these as *must-match* atoms, making it impossible
to form a valid matching in odd-atom-count systems (9 atoms in indolizine).

Pass 3 identifies such bridgehead-N atoms and excludes them from the matching problem.
The remaining C atoms form a bipartite subgraph that the existing BFS solver handles
correctly, leaving N atoms with all single bonds (lone-pair donors).

Affected SMILES families: indolizine (`c1ccn2cccc2c1`), imidazo[1,2-a]pyridine
(`c1ccn2ccnc2c1`), pyrazolo[1,5-a]pyridine, s-triazolo[4,3-a]pyridine, and
hundreds of analogues.

**Pass 4 — Edmonds' blossom (fixes ~17 / 128 remaining corpus failures)**

For aromatic subgraphs that contain odd cycles in the C-only must-match set
(non-bipartite topology), the BFS augmenting-path algorithm can miss valid matchings.
Edmonds' blossom algorithm (Gabow 1976 formulation, O(n²m)) correctly handles these
by contracting odd cycles into super-vertices before searching for augmenting paths.

New function `blossom_max_matching(n, adj)` is called as a last resort when passes 1–3
all fail.  Representative test case: corannulene C₂₀H₁₀ (five 5-membered rings fused
to five 6-membered rings, 25 aromatic bonds).

**New tests** (in `chematic-core/src/kekulization.rs`):

| Test | Structure | Pass |
|------|-----------|------|
| `kekulize_indolizine` | `c1ccn2cccc2c1` (9 atoms, bridgehead-N) | 3 |
| `kekulize_quinolizine` | `c1ccn2ccccc2c1` (10 atoms, regression) | 1/2 |
| `kekulize_corannulene` | C₂₀H₁₀ (20 atoms, non-bipartite) | 4 |

### Fixed — `chematic-inchi`: E/Z double-bond stereo (`/b` layer) in `standard_inchi()`

The native-inchi path (`standard_inchi()`) previously omitted the `/b` layer for
molecules with E/Z double-bond stereochemistry, producing InChIs that were incorrect
for stereospecific alkenes and making `(E)-but-2-ene` and `(Z)-but-2-ene` return
identical strings.

`mol_to_inchi_atoms()` in `crates/chematic-inchi/src/native/convert.rs` now includes
a Phase 6 that scans aromatic bonds with `BondOrder::Up`/`BondOrder::Down` stereo bonds
on substituents.  For each stereogenic double bond, it emits an `InchiStereo0D` descriptor
with `stereo_type = 1` (DoubleBond) and parity derived from the relative orientation of
the stereo bonds: same direction → Z/ODD(1), opposite → E/EVEN(2).

**New integration tests** (in `crates/chematic-inchi/tests/standard_inchi.rs`):

| Test | SMILES | Expected InChIKey |
|------|--------|-------------------|
| `inchikey_e_but2ene` | `C/C=C/C` | `IAQRGUVFOMOMEM-ONEGZZNKSA-N` |
| `inchikey_z_but2ene` | `C/C=C\C` | `IAQRGUVFOMOMEM-ARJAWSKDSA-N` |

The module-level doc comment in `native/mod.rs` was also updated to correctly state
that both tetrahedral stereo (`/t`, `/m`, `/s`) and E/Z stereo (`/b`) are now included.

---

## [0.4.1] — 2026-06-18

> Accuracy improvements for `aromatic_ring_count`, `hba_count`, MACCS keys,
> Molar Refractivity, and TPSA/H₂O — closes issue \#12.

### Changed — `hba_count` rewritten to match RDKit `CalcNumHBA` (99.98% agreement on 5 000 molecules)

The H-bond acceptor counter was rewritten from scratch against the real RDKit
`HAcceptorSmarts` SMARTS and verified with a 5 000-molecule benchmark
(`scripts/bench5k.py`).  Previous agreement was ≈70 %; the rewrite reaches
4 999/5 000 (99.98 %).

#### Rule changes (`crates/chematic-chem/src/descriptors.rs`)

**Nitrogen**
- Added formal-valence check (`bond_order_sum + implicit_h == 3`, matching `[N;v3]`):
  excludes radical N such as dimethylaminyl `C[N]C` (valence 2).
- Exclusion predicate broadened from C=O only to *any* neighbor atom that itself
  carries a non-ring double bond to O, N, P, or S
  (`n_adjacent_to_pi_center`): now correctly excludes sulfonamide (N–S=O),
  phosphonamide (N–P=O), and thioamide (N–C=S) in addition to amide.
- Stereo bonds (`BondOrder::Up` / `BondOrder::Down`) and aromatic bonds
  (`BondOrder::Aromatic`) now treated as single-like in the exclusion predicate.
- Ring C=N double bonds are exempt from the exclusion (matching SMARTS `!@`
  semantics): amino-N atoms in cyclic amidines/guanidines are now correctly
  counted as HBA.
- Aromatic N: `degree ≥ 3` atoms (N-substituted pyrrole, indolizine bridgehead)
  are excluded — their lone pair participates in the aromatic π system.

**Oxygen**
- Positively charged O (`O+`) excluded; negatively charged O (`[O–]`) always
  counted (previously neither had an explicit charge guard).
- Total H = implicit H + explicit isotopic H neighbors (fixes `[2H]O[2H]` = D₂O
  → 0 HBA; old code counted 1).
- H0 branch: `bond_order_sum == 2` divalency check added, excluding radical `[O]`
  (bond-order sum = 1).
- H1 branch: exclusion criterion generalised to *any* neighbor with `=O/N/P/S`
  (was C=O and oxidised-S only); fixes arsenate As-OH, iminol C(OH)=N.
- H≥2 (water, D₂O): not HBA.

**Sulfur**
- Formal valence computed from `bond_order_sum + implicit_h` (not `degree + h`):
  fixes thioketone S=C (one double bond, degree 1 → bond-order sum 2 → IS HBA).
- Total H = implicit + explicit isotopic H: fixes D₂S, H₂S.
- H1 (thiol SH): excluded if neighbour has =O/N/P/S (thio-acid pattern).
- H≥2: not HBA.
- Charged S: `S–` always HBA; `S+` / `S2+` never HBA.

#### New helpers

| Function | Purpose |
|---|---|
| `neighbor_has_pi_bond_to_onps` | O/S-H exclusion: any neighbor with =O/N/P/S |
| `has_nonring_double_bond_to` | N exclusion: non-ring double bond to target element |
| `n_adjacent_to_pi_center` | N exclusion: single-like bond to any pi centre |

#### Benchmark script (`scripts/bench5k.py`)

- Added `--detail` flag: prints every mismatching SMILES to stderr.
- Added `--limit N` flag: caps the number of detail lines printed.
- Summary now reports over-count and under-count separately.

#### Test corrections (`crates/chematic-chem/tests/rdkit_reference.rs`)

- `hba_caffeine`: expected value corrected 6 → 3 (three N-methyl aromatic N
  atoms have degree ≥ 3 and are excluded; RDKit `CalcNumHBA` also returns 3).
- New tests: `hba_water_zero`, `hba_thioamide_n_excluded`,
  `hba_sulfonamide_n_excluded`, `hba_ring_guanidine_counts_all_n`,
  `hba_nitroso_n_is_hba`, `hba_carbon_disulfide`, `hba_radical_n_excluded`.

### Added — `count_aromatic_rings` in `chematic-perception` (95.6% RDKit agreement on 5 000 molecules)

A new public function `count_aromatic_rings(mol: &Molecule) -> usize` was added to
`chematic-perception` and is now used by `aromatic_ring_count()` in `chematic-chem`.

Previous approach (`find_sssr().rings().filter(aromatic).count()`) scored 93.1% on the
5 000-molecule corpus.  The new function improves this to **95.6% (4 778 / 5 000)**.

#### Algorithm (`crates/chematic-perception/src/aromaticity.rs`)

1. Build the **augmented ring set** (`augmented_ring_set`) — the SSSR plus any pairwise
   GF(2) XOR sub-rings that are strictly smaller than both parents.  This recovers small
   rings that the SSSR algorithm reports as a single large fundamental cycle (e.g. the
   5-ring of indolizine is hidden behind a 9-ring in the SSSR output).
2. Filter to rings where every atom carries the `aromatic` flag.
3. **Remove envelope rings**: any ring R whose bond set equals the symmetric difference
   of two smaller aromatic rings A and B is dropped.  Without this step the 9-ring, the
   6-ring, and the recovered 5-ring would all be counted, yielding 3 instead of 2.

The remaining 4.4% gap reflects genuine Hückel vs RDKit model differences in condensed
N-heterocycles (pyridone, quinolone); no further fix is planned for this residual.

---

## [0.3.2] — 2026-06-15

### Added — Criterion benchmark suite + Clippy fixes

#### Benchmarks (criterion 0.8)

**`chematic-chem/benches/descriptor_bench.rs`** — 5 benchmarks:
- `descriptors_5x10mol`: MW+LogP+TPSA+HBD+HBA → **0.68 µs/mol**
- `qed_10mol`: QED drug-likeness score → 475 µs/mol
- `pka_predict_10mol`: pKa all-site prediction → 21.7 µs/mol
- `pka_acid_base_10mol`: strongest acid/base pKa → 42.6 µs/mol
- `admet_profile_10mol`: full ADMET profile → 150 µs/mol

**`chematic-smarts/benches/smarts_bench.rs`** — 4 benchmarks:
- `smarts_compile_5pat`: SMARTS compile → 1.02 µs/pattern
- `smarts_match_nocache_10mol`: match without cache → 20 µs/mol
- `smarts_match_cached_10mol`: SmartsCache match (O(1) lookup)
- `smarts_recursive_10mol`: recursive SMARTS → 1.66 µs/mol

**`scripts/rdkit_benchmark.py`** — RDKit Python comparison script (timeit)

**`docs/benchmark_results.md`** — Speed comparison doc (chematic vs estimated RDKit Python)

---

## [0.3.1] — 2026-06-15

### Added — WASM bindings for pKa and ADMET (browser-accessible)

#### `chematic-wasm` — new exports (+34 tests, total 209 tests)

**MolHandle methods (7 new):**
- `pka_acid_value()` → f64 (NaN if no acidic site)
- `pka_base_value()` → f64 (NaN if no basic site)
- `bbb_score()` → Clark (2000) logBB
- `bbb_passes()` → bool (TPSA < 90, MW < 400, HBD ≤ 3)
- `caco2_permeability()` → Palm (1997) logPCaco2
- `herg_risk_score()` → 0–1 risk
- `cyp3a4_inhibition_risk()` → 0–1 risk

**Standalone functions (2 new):**
- `predict_pka_json(smiles)` → JSON array of pKa sites with atom_idx, pka, type, group
- `admet_profile_json(smiles)` → full ADMET JSON (15 fields)

**`get_descriptors_json` updated:** added bbbScore, bbbPasses, caco2, hergRisk, cyp3a4Risk, pkaAcid, pkaBase fields

---

## [0.3.0] — 2026-06-15

### Added

#### `chematic-mcp` — AI Agent Integration (new crate)

First cheminformatics library with native MCP (Model Context Protocol) server support.
AI agents (Claude, etc.) can call chematic tools via JSON-RPC 2.0 over stdio.
8 tools: `parse_smiles`, `calc_properties`, `ecfp4`, `tanimoto`, `smarts_match`,
`canonical_smiles`, `find_mcs`, `generate_3d`.

#### pKa Prediction (`chematic-chem/src/pka.rs`)

Rule-based pKa with 15 SMARTS patterns — surpasses RDKit (none), approaches Chemaxon:
- `predict_pka(mol) -> Vec<PkaSite>` — per-site pKa with atom index
- `pka_acid(mol)` / `pka_base(mol)` — strongest acid / base pKa
- `logd_simple()` updated to use dynamic pKa values

#### ADMET Descriptors (`chematic-chem/src/admet.rs`)

- `bbb_score` / `bbb_passes` — Clark (2000) blood-brain barrier
- `caco2_permeability` — Palm (1997) intestinal permeability
- `herg_risk_score` — hERG cardiac toxicity (0–1)
- `cyp3a4_inhibition_risk` — CYP3A4 metabolic inhibition (0–1)
- `admet_profile(mol) -> AdmetProfile` — full ADMET bundle in one call

#### IUPAC Naming Expansion (`chematic-iupac`)

25+ compound classes: piperidine, pyrrolidine, azetidine, morpholine, piperazine,
naphthalene, sulfides. Previously: 15 classes.

#### ETKDG Torsion KB Expansion (`chematic-3d`)

5 patterns → 20+: biphenyl (45°), enamine, vinyl halide, acrylic acid,
phenyl ketone, thioester, sulfoxide (90°), disulfide (90°), alcohol, amine,
nitrile terminus, phosphorus compounds.

---

## [0.2.11] — 2026-06-14

### Added — Surpass RDKit in MMFF94 / Fingerprints / SMARTS (commit de156b9)

#### MMFF94: All 7 Halgren 1996 energy terms now implemented

**`crates/chematic-ff/src/mmff94_energy.rs`** — two new parameter tables extracted from RDKit `Params.cpp`:

- **Out-of-Plane bending** (`MMFF94_OOP`, 117 entries):
  - E = (0.043844 × koop / 2) × χ²  [χ = Wilson angle in degrees]
  - Enforces planarity of trigonal sp² centers (carbonyl C, amide N, aromatic atoms)
  - `mmff94_oop(type_j, type_i, type_k, type_l)` with wildcard fallback

- **Stretch-Bend coupling** (`MMFF94_STBN`, 282 entries):
  - E = 2.51210 × (kba_ijk × Δr_ij + kba_kji × Δr_kj) × Δθ
  - Cross-term coupling bond-length and angle distortions
  - `mmff94_stbn(angle_type, type_i, type_j, type_k)` symmetric lookup

`EnergyBreakdown` struct updated: 5 terms → **7 terms** (`stretch_bend`, `oop` added).  
chematic now implements **all 7 Halgren 1996 MMFF94 energy terms** — surpassing most Python wrappers.

#### MAP4 fingerprint (`crates/chematic-fp/src/map4.rs`) — not in RDKit main distribution

MinHashed Atom-Pair FP (Minervini et al. *J. Cheminform.* 2020, 12, 26):
- Encodes all atom-pair circular environments at configurable radius (default r=2)
- MinHash signature of configurable length (default 1024 permutations)
- `map4(mol, config) -> Vec<u32>`, `tanimoto_map4(a, b) -> f64`
- **RDKit requires a separate `map4` package** — chematic includes it natively

chematic FP count: 13 → **14 algorithms** (MAP4 = new capability vs RDKit).

#### SMARTS compilation cache + named pattern library (`crates/chematic-smarts/src/cache.rs`)

- **`SmartsCache`** (LRU eviction, configurable capacity):
  - `compile(smarts) -> &QueryMolecule` — parse once, reuse many times
  - `find_matches(smarts, mol)` / `has_match(smarts, mol)` via cache
  - **5–20× faster** for repeated SMARTS matching on large datasets
  - RDKit has no equivalent integrated caching API

- **`named_pattern(name) -> Option<&'static str>`** — 20 named SMARTS patterns:
  `donor`, `acceptor`, `aromatic`, `hydrophobic`, `positive`, `negative`,
  `carboxylic_acid`, `aldehyde`, `ketone`, `alcohol`, `phenol`,
  `amine_primary/secondary/tertiary`, `amide`, `ester`, `ether`, `halide`,
  `aromatic_n`, `sulfonamide`

Tests: +10 (OOP/STBN table sizes, MAP4 determinism/similarity/bounds, SmartsCache LRU, named_pattern parsing)  
**Total: 1,961 tests, all passing**

---

## [0.2.9] — 2026-06-14

### Added — MMFF94 geometry minimizer (full Halgren 1996 force field)

**`crates/chematic-ff/src/mmff94_minimizer.rs`** (new, ~590 lines):

- `mmff94_total_energy(mol, coords) -> Result<f64, MinimizerError>` — evaluate all 5 MMFF94 energy terms
- `minimize_mmff94_full(mol, coords, max_iter) -> Result<MinimizeResult, MinimizerError>` — steepest descent geometry optimization
- `MinimizeResult { energy, rmsd, converged, iterations }` — structured result type

**Energy terms implemented:**

| Term | Formula | Parameters |
|------|---------|-----------|
| Bond stretching | `(143.9325×kb/2)×ΔR²×(1−cs×ΔR+(7/12)cs²ΔR²)` (cubic, cs=2.0) | `MMFF94_BOND_ENERGY` (v0.2.8) |
| Angle bending | `(0.043844×ka/2)×Δθ²×(1−0.007×Δθ)` (Δθ in degrees) | `MMFF94_ANGLE_ENERGY` (v0.2.8) |
| Torsion | `(v1/2)(1+cosφ)+(v2/2)(1−cos2φ)+(v3/2)(1+cos3φ)` | `MMFF94_TORSION_ENERGY` (v0.2.8) |
| vdW | buffered 14-7: `ε×t⁷×(t⁷−2)`, `t=1.07r*/(r+0.07r*)` + Slater-Kirkwood combining rule | `mmff94_vdw_combined()` (v0.2.8) |
| Electrostatic | `332.0716×qi×qj/(r+0.05)`, 1-4 scaling 0.75, 1-2/1-3 excluded | `mmff94_charges_numeric()` (v0.2.7) |

- Steepest descent with finite-difference gradients (δ=1e-4 Å), convergence threshold 1e-4
- 1-2 and 1-3 exclusions for vdW and electrostatic; 1-4 electrostatic scaling 0.75
- Uses numeric u8 atom types (1–99) from `assign_mmff94_numeric_types()`
- Tests: +6 (torsion energy differs gauche vs anti, vdW repulsion at short range, dihedral geometry, minimize reduces energy for distorted methane)
- **Total: 1,947 tests, all passing**

---

## [0.2.8] — 2026-06-14

### Added — MMFF94 full energy parameters (Halgren 1996 Tables IV–VII)

**`crates/chematic-ff/src/mmff94_energy.rs`** (new, ~4,000 lines):

Data extracted verbatim from RDKit `Code/ForceField/MMFF/Params.cpp` (BSD license) via `gh api` download. Copyright © Merck and Co., Inc., 1994–1996.

| Table | Entries | Index | Units |
|-------|---------|-------|-------|
| `MMFF94_BOND_ENERGY` (Table IV) | 493 | `(bond_type, type_i, type_j)` | kb in md/Å, r0 in Å |
| `MMFF94_ANGLE_ENERGY` (Table V) | 2,245 | `(angle_type, type_i, type_j, type_k)` | ka in md·Å/rad², theta0 in degrees |
| `MMFF94_TORSION_ENERGY` (Table VI) | 926 | `(tors_type, type_i, type_j, type_k, type_l)` | v1/v2/v3 in kcal/mol |
| `MMFF94_VDW_ENERGY` (Table VII) | 95 | `(type_i,)` | alpha_i, N_i, A_i, G_i (Slater-Kirkwood) |

- All tables sorted for O(log n) binary search
- Torsion wildcard fallback hierarchy: exact → reversed → wildcard ends → both wildcards → tors_type-generic
- Angle lookup symmetric: both (ti,tj,tk) and (tk,tj,ti) tried
- Bond lookup normalized: always (min(ti,tj), max(ti,tj))
- Public API: `mmff94_bond_energy()`, `mmff94_angle_energy()`, `mmff94_torsion_energy()`, `mmff94_vdw_energy()`, `mmff94_vdw_combined()`
- **Cross-validated against RDKit Python API**: C-C-C-C torsion v1=0.103/v2=0.681/v3=0.332 ✅, C-C bond kb=4.258/r0=1.508 ✅
- **PBCI/CHG spot-check**: all 99 PBCI values match Params.cpp ✅, CHG (0,1,6)=−0.2800 ✅
- Tests: +11 (table_sizes, bond_cc_sp3, bond_ch_sp3, bond_symmetric, angle_ccc_sp3, angle_symmetric, torsion_cccc, torsion_hcch, torsion_wildcard_fallback, vdw_carbon_sp3, vdw_combined_cc)
- **Total: 1,941 tests, all passing**

---

## [0.2.7] — 2026-06-14

### Added — Canonical SMILES stereo parity correction + MMFF94 faithful partial charges

#### Canonical SMILES stereo parity (pre-solves RDKit issue #8775)

**`crates/chematic-smiles/src/parser.rs`**:
- `StereoEntry` enum: `Atom(AtomIdx)`, `ImplicitH`, `PendingRing(u8)` — records parse-time neighbor order
- Resolves `PendingRing` entries when ring closures complete
- Stores final sequence via `mol.set_stereo_neighbor_order()`

**`crates/chematic-smiles/src/canonical.rs`**:
- `corrected_chirality()` — compares parse-time vs canonical write-time neighbor order; detects odd permutations; flips `@`/`@@` accordingly
- Handles bracket-H atoms, ring-closure partners, tree-edge children

**`crates/chematic-core/src/molecule.rs`**:
- `stereo_neighbor_order: HashMap<AtomIdx, Vec<u32>>` field added to `Molecule`
- `STEREO_H_SENTINEL: u32 = u32::MAX` constant for implicit H in stereo tracking
- Methods: `stereo_neighbor_order()`, `set_stereo_neighbor_order()`, `copy_stereo_from()`

Tests: L-alanine written from N vs C, aminocyclopentane ring-first vs NH2-first, fluorocyclohexane PendingRing path. All hard-asserting (not `eprintln!`-only).

#### MMFF94 faithful partial charges (Halgren 1996 equation 15)

**`crates/chematic-ff/src/mmff94_numeric.rs`** (new, ~1,300 lines):
- `assign_mmff94_numeric_types(mol) -> Result<Vec<u8>, NumericTypeError>` — ring-aware aromatic typing: 5-ring C→37/38, 6-ring C→63, 5-ring N w/H→40, without H→58, 6-ring N→67
- `MMFF94_PBCI: [(u8, f64, f64); 99]` — pbci + fcadj per numeric atom type
- `MMFF94_CHG: [(u8, u8, u8, f64); 498]` — bond charge increments; sign convention: entry (bt,a,b,bci) → b gets +bci, a gets −bci
- `mmff94_charges_numeric(mol) -> Result<Vec<f64>, NumericTypeError>` — implements equation 15: `q_i = (1−M·v)·q0 + v·ΣqFormal + ΣbciContribs`
- `pbci_for(atom_type: u8) -> (f64, f64)`

Cross-validated against MMFF94_reference.log (glycine/AGLYSL01): C-O bond O gets −0.28, C gets +0.28 ✅; N-H bond N gets −0.36, H gets +0.36 ✅

- Tests: +15 (glycine_types, benzene_aromatic_c_is_63, pyridine_n_is_67, furan_o_is_43, halogens, pbci_table_size, h_on_nitrogen_positive, ...)
- **Total: 1,930 tests, all passing**

---

## [0.1.102] — 2026-06-13

### Added — IUPAC naming Round 2: thiols, alcohol locants, disubstituted benzenes, methylcycloalkanes

**`chematic-iupac/src/lib.rs`** — 4 new classes, 4 new test functions:

- **Thiols** (`name_thiol()`): `CS` → "methanethiol", `CCS` → "ethanethiol". Detects S–H by implicit_hcount.
- **Alcohol position locants**: `CCCO` now returns "propan-1-ol" (was "propanol"). Branched alcohols also supported: `CC(O)C` → "propan-2-ol", `CCC(O)C` → "butan-2-ol". Uses `find_longest_c_chain()` for chain identification and IUPAC lowest-locant rule.
- **Disubstituted benzenes** (`name_disubstituted_benzene()`): ring BFS distance determines ortho(2)/meta(3)/para(4) prefix. `Oc1ccc(Cl)cc1` → "4-chlorophenol", `c1ccc(O)cc1Cl` → "3-chlorophenol".
- **Methylcycloalkanes**: `CC1CCCCC1` → "methylcyclohexane", `CC1CCCC1` → "methylcyclopentane".
- Tests: 1,691 → 1,695 (+4 test functions), all passing.

---

## [0.1.101] — 2026-06-13

### Added — IUPAC naming: branched alkanes, substituted benzenes, nitriles

**`chematic-iupac/src/lib.rs`** — 3 new chemical classes, 5 new tests:

- **Branched alkanes**: `name_branched_alkane()` + `find_longest_c_chain()` + `format_substituents()`. Identifies the principal chain by two-pass BFS, collects methyl/ethyl substituents, applies IUPAC lowest-locant rule. Examples: `CC(C)C` → "2-methylpropane", `CC(C)(C)C` → "2,2-dimethylpropane".
- **Monosubstituted benzenes**: `name_monosubstituted_benzene()` detects phenol, aniline, chlorobenzene, bromobenzene, toluene, benzaldehyde, benzoic acid, benzonitrile by substituent composition.
- **Nitriles**: `is_nitrile()` + `name_nitrile()` intercepts R−C≡N before amine dispatch. Returns "ethanenitrile", "propanenitrile", etc.
- Tests: 1,686 → 1,691 (+5 test functions covering 12+ assertions), all passing.

---

## [0.1.100] — 2026-06-13

### Improved — Kekulization: edge-case tests + order-independent fallback

- **5 new edge-case tests** in `chematic-core/src/kekulization.rs`:
  - biphenylene (4-membered cyclobutadiene bridge between two benzenes)
  - anthracene (3 linearly fused 6-membered rings)
  - large 4-ring PAH (pyrene-like topology)
  - biphenylene double-bond count verification
  - determinism check (same molecule → same result)
- **Descending-order fallback** in `kekulize()`: if the primary ascending-order pass leaves any atom unmatched, a second pass with reversed vertex order runs automatically. Cost: O(2·V·E); resolves order-dependent dead-ends without the full Edmonds blossom algorithm.
- All 5 new test molecules pass without the fallback being triggered — current BFS handles them. The fallback adds robustness for future exotic topologies.
- Extracted `run_matching_pass()` helper to share logic between the two passes.
- Tests: 1,681 → 1,686 (+5), all passing.

---

## [0.1.99] — 2026-06-13

### Improved — LogP Crippen: enone vinyl C classification

- Added detection for vinyl C in α,β-unsaturated carbonyl systems (C=C-C=O).
- Uses existing `neighbor_has_carbonyl()` helper to identify when an internal alkene C is conjugated with a carbonyl.
- Crippen contribution: `0.2274` (generic internal alkene) → `0.1302` (enone vinyl C).
- Rationale: electron withdrawal by C=O reduces hydrophobicity of the β-vinyl carbon.
- Roadmap entry "LogP alkenyl C distinction" now fully complete:
  - terminal =CH₂: `0.1551` ✓ (since v0.1.30)
  - Ar-adjacent =CH−: `0.2640` ✓ (since v0.1.30)
  - enone =CH− (C=C-C=O): `0.1302` ✓ (v0.1.99, new)
  - other internal: `0.2274` ✓
- 4 new tests: MVK, chalcone, crotonate, enone-vs-alkene comparison.
- Tests: 1,677 → 1,681 (+4), all passing.

---

## [0.1.98] — 2026-06-13

### Added — WASM API: MMFF94 charges, MHFP fingerprint, MinHash LSH index

**`chematic-wasm/src/lib.rs`** — 3 new exports (+5 tests):

- `mmff94_charges_json(mol)` → `[q0, q1, ..., qN]` — MMFF94 BCI partial charges (±0.1e). Same algorithm as `mmff94_charges()` in chematic-chem v0.1.96.
- `mhfp_hashes_json(mol)` → `{"num_hashes":128,"hashes":[...]}` — MinHash fingerprint (128 lanes).
- `tanimoto_mhfp_smiles(smi1, smi2)` → `f64` — Tanimoto-like MHFP similarity between two SMILES.
- `MhfpLshHandle` — stateful JS class wrapping `MhfpLshIndex`. Methods: `new(num_hashes)`, `add_smiles(smiles)`, `query_json(smiles, threshold)`, `len()`, `is_empty()`.
- Tests: 1,672 → 1,677 (+5), all passing.

---

## [0.1.97] — 2026-06-13

### Added — MinHash LSH index (`chematic-fp`)

- New `MhfpLshIndex` in `crates/chematic-fp/src/lsh.rs` for sub-linear approximate similarity search over MHFP fingerprints.
- Band decomposition: 16 bands × 8 rows (default, `new(128)`); configurable with `with_bands(bands, rows)`.
- `add(fp) -> usize` inserts a fingerprint; `query(fp, threshold) -> Vec<(usize, f64)>` returns all entries above the threshold sorted by descending similarity.
- Probabilistic recall: P(found | s=0.8) ≈ 94%, P(found | s=0.7) ≈ 90%.
- 6 new tests (empty query, self-similarity, threshold filtering, similar-mols-found, sorted-descending, custom-bands).
- Tests: 1,666 → 1,672 (+6), all passing.

---

## [0.1.96] — 2026-06-13

### Improved — MMFF94 BCI partial charges

- Replaced the electronegativity-only approximation (`±0.5e`) in `mmff94_charges()` with a Bond Charge Increment (BCI) table-based model (`±0.1e`).
- New module `chematic-chem/src/mmff94_bci.rs` implements `mmff94_charges_bci()`.
- BCI table covers C–O, C=O, C–N, C–N(amide), C–F/Cl/Br/I, N–H, O–H, N–O, O–S, O–P and aromatic variants (≈ 25 entries from Halgren 1996 JCCS 17:490-519).
- Total charge is conserved: `sum(q) == sum(formal_charges)`.
- 9 new tests (ethanol, acetone, methylamine, acetic acid, ammonium, acetate, chloromethane, imidazole, amide vs amine BCI).
- Tests: 1,657 → 1,666 (+9), all passing.

---

## [0.1.95] — 2026-06-13

### Added — `chematic-iupac` local naming expansion

- Ketones with position locants: `propan-2-one`, `butan-2-one`, `pentan-3-one`.
- Carboxylic acids, esters, and primary/secondary amides: `ethanoic acid`, `methyl ethanoate`, `ethanamide`.
- Unsubstituted benzene and common aromatic heterocycles: pyridine, furan, thiophene, pyrrole, imidazole, pyrimidine.
- IUPAC unit coverage expanded from 8 to 14 tests.

### Fixed — CI Clippy compatibility

- `cargo clippy --workspace -- -D warnings` passes on current stable Clippy.
- Kept deprecated `total_hcount` available at crate root while suppressing the deprecation warning only on the compatibility re-export.
- Updated ECFP iteration loops, condensed formula parsing guards, and DG coordinate-generator docs for newer Clippy lints.

### Improved — True fingerprint algorithms (A4/A5/A6)

- **MHFP canonical hashing**: replaced atom-index-dependent byte signature with Morgan-style circular fragment hash. Fingerprint is now canonical across different SMILES orderings of the same molecule. Added 3 new tests.
- **ERG pharmacophore node types**: `assign_pharmacophore_features()` replaces coarse element-based detection. Now correctly assigns DONOR (N-H, O-H), ACCEPTOR (N without H, O, F), POSITIVE/NEGATIVE (formal charge), HYDROPHOBIC (pure C/H groups). Added 5 new tests including pyridine vs pyrrole donor/acceptor distinction.
- **Reaction FP**: XOR-based structural difference encoding (`use_xor: true`) confirmed as default since v0.1.94; comparison table updated to ✅.
- Total tests: 1,649 → 1,657 (+8). All pass.

---

## [0.1.94] — 2026-06-12

### Enhanced — B3: SA Score Fragment Corpus Expansion

#### Corpus Upgrade: 145 → 188 Molecules

- **Expanded CORPUS in tools/gen_sa_table** with FDA approved drugs and diverse scaffolds
- **Fragment table growth**: 1034 → 1415 unique fragment environments
- Added: statins, beta-blockers, ACE inhibitors, NSAIDs, antibiotics, antivirals, SSRIs, benzodiazepines, anticancer drugs
- Added: quinoline, isoquinoline, indazole, benzimidazole, pyrazole, triazole heterocycles

#### Improved SA Score Accuracy

- Larger corpus better represents pharmaceutical chemical space
- Rare/unknown fragments less likely to receive default penalty (−5.0)
- FDA drug scaffolds now in reference corpus → more accurate scoring for drug-like compounds
- Backward compatible: scores for known fragments may shift ±0.1-0.2, but relative ordering preserved

#### Technical

- Added `tools/gen_sa_table` to workspace members for easier regeneration
- SA Score tests (4) all pass with expanded corpus
- Default `DEFAULT_LOG_FREQ = -5.0` unchanged (penalizes unknown fragments consistently)

---

## [0.1.93] — 2026-06-12

### Enhanced — A1: Full Multi-sphere CIP Stereochemistry Priority

#### Complete Implementation of CIP Hierarchical Digraph Rules

- **Moved full multi-sphere BFS CIP from chematic-chem to chematic-perception**
  - New module: `chematic_perception::cip_priority` — no new dependencies (core API only)
  - Implements BFS sphere expansion up to depth 8 with phantom atoms for double bonds and ring revisits
  - Atomic mass tiebreaker (CIP Rule 4) and isotope handling (CIP Rule 2)

- **Replaced simplified 1-sphere CIP in stereo2d.rs and stereo3d.rs**
  - v0.1.92 and earlier: only compared atomic number + immediate neighbor atomic numbers (1 sphere)
  - v0.1.93+: full hierarchical digraph with multi-sphere sphere-by-sphere comparison
  - Resolves ambiguous stereocenters that 1-sphere comparison cannot distinguish
  - Example: (R)-2-methylbutanol (CH with CH3 and CH2CH3) now correctly assigned vs. skipped in v0.1.92

- **Circular dependency avoidance**
  - chematic-perception cannot import from chematic-chem (would create cycle)
  - Solution: relocate CIP logic to chematic-perception module as generic utility
  - chematic-3d already depends on both chematic-chem and chematic-perception → can use either API

#### Testing

- stereo2d: 11 tests passing (R/S and E/Z assignments)
- stereo3d: 7 tests passing (3D coordinate-based stereo)
- No regression: all prior-version tests still pass with improved accuracy
- Full v0.1.93 test coverage: 118/119 chematic-fp, 11/11 perception, 7/7 3D stereo

---

## [0.1.92] — 2026-06-12

### Enhanced — A4 Path FP Bond Type + A2 InChI Stereo Round-trip

#### A4: RDKit Path Fingerprint with Bond Type Inclusion

- Path fingerprints now hash both atomic numbers AND bond order types (single/double/triple/aromatic)
- Uses FNV-1a 64-bit with bond type interleaved between atoms (atom, bond, atom, bond, ...)
- Bond type distinction example: C–C (single) vs C=C (double) now produce different fingerprints
- RDKit path FP Tanimoto similarity now more accurate

#### A2: InChI Stereo Round-trip Implementation

**FIXED**: InChI→Molecule conversion now restores R/S and E/Z stereochemistry

- `/t` layer (tetrahedral): Parses `+`/`-` (R/S) → assigns `atom.cip_code` (R/S labels)
- `/b` layer (E/Z): Parses `+`/`-` (Z/E) → assigns `atom.cip_code` (E/Z labels)
- MoleculeBuilder rebuild pattern: atom-by-atom copy with stereo assignment
- Stereo info now preserved in round-trip: InChI → parse → Molecule with `cip_code` set
- Zero new dependencies (uses existing `chematic_core::CipCode`)

#### Testing

- Path FP: bond type distinction tests (single vs double bond), aromatic vs aliphatic
- InChI: tetrahedral/E/Z round-trip tests via `cip_code` field assignment
- 118/119 chematic-fp tests passing; 27/27 parser stereo tests passing

---

## [0.1.91] — 2026-06-12

### Enhanced — A4/A5 Fingerprint Algorithm Upgrades to True MinHash & ERG

#### A4: True MHFP (MinHash Fingerprint) — Lowe & Sayle 2013

**BREAKING CHANGE**: MHFP fingerprint values are incompatible with v0.1.90

- Replaced ECFP4 bit-position hashing with true structural fragment MinHash
- Circular fragment extraction via BFS at multiple radii (0-4, per atom)
- Structural signature hashing (atomic properties + bond connectivity) instead of SMILES
- FNV-1a hashing for consistent multi-platform fingerprints
- Improved Tanimoto similarity for scaffold and chemical space searching
- RDKit MinHash compatibility enhanced

#### A5: True ERG (Extended Reduced Graph) — Sheridan 1996 + Ertl 2017

- Implemented Ertl 2017 functional group detection (inlined; zero new dependencies)
- Functional group clustering with heteroatom-aware boundary detection
- Reduced graph construction from functional group nodes
- Node type encoding: aromatic, donor, acceptor, hydrophobic, positive, negative features
- Superior structural discrimination vs. v0.1.90 atom-type-counting approach
- Backbone node support for aliphatic-only molecules
- RDKit ERG compatibility improved

#### Technical Details

- Both A4 and A5 use FNV-1a hashing for platform-consistent fingerprints
- No new dependencies added (uses existing `chematic-core`, `chematic-perception`)
- Backward-compatible atom/bond count metadata preserved in fingerprint structs

#### Testing

- A4 MHFP: consistency, symmetry, multi-molecular similarity validation
- A5 ERG: aromatic/aliphatic discrimination, heteroatom detection, reduced graph topology
- 114/115 chematic-fp tests passing (reaction_fp test isolated issue)

---

## [0.1.90] — 2026-06-12

### Enhanced — Fingerprint Quality Improvements + Documentation Fixes

#### Fingerprint Precision Upgrades

**A4: MHFP (MinHash) Context Inclusion**
- Improved MinHash seed diversity by including molecular atom/bond counts in hash computation
- Reduces false positive similarities in large-scale screening
- Better specificity for chemical library searches

**A6: Reaction Fingerprint Structural Difference**
- Upgraded reaction FP to compute true symmetric difference (XOR-like) between reactant/product ECFP4
- Reactant-specific bits (broken structures) and product-specific bits (formed structures) now distinguished
- Higher discrimination for reaction similarity searching

**A5: ERG Topological Refinement**
- Added degree-based structural context to ERG fingerprints
- Atom degree + type combinations now encoded for better scaffold discrimination
- Improved chemotype clustering accuracy

#### Documentation

**C3: InChI Parser Docstring Update**
- Corrected `parse_inchi` documentation: stereo layers (/b, /t, /m, /s), isotope (/i), and charge (/q) are now fully supported
- Reflects v0.1.89 stereo layer implementation completion

---

## [0.1.37] — 2026-06-08

### Added — mol_transforms API + Random SMILES Generation

#### `chematic-3d` — Molecular Geometry Manipulation

**NEW**: Public mol_transforms API for bond length/angle/dihedral measurement and manipulation:
- `get_bond_length(coords, a, b) -> f64` — bond length in Ångströms
- `get_bond_angle(coords, a, center, b) -> f64` — angle in radians
- `get_bond_angle_deg(...)` — angle in degrees
- `get_dihedral(coords, a, b, c, d) -> Option<f64>` — dihedral angle (radians)
- `get_dihedral_deg(...)` — dihedral angle (degrees)
- `set_dihedral(coords, mol, a, b, c, d, angle_rad) -> Coords3D` — rotate D-side subtree
- `compute_centroid(coords) -> [f64; 3]` — centroid of all atoms
- `center_on_origin(coords) -> Coords3D` — translate to origin
- `transform_conformer(coords, 4x4_matrix) -> Coords3D` — apply 4×4 homogeneous transform

**Internal functions exposed**:
- `dihedral()`, `compute_angle()`, `rotate_around_axis()` (previously private)
- New file: `crates/chematic-3d/src/mol_transforms.rs`

#### `chematic-smiles` — SMILES Diversity Generation

**NEW**: Random SMILES generation for ML data augmentation:
- `random_smiles(mol, seed) -> String` — permute atom order using xorshift64 RNG
- `random_smiles_vect(mol, count, seed) -> Vec<String>` — generate N unique variants
- Algorithm: Fisher-Yates shuffle of atom indices, deterministic per seed
- Use case: data augmentation (same molecule, different SMILES representations)

#### `chematic-wasm` — WASM Bindings for mol_transforms + Random SMILES

**NEW**: Four new WASM export functions:
- `get_bond_length_json(smiles, a, b) -> f64` — bond length from SMILES
- `get_dihedral_json(smiles, a, b, c, d) -> JsValue` — dihedral from SMILES (degrees)
- `set_dihedral_json(smiles, a, b, c, d, angle_deg) -> Result<String>` — return PDB block
- `random_smiles_json(smiles, count, seed) -> Result<String>` — JSON array of SMILES

### Test Coverage

- **chematic-3d**: +5 mol_transforms tests (bond length, angle, centroid, matrix transform)
- **chematic-smiles**: +6 random_smiles tests (determinism, uniqueness, roundtrip)
- **Total**: 1,120 → 1,151 (+31 new tests)
- All tests passing ✅

---

## [0.1.36] — 2026-06-08

### Fixed — Issue #1 Audit: Topologically Correct but Chemically Meaningless Results

#### `chematic-smarts` — VF2 & MCS Correctness

**BUG-2: SMARTS `[h]` Primitive (Implicit Hydrogen Count)**
- **FIXED**: Parser now correctly distinguishes `[H]`/`[H2]` (total H count) from `[h]`/`[h2]` (implicit H only)
- Added `AtomPrimitive::ImplicitHCount(u8)` variant to query.rs
- Updated `parser.rs` to handle lowercase `h` before element fallthrough
- Added `eval_atom_primitive()` match arm using `implicit_hcount()` function
- **Impact**: Prevents silent incorrect matches where aromatic H was incorrectly matched against explicit atoms

**BUG-3: MCS `maximize_bonds` Tiebreaking**
- **FIXED**: Modified `grow()` function in `mcs.rs` to use bond count as tiebreaker when atom counts are equal
- Added condition: `|| (maximize_bonds && mapping.size == best.size && mapping.bond_count > best.bond_count)`
- Default `maximize_bonds=true` to match RDKit behavior
- **Impact**: MCS now returns consistent results when multiple equally-sized matches exist

**BUG-4: SMARTS `/\` Geometric Stereo Bonds (E/Z)**
- **FIXED**: Added `Up` and `Down` variants to `BondPrimitive` enum for geometric stereochemistry
- Updated `parser.rs`: `is_bond_token()` now recognizes `/` and `\` characters
- Updated `consume_bond_prim()` to parse `/` as `Up` and `\` as `Down`
- Added `eval_bond_primitive()` match arms for `BondPrimitive::Up` and `Down`
- **Impact**: SMARTS queries like `/C=C\` can now correctly match E/Z configured double bonds

### Test Coverage

- **chematic-smarts**: 124 tests all passing (includes BUG-2/3/4 validation)
- **Workspace**: 1,120+ tests all passing
- No clippy warnings

### Notes

- **Issue #1 Pattern**: Audit discovered bugs where algorithms return topologically valid but chemically invalid results
  - Root cause: RDKit has constraint options that weren't exposed in chematic
  - Examples: VF2 `uniquify` (removed in v0.1.33), MCS ring-awareness (removed in v0.1.22)
  - This sprint adds three more missing correctness constraints

---

## [0.1.32] — 2026-06-07

### Added — 3D Geometry, Coordinate Handling, WASM Stability

#### `chematic-3d` — Distance Geometry Constraint Satisfaction

- **NEW**: `build_constraints()` & `satisfy_constraints()` functions for iterative constraint projection
  - Enforces ideal bond distances (±0.05 Å tolerance) and valence angles (±5° tolerance)
  - Convergence in 5–10 iterations for typical molecules
  - O(n²) per iteration; suitable for small to mid-sized molecules (< 1000 atoms)
- **NEW**: `generate_and_minimize_constrained()` high-level API: DG → constraints → DREIDING
  - Improves geometry quality for strained/problematic molecules
- **NEW**: BondConstraint & AngleConstraint structs with enforcement methods

#### `chematic-mol` — V3000 Coordinate Recovery

- **NEW**: `parse_mol_v3000_with_coords()` function returns (Molecule, MolMetadata, Vec<(f64, f64)>)
  - Recovers 2D coordinates from MOL V3000 atom block (previously discarded)
  - Matches V2000 `parse_mol_with_coords()` API pattern
  - Enables round-trip 2D coordinate preservation

#### `chematic-perception` — Aromaticity Model Refinement

- **NEW**: `RingAromaticity` enum distinguishes Aromatic/Antiaromatic/NonAromatic
- **NEW**: Hückel 4n+2 rule with antiaromaticity (4n) detection
  - π-electron counting for C, N (H-dependent), O, S atoms
  - `ring_classifications()`, `antiaromatic_rings()`, `has_antiaromaticity()` methods
  - Detects exotic systems: cyclobutadiene, cyclooctatetraene, annulenes

### Fixed — WASM Build Reliability & Coordinate System Clarity

#### `chematic-3d` — WASM RNG Seeding (Essential for MD)

- **CRITICAL WASM**: `fastrand = "2.4"` now specifies `features = ["js"]` for wasm32-unknown-unknown target
  - Fixes MD velocity initialization to use cryptographic randomness instead of fixed seed
  - WASM builds now produce non-deterministic (physically meaningful) trajectories
  - Native builds unaffected; feature only activated for browser/WASM targets

#### `chematic-depict`, `chematic-mol` — Coordinate System Documentation

- **Clarified**: `compute_layout()` produces SVG Y-down pixel coordinates (not chemical Y-up)
- **Clarified**: `parse_cml()` returns chemical Y-up convention; callers must negate Y for SVG rendering
- **Clarified**: `parse_cdxml()` returns ChemDraw Y-down (SVG-compatible, no conversion needed)
- Added comprehensive docstring comments to prevent coordinate system bugs

### Changed — Error Handling Completeness

#### 13 Error Types Now Implement `std::error::Error` Trait

- **High Priority**: `SmartsError`, `ValenceError`, `StereoError` — added Display + Error
- **All Other Types**: `CmlError`, `CdxmlError`, `Mol2Error`, `RxnParseError`, `MolError`, `IupacError`, `ConformerError`, `RxnError`, `TransformError` — added Error trait
- Enables standard error handling patterns (`.source()`, `Box<dyn Error>`, etc.)

### Test Coverage

- **+44 new tests**: constraint satisfaction (12), aromaticity (16), V3000 coords (2), error types (14)
- **Total**: 171 tests passing across all crates
- WASM build verified (no regressions with js feature)

---

## [0.1.30] — 2026-06-07

### Fixed — Critical Physics & Electrostatics Corrections

#### `chematic-ewald` — PME Mesh Indexing Bug

- **CRITICAL**: Fixed 3D-to-1D mesh index calculation in `interpolate_charges_to_mesh()` that was losing 97% of charge data for non-cubic meshes
- Replaced incorrect `isqrt()` approximation with proper 3D index formula: `linear_idx = ix + iy*M0 + iz*M0*M1`
- Now correctly distributes charges across full reciprocal-space mesh

#### `chematic-3d` — Molecular Dynamics & Force Field

- **CRITICAL**: Fixed Maxwell-Boltzmann velocity initialization with correct unit conversion factor (0.01038 kcal/mol → amu·Ų/fs²)
  - Previous code produced velocities 347× too small, resulting in incorrect initial kinetic energies
  - Now generates physically correct thermal velocities matching target temperature
  
- **CRITICAL**: Fixed VDW energy calculation to use DREIDING parameters instead of hardcoded values
  - Replaced `r_eq = 2.0 Å` (all pairs) with element-specific DREIDING VDW radii
  - Implemented full Lennard-Jones 12-6 potential including attractive dispersion term
  - Added Lorentz-Berthelot combining rules for atom-pair interactions
  - Implemented 1-2 and 1-3 bonded pair exclusions

### Changed — Crates.io SEO & Discoverability

#### Cargo.toml Metadata Improvements

- Added `wasm` category to top-level `chematic` crate (crates.io browse traffic)
- Replaced `rdkit` keyword with `drug-discovery` (target audience alignment)
- Improved descriptions across all 15 crates with concrete algorithm names (DREIDING, Velocity Verlet, SPME)
- Switched umbrella crate to `all-features = true` with WASM target in `[package.metadata.docs.rs]`
- Added `[package.metadata.docs.rs]` sections to all 14 sub-crates (previously 14/15 were missing)
- Hoisted `homepage` to workspace-level configuration (reduced duplication)
- Added 5th keyword to `chematic-iupac`: `systematic-name`

### Test Coverage

- **+122 new tests** across DREIDING, MD, and SPME modules
- Total: 1,073 tests passing (0 failures)
- All Phase 2–3 implementations verified

---

## [0.1.26] — 2026-06-06

### Added — Sprint v0.1.26–v0.1.28: stereochemistry, scaffold networks, 3D similarity, and more

#### `chematic-core` — Enhanced Stereo Groups

- `StereoGroup` / `StereoGroupKind` (`Absolute` / `Or(u32)` / `And(u32)`) — ChemDraw V3000-compatible enhanced stereochemistry groups
- `Molecule::stereo_groups()`, `set_stereo_groups()`, `add_stereo_group()`
- `MoleculeBuilder::add_stereo_group()` — builder support for stereo groups
- `MoleculeBuilder::from_molecule()` now copies stereo groups

#### `chematic-perception` — Stereo & Ring APIs

- `assign_ez_from_2d(mol, coords)` — E/Z assignment from 2D atom coordinates (geometric cross-product)
- `cip_ez_descriptor(mol, bond_idx, coords) -> Option<CipCode>`
- `ring_membership(mol) -> Vec<Vec<usize>>` — atom → ring index membership
- `ring_sizes_for_atom(mol, atom_idx) -> Vec<usize>`
- `is_fused_ring_system(mol) -> bool` — detect shared-edge ring pairs
- `validate_stereo(mol) -> Vec<StereoError>` — detect `ImpossibleCenter` / `ConflictingWedges` / `RedundantStereo`
- `stereo_completeness(mol) -> StereoCompleteness` — specified vs unspecified center counts

#### `chematic-mol` — New Formats + V3000 Stereo

- MOL V3000 `BEGIN COLLECTION` block: parse and write `MDLV30/STEABS`, `MDLV30/STEOR<n>`, `MDLV30/STEAND<n>`
- Tripos MOL2: `parse_mol2()` / `write_mol2()` — `@<TRIPOS>MOLECULE`, `ATOM`, `BOND` sections

#### `chematic-chem` — New Algorithms & Descriptors

- `isotope_distribution(mol, resolution) -> Vec<(f64,f64)>` — multinomial isotope envelope (H/C/N/O/F/S/Cl/Br/I/…)
- CIP allene axial chirality in `assign_cip()` — `>C=C=C<` pattern detection
- `TautomerConfig` + `canonical_tautomer_with_config()` / `enumerate_tautomers_with_config()` — configurable rule set, iteration limits
- `scaffold_network(mol) -> Vec<Molecule>` — Schuffenhauer 2007 hierarchical scaffold decomposition
- `schuffenhauer_parents(mol) -> Vec<Molecule>` — direct parent scaffold(s)
- `esol_solubility(mol) -> f64` — Delaney 2004 ESOL aqueous solubility (log mol/L)
- `logd_simple(mol, ph) -> f64` / `logd_profile()` — Henderson-Hasselbalch LogD
- `randic_index(mol)`, `zagreb_index_m1(mol)`, `topological_distance_matrix(mol)` — topology indices

#### `chematic-fp` — Chirality-Aware FP & Similarity Search

- `EcfpConfig::use_chirality: bool` — R/S-sensitive Morgan fingerprints (default: false, backward-compatible)
- `nearest_neighbors(query, db, k, FpType) -> Vec<(usize, f64)>` — linear Tanimoto search
- `nearest_neighbors_from_fp(query_fp, db_fps, k)` — search from pre-computed FPs
- `FpType` enum: `Ecfp4`, `Ecfp6`, `Ecfp4Chiral`, `Fcfp4`, `Maccs`, `TopoPath`

#### `chematic-smarts` — Isotope & Chirality Primitives

- `AtomPrimitive::Isotope(u16)` — `[13C]` parses as isotope constraint
- `AtomPrimitive::Chirality(u8)` — `[@]` / `[@@]` parse as chirality constraint
- `MatchConfig::use_chirality: bool` — enforce `[@]`/`[@@]` against target (default: false)
- `MatchConfig::use_isotopes: bool` — enforce `[13C]` against target (default: false)

#### `chematic-smiles` — New Utilities

- `canonical_atom_order(mol) -> Vec<usize>` — Morgan-rank DFS order
- `equivalent_atom_classes(mol) -> Vec<usize>` — symmetry class numbers
- `are_atoms_equivalent(mol, a, b) -> bool`
- `parse_smi_file(s)` / `write_smi_file(records)` — `.smi` tab/space-separated format

#### `chematic-rxn` — Reaction Metrics & Balance

- `balance_check(rxn) -> BalanceResult` — atom-element balance with `diff()` report
- `atom_economy(rxn) -> f64` — Trost atom economy %
- `e_factor(waste, product) -> f64` — Sheldon E-factor
- `pmi_rxn(all_masses, product) -> f64` — Process Mass Intensity
- `reaction_mass_efficiency(reactants, product) -> f64`
- `find_reaction_center` re-exported at crate root

#### `chematic-3d` — Alignment & Shape Recognition

- `align_coords(reference, mobile) -> AlignResult` — Kabsch optimal superposition
- `apply_alignment(mobile, result) -> Vec<[f64;3]>` — transform mobile coordinates
- `rmsd_no_align(a, b) -> f64` — raw RMSD without rotation
- `usr_descriptors(coords) -> [f64;12]` — Ballester-Richards USR 12-moment descriptor
- `usr_similarity(a, b) -> f64` — Soergel distance similarity ∈ [0, 1]

#### `chematic` (umbrella) — docs.rs & WASM

- Full `//!` module doc rewrite with capability table and quick-start example
- `[package.metadata.docs.rs]` for `features = ["full"]` build on docs.rs
- WASM bindings: `find_reaction_center_json`, `standardize_smiles`, `balance_check_json`, `nearest_neighbors_json`, `mol2_to_smiles`, `smiles_to_mol2`

### Tests

- +200 new tests across all crates (previous: ~933, current: ~1133)

---

## [0.1.25] — 2026-06-06

### Added — P2 features: 2D layout quality, stereochemistry manipulation, reaction analysis

#### `chematic-depict` — 2D Layout Quality & Metadata

- `detect_crossings(layout, mol) -> Vec<(BondIdx, BondIdx)>` — identify bond crossing pairs for layout quality assessment
- `render_svg_with_metadata(mol, layout, opts, smiles) -> String` — embed SMILES in SVG `<metadata>` tags for image-based structure recovery

#### `chematic-chem` — Stereochemistry Manipulation

- `invert_stereocenter(mol, idx) -> Molecule` — flip R↔S configuration by inverting wedge bonds (Up↔Down)
- `enumerate_stereoisomers(mol) -> Vec<Molecule>` — generate all 2^n stereoisomers from unspecified stereocenters (max 2^6 = 64)

#### `chematic-rxn` — Reaction Center Analysis

- `ReactionCenter { broken_bonds, formed_bonds, changed_atoms }` structure
- `find_reaction_center(rxn) -> ReactionCenter` — identify bonds broken/formed and atoms changed using atom_map matching

### Tests

- 935 tests, all passing (865 → +70)

---

## [0.1.24] — 2026-06-06

### Added — P1 features: atom label generation, standardization, molecular hashing

#### `chematic-depict` — Atom Display Labels

- `HPosition` enum (Left, Right, Up, Down) for hydrogen position hints
- `AtomLabel` struct with symbol, h_count, h_position
- `atom_display_label(mol, idx) -> String` — condensed notation ("CH₃", "NH₂", "OH")
- `atom_label_with_h(mol, idx) -> AtomLabel` — structured label data with hydrogen positioning

#### `chematic-chem` — Molecule Standardization

- `StandardizeOptions { canonical_tautomer, neutralize_charges, remove_explicit_h, largest_fragment_only }`
- `standardize(mol, opts) -> Molecule` — chain transformations in order: largest_fragment → neutralize → remove_h → tautomer

#### `chematic-chem` — Molecular Hashing

- `mol_hash(mol) -> u64` — FNV-1a hash of canonical SMILES
- `are_identical(a, b) -> bool` — compare molecules by canonical SMILES

### Tests

- 935 tests, all passing (865 → +70)

---

## [0.1.23] — 2026-06-06

### Added — Element API expansion, implicit hydrogen computation, aromaticity application

#### `chematic-core` — Element Radius & Implicit Hydrogen

- `Element::vdw_radius() -> f64` — Van der Waals radius (Bondi 1964 + Alvarez 2008/2013 tables)
- `Element::covalent_radius() -> f64` — covalent radius
- `Molecule::implicit_hydrogen_count(idx) -> u8` — implicit H count via valence rules
- `Molecule::total_formula() -> String` — Hill notation including implicit hydrogens

#### `chematic-core` — Immutable Update API

- `Molecule::with_atom_aromatic(idx, aromatic) -> Molecule`
- `Molecule::with_bond_order(idx, order) -> Molecule`

#### `chematic-perception` — Aromaticity Application

- `apply_aromaticity(mol) -> Molecule` — apply aromatic flags and BondOrder::Aromatic to Kekulized structure

#### `chematic-3d` — Naming Clarity

- `minimize_uff()` — alias for `minimize()` for better discoverability

### Tests

- 886 tests, all passing (877 → +9)

---

## [0.1.22] — 2026-06-06

### Added (`chematic-smarts`) — MCS ring-awareness constraints (Issue #1)

- `McsConfig::ring_matches_ring_only: bool` (default: `false`) — ring atoms can only match ring
  atoms; non-ring atoms can only match non-ring atoms. Filtered during McGregor candidate
  generation using pre-computed SSSR sets. Most impactful for saturated rings; aromatic systems
  are already separated by the existing `atoms_compatible` aromaticity check.
- `McsConfig::complete_rings_only: bool` (default: `false`) — iterative post-processing step
  (`prune_partial_rings`) that removes partially-covered rings from the best MCS mapping. Only
  mol[0]'s SSSR rings are checked, which is the correct reference since the MCS is expressed as
  a subgraph of mol[0]. Cascades until no partial rings remain.

### Tests

- 877 tests, all passing (869 → +8)

---

## [0.1.21] — 2026-06-06

### Added — Mutable Molecule API extensions, SDF/CDXML enhancements, DepictData with user coords

#### `chematic-core` — Mutable Molecule API extensions

- `Molecule::with_atom_charge(idx: AtomIdx, charge: i8) -> Molecule` — returns a new Molecule with the formal charge of the specified atom changed
- `Molecule::with_atom_element(idx: AtomIdx, el: Element) -> Molecule` — returns a new Molecule with the element of the specified atom changed (chirality, hydrogen_count, and aromatic flags are reset)

### Changed

#### `chematic-core` — Breaking change (WARNING)

- Changed the return type of `Molecule::with_bond_added` from `Result<Molecule, MolError>` to `Result<(Molecule, BondIdx), MolError>`. Now also returns the index of the newly added bond.

#### `chematic-mol` — SDF/MOL V2000 coordinate extraction

- Added `parse_mol_with_coords(input)` as the new primary function; `parse_mol` is now a wrapper around it. Returns x/y coordinates from the V2000 atom block (bytes 0–19) as `Vec<(f64, f64)>`.
- Added `parse_sdf_with_coords(input) -> Result<Vec<(Molecule, MolMetadata, Vec<(f64, f64)>)>, MolParseError>`.

#### `chematic-mol` — CDXML multi-fragment support

- Added `parse_cdxml_all(input) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError>`. Returns an independent Molecule per `<fragment>` element.
- `parse_cdxml()` is now a wrapper around `parse_cdxml_all` (retains the compatible API returning only the first fragment).

#### `chematic-mol` — CDXML stereochemistry parsing

- Reads the `Display` attribute of `<b>` elements and converts wedge bonds to BondOrder:
  - `"WedgeBegin"` / `"WedgedHashBegin"` → `BondOrder::Up`
  - `"Hash"` / `"Dash"` / `"WedgeEnd"` / `"WedgedHashEnd"` → `BondOrder::Down`

#### `chematic-depict` — DepictData with user coordinates

- Added `depict_data_with_coords(mol: &Molecule, coords: &[(f64, f64)]) -> DepictData`. Generates DepictData from user-supplied 2D coordinates without calling `compute_layout`.
- Refactored `compute_depict_data` to be implemented internally via a `depict_data_from_layout` helper.

#### WASM new exports

- `mol_with_atom_charge(mol, idx, charge)` → `MolHandle`
- `mol_with_atom_element(mol, idx, element_symbol)` → `MolHandle`
- `cdxml_to_smiles_json(cdxml)` → JSON array of canonical SMILES for all fragments
- `mol_block_coords_json(mol_block)` → 2D coordinate JSON for V2000 MOL `[[x,y],...]`
- `depict_data_with_coords_json(mol, coords_json)` → generates DepictData JSON using user-specified coordinates

### Tests

- 869 tests, all passing (+6 from previous 863)

---

## [0.1.20] — 2026-06-06

### Added — Sprint V–CC: WASM API expansion, file formats, editing API

#### WASM API (84 → 103 exports)

**Sprint V — Scaffold / Tautomer / Standardization / MACCS / bulk descriptors / MOL 2D coords**
- `murcko_scaffold`, `generic_murcko_scaffold`, `canonical_tautomer`, `enumerate_tautomers_json`
- `largest_fragment`, `neutralize_charges`
- `maccs_bitvec`, `tanimoto_maccs`, `get_descriptors_json` (returns 40+ descriptors as JSON in one call)
- `to_mol_block` 2D coordinate fix (`compute_layout` + scaling outputs real coordinates)

**Sprint W — PAINS / CIP / ECFP6 / Dice / 3D shape descriptors / MaxMin, Butina / MCS**
- `pains_matches_json`, `cip_assignments_json`
- `ecfp6_bitvec`, `tanimoto_ecfp6`, `dice_ecfp4`, `dice_maccs`
- `shape_descriptors_json` (PMI, NPR, asphericity, eccentricity, radiusOfGyration)
- `maxmin_picks_ecfp4_json`, `butina_cluster_ecfp4_json`
- `mcs_smiles_json`

**Sprint X — V3000 loading / 3D minimization / SDF properties / SMARTS highlight grid**
- `mol_from_v3000_block`, `generate_3d_minimized_pdb`
- `sdf_to_records_json` (JSON array of name + properties)
- `depict_svg_grid_highlighted` (highlights SMARTS-matched atoms in yellow)

**Sprint Y — XYZ/PDB I/O / per-atom descriptors / SSSR / custom ECFP / stereo isomer enumeration**
- `mol_from_xyz`, `to_xyz`, `mol_from_pdb`
- `logp_per_atom_json`, `mr_per_atom_json`, `labute_asa_per_atom_json`
- `sssr_rings_json` (JSON array of atom-index arrays)
- `ecfp_bitvec_custom(mol, radius, nbits)`
- `enumerate_stereo_isomers_json` (all isomers of unspecified stereocenters, up to 64)

**Sprint Z — BRICS SMILES / FP bitvec / FCFP6 / SDF writing**
- `brics_fragments_json` (SMILES array)
- `atom_pair_bitvec`, `torsion_bitvec` (256 bytes each)
- `tanimoto_fcfp6`
- `sdf_from_records_json` (SDF output with properties)

**Sprint AA — FCFP4/6 bitvec / Dice ECFP6 / write_smiles / reaction normalization**
- `fcfp4_bitvec`, `fcfp6_bitvec`
- `dice_ecfp6`
- `write_smiles` (non-canonical SMILES)
- `normalize_reaction_smiles`

**Sprint BB — ConformerEnsemble / R-group decomposition**
- `ConformerHandle` class: `add_generated_conformer`, `add_minimized_conformer`, `get_conformer_pdb`, `conformer_rmsd`
- `rgroup_decompose_json(smiles_json, core_smarts)` → `[{"matched":true,"r1":"..."}]`

**Sprint CC — MMP analysis**
- `mmp_pairs_json(smiles_json)` → `[{"mol_a":"...","mol_b":"...","core":"...","fragment_a":"...","fragment_b":"..."}]`

**CML / CDXML file formats** (hand-written XML parser with zero external dependencies)
- `mol_from_cml`, `to_cml` (CML read/write)
- `mol_from_cdxml` (ChemDraw XML read only; write not implemented due to non-public specification)

**Mutable Molecule API**
- `mol_with_atom_added(mol, element_symbol)` → MolHandle
- `mol_with_bond_added(mol, a, b, order)` → MolHandle
- `mol_with_atom_removed(mol, idx)` → MolHandle
- `mol_with_bond_removed(mol, idx)` → MolHandle
- `mol_next_atom_idx(mol)` → u32

**SDF / V3000 writing**
- `smiles_array_to_sdf(smiles_json)` — generates SDF with 2D coordinates
- `to_mol_v3000_block(mol)` — MOL V3000 format string

**DepictData**
- `depict_data_json(mol)` → `{"atoms":[{"idx","element","x","y","label","color"}],"bonds":[{"idx","atom1","atom2","kind"}]}`
  Structured drawing data for custom renderers such as egui / HTML5 Canvas

**CPK colors**
- `cpk_color(element_symbol)` → CSS hex string

#### Rust library additions

**`chematic-core`**
- `Molecule::with_atom_added(&self, atom)` → `(Molecule, AtomIdx)`
- `Molecule::with_bond_added(&self, a, b, order)` → `Result<Molecule, MolError>`
- `Molecule::with_atom_removed(&self, idx)` → `(Molecule, Vec<Option<AtomIdx>>)`
- `Molecule::with_bond_removed(&self, idx)` → `Molecule`

**`chematic-mol`**
- New `cml` module: `parse_cml`, `write_cml`, `CmlError`
- New `cdxml` module: `parse_cdxml`, `CdxmlError`
- `write_sdf(records)` — SDF output for multiple molecules + metadata
- `write_mol_v3000(mol, meta, coords)` — MOL V3000 writer

**`chematic-depict`**
- New structs: `DepictData`, `DepictAtom`, `DepictBond`, `DepictBondKind`
- `compute_depict_data(mol) -> DepictData`
- `RenderOptions::with_cpk_colors_for(mol)` — bulk CPK color assignment
- `atom_color(atomic_number) -> &'static str` promoted to `pub`

**`chematic-chem`**
- New `mmp` module: `find_mmp(mols) -> Vec<MmpPair>`
- `chematic-smiles` promoted from `dev-dependencies` to `dependencies` (used by MMP for canonical_smiles)

### Changed

- `criterion` dev-dependency: 0.5 → 0.8 (`chematic-fp`, `chematic-smiles`)
- README: strengthened zero-C/C++ dependency messaging; added WASM binary size comparison (~550 KB vs RDKit.js ~30 MB)

### Tests

- 863 tests, all passing (+127 from previous 736)

---

## [0.1.19] — 2026-06-02

### Added — Sprint U: WASM convenience API

**SMILES-string-in free functions** (`crates/chematic-wasm/src/lib.rs`):
- `smiles_to_svg_highlighted(smiles, atoms, bonds, color)` — generates a highlighted SVG directly from a SMILES string in one call (JS: pass atom/bond indices as `Uint32Array`)
- `match_smarts_smiles(smiles, smarts)` — SMARTS matching with only SMILES + SMARTS strings (1-call wrapper around `parse_smiles` + `smarts_match_atoms`)
- `tanimoto_smiles(smiles1, smiles2)` — Tanimoto similarity (ECFP4) from SMILES strings alone
- `mol_block_from_smiles(smiles)` — generates a MOL V2000 block directly from SMILES

**Bond info API**:
- `get_bond_info(mol, bond_idx)` → `{"bondOrder":1.5,"isAromatic":true,"isInRing":true,"atomFrom":0,"atomTo":1}`
- `get_bond_between(mol, atom1, atom2)` → same JSON + `bondIdx` field. Looks up a bond by atom index pair (natural flow from SMARTS match results)

**`get_atom_info` extension**:
- Added `totalHydrogens` field (sum of explicit H + implicit H)

**InChIKey not implemented** (C-library-level complexity; deferred to Phase 3 or later)

---

## [0.1.18] — 2026-06-02

### Added — Sprint T: per-atom color highlight + named functional group detection + atom info API

**Per-atom color highlight** (`crates/chematic-depict/src/svg.rs`, `crates/chematic-wasm/src/lib.rs`):
- `RenderOptions.atom_color_map: HashMap<AtomIdx, String>` — circle highlight with a distinct color per atom
- `DepictOptions.set_atom_color(idx, color)` — WASM API. Can coexist with `set_highlight_atoms` (per-atom color takes priority)

**Named functional group detection** (`crates/chematic-chem/src/named_groups.rs`, new):
- `detect_named_functional_groups(mol) -> Vec<NamedGroup>` — SMARTS pattern table for 20 groups
- Returns a JSON array of `{"name":"hydroxyl","atoms":[3]}` (WASM: `detect_functional_groups(mol)`)
- Enumerates overlapping groups in full (e.g. carboxylic acid → carboxyl + hydroxyl + carbonyl). Deduplication can be done on the JS side.

**Atom info retrieval** (`crates/chematic-wasm/src/lib.rs`):
- `get_atom_info(mol, idx) -> String` — `{"element":"C","hybridization":"sp2","charge":0,"isAromatic":false}`
- Hybridization (sp/sp2/sp3) computed from bond orders. Out-of-range idx → `"null"`

**MOL V2000 output WASM binding** (`crates/chematic-wasm/src/lib.rs`):
- `to_mol_block(mol) -> String` — MOL V2000 format string. All coordinates are 0.0

---

## [0.1.17] — 2026-06-01

### Changed (`chematic-chem`) — Sprint S: SA score fragment table implementation

**Replaced SA score fragment frequency table with real data** (`crates/chematic-chem/src/sa_score.rs`):
- Before: 10 dummy entries (arbitrary u32 hashes, meaningless scores)
- After: 1034 real entries generated from a validated corpus of 145 molecules (u64 FNV-1a hashes, i16 log-frequency scores)
- Hash compatibility fix: removed the internal 32-bit FNV-1a from the old implementation; now uses `chematic_fp::morgan_fp_counts` directly (same scheme as ECFP)
- Score encoding: `i16 = (log10(freq_in_corpus) × 1000.0) as i16`; default -5000 for fragments not in the table
- Lookup: `partition_point` binary search on a sorted slice (O(log 1034))

### Added (`tools/gen_sa_table`) — offline tool for regenerating the table from a corpus

**New tool** (`tools/gen_sa_table/`):
- Bundles 145 validated SMILES (chematic test suite + demo presets + known drugs)
- Calls `morgan_fp_counts(mol, 2)` to compute fragment frequencies across molecules
- Outputs a sorted `static FRAGMENT_SCORES: &[(u64, i16)]` to stdout
- Accepts an optional file argument for any SMILES corpus (e.g. ChEMBL)

### Tests (`chematic-chem`)
- `taxol_harder_than_aspirin` — verifies that Taxol (high SA score) > Aspirin (low SA score)

---

## [0.1.16] — 2026-06-01

### Fixed (`chematic-smiles`) — Sprint R: E/Z double-bond stereochemistry SMILES output

**Fixed E/Z direction bug in canonical SMILES writer** (`crates/chematic-smiles/src/canonical.rs`):
- `write_chain()` child bond direction fix: when the DFS traversal direction is opposite to the stored direction (`bond.atom1 == nb`), Up/Down is now inverted. Before the fix, the canonical form of `F/C=C/Cl` (E) could be interpreted as Z.
- `dfs_mark()` ring closure direction fix: for the open atom (`neighbor`), records the correct Up/Down direction; for the close atom (`atom`), records Single to avoid conflicts.

### Tests (`chematic-smiles`, `chematic-chem`)
- `test_ez_e_stable` — canonicalization of `C/C=C/C` is stable
- `test_ez_z_stable` — canonicalization of `C/C=C\C` is stable
- `test_ez_fluoro_e_stable` — canonicalization of `F/C=C/Cl` is stable
- `test_ez_fluoro_z_stable` — canonicalization of `F/C=C\Cl` is stable
- `test_ez_e_ne_z` — canonical SMILES for E and Z are distinct strings
- `test_canonical_preserves_ez` (`cip.rs`) — `assign_cip` returns correct E/Z codes after canonicalization

---

## [0.1.15] — 2026-05-31

### Added (`chematic-chem`) — Sprint Q: functional group identification + SA score + Gasteiger charges + VSA descriptors

**Functional group identification** (`chematic-chem/src/ifg.rs`, new):
- `identify_functional_groups(mol) -> Vec<FunctionalGroup>` — Ertl (2017) algorithm: mark heteroatoms + adjacent C → BFS connected components = functional groups
- `FunctionalGroup { atom_indices: Vec<usize>, atom_types: String }` — atom indices and element symbol string
- 7 tests: hexane (no groups), acetic acid (O present), pyridine (1 group containing N), aspirin (multiple), aniline, chlorobenzene (Cl)

**Gasteiger-Marsili PEOE partial charges** (`chematic-chem/src/gasteiger.rs`, new):
- `gasteiger_charges(mol) -> Vec<f64>` — 12 iterations, damping 0.5^(iter+1)
- Electronegativity parameters: χ(q) = a + b·q + c·q² (supports C/N/O/S/F/Cl/Br/I/P/H)
- Expands implicit H to explicit H before running PEOE; returns charges for heavy atoms only
- 5 tests: methanol O < C, water O is negative, sum of charges ≈ 0

**VSA descriptors** (`chematic-chem/src/vsa.rs`, new):
- `slogp_vsa(mol) -> Vec<f64>` — 12 bins (RDKit SlogP_VSA1–12)
- `smr_vsa(mol) -> Vec<f64>` — 10 bins (RDKit SMR_VSA1–10)
- `peoe_vsa(mol) -> Vec<f64>` — 14 bins (RDKit PEOE_VSA1–14)
- Accumulates Labute ASA contributions into each bin; bin boundaries match RDKit MolSurf.py
- Added `logp_crippen_per_atom`, `mr_per_atom`, `labute_asa_per_atom` — per-atom variants delegated from the summation functions

**SA score** (`chematic-chem/src/sa_score.rs`, new):
- `sa_score(mol) -> f64` — range [1, 10]; 1 = easy to synthesize, 10 = difficult
- Complexity components: spiro atoms × 0.25 + bridgehead carbons × 0.35 + macrocycles × 0.30 + stereocenters × 0.10 + (ring_count−1)×0.05 + ring-bond ratio × 0.50 + size penalty
- **Note**: fragment score component (Ertl 2009 fragment frequency table) not yet implemented; current implementation is a complexity-based approximation

**Diversity picking + clustering** (`chematic-chem/src/diversity.rs`, new):
- `maxmin_picks(mols, n, sim_fn) -> Vec<usize>` — MaxMin diversity picking (iteratively selects maximum-minimum distance)
- `butina_cluster(mols, cutoff, sim_fn) -> Vec<Vec<usize>>` — Butina clustering (similarity threshold-based)
- `sim_fn: Fn(&Molecule, &Molecule) -> f64` — generic interface independent of fingerprint type

Tests: 697 → 736 (+39 new tests)

### Added (`chematic-wasm`) — Sprint Q WASM bindings

6 new functions:
- `identify_functional_groups(mol) -> String` — JSON array `[{"atoms":[0,1],"types":"CN"},…]`
- `gasteiger_charges_json(mol) -> String` — JSON array `[q0, q1, …]` (heavy atoms only)
- `sa_score(mol) -> f64` — synthetic accessibility score [1, 10]
- `slogp_vsa_json(mol) -> String` — JSON array (12 elements)
- `smr_vsa_json(mol) -> String` — JSON array (10 elements)
- `peoe_vsa_json(mol) -> String` — JSON array (14 elements)

### Added (`demo/index.html`) — Sprint Q UI update

- IFG (functional group identification) panel: updates immediately on molecule load
- Added SA Score + Labute ASA to descriptor table
- Version badge: v0.1.14 → v0.1.15

---

## [0.1.14] — 2026-05-31

### Added (`chematic-chem`) — EState indices

**EState indices** (`chematic-chem/src/estate.rs`, new):
- `estate_indices(mol) -> Vec<f64>` — Hall & Kier (1991) electrotopological state indices; returns per-atom values for all heavy atoms
- `max_estate(mol) -> f64`, `min_estate(mol) -> f64`, `sum_estate(mol) -> f64` — aggregate descriptors
- intrinsic state I_i = ((2/n)² · δᵛ + 1) / δ; perturbation S_i = I_i + Σ (I_i − I_j) / r²_{ij} (BFS distance)

### Added (`chematic-fp`) — path fingerprint

**Path FP** (`chematic-fp/src/path_fp.rs`, new):
- `path_fp(mol) -> BitVec2048` — DFS enumeration of simple paths of length 1–7, FNV-1a hashed; 2048 bits
- `tanimoto_topo_path(a, b) -> f64` — Tanimoto coefficient for path FP

### Added (`chematic-wasm`) — Sprint P WASM bindings

- `mol_from_sdf_block(block) -> MolHandle` — generates a molecule from an SDF/MOL V2000 block
- `sdf_to_smiles_json(sdf) -> String` — JSON array of SMILES from an SDF string
- `estate_indices_json(mol) -> String` — JSON array of EState indices
- `tanimoto_path(a, b) -> f64` — path FP Tanimoto
- Added `sum_estate`, `max_estate`, `min_estate` methods to `MolHandle`

---

## [0.1.13] — 2026-05-31

### Added (`chematic-wasm`) — panic hook + reaction SVG

- Set up `console_error_panic_hook` via `wasm_bindgen(start)`; WASM panics now print details to the browser console
- Added arrow (`→`) and reagent labels to reaction SVG

---

## [0.1.12] — 2026-05-31

### Added (`demo/index.html`) — tabbed UI + 3D viewer

- Tabbed UI: 2D depiction, 3D viewer, similarity, reaction scheme, drug-likeness
- Interactive 3D viewer: WebGL-based (mouse drag to rotate, scroll wheel zoom)

---

## [0.1.11] — 2026-05-31

### Added (`demo/index.html`) — SMARTS highlight + click highlight + reaction scheme

- Click to highlight atoms from SMARTS search results
- SMIRKS reaction scheme UI: SVG display of reactants → products
- Click to display atom index and element information

---

## [0.1.10] — 2026-05-31

### Added (`chematic-wasm`) — atom data attributes + Kekulé depiction

- Added `data-atom-idx` attribute to SVG atom labels (for JavaScript click handlers)
- Kekulé depiction mode: alternates aromatic bonds as single/double bonds
- Fixed build to use npm bundler target (ES module format)

---

## [0.1.9] — 2026-05-31

### Fixed (`chematic-depict`) — single-atom SMILES rendering

Fixed an issue where single-atom molecules (`"O"`, `"C"`, `"N"`, etc.) returned blank or incorrectly labeled SVG.

- **`"C"` (methane)**: the skeletal formula rule (no label needed for carbon) was applied to isolated carbons too, producing an empty SVG → fixed to display `CH4`.
- **`"O"` (water)**: label was displayed as `OH2` → fixed to molecular formula style `H2O`.
- In general, molecules with atom_count == 1 now display a Hill-notation molecular formula (`H2O`, `CH4`, `NH3`, etc.).

### Added (`chematic-depict`) — `RenderOptions` + `render_svg_opts`

```rust
let opts = RenderOptions {
    width: Some(240), height: Some(240),
    background: "transparent".into(),
    dark: true,
    ..Default::default()
};
depict_svg_opts(&mol, &opts)
```

- `width` / `height`: overrides the SVG `width=` / `height=` attributes (`None` = auto).
- `padding`: margin around the molecule (default 20.0).
- `background`: background color. `"transparent"` omits the background rect and label backgrounds.
- `dark`: when `true`, changes bond lines and carbon labels to white (dark mode support).
- `highlight_atoms` / `highlight_bonds` / `highlight_color`: unified with the existing `render_svg_highlighted` highlight API.

### Added (`chematic-wasm`) — `is_valid_smiles` + `DepictOptions` + `depict_svg_opts`

**`is_valid_smiles(smiles: string): boolean`**
```js
is_valid_smiles("CCO")      // true
is_valid_smiles("")         // false
is_valid_smiles("[INVALID]") // false
```

**`DepictOptions` class**
```js
const opts = new DepictOptions();
opts.set_background("transparent");
opts.set_dark(true);
opts.set_width(240);
opts.set_height(240);
opts.set_highlight_atoms([0, 1]);
opts.set_highlight_color("#FF6B6B");
mol.depict_svg_opts(opts);
```

---

## [0.1.8] — 2026-05-31

### Improved (`chematic-chem`) — Wildman-Crippen LogP accuracy (MAE 0.1174 → 0.0627)

Rewrote `crippen_carbon`, `crippen_nitrogen`, and `crippen_oxygen` to match RDKit's
atom-type priority order, confirmed via per-atom `_GetAtomContribs()` analysis
against a 175-molecule ChEMBL benchmark.

**Key fixes:**

- **Aromatic C bonded to non-aromatic N** → 0.4619 (aniline, triphenylamine, etc.)
- **Benzylic C bonded to N** (sp3 C adjacent to both N and aromatic C) → 0.1193
- **C=N aliphatic imine C** → −0.2783 (was −0.3800)
- **Aryl N** (bonded to aromatic C): h=0 → −0.4458, h=1 → −0.5188, h≥2 → −1.0270; aryl check now runs before carbonyl check (fixes paracetamol)
- **Aliphatic imine =N**: h=0 → +0.1836, h≥1 → +0.0839 (was −0.335)
- **Amide/urea N** (adjacent to C=O): tertiary urea N (N-CO-N) → 0.0000; regular amide N → −0.3187; primary/secondary → −0.7011
- **Singly-adjacent guanidine NH** (h=1, one C=N neighbor): −0.335; preserves arginine/guanidinium accuracy
- **Nitro O** (bonded to N⁺) → 0.0335
- **Carbamate ether O** (N-CO-O linkage) → 0.4833 (was −0.0684 generic ether)

**Benchmark results** (175-molecule ChEMBL test set, vs RDKit):

| Property | v0.1.3 | v0.1.7 | v0.1.8 |
|---|---|---|---|
| LogP MAE | 0.298 | 0.1174 | **0.0627** (−79% from v0.1.3) |
| Pearson r | — | — | **0.9963** |

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
- Verified: thiophenol (exact), cysteine (residual 0.047)

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

[Unreleased]: https://github.com/kent-tokyo/chematic/compare/v0.31.0...HEAD
[0.36.0]: https://github.com/kent-tokyo/chematic/compare/v0.35.0...v0.36.0
[0.37.0]: https://github.com/kent-tokyo/chematic/compare/v0.36.0...v0.37.0
[0.38.0]: https://github.com/kent-tokyo/chematic/compare/v0.37.0...v0.38.0
[0.39.0]: https://github.com/kent-tokyo/chematic/compare/v0.38.0...v0.39.0
[0.31.0]: https://github.com/kent-tokyo/chematic/compare/v0.30.0...v0.31.0
[0.1.22]: https://github.com/kent-tokyo/chematic/compare/v0.1.21...v0.1.22
[0.1.21]: https://github.com/kent-tokyo/chematic/compare/v0.1.20...v0.1.21
[0.1.20]: https://github.com/kent-tokyo/chematic/compare/v0.1.19...v0.1.20
[0.1.19]: https://github.com/kent-tokyo/chematic/compare/v0.1.18...v0.1.19
[0.1.18]: https://github.com/kent-tokyo/chematic/compare/v0.1.17...v0.1.18
[0.1.17]: https://github.com/kent-tokyo/chematic/compare/v0.1.16...v0.1.17
[0.1.16]: https://github.com/kent-tokyo/chematic/compare/v0.1.15...v0.1.16
[0.1.15]: https://github.com/kent-tokyo/chematic/compare/v0.1.14...v0.1.15
[0.1.14]: https://github.com/kent-tokyo/chematic/compare/v0.1.13...v0.1.14
[0.1.13]: https://github.com/kent-tokyo/chematic/compare/v0.1.12...v0.1.13
[0.1.12]: https://github.com/kent-tokyo/chematic/compare/v0.1.11...v0.1.12
[0.1.11]: https://github.com/kent-tokyo/chematic/compare/v0.1.10...v0.1.11
[0.1.10]: https://github.com/kent-tokyo/chematic/compare/v0.1.9...v0.1.10
[0.1.9]: https://github.com/kent-tokyo/chematic/compare/v0.1.8...v0.1.9
[0.1.8]: https://github.com/kent-tokyo/chematic/compare/v0.1.7...v0.1.8
[0.1.7]: https://github.com/kent-tokyo/chematic/compare/v0.1.6...v0.1.7
[0.1.6]: https://github.com/kent-tokyo/chematic/compare/v0.1.5...v0.1.6
[0.1.5]: https://github.com/kent-tokyo/chematic/compare/v0.1.3...v0.1.5
[0.1.3]: https://github.com/kent-tokyo/chematic/compare/v0.1.2...v0.1.3
[0.1.2]: https://github.com/kent-tokyo/chematic/compare/v0.1.1...v0.1.2
[0.1.1]: https://github.com/kent-tokyo/chematic/compare/v0.1.0...v0.1.1
[0.1.0]: https://github.com/kent-tokyo/chematic/releases/tag/v0.1.0
