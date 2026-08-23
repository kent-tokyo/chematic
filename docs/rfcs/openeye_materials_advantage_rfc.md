# RFC: OpenEye + Materials-Science Advantage

Status: draft, competitive audit + design round. Not merged, not marked
ready. **Zero changes under `crates/*/src/**` this round.**

## 0. Goal / non-goals

**Goal.** Grow chematic from "closing the gap with RDKit" (the existing
100-point ladder, `ROADMAP.md`) into a runtime that is:

1. OpenEye-class on 3D conformer generation, shape similarity, and
   microstate (protonation/tautomer) handling.
2. Capable of growing toward pymatgen/ASE-class materials-science
   coverage.
3. Semantically identical across Rust/Python/WASM.
4. Explainable, reproducible, and fail-closed (typed failures over silent
   wrong answers — the same discipline the 100-point ladder already
   applies to molecular work).

**Explicit non-goal**: copying OpenEye's full feature surface. OpenEye is
a mature, integrated, commercially-licensed pharma platform built over two
decades; chasing full parity is not realistic or useful. The competitive
bet this RFC lays out is narrower and sharper: build practical, real
competitiveness on **3D and shape** specifically (Track A), and win
decisively on **Rust/WASM support, explainability, and materials-science
integration** (Track B + chematic's existing strengths) — axes OpenEye
does not compete on at all.

**Hard constraints for this round and any future implementation it
recommends** (restated from the governing instructions, load-bearing for
every design choice below):

- No issues, PRs, questions, or contact of any kind with OpenEye, RDKit,
  pymatgen, ASE, spglib, or their communities. No external review
  requests.
- No reverse-engineering of OpenEye's licensed binaries; no transcription
  of proprietary algorithms or code. Only public documentation, public
  papers, public specifications, and public datasets were used to
  characterize what OpenEye/pymatgen/ASE/spglib do.
- No comparison against actual OpenEye execution output — this project
  has no OpenEye license. Every OpenEye-side claim in the gap matrix is a
  capability-shape description from public sources, not a measured
  number.
- Regular merges only (no squash/rebase/force-push); no merge without an
  explicit instruction; no version bump/tag/release/publish without a
  separate instruction; no routine full-corpus heavy remeasurement outside
  a release gate; the current zero-C/C++-dependency default build is never
  broken, and any future optional backend (e.g. a non-default spglib FFI
  bridge, if ever built) must never become the default.

## 1. Method

Three parallel code-reading passes (not assumption, not README text)
covered: (a) `chematic-3d`/`chematic-ff` for the Omega/OEFF-Szybki
comparison, (b) `chematic-chem`/`chematic-fp`/`chematic-depict`/
`chematic-mol` for the Quacpac/Shape/OEDepict/Spruce comparison, (c)
`chematic-crystal`/`chematic-mol`'s CIF-POSCAR-LAMMPS readers/
`chematic-ewald` for the materials comparison. Every claim in the two gap
matrices (`docs/competitive/openeye_gap_matrix.md`,
`docs/competitive/materials_gap_matrix.md`) traces to a specific file,
function/type name, test count, or measured percentage confirmed by
directly reading the source — not inferred from documentation or
inherited from project memory without re-verification. Three of the most
consequential claims (the 265-corpus's 241/265 figure, the `shape_
tanimoto()`-family gap already named in `ROADMAP.md`, and the confirmed
absence of any Gaussian-volume code) were independently spot-checked a
second time via direct `grep`/`CHANGELOG.md` reads before being written
into this RFC.

## 2. Two genuinely surprising findings from this audit

Worth surfacing before the gap matrices themselves, because both change
what "the gap" actually is:

1. **The public-facing "no RMSD/TFD comparison exists" claim
   (`docs/benchmark.md`, `docs/rdkit-migration.md`, and even today's
   v0.19.0 CHANGELOG) is false.** `validation/results/mmff94_bci_gap_227_
   phase2_report.md` already reports RMSD, TFD, stereo-satisfaction, and
   per-stage runtime as separate measured metrics on the 265-molecule
   corpus. This is a documentation gap, not a capability gap — chematic
   already has more 3D measurement rigor than it currently claims
   publicly. Fixing this citation is a small, separate follow-up, not part
   of Track A/B, but worth doing soon since it's an easy, free credibility
   win.
2. **The public, documented `Mol.conformer_ensemble()` API has a live,
   unfixed catastrophic-bond-tearing defect**, reachable by any caller,
   dormant since v0.14.0's CHANGELOG explicitly flagged "the default
   conformer path... is untouched." This is more concerning than an
   abstract "OpenEye ensemble quality is ahead" gap — it means chematic's
   *own* public ensemble API can silently hand back a structurally broken
   molecule today. See §3, A1.

## 3. Track A — OpenEye-class molecular modeling

Full gap grounding: `docs/competitive/openeye_gap_matrix.md`. For each
item: user value, gap vs. OpenEye, API sketch, failure taxonomy,
benchmark, fixture, holdout, CPU cost, WASM feasibility, breaking-change
risk, verdict.

### A1 — Conformer Benchmark & Failure Ledger (diagnosis only)

- **User value**: establishes ground truth for the single most-cited
  competitive weakness (`ROADMAP.md` line 48-55's own words) before
  spending engineering effort on it. Surfaces the two findings in §2 as
  formal, fixable-soon items.
- **Gap vs. OpenEye**: Omega reports success/failure with typed reasons
  and separates ensemble diversity from single-embed quality as a matter
  of course; chematic has the *data* to do this (§2, finding 1) but
  hasn't published it as a coherent ledger, and conflates a sound engine
  (`embed_pipeline_v2`) with a broken one (`conformer_ensemble()`) under
  one mental model of "3D generation."
- **API sketch**: none — this is a diagnosis round, no new API. Output is
  a written ledger (`validation/results/`-style JSON + a doc write-up),
  not code.
- **Failure taxonomy**: audit whether `PipelineV2FailureCause`'s existing
  12-stage enum is sufficient to classify 100% of the 265-corpus's
  failures, or whether the ledger surfaces failure shapes the enum
  doesn't yet name (e.g. the `conformer_ensemble()` torn-bond defect has
  no home in that enum today, since it's on a different code path
  entirely).
- **Benchmark**: re-run the *existing* `pipeline_v2_vs_rdkit_dump.rs` +
  `pipeline_v2_vs_rdkit_oracle.py` + `gen_pipeline_v2_vs_rdkit_report.py`
  pipeline against current `main` to confirm 241/265 still holds after
  v0.18.0's `assign_c_type` fix (not assumed — this RFC found the fix's
  own blast-radius note doesn't cover this specific re-validation).
- **Fixture / holdout**: `validation/openeye_advantage_fixtures.jsonl` /
  `validation/openeye_advantage_holdout.jsonl` (this RFC's own
  deliverables — see §6).
- **CPU cost**: reuses existing infrastructure; the 265-corpus
  re-measurement is a single, moderate run (minutes, not the "full
  5,000-molecule" scale this round explicitly excludes).
- **WASM feasibility**: N/A (diagnosis only).
- **Breaking-change risk**: none (no code change).
- **Verdict: GO.** See §8 for the full comparison against B1.

### A2 — Conformer Ensemble Core

- **User value**: replaces the currently-unsound `conformer_ensemble()`
  API with one built on the already-sound `embed_pipeline_v2` engine —
  directly fixes the §2 finding 2 defect, not just documents it.
- **Gap vs. OpenEye**: Omega's ensemble generation is deterministic given
  a seed, diversity-pruned, and energy-ranked. Today's chematic has none
  of these three properties on its public ensemble path (unseeded RNG,
  RMSD-only pruning with no energy step, silent minimizer failures).
- **API sketch**: a `ConformerEnsemble` type that calls `embed_pipeline_v2`
  in a loop with `derive_attempt_seed(base_seed, attempt)` (the exact
  deterministic-per-attempt pattern `distance_geometry_v2.rs` already
  uses and tests), collects successes, applies `prune_rms_threshold`
  (already an `EmbedParameters` field, currently unconsumed — this item
  is what would finally consume it), ranks survivors by MMFF94/UFF energy,
  and returns a typed result carrying every attempt's
  `PipelineV2FailureCause` for attempts that didn't produce a usable
  conformer (full provenance, not silent dropping).
- **Failure taxonomy**: extend (don't replace) `PipelineV2FailureCause`
  with an ensemble-level wrapper (e.g. "N of M attempts failed, causes:
  [...]" plus a distinct "insufficient diversity after budget exhausted"
  case) — exact shape is a **NEEDS-RESEARCH** item pending A1's audit of
  what failure shapes actually occur at ensemble scale.
- **Benchmark**: extend the 265-corpus benchmark with an ensemble-mode
  arm; compare against RDKit's existing `etkdgv3_best_of_n` (`BEST_OF_N=
  10`) arm on equal footing for the first time.
- **CPU cost**: N attempts × single-embed cost — bounded by
  `total_timeout_ms`, same cost model as today's `embed_pipeline_v2`
  multiplied by ensemble size, not a new cost category.
- **WASM feasibility**: yes in principle (no new non-WASM-safe
  dependency implied), not evaluated in depth this round.
- **Breaking-change risk**: **high** if `conformer_ensemble()`'s existing
  signature is reused with new semantics (silently different output for
  the same call) — the RFC recommends this ship as a *new* function name
  or an explicit opt-in config field, deprecating the old path openly,
  not a silent behavior swap. Exact API-compat decision deferred to the
  implementation round.
- **Verdict: NEEDS-RESEARCH** until A1 lands (the exact wiring and
  failure-taxonomy shape depend on what A1's audit finds).

**Status (2026-08-23): DONE — merged (PR #371, `cc0c0b1`; hardening pass
`669342a` after a real algorithm-ordering bug was caught in code review).**
A1's research question is answered: `embed_pipeline_v2`'s own
`max_attempts`/`timeout_ms` retry the *same* target conformer, they don't
generate different ones — a genuinely **new outer loop** was required, not
an in-place extension. What actually shipped, vs. this section's original
sketch:

- New function `embed_ensemble_v2(mol, &EnsembleV2Config) ->
  Result<EnsembleV2Result, EnsembleV2ConfigError>` in a new
  `crates/chematic-3d/src/ensemble_v2.rs` — **a new function, not a reused
  `conformer_ensemble()` signature**, resolving this section's own
  "breaking-change risk: high" concern by construction; the legacy API is
  completely untouched.
- Seed derivation reuses `derive_attempt_seed(base_seed, attempt)` exactly
  as sketched (widened to `pub(crate)`, no public API change).
- **Deviation from the sketch**: pruning is NOT done via
  `EmbedParameters::prune_rms_threshold` (that field remains an
  unconsumed, forward-compat-only field, unchanged) — it's done one layer
  up, by the ensemble loop reusing the existing
  `ConformerEnsemble::find_duplicate`/`find_duplicate_symmetric`
  (`conformer.rs`, the latter newly added this round for automorphism-
  aware pruning). In hindsight this is the more natural layering
  (pruning is inherently a multi-conformer, not single-embed, concern),
  but it means the specific mechanism this section predicted did not
  materialize — noted here rather than silently left stale.
- Energy ranking: kept conformers are grouped by `actual_force_field_used`
  and ranked by ascending energy *within* each group only — MMFF94 and UFF
  energies are never compared across a fallback boundary (a correctness
  issue this round's own code review caught before merge: an earlier
  version of the selection algorithm pruned near-duplicates in generation
  order before ranking by energy, silently keeping the first-generated
  representative of an RMSD cluster instead of the lowest-energy one).
- Failure taxonomy: full per-attempt provenance via `ConformerAttempt`
  (`Ok(ConformerSuccess)` with a typed `ConformerDisposition::Kept{..}` /
  `PrunedAsDuplicate{representative_attempt_index, rmsd, symmetric}`, or
  `Err(PipelineV2Failure)`) plus ensemble-level
  `EnsembleTermination::{Completed, TimedOut}` — no ensemble-level
  "N of M failed" wrapper enum was needed; per-attempt detail turned out
  sufficient.
- Benchmark: **not done this round** — no ensemble-mode arm was added to
  the 265-corpus benchmark, and no comparison against RDKit's
  `etkdgv3_best_of_n` was run. Verified instead with unit tests (synthetic
  selection-logic fixtures) and a small, uncommitted real-molecule spot
  check (aspirin/caffeine/cyclododecane/ibuprofen). The best-of-N parity
  gap this RFC named (chematic has zero best-of-N arm in the benchmark)
  is now addressable in principle but not yet exercised.
- **Scope, as planned**: Rust core only, `chematic-3d`. No Python/WASM
  bindings — `embed_ensemble_v2` is not yet callable from Python or WASM.
  That is the natural next step before this item delivers end-user value
  to chematic's largest current audience (Python users), per this RFC's
  own §8 reasoning for why A1/A2 were prioritized over Track B.

### A3 — Shape Runtime

- **User value**: the single largest remaining true "OpenEye is ahead"
  axis after A1/A2 (FastROCS-class shape screening has no chematic
  analog at all today).
- **Gap vs. OpenEye**: complete — no Gaussian-volume overlap, no shape/
  color Tanimoto, no rigid-body shape optimization exists (see gap
  matrix). `usr.rs`/`spectrophores.rs`/`o3a.rs` are real but
  categorically lighter-weight.
- **API sketch** (per the user's own naming, already anticipated in
  `ROADMAP.md`): `shape_tanimoto(mol_a, mol_b)`, `color_tanimoto(mol_a,
  mol_b, feature_model)`, `align_shape(query, target) -> AlignedPose`,
  `screen_shapes(query, library) -> Vec<(idx, ShapeScore)>`.
- **Failure taxonomy**: not designed this round (gated below A1/A2).
- **Benchmark/fixture/holdout**: not designed this round.
- **CPU cost**: Gaussian-volume overlap optimization is meaningfully more
  expensive per pair than 2D Tanimoto or USR — real cost, not evaluated
  in depth this round.
- **WASM feasibility**: plausible for scoring; interactive in-browser
  visualization (mentioned in the user's own item list) is a separate,
  larger UI question.
- **Breaking-change risk**: none (wholly new API surface).
- **Verdict: NOT PLANNED to start before A1/A2 land** — explicit gate,
  matching `ROADMAP.md`'s own existing B-tier framing of this exact item.

### A4 — Microstate Runtime

- **User value**: closes the single cleanest, most complete Quacpac gap
  found in this audit (protonation/ionization-state enumeration doesn't
  exist at all today, not even partially).
- **Gap vs. OpenEye**: complete on enumeration; partial on tautomer
  (existing Phase 2 machinery, still has open residuals per that RFC).
- **API sketch**: an enumeration function returning `Vec<(Molecule,
  ScoreBreakdown)>` across plausible ionization states (reusing `pka.rs`'s
  23-rule table as the *candidate-site* source, not inventing a new
  proton-placement heuristic from scratch), a dominant-microstate selector
  analogous to the existing tautomer-scoring pattern, and explicit
  abstention when no site clearly dominates (matching this project's
  fail-closed convention rather than guessing).
- **Failure taxonomy / benchmark / fixture / holdout**: not designed this
  round.
- **CPU cost**: low — combinatorial over a typically-small candidate-site
  count, similar order to existing tautomer enumeration.
- **WASM feasibility**: yes, no blocking dependency identified.
- **Breaking-change risk**: none (new API surface; does not touch
  `pka.rs`'s existing scalar-prediction functions).
- **Verdict: technically independent of A1-A3** (a molecular-graph
  problem, not a 3D problem) — could in principle run in parallel. This
  round's own instruction is to select **one** starting item across both
  tracks, so A4 is not started even though it has no *technical*
  dependency on A1/A2; this is a sequencing choice, not a blocker, and is
  recorded as such rather than silently implying a false dependency.

### A5 — Docking

- **Verdict: NOT PLANNED** this round or the next. Explicitly gated on
  A1-A4's acceptance criteria all passing. Rationale: docking quality is
  bottlenecked by conformer, shape, and microstate quality in sequence —
  building a docking engine on top of today's unsound ensemble path and
  nonexistent shape/microstate layers would just relabel the existing 3D
  gap as a docking gap, not close it. Receptor preparation (Spruce-
  equivalent, also Missing per the gap matrix) is a further prerequisite
  not even started.

## 4. Track B — Materials runtime

Full gap grounding: `docs/competitive/materials_gap_matrix.md`.

### B1 — Periodic Performance Foundation

- **User value**: removes the one real, disclosed performance ceiling
  in `chematic-crystal` (O(n²) neighbor search) for anyone working with
  structures beyond toy-scale unit cells, and closes the "diagonal-only
  supercell" gap that blocks basic materials workflows (surface slabs,
  defect supercells) even before B5 exists.
- **Gap vs. pymatgen/ASE**: both use cell-list-based neighbor search and
  accept arbitrary integer supercell transformation matrices.
- **API sketch**: replace `neighbors_within`'s internal all-pairs loop
  with a cell-list/bucket structure sized to the cutoff, while preserving
  the *exact* (non-approximate) minimum-image guarantee `periodic::
  minimum_image` already provides — the hard constraint, since a naive
  cell-list built around a cutoff sphere is not automatically equivalent
  to the existing reciprocal-vector-bounded exact search for skewed
  triclinic cells. `make_supercell` gains an arbitrary `[[i32; 3]; 3]`
  matrix overload alongside (not replacing) the existing diagonal
  `[nx,ny,nz]` entry point.
- **Failure taxonomy**: extend `CrystalError` with a typed case for
  supercell matrices that are singular or produce degenerate lattices,
  matching the existing fail-closed convention (`neighbor.rs`'s own
  `MAX_NEIGHBOR_IMAGE_CANDIDATES` cap as the template for "reject
  pathological input explicitly, don't silently truncate").
- **Benchmark**: wall-clock neighbor-search scaling on a range of
  structure sizes (small unit cell → large supercell), before/after,
  reusing the existing triclinic test structures in `lattice.rs`/
  `periodic.rs` as the correctness baseline (same neighbor results,
  faster).
- **Fixture / holdout**: `validation/materials_advantage_fixtures.jsonl`
  / `validation/materials_advantage_holdout.jsonl` (§7).
- **CPU cost**: implementation-only work is a real, non-trivial
  algorithm change; this round's own scope is diagnosis/design, not
  implementation — the cost estimate here is for a *future* round.
- **WASM feasibility**: yes — no new dependency category, pure algorithm
  change over existing types.
- **Breaking-change risk**: low if the new supercell overload is additive
  (new function/method, not a signature change to `make_supercell`);
  neighbor-search results must be bit-for-bit identical to the current
  exact search (a correctness invariant, tested via the existing
  fixtures before any performance work is trusted).
- **Verdict: GO-eligible**, but not selected as this round's recommended
  next start — see §8.

### B2 — Crystal Identity

- **User value**: enables deduplication and canonical comparison of
  crystal structures — the periodic-structure analog of chematic's
  existing canonical-SMILES/Parent-identity work for molecules.
- **Gap vs. pymatgen**: pymatgen's `StructureMatcher` + Niggli reduction +
  primitive-cell finding together let users deduplicate structures from
  different sources; chematic has none of the three.
- **Verdict: NEEDS-RESEARCH**, gated on B1. A canonical-cell definition
  (Niggli-reduced, primitive) is only meaningful once the underlying
  lattice-basis-change math is verified invariant — building B2 first
  would mean re-deriving that invariance work twice.

### B3 — Symmetry

- **User value**: the single most-requested-in-spirit materials-science
  capability (space-group determination is foundational to almost every
  downstream crystallography workflow) — but also, per this audit, the
  single hardest item in either track.
- **Gap vs. spglib**: spglib determines the space group and symmetry
  operations of an *arbitrary* structure with no prior symmetry
  declaration, via a highly-optimized, decades-refined C library.
  chematic today can only *apply* operations a CIF already states in
  text (gap matrix row 6) — it has never attempted the reverse
  (determination) direction at all.
- **Constraint**: the round's hard requirement that the **default build
  stay pure-Rust with zero C/C++ dependencies** rules out simply wrapping
  spglib via FFI as the default implementation path. A from-scratch
  pure-Rust space-group determination algorithm is a substantial,
  multi-month-class undertaking in its own right (this is not
  understated: symmetry-operation search from raw coordinates, point-
  group classification, and Hall/International-Tables-number mapping are
  each nontrivial). An *optional*, non-default FFI backend is a
  structurally different, separately-decidable question, not designed
  here.
- **Verdict: NEEDS-RESEARCH.** This RFC does not recommend starting B3
  as the next round under any framing — it is too large and too risky
  to be a first step, and gated on B1/B2's foundation existing first
  regardless.

### B4 — Materials Analysis

- **Verdict: depends on B1-B3.** Crystal fingerprint, RDF/coordination
  environment, oxidation-state inference, and XRD/diffraction pattern
  calculation all want a stable, canonical structural substrate
  (B1's performance foundation at minimum, arguably B2/B3's identity
  work too) before they're worth building. Not designed further this
  round.

### B5 — Materials Construction

- **Verdict: depends on B1** (arbitrary supercell specifically) as a
  direct prerequisite for slab/surface/defect construction. Not designed
  further this round.

### B6 — Thermodynamics & Simulation Interface

- **User value**: closest of Track B's later items to something already
  partially built — `minimize_mmff94_lbfgs` (real L-BFGS) and `run_md`
  (real velocity-Verlet NVE/NVT) already exist and are well-tested; this
  item is substantially "wire existing molecular building blocks to a
  periodic-aware Calculator abstraction," not starting from zero.
- **Gap vs. ASE**: ASE's `Calculator` trait is the load-bearing
  abstraction letting one MD/optimization driver work with any energy
  backend; chematic has no equivalent trait anywhere in the workspace.
- **API sketch**: a `Calculator`-style trait (energy + forces given a
  `PeriodicStructure`, or a `Molecule` for the non-periodic case) that
  `run_md`/`minimize_mmff94_lbfgs` could eventually be adapted to
  implement — explicitly an **adapter-design** exercise per the round's
  own instruction, not a request to build periodic MD/NEB/phonons in
  full. VASP/LAMMPS I/O already exist (materials gap matrix rows 23, 25)
  and would be reused as-is; Quantum ESPRESSO I/O (row 24) remains a
  real, separate gap this item would need to close.
- **Verdict: depends on B1-B3** for a stable structural substrate to
  build the Calculator abstraction against; adapter design (not full
  MD/NEB/phonon implementation) could reasonably start earlier than the
  rest of Track B, but is not this round's recommendation (see §8).

## 5. Acceptance criteria (design, not new test code this round)

These are the formal acceptance dimensions a future implementation round
must satisfy, each grounded in a chematic pattern that already
demonstrates the general technique — not invented from scratch.

**3D (Track A)**:
- Input-order invariance — same technique as this project's existing
  atom-permutation-invariance tests elsewhere (e.g. tautomer round 2C's
  `build_named` graph-construction helper).
- Deterministic seed reproducibility — template: `distance_geometry_v2.
  rs`'s `same_seed_reproducible`/`different_seed_gives_different_output`
  tests, already passing on the single-embed path; A2 extends this
  pattern to the ensemble level.
- Zero silent stereo inversion — template: `embed_pipeline_v2`'s existing
  `FinalStereoViolation` typed failure and the UFF-rescue-path stereo
  gate found during this audit (§3, OEFF row) — fail closed, never
  silently return inverted stereochemistry.
- Invalid geometry must never report success — template: `UffMinimizeResult.
  sound: bool` (v0.19.0), which already separates "converged" from
  "geometrically sound converged."
- 100% failure-reason classification — template: `PipelineV2FailureCause`'s
  existing 12-stage enum; A1's audit checks whether it already covers
  100% of real failures or needs extension.
- Macrocycle accounted separately from normal rings — template: the
  existing `MACROCYCLE_MIN = 9` boundary and separate macrocycle-rule
  table; A1's ledger reports macrocycle-corpus results as their own row,
  never blended into the overall pass rate.
- Best-of-N never mixed with single-conformer numbers — template: this
  RFC's own §2 finding that today's benchmark already keeps RDKit's
  best-of-10 arm distinct from single-embed arms; A2 must preserve that
  separation for chematic's own new ensemble numbers.

**Shape (Track A, A3, future)**:
- Rigid-transform invariance, symmetry invariance, identical-conformer
  score = 1, non-overlapping-conformer score ≈ 0, analytical vs. grid
  methods verified independently, atom-order independence, CPU/Rayon/
  WASM cross-backend agreement within tolerance — none of these have an
  existing chematic template today (no shape code exists yet); recorded
  here as the bar A3 must clear whenever it starts, not designed further.

**Materials (Track B)**:
- Lattice-basis-change invariance, site-order invariance, periodic-image
  invariance — template: `periodic.rs`'s existing exact minimum-image
  tests and the `cubic_half_cell_tie_uses_lexicographically_smallest_
  image` regression fixture (a real tie-break bug this project already
  found and fixed once in this exact area — B1 must not reintroduce a
  similar tie-break bug when the neighbor-search algorithm changes).
- Primitive/supercell identity — template: none exists yet (primitive
  cells aren't implemented); recorded as B2's bar.
- Skewed triclinic — template: `lattice.rs`'s existing `triclinic_volume`
  test (α=80°, β=95°, γ=110°) as the base case to extend, not invent.
- Partial occupancy — template: `Occupancy::SUM_TOLERANCE` validation
  already in `site.rs`; B2/B4's future work must respect existing
  disorder semantics, not assume full occupancy.
- Automorphism, tolerance-boundary cases, exact-oracle/independent-
  implementation cross-check — no existing chematic template (space-group
  work has never been attempted); recorded as B3's bar, acknowledged as
  the hardest to satisfy of any acceptance dimension in this RFC.

## 6. `validation/openeye_advantage_fixtures.jsonl` / `..._holdout.jsonl`

Design-driving fixtures for A1's ledger and A2's eventual acceptance
tests. Reuse the existing 265-molecule corpus and its oracle/report
scripts as the primary measurement substrate — these two files add a
small number of *new*, hand-picked rows for scenarios the existing corpus
doesn't already isolate cleanly (the live `conformer_ensemble()` defect,
a macrocycle-heavy subset, a best-of-N candidate set), plus a holdout set
checked only after A1's methodology is frozen. See the files themselves
for the row schema and exact molecules.

## 7. `validation/materials_advantage_fixtures.jsonl` / `..._holdout.jsonl`

Design-driving fixtures for B1's acceptance criteria, extending
`chematic-crystal`'s own existing triclinic/occupancy/minimum-image test
patterns (not a new methodology) to the specific invariances B1 must
preserve while changing the neighbor-search algorithm. See the files
themselves for the row schema and exact structures.

## 8. Final recommendation: A1, not B1

Comparing the two GO-eligible candidates the user named, on the three
requested dimensions:

| Dimension | A1 (Conformer Benchmark & Failure Ledger) | B1 (Periodic Performance Foundation) |
|---|---|---|
| **Point/value improvement** | Formalizes an *already-committed* roadmap phase (Phase 6, v0.24.0, already carrying a 99.75/100 target in the existing ladder) with OpenEye-specific competitive framing. Serves chematic's current largest audience (Python/WASM small-molecule users). | Opens an entirely new, currently-unscored axis with no existing phase or point target. Serves a newer, smaller audience — no WASM/MCP crystal bindings exist yet (gap matrix row 26), so B1's performance work wouldn't even be reachable from the browser runtime today. |
| **Technical risk** | Low. Diagnosis-only, reuses fully-existing tooling (265-corpus, oracle scripts, report generator). No new algorithm has to be correct on day one — this round only re-runs and writes up what's already measurable. | Medium (for the *future implementation*, not this design round). Replacing an exact all-pairs search with a cell-list while preserving bit-for-bit exact minimum-image results for skewed triclinic cells is a genuine correctness-preserving algorithm change — this exact class of bug (a triclinic tie-break error) has already bitten this project once in this file. |
| **User value** | Surfaces two concrete, fixable-soon wins independent of the larger OpenEye ambition: a stale public-docs claim (§2.1) and a live, currently-shipping catastrophic bug in a public, documented API (§2.2, dormant since v0.14.0). Both are real problems today, not speculative future value. | No known user-facing complaint or correctness bug exists today for `chematic-crystal`'s neighbor search — it is a scaling concern for large structures, not a currently-wrong-answer problem. Real value, but not urgent value. |

**Recommendation: start A1 next.** It is cheaper, lower-risk, serves the
larger current audience, continues an already-budgeted roadmap phase, and
— independent of any OpenEye-competitive framing — this audit surfaced a
real, live defect in a public API that is worth knowing about regardless
of what happens with the rest of this RFC. B1 remains the correct next
step for Track B whenever materials-science investment resumes, but
should not start in the same round as A1, per the explicit instruction
not to begin two tracks simultaneously.

## 9. What this round does and does not do

**Done this round**: this RFC; the two gap matrices; the four fixture/
holdout JSONL files; a `ROADMAP.md` update recording this round's output
and the Track A/B structure. Zero changes under `crates/*/src/**`.

**Not done this round** (explicit, matching the governing instructions):
no docking implementation; no full OpenEye-compatibility layer; no ML
potential of any kind; no large data downloads; no full 5,000-molecule 3D
remeasurement; no full MMFF94 remeasurement; no Python/WASM binding
additions; no version bump; no release. Ends as a single **draft PR** —
not marked ready, not merged, no tag, no publish.
