# Changelog

All notable changes to chematic will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

---

## [Unreleased]

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
  isn't merged yet. A related but distinct conditional trap in `chematic-smarts`'s MCS
  timeout path (already reachable from the shipped `mcs_smiles_json_with_ring_config`
  export) was found during the call-site audit and filed separately as issue #221
  rather than folded into this fix (different crate, different root cause).

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
  figures in `docs/reaction_transform_perf.md` as directional, not precise).
  Does **not** fix the larger, still-open cost on genuinely symmetric
  molecules, which needs automorphism-orbit-aware branch pruning —
  explicitly deferred as future work. Full bisect, methodology, and
  before/after numbers in `docs/reaction_transform_perf.md`.

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
  (`docs/stereo2d_reader_integration_rfc.md`): `chematic_perception::apply_local_parity_from_wedges`
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
  verified it (`docs/ecfp4_bitexact_api_rfc.md`): radius=2 (ECFP4), 2048 bits,
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
- Investigated and **ruled out** a previously-diagnosed "bridged-bicyclic ring-closure ordering" permutation-invariance bug (`docs/canonical_smiles_residual_rfc.md`'s "Root Cause 2") — the two SMILES claimed to be "two spellings of the same molecule" turned out, per independent RDKit `MolToInchi` verification, to be genuinely different constitutional isomers; chematic's differing canonical output for them is correct, not a bug. An additional 22-molecule probe using RDKit-`RenumberAtoms`-generated genuine same-molecule respellings (bridged/spiro/fused/cage systems, including stereocenters and a heteroatom bridgehead) found zero convergence failures. No production code changed for this finding; the RFC and test suite were corrected instead.

### Fixed — `scripts`

- `scripts/canonical_residual_diagnosis.py` called RDKit's `FindMolChiralCenters` with a kwarg name (`useLegacy`) the currently-pinned RDKit version no longer accepts (`useLegacyImplementation`) — a pre-existing tooling bug found while re-verifying the E/Z-marker fix above, fixed as its own small change so it doesn't need an external workaround to run.

## [0.5.0] — 2026-07-23

### Added — `chematic-perception`

- **`local_parity_from_wedges(mol, coords, center)` / `apply_local_parity_from_wedges(mol, coords)`** ("P1-S1a-core") — CIP-independent tetrahedral parity from wedge/hash bonds and 2D coordinates, producing `Atom.chirality` + `Molecule::stereo_neighbor_order` directly. Never calls CIP ranking and never touches `Atom.cip_code` — a CIP tie must not prevent a molecule from having a known local parity. Sign convention measured against RDKit's raw `CHI_TETRAHEDRAL_CW`/`CCW` tag on frame-aligned fixtures (not derived by analogy); handles 3-heavy-plus-implicit-H (no synthetic H position, mirroring RDKit's own `atomChiralTypeFromBondDirPseudo3D`) and multiple simultaneous wedges on one center (each checked in isolation for a consistent parity before the combined volume is trusted — a wedge/hash drawing with two substituents pointing opposite ways is valid notation, not automatically contradictory). Full methodology in `docs/stereo2d_local_parity_calibration.md`. Not yet called from any reader's default parse path, and the SMILES writer's own wedge-token bug (below) had to be found and separately fixed first.
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

- **Hierarchical-digraph CIP (Cahn-Ingold-Prelog) engine** — `assign_cip_accurate_experimental`, a from-scratch, provenance-carrying digraph replacement for the existing shell-pooling comparator in `chematic-chem`. Not yet wired into `chematic_chem::assign_cip()`; a separate, non-default, `publish = false` crate for now. See `docs/cip_accurate_rfc.md` for the full design and milestone history.
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
- **Milestone 5B — opt-in stabilization measurements** (no behavior change): Accurate mode is **~10× slower than legacy** (214–240µs vs 22–24µs per molecule, full 5,000-molecule corpus) — the one still-open item on the default-promotion checklist. Unresolved rate 0.392% (19/4849 stereocenters), **100% traced to the two already-known families** (17 Milestone 4A-2 cage carbons, 2 Milestone 4C-1 phosphorus atoms), zero unexplained causes. Cross-surface parity re-verified at 300 molecules (0 mismatches). Default-promotion gate criteria are listed with current status in `docs/cip_accurate_rfc.md`; promotion itself is not decided by this milestone.

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

[Unreleased]: https://github.com/kent-tokyo/chematic/compare/v0.1.22...HEAD
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
