# Roadmap snapshot before the 2026-09-13 reorganization

Historical snapshot of ROADMAP.md at `b0c24d7c`, after the v1.0.14 release.
Use the [active roadmap](../../ROADMAP.md) for current status and task counts.
The body below is preserved with relative Markdown links rebased; code paths
remain repository-relative. Its 19-item count includes the abandoned speed
stretch, which is excluded from the active backlog. This snapshot does not
establish new measurements or certify the later main merge.

---

# chematic roadmap

> Active plan, revised 2026-09-13. Current release and workspace version:
> v1.0.14. Historical measurements retain the version recorded in each
> artifact.

This file contains only the active priorities and completion boundaries.
Detailed tasks and evidence are maintained in:

- [Trust Release execution plan — T0–T6](../../docs/trust-release-plan.md)
- [RDKit accuracy plan](../../docs/rdkit-accuracy-plan.md)
- [Open-work disposition](../../docs/roadmap-open-work.md)
- [Validation guide](../../docs/validation.md)
- [Benchmark index](../../benchmarks/README.md)
- [Roadmap archive through 2026-09-12](../../docs/archive/roadmap-through-2026-09-12.md)

## North star

Make chematic a dependable, embeddable cheminformatics runtime for Rust,
Python, JavaScript/WASM, and AI agents.

Every compatibility or superiority claim must name the operation, support
domain, comparator version, corpus, configuration, failure policy, and
reproduction command. Unsupported inputs remain visible in coverage; they are
never removed from a denominator to improve a score.

## Current objective

Reach the declared exits for A0 through A6. RDKit-compatible behavior and
independently justified chemical correctness are separate outcomes: agreement
with RDKit is not evidence that chematic is more correct than RDKit.

The next 2–4 weeks prioritize a **1.x Trust Release**: trustworthy comparison
accounting, API-specific compatibility contracts, release/documentation
synchronization, browser/package usability, parser safety, and stereo
invariance. This changes execution order, not the A0–A6 acceptance targets.
Feature breadth and new 3D algorithms are not short-term release requirements.

### Status

- A3.4 adoption rerun: the same 4,500-library/500-query split passes the
  Rust/Python/Node-WASM comparison for k=1/10/100 (1,500/1,500 pairwise
  comparisons, zero mismatches, full-precision scores and deterministic order).
  Oracle, threshold, and cross-binding reports now share source/diff,
  untracked-source, lockfile, and runtime provenance, and use atomic output
  replacement to preserve the last valid report on interruption. The combined
  evidence checker also verifies common corpus/provenance and zero mismatches
  across all three reports.
  See `validation/results/rdkit-search-cross-binding-parity-v1.0.13.json`.

| Package | State | Current evidence | Remaining exit work |
|---|---|---|---|
| A0 — Accuracy evidence contract | In progress | Source-built Python + RDKit 2025.09.3: eight-field diagnostics 40,000/40,000 strict; full-field holdout 96/96 strict; native fixture 4/4; ECFP4 5,000-row Rust/source-built Python/Node-WASM gate passes with module and WASM hashes recorded; 8-field descriptor binding gates pass 5,000/5,000 and 7,737/7,737 exposed rows across Rust/Python/Node-WASM; schema-v2 descriptor diagnostics retain every input row/field, status, delta, and raw hash; raw accounting is independently recomputed for baseline/candidate; baseline/candidate packet validator passes on separately built compatibility baseline (commit 5aae7b23) and current candidate; v2 manifest now records and validates the exposed binding operation plus corpus existence, row count, and SHA-256; a reserved 4-row unused holdout was measured; reproducible dedup/split preparation produces 2,000 development and 7,737 exposed holdout rows, explicitly not sealed; release mode `check_rdkit_accuracy_manifest.py --require-sealed` now fail-closes unless sealed status/splits, frozen digests, and candidate commit/tag metadata exist | Planned genuinely unused 8,000-row sealed evaluation (the 7,737 exposed rows cannot qualify by adding 263 rows), candidate commit/tag freeze |
| A1 — Perception and descriptor coverage | In progress | Eight descriptor compatibility profile is adopted opt-in with source-built 40,000/40,000 strict and full-field holdout 96/96; Fsp3 isotope edge fixed; P=S, N+–O−–N, hypervalent sulfur, and aromatic-oxide environment boundaries are covered; the bounded chordless-aromatic-cycle candidate passes all eight fields at 7,737/7,737 strict against pinned RDKit 2025.09.3, and its promotion gate is `adopted_opt_in` with native defaults preserved; the latest 8-field Rust/Python/Node-WASM binding gate also passes 7,737/7,737; A1.1 now uses physical masses for supported isotope labels and fail-closes unknown labels instead of mass-number approximation; morphine/codeine aromatic-oxide LogP residuals are covered by the RDKit reference suite; after the P-H phosphonate, thioamide-N, and charged-sulfur TPSA fixes, fresh descriptor and ChEMBL rechecks are strict-exact for all eight fields at 5,000/5,000 each, and current-source 7,737-row exposed holdout is strict-exact for all eight fields with parse/unsupported failures 0; affected local perception/chem/fingerprint/WASM/SMARTS/3D regression suites pass (204/864/309/344/185/597, expected ignored tests only); potential centers passed 7,737/7,737 with atom-level FP=0/FN=0 before the shared-perception experiment | Unused/sealed evaluation, workflow-level candidate acceptance, broader A1.1 per-atom corpus, separately gated descriptor families |
| A2 — Stereo and identity correctness | In progress | Canonical semantic probe 200/200; an independent RDKit-InChI structural-correctness rerun passes 5,000/5,000 molecules × 4 randomized valid spellings with 0 failures and 0 oracle-invalid exclusions; current #503 K=1,024 audit reproduces 4/28 divergent components with 0 cross-correspondence failures and pinned source/corpus provenance; three current aromatic-stash residuals have explicit fail-closed and relabeling regression tests; representative diagnosis shows the shared carrier would strip the partner's sole direction while the alternate is a canonical DFS ring-close side, so the solver correctly abstains; unstable phosphorus remains fail-closed in Python and Node/WASM, with a dedicated binding regression passing and no false assignments; fresh label audit on 181 residual rows reports atom-map agreement 155/181, E/Z bond-map agreement 14/14, all 26 atom mismatches confined to P + oracle_unstable, and no assigned-label mismatch | General carrier/DFS representation that resolves the four residuals without choosing an unstable winner, phosphorus CIP adjudication, unused holdout, permutation/idempotency/collision exits |
| A3 — Fingerprint and retrieval fidelity | In progress | Rust/source-built Python/Node-WASM agree for k=1/10/100 on 500 queries × 4,500 rows; 0 pairwise mismatches with unrounded scores; independent exhaustive RDKit oracle also passes all k; 34 success + 1 explicit unsupported-bond fixture passes all bindings; inclusive threshold gate passes 96 boundary cases across Rust/source-built Python/Node-WASM; the RDKit-parity aromaticity lane resolves all 59 radius-1 raw Morgan residuals on the checked-in 5,000-row corpus (5,000/5,000 exact, RDKit 2025.09.3, no preprocessing errors) | Persist raw/provenance packet for the parity lane, re-run every binding and search gate after adoption, and complete genuinely unused evaluation |
| A4 — Workflow accuracy | In progress | Footer-verified persistent-corpus SMARTS lane covers 5,021 molecules × 31 queries (155,651 cells): parity matches 145,588 target-atom sets, 10,042 RDKit SMARTS parse errors, 21 residual cells, and RDKit atom alignment failures 0. The checked-in artifact is `validation/results/rdkit-smarts-direct-chembl-5000-footer-verified-v1.0.13.json`; corpus SHA-256 is recorded inside and `scripts/check_rdkit_smarts_evidence.py` revalidates the source hash, footer, row accounting, alignment, and bucket sum. Earlier private-corpus measurements (12 residuals, or the withdrawn 12→3 partial-output claim) are not mixed with this denominator. The current lane is not an independent frozen baseline/candidate adoption packet. Expanded 51-case/14-template SMIRKS product-set gate passes 51/51 atom-map-aware canonical structure comparisons; duplicate embeddings are sets, with a separate 1/1 incompatible-stereo refusal gate. Reaction presence is 6/6; V3000 SGROUP/COLLECTION, isotope, enhanced-stereo boundaries and opaque ENDPTS/ATTACH round trips have focused regressions | T0: complete independent baseline/candidate comparison and classify the 21 residuals under a pinned policy. Mapped embeddings, standardization, reaction breadth, typed metadata and RDKit/Indigo V3000 semantics remain open |
| A5 — Independent accuracy adjudication | Local-open / external | Four gold candidates and two placeholders; manifest integrity passes; the paired category/cluster-bootstrap evaluator is implemented and tested, but inputs are already exposed and absolute labels/review are incomplete | Independent absolute gold, unused inputs, non-maintainer review and formal adjudication |
| A6 — 3D accuracy extension | In progress | MMFF94 atom typing improved to 6,681/6,698 exact; excluding the one declared unsupported probe, 6,681/6,697 comparable (99.76%), with 16 residual atoms across 3 macrocycle molecules; BCI charge comparison is 6,665/6,693 comparable (99.58%), with 28 residual atoms. The bounded bond+angle gate passes 265/265. A row-isolated hard-timeout rerun covers all 265 pipeline rows: 215 success and 50 explicit timeout rows; all 215 successes are finite/sound with zero gross clashes, while 54 converge at the 200-step budget (Tier A 42/64 successes, Tier B 12/151). Same-heavy-coordinate RDKit MMFF94 diagnostic over the 215 successes shows median absolute energy delta 29.93 kcal/mol, p90 55.89, and max 920.36. The fixed-H/same-coordinate diagnostic corrected the buffered 14-7 implementation, explicit O-H numeric typing, explicit-H SymmSSSR representation dependence, and delocalized amine N-H typing. A narrow macrocycle-boundary candidate further reduces the maximum comparable energy delta to 8.7467 kcal/mol and raises the within-5 count to 260/262; 0009 is now exact in the macrocycle type audit, while 0029/0030 residuals remain. The latest 262/265 comparable run has median absolute energy delta 0.226 kcal/mol, p90 1.219, max 8.75, 224/262 within 1 kcal/mol, and 260/262 within 5 kcal/mol. Remaining macrocycle/drug-like residuals and convergence/gradient/stereo exits remain open. Evidence: `validation/results/mmff94-hard-timeout-pipeline-a6-v1.0.13.json`, `validation/results/mmff94-rdkit-same-heavy-energy-a6-215-v1.0.13.json`, `validation/results/mmff94-same-explicit-h-energy-a6-265-after-refined-macrocycle-boundary-v1.0.13.json` | Resolve remaining macrocycle/drug-like term residuals, eliminate timeout/non-convergence classes without weakening soundness, then complete fixed-H energy/gradient parity, stereo, and conformer-quality exits |

### Execution order

**Next: T0, the A0/A4 evidence and resource-budget audit (priority P0).**
The earlier shared-ring claim of 12 → 3 SMARTS mismatches is withdrawn: it
counted a partial dump. The current durable lane is the footer-verified
5,021-molecule corpus artifact above. The shared-ring implementation now
propagates the caller's candidate budget and reports measured exploration
counts; bounded-perception and zero-budget regressions cover that contract.
The remaining T0 work is independent baseline/candidate provenance, residual
classification, and adoption policy. See the detailed plan for counts, hashes,
matching semantics, and negative tests.

| Order | Priority | Work ID and output | Existing Phase / accuracy package |
|---|---|---|---|
| 1 | P0 | T0: complete-run accounting, measured budgets, reproducible before/after | P0/P2; A0/A4 |
| 2 | P0 | T1: API/profile Compatibility Contract and generated dashboard | P0/P2/P3; A0–A4 |
| 3 | P0 | T2: release metadata, public-channel audit, package and migration docs | P0/P6; A0 |
| 4 | P1 | T3: WASM/Worker/MCP install and 10,000-molecule usability gates | P3; A3/A4 |
| 4 | P1 | T4: sourced malicious-input corpus and isolated resource-limit CI | P1/P3; A0/A4 |
| 4 | P1 | T5: stereo permutations, round trips, confident-error and collision gates | P2/P4; A2/A5 |
| 5 | P2 | T6: independent review, maintenance, existing 3D accuracy gaps | P5/P6; A5/A6 |

T0–T6 are delivery IDs; priority P0/P1/P2 is urgency, not the product Phase
number. The [detailed plan](../../docs/trust-release-plan.md) defines subtasks,
dependencies, provisional resource budgets, and a 12-week review schedule.
Public-channel inventory, genuinely unused A0 inputs, and A5 gold preparation
start in parallel. T4 safety defects preempt other work when confirmed.

The current bounded chordless-aromatic-cycle candidate resolves the cage-like
aromatic-ring residual: all eight fields are 7,737/7,737 against pinned RDKit
2025.09.3, with native defaults preserved and promotion `adopted_opt_in`.
This is not yet a full release acceptance: the affected cross-binding and
workflow gates, plus unused evaluation, remain mandatory. The expanded
potential-center result still belongs to its earlier candidate lane and must
not be silently reused as evidence for unrelated shared-perception changes.

The exposed 7,737 rows are development/regression evidence. Adding 263 rows
would only reach a numerical total of 8,000; it would not make this corpus
unused or sealed. Obtain and freeze genuinely unused inputs with provenance,
overlap checks, and a predeclared evaluation protocol for the planned gate.

A1 additional descriptor families and A5 local preparation can run in parallel.
Extra descriptors do not delay A2. A0 and A3 retain their narrower passing
results, but are reopened against their original full acceptance scope.
A5 local evaluator/data work is open; only independent review requires an
external reviewer.

### Acceptance targets

- Core eight descriptors: zero strict mismatches and 100% coverage on existing
  and unused inputs; integer fields exact, floating fields within 1e-6.
- Stereo/identity: atom-level FP/FN=0, no wrong confident label, no structure or
  stereo loss, no false key merge, and permutation/idempotency gates passed.
- Exact retrieval: every query has recall=1.0 and identical ordered IDs for
  k=1/10/100 and thresholds; unrounded score error at most 1e-12.
- Workflows: exact mapped match/product sets and declared metadata preservation.
- Independent accuracy: reviewed gold with a predeclared paired 95% confidence
  interval; equivalence requires the whole accuracy-difference interval inside
  ±0.1 percentage point, superiority requires a positive lower bound plus
  coverage safeguards. Insufficient evidence is not equivalence.
- 3D: separately frozen typing/charge/energy/gradient and conformer-quality gates.

These are targets, not new measurements. The
[detailed plan](../../docs/rdkit-accuracy-plan.md) defines scope, task IDs,
corpus splits, failure accounting, numerical tolerances, and candidate exits.

### Next candidate boundary

The next candidate requires A0, A1 core eight-field unused evaluation, all
affected safety/binding regressions, and preserved native defaults.
It may include verified A1/A2/A3/A4 improvements without completing every
package. Full declared 2D compatibility requires A0–A4 exits; independent
chemical equivalence/superiority requires A5. 3D parity requires A6.
No version number or publication action is selected by this plan.

The Trust RC additionally requires T0, the T1 declared-operation contract,
T2 candidate metadata/docs checks, T3's 10,000-row install/Worker gate,
T4's fixed safety corpus, and T5's published suite and existing safety gates.
Week 4 is an RC **audit**, not a promised release date. A missing original
A0/A1 exit postpones the RC; a development snapshot must not be relabeled as
accepted. Public-channel synchronization is checked after publication and is
distinct from prepublication candidate readiness.

The former additional 1.10x speed stretch remains abandoned and is not a
release gate. New Playground breadth and optional performance work do not
displace the accuracy order above.

## Priority phases

P0–P6 describe product areas. A0–A6 above are the active accuracy work
packages spanning those areas.

| Phase | Purpose | Current disposition |
|---|---|---|
| P0 | Trust and measurement | T0 now propagates the shared-ring exploration budget, records measured candidate counts, rejects incomplete SMARTS dumps, and labels RDKit descriptor lanes; sealed evaluation and independent review remain open |
| P1 | Interchange throughput and safety | Bounded parser/streaming and V3000 opaque round-trip gates pass their declared local scopes; equivalent-operation breadth and external-engine fixtures remain open |
| P2 | Identity and ML primitives | Descriptor/Morgan compatibility and canonical fail-closed gates are strong; symmetry-heavy SMARTS budget accounting is hardened, while CIP/E-Z residuals and full independent correctness remain open |
| P3 | Portable production surface | Shared Rust/Python/Node/WASM contracts exist; broader browser/agent and cross-platform evidence remain open |
| P4 | Chemistry workflows | Typed reaction/document foundations exist; curated reaction/SMARTS quality remains open |
| P5 | 3D and materials | Experimental foundations and bounded gates exist; full force-field quality remains open |
| P6 | Ecosystem durability | Reproducible documentation exists; continuous maintenance and external validation remain ongoing |

## Open phase backlog

Every unchecked item is classified in
[docs/roadmap-open-work.md](../../docs/roadmap-open-work.md) as `local-open`,
`local-toolchain`, `external`, or `historical`.

- [ ] A0 — Accuracy evidence contract: finish all-field holdout and acquire/freeze the planned 8,000-row sealed evaluation; raw/provenance validation, negative tests, and paired source-built baselines are complete.
- [ ] A1 — Perception and descriptor coverage: accept the bounded A1.R candidate through affected binding/workflow regressions, then complete unused evaluation for potential stereocenters and core fields; adopt additional families separately.
- [ ] A2 — Stereo and identity correctness: close phosphorus CIP and remaining canonical/E/Z coverage without replacing fail-closed outcomes with unsupported confident labels. A2.4's current exposed corpus gate now passes standardize→canonicalize idempotency on 14,999/14,999 rows across descriptor, ChEMBL, and NCI slices; this is regression evidence, not independent correctness or full RDKit string parity. Evidence: `validation/results/canonical-idempotency-a2.4-current-v1.0.13.json`.
- [ ] A3 — Fingerprint and retrieval fidelity: gate all declared k/thresholds, unrounded scores/order, independent RDKit oracle, and unused inputs.
- [ ] A4 — Workflow accuracy: classify the 21 SMARTS residual cells and 10,042 RDKit SMARTS parse errors under a pinned supported/unsupported policy, then broaden standardization, reaction-product, and V3000 semantic evidence. The current 5,021-molecule × 31-query lane is evidence, not a frozen independent adoption packet.
- [ ] A5 — Independent accuracy adjudication: build absolute gold and the statistical evaluator, freeze genuinely unused inputs, obtain non-maintainer decisions, and meet the equivalence/superiority criteria.
- [ ] A6 — 3D accuracy extension: complete declared MMFF94/UFF typing, charge, energy, gradient, stereo, soundness, and conformer-quality gates.
- [ ] Measure only equivalent operations against installed RDKit and Open Babel versions on identical inputs; keep subprocess-only and non-equivalent lanes non-ranking.
- [ ] Exact canonical-SMILES and cross-engine V3000 parity: preserve semantic equality and unsupported richness while expanding independent RDKit/Indigo fixtures.
- [ ] Confirm the historical 1.10x performance stretch only if it is explicitly reactivated; it is currently abandoned and non-blocking.
- [ ] Replace remaining MD/UFF/MMFF94 finite-difference production paths only after analytic energy/gradient/stereo soundness gates cover the same domain.
- [ ] Replace remaining periodic neighbor all-pairs evaluation with validated cell-list/cutoff paths and exact result-set parity.
- [ ] Optimize symmetry-heavy canonical and SMARTS search only after exact-output, invariance, and budget-exhaustion gates are broad enough.
- [ ] Extend browser and agent adversarial cases for cancellation, malformed records, limits, stable JSON errors, and supported browser engines.
- [ ] Expand reaction/SMARTS/medicinal-chemistry breadth only after shared P0–P3 primitives have current evidence.
- [ ] Add curated reaction/query precision, recall, invalid-product, timeout, and ambiguity reports with independent oracles.
- [ ] Close MMFF94/UFF typing, charge, parameter, convergence, and stereo gaps using same-coordinate and soundness-aware comparisons.
- [ ] S5 independent gate: obtain a non-maintainer review or external audit of major parser, serialization, and binding boundaries.
- [ ] S6 continuous maintenance: rehearse advisory intake, fix, backport, artifact publication, and supported-version synchronization.

## Security and release gates

A candidate is releasable only when all affected scopes have:

- focused tests and the relevant workspace checks;
- format and diff checks;
- regenerated machine-readable evidence with source/comparator versions;
- documentation that states support and failure boundaries;
- no unclassified mismatch, silent fallback, non-finite metric, or hidden input
  exclusion;
- stable public defaults unless a migration and release-scoped review explicitly
  approve a change.

External gates remain open until their external artifacts exist. Local packet
preparation is not independent review.

## Compatibility boundaries

The v1 contract remains unchanged:

- `canonical_smiles()` is a representation; stable-key APIs remain fail-closed.
- Aromaticity and CIP expose explicit native/compatibility models rather than
  promising universal RDKit parity.
- Native ECFP4 and RDKit-compatible Morgan are distinct profiles.
- Python `RWMol`, CDXML editing, and Markush/polymer expansion are bounded
  subsets.
- 3D/MMFF94 remain Experimental until A6 exits pass.
- Pure-Rust InChI is approximate; exact standard InChI is opt-in native FFI.

See [docs/compatibility-scope.md](../../docs/compatibility-scope.md).

## Definition of complete

A work package is complete only when its implementation, tests, documentation,
and required measurement agree across the declared support domain. A narrow
smoke test, successful parsing, local self-review, or agreement on supported
rows alone cannot complete a broader package.

## Definition of ahead

chematic is ahead of a competitor only for a named use case and measured
dimension. “Faster”, “more accurate”, “smaller”, and “more compatible” are
hypotheses unless the claim identifies versions, corpus, configuration,
hardware, failure policy, uncertainty where applicable, and a reproduction
command.
