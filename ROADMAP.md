# chematic roadmap

> Active plan, revised 2026-09-12. Current release and workspace version:
> v1.0.13. Historical measurements retain the version recorded in each
> artifact.

This file contains only the active priorities and completion boundaries.
Detailed tasks and evidence are maintained in:

- [RDKit accuracy plan](docs/rdkit-accuracy-plan.md)
- [Open-work disposition](docs/roadmap-open-work.md)
- [Validation guide](docs/validation.md)
- [Benchmark index](benchmarks/README.md)
- [Roadmap archive through 2026-09-12](docs/archive/roadmap-through-2026-09-12.md)

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
| A4 — Workflow accuracy | In progress | Pinned historical SMARTS diagnosis: comparable-cell agreement 155,618/155,633 (99.9904%); all-cell 155,618/155,651 (99.9788%). Corrected current-source direct comparison over 16 queries × 5,000 molecules, with match order normalized and unsupported `[kN]` queries separated, reports 145,579 exact matches, 10,042 RDKit-unsupported cells, 18 explicit chematic refusals, and 12 residual mismatches; this supersedes the earlier unqualified 80,000/80,000 claim. An experimental opt-in hybrid `[RN]` selector (shared symmetrized-SSSR with a cage fallback) passes the full SMARTS regression suite and representative bicyclo/adamantane cases; its full-corpus comparison remains unmeasured after an intentionally incomplete diagnostic run. Expanded 51-case/14-template SMIRKS product-set gate passes 51/51 atom-map-aware canonical structure comparisons, with duplicate embeddings normalized as sets; incompatible stereo has a separate 1/1 fail-closed safety gate; reaction presence 6/6; V3000 SGROUP/COLLECTION, isotope, enhanced-stereo boundary, and opaque ENDPTS/ATTACH bond-attribute round trip pass focused regressions | Measure the hybrid lane to completion, triage residuals/refusals under a pinned supported/unsupported SMARTS policy, then broaden standardization, reaction, typed metadata, and RDKit/Indigo V3000 semantic coverage |
| A5 — Independent accuracy adjudication | Local-open / external | Four gold candidates and two placeholders; manifest integrity passes; the paired category/cluster-bootstrap evaluator is implemented and tested, but inputs are already exposed and absolute labels/review are incomplete | Independent absolute gold, unused inputs, non-maintainer review and formal adjudication |
| A6 — 3D accuracy extension | In progress | MMFF94 atom typing improved to 6,681/6,698 exact; excluding the one declared unsupported probe, 6,681/6,697 comparable (99.76%), with 16 residual atoms across 3 macrocycle molecules; BCI charge comparison is 6,665/6,693 comparable (99.58%), with 28 residual atoms. The bounded bond+angle gate passes 265/265. A row-isolated hard-timeout rerun covers all 265 pipeline rows: 215 success and 50 explicit timeout rows; all 215 successes are finite/sound with zero gross clashes, while 54 converge at the 200-step budget (Tier A 42/64 successes, Tier B 12/151). Same-heavy-coordinate RDKit MMFF94 diagnostic over the 215 successes shows median absolute energy delta 29.93 kcal/mol, p90 55.89, and max 920.36. The fixed-H/same-coordinate diagnostic corrected the buffered 14-7 implementation, explicit O-H numeric typing, explicit-H SymmSSSR representation dependence, and delocalized amine N-H typing. A narrow macrocycle-boundary candidate further reduces the maximum comparable energy delta to 8.7467 kcal/mol and raises the within-5 count to 260/262; 0009 is now exact in the macrocycle type audit, while 0029/0030 residuals remain. The latest 262/265 comparable run has median absolute energy delta 0.226 kcal/mol, p90 1.219, max 8.75, 224/262 within 1 kcal/mol, and 260/262 within 5 kcal/mol. Remaining macrocycle/drug-like residuals and convergence/gradient/stereo exits remain open. Evidence: `validation/results/mmff94-hard-timeout-pipeline-a6-v1.0.13.json`, `validation/results/mmff94-rdkit-same-heavy-energy-a6-215-v1.0.13.json`, `validation/results/mmff94-same-explicit-h-energy-a6-265-after-refined-macrocycle-boundary-v1.0.13.json` | Resolve remaining macrocycle/drug-like term residuals, eliminate timeout/non-convergence classes without weakening soundness, then complete fixed-H energy/gradient parity, stereo, and conformer-quality exits |

### Execution order

**Next: A0 (P0), the evidence-contract and genuinely-unused evaluation exits.**
The A1.R cage-like aromatic-ring residual is resolved in the opt-in candidate
lane; its shared-perception adoption still requires the listed cross-binding
regressions. Resume **A0 → A1 remaining exits → A2 → A3 →
A4 → A5 final adjudication → A6**.

1. Complete A0: all-field holdout, trustworthy counts/provenance, negative tests,
   separately built baseline/candidate, and a genuinely unused evaluation set.
2. A1.R acceptance is complete for the opt-in candidate (affected
   perception/fingerprint/SMARTS and binding regressions passed); finish A1 core,
   atom-level potential stereocenters, and unused exits;
   then resolve the A2 CIP/identity residuals.
3. Complete A3 remaining raw/provenance, adopt the validated aromaticity lane
   only after all binding/search reruns, and finish unused-evaluation exits while
   preserving full-k, threshold, exact-order, and independent-oracle results.
4. Finish A4 standardization → SMARTS match sets → reaction products → V3000;
   first triage the 12 current SMARTS residuals and 18 explicit refusals under
   the corrected direct-comparison accounting.
5. Complete A5 independent adjudication, then A6 as a separate 3D profile.

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
[detailed plan](docs/rdkit-accuracy-plan.md) defines scope, task IDs,
corpus splits, failure accounting, numerical tolerances, and candidate exits.

### Next candidate boundary

The next candidate requires A0, A1 core eight-field unused evaluation, all
affected safety/binding regressions, and preserved native defaults.
It may include verified A1/A2/A3/A4 improvements without completing every
package. Full declared 2D compatibility requires A0–A4 exits; independent
chemical equivalence/superiority requires A5. 3D parity requires A6.
No version number or publication action is selected by this plan.

The former additional 1.10x speed stretch remains abandoned and is not a
release gate. New Playground breadth and optional performance work do not
displace the accuracy order above.

## Priority phases

P0–P6 describe product areas. A0–A6 above are the active accuracy work
packages spanning those areas.

| Phase | Purpose | Current disposition |
|---|---|---|
| P0 | Trust and measurement | Local evidence infrastructure exists; full acceptance-contract hardening and independent review remain open |
| P1 | Interchange throughput and safety | Strong bounded parser/streaming coverage; equivalent-operation breadth remains open |
| P2 | Identity and ML primitives | Descriptor and Morgan compatibility are strong; canonical/CIP and symmetry-heavy boundaries remain open |
| P3 | Portable production surface | Shared Rust/Python/Node/WASM contracts exist; broader browser/agent and cross-platform evidence remain open |
| P4 | Chemistry workflows | Typed reaction/document foundations exist; curated reaction/SMARTS quality remains open |
| P5 | 3D and materials | Experimental foundations and bounded gates exist; full force-field quality remains open |
| P6 | Ecosystem durability | Reproducible documentation exists; continuous maintenance and external validation remain ongoing |

## Open phase backlog

Every unchecked item is classified in
[docs/roadmap-open-work.md](docs/roadmap-open-work.md) as `local-open`,
`local-toolchain`, `external`, or `historical`.

- [ ] A0 — Accuracy evidence contract: finish all-field holdout and acquire/freeze the planned 8,000-row sealed evaluation; raw/provenance validation, negative tests, and paired source-built baselines are complete.
- [ ] A1 — Perception and descriptor coverage: accept the bounded A1.R candidate through affected binding/workflow regressions, then complete unused evaluation for potential stereocenters and core fields; adopt additional families separately.
- [ ] A2 — Stereo and identity correctness: close phosphorus CIP and remaining canonical/E/Z coverage without replacing fail-closed outcomes with unsupported confident labels. A2.4's current exposed corpus gate now passes standardize→canonicalize idempotency on 14,999/14,999 rows across descriptor, ChEMBL, and NCI slices; this is regression evidence, not independent correctness or full RDKit string parity. Evidence: `validation/results/canonical-idempotency-a2.4-current-v1.0.13.json`.
- [ ] A3 — Fingerprint and retrieval fidelity: gate all declared k/thresholds, unrounded scores/order, independent RDKit oracle, and unused inputs.
- [ ] A4 — Workflow accuracy: close the 12 current SMARTS residuals and 18 explicit refusals under a pinned supported/unsupported policy, then broaden standardization, reaction-product, and V3000 semantic evidence.
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

See [docs/compatibility-scope.md](docs/compatibility-scope.md).

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
