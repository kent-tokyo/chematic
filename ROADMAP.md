# chematic roadmap

> Updated 2026-09-23. Released: **v1.0.20**. Next candidate: **v1.0.21**.
> **A0 core-eight multi-field acceptance is complete.** Preserve every exposed
> cohort and move the active correctness focus to A2 stereo/identity and the
> next version-pinned RDKit rebaseline.
> Published v1.0.20 Parse + Morgan is verified from npm: 1.398x parse-inclusive
> and 3.511x prepared, with 9,999/9,999 supported rows exact. Its PyPI wheel
> also confirms the MMFF94 stereo-safe 265/265 quality result; the remaining
> A6 exits are energy/term, timeout, conformer quality, and MMFF94 speed.

> The September 20 competitor-watch intake confirms the direction rather than
> starting a feature race: complete the Trust Release packet; stage the next
> RDKit rebaseline; turn *verified* upstream failures into narrow regressions;
> and make batch terminal accounting a binding contract. Reported upstream
> changes remain watch inputs until their primary source, released artifact, and
> affected API surface are recorded.

## Current status

- Issue #70's process-level Criterion gate is calibrated and complete. Hosted
  evidence blocks the +10% fixture, detects the +5% fixture in 8/10 independent
  runs, and turns a material one-sided runner bias into
  `environment-inconclusive` rather than a bare regression failure. The null
  control now uses ten blocks (seven cannot cross the configured Wilson
  boundary) plus a 1.04 practical-effect floor. Evidence:
  `validation/results/criterion-gate-hosted-calibration-2026-09-22.json`.
- Issue #149 is complete for the measured shared-carrier canonical E/Z scope.
  On clean `main` commit `2fe0cdd7`, all 28 coupled components converge across
  1,024 seeded atom relabelings each, with zero canonical divergence and zero
  correspondence failure. All four historical residuals are pinned as
  equivalent-spelling regressions; unmeasured coupled shapes remain fail-closed.
- Issue #618 is complete for the bounded atom-level pseudoatom contract:
  MOL V2000/V3000, SDF, and CDXML preserve real carbon, `*`, `R`, and
  `R1..R9999`; unsupported labels reject instead of becoming placeholder
  carbon, and WASM exposes index-aligned inspection plus immutable edits.
  This is interchange metadata, not full Markush expansion semantics.
- T1.5 now has an executable and CI-checked historical RDKit 2026.03.6 packet.
  Exact Python wheel and official npm tarball hashes, the measured Boost.Python
  backend, binding behavior, commands, host/toolchain, and all 10,000 exposed
  rows are retained. Public CheMatic 1.0.19 is Morgan-exact on 9,999/10,000,
  CIP-exact on 9,770/10,000 after endpoint-correspondence correction, and has
  14,306 differing SMARTS cells out of 310,000. Residual adjudication identifies
  all 18 SMILES semantic differences
  as CheMatic stereo-writer regressions ([#632](https://github.com/kent-tokyo/chematic/issues/632)) with non-isomeric graph identity
  preserved, and the one Morgan difference as a typed unsupported coordination
  contract. The 230 CIP rows and 3,364 SMARTS rows remain unresolved and are
  tracked by [#634](https://github.com/kent-tokyo/chematic/issues/634) and
  [#635](https://github.com/kent-tokyo/chematic/issues/635). The older
  9,880 CIP figure compared engine-local bond indices and is withdrawn. These are
  compatibility diagnostics, not a new-oracle adoption or complete-parity claim.
- The post-v1.0.19 candidate `5e9211a6`, frozen as annotated tag
  `trust-eval-candidate-a0-multifield-r3-20260922`, completes the A0 core-eight
  gate against RDKit 2025.09.3. A separately acquired source was overlap-audited
  against every previously exposed source; development passed 2,000/2,000 and
  the one-time sealed holdout passed all eight fields at 8,000/8,000 with zero
  mismatches, parse failures, or unsupported values. The commit-safe result is
  `validation/results/a0-core-eight-sealed-acceptance-20260922.json`.
- The published-package browser scorecard is now in `main` via
  [PR #541](https://github.com/kent-tokyo/chematic/pull/541) (`8ce8026`). It
  compares `@kent-tokyo/chematic@1.0.15` with official
  `@rdkit/rdkit@2026.3.6` (runtime `2026.03.6`) on a fixed exposed 10k corpus. It records a smaller
  WASM artifact, lower local no-store ready and parse/write timings, a
  RDKit-faster parse-inclusive fingerprint lane, one typed Fe(II)
  coordination refusal, and a three-run process-tree RSS diagnostic. Unique-memory,
  remote-network, and cross-host acceptance remain unmeasured.
  The local-ready measurement is not a CDN or cross-host claim.
- [PR #542](https://github.com/kent-tokyo/chematic/pull/542) (`3d61132`)
  updates the locked `rustls` dependency to 0.23.45. Its hosted Security Audit
  passed the local review, isolated parser-security corpus, and Cargo Audit;
  the prior RUSTSEC-2026-0285 finding is resolved for the merged lockfile.
- [PR #544](https://github.com/kent-tokyo/chematic/pull/544) (`7c5307a`)
  establishes the T1.7 ordinary-V3000 baseline: version-pinned RDKit and Indigo
  semantic round trips, SGROUP ordering/edit preservation, and bounded WASM SGROUP
  inspection. Its contract explicitly excludes coordination chemistry and ENDPTS/
  ATTACH interpretation; opaque retention is not presented as editable semantics.
- The public descriptor classification remains 52/52 strict for TPSA and its
  declared cross-descriptor fields. A separate TPSA atom-type probe is now
  63/63 strict and covers RDKit's rare three-membered-ring, protonated imine,
  ammonium, nitrilium, aromatic-cation, and Kekulé-N boundaries. Two exposed
  public 5k corpora are 10,000/10,000 strict for TPSA. These are non-sealed
  development regressions; they neither reuse nor rehabilitate any rejected
  8k candidate.
- The earlier TPSA-only pass remains historical. The completed multi-field run
  additionally exposed and fixed bounded RDKit-parity aromaticity boundaries
  for fused sulfur/nitrogen systems, a zero-aromatic fused ring system, and one
  fused cyclic-ether Crippen type. Its raw rows and per-row output are local-only
  and now exposed; no later candidate may reuse them as an unused holdout.
- The bounded V3000 E/Z query gate preserves RDKit's four-query/four-target
  `HasSubstructMatch(..., useChirality=True)` table through a CheMatic round
  trip.  Indigo 1.46.0 matches both targets for both queries in that fixture;
  record that observed behavior as an external semantic loss, not as CheMatic
  query compatibility or an Indigo correctness verdict.
- **A0 is complete; A1–A6 retain open acceptance work.** This does not complete
  the broader Trust RC, independent gold evaluation, or 3D accuracy program.
- The three T1.6 sealed candidates are historical, exposed evaluation inputs:
  `trust-eval-candidate-20260916` was rejected for molecular-weight/TPSA
  failures, and `trust-eval-candidate-20260920` was rejected at 7,998/8,000
  strict molecular-weight agreement (two isotope-table residuals), and
  `trust-eval-candidate-20260920c` was rejected with 7/8 fields at 8,000/8,000
  strict but TPSA at 7,977/8,000. None may be used for tuning or to validate a
  later candidate; public summaries deliberately exclude raw molecular input.
- PR #555 (`7d98dcd3`) merges the Parse + compatible Morgan optimization with
  supported-domain bit agreement and local/hosted browser speed evidence.
  It is a different candidate from the September 16 accuracy freeze. Published
  package timing, the original multi-session protocol, and broader runtime
  resource acceptance still require their own evidence.
- This index tracks **17 open areas**: six accuracy packages and eleven
  cross-cutting follow-ups. They overlap; they are not 18 sequential phases.
  The abandoned 1.10x speed stretch is historical and excluded from that count.

Detailed measurements retain their original version, source, corpus, and
comparator. The [2026-09-13 snapshot](docs/archive/roadmap-through-2026-09-13.md)
preserves the previous long-form status; no chemistry benchmark was rerun for
this reorganization.

## Competitive direction

Win a bounded workflow: **install a typed package, process private chemistry
data locally, and obtain correct, explainable results without blocking the
browser**. Rust and WASM alone are not differentiation. The position is a safe,
typed, local-first chemistry kernel for browser, Rust, Python, Node, and agent
workflows, with a published compatibility contract for every claimed operation.

| Competitor pressure | chematic response | Explicit non-goal for this release theme |
|---|---|---|
| RDKit npm/WASM and Python typing improve usability | Compare pinned public artifacts; publish Worker/cancellation/error contracts and version-specific wins and losses | Compete on installation friction alone, infer runtime safety from type annotations, or generalize one browser speed win |
| RDKit CIP and SMILES fixes expose order-sensitive chemistry risk | Treat atom-order and SMILES-spelling permutations as release regressions; retain uncertain cases as typed refusals | Treat agreement for one spelling as stereochemical proof |
| Indigo advances V3000, coordination, and polymers | Keep version-pinned source/output truth tables beside actual Indigo → chematic → RDKit and reverse round trips; separate opaque retention, observed external loss, and editable typed semantics | Start a haptic/polymer feature race or call retained bytes query compatibility before semantics are proven |
| CDK relaxed parsing helps interactive editors | Keep normal parsing strict; any preview mode must be a separately named incomplete-state API | Silently accept incomplete chemistry in normal APIs |
| Open Babel parser hardening raises the security bar | Publish fixed malformed-input corpus accounting with time/memory limits and refusal reasons | Treat Rust's memory model alone as parser-security proof |
| COSMolKit and other Pure-Rust parity libraries compete for the same headline | Preserve every input's index, terminal outcome and error stage across bindings; publish reproducible runtime and release evidence | Infer complete batch success from zero export failures or enter a feature-breadth race |

The September 20 review turns upstream changes into scoped regression work,
not an automatic promise of new chemistry support. The
[current execution packet](docs/trust-release-plan.md#7-2026-09-21-trust完了に向けた実行順)
records primary sources, implementation gaps, dependencies, and acceptance.
Existing published 1.0.15 measurements stay historical; new release/candidate
measurements must identify their actual package/source hashes.

### Trust Release execution sequence (September 21 revision)

The current upstream watch does not open a new feature race. It makes the
following already-scoped work more urgent. Each stage has a distinct output,
and later stages must not silently turn a rejected or exposed accuracy cohort
into fresh validation evidence.

| Stage | Work | Completion evidence | Explicit boundary |
|---|---|---|---|
| 0 | Preserve rejected A0 candidates and close the core-eight acceptance packet | Candidate provenance, declared support domain, failures, affected bindings, and the accepted one-time result are linked from the ledger | Complete. Candidate `5e9211a6` passed all eight fields on 8,000/8,000 sealed rows after a new freeze, post-freeze acquisition, and overlap audit. `scripts/check_a0_core_eight_sealed.py` validates the commit-safe summary; all raw cohorts are exposed and prohibited from reuse. |
| 1 | Prepare the RDKit rebaseline lanes before a new stable artifact is selected | Historical 2026.03.6 Python/npm artifacts, the measured Boost.Python backend, typed-error probes, fixed operation settings, complete rows, and reproducible commands are committed and CI-checked | Keep 2026.03.6 historical. Run an old/new comparison only after the next official artifact is released; do not infer nanobind or npm availability from upstream work |
| 2 | Import verified upstream failure modes as CheMatic regressions | Identity-renumber stereo, selected-atom/bond CIP, and V3000 E/Z query-round-trip truth tables pass, or a typed unsupported outcome is documented | A competitor failure is a regression source, not proof that CheMatic is more correct; do not import a report before confirming its reproducer and affected versions |
| 3 | Close A2/#503's canonical E/Z work before new 3D breadth | Remaining components converge under atom-order and spelling permutations, reparse preserves stereo, and stable-key stops only where ambiguity remains real | Issue #149's measured residual is complete. The bounded complete-slot planner includes raw directional carriers, reparses every candidate against a non-recursive E/Z signature, and selects the lexicographic minimum without raw atom/bond-index tie-breaks. A clean-main audit now passes all 28 coupled components across 1,024 relabelings each with zero divergence or correspondence failure. Stable keys remain limited to the proven aromatic-stash path; unmeasured coupled systems remain fail-closed. |
| 4 | Promote BatchResult accounting from a convention to a cross-binding contract | Rust, Python, Node, and WASM expose input count, original index, stage, terminal outcome, skipped/refused accounting, and `all_succeeded`-equivalent behavior; cancellation exposes the unprocessed range | `failed == 0` alone never means every input was exported or processed |
| 5 | Publish the browser/runtime proof after the contract is complete | Worker, cancellation, memory/time limits, structured errors, local-only behavior, and package scorecard are measured on the declared workload; performance is always split into parse/write, compatible Morgan, and native ECFP lanes | Artifact size, TypeScript declarations, or one cold-start result alone is not a browser-DX or resource-safety claim |

Stages 0–2 may proceed in parallel. Stage 3 is the first local correctness
implementation priority after immediate regressions; Stage 4 may proceed in
parallel once its public contract is fixed. Stages 5 and any new Trust RC are
gated on their preceding evidence rather than on release cadence.

### Near-term execution contract (September 21)

This is the operational order for the next one to three months. It replaces
neither the A0–A6 exits nor the T0–T6 plan; it makes their dependencies and
stop conditions explicit.

| Order | Work package | Concrete next output | Acceptance / stop condition |
|---:|---|---|---|
| 1 | **A0 evidence maintenance** | Keep the accepted core-eight packet and historical rejection ledger immutable | Reopen A0 only when the declared field/profile/oracle changes; any new adoption still requires a new source, overlap audit, annotated freeze, attestation, and one-time protocol. |
| 2 | **T1.5 RDKit rebaseline preparation** | Maintain the runnable 2026.03.6 Python/npm packet and add a native lane or explicit unavailable record; retain complete-row difference classes | Do not label a “next RDKit” comparison until its official artifact is verified. Historical 2025.09.3/2026.03.6 lanes stay visible. |
| 3 | **T5.7/T1.9 semantic regressions** | Add source-pinned reproductions for identity-renumber stereo and V3000 E/Z query round trips; include selected atom/bond CIP only after its source fixture is verified | Each fixture must compare atom/bond correspondence and semantic truth tables, or return a typed unsupported result. A competitor defect is never evidence of CheMatic superiority. |
| 4 | **T3.7 controlled batch runtime** | Specify a versioned unknown-stream envelope for Rust, Python, Node, WASM, Worker, and MCP adapters | Known-length accounting remains exact; unknown streams expose processed prefix, observed-but-unprocessed data, unread range/unknown, terminal reason, and cancellation state. Zero `failed` must not imply success. |
| 5 | **Browser proof and performance maintenance** | Measure the declared Worker/stream/cancel/export path and rerun the package scorecard only with equivalent operation contracts | Keep parse/write, compatible Morgan, native ECFP, startup, and memory as separate lanes. A local browser result is not a cross-host, CDN, or unique-memory claim. |
| 6 | **A6 correctness, not 3D breadth** | Close existing energy, term, timeout, and conformer-quality gaps | Public v1.0.20 confirms 265/265 independently sound, stereo-clean and clash-free MMFF94 stereo-safe outputs versus RDKit's 264/265. This quality lane remains slower at 0.944x paired geometric speed (0.861x 95% lower bound). Public UFF best-of-10 is 265/265 usable and 2.359x faster (lower 2.198x). Do not add embedding breadth or promote a complete 3D win before the remaining numerical and quality exits close. |

## Priority order

| Order | Priority | Delivery and next output | Product Phase / accuracy |
|---|---|---|---|
| 1 | P0 preparation; external execution | **T1.5 rebaseline**: preserve the completed 2026.03.6 Python/npm execution packet. The 18 SMILES-semantic residuals are classified as CheMatic stereo-writer regressions tracked by #632 and the one Morgan residual as a typed unsupported contract. CIP now compares proven atom order plus bond endpoints rather than engine-local bond indices, yielding 9,770/10,000 exact and 230 unresolved rows. The independent native/C++ lane is explicitly unavailable with release evidence. Product adjudication continues separately in #634/#635; prepare the same commands for the next published stable artifacts and re-measure Python boundary overhead after an actual backend change | P0/P2/P3; A1–A4 |
| 2 | P1 | **A2 and semantic interchange closure**: retain the deployed #503 aromatic-stash planner and the completed #149 28-component × 1,024-relabeling gate; expand acceptance only with new bounded semantic-reparse evidence. Resolve phosphorus adjudication before new 3D breadth. Retain the completed T5.7 identity-renumber gate and T1.9 V3000 E/Z query truth table; after primary-source verification add selected-atom/bond CIP and identity-renumber stereo cases, keep T5.6 selective CIP/isotope/atrop and T1.8 attachment boundaries, then add typed semantics before making compatibility claims | P1/P2/P4; A2/A5 |
| 3 | P1 | **T3.7 controlled runtime**: retain the merged Explorer unknown-length CSV cancellation adapter, then extend the completed row-accounting contract through Worker/stream/cancel/export/limits on 10k inputs. Every binding must report original index, stage, success/failed/refused/skipped terminal outcome, unprocessed cancellation range, and an all-inputs-succeeded predicate; publish typed install examples and a generated package scorecard only after that contract and its resource evidence exist. T3.6 maintains its scoped speed win; it does not substitute for runtime controls | P0/P3/P6; A3/A4 |
| 4 | P2; P0 for confirmed incorrectness | **T6 / A5–A6**: retain the v1.0.20 public-package 265/265 stereo-safe quality result. Resolve the two same-coordinate energy residuals, term-level adjudication, timeout, conformer-quality exits, and MMFF94 quality-lane speed deficit. Do not collapse lightweight speed and stereo-safe quality lanes into one claim | P5/P6; A5/A6 |

T3/T4/T5 can proceed in parallel once their dependencies pass. Public-channel
inventory and A5 gold preparation continue alongside frozen-packet verification.
Confirmed silent corruption, panic, or resource-limit defects take precedence.

### Parse + ECFP4/Morgan performance workstream

The [detailed T3.6 plan](docs/parse-morgan-performance-plan.md) owns the contract,
profiling hypotheses, six steps **PF0 → PF1 → PF2 → PF3 → PF4 → PF5**, and adoption gates.
PF IDs are subtasks, not new P7+ product phases. This is a new explicit performance
goal; the abandoned additional 1.10x SMILES stretch stays abandoned.

The 2026-09-23 registry-installed v1.0.20 Chromium run passes both compatible-
Morgan lanes: parse-inclusive **1.398x** geometric speedup (95% lower bound
**1.363x**) and prepared **3.511x** (lower bound **3.407x**), with
9,999/9,999 supported rows exact. The earlier v1.0.19/source measurements remain
historical rather than being relabelled.

Retained acceptance target: fixed radius-2/2048-bit compatible Morgan, same input/output,
full supported-domain bit agreement and no new refusals; paired speedup 95% CI
lower bound >1.0 and process-mean p95 no worse than RDKit. Native ECFP speed is
reported separately, never substituted for Morgan compatibility. Independent
performance inputs, cross-binding regressions and environment-specific reports
are mandatory; the sealed accuracy set is not optimization data.
The source win is recorded in PR #555. New work audits the gap between the
implemented paired-run t-interval CI and the original 20-pair/three-session
bootstrap protocol, including p95 and resource gates. Preserve both records;
do not relabel a shorter run as completion of the stronger protocol.

**T1.6 sealed cohort (2026-09-16; evaluated 2026-09-20):** the non-release annotated tag
`trust-eval-candidate-20260916` freezes `c2682e3aa75c21566c86ce1ade9cbd052838c694`
(tagged 2026-09-16T16:21:12+09:00). A ChEMBL source acquired after that freeze,
its maintainer attestation, and canonical/parent/scaffold audits against the
descriptor census, ChEMBL accuracy, and exposed browser-10k inputs are recorded
in `validation/results/sealed-cohort-preflight-trust-eval-candidate-20260916.json`.
The local-only raw input yields 2,000 development and 8,000 sealed holdout rows
from 10,239 eligible structures. The frozen candidate was evaluated once with
RDKit 2025.09.3: all 8,000 rows parsed, but six molecular-weight values were
unsupported and TPSA had 46 strict mismatches. The candidate is rejected and
this holdout is now exposed; the raw local-only result and a commit-safe summary
are recorded separately. The older 2026-09-13 preflight remains historical
`prepared_not_sealed` evidence.

**T0 result (2026-09-13):** independently built baseline/candidate arms completed
the fixed 5,021-molecule × 31-query lane with complete row accounting and zero
atom-alignment failures. The candidate regressed from **12 to 21 residual cells**
(with 10,042 pinned RDKit query-parse-error cells in each arm), so hybrid adoption
is blocked. The compact provenance record is
`validation/results/rdkit-smarts-baseline-candidate-v1.0.14.json`; the withdrawn
12 → 3 claim came from partial output and must not be reused. T0's measurement
packet is complete, while A4 residual classification and any replacement candidate
remain open.
The [Trust execution plan](docs/trust-release-plan.md) owns T0–T6 subtasks,
dependencies, resource budgets, and review dates.

## Open phase backlog

### Accuracy packages — A0–A6

These are remaining exits, not a list of all work already implemented.
See the [accuracy plan](docs/rdkit-accuracy-plan.md) for subtask acceptance and
the [disposition ledger](docs/roadmap-open-work.md) for evidence and dependencies.

- [x] **A0 — Evaluation contract:** candidate `5e9211a6` passed the frozen RDKit 2025.09.3 core-eight profile on 2,000/2,000 development and 8,000/8,000 one-time sealed rows, with full row accounting, pinned binary/evaluator/source/attestation hashes, and zero mismatches or unsupported values.
- [ ] **A1 — Perception and descriptors:** core-eight unused evaluation is complete through A0; finish potential-center evaluation and affected workflow/binding acceptance. Additional descriptor families retain separate gates.
- [ ] **A2 — Stereo and identity:** Issue #149's bounded shared-carrier residual is complete at 28/28 components × 1,024 relabelings with zero divergence and zero correspondence failure. Retain the aromatic-stash planner's scope, resolve phosphorus CIP adjudication, and meet the broader false-merge and independent-evaluation exits.
- [ ] **A3 — Fingerprints and retrieval:** retain raw/provenance evidence for the new aromaticity lane, rerun affected bindings/search after adoption, and pass unused-input evaluation. Existing k=1/10/100 and threshold evidence remains valid for its recorded builds.
- [ ] **A4 — Workflows and interchange:** finish T0 residual classification and adoption, then verify mapped SMARTS embeddings, standardization, reaction products, and typed V3000 semantics across engines.
- [ ] **A5 — Independent adjudication:** obtain absolute gold labels and unused inputs, secure non-maintainer review, and run the implemented paired evaluator against the frozen protocol.
- [ ] **A6 — 3D and force fields:** issue #337's six MMFF94 pyridinium/macrocycle typing and charge residuals are resolved on the pinned 265-molecule gate. Public v1.0.20 confirms that the production stereo-safe MMFF94 lane closes the prior 12 typed failures and four gross-clash rows: 265/265 independently sound, stereo-clean and clash-free, versus RDKit's 264/265. Its paired speed is 0.944x (95% lower 0.861x), so no quality-equivalent MMFF94 speed win is claimed. Public UFF best-of-10 reaches 265/265 usable versus RDKit's 264/265 at 2.359x (lower 2.198x). Resolve the two same-coordinate energy residuals, term-level adjudication, timeout and broader conformer-quality exits.

### Cross-cutting follow-ups

These extend the same product areas; they do not introduce new Phase numbers.

- [~] T3.6/T1.5: the v1.0.20 registry rerun closes the fixed-corpus npm/PyPI package step for prepared Morgan and bounded UFF. Complete multi-browser/host, resource, native-lane, residual-classification, and next-RDKit evidence. Keep native ECFP and T3.4 search separate.
- [ ] Make the BatchResult contract the cross-binding batch product boundary: known-length input accounting remains exact, and the public Rust/Python/WASM/Node stream schema distinguishes normal EOF from typed cancellation/time/resource/producer/consumer stops, observed-but-unprocessed rows, and an unread unknown suffix. Add retry, Worker/MCP adapters, export, clean-install, and resource proof without manufacturing success or skipped rows.
- [ ] Exact canonical-SMILES and cross-engine V3000 parity: expand semantic-identity and RDKit/Indigo fixtures beyond the completed T5.7 identity-renumber and T1.9 query-semantics gates; retain unsupported representation boundaries.
- [ ] Replace remaining MD/UFF/MMFF94 finite-difference production paths after same-domain analytic energy, gradient, and stereo soundness gates.
- [ ] Replace remaining periodic neighbor all-pairs paths after exact result-set and cutoff parity.
- [ ] Optimize symmetry-heavy canonical and SMARTS search after broader exact-output, invariance, and budget-exhaustion gates.
- [ ] Extend browser and agent adversarial coverage for cancellation, malformed records, limits, stable errors, and supported engines. T3.7 row accounting covers Rust/Python/Node/WASM; the merged Explorer unknown-length stream adapter remains browser-UI evidence only. Worker, cross-binding streaming, cancellation, export, and resource-budget exits remain.
- [ ] Expand reaction/SMARTS/medicinal-chemistry breadth after the shared P0–P3 primitives have current evidence.
- [ ] Add curated reaction/query precision, recall, invalid-product, timeout, and ambiguity reports with independent oracles.
- [ ] Close MMFF94/UFF typing, charge, parameter, convergence, and stereo gaps with same-coordinate comparisons.
- [ ] S5 independent gate: obtain non-maintainer review or external audit of parser, serialization, and binding boundaries.
- [ ] S6 continuous maintenance: rehearse advisory intake, fixes, backports, publication, and supported-version synchronization.

**Historical / abandoned:** the additional 1.10x SMILES speed stretch remains
abandoned and non-blocking. It is not an unchecked task and requires an explicit
maintainer decision to reactivate.

## v1.0.16 candidate boundary

This retained heading records the baseline planning boundary; v1.0.16 is already
released. The following requirements apply to the **next Trust RC**, not a
claim that publishing v1.0.16 completed them.

The next Trust RC requires **T0–T2**, T3's 10,000-row install/Worker gate,
T4's fixed safety corpus, and T5's public suite plus existing safety regressions.
It also retains **A0, A1 core-eight-field unused evaluation, affected binding
regressions, and native-default preservation** as mandatory exits.
The Week 4 milestone is an audit, not a promised release date; no new version
number is selected here.

The new upstream packet requires a large-ring round-trip/refusal regression and
explicit stereo/attachment capability accounting. Unsupported atrop or collapse
semantics need not become new APIs to ship, but must not silently lose information.
It also requires the T3.7 contract for the declared batch workflow and T5.7/T1.9
regressions (semantic preservation or explicit unsupported outcome). A historical
candidate's sealed result does not automatically qualify a later RC: record the
change impact and affected regression gates; a new unused-data claim needs a new
freeze/exposure audit. Future RDKit availability does not block the frozen run.
No performance claim is promoted without its matching correctness and coverage
gate; no remote-startup or unique-memory claim is required without measurement.

| Claim | Required acceptance |
|---|---|
| Declared 2D compatibility | All A0–A4 exits on the frozen support domain |
| Independent chemical equivalence or superiority | A5 absolute gold, review, and predeclared statistical criteria |
| 3D parity | A6's separately frozen numerical and conformer-quality gates |

Core descriptors require full valid-scope coverage and zero strict mismatches
(integer fields exact; floating fields within 1e-6). Retrieval requires recall
1.0, identical ordered IDs, and unrounded score error ≤1e-12. Stereo requires
zero wrong confident labels, information loss, and false merges.
Independent equivalence requires the whole paired 95% accuracy-difference
interval inside ±0.1 percentage point; superiority requires a positive lower
bound plus coverage safeguards. Full definitions live in the accuracy plan.

Release checks include affected tests, format/diff checks, versioned evidence,
and support/failure documentation. Candidate readiness and post-publication
registry/site verification are separate T2 checks.

## Product phases

The **Priority** column above expresses urgency. **Phase P0–P6** names product
areas, **A0–A6** names accuracy packages, and **T0–T6** names delivery work.

| Phase | Product area | Main active delivery |
|---|---|---|
| P0 | Trust and measurement | T0/T1/T2 |
| P1 | Interchange throughput and parser safety | T4; A4 |
| P2 | Identity, descriptors, fingerprints, and search | T0/T1/T5; A1–A3 |
| P3 | Rust/Python/Node/WASM and agent usability | T3/T4 |
| P4 | Chemistry workflows | T5; A4 |
| P5 | 3D and materials | T6; A6 |
| P6 | Ecosystem maintenance and external validation | T2/T6; A5 |

## Contracts and reference documents

The v1 boundaries remain: canonical SMILES is a representation; stable-key APIs
fail closed. Native ECFP4 and RDKit-compatible Morgan are separate profiles.
RWMol, CDXML, Markush/polymer, and interchange support have declared bounds.
3D/MMFF94 remains Experimental; pure-Rust InChI is approximate, while standard
InChI uses opt-in native FFI. See [compatibility scope](docs/compatibility-scope.md).

A package is complete when implementation, tests, documentation, and required
measurements agree across its declared domain. Safe refusal remains visible in
coverage. A superiority claim must identify its operation, versions, corpus,
configuration, hardware, failures, uncertainty, and reproduction command.

| Document | Responsibility |
|---|---|
| [Trust Release execution plan](docs/trust-release-plan.md) | T subtasks, dependencies, budgets, and schedule |
| [Parse + Morgan performance plan](docs/parse-morgan-performance-plan.md) | T3.6 / PF0–PF5, equivalent-operation speed targets and correctness gates |
| [RDKit accuracy plan](docs/rdkit-accuracy-plan.md) | A subtasks, evaluation protocol, and exits |
| [Open-work disposition](docs/roadmap-open-work.md) | Evidence ledger and local/toolchain/external dependencies |
| [Validation guide](docs/validation.md) / [benchmark index](benchmarks/README.md) | Reproduction and measurement boundaries |
| [2026-09-13 snapshot](docs/archive/roadmap-through-2026-09-13.md) | Long-form status before this reorganization |
| [Archive through 2026-09-12](docs/archive/roadmap-through-2026-09-12.md) | Earlier implementation history |
