# chematic roadmap

> Updated 2026-09-20. Released: **v1.0.17**. Next delivery theme: **1.x Trust Release**.
> **Current P0 focus: close the A0 acceptance packet, then evaluate the frozen
> T1.6 accuracy candidate once.** Sealed scores remain unrun.
> The merged Parse + Morgan browser win moves to maintenance and package verification.

## Current status

- v1.0.17 adds bounded SMILES+ large-ring labels and typed attachment-label
  inspection on the Trust/interoperability baseline. Release history belongs
  in `CHANGELOG.md`; publication does not complete the sealed accuracy gates.
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
- All **seven accuracy packages A0–A6 still have open acceptance work**.
  Publication does not complete the planned Trust RC or independent accuracy gates.
- PR #555 (`7d98dcd3`) merges the Parse + compatible Morgan optimization with
  supported-domain bit agreement and local/hosted browser speed evidence.
  It is a different candidate from the September 16 accuracy freeze. Published
  package timing, the original multi-session protocol, and broader runtime
  resource acceptance still require their own evidence.
- This index tracks **18 open areas**: seven accuracy packages and eleven
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
| RDKit npm/WASM and Python typing improve usability | Compare pinned public artifacts; publish Worker/cancellation/error contracts and version-specific wins and losses | Infer runtime safety from type annotations or generalize one browser speed win |
| RDKit CIP and SMILES fixes expose order-sensitive chemistry risk | Treat atom-order and SMILES-spelling permutations as release regressions; retain uncertain cases as typed refusals | Treat agreement for one spelling as stereochemical proof |
| Indigo advances V3000, coordination, and polymers | Test actual Indigo → chematic → RDKit and reverse round trips; separate opaque retention from editable typed semantics | Start a haptic/polymer feature race before interchange semantics are proven |
| CDK relaxed parsing helps interactive editors | Keep normal parsing strict; any preview mode must be a separately named incomplete-state API | Silently accept incomplete chemistry in normal APIs |
| Open Babel parser hardening raises the security bar | Publish fixed malformed-input corpus accounting with time/memory limits and refusal reasons | Treat Rust's memory model alone as parser-security proof |
| COSMolKit and other Pure-Rust parity libraries compete for the same headline | Preserve every input's index, terminal outcome and error stage across bindings; publish reproducible runtime and release evidence | Infer complete batch success from zero export failures |

The September 20 review turns upstream changes into scoped regression work,
not an automatic promise of new chemistry support. The
[current execution packet](docs/trust-release-plan.md#7-2026-09-20-trust完了に向けた実行順)
records primary sources, implementation gaps, dependencies, and acceptance.
Existing published 1.0.15 measurements stay historical; new release/candidate
measurements must identify their actual package/source hashes.

## Priority order

| Order | Priority | Delivery and next output | Product Phase / accuracy |
|---|---|---|---|
| 1 | P0 | **A0 recovery after frozen 8k rejection**: preserve the failed candidate's evidence, classify its declared-scope gaps without reusing the exposed holdout, and prepare a new freeze/unused cohort for a later adoption decision | P0/P2; A0–A4 |
| 2 | P0 preparation; external execution | **T1.5 rebaseline**: prepare paired old/new RDKit lanes now; run when the actual stable Python/native/npm artifacts exist. Record wrapper backend and API behavior, not just version strings | P0/P2/P3; A1–A4 |
| 3 | P1 | **A2 closure**: four #503 canonical/E/Z components and phosphorus adjudication; retain T5.6 selective CIP/isotope/atrop and T1.8 attachment boundaries. T5.7/T1.9 and T4.6 binding regressions are implemented development gates, not a broad interoperability claim | P1/P2/P4; A2/A5 |
| 4 | P1 | **T3.1–T3.4/T2 runtime and release evidence**: build on the completed T3.7 row-accounting contract with 10k Worker/cancel/offline/limits, typed install examples, current package scorecard and dashboard. T3.6 maintains the scoped speed win and audits remaining protocol conditions | P0/P3/P6; A3/A4 |
| 5 | P2; P0 for confirmed incorrectness | **T6 / A5–A6**: independent gold/review and existing 3D typing, charge, gradient, convergence, timeout and stereo gaps | P5/P6; A5/A6 |

T3/T4/T5 can proceed in parallel once their dependencies pass. Public-channel
inventory and A5 gold preparation continue alongside frozen-packet verification.
Confirmed silent corruption, panic, or resource-limit defects take precedence.

### Parse + ECFP4/Morgan performance workstream

The [detailed T3.6 plan](docs/parse-morgan-performance-plan.md) owns the contract,
profiling hypotheses, six steps **PF0 → PF1 → PF2 → PF3 → PF4 → PF5**, and adoption gates.
PF IDs are subtasks, not new P7+ product phases. This is a new explicit performance
goal; the abandoned additional 1.10x SMILES stretch stays abandoned.

Historical public v1.0.15 Chromium Parse + FP took **0.315347 vs 0.134228 ms/mol**
(CheMatic vs RDKit; medians of process means). That is about **2.35x the time**,
not a current v1.0.17 measurement. The adopted runner records its equivalent
packed output separately from the historical operation boundary.

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

- [ ] **A0 — Evaluation contract:** the September 16 frozen candidate was evaluated and rejected; prepare a new candidate freeze, unused cohort, and complete acceptance packet for any later adoption decision. The 96/96 all-field holdout and exposed-corpus validators are narrower completed slices.
- [ ] **A1 — Perception and descriptors:** complete unused core-eight-field and potential-center evaluation, plus affected workflow/binding acceptance. The opt-in profile has passing exposed-corpus evidence; additional families have separate gates.
- [ ] **A2 — Stereo and identity:** resolve the four #503 canonical/E/Z components and phosphorus CIP adjudication; meet permutation, idempotency, information-preservation, and false-merge exits.
- [ ] **A3 — Fingerprints and retrieval:** retain raw/provenance evidence for the new aromaticity lane, rerun affected bindings/search after adoption, and pass unused-input evaluation. Existing k=1/10/100 and threshold evidence remains valid for its recorded builds.
- [ ] **A4 — Workflows and interchange:** finish T0 residual classification and adoption, then verify mapped SMARTS embeddings, standardization, reaction products, and typed V3000 semantics across engines.
- [ ] **A5 — Independent adjudication:** obtain absolute gold labels and unused inputs, secure non-maintainer review, and run the implemented paired evaluator against the frozen protocol.
- [ ] **A6 — 3D and force fields:** resolve remaining MMFF94/UFF term, charge, gradient, convergence, timeout, stereo, and conformer-quality gaps. Bounded typing or energy matches do not complete this package.

### Cross-cutting follow-ups

These extend the same product areas; they do not introduce new Phase numbers.

- [ ] T3.6/T1.5: maintain the merged browser Parse + Morgan win, audit remaining statistical/resource conditions, verify the published package and rebaseline the next available RDKit artifacts; retain separate native ECFP and T3.4 search lanes.
- [ ] Exact canonical-SMILES and cross-engine V3000 parity: expand semantic-identity and RDKit/Indigo fixtures, including T5.7 identity-renumber and T1.9 query-semantics regressions and unsupported representation boundaries.
- [ ] Replace remaining MD/UFF/MMFF94 finite-difference production paths after same-domain analytic energy, gradient, and stereo soundness gates.
- [ ] Replace remaining periodic neighbor all-pairs paths after exact result-set and cutoff parity.
- [ ] Optimize symmetry-heavy canonical and SMARTS search after broader exact-output, invariance, and budget-exhaustion gates.
- [ ] Extend browser and agent adversarial coverage for cancellation, malformed records, limits, stable errors, and supported engines. T3.7 row accounting now covers Rust/Python/Node/WASM; Worker, streaming, cancellation, export and resource-budget exits remain.
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
