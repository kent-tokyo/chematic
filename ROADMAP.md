# chematic roadmap

> Updated 2026-09-16. Released: **v1.0.15**. Next delivery theme: **1.x Trust Release**.
> **Current P0 focus: freeze T1's prepared 10,000-row split with a candidate tag and
> unused-data attestation; keep the published RDKit.js comparison reproducible; then
> advance parser/stereo and V3000 contracts without widening unsupported semantics.**

## Current status

- v1.0.14 was published from `a259507d`; the release and CI follow-ups reached
  `main` in [PR #532](https://github.com/kent-tokyo/chematic/pull/532)
  (`286485be`). The tag retains the published source; the merge includes later
  formatting, native-InChI regression-test, and binding-inventory corrections.
- v1.0.15 was published from `afa61956`. Its six public channels, Pages asset
  digests, clean-install consumers, modern MCP smoke, and binary-only PyPI
  runtime smoke are recorded in
  `validation/results/release-channel-verification-v1.0.15.json` and
  `validation/results/release-package-smoke-v1.0.15-2026-09-16.json`. The
  latter verifies Linux CPython 3.9, macOS CPython 3.13, and Windows CPython
  3.13; Linux 3.9 is the actual published-wheel boundary for this release.
- The published-package browser scorecard is now in `main` via
  [PR #541](https://github.com/kent-tokyo/chematic/pull/541) (`8ce8026`). It
  compares `@kent-tokyo/chematic@1.0.15` with official
  `@rdkit/rdkit@2026.03.6` on a fixed exposed 10k corpus. It records a smaller
  WASM artifact, lower local no-store ready and parse/write timings, a
  RDKit-faster parse-inclusive fingerprint lane, one typed Fe(II)
  coordination refusal, and every unmeasured network/unique-memory dimension.
  The local-ready measurement is not a CDN or cross-host claim.
- [PR #542](https://github.com/kent-tokyo/chematic/pull/542) (`3d61132`)
  updates the locked `rustls` dependency to 0.23.45. Its hosted Security Audit
  passed the local review, isolated parser-security corpus, and Cargo Audit;
  the prior RUSTSEC-2026-0285 finding is resolved for the merged lockfile.
- All **seven accuracy packages A0–A6 still have open acceptance work**.
  Publication does not complete the planned Trust RC or independent accuracy gates.
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
| Official RDKit npm/WASM reduces installation friction | Compare pinned public artifacts on the same browser operation; publish Worker/cancellation/error contracts and losses as well as wins | Claim broad RDKit replacement or hide the RDKit-faster fingerprint lane |
| RDKit CIP and SMILES fixes expose order-sensitive chemistry risk | Treat atom-order and SMILES-spelling permutations as release regressions; retain uncertain cases as typed refusals | Treat agreement for one spelling as stereochemical proof |
| Indigo advances V3000, coordination, and polymers | Test actual Indigo → chematic → RDKit and reverse round trips; separate opaque retention from editable typed semantics | Start a haptic/polymer feature race before interchange semantics are proven |
| CDK relaxed parsing helps interactive editors | Keep normal parsing strict; any preview mode must be a separately named incomplete-state API | Silently accept incomplete chemistry in normal APIs |
| Open Babel parser hardening raises the security bar | Publish fixed malformed-input corpus accounting with time/memory limits and refusal reasons | Treat Rust's memory model alone as parser-security proof |
| Other Pure-Rust parity libraries compete for the same headline | Differentiate through reproducible cross-binding contracts, browser Worker UX, and release evidence | Chase broad 3D, protein, polymer, or format breadth without a validated workflow |

### Competitive deliverables

| Order | Deliverable | Acceptance evidence |
|---|---|---|
| 1 | **T1.6 sealed-evaluation precondition** | Frozen candidate and a legitimate unused-data attestation before consuming the prepared 2k/8k split |
| 2 | **T1.5/T3.4 compatibility and cost scorecard** | Separate RDKit oracle lanes; public-artifact bytes, local-ready, operation timings, coverage/refusals, and all missing dimensions together |
| 3 | **T4/T5 safety and stereo gates** | Hosted security audit, bounded malformed corpus, permutation/spelling suites, and typed refusal accounting |
| 4 | **T1.7 interchange contract** | Version-pinned RDKit/Indigo ordinary-V3000 round trips, SGROUP edits, COLLECTION ordering, and explicit coordination/ENDPTS exclusions |
| 5 | **T3 adoption and version truth** | `npm install` to typed Worker/free, failed-row-preserving 10k Explorer flow, and one current-version source across public surfaces |

## Execution order

| Order | Priority | Delivery and next output | Product Phase / accuracy |
|---|---|---|---|
| 1 | P0 | **T1.6 sealed-evaluation precondition**: freeze the candidate build and obtain legitimate unused-data attestation before consuming the prepared 2k/8k split. The exposed SMARTS packet is not sealed evidence | P0/P2; A0–A4 |
| 2 | P0 | **T1.5/T3.4 public compatibility-cost contract**: retain the merged v1.0.15 browser scorecard and add only equivalent, pinned RDKit lanes; do not relabel local no-store timing as internet startup | P0/P3; A3 |
| 3 | P1 | **T4/T5 safety and stereo gates**: retain hosted parser and Cargo Audit success; expand bounded malformed-input, permutation, round-trip, and refusal accounting | P1/P2/P4; A2/A4/A5 |
| 4 | P1 | **T1.7 interchange contract**: expand actual Indigo typed-SGROUP creation/editing and coordination/ENDPTS boundaries without treating opaque retention as editable semantics | P1/P4; A4 |
| 5 | P2 | **T3.1–T3.6/T2.6 adoption and version truth**: published-package install examples, diverse 10k Explorer workflow, current-version synchronization, bounded pilot evidence, independent review, and existing 3D gaps | P3/P5/P6; A5/A6 |

T3/T4/T5 can proceed in parallel once their dependencies pass. Public-channel
inventory, unused-data acquisition, and A5 gold preparation start alongside T0.
Confirmed silent corruption, panic, or resource-limit defects take precedence.

**T1.6 preflight (2026-09-13):** an independently acquired ChEMBL candidate
has 11,359 non-overlapping eligible rows after canonical/parent/scaffold audit,
which deterministically yields 2,000 development and 8,000 holdout candidates.
It is explicitly `prepared_not_sealed`: no candidate commit/tag or unused-data
attestation exists yet, so it must not be used for model selection or presented
as a sealed evaluation.

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

- [ ] A0 — Accuracy evidence contract: freeze genuinely unused 8,000-row evaluation inputs and candidate builds; complete the full acceptance packet. The 96/96 all-field holdout and exposed-corpus validators are narrower completed slices.
- [ ] A1 — Perception and descriptor coverage: complete unused core-eight-field and potential-center evaluation, plus affected workflow/binding acceptance. The opt-in profile has passing exposed-corpus evidence; additional families have separate gates.
- [ ] A2 — Stereo and identity correctness: resolve the four #503 canonical/E/Z components and phosphorus CIP adjudication; meet permutation, idempotency, information-preservation, and false-merge exits.
- [ ] A3 — Fingerprint and retrieval fidelity: retain raw/provenance evidence for the new aromaticity lane, rerun affected bindings/search after adoption, and pass unused-input evaluation. Existing k=1/10/100 and threshold evidence remains valid for its recorded builds.
- [ ] A4 — Workflow accuracy: finish T0 residual classification and adoption, then verify mapped SMARTS embeddings, standardization, reaction products, and typed V3000 semantics across engines.
- [ ] A5 — Independent accuracy adjudication: obtain absolute gold labels and unused inputs, secure non-maintainer review, and run the implemented paired evaluator against the frozen protocol.
- [ ] A6 — 3D accuracy extension: resolve remaining MMFF94/UFF term, charge, gradient, convergence, timeout, stereo, and conformer-quality gaps. Bounded typing or energy matches do not complete this package.

### Cross-cutting follow-ups

These extend the same product areas; they do not introduce new Phase numbers.

- [ ] Measure only equivalent operations against fixed RDKit/Open Babel versions on identical inputs; keep subprocess-only lanes separately labeled.
- [ ] Exact canonical-SMILES and cross-engine V3000 parity: expand semantic-identity and RDKit/Indigo fixtures, including unsupported representation boundaries.
- [ ] Replace remaining MD/UFF/MMFF94 finite-difference production paths after same-domain analytic energy, gradient, and stereo soundness gates.
- [ ] Replace remaining periodic neighbor all-pairs paths after exact result-set and cutoff parity.
- [ ] Optimize symmetry-heavy canonical and SMARTS search after broader exact-output, invariance, and budget-exhaustion gates.
- [ ] Extend browser and agent adversarial coverage for cancellation, malformed records, limits, stable errors, and supported engines.
- [ ] Expand reaction/SMARTS/medicinal-chemistry breadth after the shared P0–P3 primitives have current evidence.
- [ ] Add curated reaction/query precision, recall, invalid-product, timeout, and ambiguity reports with independent oracles.
- [ ] Close MMFF94/UFF typing, charge, parameter, convergence, and stereo gaps with same-coordinate comparisons.
- [ ] S5 independent gate: obtain non-maintainer review or external audit of parser, serialization, and binding boundaries.
- [ ] S6 continuous maintenance: rehearse advisory intake, fixes, backports, publication, and supported-version synchronization.

**Historical / abandoned:** the additional 1.10x SMILES speed stretch remains
abandoned and non-blocking. It is not an unchecked task and requires an explicit
maintainer decision to reactivate.

## Next candidate boundary

The next Trust RC requires **T0–T2**, T3's 10,000-row install/Worker gate,
T4's fixed safety corpus, and T5's public suite plus existing safety regressions.
It also retains **A0, A1 core-eight-field unused evaluation, affected binding
regressions, and native-default preservation** as mandatory exits.
The Week 4 milestone is an audit, not a promised release date; no new version
number is selected here.

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

## Product Phase numbers

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
| [RDKit accuracy plan](docs/rdkit-accuracy-plan.md) | A subtasks, evaluation protocol, and exits |
| [Open-work disposition](docs/roadmap-open-work.md) | Evidence ledger and local/toolchain/external dependencies |
| [Validation guide](docs/validation.md) / [benchmark index](benchmarks/README.md) | Reproduction and measurement boundaries |
| [2026-09-13 snapshot](docs/archive/roadmap-through-2026-09-13.md) | Long-form status before this reorganization |
| [Archive through 2026-09-12](docs/archive/roadmap-through-2026-09-12.md) | Earlier implementation history |
