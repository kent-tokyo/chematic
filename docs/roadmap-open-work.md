# Roadmap open-work disposition

This document is the current disposition of unchecked roadmap items for the
v1.0.20 candidate follow-up (updated 2026-09-22). An item is not complete merely
because a narrower local slice has evidence. Recorded v1.0.13 and earlier
measurements keep their original scope; they were not rerun for this reorganization.

2026-09-16 release-proof note: `v1.0.15` now has a complete six-channel record
and binary-only PyPI runtime artifacts from GitHub Actions run `35061960770`.
The published wheel imported and returned benzene formula `C6H6` on Linux
CPython 3.9, macOS CPython 3.13, and Windows CPython 3.13. Linux 3.9 is the
actual v1.0.15 wheel boundary; this is not a claim for newer Linux Python
versions. See `validation/results/release-package-smoke-v1.0.15-2026-09-16.json`.

The [active roadmap](../ROADMAP.md) owns execution order and the 18-area open
backlog. This ledger owns evidence and dependencies, not a second priority list.

Dependency classes used below: `local-open`, `toolchain-open`, `data-sealed`,
`external-open`, and `historical`.
The abandoned additional 1.10x SMILES speed target is historical and excluded
from the open count. The new T3.6 Parse + Morgan target is part of the existing
cross-cutting comparison area, not an additional product phase or accuracy package.
See the [pre-reorganization snapshot](archive/roadmap-through-2026-09-13.md)
for the previous long-form roadmap.

## Current delta over the historical evidence below

The long A0–A6 measurement rows below retain their September 13 provenance.
Their pending-state wording is superseded by these narrower completion records;
none completes an entire A package:

| Area | Completed slice | Remaining disposition |
|---|---|---|
| CI reliability / issue #70 | Process-level ABBA observations, reversed-order Stage 2, practical-effect gating, and an independent-build null control are now fully calibrated. Hosted fixtures block +10%, detect +5% in 8/10 independent runs, and downgrade material one-sided contamination to `environment-inconclusive`. The null uses ten blocks plus a 1.04 magnitude floor because seven unanimous blocks cannot cross the configured Wilson/0.2 boundary. | `complete`: the machine-checked contract is `validation/criterion-gate-calibration.json`; hosted evidence is `validation/results/criterion-gate-hosted-calibration-2026-09-22.json`. Temporary injection commits are evidence-only and must never be merged. |
| T0 / A4 | Independent baseline/candidate packet; 12 → 21 residuals, candidate rejected | `local-open`: classify residuals and evaluate a replacement; do not repeat the withdrawn 12 → 3 claim |
| T1.6 / A0 | Historical rejected and TPSA-only candidates remain preserved. Candidate `5e9211a6`, frozen as annotated tag `trust-eval-candidate-a0-multifield-r3-20260922` before a third source acquisition, uses the bounded RDKit-parity descriptor aromaticity path and preserves authoritative zero-aromatic results. The new source was canonical/parent/scaffold-audited against every exposed source. Development passed 2,000/2,000 and the one-time sealed holdout passed all eight declared fields at 8,000/8,000, with zero mismatches, parse failures, or unsupported values. | `complete`: `validation/results/a0-core-eight-sealed-acceptance-20260922.json` is validated by `scripts/check_a0_core_eight_sealed.py`. Raw input/output remain local-only and exposed; A0 completion does not complete A1–A6. |
| T1.7 / A4 | Ordinary V3000, external RDKit/Indigo readers and SGROUP/COLLECTION baseline (PR #544) | `local-open`: typed chemistry/edit semantics and declared refusal boundaries; general CX label retention does not complete T1.8 attachment identity/collapse |
| T1.8 / A4 | September 19: complete positive `_AP<n>` syntax parser and degree-one-wildcard identity API, with invalid/non-dummy/multi-degree rejection. September 21: RDKit 2025.09.3 gate requires 8/8 RDKit → CheMatic CLI → RDKit label/map/topology and degree-one-wildcard identity results; canonical mapped wildcard output retains `[*:n]`. | `local-open`: bond/stereo-bearing CX cases and any collapse semantics; identity is deliberately not collapse eligibility |
| T3.4 / A3/A4 | Public 1.0.15 vs official RDKit npm 10k, three browsers, 20 repetitions; separate three-run RSS diagnostic | `local-open` / `toolchain-open`: public 1.0.17 and candidate lanes, equivalent output/search, resource definitions; unavailable remote host/measurement stays explicit |
| T3.6 / A3 | Registry-installed v1.0.20 is 9,999/9,999 exact on the fixed supported corpus and passes parse-inclusive/prepared Morgan at 1.398x/3.511x geometric speedup (95% lower bounds 1.363x/3.407x). One high-coordinate Fe(II) row remains a typed refusal. | `local-open` / `external-open`: the public-package step is complete; finish the original multi-session/bootstrap, p95/resource and multi-browser/host conditions. Preserve the separate native ECFP profile. |
| T4.6 / A2/A4 | September 19 core fix plus PR #562: parser/writer use bounded `%(n)` labels ≥100; malformed/overflow legacy ambiguity rejects; a 121-closure graph round-trips and shared Rust/Python/Node/WASM validity fixtures accept `%(100)` while rejecting malformed legacy spellings | `local-open`: canonical-writer atom-order and resource-limit coverage remain; this is parser validity coverage, not external-engine round-trip evidence |
| T5.6 / A2 | Upstream CIP/atrop changes identified; existing development suite is not their acceptance | `local-open`: selection/isotope/atrop/permutation regressions; independent gold is `external-open` |
| T1.5 next-oracle rebaseline | The historical 2026.03.6 Python/npm packet is runnable and CI-checked: exact wheel/tarball hashes, Boost.Python behavior, toolchain/configuration, and 10,000 complete rows are committed. CheMatic 1.0.19 records Morgan 9,999/10,000, correspondence-correct CIP 9,770/10,000, SMILES semantic 9,982/10,000, and 14,306/310,000 differing SMARTS cells. The older 9,880 CIP figure is withdrawn because it compared engine-local bond indices. All 18 SMILES semantic differences are classified as CheMatic stereo-writer regressions with graph identity preserved and tracked by #632; the single Morgan difference is a typed unsupported coordination contract. The official RDKit GitHub release record has no attached native/C++ artifact, so that independent lane is explicitly unavailable rather than inferred from Python/WASM | T1.5 preparation is complete after CI validation; product-level adjudication continues in #634 (230 CIP rows) and #635 (3,364 SMARTS rows) without changing the frozen T1.6 oracle. `external-open`: execute the same old/new lanes only after the next released Python/native/npm artifacts are independently verified |
| T5.7 + T1.9 / A2/A4 | PR #558: atom-map potential-stereocenter invariance covers identity/reordered atoms. The follow-up T5.7 gate now exercises 32 fixed graph-storage orders through SSSR, potential-center/CIP, clone, and SMILES reparse while comparing atom-map center/CIP identities and bond-map endpoints. Its ordinary saturated tertiary-amine negative control applies the same 32-order/call-order sequence and remains at zero centers. V3000-to-V3000 preserves opaque query attributes and unmodelled query constraints refuse conversion to ordinary formats. The current 4×4 E/Z query observation gate covers symmetric 2-butene and mixed-substituent chloro/bromo alkene pairs, and preserves RDKit isomeric identity plus the chiral `HasSubstructMatch` truth table across the bounded round trip. Indigo's observed source truth table is also preserved, but Indigo 1.46.0 matches both E/Z targets for all four V3000 queries; that external semantic loss is recorded, not counted as query compatibility. | `local-open`: extend E/Z query truth tables across versions/engines and add typed semantics before compatibility claims. Unsupported paths remain refusal boundaries, not query compatibility. See `validation/results/v3000-indigo-ez-query-truth-table-current-2026-09-21.json`. |
| T3.7 / P3 | PR #559: shared Rust/Python/Node/WASM fixtures enforce `input_count = success + failed + refused + skipped`, original row indices, terminal status, parse error stage, and `all_succeeded()`. PR #584 is merged (`5678360b`): its Explorer-only `File.stream()` adapter records completed rows, observed-but-unprocessed rows, `unread input=unknown`, and a terminal reason on cancellation; the 10k Chromium smoke exercises normal CSV export/order and that cancel boundary. The public stream contract adds a shared fixture plus Rust/Python/WASM/Node `canonicalize_smiles_stream` schema: normal EOF processes every observed row, while typed stop reasons preserve the processed prefix and make the unread suffix explicitly unknown. | `local-open`: clean-install the public schema, then add retry, Worker/MCP adapters, migration examples, and measured 10k resource limits. The Explorer adapter and a core public schema are not Worker/MCP or resource-budget evidence. |

The [Trust execution packet](trust-release-plan.md#7-2026-09-21-trust完了に向けた実行順)
owns the new subtask exits. These are subdivisions of the existing 18 open areas,
not additional phases or a reason to mark A0–A6 complete.
The [performance plan](parse-morgan-performance-plan.md) owns T3.6 / PF0–PF5.
Historical 1.0.15 timings are not evidence of 1.0.17 performance; speed acceptance
requires exact supported-domain bits, preserved coverage, and paired uncertainty.
The sealed accuracy cohort must not be used for performance tuning.

## Accuracy follow-up (A0–A6)

Six packages retain open exit work; A0 core-eight is complete. The 2026-09-12 plan audit preserves
passing local results but reopens A0 and A3 against their full original scope:
the initial 12-row holdout checked MW only, but the subsequent all-field run
passes 96/96. Neither run is the planned 8,000-row sealed evaluation.
The latest A3.4 rerun supersedes the older
search description: the same 4,500-library/500-query split was checked at
k=1/10/100 using the precise WASM JSON path, with 1,500/1,500
Rust/Python/Node-WASM comparisons matching. Neither this nor the holdout proves
the broader declared contract.
P0–P6 area numbers and historical results below remain unchanged.

See the [detailed accuracy plan](rdkit-accuracy-plan.md) for stable task IDs,
API/comparator contracts, fixed denominators, unused-data splits, numerical
targets, independent adjudication, and candidate boundaries.

2026-09-13 execution note: the shared SymmSSSR path now receives the caller's
candidate cap and reports measured candidate counts; the SMARTS dump producer
emits parse-failure rows and a required completion footer; and the generated
compatibility dashboard separates RDKit compatibility from binding consistency.
These are local T0/T1 improvements, not sealed evaluation or independent review.

The 2026-09-13 A6 macrocycle-boundary experiment is explicitly rejected:
generalizing the two-large-ring condition to four large rings introduced four
new type mismatches in each of `chembl_tier_b_0029` and `_0030`. The production
condition was restored and the full 230-test `chematic-ff` suite passed.
Evidence: `validation/results/mmff94-macrocycle-boundary-experiment-a6-2026-09-13.json`.

> **2026-09-22 disposition override:** A0 is complete on the core-eight frozen
> profile. The retained long-form A0 row below is historical context; the
> authoritative current evidence is
> `validation/results/a0-core-eight-sealed-acceptance-20260922.json` and its
> fail-closed checker. A1–A6 remain independently open.

| Work package | Remaining work | Dependency class and current evidence |
|---|---|---|
| A0 / P0 | A0.1–A0.5: all-field strict holdout, raw-result accounting, pinned source/binary provenance, adversarial gate tests, separately built baseline/candidate | `local-open`. Source-built Python with RDKit 2025.09.3 passes all 40,000 field checks strictly; the 12-row all-field holdout passes 96/96. The 5,000-row ECFP4 binding gate passes 5,000/5,000 on all three bindings and records Python extension and Node/WASM hashes. The 8-field descriptor binding gate now passes 5,000/5,000 and 7,737/7,737 exposed rows across Rust/Python/Node/WASM; this is binding parity, not an RDKit oracle. Negative coverage rejects bool counts, empty successes, extra fields, missing provenance, duplicate IDs, and non-finite values. On main, the RDKit descriptor profile now looks up every one of the 3,111 isotope masses recognized by the pinned public periodic table; generated-table lookup and a source-built Python carbon-11 check are development evidence only. The September 20 candidate had fresh provenance, an annotated tag, unused-data attestation, and 8,000/8,000 parsed rows, but was rejected at 7,998/8,000 strict molecular-weight agreement; its input/result are exposed/local-only. A later adoption decision requires a new freeze, genuinely unused cohort, provenance/overlap checks, and full scorecard automation. |
| A1 / P1/P2 | A1.1–A1.4: common perception and atom-level center FP/FN; A1.5: additional descriptor families adopted separately | `local-open`. The eight-descriptor compatibility profile is now `adopted_opt_in`: source-built core fields pass 5,000/5,000 strictly, the all-field holdout passes 96/96, and native-default preservation passes. A provenance-preserving carbon-branch fallback plus a narrow RDKit-compatible ring-constrained sp3-N rule yields potential-center parity 5,000/5,000; the expanded 7,737-row exposed holdout also passes 7,737/7,737 with atom-level FP=0/FN=0. The bounded chordless-aromatic-cycle candidate now passes all eight descriptor fields 7,737/7,737 strictly against pinned RDKit 2025.09.3, with promotion decision `adopted_opt_in` and native defaults preserved. A1.1 uses physical nuclide masses for known isotope labels and RDKit's mass-number fallback for unlisted valid labels; this behavior is independently covered by the descriptor unit gate. Morphine/codeine aromatic-oxide LogP residuals are resolved while the 6-member cyclic-ether negative remains protected. After the P-H phosphonate, thioamide-N, and charged-sulfur fixes, fresh descriptor and ChEMBL rechecks each pass all eight fields 5,000/5,000 strict; the current-source 7,737-row exposed holdout is also strict-exact for all eight fields with parse/unsupported failures 0. Evidence: `validation/results/descriptor-rdkit-a1.1-tpsa-fixed-v1.0.13.json` and `validation/results/descriptor-rdkit-a1.4-recheck-v1.0.13.json`. The latest 8-field Rust/Python/Node-WASM binding regression also passes 7,737/7,737; affected local perception/chem/fingerprint/WASM/SMARTS/3D suites pass 204/864/309/344/185/597. These are binding and regression consistency results, not independent chemical correctness. Potential-center results above predate this shared-perception change (fixes3). The center runner records five explicit RDKit atom-set fallbacks and zero unresolved oracle rows. The implementation excludes ordinary tertiary amines, fused degree-3 junctions, amide/sulfonamide environments, and aromatic N. Unused/sealed evaluation, workflow-level candidate acceptance, broader A1.1 per-atom corpus, and additional descriptor families remain. Count equality alone does not prove center-set equality |
| A2 / P2 | A2.1–A2.5: phosphorus CIP, independent label checks, current K=1,024, semantic preservation, permutations, idempotency and key collisions | `local-open`, with independent adjudication in A5. Issue #149's measured shared-carrier residual is complete: the clean-main K=1,024 audit passes all 28 coupled components with zero divergence and zero cross-correspondence failure, and all four historical residuals have explicit equivalent-spelling regressions. Historical 4/28 evidence and the rejected 7/28 close-side experiment remain preserved. `scripts/check_ez_residual_evidence.py` validates both histories and the current 0/28 result. Stable-key admission remains limited to the semantically reparsed aromatic-stash path; unmeasured shapes fail closed. Existing phosphorus remains fail-closed. The 181-row label audit reports atom-map agreement 155/181, E/Z bond-map agreement 14/14, all 26 atom mismatches confined to P + `oracle_unstable`, and no assigned-label mismatch. Standardize→canonicalize idempotency passes 14,999/14,999 exposed rows. Remaining A2 work is phosphorus adjudication, independent labels, false-merge analysis, and unused evaluation. Evidence: `validation/results/ez_shared_carrier_coupling_mechanism_audit_summary_1024_2026-09-22.json`, `validation/results/cip-residual-stereocenter-v1.0.13.json`, `validation/results/cip-label-parity-a2-v1.0.13.json`, and `validation/results/canonical-idempotency-a2.4-current-v1.0.13.json`. |
| A3 / P2/P3 | A3.1–A3.5: raw FP oracle, k=1/10/100 and threshold edges, unrounded scores/order, shared input IDs, independent RDKit plus every binding, unused evaluation | `local-open`. Rust/source-built Python/Node-WASM agree at k=1/10/100 for 500 queries × 4,500 rows with 0 pairwise mismatches and unrounded f64 scores; an independent exhaustive RDKit oracle also passes all k. The shared fixture has 34 success cases plus one explicit unsupported-bond error case, passing Rust/Python/WASM. The inclusive threshold gate passes 96 RDKit-derived boundary cases across all three bindings; finite `[0, 1]` validation and one-ulp equality handling are explicit. A separate RDKit-parity aromaticity lane resolves all 59 Hueckel radius-1 raw Morgan residuals on the checked-in 5,000-row corpus: 5,000/5,000 exact against pinned RDKit 2025.09.3, with zero parity preprocessing errors. Evidence: `validation/results/ecfp-rdkit-raw-identifier-parity-aromaticity-variant-2025.09.3-5000.json`. This is candidate evidence, not yet an adopted default: raw trace/provenance retention, post-adoption all-binding/search reruns, and genuinely unused-input evaluation remain open |
| A4 / P1/P4 | A4.1–A4.5: standardization semantics, SMARTS mapped sets, reaction products and typed V3000 metadata in both directions | `local-open`; uninstalled Indigo is `local-toolchain`, data permissions or external fixtures may be `external`. Footer-verified persistent-corpus comparison over 5,021 molecules × 31 queries covers 155,651 cells: parity target-atom-set matches 145,588, RDKit SMARTS parse errors 10,042, residuals 21, and atom-alignment failures 0. Evidence: `validation/results/rdkit-smarts-direct-chembl-5000-footer-verified-v1.0.13.json`, including input/dump SHA-256 and comparison configuration; `scripts/check_rdkit_smarts_evidence.py` validates the source corpus hash, footer, row accounting, alignment, and bucket sum. Earlier private-corpus numbers are not mixed into this denominator; the 12→3 claim remains withdrawn. Matching compares target atom sets, not all mapped embeddings, and is not a frozen-build adoption packet. The shared path now propagates the caller budget and reports measured candidate counts, but independent baseline/candidate acceptance remains open. Expanded 51-case/14-template atom-map-aware canonical product-set gate passes 51/51; the opposite-enantiomer case remains a separate 1/1 fail-closed safety gate. V3000 SGROUP/COLLECTION, isotope, enhanced-stereo boundary, and opaque `ENDPTS=`/`ATTACH=` bond attributes have focused parse/write regressions; typed semantics and RDKit/Indigo bidirectional fixtures remain open |
| A5 / P0/P6 | A5.1–A5.3: absolute gold, genuinely unused data, sample-size/protocol and statistical evaluator; A5.4–A5.5: independent review and formal comparison | Preparation is `local-open`; non-maintainer adjudication is `external`. The four gold candidates and two blind placeholders still reference exposed inputs and lack absolute labels/review. `scripts/evaluate_a5_paired.py` now provides fail-closed category accounting and cluster bootstrap with unresolved rows retained; it is evaluator infrastructure, not independent accuracy evidence. Manifest integrity is not independent accuracy evidence |
| A6 / P5 | A6.1–A6.5: representative rings, typing/charge/all terms/gradients, soundness, stereo and conformer-quality non-inferiority | `local-open` / `local-toolchain`; experimental references/review may be `external`. Issue #337's typing/charge residuals remain closed. Public v1.0.20 closes the historical production-quality cohort of 12 typed failures and four clash rows: stereo-safe MMFF94 is 265/265 independently sound, stereo-clean, and clash-free versus RDKit's 264/265. Its paired speed is 0.944x (95% lower 0.861x), so it is not an MMFF94 speed win. Public UFF best-of-10 is 265/265 usable at 2.359x (lower 2.198x). Rows 166 and 231 remain above 5 kcal/mol in the separate same-coordinate energy gate; term parity, timeout, broader conformer quality, and MMFF94 speed remain open under #637. Evidence: `benchmarks/2026-09-23-public-package-fingerprint-3d-v1.0.20.md`, `validation/results/mmff94-issue337-resolution-v1.0.19.json`, and `validation/results/mmff94-same-explicit-h-energy-current-main-v1.0.19-2026-09-23.json`. |

The original dependency chain is T0 evidence/budget repair → T1 Compatibility Contract
→ T2 release/docs synchronization → T3 WASM/T4 parser safety/T5 stereo → T6
independent review, maintenance and existing 3D gaps. The
[Trust Release plan](trust-release-plan.md) defines subtasks and review dates;
the active roadmap now orders work by the remaining exits, not completed stages.
All A1–A6 exits remain open where stated; A0 is complete and A5 local
preparation continues in parallel. No new 3D feature race is required, but known
3D incorrectness is not waived.
A confirmed silent-corruption regression takes precedence in any area.

| Delivery ID | Product Phase / accuracy | Dependency classification |
|---|---|---|
| T0 | P0/P2; A0/A4 | `local-open`: independent baseline/candidate packet, residual classification, and hybrid adoption audit; measured exploration-budget propagation and footer accounting are implemented and regression-tested |
| T1 | P0/P2/P3; A0–A4 | `local-open`: profile manifest/dashboard; oracle installation is `local-toolchain`; unavailable corpus rights are `external` |
| T2 | P0/P6; A0 | `local-open`: metadata/package/docs checks; registry credentials, online publication and separate-site deployment are `external` |
| T3 | P3; A3/A4 | `local-open`: Worker/MCP/Explorer contracts; missing browsers/build targets are `local-toolchain` |
| T4 | P1/P3; A0/A4 | `local-open`: sourced corpus and isolated runner; unavailable engines/sanitizers are `local-toolchain`; redistribution rights may be `external` |
| T5 | P2/P4; A2/A5 | `local-open`: permutations, round trips and existing residuals; independent absolute-label review is `external` |
| T6 | P5/P6; A5/A6 | `local-open`: maintenance and existing 3D fixes; non-maintainer adjudication/publication is `external` |

These are planning dispositions, not assertions that a tool or credential is
currently missing. Assess each dependency when executing the corresponding task.

The September 16 preflight records a post-freeze 10,000-molecule cohort split
into 2,000 development and 8,000 sealed evaluation inputs:
`validation/results/sealed-cohort-preflight-trust-eval-candidate-20260916.json`.
Its annotated tag is `trust-eval-candidate-20260916`; the raw inputs are local-only
and scores remain unrun. Acquisition/attestation/overlap auditing and that
candidate freeze are complete, not the full acceptance packet. Newer builds need
their own freeze/exposure audit. The 300-case development stereo suite is not
the separate unused challenge or independent gold set.

Declared compatibility requires zero mismatches and 100% coverage on the
frozen valid scope. Independent equivalence requires a reviewed gold set and
the entire paired 95% accuracy-difference interval within ±0.1 percentage
point, with coverage safeguards; superiority requires a positive lower bound.
Insufficient evidence, safe refusal, or a successful local packet alone cannot
complete those exits. A6 uses separate predeclared numerical and quality gates.
The former 1.10x speed stretch is not reactivated.

## Existing area disposition

| Roadmap area | Open work | Disposition |
|---|---|---|
| P0 | Official `@rdkit/rdkit` browser comparison gate | The v1.0.12 Node/WASM rerun against `@rdkit/rdkit@2025.3.4-1.0.0` records 1,000/1,000 RDKit-compatible ECFP4/Morgan matches plus same-process initialization, parse, write, timing, artifact-size, digest, and RSS values in `validation/results/competitive-benchmark-rdkitjs-2026-09-11-v1.0.12.json`. The existing Node/Chromium evidence remains a bounded historical slice; later stable official releases, cross-platform memory, and broader feature/API parity remain open |
| P1 | CIP order-invariance corpus for negative-charge resonance | A checked-in four-case corpus (`validation/cip_negative_resonance_order_invariance.jsonl`) now covers oxygen, nitrogen, fluorine, and sulfur branches attached to the anionic cyclohexadienyl neighborhood; each case passes sixteen deterministic atom-order permutations (three fixed shapes plus thirteen seeded shuffles) with the mapped center label preserved. Broader generated permutations, independent oracle checks, and explicit fail-closed cases remain open |
| P1 | Browser API and package-experience comparison | The v1.0.12 Node/WASM artifact is regenerated and the shared 57-operation contract passes locally. The current Chromium Playground smoke also covers the guided launcher, report JSON action, share-link restoration, malformed/oversized URL fragments, and existing chemistry flows. Exact canonical-string parity, equivalent RDKit batch/cancellation APIs, broader serialization coverage, cross-platform runs, and later package versions remain open |
| P2 | V3000 complex-structure interoperability | The core parser/writer preserves V3000 `SGROUP` logical lines as opaque metadata, accepts either `SGROUP`/`COLLECTION` order, and emits SGROUP before COLLECTION. The current same-process v1.0.12 V3000 contract passes 20 repetitions with zero structural-signature mismatches. Typed SGROUP semantics, Indigo/RDKit cross-engine fixtures, and relaxed editor-preview parsing remain open |
| P1 | Exhaustive malformed corpus across all interchange formats | The current v1.0.13 rerun covers 800 malformed-input attempts (800 unique payloads; 80 attempts per format), 10 oversized cases, and 20 gzip cases with zero failures. The twelve checked-in base cases per format and all eight generated parser-entry families are represented; 480 generated cases have zero duplicate reuses and non-empty typed failure kinds for all 80 format/category combinations. Current evidence is `benchmarks/2026-09-11-streaming-safety-v1.0.13.json`, `validation/results/streaming-parser-entry-categories-v1.0.13.json`, and `validation/results/streaming-parser-entry-failure-kinds-v1.0.13.json`; the v1.0.12 taxonomy and supplemental fixtures remain historical context. The common PDB runner remains intentionally lenient and exposes only its resource-limit failure boundary; the strict PDB API now has Rust, Python, and generated Node/WASM coverage for serial, residue sequence, and coordinate fixed-column parse errors. CML now has a separate opt-in strict parser boundary while the historical lenient path remains unchanged; its quote-aware structural scanner preserves `>` inside attributes, treats comments/CDATA as opaque, rejects multiple root elements and unterminated quoted attributes/comments/CDATA sections, with focused regression coverage, while the current regenerated Node/WASM artifact executes the strict valid case and rejects missing, empty, unbalanced, and mismatched-tag cases. Broader strict fixtures and exhaustive parser-state coverage remain local open work |
| P1 | Equivalent-operation and same-process cross-engine comparison | Same-process Python contracts cover SDF, V2000 MOL, V3000 MOL, MOL2, XYZ, Extended XYZ, PDB, and CDXML: 20 identical-input repetitions each with zero valid-record/failure or applicable signature mismatches. The v1.0.13 ten-format matrix (`validation/results/cross-engine-matrix-v1.0.13.json`) records matching record counts and zero failures for chematic plus each available RDKit/Open Babel lane across SDF, MOL, XYZ, Extended XYZ, V3000, MOL2, CML, CDXML, mmCIF, and PDB. Open Babel remains CLI-only, so its rows are subprocess boundary evidence and not same-process speed evidence. Equivalent timing across all ten formats, full malformed/signature parity, and broader binding-level strict-mode coverage remain open |
| P2 | Full canonical atom-order and E/Z invariance | The focused local gate and full `chematic-smiles` suite pass. The corrected independent probe covers 200 descriptor-census molecules and eight randomized traversals per molecule with 0/200 structural-correctness failures; one RDKit-generated traversal was excluded because its own InChI differed from the source. The clean-main 1,024-relabeling Issue #149/#503 audit now covers all 28 coupled components with 0/28 divergence and zero correspondence failures. The four components from the historical 4/28 result have permanent equivalent-spelling regressions. Historical 4/28 evidence and the rejected 7/28 close-side experiment remain available for traceability. This closes the measured Issue #149 residual, not the broader A2 identity program; unmeasured or over-budget coupled shapes retain the fail-closed stable-key boundary. Current evidence: `validation/results/ez_shared_carrier_coupling_mechanism_audit_summary_1024_2026-09-22.json` and matching JSONL detail. |
| P2 | Analytic force gradients, prepared topology, neighbor lists, periodic/Ewald paths, symmetry-heavy optimization | The bounded Ewald real-space, non-periodic Coulomb, periodic-neighbor, UFF, and prepared MMFF94 analytic-gradient slices retain their recorded parity evidence. Candidate `af7c0c44` promotes the gated analytic MMFF94 path with 80-step history and a 300-iteration budget. Candidate `dd7fe3e9` adds an admissibility predicate to bounded analytic line search and uses it only when expanded and heavy-only stereo evaluation disagree, preserving 265/265 stereo-safe outputs without committing rejected trial steps. The quality lane is slower than RDKit, and two same-coordinate energy residuals plus broader term/timeout/conformer exits remain; no general 3D-performance or quality win is claimed. |
| P3 | Broader public API manifest coverage | The declared shared surface is tracked by one versioned manifest with 57 four-binding operations and zero source-only operations. Strict CML is covered by regenerated Rust/Python/Node/WASM artifacts and shared fixture tests. Broader public-API coverage remains a separate scope-expansion task and is not part of the current unchecked roadmap count |
| P3 | Full clean-install/cross-platform/browser-memory evidence | The local macOS arm64 clean-install, cold-start, throughput, peak-RSS, and v1.0.12 Node/WASM artifact regeneration slice is complete and linked from `ROADMAP.md`; cross-platform/browser-memory evidence remains open. The current artifact report `validation/results/wasm-node-artifact-rebuild-v1.0.12.json` records matching 263/263 source/artifact exports, wasm-pack 0.13.1, wasm-bindgen 0.2.121, and wasm-opt 130. `scripts/bench_browser_wasm_vs_rdkit.mjs` now records Chromium `performance.memory` checkpoints when exposed, explicitly excluding WASM linear memory; no browser-memory result is claimed because the local Chromium DOM-dump runner did not terminate reliably |
| P3 | Broader browser/agent adversarial matrix | Current Node-hosted agent-side pipeline/MCS adversarial rerun passes both suites and is recorded in `validation/results/wasm-agent-adversarial-v1.0.13.json`. The current Chromium, Firefox, and WebKit browser smoke also pass the Playground launcher, Explorer malformed-record/cancellation flows, share-link restoration, report JSON action, malformed/oversized URL fragments, and existing chemistry flows; broader agent coverage and browser-memory measurements remain open |
| P4 | Reaction/SMARTS/medchem breadth and curated quality reports | The local `chematic-rxn` v1.0.13 regression gate passes 224/224 unit tests plus the shared five-test contract. The bounded application contract covers the prior fifteen cases and now adds no-match-empty, reactant-count-mismatch, and invalid-template boundaries across Rust/Python/Node/WASM; evidence is recorded in `validation/results/reaction-quality-boundary-v1.0.13.json`. Map-number mismatch diagnostics now sort offending IDs and report the smallest symmetric-difference ID, removing process-dependent `HashSet` ordering from machine-readable errors. Full query-aware application semantics, broader medicinal-chemistry coverage, exhaustive matcher-equivalence testing, and curated precision/recall/timeout quality reports remain open |
| P5 | MMFF94/UFF typing, parameters, convergence, stereo, and independent-oracle coverage | Public v1.0.20 confirms 265/265 independently sound, stereo-clean and clash-free MMFF94 stereo-safe outputs; RDKit has 264/265 usable rows in the same scorer. The paired speed is 0.944x (95% lower 0.861x), so quality-equivalent MMFF94 remains a speed gap. Public UFF best-of-10 is 265/265 usable and 2.359x faster (lower 2.198x). Energy/term, timeout, and broader conformer-quality gates remain open. Evidence is recorded in `benchmarks/2026-09-23-public-package-fingerprint-3d-v1.0.20.md`, `validation/results/mmff94-issue337-resolution-v1.0.19.json`, `validation/results/uff-element-finiteness-v1.0.10.json`, and `validation/results/mmff94-rdkit-availability-oracle-v1.0.10.json`. |
| Performance stretch target | Additional 1.10x SMILES parse speed | Explicitly abandoned by the maintainer on 2026-09-08; historical context only, not a release gate |
| S5 | Non-maintainer review or external audit | External coordination; local review packet is complete, but this remains open |
| S6 | Advisory intake, backport, publication, and supported-version rehearsal | External coordination and future publication; remains open |

Evidence note: the current `chematic-ff` suite count above includes the
reversed-endpoint UFF torsion regression. It is a regression count, not a claim
that the broader P5 force-field parity item is complete.

The direct UFF convergence contract was tightened during this disposition:
low-gradient unsound stationary points are now reported with
`converged: false`, while `sound`, `worst_bond_length`, and
`rejected_unsound_step` remain available for diagnosis. This is a bounded
false-success fix and does not close the broader P5 force-field item.

The active autonomous scope is the local implementation and evidence work,
subject to available local toolchains. Environment-dependent evidence is
recorded as such and is never promoted to a release or parity claim. S5/S6 are
intentionally excluded from autonomous completion because they require external
parties or future publication actions.

## Dependency classification

Active work has three dependency classes; a fourth retains abandoned history:

| Class | Meaning | Current items |
|---|---|---|
| `local-open` | Can be implemented and verified in this repository without an external service, reviewer, future release, or user decision | exhaustive malformed/parser-state coverage; full canonical atom-order/E/Z invariance; remaining analytic force terms and electrostatic/full coordinate-dependent neighbor-list integration; broader triclinic periodic cell-list proof; broader stable-operation fixture coverage; broader browser/agent cases when the local runner is available; reaction/SMARTS breadth after the P3 gate |
| `local-toolchain` | Repository work is local, but its required runtime or package must be available and verified for that lane | candidate-specific Node/WASM and Python artifact regeneration; browser-memory and cross-platform execution; stable browser runner termination. Prior successful builds do not establish availability or evidence for every future lane |
| `external` | Requires a non-maintainer review, external audit, future maintenance event, publication, or external coordination | S5 independent review/audit; S6 advisory/backport/publication rehearsal; curated quality reports that require an agreed external/reference corpus |
| `historical` | Intentionally retained for traceability but abandoned or not a current completion gate | the additional 1.10x SMILES parse-throughput stretch target |

Only `local-open` and `local-toolchain` work is in the autonomous scope. A
`local-toolchain` item is not complete until the required local toolchain is
available and its artifact or runtime evidence passes; an `external` item is
never marked complete by local tests alone.

## Historical implementation notes

The following versioned results explain completed slices and earlier blockers;
they are not a current toolchain inventory or new v1.0.14 measurements.

The prepared MMFF94 slice also has an opt-in analytic L-BFGS reduction smoke;
its boundary is recorded in `validation/results/mmff94-bounded-analytic-lbfgs-v1.0.10.json`.
The vdW and opt-in electrostatic terms rebuild coordinate-dependent cutoff
neighbor lists and match their all-pair cutoff references on fixed fixtures;
see `validation/results/mmff94-vdw-neighbor-list-v1.0.10.json` and
`validation/results/mmff94-electrostatic-neighbor-list-v1.0.10.json`.
The opt-in analytic L-BFGS objective switches the cutoff energy and gradient
together; the production-default minimizer remains unchanged.

The current UFF implementation now includes a bounded common-organic torsion
slice for sp3–sp3, sp2–sp2, and mixed central bonds. The butane energy and
analytic-gradient regression passes, and degenerate dihedrals fail closed.
Full metal-class parameter coverage and independent RDKit energy parity remain
open.

The distance-descriptor audit also fixed a compact-indexing defect for explicit
H/D/T graphs: AutoCorr2D, Moran, Geary, and IPC now translate the heavy-atom
distance-matrix positions back to original atom indices. The regression case
`[2H]C([2H])([2H])NC=O` passes the full `chematic-chem` suite (854 passed,
one ignored) and clippy; this is a correctness fix outside the seven tracked
GitHub Issue acceptance gates.

The torsion builder no longer skips bonds whose stored endpoint order is
`atom1 > atom2`; bond iteration is already unique, and a reversed-endpoint
C-C-C-C regression confirms that torsion coverage is independent of graph
construction order.

The full local gate was rerun on 2026-09-11 with `TMPDIR=/private/tmp bash
scripts/check.sh` for the static and unit lanes, followed by
`TMPDIR=/private/tmp CARGO_INCREMENTAL=0 cargo test --workspace --tests --quiet`
for the integration lane after reclaiming the repository's regenerable
`target/debug` cache. All Rust unit and integration targets completed with zero
failures, and the static,
format, clippy, binding-manifest, strict-PDB-binding, benchmark-index,
streaming-safety-manifest,
cross-engine streaming matrix, MMFF94 RDKit availability oracle,
streaming failure taxonomy, Issue #503 E/Z residual evidence,
roadmap-disposition, and workflow-pin checks
passed. The default advisory `deny` step cannot acquire its lock because the
Cargo advisory database is under a read-only path, but a writable-copy rerun
against advisory-db revision `5a0ebedfe8bdd2e295b171f4162f8c977bcad9a5`
passes advisories, bans, licenses, and sources; see
`validation/results/cargo-deny-v1.0.10.json`. A source-built v1.0.12 wheel is
also reproducible offline in the project virtual environment; the strict CML
shared contract passes after rebuilding the wheel, while the system Python
installation is not treated as workspace evidence. The wheel is a local
candidate only and is not published or installed into the normal user
environment. Neither condition is promoted to a completed roadmap item.

The P1 streaming safety runner was also rerun against the then-current v1.0.13 workspace:
800 malformed-input attempts (800 unique payloads; 80 per format), 10
oversized-input cases, and 20 gzip cases passed with zero
failures across SDF, MOL, XYZ, Extended XYZ, V3000,
MOL2, CML, CDXML, mmCIF, and PDB. This remains bounded safety evidence;
exhaustive malformed
corpus coverage and same-condition cross-engine equivalence remain open.
