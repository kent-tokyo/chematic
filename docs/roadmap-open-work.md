# Roadmap open-work disposition

This document is the current disposition of unchecked roadmap items for the
v1.0.17 release follow-up (updated 2026-09-20). An item is not complete merely
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
| T0 / A4 | Independent baseline/candidate packet; 12 → 21 residuals, candidate rejected | `local-open`: classify residuals and evaluate a replacement; do not repeat the withdrawn 12 → 3 claim |
| T1.6 / A0 | Two annotated 2k/8k candidates have now been evaluated and rejected. The September 16 candidate had 6 molecular-weight unsupported values and 46 TPSA strict mismatches. The `trust-eval-candidate-20260920` candidate (commit `8a0b28a`, RDKit 2025.09.3) parsed 8,000/8,000 with 0 unsupported values and had 7 fields at 8,000/8,000 strict, but molecular weight was 7,998/8,000 strict (two isotope-table residuals; max error 0.0114336 Da). | `historical`: both raw results remain local-only and exposed; neither holdout may tune or validate a later candidate. See `validation/results/sealed-candidate-trust-eval-20260920-summary.json`. A later adoption claim requires another unused freeze and input. |
| T1.7 / A4 | Ordinary V3000, external RDKit/Indigo readers and SGROUP/COLLECTION baseline (PR #544) | `local-open`: typed chemistry/edit semantics and declared refusal boundaries; general CX label retention does not complete T1.8 attachment identity/collapse |
| T1.8 / A4 | September 19: complete positive `_AP<n>` syntax parser and degree-one-wildcard identity API, with invalid/non-dummy/multi-degree rejection | `local-open`: external RDKit round-trip and any collapse semantics; identity is deliberately not collapse eligibility |
| T3.4 / A3/A4 | Public 1.0.15 vs official RDKit npm 10k, three browsers, 20 repetitions; separate three-run RSS diagnostic | `local-open` / `toolchain-open`: public 1.0.17 and candidate lanes, equivalent output/search, resource definitions; unavailable remote host/measurement stays explicit |
| T3.6 / A3 | PR #555 source: exact on all 9,999 supported fixed-corpus rows in Chrome/Firefox/WebKit, exact on 5,000 independent ChEMBL rows in Chrome, and faster than pinned RDKit.js in all local browser lanes (2.06–6.17x paired median). Hosted run 35483926941 passes the paired-95% gate (lower bounds 2.37–3.65x); Rust/Python/Node-WASM cross-binding is 5,000/5,000 | `local-open` / `toolchain-open`: audit original 20-pair/three-session, bootstrap, p95 and resource conditions against shorter t-interval runs. `external-open`: publication then actual-package rerun. Preserve search regressions and native ECFP profile |
| T4.6 / A2/A4 | September 19 core fix plus PR #562: parser/writer use bounded `%(n)` labels ≥100; malformed/overflow legacy ambiguity rejects; a 121-closure graph round-trips and shared Rust/Python/Node/WASM validity fixtures accept `%(100)` while rejecting malformed legacy spellings | `local-open`: canonical-writer atom-order and resource-limit coverage remain; this is parser validity coverage, not external-engine round-trip evidence |
| T5.6 / A2 | Upstream CIP/atrop changes identified; existing development suite is not their acceptance | `local-open`: selection/isotope/atrop/permutation regressions; independent gold is `external-open` |
| T1.5 next-oracle rebaseline | RDKit 2026.03.6 is the latest GitHub stable release verified September 20; `validation/rdkit_rebaseline_manifest.json` now fixes Python/native/npm lanes and an explicit unavailable next-stable lane | `local-open`: add per-artifact backend/config/corpus/toolchain provenance and old/new diff classification. `external-open`: wait for each actual artifact; no change to the frozen T1.6 oracle |
| T5.7 + T1.9 / A2/A4 | PR #558: atom-map potential-stereocenter invariance covers identity/reordered atoms; the follow-up bicyclic-amine gate exercises 32 fixed graph-storage orders through SSSR, potential-center/CIP, clone, canonicalize and reparse without inventing the tertiary-N center. V3000-to-V3000 preserves opaque query attributes and unmodelled query constraints refuse conversion to ordinary formats | `local-open`: nearby negative cases and label-level stereochemistry checks; E/Z 2×2 query truth tables across serialization; unsupported paths remain refusal boundaries, not query compatibility |
| T3.7 / P3 | PR #559: shared Rust/Python/Node/WASM fixtures now enforce `input_count = success + failed + refused + skipped`, original row indices, terminal status, parse error stage, and `all_succeeded()` | `local-open`: cancellation, streaming/unknown-length input, retry, Worker/MCP adapters, migration examples and measured resource limits remain |

The [Trust execution packet](trust-release-plan.md#7-2026-09-20-trust完了に向けた実行順)
owns the new subtask exits. These are subdivisions of the existing 18 open areas,
not additional phases or a reason to mark A0–A6 complete.
The [performance plan](parse-morgan-performance-plan.md) owns T3.6 / PF0–PF5.
Historical 1.0.15 timings are not evidence of 1.0.17 performance; speed acceptance
requires exact supported-domain bits, preserved coverage, and paired uncertainty.
The sealed accuracy cohort must not be used for performance tuning.

## Accuracy follow-up (A0–A6)

All seven packages retain open exit work. The 2026-09-12 plan audit preserves
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

| Work package | Remaining work | Dependency class and current evidence |
|---|---|---|
| A0 / P0 | A0.1–A0.5: all-field strict holdout, raw-result accounting, pinned source/binary provenance, adversarial gate tests, separately built baseline/candidate | `local-open`. Source-built Python with RDKit 2025.09.3 passes all 40,000 field checks strictly; the 12-row all-field holdout passes 96/96. The 5,000-row ECFP4 binding gate passes 5,000/5,000 on all three bindings and records Python extension and Node/WASM hashes. The 8-field descriptor binding gate now passes 5,000/5,000 and 7,737/7,737 exposed rows across Rust/Python/Node/WASM; this is binding parity, not an RDKit oracle. Negative coverage rejects bool counts, empty successes, extra fields, missing provenance, duplicate IDs, and non-finite values. The September 20 candidate had fresh provenance, an annotated tag, unused-data attestation, and 8,000/8,000 parsed rows, but was rejected at 7,998/8,000 strict molecular-weight agreement; its input/result are exposed/local-only. A later adoption decision requires a new freeze, genuinely unused cohort, provenance/overlap checks, and full scorecard automation. |
| A1 / P1/P2 | A1.1–A1.4: common perception and atom-level center FP/FN; A1.5: additional descriptor families adopted separately | `local-open`. The eight-descriptor compatibility profile is now `adopted_opt_in`: source-built core fields pass 5,000/5,000 strictly, the all-field holdout passes 96/96, and native-default preservation passes. A provenance-preserving carbon-branch fallback plus a narrow RDKit-compatible ring-constrained sp3-N rule yields potential-center parity 5,000/5,000; the expanded 7,737-row exposed holdout also passes 7,737/7,737 with atom-level FP=0/FN=0. The bounded chordless-aromatic-cycle candidate now passes all eight descriptor fields 7,737/7,737 strictly against pinned RDKit 2025.09.3, with promotion decision `adopted_opt_in` and native defaults preserved. A1.1 uses physical nuclide masses for known isotope labels and RDKit's mass-number fallback for unlisted valid labels; this behavior is independently covered by the descriptor unit gate. Morphine/codeine aromatic-oxide LogP residuals are resolved while the 6-member cyclic-ether negative remains protected. After the P-H phosphonate, thioamide-N, and charged-sulfur fixes, fresh descriptor and ChEMBL rechecks each pass all eight fields 5,000/5,000 strict; the current-source 7,737-row exposed holdout is also strict-exact for all eight fields with parse/unsupported failures 0. Evidence: `validation/results/descriptor-rdkit-a1.1-tpsa-fixed-v1.0.13.json` and `validation/results/descriptor-rdkit-a1.4-recheck-v1.0.13.json`. The latest 8-field Rust/Python/Node-WASM binding regression also passes 7,737/7,737; affected local perception/chem/fingerprint/WASM/SMARTS/3D suites pass 204/864/309/344/185/597. These are binding and regression consistency results, not independent chemical correctness. Potential-center results above predate this shared-perception change (fixes3). The center runner records five explicit RDKit atom-set fallbacks and zero unresolved oracle rows. The implementation excludes ordinary tertiary amines, fused degree-3 junctions, amide/sulfonamide environments, and aromatic N. Unused/sealed evaluation, workflow-level candidate acceptance, broader A1.1 per-atom corpus, and additional descriptor families remain. Count equality alone does not prove center-set equality |
| A2 / P2 | A2.1–A2.5: phosphorus CIP, independent label checks, current K=1,024, semantic preservation, permutations, idempotency and key collisions | `local-open`, with independent adjudication in A5. Existing phosphorus remains fail-closed; canonical RDKit recanonicalization probe is 200/200. The current source-provenance-pinned #503 K=1,024 audit reproduces 4/28 divergent components with 0 cross-correspondence failures; three current aromatic-stash residuals have explicit fail-closed and relabeling regression tests, while the common traversal mechanism remains open. `num_unspecified_stereocenters` now reuses Rule-5-aware potential-center selection and excludes ordinary CH2/CH3 false positives. The fail-closed result is validated by `scripts/check_ez_residual_evidence.py` against `validation/results/ez_shared_carrier_coupling_mechanism_audit_summary_1024_2026-09-12.json`. A dedicated Node/WASM regression also confirms the representation-unstable phosphorus case returns zero accurate assignments plus `oracleUnstable` entries, matching the Python boundary. The runner now accepts provenance-carrying JSONL rows; the 181-row residual classification corpus recheck has zero potential-center count mismatches and zero parse failures, but does not close CIP label adjudication. The new atom/bond label audit reports atom-map agreement 155/181, E/Z bond-map agreement 14/14, all 26 atom mismatches confined to P + `oracle_unstable`, and no assigned-label mismatch. The exposed descriptor, ChEMBL, and NCI corpus gate now passes standardize→canonicalize idempotency 14,999/14,999; this is regression evidence only. Evidence: `validation/results/cip-residual-stereocenter-v1.0.13.json`, `validation/results/cip-label-parity-a2-v1.0.13.json`, `validation/results/canonical-idempotency-a2.4-current-v1.0.13.json`. |
| A3 / P2/P3 | A3.1–A3.5: raw FP oracle, k=1/10/100 and threshold edges, unrounded scores/order, shared input IDs, independent RDKit plus every binding, unused evaluation | `local-open`. Rust/source-built Python/Node-WASM agree at k=1/10/100 for 500 queries × 4,500 rows with 0 pairwise mismatches and unrounded f64 scores; an independent exhaustive RDKit oracle also passes all k. The shared fixture has 34 success cases plus one explicit unsupported-bond error case, passing Rust/Python/WASM. The inclusive threshold gate passes 96 RDKit-derived boundary cases across all three bindings; finite `[0, 1]` validation and one-ulp equality handling are explicit. A separate RDKit-parity aromaticity lane resolves all 59 Hueckel radius-1 raw Morgan residuals on the checked-in 5,000-row corpus: 5,000/5,000 exact against pinned RDKit 2025.09.3, with zero parity preprocessing errors. Evidence: `validation/results/ecfp-rdkit-raw-identifier-parity-aromaticity-variant-2025.09.3-5000.json`. This is candidate evidence, not yet an adopted default: raw trace/provenance retention, post-adoption all-binding/search reruns, and genuinely unused-input evaluation remain open |
| A4 / P1/P4 | A4.1–A4.5: standardization semantics, SMARTS mapped sets, reaction products and typed V3000 metadata in both directions | `local-open`; uninstalled Indigo is `local-toolchain`, data permissions or external fixtures may be `external`. Footer-verified persistent-corpus comparison over 5,021 molecules × 31 queries covers 155,651 cells: parity target-atom-set matches 145,588, RDKit SMARTS parse errors 10,042, residuals 21, and atom-alignment failures 0. Evidence: `validation/results/rdkit-smarts-direct-chembl-5000-footer-verified-v1.0.13.json`, including input/dump SHA-256 and comparison configuration; `scripts/check_rdkit_smarts_evidence.py` validates the source corpus hash, footer, row accounting, alignment, and bucket sum. Earlier private-corpus numbers are not mixed into this denominator; the 12→3 claim remains withdrawn. Matching compares target atom sets, not all mapped embeddings, and is not a frozen-build adoption packet. The shared path now propagates the caller budget and reports measured candidate counts, but independent baseline/candidate acceptance remains open. Expanded 51-case/14-template atom-map-aware canonical product-set gate passes 51/51; the opposite-enantiomer case remains a separate 1/1 fail-closed safety gate. V3000 SGROUP/COLLECTION, isotope, enhanced-stereo boundary, and opaque `ENDPTS=`/`ATTACH=` bond attributes have focused parse/write regressions; typed semantics and RDKit/Indigo bidirectional fixtures remain open |
| A5 / P0/P6 | A5.1–A5.3: absolute gold, genuinely unused data, sample-size/protocol and statistical evaluator; A5.4–A5.5: independent review and formal comparison | Preparation is `local-open`; non-maintainer adjudication is `external`. The four gold candidates and two blind placeholders still reference exposed inputs and lack absolute labels/review. `scripts/evaluate_a5_paired.py` now provides fail-closed category accounting and cluster bootstrap with unresolved rows retained; it is evaluator infrastructure, not independent accuracy evidence. Manifest integrity is not independent accuracy evidence |
| A6 / P5 | A6.1–A6.5: representative rings, typing/charge/all terms/gradients, soundness, stereo and conformer-quality non-inferiority | `local-open` / `local-toolchain`; experimental references/review may be `external`. Current-source MMFF94 type IDs match 6,681/6,698 (99.76%); excluding the declared unsupported probe, 16 comparable residual atoms remain. Current BCI charge comparison matches 6,665/6,693 (99.58%) on comparable rows, with 28 residual atoms across three macrocycles. On current main, the low-level production strict bond+angle harness now succeeds 265/265 (missing parameter/type/minimization failures all 0), closing Issue #227's stated coverage gap; this is not full force-field or pipeline parity. The row-isolated hard-timeout rerun covers all 265 rows with 215 successes and 50 explicit timeouts; all successes are finite/sound with zero gross clashes, while convergence is 54/215 successes at 200 steps (Tier A 42/64 and Tier B 12/151). Same-heavy-coordinate RDKit diagnostic over 215 successes reports median absolute energy delta 29.93 kcal/mol, p90 55.89, max 920.36. The fixed-H/same-coordinate diagnostic corrected the buffered 14-7 formula, explicit O-H numeric typing, explicit-H SymmSSSR representation dependence, and delocalized amine N-H typing. A narrow macrocycle-boundary candidate reduces the maximum comparable energy delta to 8.7467 kcal/mol and raises within-5 to 260/262; 0009 is exact in the macrocycle type audit, while 0029/0030 remain residuals. The latest 262/265 run reports median absolute energy delta 0.226 kcal/mol, p90 1.219, max 8.75, 224/262 within 1 kcal/mol, and 260/262 within 5 kcal/mol. Remaining macrocycle/drug-like term residuals and convergence/gradient/stereo exits remain open. Evidence: `validation/results/mmff94-strict-gate-current-main-v1.0.17-2026-09-20.json`, `validation/results/mmff94-hard-timeout-pipeline-a6-v1.0.13.json`, `validation/results/mmff94-rdkit-same-heavy-energy-a6-215-v1.0.13.json`, and `validation/results/mmff94-same-explicit-h-energy-a6-265-after-refined-macrocycle-boundary-v1.0.13.json`. Itemized energy and term-gradient diagnostics expose cancellation and high bond/nonbonded residuals in remaining stress cases |

The original dependency chain is T0 evidence/budget repair → T1 Compatibility Contract
→ T2 release/docs synchronization → T3 WASM/T4 parser safety/T5 stereo → T6
independent review, maintenance and existing 3D gaps. The
[Trust Release plan](trust-release-plan.md) defines subtasks and review dates;
the active roadmap now orders work by the remaining exits, not completed stages.
All A0–A6 exits remain open where stated; A0 unused-data work and A5 local
preparation start in parallel. No new 3D feature race is required, but known
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
| P2 | Full canonical atom-order and E/Z invariance | The focused local gate and full `chematic-smiles` suite pass, including the 2026-09-11 aromatic direction-anchor hardening. The corrected independent probe covers 200 descriptor-census molecules and eight randomized traversals per molecule: 0/200 structural-correctness failures, with one RDKit-generated traversal excluded because its own InChI differed from the source molecule. The current source-provenance-pinned 1024-relabeling Issue #503 audit covers all 28 coupled components and reproduces 4/28 divergent outputs with zero correspondence failures; the summary is `validation/results/ez_shared_carrier_coupling_mechanism_audit_summary_1024_2026-09-12.json`, and the detailed rows are in the matching JSONL artifact. A close-side guard relaxation experiment worsened the residual to 7/28 and remains rejected in `validation/results/ez_shared_carrier_coupling_mechanism_close_side_experiment_2026-09-11.json`. `scripts/check_ez_residual_evidence.py` validates both the current provenance boundary and the rejected experiment. This is stable residual evidence, not a complete atom-order-invariance proof; affected aromatic-stash systems still intentionally fail closed |
| P2 | Analytic force gradients, prepared topology, neighbor lists, periodic/Ewald paths, symmetry-heavy optimization | The bounded Ewald real-space slice is complete: PME now uses a validated uniform cell list with the configured cutoff, and package tests cover direct-reference parity, negative coordinates, and the strict cell-boundary rule. The public non-periodic direct Coulomb cutoff path now has the same bounded cell-list treatment, with independent all-pair parity recorded in `validation/results/ewald-coulomb-cutoff-neighbor-list-v1.0.10.json`. The periodic-neighbor cell-list slice now covers orthorhombic plus 4 hand-selected and 12 fixed-seed generated triclinic cells, with exact result-key parity and deterministic ordering against the existing search; the repaired triclinic evidence is in `validation/results/periodic-neighbor-triclinic-cell-list-v1.0.10.json`, while the earlier cutoff-boundary rejection remains preserved as historical evidence. UFF now prepares bond/angle/vdW topology and uses analytic gradients verified against finite differences. MMFF94 now has bounded prepared bond/angle, stretch-bend, nonbonded, non-singular torsion, and non-singular Wilson out-of-plane analytic-gradient building blocks, plus an aggregate prepared-gradient gate matching total energy and an opt-in analytic L-BFGS reduction smoke; the coordinate-dependent MMFF94 vdW and opt-in electrostatic lists match their all-pair 10 Å cutoff references, and the opt-in analytic L-BFGS objective switches both together. Central-difference checks are recorded in `validation/results/mmff94-stretch-bend-gradient-v1.0.10.json`, `validation/results/mmff94-nonbonded-gradient-v1.0.10.json`, `validation/results/mmff94-torsion-gradient-v1.0.10.json`, `validation/results/mmff94-oop-gradient-v1.0.10.json`, `validation/results/mmff94-bounded-analytic-gradient-v1.0.10.json`, `validation/results/mmff94-vdw-neighbor-list-v1.0.10.json`, and `validation/results/mmff94-electrostatic-neighbor-list-v1.0.10.json`. OOP gradients fail closed at zero-area and asin-branch singularities, and the public finite-difference minimization path remains unchanged until a full soundness gate exists. Full production integration and the broader numerical/performance gate remain local open work |
| P3 | Broader public API manifest coverage | The declared shared surface is tracked by one versioned manifest with 57 four-binding operations and zero source-only operations. Strict CML is covered by regenerated Rust/Python/Node/WASM artifacts and shared fixture tests. Broader public-API coverage remains a separate scope-expansion task and is not part of the current unchecked roadmap count |
| P3 | Full clean-install/cross-platform/browser-memory evidence | The local macOS arm64 clean-install, cold-start, throughput, peak-RSS, and v1.0.12 Node/WASM artifact regeneration slice is complete and linked from `ROADMAP.md`; cross-platform/browser-memory evidence remains open. The current artifact report `validation/results/wasm-node-artifact-rebuild-v1.0.12.json` records matching 263/263 source/artifact exports, wasm-pack 0.13.1, wasm-bindgen 0.2.121, and wasm-opt 130. `scripts/bench_browser_wasm_vs_rdkit.mjs` now records Chromium `performance.memory` checkpoints when exposed, explicitly excluding WASM linear memory; no browser-memory result is claimed because the local Chromium DOM-dump runner did not terminate reliably |
| P3 | Broader browser/agent adversarial matrix | Current Node-hosted agent-side pipeline/MCS adversarial rerun passes both suites and is recorded in `validation/results/wasm-agent-adversarial-v1.0.13.json`. The current Chromium, Firefox, and WebKit browser smoke also pass the Playground launcher, Explorer malformed-record/cancellation flows, share-link restoration, report JSON action, malformed/oversized URL fragments, and existing chemistry flows; broader agent coverage and browser-memory measurements remain open |
| P4 | Reaction/SMARTS/medchem breadth and curated quality reports | The local `chematic-rxn` v1.0.13 regression gate passes 224/224 unit tests plus the shared five-test contract. The bounded application contract covers the prior fifteen cases and now adds no-match-empty, reactant-count-mismatch, and invalid-template boundaries across Rust/Python/Node/WASM; evidence is recorded in `validation/results/reaction-quality-boundary-v1.0.13.json`. Map-number mismatch diagnostics now sort offending IDs and report the smallest symmetric-difference ID, removing process-dependent `HashSet` ordering from machine-readable errors. Full query-aware application semantics, broader medicinal-chemistry coverage, exhaustive matcher-equivalence testing, and curated precision/recall/timeout quality reports remain open |
| P5 | MMFF94/UFF typing, parameters, convergence, stereo, and independent-oracle coverage | The bounded diversity/failure-rate/RMSD/TFD/energy measurement task is now complete for its declared scopes: local Rust regression 596/596 executed tests, current 250-molecule RDKit TFD boundary (201 finite paired values plus 49 rigid/not-applicable probes), and the existing 58-molecule quality bundle. A local UFF dispatch soundness probe covers 18 explicit metal/halogen types, and a second prepared-parameter fixture covers all 45 declared `UffType` variants with finite energy and central-difference-consistent analytic gradients, plus a genuine 1-4 vdW gradient fixture; the current `chematic-ff` library suite passes 225/225. The independent RDKit MMFF94 availability oracle also passes 265/265 parse, embed, property, force-field construction, and finite-energy checks; its 142/123 minimize return-code split is retained without claiming convergence parity. The v1.0.13 low-level production strict-gate remeasurement now reaches 265/265 on the pinned corpus with zero missing parameters, unsupported atom types, or minimization failures after adding the registry-backed C5A equivalence row and routing the concrete `(43,18,63)` residual through the documented Halgren empirical angle fallback. The follow-up gap diagnosis now agrees with the live RDKit oracle for all 10 empirical angle tuples and both empirical bond tuples; the former 9/10 angle result was a checker omission of type 63, not a Rust parameter mismatch. This is a low-level `generate_coords` boundary and not full pipeline_v2 parity; the current audit has zero Torsion `missing` rows, with five explicit `no_term_by_design` rows retained and broader table/parity coverage still open. The #337 D2 candidate pass now gathers roots from one immutable bond graph before duplicate grouping, and the 2026-09-11 follow-up makes duplicate-group and root expansion ordering explicit instead of relying on `FxHashMap` iteration; its machine-readable local determinism evidence records six fixtures and 1,536 seeded checks, while the representative-family parity mismatch remains open. Candidate ordering now also compares a rotation/reversal-normalized ring traversal key (atom rank plus bond order) before the compact rank-multiset key, removing one remaining tie-break ambiguity without claiming RDKit representative-family parity. The diagnostic all-root replay confirms zero candidates missing from the D2 population for all six fixtures; exact-cycle enumeration still exposes 3, 15, 15, 15, 15, and 3 non-basis alternatives for the 28/30/32/32/31/29-atom macrocycles, so the remaining work is representative-family selection rather than root coverage. UFF's bounded soundness diagnostics now propagate through Rust, Python, and WASM, including whether line search rejected an energy-decreasing but geometrically unsound proposal (`rejected_unsound_step`). The current development line also adds a smooth aromatic trigonal-center inversion term with an analytic-gradient regression; general UFF torsion, independent energy parity, and the broader P5 typing/parameter/convergence/stereo item remain open. Evidence is recorded in `validation/results/uff-element-finiteness-v1.0.10.json`, `validation/results/mmff94-issue337-determinism-v1.0.10.json`, `benchmarks/2026-09-09-reaction-3d-focus-v1.0.10.json` / `.md`, `validation/results/tfd-oracle-evidence-v1.0.10.json`, `validation/results/mmff94-rdkit-availability-oracle-v1.0.10.json`, `validation/results/mmff94_strict_gate_remeasure_227_v1.0.12.json`, `validation/results/mmff94_strict_gate_remeasure_227_v1.0.13.json`, `validation/results/mmff94-angle-gap-diagnosis-v1.0.13.json`, and `validation/results/3d-quality-evidence-v1.0.10.json` |
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
