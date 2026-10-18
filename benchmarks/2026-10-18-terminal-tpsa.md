# Terminal N/O and zero-model TPSA source follow-up

Source `f8e9d7dee5a8ece68a61bb9197f616bc86f6bfeb` → `add2b6b445c240b6c313927a98c551587239aff4`: the isolated RDKit-default TPSA change passes both first-use improvement gates and every source median floor. This is one-host source-wheel evidence on two exposed 1,000-row windows, not a published-package or toolkit-wide claim.

## Result

| Window offset | Operation | Source speedup | Paired 95% CI | Faster blocks |
|---:|---|---:|---|---:|
| 0 | first-use `rdkit_tpsa` | 1.036887x | [1.024228, 1.049287] | 20/24 |
| 5000 | first-use `rdkit_tpsa` | 1.270154x | [1.260643, 1.280864] | 24/24 |

All 42 source medians meet the predeclared 0.99 ratio floor. The minimum is 0.991554258 at offset 0 `labute_asa/parse_inclusive`. 19 source medians are below 1; all are retained in the [78-lane CSV](2026-10-18-terminal-tpsa-lanes.csv) and [numeric packet](2026-10-18-terminal-tpsa-numeric.json.gz). The largest source median speed decrease is 0.845%.

Source lanes with their entire interval below 1:

| Window | Lane | Speedup | 95% CI |
|---:|---|---:|---|
| 0 | `morgan2_chiral/parse_inclusive` | 0.995243956 | [0.988897257, 0.999000090] |
| 0 | `tpsa/parse_inclusive` | 0.995302018 | [0.983541747, 0.999175936] |
| 0 | `labute_asa/first_use` | 0.994912874 | [0.987665785, 0.999137159] |
| 0 | `labute_asa/parse_inclusive` | 0.991554258 | [0.983802614, 0.998558949] |
| 5000 | `morgan2_chiral/parse_inclusive` | 0.993387260 | [0.990768653, 0.998404307] |
| 5000 | `maccs/repeated_use` | 0.998439242 | [0.991038296, 0.999919713] |

## What changed

- RDKit-default TPSA sums N/O contributions. Model cleanup can remove graph H atoms and rewrite charges/bonds, but cannot create N/O atoms. A graph without N/O returns positive zero before the model callback. The scalar remains cached on its original molecule; clones, edits and nested model callbacks retain the existing cache rules.
- When every N/O is nonaromatic, has degree at most 1, has no aromatic incident bond, and the RDKit model gate agrees, global aromaticity cannot change those local polar environments. Use the original graph for this case. Already-sanitized callback views keep their prepared global-view route.
- The final runtime delta is limited to TPSA, the additive hidden scalar-shortcut helper, and its two reexports. Unconfirmed H-removal, explicit-valence, Labute coefficient and model-gate allocation experiments are excluded. Parsing behavior follows the baseline. Normal parsing may preseed the model gate; tests also cover unseeded template parsing.

## Correctness

- 1380 related Rust library tests pass; 2 existing ignores remain. Workspace/all-target Clippy passes. 329 selected Python compatibility/stable-ECFP4/tetrahedral tests pass.
- The existing frozen TPSA implementation matches float bits over all 14,072 normal/stress rows in raw and Kekulé forms, including repeated evaluation. New tests cover terminal-environment guards, N/O-free graphs, parser-preseeded versus unseeded caches, cloning, edits and nested callbacks.
- All 56,288 Morgan cells equal both the source baseline and pinned RDKit. All 42,216 scalar values/refusals (`tpsa`, `labute_asa`, `num_rings`) remain source-exact. The normal 10,000-row scalar outputs equal RDKit. Stress ASA has 3 shared unavailable cases; ring count keeps 4 pre-existing RDKit differences. COSMolKit keeps 6 stress parse refusals and its existing scalar/ring/stereo differences; none are filtered.

## COSMolKit deficit screen

The initial seven-operation screen used `bc36b398` and identified first-use TPSA/ASA in both windows, chiral Morgan in the first window, and ring count in the second. First-window ring count differs on 6 COSMolKit rows and is unranked. Its original 50ms exploratory attempt stopped before any complete block; the full 5ms replacement is retained.

The following table is a fresh post-change diagnostic with six balanced three-arm blocks per window, 5ms accumulated operation time per cell, and full output checks. Values are paired COSMolKit/source medians; values above 1 favor this source. These six-block diagnostics do not decide adoption or establish strict external victories. The [84-lane diagnostic CSV](2026-10-18-terminal-tpsa-screen-lanes.csv) and [numeric diagnostics](2026-10-18-terminal-tpsa-screen.json.gz) retain both epochs with their distinct source commits. Do not infer causal gains by pooling these two diagnostic epochs.

| Operation | First0 | First5000 | Reused0 | Reused5000 | Parse+op0 | Parse+op5000 |
|---|---:|---:|---:|---:|---:|---:|
| `tpsa` | 0.553 | 0.139 | 1.208 | 1.187 | 12.647 | 6.100 |
| `labute_asa` | 0.334 | 0.072 | 1.073 | 1.054 | 11.373 | 5.012 |
| `num_rings` | unranked | 0.084 | unranked | 1.016 | unranked | 13.190 |
| `morgan2_2048` | 2.467 | 1.052 | 3.218 | 3.001 | 8.557 | 4.637 |
| `morgan2_chiral` | 0.886 | 1.065 | 3.154 | 2.949 | 4.149 | 4.570 |
| `morgan3_2048` | 2.530 | 1.217 | 3.139 | 2.933 | 7.784 | 4.487 |
| `maccs` | 2.070 | 1.678 | 2.144 | 1.773 | 2.277 | 1.851 |

Default parsers prepare different amounts of work. COSMolKit uses sanitization/removal defaults; chematic keeps descriptor preparation lazy while seeding its model gate. First-use means the first operation on fresh default-parser objects, not empty global caches or cold process startup. Parse-inclusive and repeated-use remain separate.

The formal external timing scope is the three Morgan operations only; scalar/MACCS formal lanes compare source arms. Every strict external Morgan win additionally requires equivalent outputs, an unrounded interval lower bound above 1 and 24/24 faster blocks.
- `rdkit`: 18/18 strict Morgan wins in this source epoch.
- `cosmolkit`: 17/18 strict Morgan wins in this source epoch.

## Retained rejected experiments

| Candidate | Local source | Improvement gates passed | Source floor misses | Decision |
|---|---|---:|---:|---|
| Allocation/H-preparation | `09d9f70cdd23` | 1/6 | 0 | rejected; [all numbers](2026-10-18-cosmolkit-allocation-candidate-numeric.json.gz) |
| Valence/constant/local-TPSA combination | `36846e510457` | 3/6 | 0 | rejected; [all numbers](2026-10-18-cosmolkit-cold-terms-candidate-numeric.json.gz) |

Both earlier epochs remain rejected under their own six-primary rules. The final variant is genuinely different: remove those unconfirmed components, isolate terminal TPSA, and add a no-N/O shortcut before model preparation. Its two TPSA primaries and all 42 median floors were fixed locally before its first timing. No earlier elapsed samples are reused. Source archives, wheels, full canonical values, refusals and journals remain local; only numeric timing/accounting/digest packets are exported.

## Method and provenance

- Baseline: `f8e9d7dee5a8ece68a61bb9197f616bc86f6bfeb` (measurement baseline includes Cleanup PR #1034). Candidate: `add2b6b445c240b6c313927a98c551587239aff4`. Earlier rejected epochs use their explicitly recorded `bc36b398` baseline.
- Apple M4, 10 CPUs, 16GiB RAM; macOS 27.0.1 arm64; CPython 3.13.6. Pinned RDKit 2026.09.1 and COSMolKit 0.5.0. Same interpreter and native-inchi/opt3/full-LTO/one-CGU/one-job/no-incremental/no-strip settings; both final wheels built from empty targets. Loaded native hashes are checked against isolated installed wheels.
- Two disjoint 1,000-row windows at offsets 0/5000 of `rdkit-js-browser-10k-v1.smi`. Every permutation of baseline/candidate/RDKit/COSMolKit occurs once per window; each arm occupies each position six times. Warmup objects are disjoint from timed objects. Morgan cells accumulate 50ms; scalar/MACCS cells accumulate 20ms, including the TPSA adoption primaries.
- Reparse outside every timed full batch for first-use; prime and reuse one batch for repeated-use; put parsing inside the loop for parse-inclusive. Return native API objects and normalize outside each timed batch. Every full batch has raw duration/count/row/output-digest accounting. Per-batch means are summed time divided by actual batch count, not a contiguous steady-state run.
- Pointwise paired 95% bootstrap intervals use 10,000 fixed-seed resamples. Unrounded ratios decide adoption; no block, error, refusal or negative observation is removed. Full host-load observations and preregistration hashes are retained. Other host activity and this exposed development corpus limit transfer to other workloads, hosts and packages.
- Recomputed accounting: 192 journals, 2880 cells, 188696 timed batches and 188696000 processed rows. This row total counts repeated work, not independent chemical structures.
- Required benchmark CI recomputes the legacy and new packets, the selected TPSA primaries, all 42 floors, scalar accounting and rejected decisions, with corruption controls. Hosted CI validates the source; it does not reproduce this local timing experiment.
