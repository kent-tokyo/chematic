# Morgan: shared H inference and selective CIP refinement

**Combined source measurement gate: FAIL.** The gate and the separate CIP-only failure are computed from the original unrounded paired observations. This report does not establish published-package, current-main or toolkit-wide superiority.

The combined candidate computes the common implicit/bracket hydrogen count once for Morgan total-degree and total-H components, folds exactly the same invariant components in order without a temporary SmallVec, and rebuilds CIP comparison entries only in the union of actual resort sections. Shared boundaries are rebuilt once. Wildcards, explicit-H exclusion, integer casts, saturation, isotope/charge bits, optional ring suffix, rank updates and stereo/sanitization rules are preserved. No new molecule cache, public API or dependency.

## Fixed protocol and scope

- Baseline: source `a07a1294c550d5fa55e5fdac66c5755556d76ff2`, freshly built from an empty target before either epoch. Candidate wheels each use a new empty target. Identical Rust 1.97.0, opt-level=3, LTO=true, codegen-units=1, strip=none, native-inchi, one Cargo job, incremental off, and CPython 3.13.6 interpreter. Baseline ZIP/installed-native hashes are reverified in each epoch.
- One Apple M4 arm64 host (Mac16,12, 10 CPUs, 16 GiB), macOS 27.0.1; pinned installed RDKit 2026.09.1 and COSMolKit 0.5.0. Source wheels and an exposed development corpus; no cross-host, released-wheel or held-out generalization claim.
- Two unfiltered 1,000-row windows, offsets 0 and 5,000, of `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`. All 24 permutations of baseline/candidate/RDKit/COS once per profile, seeds 20261017+offset. Every arm occurs in every position six times.
- First use: parse fresh molecules outside each timed full batch, one operation per object. Repeated use: fresh objects parsed and primed outside timing, then reused. Parse inclusive: fresh parsing and operation both timed.
- Accumulate independently timed complete batches for at least 50 ms (Morgan) or 20 ms (controls); each actual duration, batch/row count and output hash retained. Preparation and normalization gaps are excluded; this is neither continuous steady state nor process cold start.
- Time public native results (chematic bytes versus comparator fingerprint objects); canonical bit conversion and output checking are outside timing. Validate every batch, including all iterations/modes. Whole-process RSS is not ranked.
- Source controls: TPSA, Labute ASA, ring count, MACCS. External timing covers R2 nonchiral, R2 chiral and R3 nonchiral Morgan only. Detailed/sparse/provenance APIs and other toolkit operations are outside these timing lanes.
- Adoption: both source first-use chiral Morgan paired 95% CI lower bounds >1; all 42 unrounded source medians >=0.99; complete experiments and exact source outputs including broad validation. No exception.
- Each external strict win requires equal outputs, paired 95% lower bound >1 and 24/24 faster blocks. Statistics: 10,000 paired bootstrap samples, fixed comparator-index seeds. Ratios are comparator/candidate, so >1 favors the candidate.
- Each fixed candidate has exactly one fresh formal epoch. The combined candidate was developed after the CIP-only epoch failed; no earlier elapsed samples are reused, pooled, or substituted. No own build, profiler, source edit or heavy audit runs during formal timing. The baseline is pinned despite independent main updates.

## Correctness and export boundary

CIP references contain the old invariant, segmentation and refinement bodies with only function/pool-name and formatting changes; all 14,072 rows compare in both seed modes and four invariant families, before/after supported sanitization. Shared-boundary, generated graph, isotope/map, empty/shrinking and thread-local cases pass. Combined candidate: 551 smiles tests passed (3 existing ignored), 388 fingerprint tests passed, all-target Clippy on both crates, fmt/refactor checks, and 329 selected Python binding tests. Atom-invariant references preserve the exact old construction and unchanged helpers; all corpus atoms and explicit-H/wildcard cast/saturation builder boundaries match.

In each epoch all 14,072 normal/stress rows × four R2/R3 chiral/nonchiral fingerprints = 56,288 cells equal baseline/candidate/RDKit exactly. COSMolKit normal 40,000 cells match; stress retains six refusals per operation and seven additional chiral mismatches per radius. Refusals and mismatches are not agreement.

Numeric exports retain every batch duration/count/output hash, all 78 lanes, source/wheel hashes and full-journal commitments, while excluding raw molecular outputs, binary payloads, local paths and SDK/process inventories. Complete raw outputs, source/runtime manifests, wheels and logs remain local. The [stdlib auditor](../scripts/check_morgan_numeric_packet.py) verifies numeric consistency, not the chemistry of unexported outputs.

## CIP sections only

Measured source `29cb4b5a889565b82820904e4530ac8ff0820b9e` versus `a07a1294c550d5fa55e5fdac66c5755556d76ff2`. [Numeric packet](2026-10-18-morgan-tied-cip-numeric.json.gz). Measurement gate: **FAIL**. 0 source median-floor violations; 2 source intervals have upper bound below one. 33/36 external strict wins.

Independent audit: 192 journals, 2,880 accounting cells, 187,269 timed full batches, 187,269,000 processed rows, no missing/partial journals. Source/runtime/wheel seals checked before and after both profiles.

| Offset | First-use chiral source/candidate | Paired 95% CI | Wins/24 |
|---:|---:|---|---:|
| 0 | 1.031754 | [1.024248, 1.035479] | 24 |
| 5000 | 1.002092 | [0.994040, 1.006608] | 13 |

All source conditions, including every negative interval:

| Offset | Operation | Mode | Source/candidate | Paired 95% CI | Wins/24 |
|---:|---|---|---:|---|---:|
| 0 | morgan2_2048 | first_use | 1.000050 | [0.992228, 1.005508] | 12 |
| 0 | morgan2_2048 | repeated_use | 1.000152 | [0.994014, 1.007913] | 12 |
| 0 | morgan2_2048 | parse_inclusive | 0.999684 | [0.997291, 1.004757] | 12 |
| 0 | morgan2_chiral | first_use | 1.031754 | [1.024248, 1.035479] | 24 |
| 0 | morgan2_chiral | repeated_use | 0.996300 | [0.988858, 1.011948] | 11 |
| 0 | morgan2_chiral | parse_inclusive | 1.025938 | [1.018438, 1.029873] | 24 |
| 0 | morgan3_2048 | first_use | 0.999137 | [0.994164, 1.006195] | 11 |
| 0 | morgan3_2048 | repeated_use | 0.997605 | [0.993623, 1.001720] | 8 |
| 0 | morgan3_2048 | parse_inclusive | 0.998295 | [0.991116, 0.999518] | 7 |
| 0 | tpsa | first_use | 1.001795 | [0.991315, 1.010519] | 14 |
| 0 | tpsa | repeated_use | 0.995486 | [0.992281, 1.003244] | 8 |
| 0 | tpsa | parse_inclusive | 0.991211 | [0.984338, 1.006685] | 9 |
| 0 | labute_asa | first_use | 0.999523 | [0.984326, 1.007706] | 12 |
| 0 | labute_asa | repeated_use | 0.993433 | [0.990741, 1.001582] | 8 |
| 0 | labute_asa | parse_inclusive | 1.005183 | [0.998641, 1.008742] | 15 |
| 0 | num_rings | first_use | 0.998145 | [0.989420, 1.003137] | 10 |
| 0 | num_rings | repeated_use | 0.995919 | [0.989747, 1.009059] | 11 |
| 0 | num_rings | parse_inclusive | 0.994758 | [0.986979, 1.007990] | 10 |
| 0 | maccs | first_use | 1.001544 | [0.998738, 1.004925] | 14 |
| 0 | maccs | repeated_use | 0.999939 | [0.996840, 1.002182] | 11 |
| 0 | maccs | parse_inclusive | 0.997653 | [0.995722, 0.999940] | 7 |
| 5000 | morgan2_2048 | first_use | 1.004222 | [0.998356, 1.010601] | 15 |
| 5000 | morgan2_2048 | repeated_use | 1.003813 | [0.986627, 1.008541] | 14 |
| 5000 | morgan2_2048 | parse_inclusive | 1.001707 | [0.994984, 1.007333] | 13 |
| 5000 | morgan2_chiral | first_use | 1.002092 | [0.994040, 1.006608] | 13 |
| 5000 | morgan2_chiral | repeated_use | 0.995341 | [0.981827, 1.004541] | 10 |
| 5000 | morgan2_chiral | parse_inclusive | 1.001606 | [0.995215, 1.013616] | 13 |
| 5000 | morgan3_2048 | first_use | 0.999518 | [0.995138, 1.005664] | 11 |
| 5000 | morgan3_2048 | repeated_use | 0.996935 | [0.985096, 1.002715] | 9 |
| 5000 | morgan3_2048 | parse_inclusive | 0.999561 | [0.990955, 1.003779] | 12 |
| 5000 | tpsa | first_use | 1.000053 | [0.993769, 1.015092] | 12 |
| 5000 | tpsa | repeated_use | 1.005396 | [1.000340, 1.013510] | 17 |
| 5000 | tpsa | parse_inclusive | 1.002663 | [0.994424, 1.011274] | 15 |
| 5000 | labute_asa | first_use | 0.994401 | [0.989464, 1.000623] | 8 |
| 5000 | labute_asa | repeated_use | 1.001889 | [0.995152, 1.009805] | 13 |
| 5000 | labute_asa | parse_inclusive | 1.011004 | [0.997878, 1.017904] | 16 |
| 5000 | num_rings | first_use | 0.991107 | [0.980487, 1.010822] | 10 |
| 5000 | num_rings | repeated_use | 0.993120 | [0.980628, 1.001407] | 9 |
| 5000 | num_rings | parse_inclusive | 1.008743 | [0.991259, 1.019713] | 13 |
| 5000 | maccs | first_use | 1.002235 | [1.000405, 1.005071] | 17 |
| 5000 | maccs | repeated_use | 1.000847 | [0.995696, 1.006910] | 14 |
| 5000 | maccs | parse_inclusive | 1.002505 | [0.999611, 1.011245] | 16 |

All external Morgan conditions:

| Offset | Comparator | Operation | Mode | Comparator/candidate | Paired 95% CI | Wins/24 | Strict win |
|---:|---|---|---|---:|---|---:|---|
| 0 | rdkit | morgan2_2048 | first_use | 3.364670 | [3.343272, 3.374314] | 24 | yes |
| 0 | cosmolkit | morgan2_2048 | first_use | 2.114647 | [2.098088, 2.128121] | 24 | yes |
| 0 | rdkit | morgan2_2048 | repeated_use | 5.462069 | [5.435625, 5.497309] | 24 | yes |
| 0 | cosmolkit | morgan2_2048 | repeated_use | 3.286080 | [3.260520, 3.319438] | 24 | yes |
| 0 | rdkit | morgan2_2048 | parse_inclusive | 8.713808 | [8.652259, 8.751984] | 24 | yes |
| 0 | cosmolkit | morgan2_2048 | parse_inclusive | 7.646167 | [7.620705, 7.679572] | 24 | yes |
| 0 | rdkit | morgan2_chiral | first_use | 1.338561 | [1.331377, 1.350822] | 24 | yes |
| 0 | cosmolkit | morgan2_chiral | first_use | 0.850357 | [0.848815, 0.853271] | 0 | no |
| 0 | rdkit | morgan2_chiral | repeated_use | 5.321590 | [5.275896, 5.388758] | 24 | yes |
| 0 | cosmolkit | morgan2_chiral | repeated_use | 3.245106 | [3.224543, 3.288568] | 24 | yes |
| 0 | rdkit | morgan2_chiral | parse_inclusive | 4.564693 | [4.544256, 4.604922] | 24 | yes |
| 0 | cosmolkit | morgan2_chiral | parse_inclusive | 4.024064 | [3.996171, 4.060343] | 24 | yes |
| 0 | rdkit | morgan3_2048 | first_use | 3.467163 | [3.435446, 3.488547] | 24 | yes |
| 0 | cosmolkit | morgan3_2048 | first_use | 2.274070 | [2.250229, 2.283815] | 24 | yes |
| 0 | rdkit | morgan3_2048 | repeated_use | 5.062427 | [5.003820, 5.127542] | 24 | yes |
| 0 | cosmolkit | morgan3_2048 | repeated_use | 3.200693 | [3.186863, 3.223705] | 24 | yes |
| 0 | rdkit | morgan3_2048 | parse_inclusive | 8.123868 | [8.073395, 8.191414] | 24 | yes |
| 0 | cosmolkit | morgan3_2048 | parse_inclusive | 7.067004 | [7.022883, 7.115927] | 24 | yes |
| 5000 | rdkit | morgan2_2048 | first_use | 1.609183 | [1.599930, 1.617517] | 24 | yes |
| 5000 | cosmolkit | morgan2_2048 | first_use | 1.009561 | [0.994775, 1.025772] | 14 | no |
| 5000 | rdkit | morgan2_2048 | repeated_use | 5.133077 | [5.089667, 5.160692] | 24 | yes |
| 5000 | cosmolkit | morgan2_2048 | repeated_use | 3.047177 | [3.023400, 3.085569] | 24 | yes |
| 5000 | rdkit | morgan2_2048 | parse_inclusive | 4.645223 | [4.613916, 4.674062] | 24 | yes |
| 5000 | cosmolkit | morgan2_2048 | parse_inclusive | 4.475485 | [4.420905, 4.521482] | 24 | yes |
| 5000 | rdkit | morgan2_chiral | first_use | 1.634777 | [1.613851, 1.642829] | 24 | yes |
| 5000 | cosmolkit | morgan2_chiral | first_use | 1.029894 | [1.016617, 1.042746] | 23 | no |
| 5000 | rdkit | morgan2_chiral | repeated_use | 5.066309 | [4.984046, 5.118565] | 24 | yes |
| 5000 | cosmolkit | morgan2_chiral | repeated_use | 3.009972 | [2.980146, 3.037048] | 24 | yes |
| 5000 | rdkit | morgan2_chiral | parse_inclusive | 4.631947 | [4.576664, 4.658979] | 24 | yes |
| 5000 | cosmolkit | morgan2_chiral | parse_inclusive | 4.468147 | [4.432554, 4.487789] | 24 | yes |
| 5000 | rdkit | morgan3_2048 | first_use | 1.828469 | [1.810455, 1.844697] | 24 | yes |
| 5000 | cosmolkit | morgan3_2048 | first_use | 1.183517 | [1.172013, 1.190569] | 24 | yes |
| 5000 | rdkit | morgan3_2048 | repeated_use | 4.830204 | [4.776713, 4.867810] | 24 | yes |
| 5000 | cosmolkit | morgan3_2048 | repeated_use | 3.016904 | [2.994644, 3.038577] | 24 | yes |
| 5000 | rdkit | morgan3_2048 | parse_inclusive | 4.566046 | [4.546741, 4.586555] | 24 | yes |
| 5000 | cosmolkit | morgan3_2048 | parse_inclusive | 4.341548 | [4.307803, 4.379497] | 24 | yes |

## CIP sections + shared H inference

Measured source `e781de72a721e1903e8fad455fe253acacf66c0a` versus `a07a1294c550d5fa55e5fdac66c5755556d76ff2`. [Numeric packet](2026-10-18-morgan-shared-h-numeric.json.gz). Measurement gate: **FAIL**. 2 source median-floor violations; 0 source intervals have upper bound below one. 33/36 external strict wins.

Independent audit: 192 journals, 2,880 accounting cells, 152,817 timed full batches, 152,817,000 processed rows, no missing/partial journals. Source/runtime/wheel seals checked before and after both profiles.

| Offset | First-use chiral source/candidate | Paired 95% CI | Wins/24 |
|---:|---:|---|---:|
| 0 | 1.055363 | [1.048913, 1.061574] | 24 |
| 5000 | 1.026339 | [1.019458, 1.040012] | 21 |

All source conditions, including every negative interval:

| Offset | Operation | Mode | Source/candidate | Paired 95% CI | Wins/24 |
|---:|---|---|---:|---|---:|
| 0 | morgan2_2048 | first_use | 1.068162 | [1.057896, 1.083483] | 24 |
| 0 | morgan2_2048 | repeated_use | 1.137044 | [1.110064, 1.154848] | 24 |
| 0 | morgan2_2048 | parse_inclusive | 1.039951 | [1.022056, 1.061122] | 22 |
| 0 | morgan2_chiral | first_use | 1.055363 | [1.048913, 1.061574] | 24 |
| 0 | morgan2_chiral | repeated_use | 1.102650 | [1.083078, 1.118365] | 24 |
| 0 | morgan2_chiral | parse_inclusive | 1.035132 | [1.020935, 1.048552] | 22 |
| 0 | morgan3_2048 | first_use | 1.056115 | [1.040825, 1.061039] | 24 |
| 0 | morgan3_2048 | repeated_use | 1.085237 | [1.068746, 1.091245] | 23 |
| 0 | morgan3_2048 | parse_inclusive | 1.027236 | [1.018505, 1.040334] | 21 |
| 0 | tpsa | first_use | 0.991812 | [0.975509, 1.006836] | 9 |
| 0 | tpsa | repeated_use | 1.000600 | [0.988515, 1.018665] | 12 |
| 0 | tpsa | parse_inclusive | 0.999263 | [0.987165, 1.015217] | 12 |
| 0 | labute_asa | first_use | 0.995740 | [0.972391, 1.012618] | 12 |
| 0 | labute_asa | repeated_use | 0.997863 | [0.987384, 1.009967] | 11 |
| 0 | labute_asa | parse_inclusive | 1.005486 | [0.981389, 1.015461] | 13 |
| 0 | num_rings | first_use | 1.001133 | [0.983466, 1.014161] | 12 |
| 0 | num_rings | repeated_use | 1.002992 | [0.994982, 1.017196] | 14 |
| 0 | num_rings | parse_inclusive | 0.999239 | [0.978344, 1.018367] | 12 |
| 0 | maccs | first_use | 0.998743 | [0.991114, 1.005372] | 10 |
| 0 | maccs | repeated_use | 1.001375 | [0.992991, 1.007629] | 13 |
| 0 | maccs | parse_inclusive | 0.997747 | [0.992362, 1.003377] | 10 |
| 5000 | morgan2_2048 | first_use | 1.029939 | [1.024146, 1.044361] | 21 |
| 5000 | morgan2_2048 | repeated_use | 1.110650 | [1.096063, 1.140057] | 23 |
| 5000 | morgan2_2048 | parse_inclusive | 1.020511 | [1.006740, 1.035372] | 19 |
| 5000 | morgan2_chiral | first_use | 1.026339 | [1.019458, 1.040012] | 21 |
| 5000 | morgan2_chiral | repeated_use | 1.117238 | [1.089161, 1.137584] | 23 |
| 5000 | morgan2_chiral | parse_inclusive | 1.018314 | [1.008043, 1.036923] | 21 |
| 5000 | morgan3_2048 | first_use | 1.026283 | [1.014109, 1.038029] | 19 |
| 5000 | morgan3_2048 | repeated_use | 1.094057 | [1.075959, 1.107871] | 23 |
| 5000 | morgan3_2048 | parse_inclusive | 1.022613 | [1.013690, 1.033748] | 19 |
| 5000 | tpsa | first_use | 0.989602 | [0.977494, 1.009787] | 9 |
| 5000 | tpsa | repeated_use | 1.004680 | [0.998597, 1.014802] | 16 |
| 5000 | tpsa | parse_inclusive | 0.998152 | [0.982563, 1.016388] | 11 |
| 5000 | labute_asa | first_use | 0.996258 | [0.990618, 1.009112] | 10 |
| 5000 | labute_asa | repeated_use | 1.005552 | [0.993263, 1.019732] | 14 |
| 5000 | labute_asa | parse_inclusive | 0.989197 | [0.981064, 1.018999] | 10 |
| 5000 | num_rings | first_use | 0.992488 | [0.976606, 1.011771] | 10 |
| 5000 | num_rings | repeated_use | 1.000629 | [0.992847, 1.007506] | 13 |
| 5000 | num_rings | parse_inclusive | 1.013973 | [0.980412, 1.036044] | 13 |
| 5000 | maccs | first_use | 0.996316 | [0.990694, 1.005642] | 9 |
| 5000 | maccs | repeated_use | 0.999300 | [0.996437, 1.004987] | 11 |
| 5000 | maccs | parse_inclusive | 1.002249 | [0.994677, 1.009516] | 13 |

All external Morgan conditions:

| Offset | Comparator | Operation | Mode | Comparator/candidate | Paired 95% CI | Wins/24 | Strict win |
|---:|---|---|---|---:|---|---:|---|
| 0 | rdkit | morgan2_2048 | first_use | 3.469446 | [3.380638, 3.505783] | 24 | yes |
| 0 | cosmolkit | morgan2_2048 | first_use | 2.208196 | [2.156906, 2.222550] | 24 | yes |
| 0 | rdkit | morgan2_2048 | repeated_use | 5.929852 | [5.776862, 6.033222] | 24 | yes |
| 0 | cosmolkit | morgan2_2048 | repeated_use | 3.634059 | [3.540477, 3.712033] | 24 | yes |
| 0 | rdkit | morgan2_2048 | parse_inclusive | 8.790540 | [8.696836, 8.880971] | 24 | yes |
| 0 | cosmolkit | morgan2_2048 | parse_inclusive | 7.762934 | [7.698606, 7.875700] | 24 | yes |
| 0 | rdkit | morgan2_chiral | first_use | 1.339730 | [1.319672, 1.376986] | 24 | yes |
| 0 | cosmolkit | morgan2_chiral | first_use | 0.858715 | [0.849144, 0.868872] | 1 | no |
| 0 | rdkit | morgan2_chiral | repeated_use | 5.689559 | [5.594343, 5.806008] | 24 | yes |
| 0 | cosmolkit | morgan2_chiral | repeated_use | 3.549081 | [3.481125, 3.580866] | 24 | yes |
| 0 | rdkit | morgan2_chiral | parse_inclusive | 4.546589 | [4.448359, 4.589995] | 24 | yes |
| 0 | cosmolkit | morgan2_chiral | parse_inclusive | 4.006301 | [3.962171, 4.052323] | 24 | yes |
| 0 | rdkit | morgan3_2048 | first_use | 3.516608 | [3.492314, 3.559742] | 24 | yes |
| 0 | cosmolkit | morgan3_2048 | first_use | 2.283418 | [2.261768, 2.297213] | 24 | yes |
| 0 | rdkit | morgan3_2048 | repeated_use | 5.188778 | [5.148372, 5.246383] | 24 | yes |
| 0 | cosmolkit | morgan3_2048 | repeated_use | 3.290945 | [3.243166, 3.418495] | 24 | yes |
| 0 | rdkit | morgan3_2048 | parse_inclusive | 8.101123 | [8.025969, 8.225364] | 24 | yes |
| 0 | cosmolkit | morgan3_2048 | parse_inclusive | 7.008550 | [6.889190, 7.140982] | 24 | yes |
| 5000 | rdkit | morgan2_2048 | first_use | 1.626349 | [1.614458, 1.647254] | 24 | yes |
| 5000 | cosmolkit | morgan2_2048 | first_use | 1.048297 | [1.032603, 1.056836] | 23 | no |
| 5000 | rdkit | morgan2_2048 | repeated_use | 5.556805 | [5.454833, 5.666898] | 24 | yes |
| 5000 | cosmolkit | morgan2_2048 | repeated_use | 3.375473 | [3.327521, 3.406491] | 24 | yes |
| 5000 | rdkit | morgan2_2048 | parse_inclusive | 4.646813 | [4.595532, 4.733317] | 24 | yes |
| 5000 | cosmolkit | morgan2_2048 | parse_inclusive | 4.502740 | [4.475135, 4.551485] | 24 | yes |
| 5000 | rdkit | morgan2_chiral | first_use | 1.645758 | [1.628285, 1.686155] | 24 | yes |
| 5000 | cosmolkit | morgan2_chiral | first_use | 1.066497 | [1.055195, 1.071910] | 23 | no |
| 5000 | rdkit | morgan2_chiral | repeated_use | 5.471632 | [5.315211, 5.567605] | 24 | yes |
| 5000 | cosmolkit | morgan2_chiral | repeated_use | 3.299007 | [3.258501, 3.361127] | 24 | yes |
| 5000 | rdkit | morgan2_chiral | parse_inclusive | 4.652863 | [4.553357, 4.688322] | 24 | yes |
| 5000 | cosmolkit | morgan2_chiral | parse_inclusive | 4.517016 | [4.492727, 4.556373] | 24 | yes |
| 5000 | rdkit | morgan3_2048 | first_use | 1.827701 | [1.819502, 1.841624] | 24 | yes |
| 5000 | cosmolkit | morgan3_2048 | first_use | 1.217254 | [1.202243, 1.221859] | 24 | yes |
| 5000 | rdkit | morgan3_2048 | repeated_use | 5.065410 | [5.003689, 5.179853] | 24 | yes |
| 5000 | cosmolkit | morgan3_2048 | repeated_use | 3.234347 | [3.179274, 3.275776] | 24 | yes |
| 5000 | rdkit | morgan3_2048 | parse_inclusive | 4.552718 | [4.523071, 4.573430] | 24 | yes |
| 5000 | cosmolkit | morgan3_2048 | parse_inclusive | 4.363944 | [4.325509, 4.381237] | 24 | yes |

## Recompute numeric consistency

```sh
python3 scripts/check_morgan_numeric_packet.py benchmarks/2026-10-18-morgan-tied-cip-numeric.json.gz
python3 scripts/check_morgan_numeric_packet.py benchmarks/2026-10-18-morgan-shared-h-numeric.json.gz
```

The scalar and MACCS controls do not establish external descriptor wins. The chemistry gate is distinct from speed, and the exact measured source remains distinct from later main integration and hosted/public artifacts.

## Integration status

Both fixed runtime candidates missed the preregistered adoption gate and are excluded from this evidence-only change. The combined [candidate source](https://github.com/kent-tokyo/chematic/commit/e781de72a721e1903e8fad455fe253acacf66c0a) is retained on its own branch for inspection; the [source diff](https://github.com/kent-tokyo/chematic/compare/a07a1294c550d5fa55e5fdac66c5755556d76ff2...e781de72a721e1903e8fad455fe253acacf66c0a) contains no private output/environment payload.

The combined primary improvements are established within their intervals, but TPSA first use and Labute ASA parse inclusive at offset 5,000 have paired medians 0.989602 and 0.989197, below the fixed 0.99 floor. Their intervals cross one; this does not erase the median-budget failure. COSMolKit remains faster in offset-0 first-use chiral Morgan (COS/candidate 0.858715, 1/24 faster candidate blocks). No failed gate is relabelled as passed.

Evidence branch base: `0e7a6d9f8d50fd75a9b1bfada00336f8e320f8ec`. No crate source changes are proposed. Main updates and published artifacts are distinct from the historical measured revisions above.
