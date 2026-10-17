# Legacy stereo code projection: rejected source experiment — 2026-10-18

**Adoption gate failed; runtime implementation excluded.** The first-use chiral Morgan primary improves by **1.044867464×** in the offset-0 window but regresses to **0.992939314×** in offset-5,000. Both source confidence intervals must have lower bound above one. The second does not, so the candidate is rejected even though every one of the 42 source medians meets the separate 0.99 floor. Seven source intervals have upper bound below one and remain visible below.

Measured baseline: `01b874e4fd5f0f59d0a3ff565c52f5cc7f1de384`. Measured rejected candidate: [`4f5e14678c7d0b6c09333877ad8f96dc30010362`](https://github.com/kent-tokyo/chematic/commit/4f5e14678c7d0b6c09333877ad8f96dc30010362) ([source diff](https://github.com/kent-tokyo/chematic/compare/01b874e4fd5f0f59d0a3ff565c52f5cc7f1de384...4f5e14678c7d0b6c09333877ad8f96dc30010362)), retained on `codex/morgan-stereo-projection-candidate-20261018`. The evidence branch contains the report, numeric packet, index link and numeric-corruption test extension only. No part of this failed runtime candidate is promoted. Main advanced during measurement (report-preparation snapshot `229881de8d1cc2639874eba3c5d416d2ceccc469`, including separate aromaticity edits); this is a pinned historical source comparison, not a benchmark of that newer main or published packages.

Ratios mean comparator time / candidate time, computed per paired block and summarized from those ratios. Above one favors the candidate. Unrounded medians and bootstrap intervals determine adoption, not ratios of aggregate medians.

| Offset | Operation | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks | Strict external win |
| ---: | --- | --- | --- | ---: | --- | ---: | --- |
| 0 | morgan2_chiral | first_use | baseline | 1.044867464 | [1.042886612, 1.050407089] | 24/24 | — |
| 5000 | morgan2_chiral | first_use | baseline | 0.992939314 | [0.988321356, 0.998530433] | 6/24 | — |

The measured Morgan lanes meet the strict external win criterion in **18/18 RDKit** and **15/18 COSMolKit** comparisons. That is bounded to the three selected Morgan configurations, two windows and three usage modes. Complete toolkit superiority is not established. All three COSMolKit non-wins follow; an interval above one with 21/24 faster blocks still fails the predeclared 24/24 requirement. The slow offset-0 chiral first-use lane is 0.863248328× (candidate latency roughly 16% greater).

| Offset | Operation | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks | Strict external win |
| ---: | --- | --- | --- | ---: | --- | ---: | --- |
| 0 | morgan2_chiral | first_use | cosmolkit | 0.863248328 | [0.860401706, 0.867222351] | 0/24 | no |
| 5000 | morgan2_2048 | first_use | cosmolkit | 0.994213707 | [0.990517626, 1.003244410] | 9/24 | no |
| 5000 | morgan2_chiral | first_use | cosmolkit | 1.013368568 | [1.009409602, 1.022563447] | 21/24 | no |

The candidate specializes legacy stereo processing for consumers that read CIP codes, atom tags and bond stereo, while discarding `_ChiralityPossible` and final CIP ranks. Untagged atoms still initialize ranks at the same point and contribute the same legal-center/unassigned counts; only their unused flag writes and duplicate-neighbor work are omitted. When every tagged atom already has a CIP code and no unresolved stereo bond remains, a final rerank/flag-only iteration has no consumer, including the ring-special-case cleanup, and is omitted. Any unresolved tagged atom or stereo bond retains the ranking path. APIs that return potential-center flags or ranks keep full behavior.

A reused neighbor Vec, four-rank stack storage and inline permutation indices remove local allocations; the permutation buffer retains a heap fallback for high physical degree. Outgoing dative bonds still enter neighbor order but are excluded before bounded rank storage, exactly as before. Ring-case cleanup is factored into a small private helper to respect the existing refactor ratchet. CIP invariant/refinement functions, bond-stereo assignment and reranking bodies are byte-identical to the baseline; no graph representation, molecule cache, Morgan invariant, public API or dependency changes. The earlier failed CIP-pool, tied-section, shared-H and inline-adjacency deltas are absent.

Correctness before timing: **950 Rust tests passed, 0 failed, 3 existing ignored, 40 suites**; all-target Clippy on both crates, fmt, refactor and normal commit hook pass. Frozen old potential-center, atom-code, legacy-perception and public API bodies match after only names and Rustfmt changes. All **14,072 input rows** match the old public stereo API, including errors, and all **14,072 raw + 14,072 sanitized model views** match the full old internal state after excluding only the intentionally unused possible-center flag. Three rank families check neighbor append order and duplicate stopping. Another **384 generated cases** cover high degree, both dative orientations, bracket H, six elements, and full clean/flag modes; empty input also matches. Seven CIP/bond function bodies are independently checked unchanged.

Selected Python binding tests: **329 passed**. Across the normal 10k and 4,072-row stress corpora, all **56,288 R2/R3 chiral/nonchiral fingerprint cells** match baseline/candidate/RDKit exactly. COSMolKit matches all 40,000 normal cells; stress retains six refusals per operation and seven additional chiral mismatches per radius. Refusals/mismatches are retained and do not count as agreement. On the two formal windows all four arms' Morgan outputs agree; source TPSA/ASA/ring/MACCS outputs also agree in every timed batch. This evidence PR's hosted checks concern integration, separately from measured-source correctness.

The plan is fixed before candidate timing. During implementation, before any executable was timed, a recorded plan amendment added the proof-limited unused-final-rerank elimination; the earlier plan is retained. Each source arm has its own empty-target native-inchi wheel build with Rust 1.97.0, opt3/full LTO/one CGU/one job, no incremental compilation or stripping, and CPython 3.13.6. Installed native hashes match each wheel. RDKit 2026.09.1 and COSMolKit 0.5.0 are the pinned, reverified installed comparators on one Apple M4 arm64 macOS 27.0.1 host (10 CPUs, 16GiB RAM). No earlier elapsed observations are reused or pooled.

Two unfiltered 1,000-row windows at offsets 0 and 5,000 use all 24 permutations of baseline/candidate/RDKit/COSMolKit exactly once, shuffled with seeds 20261017+offset. Every arm appears in every position six times. First use parses fresh objects outside each timed full batch and performs one operation per object. Repeated use primes fresh objects outside timing then reuses them. Parse inclusive times both fresh parse and operation. Each operation accumulates independent full batches for at least 50ms (Morgan) or 20ms (controls); preparation, canonical bit conversion, and every-batch output normalization/checking are outside timed operation intervals. Public native results are timed (chematic bytes and comparator fingerprint objects). This is not continuous steady state or a process cold start; whole-process RSS is unranked.

Adoption requires both first-use chiral source intervals' lower bounds>1, all 42 source unrounded medians>=0.99, complete experiments and exact source outputs including broad checks. Each external strict win additionally requires equal outputs and 24/24 faster paired blocks. Statistics use 10,000 paired bootstrap samples with fixed comparator-index seeds. External descriptors, other fingerprint configurations, detailed/sparse/provenance APIs and other toolkit operations have no new external speed ranking.

The final formal epoch completed all **48 blocks / 192 arm journals / 2,880 accounting cells / 189,469 timed full batches / 189,469,000 processed rows**, without interruption, filtering or reruns. Before the epoch began, capacity checks blocked only operational preparation with zero formal samples; verified unused compilation outputs were cleaned while retaining source/wheels/binaries/logs/evidence. No own build, source edit, profiler or heavy audit ran during formal timing. Source, runtime, wrapper, wheel and input seals were verified before and after both profiles. This shared host has no guarantee of global quiescence; the observations do not isolate the cause of sub-percent regressions. Original dirty documentation and unrelated worktrees were preserved.

The [2026-10-18-morgan-stereo-projection-numeric.json.gz](2026-10-18-morgan-stereo-projection-numeric.json.gz) packet retains every batch duration, actual count, output hash, all 78 lanes, wheel/native hashes and raw-journal commitments. It exports no compound streams, SDK/process/environment inventories, source archives, binaries or absolute local paths. Full raw journals, molecular streams, build/runtime manifests, both wheels and source archives remain local. Independent full-data and stdlib numeric auditors agree the gate failed. The numeric audit verifies consistency, not the chemistry of unexported molecular streams. Corrupt elapsed sums, missing arms and changed output digests are tested separately.

```sh
python3 scripts/check_morgan_numeric_packet.py benchmarks/2026-10-18-morgan-stereo-projection-numeric.json.gz
```

All seven source intervals wholly below one:

| Offset | Operation | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks | Strict external win |
| ---: | --- | --- | --- | ---: | --- | ---: | --- |
| 0 | maccs | first_use | baseline | 0.996071966 | [0.994279931, 0.998370414] | 5/24 | — |
| 0 | maccs | repeated_use | baseline | 0.994447071 | [0.989885953, 0.997934976] | 6/24 | — |
| 0 | maccs | parse_inclusive | baseline | 0.995658132 | [0.993053915, 0.999972233] | 7/24 | — |
| 5000 | morgan2_2048 | parse_inclusive | baseline | 0.994682854 | [0.990599723, 0.997909069] | 6/24 | — |
| 5000 | morgan2_chiral | first_use | baseline | 0.992939314 | [0.988321356, 0.998530433] | 6/24 | — |
| 5000 | maccs | first_use | baseline | 0.994801024 | [0.992332936, 0.999126640] | 7/24 | — |
| 5000 | maccs | parse_inclusive | baseline | 0.994960988 | [0.991828282, 0.998483291] | 5/24 | — |

All 42 source comparisons, including every negative:

| Offset | Operation | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks | Strict external win |
| ---: | --- | --- | --- | ---: | --- | ---: | --- |
| 0 | morgan2_2048 | first_use | baseline | 1.005217266 | [1.000846517, 1.006543015] | 18/24 | — |
| 0 | morgan2_2048 | repeated_use | baseline | 1.011848186 | [1.001252115, 1.016803005] | 18/24 | — |
| 0 | morgan2_2048 | parse_inclusive | baseline | 1.002516665 | [0.998497727, 1.009787382] | 14/24 | — |
| 0 | morgan2_chiral | first_use | baseline | 1.044867464 | [1.042886612, 1.050407089] | 24/24 | — |
| 0 | morgan2_chiral | repeated_use | baseline | 0.998355785 | [0.992865989, 1.007069276] | 11/24 | — |
| 0 | morgan2_chiral | parse_inclusive | baseline | 1.035196917 | [1.033813800, 1.040012149] | 24/24 | — |
| 0 | morgan3_2048 | first_use | baseline | 1.002788312 | [0.999467711, 1.006802907] | 15/24 | — |
| 0 | morgan3_2048 | repeated_use | baseline | 1.003548596 | [0.995442060, 1.007312082] | 14/24 | — |
| 0 | morgan3_2048 | parse_inclusive | baseline | 1.000669099 | [0.996259327, 1.004248550] | 14/24 | — |
| 0 | tpsa | first_use | baseline | 1.010950533 | [0.995364795, 1.017343493] | 14/24 | — |
| 0 | tpsa | repeated_use | baseline | 1.001950075 | [0.994537212, 1.005858049] | 13/24 | — |
| 0 | tpsa | parse_inclusive | baseline | 1.004048811 | [0.993561218, 1.010062442] | 13/24 | — |
| 0 | labute_asa | first_use | baseline | 1.007614529 | [0.997667040, 1.012026584] | 16/24 | — |
| 0 | labute_asa | repeated_use | baseline | 1.005324708 | [1.002796452, 1.009032038] | 20/24 | — |
| 0 | labute_asa | parse_inclusive | baseline | 0.997756302 | [0.986159885, 1.012800872] | 11/24 | — |
| 0 | num_rings | first_use | baseline | 1.001912930 | [0.994342189, 1.018432632] | 12/24 | — |
| 0 | num_rings | repeated_use | baseline | 1.001400892 | [0.998339943, 1.007577169] | 16/24 | — |
| 0 | num_rings | parse_inclusive | baseline | 1.001913327 | [0.989092778, 1.013270871] | 13/24 | — |
| 0 | maccs | first_use | baseline | 0.996071966 | [0.994279931, 0.998370414] | 5/24 | — |
| 0 | maccs | repeated_use | baseline | 0.994447071 | [0.989885953, 0.997934976] | 6/24 | — |
| 0 | maccs | parse_inclusive | baseline | 0.995658132 | [0.993053915, 0.999972233] | 7/24 | — |
| 5000 | morgan2_2048 | first_use | baseline | 0.999568038 | [0.992156527, 1.004808727] | 11/24 | — |
| 5000 | morgan2_2048 | repeated_use | baseline | 0.999503061 | [0.994046080, 1.007005064] | 11/24 | — |
| 5000 | morgan2_2048 | parse_inclusive | baseline | 0.994682854 | [0.990599723, 0.997909069] | 6/24 | — |
| 5000 | morgan2_chiral | first_use | baseline | 0.992939314 | [0.988321356, 0.998530433] | 6/24 | — |
| 5000 | morgan2_chiral | repeated_use | baseline | 1.005744680 | [0.999368513, 1.012620015] | 16/24 | — |
| 5000 | morgan2_chiral | parse_inclusive | baseline | 0.998073088 | [0.994338076, 1.004582222] | 10/24 | — |
| 5000 | morgan3_2048 | first_use | baseline | 0.999233334 | [0.996165965, 1.006474115] | 11/24 | — |
| 5000 | morgan3_2048 | repeated_use | baseline | 1.008958819 | [1.000983791, 1.014283799] | 17/24 | — |
| 5000 | morgan3_2048 | parse_inclusive | baseline | 1.001350938 | [0.998284583, 1.004838716] | 17/24 | — |
| 5000 | tpsa | first_use | baseline | 0.994613924 | [0.983043616, 1.005527755] | 9/24 | — |
| 5000 | tpsa | repeated_use | baseline | 1.001461016 | [0.995987907, 1.003965687] | 13/24 | — |
| 5000 | tpsa | parse_inclusive | baseline | 1.007221245 | [0.985740297, 1.011527652] | 14/24 | — |
| 5000 | labute_asa | first_use | baseline | 0.998555764 | [0.988211089, 1.004966599] | 12/24 | — |
| 5000 | labute_asa | repeated_use | baseline | 1.006653039 | [0.997547361, 1.010132834] | 15/24 | — |
| 5000 | labute_asa | parse_inclusive | baseline | 1.002166614 | [0.992552177, 1.009734177] | 14/24 | — |
| 5000 | num_rings | first_use | baseline | 0.999678248 | [0.988348538, 1.010635805] | 12/24 | — |
| 5000 | num_rings | repeated_use | baseline | 1.001740535 | [0.994600298, 1.012756932] | 13/24 | — |
| 5000 | num_rings | parse_inclusive | baseline | 0.997552058 | [0.987931695, 1.005752348] | 11/24 | — |
| 5000 | maccs | first_use | baseline | 0.994801024 | [0.992332936, 0.999126640] | 7/24 | — |
| 5000 | maccs | repeated_use | baseline | 0.996413499 | [0.991575526, 1.000094205] | 8/24 | — |
| 5000 | maccs | parse_inclusive | baseline | 0.994960988 | [0.991828282, 0.998483291] | 5/24 | — |

All 36 external Morgan comparisons:

| Offset | Operation | Mode | Comparator | Paired median ratio | Paired 95% CI | Faster blocks | Strict external win |
| ---: | --- | --- | --- | ---: | --- | ---: | --- |
| 0 | morgan2_2048 | first_use | rdkit | 3.323045419 | [3.311872701, 3.332223952] | 24/24 | yes |
| 0 | morgan2_2048 | first_use | cosmolkit | 2.106733600 | [2.093281286, 2.113145349] | 24/24 | yes |
| 0 | morgan2_2048 | repeated_use | rdkit | 5.422971191 | [5.395041062, 5.454817573] | 24/24 | yes |
| 0 | morgan2_2048 | repeated_use | cosmolkit | 3.275026421 | [3.269679965, 3.295594665] | 24/24 | yes |
| 0 | morgan2_2048 | parse_inclusive | rdkit | 8.698756761 | [8.645779776, 8.771903350] | 24/24 | yes |
| 0 | morgan2_2048 | parse_inclusive | cosmolkit | 7.657012641 | [7.621190901, 7.720097686] | 24/24 | yes |
| 0 | morgan2_chiral | first_use | rdkit | 1.353153319 | [1.343643340, 1.358430975] | 24/24 | yes |
| 0 | morgan2_chiral | first_use | cosmolkit | 0.863248328 | [0.860401706, 0.867222351] | 0/24 | no |
| 0 | morgan2_chiral | repeated_use | rdkit | 5.330759386 | [5.253745265, 5.357776079] | 24/24 | yes |
| 0 | morgan2_chiral | repeated_use | cosmolkit | 3.226450735 | [3.207221408, 3.264768036] | 24/24 | yes |
| 0 | morgan2_chiral | parse_inclusive | rdkit | 4.650026852 | [4.637948903, 4.658645251] | 24/24 | yes |
| 0 | morgan2_chiral | parse_inclusive | cosmolkit | 4.086913719 | [4.074529563, 4.108705612] | 24/24 | yes |
| 0 | morgan3_2048 | first_use | rdkit | 3.457759072 | [3.425111117, 3.473671903] | 24/24 | yes |
| 0 | morgan3_2048 | first_use | cosmolkit | 2.251969459 | [2.233001195, 2.262060883] | 24/24 | yes |
| 0 | morgan3_2048 | repeated_use | rdkit | 5.023076352 | [4.971563675, 5.078413236] | 24/24 | yes |
| 0 | morgan3_2048 | repeated_use | cosmolkit | 3.174650876 | [3.156519884, 3.198959904] | 24/24 | yes |
| 0 | morgan3_2048 | parse_inclusive | rdkit | 8.125570995 | [8.072129083, 8.160860950] | 24/24 | yes |
| 0 | morgan3_2048 | parse_inclusive | cosmolkit | 7.066427376 | [7.040897783, 7.107480540] | 24/24 | yes |
| 5000 | morgan2_2048 | first_use | rdkit | 1.593558020 | [1.585819766, 1.606191684] | 24/24 | yes |
| 5000 | morgan2_2048 | first_use | cosmolkit | 0.994213707 | [0.990517626, 1.003244410] | 9/24 | no |
| 5000 | morgan2_2048 | repeated_use | rdkit | 5.092562577 | [5.051157359, 5.134217371] | 24/24 | yes |
| 5000 | morgan2_2048 | repeated_use | cosmolkit | 3.012849857 | [2.994931559, 3.042233625] | 24/24 | yes |
| 5000 | morgan2_2048 | parse_inclusive | rdkit | 4.625421231 | [4.610290413, 4.656436915] | 24/24 | yes |
| 5000 | morgan2_2048 | parse_inclusive | cosmolkit | 4.465861473 | [4.440712817, 4.494541270] | 24/24 | yes |
| 5000 | morgan2_chiral | first_use | rdkit | 1.610158637 | [1.597577485, 1.616816002] | 24/24 | yes |
| 5000 | morgan2_chiral | first_use | cosmolkit | 1.013368568 | [1.009409602, 1.022563447] | 21/24 | no |
| 5000 | morgan2_chiral | repeated_use | rdkit | 5.060807042 | [5.023167011, 5.106885976] | 24/24 | yes |
| 5000 | morgan2_chiral | repeated_use | cosmolkit | 3.018887460 | [3.004519829, 3.034005040] | 24/24 | yes |
| 5000 | morgan2_chiral | parse_inclusive | rdkit | 4.617788615 | [4.604076056, 4.649089757] | 24/24 | yes |
| 5000 | morgan2_chiral | parse_inclusive | cosmolkit | 4.460085707 | [4.434480534, 4.475807352] | 24/24 | yes |
| 5000 | morgan3_2048 | first_use | rdkit | 1.810432108 | [1.800414640, 1.815329526] | 24/24 | yes |
| 5000 | morgan3_2048 | first_use | cosmolkit | 1.170340485 | [1.163129970, 1.173970820] | 24/24 | yes |
| 5000 | morgan3_2048 | repeated_use | rdkit | 4.826753572 | [4.800279092, 4.854209030] | 24/24 | yes |
| 5000 | morgan3_2048 | repeated_use | cosmolkit | 3.015556089 | [3.003455609, 3.031773297] | 24/24 | yes |
| 5000 | morgan3_2048 | parse_inclusive | rdkit | 4.577858978 | [4.560544483, 4.590294075] | 24/24 | yes |
| 5000 | morgan3_2048 | parse_inclusive | cosmolkit | 4.352199821 | [4.324443377, 4.378987154] | 24/24 | yes |
