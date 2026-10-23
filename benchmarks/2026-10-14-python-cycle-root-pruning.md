# Python relevant-cycle root pruning and fixed-path reuse (2026-10-14)

First-use ring counting still trails eager RDKit ring-info lookup. This follow-up removes repeated shortest-path construction in relevant-cycle candidate generation; it preserves the ring definition, relevance reduction, family expansion, cap errors, and descriptor fallback. It moves no computation into parsing.

Baseline: source wheel built from main `cac480fb23cd8adeec6e5e9bed1f18677f87b392` (1.1.4). Candidate: the same source plus the [candidate patch](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-14-python-cycle-root-pruning-candidate.patch). Both are local source builds, not registry artifacts. See the [source manifest](2026-10-14-python-cycle-root-pruning-source.json) for hashes.

## Change and correctness

For each BFS root, the first predecessor defines a rooted tree. Two fixed paths intersect only at the root exactly when their first branches differ. Cache their edge bitsets once, compare the first branches, and XOR cached bitsets to construct a candidate. A candidate contains its highest-index root, which needs two eligible lower-index neighbours. Roots that fail this necessary condition are skipped before BFS. Alternative shortest paths are still expanded by the unchanged family code. Bitset rows are allocated only for reached 2-core vertices. The cache exists only during candidate generation for one root; it is not a molecule-level cache.

- Exact candidate sequence comparison against the frozen prior collector: all 1,024 labelled five-vertex graphs, 256 deterministic relabelled sparse/dense graphs, a 66-edge complete graph, disconnected/ineligible-edge cases, long 40-edge bridges crossing bitset word boundaries, and bounded path rows on graphs with long acyclic tails.
- Perception library: 254 passed, 1 ignored. Python suite: 1,389 passed, 1 skipped. Initial Rust build ran out of disk space; the complete retry passed.
- 14,072 inputs (10,000 fixed ordinary inputs plus 4,072 stress inputs): 112,576 first/repeated descriptor outcomes and 56,288 Morgan outcomes retain byte-identical dump files. This includes retained errors; equality does not turn error/nonfinite outcomes into valid chemistry. See [parity hashes](2026-10-14-python-cycle-root-pruning-parity.json).

## Measurement

CPython 3.13.6, rustc 1.97.0, maturin 1.13.3, macOS arm64, release/LTO defaults. Competitors: RDKit 2026.09.1 and COSMolKit 0.5.0rc22. The same unfiltered 1,000-input slices at offsets 0 and 5,000 of `rdkit-js-browser-10k-v1.smi` are used. Corpus SHA-256: `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`. Existing corpus provenance/licenses are unchanged. Both slices are exposed development data, not a sealed holdout.

Each slice runs 20 rotating four-arm blocks in fresh processes. First use times only the operation on a fresh parsed object; repeated use follows one untimed operation; parse-inclusive times actual parsing and the operation in the same loop. Native return objects are retained and normalization happens outside timing. Every block verifies outputs and installed extension hashes; all rows/blocks are retained. Whole-process RSS includes imports/loading, not just operation allocation.

Tables show comparator time / candidate time, paired-median 95% bootstrap interval (10,000 deterministic resamples), and faster blocks. Above 1 means candidate faster. A strict demonstrated win requires the lower interval endpoint above 1 and all 20 blocks faster. One host and seven operations do not establish whole-toolkit superiority. The pre-existing macrocycle threshold/bounded fallback is unchanged and remains a chemistry limitation.

## Ring counting

| Offset | Mode | Previous main / candidate | RDKit / candidate | COSMolKit / candidate |
|---:|---|---|---|---|
| 0 | first_use | 1.765 [1.662, 1.784]; 20/20 | 0.033 [0.030, 0.035]; 0/20 | unranked |
| 0 | repeated_use | 1.012 [0.973, 1.026]; 11/20 | 1.666 [1.578, 1.772]; 18/20 | unranked |
| 0 | parse_inclusive | 1.284 [1.266, 1.308]; 20/20 | 11.093 [9.945, 11.441]; 20/20 | unranked |
| 5000 | first_use | 0.981 [0.941, 1.018]; 7/20 | 0.099 [0.094, 0.108]; 0/20 | 0.098 [0.089, 0.110]; 0/20 |
| 5000 | repeated_use | 0.995 [0.906, 1.075]; 10/20 | 1.703 [1.576, 1.833]; 18/20 | 1.085 [0.999, 1.142]; 13/20 |
| 5000 | parse_inclusive | 0.973 [0.928, 0.999]; 6/20 | 11.413 [11.057, 12.989]; 20/20 | 12.729 [12.134, 13.838]; 20/20 |

## Strict comparator wins across seven operations

| Mode | Demonstrated / eligible comparisons |
|---|---:|
| first_use | 8/27 |
| repeated_use | 19/27 |
| parse_inclusive | 26/27 |

COSMolKit ring-count comparisons are suppressed for any slice whose outputs disagree with RDKit. See every mismatch count in the raw agreement records. No mismatching rows are removed.

## Complete matrix

| Offset | Operation / mode | Previous main / candidate | RDKit / candidate | COSMolKit / candidate |
|---:|---|---|---|---|
| 0 | tpsa/first_use | 0.999 [0.984, 1.024]; 10/20 | 1.401 [1.338, 1.435]; 18/20 | 0.455 [0.400, 0.486]; 0/20 |
| 0 | tpsa/repeated_use | 0.994 [0.965, 1.105]; 9/20 | 9.585 [8.918, 9.805]; 20/20 | 1.507 [1.298, 1.605]; 17/20 |
| 0 | tpsa/parse_inclusive | 0.996 [0.973, 1.007]; 6/20 | 12.855 [12.087, 13.042]; 20/20 | 11.948 [10.950, 12.227]; 20/20 |
| 0 | labute_asa/first_use | 0.998 [0.993, 1.007]; 9/20 | 0.681 [0.652, 0.699]; 0/20 | 0.314 [0.290, 0.330]; 0/20 |
| 0 | labute_asa/repeated_use | 0.988 [0.974, 1.008]; 7/20 | 2.596 [2.410, 2.715]; 19/20 | 1.362 [1.313, 1.454]; 16/20 |
| 0 | labute_asa/parse_inclusive | 0.995 [0.970, 1.005]; 9/20 | 12.169 [10.821, 12.322]; 20/20 | 11.537 [10.570, 11.877]; 20/20 |
| 0 | num_rings/first_use | 1.765 [1.662, 1.784]; 20/20 | 0.033 [0.030, 0.035]; 0/20 | unranked |
| 0 | num_rings/repeated_use | 1.012 [0.973, 1.026]; 11/20 | 1.666 [1.578, 1.772]; 18/20 | unranked |
| 0 | num_rings/parse_inclusive | 1.284 [1.266, 1.308]; 20/20 | 11.093 [9.945, 11.441]; 20/20 | unranked |
| 0 | morgan2_2048/first_use | 1.002 [0.972, 1.015]; 11/20 | 3.072 [2.746, 3.137]; 20/20 | 1.982 [1.901, 2.032]; 20/20 |
| 0 | morgan2_2048/repeated_use | 0.989 [0.983, 1.004]; 7/20 | 5.218 [4.602, 5.250]; 20/20 | 3.147 [3.030, 3.239]; 20/20 |
| 0 | morgan2_2048/parse_inclusive | 0.990 [0.966, 1.002]; 7/20 | 8.486 [7.653, 8.651]; 20/20 | 7.429 [7.133, 7.658]; 20/20 |
| 0 | morgan2_chiral/first_use | 0.986 [0.975, 1.003]; 8/20 | 1.245 [1.089, 1.292]; 16/20 | 0.792 [0.763, 0.814]; 2/20 |
| 0 | morgan2_chiral/repeated_use | 1.004 [0.976, 1.046]; 11/20 | 5.018 [4.448, 5.295]; 20/20 | 3.122 [2.999, 3.308]; 20/20 |
| 0 | morgan2_chiral/parse_inclusive | 1.009 [0.991, 1.058]; 12/20 | 4.434 [3.916, 4.555]; 20/20 | 3.977 [3.874, 4.213]; 20/20 |
| 0 | morgan3_2048/first_use | 1.009 [0.984, 1.042]; 14/20 | 3.204 [2.927, 3.368]; 20/20 | 2.249 [2.190, 2.318]; 20/20 |
| 0 | morgan3_2048/repeated_use | 0.997 [0.963, 1.015]; 10/20 | 4.731 [4.417, 4.937]; 20/20 | 3.146 [2.979, 3.208]; 20/20 |
| 0 | morgan3_2048/parse_inclusive | 0.999 [0.958, 1.025]; 10/20 | 7.913 [6.998, 8.217]; 20/20 | 6.885 [6.550, 7.010]; 20/20 |
| 0 | maccs/first_use | 1.004 [0.977, 1.029]; 12/20 | 1.873 [1.714, 1.895]; 20/20 | 2.177 [2.142, 2.209]; 20/20 |
| 0 | maccs/repeated_use | 0.998 [0.959, 1.028]; 10/20 | 1.940 [1.809, 1.967]; 20/20 | 2.256 [2.225, 2.309]; 20/20 |
| 0 | maccs/parse_inclusive | 0.999 [0.975, 1.020]; 10/20 | 2.111 [2.047, 2.131]; 20/20 | 2.383 [2.322, 2.411]; 20/20 |
| 5000 | tpsa/first_use | 0.941 [0.870, 1.007]; 8/20 | 0.348 [0.335, 0.370]; 0/20 | 0.089 [0.083, 0.097]; 0/20 |
| 5000 | tpsa/repeated_use | 0.927 [0.874, 1.045]; 7/20 | 7.833 [6.956, 9.389]; 20/20 | 1.227 [1.099, 1.316]; 18/20 |
| 5000 | tpsa/parse_inclusive | 0.982 [0.953, 1.002]; 7/20 | 4.453 [4.167, 5.048]; 20/20 | 4.404 [4.299, 4.791]; 20/20 |
| 5000 | labute_asa/first_use | 0.983 [0.946, 1.035]; 7/20 | 0.172 [0.162, 0.204]; 0/20 | 0.070 [0.066, 0.075]; 0/20 |
| 5000 | labute_asa/repeated_use | 0.970 [0.886, 1.031]; 8/20 | 2.151 [1.904, 2.544]; 20/20 | 1.134 [0.918, 1.251]; 13/20 |
| 5000 | labute_asa/parse_inclusive | 0.990 [0.956, 1.004]; 8/20 | 4.109 [4.006, 4.667]; 20/20 | 4.413 [4.037, 4.727]; 20/20 |
| 5000 | num_rings/first_use | 0.981 [0.941, 1.018]; 7/20 | 0.099 [0.094, 0.108]; 0/20 | 0.098 [0.089, 0.110]; 0/20 |
| 5000 | num_rings/repeated_use | 0.995 [0.906, 1.075]; 10/20 | 1.703 [1.576, 1.833]; 18/20 | 1.085 [0.999, 1.142]; 13/20 |
| 5000 | num_rings/parse_inclusive | 0.973 [0.928, 0.999]; 6/20 | 11.413 [11.057, 12.989]; 20/20 | 12.729 [12.134, 13.838]; 20/20 |
| 5000 | morgan2_2048/first_use | 0.990 [0.962, 1.011]; 8/20 | 1.400 [1.361, 1.579]; 19/20 | 0.892 [0.871, 0.977]; 4/20 |
| 5000 | morgan2_2048/repeated_use | 0.974 [0.897, 1.012]; 7/20 | 4.757 [4.307, 5.217]; 20/20 | 2.882 [2.624, 3.237]; 20/20 |
| 5000 | morgan2_2048/parse_inclusive | 0.982 [0.928, 1.001]; 6/20 | 4.192 [3.913, 4.346]; 20/20 | 4.076 [3.918, 4.409]; 20/20 |
| 5000 | morgan2_chiral/first_use | 0.994 [0.911, 1.019]; 8/20 | 1.417 [1.281, 1.467]; 19/20 | 0.920 [0.819, 0.989]; 5/20 |
| 5000 | morgan2_chiral/repeated_use | 0.991 [0.958, 1.009]; 8/20 | 4.625 [4.384, 4.754]; 20/20 | 2.883 [2.751, 3.241]; 20/20 |
| 5000 | morgan2_chiral/parse_inclusive | 0.988 [0.954, 0.993]; 5/20 | 4.206 [4.080, 4.311]; 20/20 | 4.226 [3.980, 4.377]; 20/20 |
| 5000 | morgan3_2048/first_use | 1.002 [0.938, 1.008]; 11/20 | 1.585 [1.509, 1.621]; 19/20 | 1.077 [1.048, 1.164]; 16/20 |
| 5000 | morgan3_2048/repeated_use | 1.008 [0.994, 1.082]; 13/20 | 4.567 [4.464, 4.847]; 20/20 | 3.030 [2.907, 3.201]; 20/20 |
| 5000 | morgan3_2048/parse_inclusive | 0.971 [0.931, 1.021]; 6/20 | 4.162 [4.046, 4.363]; 20/20 | 3.981 [3.842, 4.248]; 20/20 |
| 5000 | maccs/first_use | 0.989 [0.958, 1.002]; 7/20 | 1.625 [1.592, 1.651]; 20/20 | 1.749 [1.676, 1.817]; 20/20 |
| 5000 | maccs/repeated_use | 0.999 [0.960, 1.012]; 10/20 | 1.746 [1.664, 1.773]; 20/20 | 1.849 [1.809, 1.902]; 20/20 |
| 5000 | maccs/parse_inclusive | 0.989 [0.942, 1.014]; 9/20 | 1.806 [1.694, 1.970]; 20/20 | 1.940 [1.870, 1.987]; 19/20 |

## Raw evidence and reproduction

- [Slice 0, all 20 raw blocks](2026-10-14-python-cycle-root-pruning-slice0-paired20.json)
- [Slice 5000, all 20 raw blocks](2026-10-14-python-cycle-root-pruning-slice5000-paired20.json)
- [Source manifest](2026-10-14-python-cycle-root-pruning-source.json) and [parity hashes](2026-10-14-python-cycle-root-pruning-parity.json).
- Compressed baseline/candidate dumps: `validation/results/python-cycle-root-pruning-{baseline,candidate}-{descriptors,morgan}-2026-10-14.jsonl.gz`.

```sh
python scripts/bench_python_speed_frontier.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --offset 0 --limit 1000 --blocks 20 \
  --base-revision cac480fb23cd8adeec6e5e9bed1f18677f87b392 \
  --source-manifest benchmarks/2026-10-14-python-cycle-root-pruning-source.json \
  --baseline-wheel <baseline.whl> --baseline-path <baseline-install> \
  --candidate-wheel <candidate.whl> --candidate-path <candidate-install> \
  --output <slice0.json>
```

Repeat with `--offset 5000`; use the pinned competitor versions and same build flags. Follow-up work remains in first-use compatible-model construction, ring finding, and any fingerprint lanes that do not meet the strict criterion.

## Checks and retained intermediate evidence

- [Validation logs and hashes](../validation/results/python-cycle-root-pruning-checks-2026-10-14.json)
- [Completed one-slice intermediate fixed-path candidate](2026-10-14-python-cycle-path-reuse.md). Its source patch, exact dumps and all 20 raw blocks remain available; it was measured before root pruning and reached-row allocation.

Baseline comparisons whose entire 95% interval is below 1 (all matrix rows retained):

- Offset 5000, `num_rings/parse_inclusive`: 0.973 [0.928, 0.999]; 6/20
- Offset 5000, `morgan2_chiral/parse_inclusive`: 0.988 [0.954, 0.993]; 5/20


## Adaptive focused repeat on offset 5,000

The two baseline-relative parse-inclusive intervals below 1 in the complete matrix prompted a separate 20-block run of `num_rings,morgan2_chiral`. This is an adaptive follow-up on the same exposed data, not an independent corpus. The original matrix and every follow-up block are retained. Source and wheel hashes are unchanged. No samples are pooled or discarded, and the strict-win counts above still refer only to the complete seven-operation matrix.

| Operation / mode | Previous main / candidate | Faster blocks |
|---|---|---:|
| num_rings/first_use | 1.000 [0.977, 1.012] | 10/20 |
| num_rings/repeated_use | 1.031 [0.903, 1.092] | 11/20 |
| num_rings/parse_inclusive | 0.980 [0.946, 1.007] | 8/20 |
| morgan2_chiral/first_use | 0.994 [0.921, 1.037] | 9/20 |
| morgan2_chiral/repeated_use | 0.974 [0.914, 1.007] | 9/20 |
| morgan2_chiral/parse_inclusive | 0.961 [0.902, 0.998] | 5/20 |

[All focused-repeat blocks and output gates](2026-10-14-python-cycle-root-pruning-slice5000-focused-repeat20.json). Reproduce with the earlier command plus `--offset 5000 --ops num_rings,morgan2_chiral` and a new output path.

The focused ring parse-inclusive interval crosses 1, so it did not independently confirm the initial ring slowdown. The chiral Morgan parse-inclusive interval remains below 1 (`0.961 [0.902, 0.998]`, 5/20 faster blocks): this observed regression remains unresolved. These timings do not identify its cause, and this candidate is not a uniformly faster replacement for main. The source changes remain isolated for review; no new merge or release is asserted.
