# Python aromatic component storage and fused-ring FIFO — 2026-10-15

The production experiment is rejected for adoption; only its report and evidence packet are delivered. It is unmerged and does not establish complete superiority over RDKit or COSMolKit, or verify current main.

The candidate reuses the component member vector across aromatic connected components and uses a contiguous Vec with a monotonically increasing head as the fused-ring BFS FIFO. Parent initialization, discovery/adjacency order, early termination, path reconstruction, electron arithmetic and all Hückel predicates remain unchanged. No aromaticity shortcut, public API, dependency or unsafe code is added.

The immediate baseline includes all ten previously accepted source candidates, including aromatic normalization storage and the V3000 bounded-capacity/zero-map writer changes. All rejected experiments remain excluded. Baseline and candidate Python wheels are built freshly before and after this storage-only change; all source/input and extension hashes are verified.

Decision: Reject component reuse/fused FIFO for adoption as a cold-descriptor optimization. Initial first-use TPSA/ASA intervals cross 1 on both profiles: ChEMBL medians 1.000x and 1.002x, NCI 0.998x and 1.006x. ChEMBL parse-inclusive TPSA 1.007x [1.003,1.019] (15/20) and ASA 1.015x [1.007,1.035] (17/20) are small non-strict improvements, retained without reclassifying them as cold gains. No source-relative lane wins all20; no initial CI wholly below 1, hence no adaptive selection. Exact 14,072-input outputs, 1,220 native tests and 1,389 Python tests pass. Earlier sampling mixed parsing and getters; getter-ancestry review shows this gate is mostly parsing-side work. Entire experiment archived, not applied to production.

The cold target is not achieved; earlier sampling mixed parser/getter work. The getter-ancestry diagnostic review is retained in previous-getter-profile-review.json and is explicitly historical and non-additive. All 14,072 fixed inputs have byte-identical descriptor and Morgan dumps: 112,576 descriptor outcomes (first/repeated, exact float bits), 56,288 Morgan outputs (R2/R3 with both chirality settings), zero parse/evaluation errors. The added differential test freezes the original model gate and aromatic/fused-system implementations. On all 14,072 existing corpus inputs it compares direct aromatic-system decisions, the uncached gate and cold/warm memoized gates. Fused aromatic inputs are explicitly required to be present. No new chemical input records were added.

Scoped release Rust tests: 1,220 passed, 0 failed, 2 existing ignored. Full Python suite: 1,389 passed, one skipped. Release all-target Clippy for perception and chemistry, formatting and the refactor budget passed. The complete logs are retained. The initial baseline wheel build and first native-test compilation stopped on ENOSPC; their failure logs and unchanged-source successful retries are retained. Candidate Morgan output persistence also hit ENOSPC after both wheel builds completed; its exact partial bytes and failure log are retained. Only the completed owned build cache was then reclaimed, and full parity/Python validation completed before any timing. Only byte-verified duplicate staging/raw outcomes or owned completed Cargo targets are reclaimed, with cleanup records; primary source/evidence remains. Only this run owned regenerable targets were removed between phases. Timing starts after all builds and validation finish.

## Protocol

Pinned measured environment: RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, CPython 3.13.6, Rust 1.97.0, maturin 1.13.3, macOS arm64. The fixed, unfiltered slices contain 1,000 inputs each, offsets 0 and 5,000. Each uses 20 rotating fresh-process blocks across baseline, candidate and the two comparators. This is an exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.

Ratios are reference elapsed time / candidate elapsed time; greater than 1 means faster. Brackets are paired median-ratio 95% bootstrap intervals. Strict wins require semantic agreement, CI lower > 1 and all 20 blocks faster. Every operation, mode, block, error and mismatch is retained. Unchanged-operation movement is not evidence of a new algorithm.

First-use parses fresh objects outside the timer and may compare different eager/lazy work; repeated-use warms the operation once per object; parse-inclusive times actual fresh parsing and evaluation. Native return objects remain retained during timing; normalization and semantic checks occur outside timing. Whole-process peak RSS includes imports and preparation and is not per-operation allocation. Source, corpus, wheels and installed extension hashes are verified.

Initial comparator strict wins: first_use: 14/27, repeated_use: 24/27, parse_inclusive: 27/27. COSMolKit produces six unequal ring counts at offset 0, so all three corresponding modes remain unranked. The full remaining descriptor deficits are retained below.

## Initial matrix: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.000 [0.986, 1.010]; 10/20 | 1.469 [1.459, 1.494]; 20/20 **strict win** | 0.516 [0.503, 0.520]; 0/20 |
| tpsa/repeated_use | 1.030 [1.017, 1.069]; 17/20 | 10.073 [9.882, 10.211]; 20/20 **strict win** | 1.843 [1.800, 1.907]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.007 [1.003, 1.019]; 15/20 | 13.222 [13.184, 13.293]; 20/20 **strict win** | 12.405 [12.337, 12.459]; 20/20 **strict win** |
| labute_asa/first_use | 1.002 [0.997, 1.010]; 12/20 | 0.722 [0.717, 0.733]; 0/20 | 0.347 [0.336, 0.354]; 0/20 |
| labute_asa/repeated_use | 1.008 [0.994, 1.041]; 12/20 | 2.675 [2.630, 2.745]; 20/20 **strict win** | 1.438 [1.354, 1.544]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.015 [1.007, 1.035]; 17/20 | 12.462 [12.392, 12.570]; 20/20 **strict win** | 11.857 [11.764, 11.976]; 20/20 **strict win** |
| num_rings/first_use | 1.007 [0.995, 1.036]; 12/20 | 0.039 [0.038, 0.040]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 1.026 [1.009, 1.054]; 16/20 | 1.764 [1.689, 1.797]; 19/20 | unranked (output mismatch) |
| num_rings/parse_inclusive | 1.013 [1.001, 1.032]; 14/20 | 11.772 [11.668, 11.837]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 1.001 [0.993, 1.011]; 10/20 | 3.126 [3.098, 3.166]; 20/20 **strict win** | 2.040 [2.015, 2.063]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.996 [0.989, 1.020]; 8/20 | 5.264 [5.180, 5.296]; 20/20 **strict win** | 3.243 [3.202, 3.359]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.013 [0.998, 1.035]; 13/20 | 8.716 [8.617, 8.768]; 20/20 **strict win** | 7.688 [7.660, 7.745]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.006 [0.998, 1.020]; 13/20 | 1.271 [1.253, 1.288]; 20/20 **strict win** | 0.831 [0.825, 0.835]; 0/20 |
| morgan2_chiral/repeated_use | 1.001 [0.991, 1.011]; 10/20 | 5.005 [4.975, 5.109]; 20/20 **strict win** | 3.136 [3.116, 3.210]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.006 [1.000, 1.011]; 14/20 | 4.491 [4.463, 4.510]; 20/20 **strict win** | 3.951 [3.935, 3.998]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.000 [0.983, 1.011]; 11/20 | 3.349 [3.331, 3.387]; 20/20 **strict win** | 2.236 [2.221, 2.252]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.002 [0.991, 1.009]; 10/20 | 4.905 [4.851, 4.953]; 20/20 **strict win** | 3.154 [3.127, 3.167]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.008 [0.999, 1.010]; 14/20 | 8.068 [8.007, 8.089]; 20/20 **strict win** | 7.037 [7.011, 7.085]; 20/20 **strict win** |
| maccs/first_use | 1.000 [0.996, 1.007]; 10/20 | 1.886 [1.879, 1.894]; 20/20 **strict win** | 2.174 [2.169, 2.183]; 20/20 **strict win** |
| maccs/repeated_use | 0.999 [0.995, 1.001]; 7/20 | 1.946 [1.934, 1.954]; 20/20 **strict win** | 2.242 [2.232, 2.249]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.999 [0.995, 1.007]; 10/20 | 2.109 [2.099, 2.124]; 20/20 **strict win** | 2.372 [2.365, 2.430]; 20/20 **strict win** |

## Initial matrix: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.998 [0.996, 1.011]; 9/20 | 0.428 [0.420, 0.435]; 0/20 | 0.114 [0.111, 0.117]; 0/20 |
| tpsa/repeated_use | 1.004 [0.990, 1.030]; 10/20 | 9.750 [9.595, 9.937]; 20/20 **strict win** | 1.315 [1.270, 1.392]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.004 [0.993, 1.013]; 11/20 | 5.100 [5.040, 5.147]; 20/20 **strict win** | 5.253 [5.222, 5.308]; 20/20 **strict win** |
| labute_asa/first_use | 1.006 [0.994, 1.011]; 12/20 | 0.203 [0.200, 0.213]; 0/20 | 0.084 [0.080, 0.086]; 0/20 |
| labute_asa/repeated_use | 1.005 [0.990, 1.034]; 11/20 | 2.596 [2.546, 2.691]; 20/20 **strict win** | 1.150 [1.116, 1.212]; 19/20 |
| labute_asa/parse_inclusive | 1.003 [1.000, 1.006]; 14/20 | 4.843 [4.768, 4.901]; 20/20 **strict win** | 5.212 [5.175, 5.246]; 20/20 **strict win** |
| num_rings/first_use | 0.990 [0.980, 1.000]; 6/20 | 0.104 [0.099, 0.107]; 0/20 | 0.077 [0.072, 0.083]; 0/20 |
| num_rings/repeated_use | 1.003 [0.989, 1.023]; 11/20 | 1.727 [1.686, 1.763]; 20/20 **strict win** | 1.058 [1.024, 1.070]; 16/20 |
| num_rings/parse_inclusive | 1.007 [0.997, 1.019]; 14/20 | 11.893 [11.683, 12.000]; 20/20 **strict win** | 12.939 [12.900, 13.017]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.005 [0.993, 1.022]; 11/20 | 1.564 [1.528, 1.572]; 20/20 **strict win** | 1.005 [0.996, 1.019]; 12/20 |
| morgan2_2048/repeated_use | 1.005 [0.997, 1.019]; 12/20 | 4.735 [4.646, 4.773]; 20/20 **strict win** | 2.888 [2.857, 2.905]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.005 [0.998, 1.010]; 14/20 | 4.630 [4.582, 4.659]; 20/20 **strict win** | 4.482 [4.450, 4.514]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.003 [0.991, 1.006]; 11/20 | 1.596 [1.571, 1.612]; 20/20 **strict win** | 1.021 [1.016, 1.053]; 16/20 |
| morgan2_chiral/repeated_use | 0.998 [0.991, 1.007]; 9/20 | 4.713 [4.635, 4.772]; 20/20 **strict win** | 2.900 [2.874, 2.963]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.004 [1.001, 1.012]; 15/20 | 4.605 [4.569, 4.656]; 20/20 **strict win** | 4.459 [4.409, 4.488]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.001 [0.987, 1.013]; 11/20 | 1.779 [1.749, 1.788]; 20/20 **strict win** | 1.181 [1.161, 1.193]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.000 [0.995, 1.005]; 10/20 | 4.561 [4.482, 4.602]; 20/20 **strict win** | 2.911 [2.895, 2.944]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.000 [0.992, 1.011]; 10/20 | 4.531 [4.474, 4.559]; 20/20 **strict win** | 4.317 [4.274, 4.327]; 20/20 **strict win** |
| maccs/first_use | 0.996 [0.990, 1.000]; 6/20 | 1.629 [1.624, 1.638]; 20/20 **strict win** | 1.747 [1.742, 1.764]; 20/20 **strict win** |
| maccs/repeated_use | 0.993 [0.990, 1.000]; 6/20 | 1.729 [1.725, 1.736]; 20/20 **strict win** | 1.853 [1.846, 1.859]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.995 [0.990, 1.000]; 6/20 | 1.789 [1.782, 1.796]; 20/20 **strict win** | 1.923 [1.918, 1.933]; 20/20 **strict win** |

## Adaptive negative observations

The plan was saved before timing: any initial baseline-relative lane with CI upper < 1 selects its entire operation for one separate 20-block follow-up in all modes. Initial observations remain above; follow-ups are not pooled or substituted. New follow-up negatives remain recorded without recursive retries.

```json
{
  "5000": {
    "negative_initial_lanes": {},
    "operations": [],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  },
  "0": {
    "negative_initial_lanes": {},
    "operations": [],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  }
}
```

## Reproduction and status

[The evidence packet](2026-10-15-python-huckel-fifo-artifacts/README.md) contains both measured wheels, every raw matrix/follow-up, exact compressed output dumps, build/test logs, the original implementation, the patch and all frozen source/input files. The archive is read back and all 10058 candidate file hashes verified; baseline overrides and removal of the added test reconstruct all 10057 baseline hashes. The patch also reconstructs the changed files in an independent tiny Git repository. Hashes verify source reconstruction, not bit-reproducibility of rebuilt wheels.

Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared metadata was not refreshed; this is source-version 1.1.4 evidence with no hosted CI, merge or publication. Previous unresolved regressions against older candidates, including chiral Morgan parse-inclusive, ring/Morgan repeated-use and MACCS observations, are not resolved by an immediate-baseline comparison. Separate run win counts must not be treated as attributable code improvements.
