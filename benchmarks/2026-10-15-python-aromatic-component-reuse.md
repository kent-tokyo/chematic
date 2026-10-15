# Python aromaticity component reuse — 2026-10-15

The production change is adopted in the attached source candidate. It is unmerged and does not establish complete superiority over RDKit or COSMolKit, or verify current main.

During RDKit-compatible aromaticity typing, the bridge pass already provides the connected components of the cyclic subgraph. The candidate reuses those cached labels instead of building an identical union-find and repeatedly finding its roots. The separate candidate-only cycle-detection forest remains unchanged. Donor typing, candidate predicates, ring order and selection, Hückel evaluation and aromaticity application are unchanged. No public API, dependency or unsafe code is added.

The immediate baseline includes all six previous accepted source candidates: cycle-root pruning, V3000 writer buffers, shared Kekulé descriptor views, V3000 parser token reuse, ring-component counting and the no-double-bond E/Z skip. It excludes the rejected direct-line experiment. Both Python wheels were built freshly before/after this single production change; this avoids attributing inherited wins to the new optimization.

Decision: Adopt as a small unmerged source optimization: exact outputs and all scoped gates pass; NCI first-use TPSA is 1.019x [1.007, 1.030], 17/20 blocks faster, and parse-inclusive TPSA is 1.013x [1.004, 1.022], 15/20. Neither satisfies the all-block strict criterion. Labute ASA improvement is unproved. Five initial offset-0 Morgan negative lanes select one separate follow-up; none has a wholly-below-1 interval there. Initial negatives remain visible and the follow-up does not prove absence of regression. No claim of toolkit-wide or main/release superiority.

All 14,072 fixed inputs have byte-identical descriptor and Morgan dumps: 112,576 descriptor outcomes (first/repeated, exact float bits), 56,288 Morgan outputs (R2/R3 with both chirality settings), zero parse/evaluation errors. Added tests compare donor types, candidate flags and kept component masks with the frozen original implementation on all 10,000 development molecules, with cold/warm calls and separately seeded valid bond-flag/component caches. No new chemical input records were added.

Scoped release Rust tests: 1,218 passed, 0 failed, 2 existing ignored. Full Python suite: 1,389 passed, one skipped. Release all-target Clippy for perception and chemistry, formatting and the refactor budget passed. The complete logs are retained, including the initial native-test ENOSPC link failure and successful retry. Only this run owned regenerable targets were removed between phases. Timing starts after all builds and validation finish.

## Protocol

Pinned measured environment: RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, CPython 3.13.6, Rust 1.97.0, maturin 1.13.3, macOS arm64. The fixed, unfiltered slices contain 1,000 inputs each, offsets 0 and 5,000. Each uses 20 rotating fresh-process blocks across baseline, candidate and the two comparators. This is an exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.

Ratios are reference elapsed time / candidate elapsed time; greater than 1 means faster. Brackets are paired median-ratio 95% bootstrap intervals. Strict wins require semantic agreement, CI lower > 1 and all 20 blocks faster. Every operation, mode, block, error and mismatch is retained. Unchanged-operation movement is not evidence of a new algorithm.

First-use parses fresh objects outside the timer and may compare different eager/lazy work; repeated-use warms the operation once per object; parse-inclusive times actual fresh parsing and evaluation. Native return objects remain retained during timing; normalization and semantic checks occur outside timing. Whole-process peak RSS includes imports and preparation and is not per-operation allocation. Source, corpus, wheels and installed extension hashes are verified.

Initial comparator strict wins: first_use: 11/27, repeated_use: 21/27, parse_inclusive: 27/27. COSMolKit produces six unequal ring counts at offset 0, so all three corresponding modes remain unranked. The full remaining descriptor deficits are retained below.

## Initial matrix: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.002 [0.985, 1.008]; 11/20 | 1.449 [1.280, 1.512]; 19/20 | 0.514 [0.440, 0.521]; 0/20 |
| tpsa/repeated_use | 0.984 [0.957, 1.003]; 7/20 | 9.415 [8.235, 9.838]; 20/20 **strict win** | 1.599 [1.436, 1.652]; 19/20 |
| tpsa/parse_inclusive | 1.000 [0.976, 1.015]; 10/20 | 12.910 [11.574, 13.138]; 20/20 **strict win** | 12.169 [10.460, 12.268]; 20/20 **strict win** |
| labute_asa/first_use | 0.987 [0.917, 1.006]; 9/20 | 0.673 [0.607, 0.720]; 0/20 | 0.315 [0.288, 0.345]; 0/20 |
| labute_asa/repeated_use | 0.976 [0.937, 1.002]; 7/20 | 2.479 [2.325, 2.609]; 20/20 **strict win** | 1.338 [1.262, 1.433]; 18/20 |
| labute_asa/parse_inclusive | 0.990 [0.926, 1.004]; 6/20 | 12.060 [10.538, 12.324]; 20/20 **strict win** | 11.445 [10.046, 11.677]; 20/20 **strict win** |
| num_rings/first_use | 0.980 [0.958, 1.001]; 6/20 | 0.036 [0.034, 0.038]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 0.989 [0.947, 1.000]; 6/20 | 1.686 [1.520, 1.760]; 18/20 | unranked (output mismatch) |
| num_rings/parse_inclusive | 1.000 [0.950, 1.004]; 9/20 | 11.535 [10.275, 11.639]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 0.987 [0.973, 0.999]; 5/20 | 3.071 [2.864, 3.117]; 20/20 **strict win** | 2.009 [1.741, 2.028]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.992 [0.980, 0.998]; 5/20 | 5.060 [4.638, 5.167]; 20/20 **strict win** | 3.145 [2.917, 3.188]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.979 [0.922, 1.003]; 7/20 | 8.274 [7.406, 8.619]; 20/20 **strict win** | 7.377 [6.468, 7.620]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.999 [0.956, 1.011]; 9/20 | 1.245 [1.166, 1.260]; 18/20 | 0.818 [0.733, 0.836]; 0/20 |
| morgan2_chiral/repeated_use | 0.993 [0.968, 0.998]; 5/20 | 4.928 [4.531, 5.048]; 20/20 **strict win** | 3.091 [2.826, 3.123]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.975 [0.942, 0.997]; 5/20 | 4.337 [4.021, 4.433]; 20/20 **strict win** | 3.857 [3.568, 3.909]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.999 [0.987, 1.013]; 8/20 | 3.327 [3.130, 3.349]; 20/20 **strict win** | 2.198 [1.998, 2.235]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.986 [0.967, 0.999]; 6/20 | 4.772 [4.566, 4.845]; 20/20 **strict win** | 3.028 [2.941, 3.146]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.996 [0.989, 1.015]; 8/20 | 7.833 [7.509, 7.957]; 20/20 **strict win** | 6.926 [6.771, 6.996]; 20/20 **strict win** |
| maccs/first_use | 1.001 [0.994, 1.012]; 11/20 | 1.860 [1.800, 1.877]; 20/20 **strict win** | 2.165 [2.141, 2.183]; 20/20 **strict win** |
| maccs/repeated_use | 1.002 [0.996, 1.012]; 12/20 | 1.936 [1.825, 1.945]; 20/20 **strict win** | 2.236 [2.214, 2.248]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.998 [0.993, 1.002]; 7/20 | 2.098 [1.998, 2.112]; 20/20 **strict win** | 2.382 [2.363, 2.390]; 20/20 **strict win** |

## Initial matrix: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.019 [1.007, 1.030]; 17/20 | 0.410 [0.393, 0.422]; 0/20 | 0.106 [0.102, 0.119]; 0/20 |
| tpsa/repeated_use | 0.986 [0.943, 1.021]; 7/20 | 9.665 [8.605, 9.839]; 20/20 **strict win** | 1.313 [1.163, 1.455]; 19/20 |
| tpsa/parse_inclusive | 1.013 [1.004, 1.022]; 15/20 | 4.921 [4.661, 5.087]; 20/20 **strict win** | 5.065 [4.879, 5.354]; 20/20 **strict win** |
| labute_asa/first_use | 1.008 [0.999, 1.064]; 14/20 | 0.198 [0.184, 0.205]; 0/20 | 0.078 [0.075, 0.085]; 0/20 |
| labute_asa/repeated_use | 0.993 [0.939, 1.025]; 8/20 | 2.613 [2.442, 2.667]; 20/20 **strict win** | 1.243 [1.116, 1.315]; 18/20 |
| labute_asa/parse_inclusive | 1.011 [0.997, 1.060]; 12/20 | 4.700 [4.630, 4.925]; 20/20 **strict win** | 4.987 [4.840, 5.200]; 20/20 **strict win** |
| num_rings/first_use | 1.001 [0.955, 1.042]; 10/20 | 0.104 [0.094, 0.109]; 0/20 | 0.083 [0.075, 0.091]; 0/20 |
| num_rings/repeated_use | 1.000 [0.952, 1.018]; 10/20 | 1.756 [1.652, 1.867]; 20/20 **strict win** | 1.070 [1.011, 1.134]; 15/20 |
| num_rings/parse_inclusive | 1.000 [0.985, 1.049]; 11/20 | 11.848 [11.637, 12.289]; 20/20 **strict win** | 12.750 [12.482, 12.888]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.008 [0.984, 1.058]; 11/20 | 1.512 [1.464, 1.554]; 20/20 **strict win** | 0.950 [0.925, 0.977]; 4/20 |
| morgan2_2048/repeated_use | 1.009 [0.984, 1.019]; 13/20 | 4.744 [4.450, 4.862]; 20/20 **strict win** | 2.880 [2.794, 2.956]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.009 [0.974, 1.033]; 11/20 | 4.440 [4.252, 4.509]; 20/20 **strict win** | 4.329 [4.138, 4.373]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.018 [0.978, 1.048]; 11/20 | 1.527 [1.429, 1.583]; 20/20 **strict win** | 0.968 [0.939, 0.989]; 5/20 |
| morgan2_chiral/repeated_use | 1.002 [0.988, 1.067]; 11/20 | 4.705 [4.494, 4.811]; 20/20 **strict win** | 2.872 [2.735, 3.014]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.013 [0.989, 1.059]; 12/20 | 4.442 [4.215, 4.518]; 20/20 **strict win** | 4.273 [4.058, 4.490]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.014 [0.996, 1.031]; 13/20 | 1.699 [1.600, 1.757]; 20/20 **strict win** | 1.118 [1.068, 1.155]; 16/20 |
| morgan3_2048/repeated_use | 1.009 [0.963, 1.033]; 11/20 | 4.487 [4.179, 4.555]; 20/20 **strict win** | 2.891 [2.678, 3.003]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.007 [0.994, 1.077]; 13/20 | 4.314 [4.066, 4.477]; 20/20 **strict win** | 4.179 [3.958, 4.225]; 20/20 **strict win** |
| maccs/first_use | 0.998 [0.979, 1.029]; 10/20 | 1.630 [1.467, 1.649]; 20/20 **strict win** | 1.753 [1.679, 1.784]; 20/20 **strict win** |
| maccs/repeated_use | 1.001 [0.991, 1.020]; 11/20 | 1.707 [1.590, 1.756]; 20/20 **strict win** | 1.902 [1.854, 1.962]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.994 [0.968, 1.009]; 9/20 | 1.762 [1.660, 1.812]; 20/20 **strict win** | 1.957 [1.926, 2.025]; 20/20 **strict win** |

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
    "negative_initial_lanes": {
      "morgan2_2048/first_use": {
        "median": 0.987,
        "ci95": [
          0.973,
          0.999
        ],
        "winning_blocks": 5,
        "blocks": 20
      },
      "morgan2_2048/repeated_use": {
        "median": 0.992,
        "ci95": [
          0.98,
          0.998
        ],
        "winning_blocks": 5,
        "blocks": 20
      },
      "morgan2_chiral/repeated_use": {
        "median": 0.993,
        "ci95": [
          0.968,
          0.998
        ],
        "winning_blocks": 5,
        "blocks": 20
      },
      "morgan2_chiral/parse_inclusive": {
        "median": 0.975,
        "ci95": [
          0.942,
          0.997
        ],
        "winning_blocks": 5,
        "blocks": 20
      },
      "morgan3_2048/repeated_use": {
        "median": 0.986,
        "ci95": [
          0.967,
          0.999
        ],
        "winning_blocks": 6,
        "blocks": 20
      }
    },
    "operations": [
      "morgan2_2048",
      "morgan2_chiral",
      "morgan3_2048"
    ],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  }
}
```

### Separate follow-up: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| morgan2_2048/first_use | 1.006 [1.002, 1.013]; 15/20 | 3.073 [3.059, 3.096]; 20/20 **strict win** | 2.006 [1.990, 2.025]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 1.001 [0.997, 1.005]; 11/20 | 5.306 [5.258, 5.374]; 20/20 **strict win** | 3.222 [3.184, 3.237]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.003 [0.996, 1.015]; 12/20 | 8.572 [8.449, 8.625]; 20/20 **strict win** | 7.540 [7.514, 7.574]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.001 [0.993, 1.007]; 11/20 | 1.263 [1.244, 1.271]; 20/20 **strict win** | 0.828 [0.814, 0.833]; 0/20 |
| morgan2_chiral/repeated_use | 1.001 [0.989, 1.007]; 11/20 | 4.989 [4.939, 5.019]; 20/20 **strict win** | 3.089 [3.054, 3.137]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.997 [0.995, 1.001]; 6/20 | 4.434 [4.401, 4.445]; 20/20 **strict win** | 3.900 [3.878, 3.915]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.007 [0.991, 1.012]; 12/20 | 3.360 [3.336, 3.392]; 20/20 **strict win** | 2.225 [2.205, 2.244]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.997 [0.976, 1.002]; 7/20 | 4.830 [4.730, 4.911]; 20/20 **strict win** | 3.113 [3.063, 3.169]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.998 [0.988, 1.008]; 9/20 | 7.942 [7.843, 8.014]; 20/20 **strict win** | 6.937 [6.896, 6.976]; 20/20 **strict win** |

## Reproduction and status

[The evidence packet](2026-10-15-python-aromatic-component-reuse-artifacts/README.md) contains both measured wheels, every raw matrix/follow-up, exact compressed output dumps, build/test logs, the original implementation, the patch and all frozen source/input files. The archive is read back and all 10050 candidate file hashes verified; baseline overrides and removal of the added test reconstruct all 10049 baseline hashes. The patch also reconstructs the changed files in an independent tiny Git repository. Hashes verify source reconstruction, not bit-reproducibility of rebuilt wheels.

Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared metadata was not refreshed; this is source-version 1.1.4 evidence with no hosted CI, merge or publication. Previous unresolved regressions against older candidates, including chiral Morgan parse-inclusive, ring/Morgan repeated-use and MACCS observations, are not resolved by an immediate-baseline comparison. Separate run win counts must not be treated as attributable code improvements.
