# Python inline bridge and aromatic-component scratch — 2026-10-15

The production experiment is rejected for adoption; only its report and evidence packet are delivered. It is unmerged and does not establish complete superiority over RDKit or COSMolKit, or verify current main.

The candidate replaces three temporary bridge-search vectors (discovery/low-link pairs, DFS frames, visit stack) and four private aromatic-component fact vectors (roots, kinds, edge and atom counts) with existing smallvec inline storage for up to 64 elements. Larger buffers spill to dynamic storage. Public cached bond flags and RingComponents labels remain ordinary Vec values. DFS traversal, component numbering, donor/candidate/Hückel predicates and arithmetic are unchanged. No public API, dependency or unsafe code is added. The preceding FIFO trial is rejected, archived separately and excluded from this source candidate.

The immediate baseline includes all ten previously accepted source candidates, including aromatic normalization storage and the V3000 bounded-capacity/zero-map writer changes. All rejected experiments remain excluded. The exact baseline Python wheel is reused after matching all 10,057 baseline source/input hashes and its wheel/installed-extension hashes. It was built before either trial; the candidate is built freshly. New timing observations are separate and never pooled. The origin is recorded in baseline-artifact-origin.json.

Decision: Reject inline SmallVec scratch. Initial first-use ChEMBL TPSA 0.914x [0.909,0.922], ASA 0.917x [0.903,0.928], and ring count 0.939x [0.935,0.953] are slower. Separate predeclared follow-ups retain the TPSA/ASA regressions at 0.912x/0.916x and NCI ring count at 0.937x. Both full initial matrices and all selected follow-ups remain recorded, including new negatives. No source-relative lane is a strict win. Exact corpus parity and all tests pass, but speed does not justify adoption. No production change is delivered.

All 14,072 fixed inputs have byte-identical descriptor and Morgan dumps: 112,576 descriptor outcomes (first/repeated, exact float bits), 56,288 Morgan outputs (R2/R3 with both chirality settings), zero parse/evaluation errors. Two differential tests freeze the original bridge search and ComponentInfo constructor. All 14,072 inputs compare exact cyclic-bond flags, component label numbers/counts and every component fact field, including cold/warm calls. Both inline and greater-than-64-atom spill cases are required to be present. No new chemical input records were added.

Scoped release Rust tests: 1,221 passed, 0 failed, 2 existing ignored. Full Python suite: 1,389 passed, one skipped. Release all-target Clippy for perception and chemistry, formatting and the refactor budget passed. The complete logs are retained. The reused baseline build originally hit ENOSPC; its failure and successful retry logs are retained with the baseline origin record. The first current candidate wheel build also hit ENOSPC while creating the PyO3 archive; its unchanged-source retry and verified prior wheel-staging cleanup record are retained. Successful current candidate validation logs are retained. Only this run owned regenerable targets were removed between phases. Timing starts after all builds and validation finish.

## Protocol

Pinned measured environment: RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, CPython 3.13.6, Rust 1.97.0, maturin 1.13.3, macOS arm64. The fixed, unfiltered slices contain 1,000 inputs each, offsets 0 and 5,000. Each uses 20 rotating fresh-process blocks across baseline, candidate and the two comparators. This is an exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.

Ratios are reference elapsed time / candidate elapsed time; greater than 1 means faster. Brackets are paired median-ratio 95% bootstrap intervals. Strict wins require semantic agreement, CI lower > 1 and all 20 blocks faster. Every operation, mode, block, error and mismatch is retained. Unchanged-operation movement is not evidence of a new algorithm.

First-use parses fresh objects outside the timer and may compare different eager/lazy work; repeated-use warms the operation once per object; parse-inclusive times actual fresh parsing and evaluation. Native return objects remain retained during timing; normalization and semantic checks occur outside timing. Whole-process peak RSS includes imports and preparation and is not per-operation allocation. Source, corpus, wheels and installed extension hashes are verified.

Initial comparator strict wins: first_use: 14/27, repeated_use: 24/27, parse_inclusive: 27/27. COSMolKit produces six unequal ring counts at offset 0, so all three corresponding modes remain unranked. The full remaining descriptor deficits are retained below.

## Initial matrix: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.914 [0.909, 0.922]; 2/20 | 1.357 [1.317, 1.376]; 20/20 **strict win** | 0.469 [0.458, 0.483]; 0/20 |
| tpsa/repeated_use | 0.989 [0.960, 1.016]; 9/20 | 9.780 [9.414, 9.998]; 20/20 **strict win** | 1.716 [1.617, 1.796]; 19/20 |
| tpsa/parse_inclusive | 0.965 [0.955, 0.975]; 3/20 | 12.681 [12.510, 12.820]; 20/20 **strict win** | 11.821 [11.755, 11.958]; 20/20 **strict win** |
| labute_asa/first_use | 0.917 [0.903, 0.928]; 0/20 | 0.671 [0.641, 0.684]; 0/20 | 0.313 [0.306, 0.320]; 0/20 |
| labute_asa/repeated_use | 0.997 [0.983, 1.019]; 9/20 | 2.701 [2.608, 2.724]; 20/20 **strict win** | 1.385 [1.300, 1.457]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.967 [0.959, 0.982]; 4/20 | 11.934 [11.803, 12.086]; 20/20 **strict win** | 11.397 [11.189, 11.497]; 20/20 **strict win** |
| num_rings/first_use | 0.939 [0.935, 0.953]; 3/20 | 0.035 [0.034, 0.037]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 1.005 [0.988, 1.033]; 10/20 | 1.741 [1.667, 1.826]; 19/20 | unranked (output mismatch) |
| num_rings/parse_inclusive | 0.978 [0.971, 0.987]; 4/20 | 11.246 [11.177, 11.359]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 0.978 [0.966, 0.987]; 3/20 | 3.059 [3.043, 3.070]; 20/20 **strict win** | 1.978 [1.963, 2.001]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 1.002 [0.996, 1.012]; 12/20 | 5.247 [5.221, 5.310]; 20/20 **strict win** | 3.230 [3.188, 3.268]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.993 [0.986, 0.998]; 6/20 | 8.532 [8.441, 8.638]; 20/20 **strict win** | 7.507 [7.455, 7.551]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.995 [0.980, 1.004]; 7/20 | 1.266 [1.240, 1.289]; 20/20 **strict win** | 0.821 [0.815, 0.829]; 0/20 |
| morgan2_chiral/repeated_use | 1.004 [0.999, 1.007]; 12/20 | 5.103 [4.991, 5.122]; 20/20 **strict win** | 3.131 [3.115, 3.201]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.996 [0.992, 1.003]; 7/20 | 4.402 [4.363, 4.447]; 20/20 **strict win** | 3.906 [3.846, 3.938]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.984 [0.970, 0.993]; 5/20 | 3.273 [3.242, 3.328]; 20/20 **strict win** | 2.188 [2.177, 2.210]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.999 [0.992, 1.010]; 10/20 | 4.932 [4.834, 4.989]; 20/20 **strict win** | 3.138 [3.100, 3.177]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.990 [0.985, 1.001]; 6/20 | 7.845 [7.751, 7.900]; 20/20 **strict win** | 6.892 [6.836, 6.956]; 20/20 **strict win** |
| maccs/first_use | 1.002 [0.996, 1.005]; 11/20 | 1.882 [1.869, 1.887]; 20/20 **strict win** | 2.174 [2.166, 2.196]; 20/20 **strict win** |
| maccs/repeated_use | 0.999 [0.996, 1.004]; 8/20 | 1.940 [1.934, 1.950]; 20/20 **strict win** | 2.247 [2.229, 2.259]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.001 [0.991, 1.007]; 10/20 | 2.119 [2.104, 2.135]; 20/20 **strict win** | 2.390 [2.372, 2.402]; 20/20 **strict win** |

## Initial matrix: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.989 [0.977, 0.995]; 3/20 | 0.424 [0.409, 0.427]; 0/20 | 0.112 [0.110, 0.114]; 0/20 |
| tpsa/repeated_use | 1.008 [0.995, 1.022]; 12/20 | 9.981 [9.752, 10.031]; 20/20 **strict win** | 1.312 [1.276, 1.339]; 20/20 **strict win** |
| tpsa/parse_inclusive | 0.995 [0.987, 1.003]; 8/20 | 5.049 [5.001, 5.061]; 20/20 **strict win** | 5.222 [5.194, 5.269]; 20/20 **strict win** |
| labute_asa/first_use | 0.990 [0.984, 0.992]; 3/20 | 0.204 [0.196, 0.207]; 0/20 | 0.084 [0.081, 0.088]; 0/20 |
| labute_asa/repeated_use | 1.001 [0.989, 1.014]; 10/20 | 2.647 [2.593, 2.669]; 20/20 **strict win** | 1.171 [1.127, 1.245]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.999 [0.991, 1.003]; 8/20 | 4.806 [4.772, 4.845]; 20/20 **strict win** | 5.157 [5.126, 5.192]; 20/20 **strict win** |
| num_rings/first_use | 0.946 [0.934, 0.953]; 1/20 | 0.096 [0.092, 0.099]; 0/20 | 0.073 [0.071, 0.077]; 0/20 |
| num_rings/repeated_use | 0.988 [0.978, 0.999]; 6/20 | 1.689 [1.659, 1.743]; 20/20 **strict win** | 1.055 [1.025, 1.090]; 17/20 |
| num_rings/parse_inclusive | 0.972 [0.966, 0.995]; 5/20 | 11.545 [11.368, 11.577]; 20/20 **strict win** | 12.576 [12.437, 12.705]; 20/20 **strict win** |
| morgan2_2048/first_use | 0.989 [0.985, 0.997]; 4/20 | 1.545 [1.537, 1.559]; 20/20 **strict win** | 0.990 [0.981, 1.000]; 7/20 |
| morgan2_2048/repeated_use | 0.999 [0.986, 1.004]; 8/20 | 4.735 [4.640, 4.773]; 20/20 **strict win** | 2.881 [2.837, 2.892]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.993 [0.985, 0.998]; 5/20 | 4.573 [4.550, 4.599]; 20/20 **strict win** | 4.417 [4.397, 4.460]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.995 [0.988, 0.997]; 5/20 | 1.589 [1.574, 1.599]; 20/20 **strict win** | 1.012 [1.006, 1.020]; 17/20 |
| morgan2_chiral/repeated_use | 1.000 [0.997, 1.004]; 9/20 | 4.761 [4.721, 4.811]; 20/20 **strict win** | 2.891 [2.881, 2.910]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.989 [0.983, 0.992]; 4/20 | 4.559 [4.538, 4.580]; 20/20 **strict win** | 4.432 [4.415, 4.453]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.992 [0.986, 0.996]; 2/20 | 1.775 [1.763, 1.789]; 20/20 **strict win** | 1.157 [1.146, 1.170]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.995 [0.993, 1.007]; 9/20 | 4.553 [4.537, 4.581]; 20/20 **strict win** | 2.898 [2.882, 2.941]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.993 [0.989, 0.995]; 3/20 | 4.469 [4.451, 4.496]; 20/20 **strict win** | 4.291 [4.267, 4.337]; 20/20 **strict win** |
| maccs/first_use | 1.001 [0.996, 1.007]; 11/20 | 1.640 [1.636, 1.644]; 20/20 **strict win** | 1.753 [1.747, 1.760]; 20/20 **strict win** |
| maccs/repeated_use | 0.998 [0.995, 1.002]; 8/20 | 1.732 [1.729, 1.739]; 20/20 **strict win** | 1.854 [1.845, 1.859]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.003 [0.999, 1.004]; 14/20 | 1.805 [1.799, 1.810]; 20/20 **strict win** | 1.932 [1.927, 1.940]; 20/20 **strict win** |

## Adaptive negative observations

The plan was saved before timing: any initial baseline-relative lane with CI upper < 1 selects its entire operation for one separate 20-block follow-up in all modes. Initial observations remain above; follow-ups are not pooled or substituted. New follow-up negatives remain recorded without recursive retries.

```json
{
  "5000": {
    "negative_initial_lanes": {
      "tpsa/first_use": {
        "median": 0.989,
        "ci95": [
          0.977,
          0.995
        ],
        "winning_blocks": 3,
        "blocks": 20
      },
      "labute_asa/first_use": {
        "median": 0.99,
        "ci95": [
          0.984,
          0.992
        ],
        "winning_blocks": 3,
        "blocks": 20
      },
      "num_rings/first_use": {
        "median": 0.946,
        "ci95": [
          0.934,
          0.953
        ],
        "winning_blocks": 1,
        "blocks": 20
      },
      "num_rings/repeated_use": {
        "median": 0.988,
        "ci95": [
          0.978,
          0.999
        ],
        "winning_blocks": 6,
        "blocks": 20
      },
      "num_rings/parse_inclusive": {
        "median": 0.972,
        "ci95": [
          0.966,
          0.995
        ],
        "winning_blocks": 5,
        "blocks": 20
      },
      "morgan2_2048/first_use": {
        "median": 0.989,
        "ci95": [
          0.985,
          0.997
        ],
        "winning_blocks": 4,
        "blocks": 20
      },
      "morgan2_2048/parse_inclusive": {
        "median": 0.993,
        "ci95": [
          0.985,
          0.998
        ],
        "winning_blocks": 5,
        "blocks": 20
      },
      "morgan2_chiral/first_use": {
        "median": 0.995,
        "ci95": [
          0.988,
          0.997
        ],
        "winning_blocks": 5,
        "blocks": 20
      },
      "morgan2_chiral/parse_inclusive": {
        "median": 0.989,
        "ci95": [
          0.983,
          0.992
        ],
        "winning_blocks": 4,
        "blocks": 20
      },
      "morgan3_2048/first_use": {
        "median": 0.992,
        "ci95": [
          0.986,
          0.996
        ],
        "winning_blocks": 2,
        "blocks": 20
      },
      "morgan3_2048/parse_inclusive": {
        "median": 0.993,
        "ci95": [
          0.989,
          0.995
        ],
        "winning_blocks": 3,
        "blocks": 20
      }
    },
    "operations": [
      "tpsa",
      "labute_asa",
      "num_rings",
      "morgan2_2048",
      "morgan2_chiral",
      "morgan3_2048"
    ],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  },
  "0": {
    "negative_initial_lanes": {
      "tpsa/first_use": {
        "median": 0.914,
        "ci95": [
          0.909,
          0.922
        ],
        "winning_blocks": 2,
        "blocks": 20
      },
      "tpsa/parse_inclusive": {
        "median": 0.965,
        "ci95": [
          0.955,
          0.975
        ],
        "winning_blocks": 3,
        "blocks": 20
      },
      "labute_asa/first_use": {
        "median": 0.917,
        "ci95": [
          0.903,
          0.928
        ],
        "winning_blocks": 0,
        "blocks": 20
      },
      "labute_asa/parse_inclusive": {
        "median": 0.967,
        "ci95": [
          0.959,
          0.982
        ],
        "winning_blocks": 4,
        "blocks": 20
      },
      "num_rings/first_use": {
        "median": 0.939,
        "ci95": [
          0.935,
          0.953
        ],
        "winning_blocks": 3,
        "blocks": 20
      },
      "num_rings/parse_inclusive": {
        "median": 0.978,
        "ci95": [
          0.971,
          0.987
        ],
        "winning_blocks": 4,
        "blocks": 20
      },
      "morgan2_2048/first_use": {
        "median": 0.978,
        "ci95": [
          0.966,
          0.987
        ],
        "winning_blocks": 3,
        "blocks": 20
      },
      "morgan2_2048/parse_inclusive": {
        "median": 0.993,
        "ci95": [
          0.986,
          0.998
        ],
        "winning_blocks": 6,
        "blocks": 20
      },
      "morgan3_2048/first_use": {
        "median": 0.984,
        "ci95": [
          0.97,
          0.993
        ],
        "winning_blocks": 5,
        "blocks": 20
      }
    },
    "operations": [
      "tpsa",
      "labute_asa",
      "num_rings",
      "morgan2_2048",
      "morgan3_2048"
    ],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  }
}
```

### Separate follow-up: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.912 [0.909, 0.918]; 1/20 | 1.333 [1.304, 1.353]; 20/20 **strict win** | 0.461 [0.457, 0.466]; 0/20 |
| tpsa/repeated_use | 1.007 [0.984, 1.017]; 12/20 | 9.866 [9.765, 10.119]; 20/20 **strict win** | 1.683 [1.636, 1.818]; 20/20 **strict win** |
| tpsa/parse_inclusive | 0.964 [0.955, 0.971]; 1/20 | 12.393 [12.330, 12.440]; 20/20 **strict win** | 11.645 [11.581, 11.780]; 20/20 **strict win** |
| labute_asa/first_use | 0.916 [0.912, 0.921]; 1/20 | 0.641 [0.631, 0.657]; 0/20 | 0.308 [0.302, 0.313]; 0/20 |
| labute_asa/repeated_use | 1.001 [0.983, 1.019]; 11/20 | 2.684 [2.629, 2.700]; 20/20 **strict win** | 1.361 [1.325, 1.431]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.965 [0.956, 0.970]; 0/20 | 11.687 [11.617, 11.752]; 20/20 **strict win** | 11.179 [11.122, 11.231]; 20/20 **strict win** |
| num_rings/first_use | 0.939 [0.935, 0.946]; 0/20 | 0.035 [0.034, 0.035]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 1.005 [0.996, 1.014]; 12/20 | 1.735 [1.686, 1.764]; 20/20 **strict win** | unranked (output mismatch) |
| num_rings/parse_inclusive | 0.970 [0.964, 0.977]; 0/20 | 11.065 [10.984, 11.197]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 0.974 [0.966, 0.978]; 1/20 | 3.045 [3.025, 3.097]; 20/20 **strict win** | 1.985 [1.973, 1.995]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.997 [0.989, 1.001]; 8/20 | 5.223 [5.189, 5.320]; 20/20 **strict win** | 3.242 [3.198, 3.271]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.981 [0.978, 0.985]; 1/20 | 8.451 [8.371, 8.483]; 20/20 **strict win** | 7.428 [7.397, 7.473]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.976 [0.975, 0.980]; 2/20 | 3.351 [3.333, 3.366]; 20/20 **strict win** | 2.204 [2.195, 2.217]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.993 [0.990, 0.998]; 5/20 | 4.878 [4.847, 4.929]; 20/20 **strict win** | 3.117 [3.099, 3.141]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.986 [0.980, 0.989]; 3/20 | 7.868 [7.827, 7.890]; 20/20 **strict win** | 6.914 [6.865, 6.951]; 20/20 **strict win** |

### Separate follow-up: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.983 [0.980, 0.992]; 2/20 | 0.419 [0.418, 0.427]; 0/20 | 0.115 [0.112, 0.119]; 0/20 |
| tpsa/repeated_use | 1.008 [0.984, 1.022]; 11/20 | 9.880 [9.628, 9.959]; 20/20 **strict win** | 1.332 [1.249, 1.429]; 19/20 |
| tpsa/parse_inclusive | 0.995 [0.989, 1.005]; 8/20 | 4.952 [4.912, 5.001]; 20/20 **strict win** | 5.158 [5.115, 5.186]; 20/20 **strict win** |
| labute_asa/first_use | 0.988 [0.978, 0.996]; 5/20 | 0.201 [0.198, 0.204]; 0/20 | 0.084 [0.080, 0.088]; 0/20 |
| labute_asa/repeated_use | 0.991 [0.983, 1.012]; 9/20 | 2.608 [2.556, 2.692]; 20/20 **strict win** | 1.253 [1.178, 1.349]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.990 [0.984, 0.997]; 4/20 | 4.741 [4.706, 4.754]; 20/20 **strict win** | 5.075 [5.058, 5.130]; 20/20 **strict win** |
| num_rings/first_use | 0.937 [0.931, 0.942]; 0/20 | 0.096 [0.093, 0.100]; 0/20 | 0.077 [0.068, 0.084]; 0/20 |
| num_rings/repeated_use | 0.997 [0.979, 1.009]; 8/20 | 1.690 [1.643, 1.760]; 20/20 **strict win** | 1.057 [1.028, 1.087]; 18/20 |
| num_rings/parse_inclusive | 0.977 [0.972, 0.989]; 3/20 | 11.224 [11.169, 11.357]; 20/20 **strict win** | 12.396 [12.251, 12.546]; 20/20 **strict win** |
| morgan2_2048/first_use | 0.990 [0.986, 0.996]; 4/20 | 1.544 [1.536, 1.559]; 20/20 **strict win** | 0.987 [0.978, 1.004]; 7/20 |
| morgan2_2048/repeated_use | 0.992 [0.984, 0.999]; 5/20 | 4.702 [4.648, 4.752]; 20/20 **strict win** | 2.873 [2.840, 2.901]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.987 [0.979, 0.997]; 6/20 | 4.522 [4.459, 4.546]; 20/20 **strict win** | 4.393 [4.349, 4.431]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.992 [0.986, 1.005]; 6/20 | 1.583 [1.567, 1.601]; 20/20 **strict win** | 1.019 [0.996, 1.026]; 13/20 |
| morgan2_chiral/repeated_use | 0.995 [0.991, 1.000]; 7/20 | 4.757 [4.717, 4.783]; 20/20 **strict win** | 2.924 [2.867, 2.949]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.995 [0.983, 1.010]; 10/20 | 4.514 [4.493, 4.544]; 20/20 **strict win** | 4.369 [4.335, 4.389]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.994 [0.991, 1.005]; 8/20 | 1.769 [1.754, 1.773]; 20/20 **strict win** | 1.164 [1.151, 1.174]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.000 [0.994, 1.011]; 10/20 | 4.518 [4.495, 4.551]; 20/20 **strict win** | 2.901 [2.890, 2.960]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.988 [0.985, 1.000]; 6/20 | 4.409 [4.373, 4.445]; 20/20 **strict win** | 4.254 [4.229, 4.289]; 20/20 **strict win** |

## Reproduction and status

[The evidence packet](2026-10-15-python-ring-scratch-artifacts/README.md) contains both measured wheels, every raw matrix/follow-up, exact compressed output dumps, build/test logs, the original implementation, the patch and all frozen source/input files. The archive is read back and all 10059 candidate file hashes verified; baseline overrides and removal of the added test reconstruct all 10057 baseline hashes. The patch also reconstructs the changed files in an independent tiny Git repository. Hashes verify source reconstruction, not bit-reproducibility of rebuilt wheels.

Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared metadata was not refreshed; this is source-version 1.1.4 evidence with no hosted CI, merge or publication. Previous unresolved regressions against older candidates, including chiral Morgan parse-inclusive, ring/Morgan repeated-use and MACCS observations, are not resolved by an immediate-baseline comparison. Separate run win counts must not be treated as attributable code improvements.
