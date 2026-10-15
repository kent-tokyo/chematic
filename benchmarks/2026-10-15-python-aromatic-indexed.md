# Python indexed aromatic normalization — 2026-10-15

The production experiment is rejected for adoption; only its report and evidence are delivered. This is unmerged source-version 1.1.4 evidence, without current-main verification or hosted CI. Complete superiority over RDKit or COSMolKit is not established.

Decision: Neither fixed profile establishes the intended cold TPSA/Labute ASA improvement: all four first-use source-relative intervals include 1. NCI median gains of 3.0% and 5.5% remain uncertain; ChEMBL medians show no gain. Initial and adaptive regressions are retained separately. Reject production adoption rather than infer a benefit from temporary hash-table removal.

## Change and correctness

The candidate stores temporary atom flags and bond orders in dense index-aligned arrays instead of a second atom hash set and bond hash map. Molecule iterators enumerate dense indices in order, including after removal. Model membership predicates, promote-only versus authoritative atom semantics, aromatic-bond endpoint fallback, traversal order, implicit-H preservation, direction stashing, metadata copying and the existing patched-rebuild behavior are unchanged. The model itself still uses its original sets; no ring or aromaticity algorithm is changed. No public API, dependency or unsafe code is added.

The immediate baseline includes all fourteen previously accepted source candidates, including the inactive E/Z guard and inline small-ring queue. Rejected prototypes are excluded. Both current wheels were built freshly; the baseline was built before edits and all 10,061 baseline source/input hashes were checked. No previous wheel or timing was reused.

The frozen-oracle proof confirms that the complete immediate-baseline builders, bond-order map and finish function differ only by rustfmt whitespace. This independent copy checks 14,072 pinned inputs in parsed and attempted Kekule forms (original clone on Kekule error), with both aromaticity algorithms and both promote-only/authoritative builders: 112,576 full graph/metadata comparisons. Additional empty, reordered-adjacency, stale-flag, arbitrary-model and metadata cases cover authoritative aromatic-bond endpoint fallback and ignored out-of-range model entries. The older normalization storage oracle remains intact; its call adapter uses the new production builder. Independent baseline/candidate processes additionally produce byte-identical full descriptor and Morgan dumps: 112,576 descriptor outcomes including exact floating bits and 56,288 Morgan outputs, with zero parse/evaluation errors. Existing licensed input corpora are reused; no new molecular records were added.

Scoped release Rust tests: 1,224 passed, 0 failed, 2 existing ignored. Python: 1,389 passed, one skipped. Release all-target Clippy for chemistry/perception, formatting and the unchanged refactor budget (487 production files, 8,063 functions, 18 file allowances, 88 function allowances) passed. Full logs are retained. Builds and CPU sampling completed before formal timing; only owned completed build targets were removed.

Fresh baseline CPU samples are retained in sample-*.txt and getter-profile-review.json. At offset 5,000, the normalization builder appears beneath 468/2,084 TPSA getter samples and 462/2,194 Labute ASA getter samples. At offset 0 it has no separately named samples. Counts are filtered to the descriptor getter ancestry; inclusive function counts overlap and are not additive, and LTO inlining can hide names. Sampling identifies a plausible target, not a predicted speedup. Initial sandbox sampling, the corrected shared-import compilation failure and a separate SIGTERM-interrupted validation compilation are retained. The interruption had no completed test results; the same source/tests resumed from Cargo fingerprints before any formal timing.

## Protocol

Pinned environment: CPython 3.13.6, RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, Rust 1.97.0, maturin 1.13.3, macOS arm64. Two fixed unfiltered 1,000-input slices, offsets 5,000 then 0, use 20 rotating fresh-process four-arm blocks, seven operations and all three preparation modes. This is one exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.

Ratios are reference elapsed / candidate elapsed; greater than 1 is faster. Brackets are paired median-ratio 95% bootstrap intervals. A strict win requires semantic agreement, CI lower > 1 and all 20 blocks faster. Output mismatches suppress ranking; all observations and failures are retained. Whole-process RSS includes imports and preparation.

First-use prepares fresh parsed objects outside timing and can compare different eager/lazy preparation. Repeated-use warms the operation once on each object; parse-inclusive times fresh parsing and evaluation. Return objects are retained during timing; normalization and semantic checks occur outside timing. Every block validates output determinism and wheel/extension identity.

An external checkpoint wrapper leaves frozen benchmark workers, timers, output checks and selection policy unchanged. It reserves 16 MiB of owned storage and on ENOSPC releases that reserve and retries writing the same computed checkpoint for at most 60 seconds. It never removes or remeasures a performance observation. Wrapper logs retain any such event. Matrix completion and any checkpoint events are recorded in the packet; all initial and selected adaptive matrices are retained.

Initial comparator strict wins: first_use: 13/27, repeated_use: 22/27, parse_inclusive: 26/27. COSMolKit has six unequal ring counts at offset 0; the corresponding modes remain unranked. Every remaining deficit is shown below.

## Initial matrix: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.999 [0.968, 1.008]; 9/20 | 1.496 [1.478, 1.506]; 20/20 **strict win** | 0.512 [0.500, 0.522]; 0/20 |
| tpsa/repeated_use | 1.006 [0.971, 1.019]; 11/20 | 9.547 [9.417, 9.848]; 20/20 **strict win** | 1.690 [1.637, 1.790]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.000 [0.977, 1.009]; 10/20 | 13.109 [12.784, 13.399]; 20/20 **strict win** | 12.282 [11.723, 12.358]; 20/20 **strict win** |
| labute_asa/first_use | 0.992 [0.951, 1.000]; 6/20 | 0.703 [0.683, 0.745]; 0/20 | 0.340 [0.327, 0.348]; 0/20 |
| labute_asa/repeated_use | 1.000 [0.970, 1.033]; 10/20 | 2.630 [2.451, 2.671]; 20/20 **strict win** | 1.371 [1.336, 1.424]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 0.986 [0.946, 1.003]; 7/20 | 12.121 [11.761, 12.550]; 20/20 **strict win** | 11.646 [11.043, 11.782]; 20/20 **strict win** |
| num_rings/first_use | 0.996 [0.959, 1.008]; 9/20 | 0.036 [0.035, 0.038]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 1.000 [0.980, 1.041]; 10/20 | 1.698 [1.605, 1.742]; 19/20 | unranked (output mismatch) |
| num_rings/parse_inclusive | 0.989 [0.945, 1.006]; 8/20 | 11.374 [11.067, 11.564]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 0.991 [0.952, 1.003]; 7/20 | 3.087 [2.963, 3.217]; 20/20 **strict win** | 1.977 [1.909, 2.019]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.994 [0.939, 1.030]; 8/20 | 5.117 [4.942, 5.336]; 20/20 **strict win** | 3.208 [3.069, 3.299]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.983 [0.954, 0.990]; 4/20 | 8.314 [8.085, 8.612]; 20/20 **strict win** | 7.504 [7.175, 7.690]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.996 [0.952, 1.006]; 8/20 | 1.235 [1.188, 1.275]; 20/20 **strict win** | 0.816 [0.785, 0.822]; 0/20 |
| morgan2_chiral/repeated_use | 0.995 [0.985, 1.000]; 6/20 | 4.877 [4.659, 5.016]; 20/20 **strict win** | 3.112 [3.079, 3.192]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.999 [0.994, 1.007]; 9/20 | 4.472 [4.442, 4.553]; 20/20 **strict win** | 3.912 [3.892, 3.949]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.996 [0.982, 1.027]; 7/20 | 3.291 [3.129, 3.356]; 20/20 **strict win** | 2.189 [2.172, 2.219]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.987 [0.963, 1.001]; 6/20 | 4.749 [4.486, 4.865]; 20/20 **strict win** | 3.094 [2.991, 3.130]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.992 [0.963, 1.010]; 7/20 | 7.886 [7.548, 8.000]; 20/20 **strict win** | 6.878 [6.543, 7.024]; 20/20 **strict win** |
| maccs/first_use | 1.000 [0.986, 1.019]; 10/20 | 1.861 [1.842, 1.935]; 20/20 **strict win** | 2.189 [2.171, 2.218]; 20/20 **strict win** |
| maccs/repeated_use | 1.001 [0.984, 1.016]; 10/20 | 1.913 [1.866, 1.959]; 20/20 **strict win** | 2.281 [2.266, 2.303]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.999 [0.970, 1.022]; 10/20 | 2.097 [2.045, 2.122]; 20/20 **strict win** | 2.405 [2.381, 2.434]; 20/20 **strict win** |

## Initial matrix: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.030 [0.983, 1.053]; 12/20 | 0.466 [0.456, 0.493]; 0/20 | 0.128 [0.120, 0.137]; 0/20 |
| tpsa/repeated_use | 0.949 [0.870, 1.083]; 8/20 | 9.550 [9.155, 10.517]; 20/20 **strict win** | 1.358 [1.293, 1.437]; 19/20 |
| tpsa/parse_inclusive | 1.017 [0.968, 1.039]; 12/20 | 5.284 [5.160, 5.637]; 20/20 **strict win** | 5.457 [5.171, 5.788]; 20/20 **strict win** |
| labute_asa/first_use | 1.055 [0.969, 1.068]; 13/20 | 0.237 [0.221, 0.247]; 1/20 | 0.094 [0.088, 0.105]; 0/20 |
| labute_asa/repeated_use | 1.018 [0.901, 1.087]; 12/20 | 2.545 [2.472, 2.669]; 20/20 **strict win** | 1.157 [1.030, 1.265]; 15/20 |
| labute_asa/parse_inclusive | 1.029 [0.977, 1.045]; 14/20 | 5.043 [4.923, 5.311]; 20/20 **strict win** | 5.379 [5.130, 5.694]; 20/20 **strict win** |
| num_rings/first_use | 1.007 [0.959, 1.053]; 11/20 | 0.109 [0.103, 0.112]; 0/20 | 0.105 [0.095, 0.111]; 0/20 |
| num_rings/repeated_use | 1.025 [0.952, 1.060]; 13/20 | 1.767 [1.640, 1.872]; 19/20 | 1.082 [1.016, 1.141]; 14/20 |
| num_rings/parse_inclusive | 0.997 [0.940, 1.018]; 9/20 | 11.875 [11.406, 12.776]; 20/20 **strict win** | 13.403 [12.301, 13.758]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.038 [0.980, 1.065]; 12/20 | 1.626 [1.593, 1.728]; 20/20 **strict win** | 1.069 [1.000, 1.108]; 14/20 |
| morgan2_2048/repeated_use | 1.004 [0.978, 1.016]; 12/20 | 4.901 [4.488, 5.013]; 20/20 **strict win** | 2.851 [2.750, 2.898]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.019 [0.959, 1.072]; 11/20 | 4.793 [4.286, 4.927]; 20/20 **strict win** | 4.584 [4.240, 4.770]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.031 [0.958, 1.078]; 12/20 | 1.691 [1.542, 1.782]; 20/20 **strict win** | 1.069 [1.015, 1.117]; 15/20 |
| morgan2_chiral/repeated_use | 1.006 [0.935, 1.055]; 11/20 | 4.691 [4.376, 4.919]; 20/20 **strict win** | 2.985 [2.751, 3.099]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.017 [0.977, 1.047]; 12/20 | 4.660 [4.459, 4.780]; 20/20 **strict win** | 4.521 [4.336, 4.763]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.028 [0.975, 1.040]; 13/20 | 1.799 [1.661, 1.864]; 20/20 **strict win** | 1.217 [1.160, 1.266]; 17/20 |
| morgan3_2048/repeated_use | 0.986 [0.927, 1.010]; 7/20 | 4.388 [4.044, 4.522]; 20/20 **strict win** | 2.826 [2.620, 2.978]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.014 [0.955, 1.045]; 12/20 | 4.547 [4.154, 4.689]; 20/20 **strict win** | 4.330 [4.129, 4.579]; 20/20 **strict win** |
| maccs/first_use | 0.985 [0.958, 1.002]; 6/20 | 1.590 [1.459, 1.633]; 20/20 **strict win** | 1.746 [1.671, 1.812]; 20/20 **strict win** |
| maccs/repeated_use | 1.001 [0.954, 1.018]; 10/20 | 1.686 [1.526, 1.717]; 20/20 **strict win** | 1.839 [1.627, 1.976]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.969 [0.923, 0.999]; 6/20 | 1.747 [1.504, 1.919]; 19/20 | 1.862 [1.767, 1.969]; 20/20 **strict win** |

## Adaptive observations

The plan was saved before timing: an initial source-relative CI upper < 1 selects its entire operation for one separate 20-block follow-up in all modes. Initial and follow-up observations remain separate, without pooling or recursive retries. All new follow-up negatives remain recorded.

```json
{
  "5000": {
    "negative_initial_lanes": {
      "maccs/parse_inclusive": {
        "median": 0.969,
        "ci95": [
          0.923,
          0.999
        ],
        "winning_blocks": 6,
        "blocks": 20
      }
    },
    "operations": [
      "maccs"
    ],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  },
  "0": {
    "negative_initial_lanes": {
      "morgan2_2048/parse_inclusive": {
        "median": 0.983,
        "ci95": [
          0.954,
          0.99
        ],
        "winning_blocks": 4,
        "blocks": 20
      }
    },
    "operations": [
      "morgan2_2048"
    ],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  }
}
```

### Separate follow-up: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| morgan2_2048/first_use | 0.995 [0.971, 1.001]; 7/20 | 3.040 [3.010, 3.121]; 19/20 | 1.986 [1.946, 2.022]; 19/20 |
| morgan2_2048/repeated_use | 1.001 [0.998, 1.014]; 11/20 | 5.197 [5.152, 5.260]; 20/20 **strict win** | 3.176 [3.128, 3.221]; 19/20 |
| morgan2_2048/parse_inclusive | 0.994 [0.971, 1.012]; 9/20 | 8.468 [8.297, 8.595]; 19/20 | 7.507 [7.296, 7.554]; 19/20 |

### Separate follow-up: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| maccs/first_use | 1.001 [0.995, 1.008]; 12/20 | 1.642 [1.637, 1.655]; 20/20 **strict win** | 1.755 [1.747, 1.769]; 20/20 **strict win** |
| maccs/repeated_use | 0.999 [0.992, 1.001]; 9/20 | 1.736 [1.728, 1.740]; 20/20 **strict win** | 1.852 [1.843, 1.854]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.003 [0.997, 1.006]; 13/20 | 1.808 [1.793, 1.813]; 20/20 **strict win** | 1.931 [1.912, 1.940]; 20/20 **strict win** |

## Reproduction and status

[The evidence packet](2026-10-15-python-aromatic-indexed-artifacts/README.md) contains both measured wheels, every initial/adaptive matrix, compressed exact output dumps, build/test and profile logs, original implementation, exact patch, and all frozen source/input files. Archive readback checks all 10062 candidate hashes; baseline overrides and removal of the added test reconstruct all 10061 baseline hashes. The patch independently reconstructs changed files in a tiny Git repository. Hashes demonstrate source reconstruction, not bit-reproducibility of rebuilt wheels.

Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared Git metadata was not refreshed; this is not merged, published or hosted-CI-validated. Prior unresolved regressions against older source candidates, including V3000 reading and Morgan/MACCS observations, are not resolved by an immediate-baseline comparison. OpenBabel format lanes are outside this descriptor experiment. Separate-run win counts do not establish attributable algorithm gains.
