# Python selected-component roots — 2026-10-15

The production experiment is rejected for adoption; only its report and evidence are delivered. This is unmerged source-version 1.1.4 evidence, without current-main verification or hosted CI. Complete superiority over RDKit or COSMolKit is not established.

Decision: Reject the selected-component-root prototype: none of the four intended first-use TPSA/ASA lanes demonstrates an improvement. NCI TPSA 1.006x [0.992,1.012] (12/20 faster), ASA 0.997x [0.978,1.008] (9/20); ChEMBL TPSA 0.995x [0.981,0.999] (6/20), ASA 0.999x [0.989,1.005] (9/20). All 42 initial source-relative lanes lack a strict win. The predeclared separate follow-up reproduces slower ChEMBL first-use TPSA: 0.991x [0.986,0.998], 5/20 faster; NCI parse-inclusive MACCS also remains below baseline: 0.997x [0.993,0.999], 4/20. No pooling, substitutions or recursive retries were used. MACCS timing movement in an unchanged path is not evidence of a new algorithm; binary layout and host noise may contribute. Exact ring sequences/order match in 168864 differential comparisons; all 14072-input descriptor and Morgan outputs, 1221 native tests and 1389 Python tests pass. Correctness does not justify a speed optimization without benefit. Existing twelve accepted source candidates remain unchanged. Complete superiority is still not established: NCI first-use ASA reference/candidate ratios are RDKit 0.211x and COSMolKit 0.088x; preparation differs, and COSMolKit ring-count mismatches remain unranked.

## Change and correctness

The candidate reuses the final roots from the existing ordered union pass. It counts non-root atoms and each final root once, then counts selected edges by direct root lookup. This removes a temporary ring-atom flag allocation and avoids repeating root searches for every selected edge. Singleton selected self-edges retain their existing counting behavior. Root labels, selected-edge membership, component iteration order, cycle candidate generation, rank/tie selection and arithmetic are unchanged. The ordered and single-cycle fast paths are unchanged. No public API, dependency or unsafe code is added.

The immediate baseline includes all twelve previously accepted source candidates, including packed private component facts and direct V3000 integer output. Rejected prototypes are excluded. Both current wheels were built freshly; the baseline was built before edits and all 10,059 baseline source/input hashes were checked. No previous wheel or timing was reused.

A frozen copy of the entire original component-basis function checks exact selected ring sequences and order on 14,072 pinned inputs. Three whole-component masks (all, even component labels, odd component labels), both rank settings and cold/repeated calls give 168,864 differential comparisons. Independent baseline/candidate processes additionally produce byte-identical full descriptor and Morgan dumps: 112,576 descriptor outcomes including exact floating bits and 56,288 Morgan outputs, with zero parse/evaluation errors. Existing licensed input corpora are reused; no new molecular records were added.

Scoped release Rust tests: 1,221 passed, 0 failed, 2 existing ignored. Python: 1,389 passed, one skipped. Release all-target Clippy for chemistry/perception, formatting and the unchanged refactor budget (487 production files, 8,063 functions, 18 file allowances, 88 function allowances) passed. Full logs are retained. Builds and CPU sampling completed before formal timing; only owned completed build targets were removed.

Fresh baseline CPU samples are retained in sample-*.txt and getter-profile-review.json. Counts are filtered to the descriptor getter ancestry; inclusive function counts overlap and are not additive. Sampling identifies a plausible target, not a predicted speedup.

## Protocol

Pinned environment: CPython 3.13.6, RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, Rust 1.97.0, maturin 1.13.3, macOS arm64. Two fixed unfiltered 1,000-input slices, offsets 5,000 then 0, use 20 rotating fresh-process four-arm blocks, seven operations and all three preparation modes. This is one exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.

Ratios are reference elapsed / candidate elapsed; greater than 1 is faster. Brackets are paired median-ratio 95% bootstrap intervals. A strict win requires semantic agreement, CI lower > 1 and all 20 blocks faster. Output mismatches suppress ranking; all observations and failures are retained. Whole-process RSS includes imports and preparation.

First-use prepares fresh parsed objects outside timing and can compare different eager/lazy preparation. Repeated-use warms the operation once on each object; parse-inclusive times fresh parsing and evaluation. Return objects are retained during timing; normalization and semantic checks occur outside timing. Every block validates output determinism and wheel/extension identity.

An external checkpoint wrapper leaves frozen benchmark workers, timers, output checks and selection policy unchanged. It reserves 16 MiB of owned storage and on ENOSPC releases that reserve and retries writing the same computed checkpoint for at most 60 seconds. It never removes or remeasures a performance observation. Wrapper logs retain any such event. This run completed all four matrices without checkpoint interruption or reserve release.

Initial comparator strict wins: first_use: 14/27, repeated_use: 24/27, parse_inclusive: 27/27. COSMolKit has six unequal ring counts at offset 0; the corresponding modes remain unranked. Every remaining deficit is shown below.

## Initial matrix: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.995 [0.981, 0.999]; 6/20 | 1.510 [1.484, 1.535]; 20/20 **strict win** | 0.524 [0.512, 0.536]; 0/20 |
| tpsa/repeated_use | 0.994 [0.984, 1.004]; 7/20 | 9.885 [9.772, 10.014]; 20/20 **strict win** | 1.743 [1.662, 1.804]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.003 [0.991, 1.010]; 12/20 | 13.114 [13.000, 13.220]; 20/20 **strict win** | 12.212 [12.143, 12.276]; 20/20 **strict win** |
| labute_asa/first_use | 0.999 [0.989, 1.005]; 9/20 | 0.731 [0.720, 0.749]; 0/20 | 0.343 [0.341, 0.346]; 0/20 |
| labute_asa/repeated_use | 0.989 [0.984, 1.008]; 7/20 | 2.715 [2.687, 2.743]; 20/20 **strict win** | 1.456 [1.334, 1.534]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.001 [0.993, 1.011]; 11/20 | 12.388 [12.297, 12.522]; 20/20 **strict win** | 11.700 [11.666, 11.781]; 20/20 **strict win** |
| num_rings/first_use | 1.000 [0.995, 1.007]; 10/20 | 0.039 [0.036, 0.040]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 0.996 [0.984, 1.003]; 6/20 | 1.779 [1.750, 1.812]; 20/20 **strict win** | unranked (output mismatch) |
| num_rings/parse_inclusive | 1.004 [0.996, 1.011]; 13/20 | 11.636 [11.560, 11.723]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 1.004 [0.999, 1.011]; 13/20 | 3.173 [3.157, 3.184]; 20/20 **strict win** | 2.049 [2.038, 2.075]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.999 [0.991, 1.006]; 10/20 | 5.229 [5.175, 5.286]; 20/20 **strict win** | 3.209 [3.186, 3.235]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.000 [0.991, 1.008]; 10/20 | 8.675 [8.647, 8.719]; 20/20 **strict win** | 7.597 [7.570, 7.645]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.003 [0.997, 1.005]; 12/20 | 1.273 [1.269, 1.280]; 20/20 **strict win** | 0.832 [0.828, 0.835]; 0/20 |
| morgan2_chiral/repeated_use | 1.002 [0.992, 1.010]; 10/20 | 5.037 [4.991, 5.079]; 20/20 **strict win** | 3.114 [3.098, 3.145]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.001 [0.997, 1.005]; 11/20 | 4.497 [4.482, 4.513]; 20/20 **strict win** | 3.950 [3.937, 3.966]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.996 [0.992, 1.002]; 6/20 | 3.381 [3.356, 3.402]; 20/20 **strict win** | 2.249 [2.224, 2.261]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.001 [0.994, 1.006]; 12/20 | 4.890 [4.841, 4.921]; 20/20 **strict win** | 3.118 [3.105, 3.141]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.002 [0.996, 1.007]; 12/20 | 8.000 [7.962, 8.031]; 20/20 **strict win** | 7.001 [6.927, 7.031]; 20/20 **strict win** |
| maccs/first_use | 1.000 [0.996, 1.004]; 10/20 | 1.875 [1.869, 1.882]; 20/20 **strict win** | 2.163 [2.157, 2.169]; 20/20 **strict win** |
| maccs/repeated_use | 1.002 [0.998, 1.003]; 13/20 | 1.939 [1.933, 1.947]; 20/20 **strict win** | 2.233 [2.222, 2.239]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.999 [0.996, 1.004]; 8/20 | 2.100 [2.097, 2.105]; 20/20 **strict win** | 2.371 [2.361, 2.376]; 20/20 **strict win** |

## Initial matrix: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.006 [0.992, 1.012]; 12/20 | 0.445 [0.431, 0.460]; 0/20 | 0.120 [0.117, 0.123]; 0/20 |
| tpsa/repeated_use | 0.990 [0.978, 1.024]; 8/20 | 9.616 [9.366, 9.913]; 20/20 **strict win** | 1.331 [1.226, 1.383]; 18/20 |
| tpsa/parse_inclusive | 1.003 [0.987, 1.014]; 11/20 | 5.173 [5.107, 5.252]; 20/20 **strict win** | 5.278 [5.225, 5.315]; 20/20 **strict win** |
| labute_asa/first_use | 0.997 [0.978, 1.008]; 9/20 | 0.211 [0.206, 0.224]; 0/20 | 0.088 [0.084, 0.091]; 0/20 |
| labute_asa/repeated_use | 1.001 [0.981, 1.011]; 11/20 | 2.605 [2.547, 2.642]; 20/20 **strict win** | 1.166 [1.108, 1.246]; 19/20 |
| labute_asa/parse_inclusive | 1.000 [0.995, 1.012]; 10/20 | 4.961 [4.923, 4.977]; 20/20 **strict win** | 5.243 [5.192, 5.277]; 20/20 **strict win** |
| num_rings/first_use | 0.993 [0.974, 1.010]; 8/20 | 0.100 [0.098, 0.106]; 0/20 | 0.077 [0.075, 0.085]; 0/20 |
| num_rings/repeated_use | 0.995 [0.985, 1.010]; 9/20 | 1.701 [1.668, 1.733]; 20/20 **strict win** | 1.073 [1.038, 1.121]; 16/20 |
| num_rings/parse_inclusive | 1.003 [0.995, 1.018]; 13/20 | 11.989 [11.861, 12.097]; 20/20 **strict win** | 13.051 [12.883, 13.151]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.003 [0.997, 1.009]; 12/20 | 1.557 [1.538, 1.577]; 20/20 **strict win** | 1.008 [0.994, 1.032]; 12/20 |
| morgan2_2048/repeated_use | 0.993 [0.986, 1.008]; 7/20 | 4.693 [4.619, 4.733]; 20/20 **strict win** | 2.867 [2.808, 2.943]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.003 [0.989, 1.022]; 10/20 | 4.687 [4.641, 4.724]; 20/20 **strict win** | 4.516 [4.434, 4.544]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.000 [0.989, 1.008]; 10/20 | 1.604 [1.597, 1.628]; 20/20 **strict win** | 1.041 [1.026, 1.056]; 17/20 |
| morgan2_chiral/repeated_use | 1.000 [0.990, 1.014]; 9/20 | 4.774 [4.736, 4.831]; 20/20 **strict win** | 2.897 [2.877, 2.996]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.001 [0.992, 1.009]; 11/20 | 4.652 [4.616, 4.717]; 20/20 **strict win** | 4.512 [4.438, 4.535]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.995 [0.976, 1.008]; 8/20 | 1.771 [1.757, 1.801]; 20/20 **strict win** | 1.181 [1.155, 1.197]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.001 [0.978, 1.021]; 10/20 | 4.486 [4.435, 4.560]; 20/20 **strict win** | 2.889 [2.848, 2.939]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.997 [0.993, 1.008]; 8/20 | 4.549 [4.488, 4.601]; 20/20 **strict win** | 4.288 [4.252, 4.343]; 20/20 **strict win** |
| maccs/first_use | 0.997 [0.993, 1.005]; 9/20 | 1.635 [1.622, 1.645]; 20/20 **strict win** | 1.744 [1.737, 1.752]; 20/20 **strict win** |
| maccs/repeated_use | 0.996 [0.983, 1.003]; 7/20 | 1.719 [1.705, 1.727]; 20/20 **strict win** | 1.846 [1.835, 1.857]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.989 [0.975, 0.998]; 5/20 | 1.793 [1.772, 1.806]; 20/20 **strict win** | 1.908 [1.900, 1.937]; 20/20 **strict win** |

## Adaptive observations

The plan was saved before timing: an initial source-relative CI upper < 1 selects its entire operation for one separate 20-block follow-up in all modes. Initial and follow-up observations remain separate, without pooling or recursive retries. All new follow-up negatives remain recorded.

```json
{
  "5000": {
    "negative_initial_lanes": {
      "maccs/parse_inclusive": {
        "median": 0.989,
        "ci95": [
          0.975,
          0.998
        ],
        "winning_blocks": 5,
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
      "tpsa/first_use": {
        "median": 0.995,
        "ci95": [
          0.981,
          0.999
        ],
        "winning_blocks": 6,
        "blocks": 20
      }
    },
    "operations": [
      "tpsa"
    ],
    "policy": "One separate adaptive follow-up of 20 blocks. Retain all three modes and all outcomes. Do not pool or iterate."
  }
}
```

### Separate follow-up: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 0.991 [0.986, 0.998]; 5/20 | 1.501 [1.477, 1.519]; 20/20 **strict win** | 0.527 [0.516, 0.532]; 0/20 |
| tpsa/repeated_use | 0.998 [0.982, 1.009]; 9/20 | 9.907 [9.750, 10.068]; 20/20 **strict win** | 1.705 [1.624, 1.781]; 20/20 **strict win** |
| tpsa/parse_inclusive | 0.996 [0.989, 1.003]; 6/20 | 12.931 [12.839, 13.052]; 20/20 **strict win** | 12.129 [12.052, 12.169]; 20/20 **strict win** |

### Separate follow-up: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| maccs/first_use | 0.997 [0.995, 1.002]; 7/20 | 1.630 [1.625, 1.636]; 20/20 **strict win** | 1.740 [1.735, 1.750]; 20/20 **strict win** |
| maccs/repeated_use | 0.995 [0.993, 1.002]; 7/20 | 1.725 [1.719, 1.730]; 20/20 **strict win** | 1.845 [1.841, 1.849]; 20/20 **strict win** |
| maccs/parse_inclusive | 0.997 [0.993, 0.999]; 4/20 | 1.794 [1.792, 1.802]; 20/20 **strict win** | 1.923 [1.919, 1.932]; 20/20 **strict win** |

## Reproduction and status

[The evidence packet](2026-10-15-python-basis-roots-artifacts/README.md) contains both measured wheels, every initial/adaptive matrix, compressed exact output dumps, build/test and profile logs, original implementation, exact patch, and all frozen source/input files. Archive readback checks all 10060 candidate hashes; baseline overrides and removal of the added test reconstruct all 10059 baseline hashes. The patch independently reconstructs changed files in a tiny Git repository. Hashes demonstrate source reconstruction, not bit-reproducibility of rebuilt wheels.

Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared Git metadata was not refreshed; this is not merged, published or hosted-CI-validated. Prior unresolved regressions against older source candidates, including V3000 reading and Morgan/MACCS observations, are not resolved by an immediate-baseline comparison. OpenBabel format lanes are outside this descriptor experiment. Separate-run win counts do not establish attributable algorithm gains.
