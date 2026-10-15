# Python aromatic normalization storage — 2026-10-15

The production change is adopted in the attached source candidate. It is unmerged and does not establish complete superiority over RDKit or COSMolKit, or verify current main.

CPU sampling identified aromatic normalization as a material part of cold descriptors. The candidate reserves known atom/bond capacities in both normalization builders and reuses the original per-atom implicit-H vector as an indexed correction map. This removes a second correction-vector allocation and per-atom linear searches. The implicit-H calls, aromaticity model, bond-order decisions, atom/bond enumeration and metadata-copy sequence remain unchanged, including historical behavior in the patched rebuild. Two pre-existing formatting fragments in the same included source file are normalized. No public API, dependency or unsafe code is added.

The immediate baseline includes all eight previously accepted source candidates: cycle-root pruning, V3000 writer buffers, shared Kekulé descriptor views, V3000 parser token reuse, ring-component counting, the no-double-bond E/Z skip, aromaticity component reuse and V3000 positive-zero parsing. It excludes all rejected experiments. Both Python wheels were built freshly before/after this single production change; this avoids attributing inherited wins to the new optimization.

Decision: NCI first-use TPSA is 1.063x [1.054,1.084] and Labute ASA 1.059x [1.050,1.069] against the fresh eight-candidate baseline, both faster in all 20 blocks. NCI parse-inclusive medians improve 1.032x and 1.033x with intervals above 1, but only 18/20 blocks faster, so these are not strict lane wins. ChEMBL TPSA/ASA improvements are not established. No initial baseline-relative CI is wholly below 1 on either profile; the predeclared adaptive selection is empty. All native/Python gates and exact descriptor/fingerprint dumps passed. Complete superiority is not achieved: NCI first-use TPSA is only 0.428x RDKit and 0.119x COSMolKit, ASA 0.206x and 0.088x respectively. This is the ninth unmerged source candidate; no current main, hosted CI or publication is verified. Older baseline-relative unresolved negatives remain unresolved.

All 14,072 fixed inputs have byte-identical descriptor and Morgan dumps: 112,576 descriptor outcomes (first/repeated, exact float bits), 56,288 Morgan outputs (R2/R3 with both chirality settings), zero parse/evaluation errors. The added differential test freezes the original normalization finish function and compares every one of the 10,000 existing development inputs under both Hückel and RDKit-like models, on parsed graphs and their successful Kekulé forms. It compares atoms, bonds, atom tags, R-group labels, adjacency order, implicit H, stereo groups with bond memberships, neighbor order, bond-direction markers/anchors and SMILES closure metadata. It asserts that hydrogen corrections are exercised. Existing native/Python suites separately cover tags and other metadata. No new chemical input records were added.

Scoped release Rust tests: 1,219 passed, 0 failed, 2 existing ignored. Full Python suite: 1,389 passed, one skipped. Release all-target Clippy for perception and chemistry, formatting and the refactor budget passed. The complete logs are retained. The initial native-test build stopped on an incorrect include-module path; that failure log is retained with the corrected successful rerun. Only this run owned regenerable targets were removed between phases. Timing starts after all builds and validation finish.

## Diagnostic profile

Four separate 8-second workers repeatedly parse all 1,000 rows and evaluate TPSA or Labute ASA, with fresh descriptor state. macOS sample records five seconds at a one-millisecond requested interval. These are diagnostic CPU samples, not benchmark timings or attribution of wall time. NCI Labute ASA has 2,063 of 3,601 main-thread samples under descriptor aromaticity and 585 under its normalized molecule build, summed across nonrecursive call sites; full trees and worker checksums are retained. The no-N/O input audit uses pinned RDKit over all unfiltered rows and is retained in profile-summary.json. No chemical rows were added or filtered.

## Protocol

Pinned measured environment: RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, CPython 3.13.6, Rust 1.97.0, maturin 1.13.3, macOS arm64. The fixed, unfiltered slices contain 1,000 inputs each, offsets 0 and 5,000. Each uses 20 rotating fresh-process blocks across baseline, candidate and the two comparators. This is an exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.

Ratios are reference elapsed time / candidate elapsed time; greater than 1 means faster. Brackets are paired median-ratio 95% bootstrap intervals. Strict wins require semantic agreement, CI lower > 1 and all 20 blocks faster. Every operation, mode, block, error and mismatch is retained. Unchanged-operation movement is not evidence of a new algorithm.

First-use parses fresh objects outside the timer and may compare different eager/lazy work; repeated-use warms the operation once per object; parse-inclusive times actual fresh parsing and evaluation. Native return objects remain retained during timing; normalization and semantic checks occur outside timing. Whole-process peak RSS includes imports and preparation and is not per-operation allocation. Source, corpus, wheels and installed extension hashes are verified.

Initial comparator strict wins: first_use: 14/27, repeated_use: 24/27, parse_inclusive: 27/27. COSMolKit produces six unequal ring counts at offset 0, so all three corresponding modes remain unranked. The full remaining descriptor deficits are retained below.

## Initial matrix: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.006 [0.991, 1.029]; 11/20 | 1.477 [1.451, 1.505]; 20/20 **strict win** | 0.518 [0.512, 0.533]; 0/20 |
| tpsa/repeated_use | 0.994 [0.972, 1.013]; 8/20 | 9.764 [9.569, 9.973]; 20/20 **strict win** | 1.780 [1.670, 1.833]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.006 [0.995, 1.009]; 12/20 | 13.228 [13.063, 13.361]; 20/20 **strict win** | 12.337 [12.181, 12.430]; 20/20 **strict win** |
| labute_asa/first_use | 1.005 [0.976, 1.014]; 12/20 | 0.714 [0.699, 0.721]; 0/20 | 0.341 [0.331, 0.352]; 0/20 |
| labute_asa/repeated_use | 1.009 [0.976, 1.022]; 11/20 | 2.700 [2.569, 2.806]; 20/20 **strict win** | 1.387 [1.288, 1.422]; 19/20 |
| labute_asa/parse_inclusive | 1.003 [0.998, 1.014]; 12/20 | 12.537 [12.366, 12.638]; 20/20 **strict win** | 11.791 [11.635, 11.830]; 20/20 **strict win** |
| num_rings/first_use | 1.005 [0.992, 1.016]; 11/20 | 0.037 [0.035, 0.039]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 1.012 [0.989, 1.034]; 12/20 | 1.714 [1.664, 1.768]; 20/20 **strict win** | unranked (output mismatch) |
| num_rings/parse_inclusive | 1.008 [0.993, 1.013]; 13/20 | 11.714 [11.554, 11.856]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 1.008 [0.994, 1.017]; 12/20 | 3.119 [3.061, 3.163]; 20/20 **strict win** | 2.044 [2.008, 2.050]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.999 [0.994, 1.002]; 8/20 | 5.194 [5.107, 5.292]; 20/20 **strict win** | 3.201 [3.131, 3.236]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.998 [0.987, 1.008]; 9/20 | 8.671 [8.511, 8.722]; 20/20 **strict win** | 7.575 [7.520, 7.639]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.002 [0.993, 1.008]; 11/20 | 1.248 [1.236, 1.269]; 20/20 **strict win** | 0.822 [0.807, 0.836]; 0/20 |
| morgan2_chiral/repeated_use | 0.992 [0.980, 1.000]; 5/20 | 5.002 [4.941, 5.045]; 20/20 **strict win** | 3.107 [3.070, 3.197]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.001 [0.994, 1.004]; 10/20 | 4.485 [4.447, 4.511]; 20/20 **strict win** | 3.939 [3.897, 3.950]; 20/20 **strict win** |
| morgan3_2048/first_use | 0.997 [0.992, 1.004]; 8/20 | 3.349 [3.287, 3.394]; 20/20 **strict win** | 2.232 [2.196, 2.249]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 0.999 [0.997, 1.001]; 8/20 | 4.828 [4.736, 4.899]; 20/20 **strict win** | 3.122 [3.065, 3.140]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.999 [0.986, 1.004]; 9/20 | 7.927 [7.871, 8.067]; 20/20 **strict win** | 6.952 [6.921, 6.990]; 20/20 **strict win** |
| maccs/first_use | 1.004 [0.998, 1.007]; 13/20 | 1.885 [1.871, 1.893]; 20/20 **strict win** | 2.169 [2.159, 2.177]; 20/20 **strict win** |
| maccs/repeated_use | 1.001 [0.997, 1.005]; 12/20 | 1.942 [1.933, 1.955]; 20/20 **strict win** | 2.237 [2.214, 2.246]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.004 [1.002, 1.008]; 15/20 | 2.110 [2.086, 2.122]; 20/20 **strict win** | 2.385 [2.369, 2.391]; 20/20 **strict win** |

## Initial matrix: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.063 [1.054, 1.084]; 20/20 **strict win** | 0.428 [0.423, 0.444]; 0/20 | 0.119 [0.114, 0.126]; 0/20 |
| tpsa/repeated_use | 1.007 [0.973, 1.057]; 13/20 | 9.719 [9.490, 10.125]; 20/20 **strict win** | 1.363 [1.283, 1.542]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.032 [1.016, 1.050]; 18/20 | 5.077 [5.044, 5.112]; 20/20 **strict win** | 5.259 [5.217, 5.287]; 20/20 **strict win** |
| labute_asa/first_use | 1.059 [1.050, 1.069]; 20/20 **strict win** | 0.206 [0.201, 0.215]; 0/20 | 0.088 [0.085, 0.094]; 0/20 |
| labute_asa/repeated_use | 1.014 [0.990, 1.046]; 11/20 | 2.595 [2.561, 2.648]; 20/20 **strict win** | 1.203 [1.165, 1.384]; 19/20 |
| labute_asa/parse_inclusive | 1.033 [1.018, 1.037]; 18/20 | 4.836 [4.808, 4.870]; 20/20 **strict win** | 5.205 [5.165, 5.238]; 20/20 **strict win** |
| num_rings/first_use | 1.010 [0.995, 1.024]; 12/20 | 0.100 [0.098, 0.106]; 0/20 | 0.084 [0.078, 0.107]; 0/20 |
| num_rings/repeated_use | 0.997 [0.982, 1.017]; 7/20 | 1.692 [1.669, 1.730]; 20/20 **strict win** | 1.075 [1.039, 1.105]; 17/20 |
| num_rings/parse_inclusive | 1.000 [0.987, 1.015]; 10/20 | 11.830 [11.746, 11.875]; 20/20 **strict win** | 13.043 [12.962, 13.094]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.037 [1.033, 1.053]; 19/20 | 1.557 [1.537, 1.568]; 20/20 **strict win** | 1.002 [0.984, 1.023]; 10/20 |
| morgan2_2048/repeated_use | 1.005 [0.999, 1.010]; 13/20 | 4.711 [4.660, 4.755]; 20/20 **strict win** | 2.890 [2.852, 2.930]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.024 [1.020, 1.037]; 19/20 | 4.578 [4.558, 4.635]; 20/20 **strict win** | 4.469 [4.432, 4.515]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.039 [1.027, 1.052]; 17/20 | 1.592 [1.559, 1.611]; 20/20 **strict win** | 1.035 [1.020, 1.049]; 17/20 |
| morgan2_chiral/repeated_use | 1.000 [0.995, 1.012]; 10/20 | 4.751 [4.675, 4.802]; 20/20 **strict win** | 2.906 [2.893, 2.953]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.021 [1.015, 1.040]; 18/20 | 4.576 [4.544, 4.633]; 20/20 **strict win** | 4.470 [4.420, 4.497]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.039 [1.029, 1.054]; 18/20 | 1.785 [1.765, 1.798]; 20/20 **strict win** | 1.178 [1.167, 1.204]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.002 [0.994, 1.011]; 12/20 | 4.549 [4.475, 4.581]; 20/20 **strict win** | 2.900 [2.877, 2.933]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.025 [1.012, 1.029]; 17/20 | 4.478 [4.440, 4.526]; 20/20 **strict win** | 4.297 [4.271, 4.326]; 20/20 **strict win** |
| maccs/first_use | 1.008 [1.005, 1.015]; 18/20 | 1.642 [1.638, 1.650]; 20/20 **strict win** | 1.759 [1.752, 1.776]; 20/20 **strict win** |
| maccs/repeated_use | 1.013 [1.005, 1.021]; 17/20 | 1.738 [1.729, 1.752]; 20/20 **strict win** | 1.867 [1.853, 1.872]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.010 [1.008, 1.018]; 18/20 | 1.813 [1.806, 1.817]; 20/20 **strict win** | 1.943 [1.936, 1.961]; 20/20 **strict win** |

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

[The evidence packet](2026-10-15-python-normalization-storage-artifacts/README.md) contains both measured wheels, every raw matrix/follow-up, exact compressed output dumps, build/test logs, the original implementation, the patch and all frozen source/input files. The archive is read back and all 10054 candidate file hashes verified; baseline overrides and removal of the added test reconstruct all 10053 baseline hashes. The patch also reconstructs the changed files in an independent tiny Git repository. Hashes verify source reconstruction, not bit-reproducibility of rebuilt wheels.

Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared metadata was not refreshed; this is source-version 1.1.4 evidence with no hosted CI, merge or publication. Previous unresolved regressions against older candidates, including chiral Morgan parse-inclusive, ring/Morgan repeated-use and MACCS observations, are not resolved by an immediate-baseline comparison. Separate run win counts must not be treated as attributable code improvements.
