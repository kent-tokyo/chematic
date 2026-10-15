# Python ring-count cyclic-component reuse — 2026-10-15

This is an unmerged source candidate. It does not establish complete superiority over RDKit or COSMolKit, or verify the current state of `main`. The immediate baseline includes all four preceding candidates: cycle-root pruning, V3000 writer buffers, shared Kekulé descriptor views and V3000 parser token storage. Both wheels were built freshly from that snapshot before and after this change.

Immediate-baseline ring first-use medians: offset 0: 1.181 [1.080, 1.367]; 17/20, offset 5000: 1.101 [1.042, 1.150]; 15/20. Neither meets the all-20-block strict-win criterion. Parse-inclusive baseline intervals and all other lanes are retained below. First-use NCI TPSA and Labute ASA remain slower than both comparators; those descriptor algorithms are not changed by this patch.

The ring bridge pass already computes cyclic-component labels. Ring counting now reuses those labels instead of constructing another union-find, scanning the cyclic edges twice and repeatedly finding roots. Query-bond and coordination-bond guards, single-closure proofs, symmetrized-SSSR selection and bounded macrocycle fallbacks are unchanged. `chematic_perception::sssr::RingComponents` and `ring_components_shared` become documented public APIs; the struct is non-exhaustive. Component labels are topology-specific, not persistent atom identifiers.

All 14,072 inputs retain byte-identical descriptor and Morgan dumps: 112,576 descriptor outcomes (first/repeated, exact float bits) and 56,288 Morgan outputs (R2/R3, both chirality settings), with no parse or evaluation failures. Rust: 955 chem + 257 perception tests passed, with two existing ignored tests; Python: 1,389 passed, one skipped. Added tests compare the frozen pre-change implementation on 10,000 molecules, every one of the 1,024 five-atom simple graphs, graph mutations, query/coordination bonds and valid externally seeded bond flags. Release Clippy for all targets, formatting and the refactor budget passed.

Ratios below are reference elapsed time / candidate elapsed time; greater than 1 is faster. Brackets are the paired median-ratio 95% bootstrap interval. Every row retains all 20 blocks; a strict win requires the lower bound > 1 and all 20 blocks faster. No operations or observations were removed. First-use excludes parsing and therefore includes different eager/lazy preparation costs; repeated-use has a warm per-object operation; parse-inclusive times actual SMILES parsing and evaluation together. Whole-process peak RSS includes imports and all operations and is not an operation-allocation measurement.

Pinned comparators: RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, CPython 3.13.6, macOS arm64. Two exposed development slices of 1,000 inputs each from the fixed 10,000-input corpus; this is one host and these seven operations, without a sealed holdout or toolkit-wide claim. Source/wheel/corpus hashes, toolchain, load averages, every rotating fresh-process arm, errors and output gates are retained.

Initial full-matrix strict comparator wins: first_use 5/27, repeated_use 17/27, parse_inclusive 26/27. COSMolKit has six unequal ring counts on offset 0, so that comparison is unranked in every preparation mode. Baseline-relative movement in unchanged operations is not attributed to a new algorithm.

## Initial matrix: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.029 [0.983, 1.262]; 13/20 | 1.424 [1.367, 1.488]; 17/20 | 0.467 [0.404, 0.570]; 1/20 |
| tpsa/repeated_use | 1.076 [0.971, 1.327]; 14/20 | 8.597 [7.558, 9.803]; 20/20 **strict win** | 1.784 [1.336, 2.177]; 17/20 |
| tpsa/parse_inclusive | 1.087 [0.989, 1.317]; 12/20 | 12.854 [12.173, 13.300]; 20/20 **strict win** | 12.284 [11.943, 15.643]; 20/20 **strict win** |
| labute_asa/first_use | 1.055 [1.011, 1.372]; 15/20 | 0.676 [0.640, 0.723]; 1/20 | 0.334 [0.305, 0.391]; 1/20 |
| labute_asa/repeated_use | 1.190 [1.030, 1.411]; 16/20 | 2.556 [2.355, 2.716]; 19/20 | 1.460 [1.318, 1.792]; 17/20 |
| labute_asa/parse_inclusive | 1.081 [1.002, 1.344]; 15/20 | 12.403 [11.390, 13.456]; 20/20 **strict win** | 11.870 [11.147, 15.903]; 20/20 **strict win** |
| num_rings/first_use | 1.181 [1.080, 1.367]; 17/20 | 0.036 [0.033, 0.039]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 1.089 [0.971, 1.462]; 13/20 | 1.550 [1.377, 1.634]; 19/20 | unranked (output mismatch) |
| num_rings/parse_inclusive | 1.101 [1.008, 1.346]; 15/20 | 11.388 [10.437, 12.004]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 1.074 [0.998, 1.308]; 13/20 | 3.091 [2.828, 3.215]; 20/20 **strict win** | 2.160 [1.902, 2.526]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 1.055 [1.004, 1.214]; 14/20 | 5.036 [4.677, 5.289]; 20/20 **strict win** | 3.340 [3.061, 3.895]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.064 [0.973, 1.213]; 13/20 | 8.531 [7.516, 9.270]; 20/20 **strict win** | 7.812 [7.254, 9.190]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.084 [0.996, 1.240]; 13/20 | 1.214 [1.134, 1.276]; 17/20 | 0.854 [0.796, 0.982]; 6/20 |
| morgan2_chiral/repeated_use | 1.075 [1.002, 1.161]; 14/20 | 4.966 [4.620, 5.434]; 20/20 **strict win** | 3.199 [2.998, 3.862]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.065 [1.002, 1.144]; 15/20 | 4.356 [4.174, 4.924]; 20/20 **strict win** | 4.017 [3.794, 4.791]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.047 [0.997, 1.169]; 13/20 | 3.311 [3.066, 3.573]; 20/20 **strict win** | 2.343 [2.200, 2.624]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.071 [1.014, 1.217]; 15/20 | 4.826 [4.395, 5.448]; 20/20 **strict win** | 3.392 [3.151, 3.731]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.046 [0.983, 1.177]; 12/20 | 7.900 [7.328, 8.093]; 20/20 **strict win** | 7.755 [6.929, 8.735]; 20/20 **strict win** |
| maccs/first_use | 1.054 [1.001, 1.128]; 14/20 | 1.877 [1.750, 2.041]; 19/20 | 2.409 [2.165, 2.630]; 19/20 |
| maccs/repeated_use | 1.046 [0.992, 1.083]; 13/20 | 2.008 [1.882, 2.414]; 20/20 **strict win** | 2.414 [2.233, 2.668]; 19/20 |
| maccs/parse_inclusive | 1.029 [0.978, 1.077]; 12/20 | 2.132 [2.028, 2.491]; 20/20 **strict win** | 2.416 [2.325, 2.656]; 19/20 |

## Initial matrix: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.001 [0.971, 1.041]; 10/20 | 0.398 [0.368, 0.409]; 0/20 | 0.102 [0.078, 0.107]; 0/20 |
| tpsa/repeated_use | 0.994 [0.839, 1.049]; 10/20 | 8.270 [7.336, 9.903]; 20/20 **strict win** | 1.079 [0.942, 1.325]; 11/20 |
| tpsa/parse_inclusive | 1.010 [0.931, 1.036]; 12/20 | 4.846 [4.371, 4.954]; 20/20 **strict win** | 4.871 [3.966, 5.107]; 20/20 **strict win** |
| labute_asa/first_use | 1.009 [0.962, 1.054]; 11/20 | 0.193 [0.172, 0.199]; 0/20 | 0.075 [0.057, 0.081]; 0/20 |
| labute_asa/repeated_use | 1.010 [0.899, 1.049]; 10/20 | 2.288 [1.956, 2.732]; 19/20 | 1.003 [0.848, 1.309]; 10/20 |
| labute_asa/parse_inclusive | 1.003 [0.933, 1.040]; 11/20 | 4.583 [4.012, 4.731]; 20/20 **strict win** | 4.932 [3.902, 5.092]; 20/20 **strict win** |
| num_rings/first_use | 1.101 [1.042, 1.150]; 15/20 | 0.103 [0.097, 0.106]; 0/20 | 0.082 [0.069, 0.098]; 0/20 |
| num_rings/repeated_use | 1.023 [0.963, 1.095]; 11/20 | 1.692 [1.477, 1.757]; 19/20 | 1.090 [0.860, 1.143]; 12/20 |
| num_rings/parse_inclusive | 1.009 [0.961, 1.047]; 11/20 | 11.754 [11.218, 12.015]; 20/20 **strict win** | 12.553 [10.872, 12.996]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.012 [0.962, 1.027]; 12/20 | 1.474 [1.298, 1.511]; 19/20 | 0.956 [0.760, 0.997]; 6/20 |
| morgan2_2048/repeated_use | 0.994 [0.916, 1.032]; 8/20 | 4.601 [4.325, 4.721]; 20/20 **strict win** | 2.787 [2.417, 3.031]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 0.984 [0.909, 1.026]; 9/20 | 4.453 [4.055, 4.603]; 20/20 **strict win** | 4.221 [3.449, 4.495]; 20/20 **strict win** |
| morgan2_chiral/first_use | 0.988 [0.916, 1.034]; 9/20 | 1.466 [1.250, 1.514]; 19/20 | 0.956 [0.781, 1.001]; 6/20 |
| morgan2_chiral/repeated_use | 0.992 [0.927, 1.013]; 9/20 | 4.623 [4.278, 4.712]; 20/20 **strict win** | 2.869 [2.380, 3.138]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 0.992 [0.917, 1.018]; 8/20 | 4.327 [3.915, 4.484]; 20/20 **strict win** | 4.195 [3.396, 4.424]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.001 [0.972, 1.040]; 10/20 | 1.649 [1.422, 1.701]; 20/20 **strict win** | 1.101 [0.888, 1.215]; 13/20 |
| morgan3_2048/repeated_use | 0.995 [0.957, 1.016]; 8/20 | 4.396 [3.895, 4.571]; 20/20 **strict win** | 2.814 [2.214, 2.906]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 0.999 [0.952, 1.040]; 10/20 | 4.247 [3.976, 4.398]; 20/20 **strict win** | 4.077 [3.161, 4.236]; 20/20 **strict win** |
| maccs/first_use | 0.991 [0.944, 1.003]; 8/20 | 1.599 [1.267, 1.631]; 18/20 | 1.718 [1.383, 1.759]; 18/20 |
| maccs/repeated_use | 0.994 [0.961, 1.023]; 8/20 | 1.656 [1.495, 1.716]; 20/20 **strict win** | 1.834 [1.642, 1.857]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.005 [0.979, 1.035]; 11/20 | 1.690 [1.590, 1.811]; 20/20 **strict win** | 1.921 [1.648, 1.982]; 20/20 **strict win** |

## Adaptive negative observations

Every initial baseline-relative lane whose reported 95% interval is wholly below 1 selects its operation for one separate 20-block follow-up, retaining all three modes. Initial observations remain above. Follow-ups are adaptive diagnostics, not replacements or pooled confirmation; newly appearing follow-up negatives are retained without recursive retries.

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

## Historical ring regression audit

This additional audit was declared before the current timing runs. Its baseline is the preceding pre-Kekulé-view-sharing wheel, already containing cycle-root pruning and V3000 writer buffers. The candidate also contains view sharing, parser token reuse and this ring change, so the audit measures their cumulative result rather than isolating one optimization. It examines the previously reproduced roughly 2% offset-5000 ring parse-inclusive regression; results are separate from the immediate baseline and adaptive observations.

| Operation / preparation | Historical baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| num_rings/first_use | 1.067 [1.019, 1.109]; 16/20 | 0.097 [0.090, 0.105]; 0/20 | 0.094 [0.086, 0.102]; 0/20 |
| num_rings/repeated_use | 1.006 [0.905, 1.053]; 11/20 | 1.492 [1.406, 1.811]; 18/20 | 0.942 [0.873, 1.149]; 9/20 |
| num_rings/parse_inclusive | 1.018 [1.004, 1.039]; 14/20 | 11.477 [10.552, 12.299]; 20/20 **strict win** | 12.583 [11.658, 13.453]; 20/20 **strict win** |

Historical parse-inclusive result against the pre-view-sharing checkpoint: 1.018 [1.004, 1.039]; 14/20. The separate audit does not establish a strict recovery of the earlier regression.

## Reproduction and status

All raw reports, logs, exact output dumps, measured wheels, reconstructable patch, baseline overrides and the full frozen build/input snapshot are in [the evidence packet](2026-10-15-python-ring-component-reuse-artifacts/README.md). The archive was read back and all 10,048 candidate hashes verified; applying baseline overrides and removing the added test reconstructs all 10,047 baseline hashes. The patch was applied in an independent tiny Git repository and verified against the candidate. Build logs preserve initial ENOSPC failures and successful retries; only working files belonging to these prior runs were deduplicated without changing their bytes or paths. Timing began after builds and validation completed.

Last verified shared Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared Git metadata remains inaccessible; the attached source snapshot and prior candidate files are guarded by hashes. This candidate has no hosted CI, merge, release or publication. Previous candidate regressions against older baselines, including the earlier cycle-root chiral-Morgan parse-inclusive observation, are not resolved by an incremental comparison unless explicitly tested here.
