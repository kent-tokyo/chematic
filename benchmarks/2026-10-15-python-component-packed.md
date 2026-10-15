# Python contiguous aromatic-component scratch — 2026-10-15

The production change is adopted in the attached source candidate. It is unmerged and does not establish complete superiority over RDKit or COSMolKit, or verify current main.

The candidate combines three private u32 component-fact vectors into one allocation with adjacent root, cyclic-edge-count and atom-count slices. Their total element count stays n + 2c; the kind vector and public cached flags/components are unchanged. Traversal, component numbering, donor/candidate/Hückel predicates and arithmetic are unchanged. No public API, dependency or unsafe code is added. Both preceding FIFO and inline SmallVec trials are rejected, archived separately and excluded from this candidate.

The immediate baseline includes all ten previously accepted source candidates, including aromatic normalization storage and the V3000 bounded-capacity/zero-map writer changes. All rejected experiments remain excluded. The exact baseline Python wheel is reused after matching all 10,057 baseline source/input hashes and its wheel/installed-extension hashes. It was built before either trial; the candidate is built freshly. New timing observations are separate and never pooled. The origin is recorded in baseline-artifact-origin.json.

Decision: Adopt the contiguous private ComponentInfo allocation as a modest source-candidate improvement. ChEMBL first-use TPSA 1.019x [1.012,1.032] and ASA 1.018x [1.013,1.028]; uninterrupted NCI TPSA 1.010x [1.006,1.017] and ASA 1.011x [1.004,1.016]. Each cold lane wins 17/20, so none is labeled a strict cold win. NCI parse-inclusive TPSA 1.026x (19/20), ASA 1.029x (19/20), ring count 1.048x [1.039,1.056] (20/20, strict). No initial source-relative CI is wholly below 1, so the predeclared adaptive policy selects no operation. Cold gains are small and interpretation is limited by the documented ChEMBL ENOSPC checkpoint interruption; all 16 saved original blocks remain unchanged, the unavailable block is disclosed, and only missing checkpoint blocks are recollected. Exact 14,072-input parity, 1,220 native tests, 1,389 Python tests, formatting, Clippy and unchanged refactor budget pass. Three allocations become one with the same n+2c u32 payload. This does not establish complete superiority: NCI first-use ASA reference/candidate ratios are RDKit 0.206x and COSMolKit 0.086x; eager/lazy preparation differs. Prior unresolved regressions against older candidates are unchanged. Both rejected predecessor experiments remain excluded.

All 14,072 fixed inputs have byte-identical descriptor and Morgan dumps: 112,576 descriptor outcomes (first/repeated, exact float bits), 56,288 Morgan outputs (R2/R3 with both chirality settings), zero parse/evaluation errors. A differential test freezes the original ComponentInfo constructor. All 14,072 inputs compare exact root labels, component kinds, edge/atom counts, cyclic flags, booleans and root-iteration order on repeated constructor calls. The frozen reference constructor populates shared ring caches before the candidate constructor; fresh/cold versus repeated descriptor output parity is separately checked by the two independent Python dump processes. No new chemical input records were added.

Scoped release Rust tests: 1,220 passed, 0 failed, 2 existing ignored. Full Python suite: 1,389 passed, one skipped. Release all-target Clippy for perception and chemistry, formatting and the refactor budget passed. The complete logs are retained. The reused baseline build originally hit ENOSPC; its failure and successful retry logs are retained with the baseline origin record. Successful current candidate validation logs are retained. Only this run owned regenerable targets were removed between phases. Timing starts after all builds and validation finish.

## Measurement interruption

ChEMBL initial collection stopped after 16 saved blocks when opening the block-17 checkpoint failed with transient ENOSPC. The just-computed four-arm block 17 was not persisted and its timings could not be recovered. Both the 16-block incomplete matrix and error logs are retained. The existing runner --resume revalidated the exact source, wheels, corpus, protocol and saved output hashes in untimed fresh workers, then recollected block 17 and continued through block 20. The original 16 saved blocks are unchanged. NCI completed without interruption and was not repeated. This collection boundary limits interpretation of small differences; no performance-based sample was removed and no measurements from another candidate were pooled. See measurement-recovery.json and frontier-0-resume.log.

## Protocol

Pinned measured environment: RDKit 2026.09.1, COSMolKit 0.5.0-rc.22, CPython 3.13.6, Rust 1.97.0, maturin 1.13.3, macOS arm64. The fixed, unfiltered slices contain 1,000 inputs each, offsets 0 and 5,000. Each uses 20 rotating fresh-process blocks across baseline, candidate and the two comparators. This is an exposed development corpus on one host, without a sealed holdout or toolkit-wide claim.

Ratios are reference elapsed time / candidate elapsed time; greater than 1 means faster. Brackets are paired median-ratio 95% bootstrap intervals. Strict wins require semantic agreement, CI lower > 1 and all 20 blocks faster. Every successfully checkpointed operation, mode, block, error and mismatch is retained; the unavailable observation is disclosed below. Unchanged-operation movement is not evidence of a new algorithm.

First-use parses fresh objects outside the timer and may compare different eager/lazy work; repeated-use warms the operation once per object; parse-inclusive times actual fresh parsing and evaluation. Native return objects remain retained during timing; normalization and semantic checks occur outside timing. Whole-process peak RSS includes imports and preparation and is not per-operation allocation. Source, corpus, wheels and installed extension hashes are verified.

Initial comparator strict wins: first_use: 13/27, repeated_use: 26/27, parse_inclusive: 27/27. COSMolKit produces six unequal ring counts at offset 0, so all three corresponding modes remain unranked. The full remaining descriptor deficits are retained below.

## Initial matrix: offset 0

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.019 [1.012, 1.032]; 17/20 | 1.524 [1.490, 1.540]; 19/20 | 0.538 [0.516, 0.555]; 0/20 |
| tpsa/repeated_use | 1.005 [0.990, 1.050]; 11/20 | 9.903 [9.668, 10.275]; 20/20 **strict win** | 1.798 [1.699, 1.914]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.009 [1.005, 1.012]; 16/20 | 12.684 [12.604, 12.819]; 20/20 **strict win** | 11.861 [11.802, 11.943]; 20/20 **strict win** |
| labute_asa/first_use | 1.018 [1.013, 1.028]; 17/20 | 0.735 [0.716, 0.771]; 0/20 | 0.349 [0.342, 0.362]; 0/20 |
| labute_asa/repeated_use | 1.018 [1.001, 1.029]; 14/20 | 2.739 [2.655, 2.862]; 20/20 **strict win** | 1.458 [1.366, 1.521]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.003 [0.996, 1.012]; 10/20 | 11.977 [11.826, 12.098]; 20/20 **strict win** | 11.329 [11.263, 11.402]; 20/20 **strict win** |
| num_rings/first_use | 1.005 [0.997, 1.010]; 14/20 | 0.037 [0.036, 0.039]; 0/20 | unranked (output mismatch) |
| num_rings/repeated_use | 1.019 [0.998, 1.047]; 12/20 | 1.737 [1.713, 1.876]; 20/20 **strict win** | unranked (output mismatch) |
| num_rings/parse_inclusive | 1.004 [0.986, 1.020]; 12/20 | 11.124 [11.022, 11.306]; 20/20 **strict win** | unranked (output mismatch) |
| morgan2_2048/first_use | 1.010 [0.995, 1.033]; 13/20 | 3.150 [3.130, 3.170]; 20/20 **strict win** | 2.037 [2.027, 2.056]; 20/20 **strict win** |
| morgan2_2048/repeated_use | 0.999 [0.995, 1.011]; 9/20 | 5.254 [5.191, 5.293]; 20/20 **strict win** | 3.206 [3.201, 3.238]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.013 [1.007, 1.021]; 19/20 | 8.577 [8.482, 8.632]; 20/20 **strict win** | 7.519 [7.475, 7.543]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.016 [1.005, 1.032]; 16/20 | 1.276 [1.251, 1.285]; 20/20 **strict win** | 0.830 [0.823, 0.843]; 0/20 |
| morgan2_chiral/repeated_use | 1.002 [0.998, 1.020]; 12/20 | 5.016 [4.932, 5.079]; 20/20 **strict win** | 3.115 [3.068, 3.165]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.003 [0.987, 1.012]; 11/20 | 4.401 [4.338, 4.481]; 20/20 **strict win** | 3.874 [3.804, 3.920]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.005 [0.993, 1.019]; 13/20 | 3.378 [3.331, 3.440]; 20/20 **strict win** | 2.261 [2.213, 2.322]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.005 [0.996, 1.011]; 13/20 | 4.887 [4.825, 4.944]; 20/20 **strict win** | 3.145 [3.109, 3.161]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.004 [0.999, 1.008]; 13/20 | 7.906 [7.880, 8.001]; 20/20 **strict win** | 6.889 [6.809, 6.951]; 20/20 **strict win** |
| maccs/first_use | 1.001 [0.993, 1.006]; 11/20 | 1.888 [1.880, 1.899]; 20/20 **strict win** | 2.174 [2.168, 2.182]; 20/20 **strict win** |
| maccs/repeated_use | 1.000 [0.994, 1.007]; 9/20 | 1.947 [1.938, 1.959]; 20/20 **strict win** | 2.247 [2.238, 2.253]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.000 [0.997, 1.007]; 10/20 | 2.109 [2.100, 2.121]; 20/20 **strict win** | 2.368 [2.360, 2.405]; 20/20 **strict win** |

## Initial matrix: offset 5000

| Operation / preparation | Immediate baseline | RDKit | COSMolKit |
|---|---:|---:|---:|
| tpsa/first_use | 1.010 [1.006, 1.017]; 17/20 | 0.433 [0.430, 0.437]; 0/20 | 0.114 [0.113, 0.115]; 0/20 |
| tpsa/repeated_use | 1.001 [0.993, 1.017]; 10/20 | 9.941 [9.693, 10.061]; 20/20 **strict win** | 1.317 [1.305, 1.378]; 20/20 **strict win** |
| tpsa/parse_inclusive | 1.026 [1.018, 1.036]; 19/20 | 5.108 [5.063, 5.138]; 20/20 **strict win** | 5.267 [5.241, 5.304]; 20/20 **strict win** |
| labute_asa/first_use | 1.011 [1.004, 1.016]; 17/20 | 0.206 [0.202, 0.208]; 0/20 | 0.086 [0.083, 0.087]; 0/20 |
| labute_asa/repeated_use | 0.999 [0.991, 1.028]; 10/20 | 2.642 [2.590, 2.709]; 20/20 **strict win** | 1.154 [1.132, 1.177]; 20/20 **strict win** |
| labute_asa/parse_inclusive | 1.029 [1.024, 1.033]; 19/20 | 4.872 [4.822, 4.906]; 20/20 **strict win** | 5.207 [5.192, 5.242]; 20/20 **strict win** |
| num_rings/first_use | 0.998 [0.988, 1.002]; 9/20 | 0.099 [0.097, 0.102]; 0/20 | 0.079 [0.074, 0.082]; 0/20 |
| num_rings/repeated_use | 1.009 [0.997, 1.030]; 13/20 | 1.700 [1.670, 1.786]; 20/20 **strict win** | 1.068 [1.032, 1.099]; 19/20 |
| num_rings/parse_inclusive | 1.048 [1.039, 1.056]; 20/20 **strict win** | 11.881 [11.738, 11.910]; 20/20 **strict win** | 12.837 [12.804, 12.993]; 20/20 **strict win** |
| morgan2_2048/first_use | 1.006 [1.002, 1.014]; 14/20 | 1.553 [1.543, 1.574]; 20/20 **strict win** | 1.001 [0.992, 1.007]; 10/20 |
| morgan2_2048/repeated_use | 1.004 [1.000, 1.009]; 14/20 | 4.742 [4.689, 4.770]; 20/20 **strict win** | 2.885 [2.869, 2.936]; 20/20 **strict win** |
| morgan2_2048/parse_inclusive | 1.020 [1.018, 1.025]; 19/20 | 4.635 [4.627, 4.656]; 20/20 **strict win** | 4.499 [4.473, 4.521]; 20/20 **strict win** |
| morgan2_chiral/first_use | 1.011 [1.007, 1.019]; 16/20 | 1.596 [1.588, 1.613]; 20/20 **strict win** | 1.024 [1.019, 1.032]; 19/20 |
| morgan2_chiral/repeated_use | 1.007 [1.000, 1.014]; 14/20 | 4.739 [4.704, 4.778]; 20/20 **strict win** | 2.915 [2.884, 2.946]; 20/20 **strict win** |
| morgan2_chiral/parse_inclusive | 1.024 [1.015, 1.032]; 17/20 | 4.644 [4.629, 4.674]; 20/20 **strict win** | 4.504 [4.465, 4.539]; 20/20 **strict win** |
| morgan3_2048/first_use | 1.008 [0.996, 1.015]; 13/20 | 1.790 [1.785, 1.796]; 20/20 **strict win** | 1.176 [1.167, 1.181]; 20/20 **strict win** |
| morgan3_2048/repeated_use | 1.002 [0.998, 1.010]; 12/20 | 4.550 [4.524, 4.570]; 20/20 **strict win** | 2.903 [2.890, 2.920]; 20/20 **strict win** |
| morgan3_2048/parse_inclusive | 1.016 [1.012, 1.022]; 19/20 | 4.538 [4.515, 4.569]; 20/20 **strict win** | 4.323 [4.293, 4.347]; 20/20 **strict win** |
| maccs/first_use | 1.002 [0.999, 1.004]; 13/20 | 1.643 [1.636, 1.648]; 20/20 **strict win** | 1.755 [1.752, 1.764]; 20/20 **strict win** |
| maccs/repeated_use | 1.003 [0.999, 1.006]; 13/20 | 1.734 [1.731, 1.738]; 20/20 **strict win** | 1.860 [1.853, 1.865]; 20/20 **strict win** |
| maccs/parse_inclusive | 1.001 [0.997, 1.005]; 11/20 | 1.800 [1.798, 1.810]; 20/20 **strict win** | 1.933 [1.930, 1.939]; 20/20 **strict win** |

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

[The evidence packet](2026-10-15-python-component-packed-artifacts/README.md) contains both measured wheels, every raw matrix/follow-up, exact compressed output dumps, build/test logs, the original implementation, the patch and all frozen source/input files. The archive is read back and all 10058 candidate file hashes verified; baseline overrides and removal of the added test reconstruct all 10057 baseline hashes. The patch also reconstructs the changed files in an independent tiny Git repository. Hashes verify source reconstruction, not bit-reproducibility of rebuilt wheels.

Last verified Git base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared metadata was not refreshed; this is source-version 1.1.4 evidence with no hosted CI, merge or publication. Previous unresolved regressions against older candidates, including chiral Morgan parse-inclusive, ring/Morgan repeated-use and MACCS observations, are not resolved by an immediate-baseline comparison. Separate run win counts must not be treated as attributable code improvements.
