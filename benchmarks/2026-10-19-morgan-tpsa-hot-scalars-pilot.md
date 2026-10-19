# Numeric cache placement: rejected pilot — 2026-10-19

**Rejected before formal measurement.** This third candidate combines direct Morgan hashing, ring-sensitive TPSA typing and a private `repr(C)` layout that moves the ring atomic before the inline type-erased array. It preserves the same scalar offsets and 528-byte cache size, but the eight-block diagnostic showed TPSA repeated-use 0.928339526 and ring repeated-use 0.975118018. Both miss the unchanged 0.99 source median floor. It did not advance to a 24-block formal experiment and is not adopted.

All eight complete blocks, sixteen fresh-process arms, all 21 diagnostic lane medians, native hashes and output accounting remain in the [numeric packet](2026-10-19-morgan-tpsa-hot-scalars-pilot-numeric.json.gz) and [CSV](2026-10-19-morgan-tpsa-hot-scalars-pilot-lanes.csv). No pilot sample is reused in another epoch. These are diagnostic observations, not paired 24-block confidence intervals or external rankings.

The fixed full 10,000-row corpus and seven operations follow the [preceding rejected combined experiment](2026-10-19-morgan-tpsa-environments.md). Baseline/candidate orders alternate AB/BA four times each, with PYTHONHASHSEED=20261019. The primary plan required improvements in chiral Morgan first-use, TPSA first-use and ring repeated-use plus all 21 source medians at least 0.99. The original baseline wheel is reused as an immutable artifact; all elapsed samples are fresh.

Candidate source `fd96f67c`; baseline source `00596fb1ea5397524f98e3789c41f17ee8616fdc`. Each native module originated in an independent empty target, release opt3/full LTO/one codegen unit/native-inchi. Loaded extension hashes match their wheel builds. All 56,288 fingerprint cells match source baseline and RDKit; all 42,216 scalar value/refusal states remain source-exact; all 14,072 TPSA values match RDKit within 1e-9. Native tests: 989 chemistry, 163 core, 352 fingerprint passes, one existing chemistry ignore; Clippy passes.

| Diagnostic lane | Source median ratio |
|---|---:|
| `tpsa/first_use` | 1.036478615 |
| `tpsa/repeated_use` | 0.928339526 |
| `tpsa/parse_inclusive` | 1.002289999 |
| `labute_asa/first_use` | 0.986189808 |
| `labute_asa/repeated_use` | 0.997309864 |
| `labute_asa/parse_inclusive` | 1.002335536 |
| `num_rings/first_use` | 1.186164216 |
| `num_rings/repeated_use` | 0.975118018 |
| `num_rings/parse_inclusive` | 0.994507441 |
| `morgan2_2048/first_use` | 1.036347769 |
| `morgan2_2048/repeated_use` | 1.067568548 |
| `morgan2_2048/parse_inclusive` | 1.013671067 |
| `morgan2_chiral/first_use` | 1.015730319 |
| `morgan2_chiral/repeated_use` | 1.051212045 |
| `morgan2_chiral/parse_inclusive` | 1.008025683 |
| `morgan3_2048/first_use` | 1.028217081 |
| `morgan3_2048/repeated_use` | 1.037983421 |
| `morgan3_2048/parse_inclusive` | 1.009152595 |
| `maccs/first_use` | 1.004068585 |
| `maccs/repeated_use` | 1.002967759 |
| `maccs/parse_inclusive` | 1.001870966 |

The next distinct candidate separates the type-erased payload allocation from inline numeric cache fields. It needs its own fresh build, correctness audit, pilot and formal source gates. This rejected layout-only experiment remains separate.

A later [combined source candidate](2026-10-19-python-ring-cache.md) passes its own registered gates after additional runtime changes. This earlier record remains rejected; its elapsed observations are not pooled into the later result.
