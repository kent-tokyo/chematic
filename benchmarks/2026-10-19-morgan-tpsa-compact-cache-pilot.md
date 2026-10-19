# Compact core cache: rejected pilot — 2026-10-19

**Rejected before formal measurement.** The fourth source candidate combines direct Morgan hashing, ring-sensitive TPSA typing and a lazily allocated boxed type-erased cache. Inline numeric/cache storage shrinks from 528 to 88 bytes on arm64. The eight-block pilot improves first-use TPSA 1.175433203x, ASA 1.204863732x and ring count 1.586285207x, but ring repeated-use remains 0.997178657 and ring parse-inclusive is 0.989811712. The latter misses the unchanged 0.99 source median floor. The candidate did not advance to formal measurement; favorable pilot medians do not establish external superiority.

All eight complete diagnostic blocks and 21 lane medians remain in the [numeric packet](2026-10-19-morgan-tpsa-compact-cache-pilot-numeric.json.gz) and [CSV](2026-10-19-morgan-tpsa-compact-cache-pilot-lanes.csv). Baseline/candidate orders alternate AB/BA four times each in fresh processes, all 10,000 normal-corpus rows and seven operations, PYTHONHASHSEED=20261019. No elapsed sample is reused in any other epoch. The original baseline wheel is an immutable artifact reused with fresh elapsed observations.

Source baseline `00596fb1ea5397524f98e3789c41f17ee8616fdc`; candidate `d6b87215`. Release opt3/full LTO/one codegen unit/native-inchi; candidate native build starts in an independent empty target. All 56,288 fingerprint cells equal source baseline and RDKit; all 42,216 scalar value/refusal states equal source baseline; all 14,072 TPSA values agree with RDKit within 1e-9. Native tests: chemistry 989 passes / 1 existing ignore, core 165 passes, fingerprints 352 passes; modified-library Clippy passes.

The [layout-only pilot](2026-10-19-morgan-tpsa-hot-scalars-pilot.md) and both preceding full-24-block rejected experiments remain separate. The next distinct candidate adds Python-wrapper ring-count memoization for its immutable graph, seeded from the already existing bounded parser preparation. It requires fresh build, correctness audit, pilot and formal gates.

| Diagnostic lane | Source median ratio |
|---|---:|
| `tpsa/first_use` | 1.175433203 |
| `tpsa/repeated_use` | 1.213741335 |
| `tpsa/parse_inclusive` | 0.999789468 |
| `labute_asa/first_use` | 1.204863732 |
| `labute_asa/repeated_use` | 1.670394965 |
| `labute_asa/parse_inclusive` | 1.000722425 |
| `num_rings/first_use` | 1.586285207 |
| `num_rings/repeated_use` | 0.997178657 |
| `num_rings/parse_inclusive` | 0.989811712 |
| `morgan2_2048/first_use` | 1.030271826 |
| `morgan2_2048/repeated_use` | 1.055775683 |
| `morgan2_2048/parse_inclusive` | 1.010238818 |
| `morgan2_chiral/first_use` | 1.015921618 |
| `morgan2_chiral/repeated_use` | 1.061057775 |
| `morgan2_chiral/parse_inclusive` | 1.011041629 |
| `morgan3_2048/first_use` | 1.032573853 |
| `morgan3_2048/repeated_use` | 1.037846195 |
| `morgan3_2048/parse_inclusive` | 1.011559638 |
| `maccs/first_use` | 0.994690647 |
| `maccs/repeated_use` | 0.995232004 |
| `maccs/parse_inclusive` | 0.996018844 |

A later [combined source candidate](2026-10-19-python-ring-cache.md) passes its own registered gates after additional runtime changes. This earlier record remains rejected; its elapsed observations are not pooled into the later result.
