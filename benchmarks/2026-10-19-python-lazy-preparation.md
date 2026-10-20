# Python SMILES preparation boundary

Date: 2026-10-19. This source experiment compares the published macOS arm64
CPython 3.13 v1.1.8 wheel with a candidate built from the same revision after
removing descriptor, aromatic-view and ring preparation from `from_smiles()`.
The fixed input is rows 0–999 of the 10k browser corpus. Each arm runs in a
fresh process in rotating order for 20 paired blocks. Candidate and baseline
outputs are identical for every measured operation.

## Result

Parse-only median falls from **5.715 ms to 3.185 ms** for 1,000 molecules
(1.79x). Peak process RSS is effectively unchanged: 55.77 versus 55.89 MiB.

The table reports baseline time divided by candidate time. Values above one
favor the lazy candidate.

| Work boundary | TPSA | Labute ASA | QED | Ring count | Morgan R2 |
|---|---:|---:|---:|---:|---:|
| Parse + one operation | **1.329x** | **1.318x** | **1.010x** | **1.155x** | **1.203x** |
| First operation after parsing | 0.209x | 0.251x | 0.993x | 0.025x | 0.749x |
| Same operation after preparation | 1.004x | 1.001x | 0.963x | 0.999x | 1.002x |

Every parse-inclusive lane wins 19 or 20 of 20 blocks; TPSA, ASA, ring count
and Morgan win all 20 with a paired 95% lower bound above 1.0. The first-use
row is not a regression in total work: it shows that preparation moved out of
parse and into the first consumer. It is nevertheless a real API-boundary
trade-off and is published explicitly.

For two descriptors on the same molecule, parse + TPSA + ASA improves
**1.290x** with 20/20 wins. A prepared 20-molecule, 20,000-call reuse lane is
statistically unchanged for TPSA (1.003x; CI includes 1.0) and 0.992x for ASA
(95% CI 0.984–0.998). The candidate therefore helps one-shot and small-bundle
workloads without claiming a hot-loop speedup.

## Evidence boundary

- [single-operation raw data](2026-10-19-python-lazy-preparation-frontier.json)
- [two-descriptor raw data](2026-10-19-python-lazy-preparation-bundle.json)
- [prepared-reuse raw data](2026-10-19-python-lazy-preparation-reuse.json)
- [source and artifact manifest](2026-10-19-python-lazy-preparation-source-manifest.json)

The candidate wheel required four bytes of `__LINKEDIT` alignment padding and
local ad-hoc re-signing because the host linker emitted a misaligned string
table. Executable code was unchanged, and the repaired wheel hash is recorded.
This is one-host source evidence, not a published-package or cross-platform
claim. A hosted clean-wheel rerun remains required before a release claim.

