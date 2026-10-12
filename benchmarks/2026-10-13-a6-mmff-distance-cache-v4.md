# A6 MMFF distance-cache holdout v4

Date: 2026-10-13

This record is the once-only outcome-blind v4 evaluation of candidate commit
`43619387545cc01f207d7bba8f69c06240a4a92b`. The candidate was frozen before
the manifest was selected. The result is useful quality evidence, but it does
not establish robust or universal speed superiority over RDKit.

## Provenance

- CheMatic source candidate: `43619387545cc01f207d7bba8f69c06240a4a92b`
- Candidate wheel: `chematic-1.1.1-cp313-cp313-macosx_11_0_arm64.whl`
- Candidate wheel SHA-256:
  `3683423dc908a2da2f04e95ac61b41e3a04b01de5de871b3855c01e417c5c381`
- Oracle: RDKit 2026.09.1
- Host: macOS arm64, Python 3.13.6
- Manifest:
  [`mmff94_stereo_safe_post_freeze_holdout_v4.json`](../validation/manifests/mmff94_stereo_safe_post_freeze_holdout_v4.json)
- Manifest SHA-256:
  `2c97a2e07eca4718a04af357d1883611699caf9cea2dd0b95e97a6fd4ab2afe9`
- Canonical corpus SHA-256:
  `e9e859b2de426470d667d60eb155330a1f2cd6f919d8dba68b620bf68f66c83f`
- Raw result:
  [`a6-mmff-distance-cache-holdout-v4-r3-2026-10-13.json`](../validation/results/a6-mmff-distance-cache-holdout-v4-r3-2026-10-13.json)
- Raw-result SHA-256:
  `9e5a52e5e6694d95b59ccab3d50c12dbcc07dd677eb72c439f068d56a57ddd17`

The 100 rows contain 50 structures with declared stereochemistry and 50
general drug-like structures. Selection excluded Tier A, Tier B, and inspected
holdouts v1-v3. CheMatic was not run while selecting the rows.

## Preregistered workload

- ETKDGv3, seed 42
- 20 candidates per molecule
- 200 MMFF94 iterations per initial candidate
- only when none converged, retry the lowest-energy non-converged candidate
  once for 100 iterations
- three alternating CheMatic/RDKit timing blocks

## Result

| Metric | CheMatic | RDKit |
|---|---:|---:|
| Seconds, block 1 | 29.118 | 29.451 |
| Seconds, block 2 | 29.912 | 30.506 |
| Seconds, block 3 | 30.274 | 30.026 |
| Median seconds | 29.912 | 30.026 |
| Generated candidates | 1,997 | 1,997 |
| Converged candidates | 912 | 913 |

The median RDKit/CheMatic ratio was `1.0038x`. CheMatic was faster in two of
three blocks, so this is an effective tie rather than a robust speed win.

| Quality metric | Result |
|---|---:|
| Compared rows | 100 |
| Failures | 0 |
| Best energy at or below RDKit | 72/100 |
| Within 1 kcal/mol | 95/100 |
| Within 2 kcal/mol | 98/100 |
| Median energy delta | -0.000000013 kcal/mol |
| p90 energy delta | +0.192 kcal/mol |
| Maximum positive delta | +2.814 kcal/mol |

The quality profile improves the once-only v3 within-1-kcal/mol count from
94/100 to 95/100 while keeping timing effectively tied. It still misses the
predeclared robust speed condition: every timing block must favor CheMatic and
the paired interval must be favorable. v4 is now inspected data and must not
be reused as a sealed final holdout.

## Remaining acceptance work

- diagnose quality-tail rows without tuning against v4;
- add minimizer diagnostics that can select extra work without using oracle
  outcomes;
- freeze a new candidate before any further once-only evaluation;
- run at least 21 alternating blocks, peak memory, and a second host;
- repeat from exact published artifacts before making a release claim.
