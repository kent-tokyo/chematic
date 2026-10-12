# MMFF final-energy cache

Date: 2026-10-13

## Change

RDKit-compatible BFGS already evaluates the energy at every accepted point.
MMFF94 and UFF previously evaluated the complete force field once more after
minimization only to report that same value. The candidate reuses the accepted
energy. RDKit's stale-distance termination remains on the old path because its
reported value deliberately combines restored coordinates with the last trial
point's distance cache.

The public `BfgsOutcome` shape and minimizer status are unchanged.

## Verification

- Oracle: RDKit 2026.09.1
- Host: macOS arm64, Python 3.13.6
- Work: first 100 A6 development rows, 20 ETKDGv3 candidates, 200 MMFF94
  iterations, three alternating repeats
- Before wheel SHA-256:
  `ade612cda0ef976b95b8847caa1488fd83c9338fa517a85fa149718c0de3e04d`
- Candidate wheel SHA-256:
  `bcb57336d4fcd09af4cfbfc7ffc746bd0bbe67a37f54bc676d632a2b885cc71f`

All 100 per-molecule rows are identical before and after, including generated
and converged counts, best energy, and RDKit delta. The wider 265-row
development result is also identical in every row.

| Lane | Before | Candidate |
|---|---:|---:|
| CheMatic median | 10.834 s | 10.821 s |
| RDKit median in the same run | 11.522 s | 11.557 s |
| RDKit / CheMatic ratio | 1.063x | 1.068x |

The removed force-field evaluation is deterministic work reduction, but this
three-repeat change is small relative to host variation. It is not a universal
speed claim and does not address the remaining conformer-quality tail.
