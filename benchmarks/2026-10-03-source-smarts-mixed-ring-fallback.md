# Source-only SMARTS mixed-size bridge-ring follow-up

This follows the [opt-in `[kN]` packet](2026-10-03-source-smarts-ring-size-follow-up.md).
It changes only the opt-in RDKit-style `[RN]` ring-count model; the native
SSSR matcher and its v1.0.30 published 200 match-set / 43 Boolean residuals
are unchanged. The pinned comparator is RDKit 2026.03.6 on the same 10,000
molecules × 31 queries. This is **source-only**, not a Python/npm release
artifact comparison.

In corpus row 3498, a four-membered basis ring shares a three-atom path
with a six-membered basis ring. RDKit includes a second, symmetry-equivalent
six-membered ring, giving five rings in its ring information. CheMatic's
previous opt-in fallback searched overlapping basis rings only when the two
basis rings had equal size; it therefore missed this valid same-size
*candidate* ring. The fallback now considers any overlapping basis-ring
pair, while the candidate acceptance rule still requires a same-size
replacement, shared bond and preservation of unique basis bonds, and retains
the existing candidate budget. A row-specific test pins all 27 RDKit ring
membership counts.

| Opt-in outcome | Before this follow-up | After |
|---|---:|---:|
| Published/default match-set residuals corrected | 164 | 170 |
| Remaining wrong-confident match sets | 19 | 13 |
| Published/default Boolean residuals corrected | 37 | 38 |
| Remaining wrong-confident Booleans | 6 | 5 |
| Typed refusal cells | 18 | 18 |
| Newly wrong-confident cells | 0 | 0 |

The six newly corrected match sets are `[R1]`, `[R2]`, `[R3]` on rows 3498
and 3995. Only row 3498 `[R3]` changes its hit/no-hit Boolean. The remaining
Boolean errors are `[R1]` (one), `[R2]` (three), and organometallic `[x2]`
(one). One `[k6]` matched-atom-set difference remains on row 3498 because
the separate bounded symmetrized ring-size selector still lacks this
alternative ring. The 18 existing typed refusals are unchanged: 17 formerly
wrong match sets and one formerly correct match set, all with formerly
correct Booleans. Do not report the 170 corrections as unconditional parity.

Reproduce using `rdkit_parity_dump` followed by
`scripts/check_smarts_310k_profiles.py` as in the parent packet. The
external source dump SHA-256 was
`c6acbe45d8f7fdec4ec4bb662a4123794e62e6d5f2d5368454f2318d4779e1db`;
the checker summary SHA-256 was
`9eba96a98f90f6ddcba1acb3fe669ba078ddfe153218ae504549712753e2dae4`.
The checker verified all 310,000 cells and the exact source-default baseline
identities. The full source dump is not checked in. No new typed refusal,
wrong-confident result or publication-level claim is implied.
