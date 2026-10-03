# SMARTS `[RN]`: six refused polycycles, ring-family diagnosis

This is a **source-only diagnostic** on main `5ac44077` with pinned RDKit
2026.03.6, not a matcher fix or a published-package comparison. It follows
the [18-cell refusal census](2026-10-03-smarts-ring-refusal-census.md).
The six SMILES come from the 10,000-row corpus (SHA-256
`f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`).
The checker validates atom count, element, charge and input row alignment
before comparing *ring atom sets, ring bond sets* and per-atom ring counts.

| Row | RDKit rings | CheMatic basis / shared rings | Shared sets missing / extra vs original RDKit | Atoms with wrong shared ring count | Distinct RDKit ring families in 1,000 renumberings |
| ---: | ---: | ---: | ---: | ---: | ---: |
| 9 | 8 | 7 / 8 | 1 / 1 | 4 | 3 |
| 23 | 10 | 7 / 10 | 0 / 0 | 0 | 3 |
| 28 | 9 | 7 / 9 | 0 / 0 | 0 | 1 |
| 29 | 10 | 7 / 10 | 1 / 1 | 4 | 3 |
| 30 | 10 | 7 / 10 | 1 / 1 | 4 | 3 |
| 34 | 9 | 8 / 9 | 0 / 0 | 0 | 2 |

The opt-in SMARTS-specific approximate model agrees with the original RDKit
per-atom ring counts on only row 34. The alternate shared symmetrized model
agrees exactly on rows 23, 28 and 34, but not 9, 29 or 30. On those three,
one large equal-size ring is selected along a different path; the affected
atom indices are 14,15,17,22 (row 9), 27,28,30,35 (row 29), and
26,27,29,34 (row 30). The 1,000 deterministic `RenumberAtoms` checks
preserved canonical isomeric SMILES and observed multiple RDKit ring families
on five of the six graphs. This is **observed order dependence**, not an
exhaustive enumeration of all ring families; the single observed family on
row 28 is not proof of uniqueness. Bond-set comparisons give the same 3/3
exact/different split as atom-set comparisons, including under permutation.

For rows 9, 29 and 30, the current CheMatic shared ring family did **not**
occur in the 1,000 observed RDKit families. The supported-domain claim thus
cannot be widened by merely swapping the two CheMatic selectors or removing
the `RingModelAmbiguous` guard. It would create wrong-confident atom sets on
three rows. Preserve typed refusal while deriving a principled ring-family
selection or an explicit, order-independent semantic contract. Even the
three original-input exact cases need a runtime-applicable safety criterion;
do not special-case the six corpus IDs.

Reproduce from this source revision with RDKit 2026.03.6 installed:

```sh
set -o pipefail
cargo run -p chematic-smarts --example ring_refusal_probe --locked -- \
  validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  | python scripts/check_smarts_ring_refusal_probe.py --permutations 1000
```

The Rust probe is diagnostic-only and does not expose atom-map or production
SMARTS API changes. Its checked source-profile input is
`validation/results/v1.0.31-source-smarts-310k-profiles.json` (SHA-256
`8797e777c6b203d5994fedae80a254e10284f5dd721b110e43e48bf7e29624b6`).
The original 200 match-set and 43 Boolean published residuals, the 18 opt-in
refusals, and the full 310,000-cell publication gate remain open.
