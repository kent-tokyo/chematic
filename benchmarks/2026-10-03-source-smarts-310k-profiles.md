# Source SMARTS profile on the pinned 310,000-cell corpus

This is a **source-only diagnostic**, not a rerun of a published wheel, crate,
or npm package. The source branch starts at `2434d8b7` (v1.0.31 line); the
only Rust change in this packet is to the diagnostic example, which now mirrors
the Python binding's aromaticity preprocessing. It does not change production
SMARTS behavior. The pinned comparator is RDKit 2026.03.6.

The v1.0.30 published Python archive has **200 differing match sets** across
64 of 10,000 molecules and 31 queries (310,000 cells). Only **43** of those
cells differ in hit/no-hit Boolean; 157 retain the same Boolean. The current
source default reproduces the exact 200 and 43 cell identities. That gate is
necessary before attributing any change to the opt-in profile.

| Opt-in profile outcome | Match sets | Hit/no-hit |
|---|---:|---:|
| Published/default residual cells corrected | 154 | 37 |
| Wrong confident result remains | 29 | 6 |
| Typed refusal among formerly wrong cells | 17 | 0 |
| Typed refusal among formerly correct cells | 1 | 18 |
| Newly wrong confident cells | 0 | 0 |

The opt-in profile refuses 18 ring-model cells in six charged-polycycle rows
(`RingModelAmbiguous`, three ring-count queries per row). These are the **same
18 cells**, not 18 additional to each column. In match-set terms, 154 + 29 +
17 = 200 original residuals; the extra refusal removes one formerly correct
match set. In Boolean terms, 37 + 6 = 43 original residuals, while all 18
refusals remove formerly correct Boolean answers. Thus this is not a simple
“200 down to 29” compatibility claim: applicability narrows.

The two opt-in configurations measured here (`parity` and
`shared_symmetrized`) produce identical counts on this corpus. Remaining
wrong-confident match sets are concentrated in `[k6]` (10), `[R2]` (6), and
`[R3]` (5); the remaining Boolean errors are `[R2]` (3), `[R1]` (1), `[R3]`
(1), and `[x2]` (1). This does not establish that RDKit's ring model is the
chemical specification. Native SSSR semantics remain a separate contract.

## Reproduction and scope

The checked [machine summary](../validation/results/v1.0.31-source-smarts-310k-profiles.json)
pins the source dump, archived oracle, published rows, query list, and
adjudicated baseline by SHA-256 and includes row IDs for the remaining
Boolean errors and typed refusals. The large, intermediate source dump was
kept outside Git. Recreate it with:

```sh
CARGO_TARGET_DIR=/tmp/chematic-smarts-target CARGO_INCREMENTAL=0 CARGO_BUILD_JOBS=2 \
  cargo run -p chematic-smarts --example rdkit_parity_dump -- \
  validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  > /tmp/chematic-smarts-10k-perceived.jsonl
python3 scripts/check_smarts_310k_profiles.py \
  --dump /tmp/chematic-smarts-10k-perceived.jsonl \
  --output /tmp/chematic-smarts-10k-perceived-full-profile.json
```

The checker requires the exact archived source-default 200 match-set and 43
Boolean cell IDs, corpus/query hashes, row order, complete footer, and no
parse refusals. Missing opt-in results without an explicit error or budget
reason fail closed. A 21-structure hand corpus is in the dump, but the table
above is only the pinned 10,000-molecule corpus. No Python/npm opt-in binding
or published-artifact equivalence is claimed. The next implementation target
is the six Boolean wrong-confident cells, then the remaining `[k6]`/ring
match-set differences; retain the refusal boundary until independently
adjudicated.
