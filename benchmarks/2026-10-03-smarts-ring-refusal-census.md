# SMARTS opt-in ring-count refusals: six rows, 18 cells

This is a classification of the **source-only** opt-in profile, not a reduction
of the published v1.0.30 SMARTS residuals. The pinned inputs are the
[source profile](../validation/results/v1.0.31-source-smarts-310k-profiles.json),
RDKit 2026.03.6 oracle
`validation/results/v1.0.30-rdkit-smarts-all-cells.jsonl.gz`, and published
v1.0.30 Python rows
`validation/results/v1.0.30-published-python-chemistry-rows.jsonl.gz`.
The profile identifies exactly six refused corpus rows; each refuses `[R1]`,
`[R2]`, and `[R3]` with `RingModelAmbiguous`. The counts below are matched
atom counts, **not** proof that equal-size match sets are identical.

| Corpus row | RDKit `[R1]/[R2]/[R3]` | Published-default differing counts | Refused cells |
| ---: | ---: | --- | ---: |
| 9 | 14 / 18 / 14 | 32 / 10 / 4 | 3 |
| 23 | 8 / 6 / 2 | 28 / 16 / exact `[R3]` | 3 |
| 28 | 8 / 8 / 20 | 30 / 16 / 2 | 3 |
| 29 | 8 / 6 / 2 | 32 / 12 / 4 | 3 |
| 30 | 8 / 6 / 2 | 31 / 12 / 4 | 3 |
| 34 | 12 / 20 / 13 | 28 / 15 / 4 | 3 |

Thus 17 of the 18 refused cells replace a published-default **match-set**
error; row 23 `[R3]` replaces a previously correct match set. The later
[Fe ring-view source check](2026-10-03-source-smarts-organometallic-ring-view.md)
reports zero wrong-confident cells elsewhere, but these 18 refusals are not
exact matches. All six structures are
large, doubly charged fused polycycles with aromatic `[n+]` atoms. The
guard in `rdkit_parity_match.rs` triggers when at least two such cations and
one extra candidate ring occur; it rejects the entire `[RN]` query instead
of guessing among competing replacement-ring families.

The next implementation step is to compare the *actual ring sets* and their
atom/bond memberships, not merely their counts, for all six structures
against the pinned RDKit oracle. Only relax the guard for a chemically and
algorithmically justified subset, then rerun all 310,000 pinned cells and
the original 200 residual identities. Preserve typed refusal for unresolved
families. A published opt-in binding and artifact rerun remain separate gates.
