# Published v1.0.31 MMFF94 iteration-limit diagnostic

This changes **only the runtime configuration** of the already-published
v1.0.31 wheel (SHA-256
`b481316bb5a49cb572b035c805316e38a381f4f0e6fa141a05d1da3dc6f05c1f`).
It does not change the package default, the force-field implementation, or
the release's [200-iteration quality result](2026-10-03-a6-published-v1031-mmff94-quality.md).
Both runs used the same 265 Tier A+B molecules, seed 20260801, eight embed
attempts, 20-second row limit, strict stereo-safe MMFF94 arm, diagnostic-only
ring torsions, and the same RDKit 2026.03.6 reference rows. The runner's
metadata confirms that `force_field_max_iterations` was the only changed
configuration field: **200 → 400**.

| Quality dimension | 200 iterations | 400 iterations |
|---|---:|---:|
| Completed CheMatic rows | 265/265 | 265/265 |
| Force-field converged | 100 | 164 |
| Force-field not converged | 165 | 101 |
| Independent geometry sound | 265 | 265 |
| Independent stereo satisfied / violated / unevaluable | 170 / 0 / 0 | 170 / 0 / 0 |
| Independent gross clashes | 0 | 0 |
| Rows with a >15% bond-length deviation | 26 | 26 |
| Rows with a >50% bond-length deviation | 0 | 0 |

Of the 165 previously non-converged rows, 161 had hit the 200-iteration cap.
With 400 iterations, 64 gain convergence and none lose it. Of the remaining
101, 97 hit the new 400-iteration cap and four stopped earlier. Coordinates
changed on 161 rows; result statuses and parameter coverage did not change.
The independent scorer found no new geometry, stereo or gross-clash failures,
but this single cohort does **not** justify changing the shipped default.
The 101 non-converged rows still need residual-force/optimizer-failure
classification, and independent conformer-quality non-inferiority is open.

The external 400-iteration runner JSONL has SHA-256
`f5a12d84230ffe078dce14e9571d9d9f2b91b7db7bd1ada960212ceea0c67155`;
its metadata has SHA-256
`a91c1947522874a48f38ea7d72e2d6e306afaec8469557b60aba7095729f8239`;
the common-scored JSONL has SHA-256
`9ab3b1b8328a435cb6b480394994fbca20dde6e9d0b13cdaaa6745c4ff0ac07d`.
These large intermediates remain outside Git. To reproduce without editing
the runner file, import `scripts.public_package_3d_chematic`, set its
`FORCE_FIELD_MAX_ITERATIONS` constant to 400 in memory, and call its `main()`
with the same wheel, tiers, arm, and output arguments shown in the 200-iteration
record. The metadata must record 400 before scoring with
`pipeline_v2_vs_rdkit_common_scorer` against the pinned RDKit rows.
