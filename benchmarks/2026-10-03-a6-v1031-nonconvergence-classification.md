# A6: classify the 101 MMFF94 non-converged rows

This pairs the **published v1.0.31** wheel's already-recorded 265-row runs
at 200 and 400 maximum iterations. It is an optimizer diagnostic, not a
source-wheel result, a new conformer-quality result, or permission to raise
the default iteration limit. The raw 200-row and 400-row JSONL SHA-256 values
are `c2525cf7a12f7fd77b1881ee868b49785b7b62b44e006eb708d2c630b6e0d850`
and `f5a12d84230ffe078dce14e9571d9d9f2b91b7db7bd1ada960212ceea0c67155`.
The pinned [row-level classification](../validation/results/a6-mmff94-v1031-nonconvergence-paired-200-400.json)
has SHA-256 `7c6aa135a55cd8b63b4a04c8187e95a1f4fdb8946442ff1f0140e24d418de0f6`.

| Outcome | Rows |
|---|---:|
| Converged at 200 iterations | 100 / 265 |
| Converged at 400 iterations | 164 / 265 |
| Gained convergence / lost convergence | 64 / 0 |
| Still not converged at 400 | 101 |
| Of those: hit 400-iteration cap | 97 |
| Of those: stopped before cap | 4 |

The four early stops are rows **170, 191, 227 and 234**, at 136, 72, 74
and 74 iterations. Each stops at the *same* iteration in both runs, with
recorded maximum residual forces of approximately 23.65, 8.05, 23.25 and
14.69, respectively. The raw result does not contain an optimizer termination
reason, so calling these “line-search failures” or “converged” would be an
unsupported inference. Add a typed termination reason before diagnosing them.

Among the 97 cap-hit rows, residual-force bands are 11 below 0.01, 20 from
0.01 to 0.1, 50 from 0.1 to 1, and 16 from 1 to 10. Six non-converged rows
(74, 88, 93, 94, 95, 99) also record missing torsion parameters; this overlaps
with the cap-hit group and is not a mutually exclusive cause of non-convergence.
All 101 have finite energy that did not increase from before to after the
optimizer, but this does not prove a stationary point or correct conformer.

Reproduce from the pinned raw files described in the
[400-iteration record](2026-10-03-a6-v1031-iteration400-diagnostic.md):

```sh
python3 scripts/classify_a6_mmff94_nonconvergence.py \
  --baseline /path/to/chematic-3d.jsonl \
  --extended /path/to/chematic-3d-iter400.jsonl \
  --output validation/results/a6-mmff94-v1031-nonconvergence-paired-200-400.json \
  --gate-v1031
```

Next: instrument optimizer termination and residual-gradient checks for the
four early stops; stratify the 97 cap-hit rows by residual force and missing
parameters; rerun independent stereo/clash/geometry and conformer-quality
gates before changing defaults. This is separate from the 451-heavy-atom
source-wheel typing result and the published 2,908-heavy-atom baseline.
