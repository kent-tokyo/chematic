# A6: independent MMFF94 quality dimensions

This is an **audit of existing evidence**, not a v1.0.31 rerun or a new
force-field result. `scripts/a6_quality_matrix.py` validates the two reversed
v1.0.26 published-wheel runs with their existing hash/row/scorer checker,
and the separate v1.0.19-source same-coordinate packet with its own checker.
The 10k atom-type census comes from a v1.0.25 wheel. These artifacts have
different code versions and corpora; the matrix never merges their success
counts into a single parity percentage.

The subsequent [published v1.0.31 rerun](2026-10-03-a6-published-v1031-mmff94-quality.md)
confirms 265/265 geometry/stereo/clash and the same 100/265 convergence count.
Its separate same-coordinate lane closes the two historical >5 kcal/mol
total-energy residuals on 262 comparable rows; per-term and coverage exits
remain open. This historical matrix remains unchanged so its provenance is
unambiguous.

| Dimension | Measured evidence | A6 status |
|---|---|---|
| Atom types and parameters | v1.0.26 strict bond-angle arm has no missing bond/angle/OOP/stretch-bend terms, but 32 missing torsion terms across six rows; separate v1.0.25 10k census has 2,908 differing heavy-atom types | Open; strict arm permits those torsion gaps |
| Convergence | v1.0.26 published: 100/265 `converged`, 165/265 not converged; identical row IDs in both run orders | Open; a usable geometry is not proof of optimization convergence |
| Stereo | v1.0.26: 265/265 independently stereo-clean on both runs | Passed on measured artifact only |
| Clash/geometry | v1.0.26: 265/265 independently sound and gross-clash-free on both runs | Passed on measured artifact only |
| Same-coordinate energy | v1.0.19 source: 262 comparable, 260 within 5 kcal/mol; rows 166 and 231 exceed 5 kcal/mol | Open; different source/artifact from the v1.0.26 wheel |
| Independent conformer quality | No non-inferiority gate; between-engine RMSD is not one | Not measured |

The [machine-readable matrix](../validation/results/a6-separated-quality-matrix-v1.0.31-audit.json)
pins its source files by SHA-256 and records all 165 non-converged row IDs.
Its `ready_for_target_version` is **false** for v1.0.31. The `--require-ready`
mode exits 2, which prevents a historical 265/265 geometry result from being
mistaken for a current-release A6 acceptance. Speed is explicitly excluded.

```sh
python3 scripts/a6_quality_matrix.py --target-version 1.0.31 \
  --output /tmp/chematic-a6-quality-matrix.json
python3 scripts/a6_quality_matrix.py --target-version 1.0.31 --require-ready \
  --output /tmp/chematic-a6-quality-matrix.json
```

Next: rerun each dimension on one pinned current artifact and matched cohort;
resolve the heavy-atom typing classes and same-coordinate energy residuals;
make non-convergence a typed, visible outcome; and add an independent
conformer-quality criterion. Until then 3D/MMFF94 remains Experimental.
