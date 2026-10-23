# 3D alignment fix and ETKDG quality follow-up

Date: 2026-10-24. The source candidate starts from main
`38b4c92c430d9b3bfe70e47e97293b9017ced990`. The measured ETKDG code is
unchanged from that base; this candidate fixes Kabsch transform application
and O3A partial-correspondence refinement. The local comparator is RDKit
2026.09.1 on Darwin arm64 with Python 3.13.6. The binary and input hashes are
in the [machine record](../validation/results/etkdg-quality-tier-b-followup-2026-10-24.json).

## Alignment correction

`align_coords` now returns an ordinary rigid transform, with translation
computed from the centroids used for the fit. `apply_alignment` applies this
same transform to every supplied atom. Previously it recomputed a centroid
from the supplied atoms, which shifted a full molecule incorrectly after
fitting only a subset of its atoms. Kabsch's rank-deficient planar and linear
fits now also return a proper rotation instead of a projection matrix.

A fixed O3A regression has two corresponding atoms and one unmatched atom far
from them. Its iterative overlap score reaches 2.0 after the subset fit;
the old full-molecule centroid shift prevented that exact match. Planar and
linear regression tests check the rotation determinant and transformation of
an atom excluded from the fit. This is a geometry-correctness fix, not a
measured win in experimental conformer accuracy.

## Additional RDKit comparison

The pre-existing Tier B manifest has 33 eligible single-fragment molecules
with at most 30 heavy atoms, 3–10 rotatable bonds, available MMFF94 parameters,
and no fifth- or sixth-order experimental torsion coefficient. The first 12
were used in the earlier [third-seed batch](2026-10-22-etkdg-third-seed-batch.md).
This follow-up takes the next 21 by SHA-256 of SMILES. The manifest was already
used in development, and row `chembl_tier_b_0022` was inspected before this
record was finalized. These 21 are **development examples**, not a holdout.

CheMatic's quality profile embeds up to 60 candidates from three 20-seed
batches, optimizes each with MMFF94 for up to 200 steps, may make two further
400-step retries, and returns at most 20. RDKit embeds 60 ETKDGv3 candidates
with seed 42 and optimizes each with MMFF94 for up to 200 steps. The embedded
candidate count matches; optimization work and returned ensemble sizes do not.
The metric is the lowest final MMFF94 energy, including unfinished candidates,
with a 1e-6 kcal/mol equality tolerance.

| Metric, 21 molecules | Result |
| --- | ---: |
| CheMatic minimum at or below RDKit within tolerance | 20/21 |
| CheMatic strictly lower by more than tolerance | 8/21 |
| CheMatic higher by more than tolerance | 1/21 |
| Largest deficit: `chembl_tier_b_0022` | +0.433597 kcal/mol |
| Largest same-coordinate energy difference, independent RDKit re-score | 1.14e-13 kcal/mol |
| Returned candidate CIP / double-bond stereo violations | 0 / 0 |

The scorer verifies identical atom-indexed graphs and bond orders before
re-scoring all 420 returned CheMatic conformers with RDKit MMFF94. It checks
declared tetrahedral CIP and double-bond stereo on each returned conformer.
The remaining energy deficit is a conformer-generation problem on this sample;
same-coordinate MMFF94 energy agreement alone cannot close it. A fourth
20-seed batch was tried locally on row `0022`: it raised the number of
converged returned conformers from 16 to 20 but did not lower the minimum,
so it was not retained in the quality profile.

Reproduce the machine record with:

```sh
cargo build --locked -p chematic-3d --example etkdg_quality_seed_budget_probe
PYTHONDONTWRITEBYTECODE=1 PYTHONPATH=/path/to/rdkit-2026.09.1 \
  python3 scripts/benchmarks/etkdg_quality_followup_probe.py \
  --binary target/debug/examples/etkdg_quality_seed_budget_probe \
  --output validation/results/etkdg-quality-tier-b-followup-2026-10-24.json
```

The Rust `chematic-3d` test suites, Clippy with warnings denied, formatting,
and refactor-budget check passed locally. This comparison does not measure
experimental geometry, conformer diversity or population coverage, memory,
elapsed time, equal compute, or a published package. Broad RDKit 3D, MMFF,
or conformer superiority is not established.
