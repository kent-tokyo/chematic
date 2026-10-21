# ETKDG quality profile: bounded convergence retry

Date: 2026-10-22. This is a local macOS arm64 source candidate on an already used
12-molecule development slice, not a held-out or released-package result.

## Change

The quality profile embeds 40 candidate conformers in two batches, optimizes
each with MMFF94 for up to 200 iterations, and returns at most 20. When none
of the returned candidates has converged, it now retries the lowest finite
energy candidate for up to 400 additional iterations. The retry is retained
only if it converges, does not raise the energy, has finite coordinates, and
passes the explicit stereochemistry check. The compatibility API retains its
existing budget and behavior. This retry adds optimization work in the
affected case; the comparison below is not an equal-compute speed or quality
ranking.

## Fixed development comparison

RDKit 2026.09.1 used ETKDGv3, seed 42, 40 candidates and up to 200 MMFF94
iterations per candidate. CheMatic used two 20-candidate batches with the
same initial limit, plus the conditional retry described above.

| Metric | Previous pooled profile | Retry profile |
|---|---:|---:|
| Returned candidates across 12 molecules | 240 | 240 |
| Converged returned candidates | 187 | 188 |
| Molecules with lower raw minimum energy | — | 1 |
| Molecules with higher raw minimum energy | — | 0 |
| Molecules with a new converged candidate | — | 1 |
| Returned candidate stereo violations | 0 | 0 |

For `chembl_tier_b_0154`, the previous pooled result had no converged
candidate and a 200.903754 kcal/mol raw minimum. The retry yielded one
converged candidate at 198.265822 kcal/mol. RDKit's 40 candidates had no
converged result and a 201.431118 kcal/mol raw minimum under its 200-iteration
limit. These are different optimization budgets. The result demonstrates a
useful bounded rescue for this molecule, not superior general 3D accuracy.

The quality profile's raw minimum was at or below RDKit's on 9/12 molecules.
The remaining converged-energy gap for `chembl_tier_b_0146` is 0.308065
kcal/mol. All 240 returned conformers passed the input CIP and double-bond
stereo checks. Independent RDKit MMFF94 scoring of the same atom-indexed
coordinates differed by at most 9.95e-14 kcal/mol after graph and bond-order
checks. The [machine record](../validation/results/etkdg-quality-convergence-retry-2026-10-22.json)
pins both source files, the binary, script, manifest, and prior report hashes.

## Reproduce and limits

```sh
cargo build --locked -p chematic-3d --example etkdg_quality_seed_budget_probe
/path/to/rdkit-2026.09.1/bin/python scripts/benchmarks/etkdg_quality_seed_budget_probe.py \
  --binary target/debug/examples/etkdg_quality_seed_budget_probe \
  --output validation/results/etkdg-quality-convergence-retry-2026-10-22.json \
  --prior-report validation/results/etkdg-quality-pooled-2026-10-22.json \
  --profile etkdg_quality_pooled_development_v3
```

Rust 3D library tests: 710 passed, 4 ignored. Clippy and formatting passed.
This report does not measure
elapsed time, memory, conformer diversity, experimental geometry or
populations, cross-platform numeric agreement, a fresh held-out corpus, or
public artifacts. The Python binding was not rebuilt in this measurement.
