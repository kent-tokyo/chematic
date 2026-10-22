# ETKDG quality profile: third seed batch

Date: 2026-10-22. This source comparison uses RDKit 2026.09.1 and the same
already used 12-molecule Tier B development slice as the previous ETKDG
quality reports. It is not a fresh holdout or a published-package result.

For molecules without fifth- or sixth-order torsion coefficients, the quality
profile now searches three distinct batches of `num_confs` seeds. The existing
first two batches and their bounded MMFF94 retry are preserved. The third batch
adds new seeds; the lowest unfinished returned candidate may receive one more
bounded MMFF94 retry. Each retry is retained only if it converges, does not
raise energy, has finite coordinates, and preserves stereochemistry. The
profile still returns at most `num_confs` conformers. Molecules whose gradient
modes differ keep the existing two-batch path.

The development comparison uses `num_confs=20`, hence 60 CheMatic candidates.
RDKit embeds 60 ETKDGv3 candidates with seed 42. Both optimize each candidate
with MMFF94 for at most 200 iterations. CheMatic can spend up to two additional
400-iteration retries. This is an equal embedding-candidate count, **not** an
equal compute budget; CheMatic returns at most 20 versus RDKit's 60.

| Metric on 12 molecules | Result |
| --- | ---: |
| CheMatic minima improved versus previous 40-candidate quality profile | 5/12 |
| CheMatic minima regressed versus previous profile | 0/12 |
| CheMatic minima at or below RDKit 60 within 1e-6 kcal/mol | 12/12 |
| Maximum CheMatic minus RDKit 60 minimum | +0.000000787 kcal/mol |
| Returned candidate CIP / double-bond stereo violations | 0 / 0 |
| Maximum independent RDKit MMFF94 re-score difference on CheMatic coordinates | 9.95e-14 kcal/mol |

The former 0.308065 kcal/mol deficit on `chembl_tier_b_0146` versus RDKit 40
closes: CheMatic reaches 93.677174 kcal/mol, within 1e-6 of RDKit 60's
93.677173 kcal/mol. `chembl_tier_b_0051` reaches -46.669373 versus RDKit 60's
-46.666913 kcal/mol. On `chembl_tier_b_0154`, CheMatic has a converged
candidate at 198.265822 kcal/mol while RDKit 60 has no converged candidate at
the 200-iteration limit. These are bounded development comparisons, not
experimental geometry or population measurements.

The [machine record](../validation/results/etkdg-quality-third-seed-batch-2026-10-22.json)
contains every row and hashes of the source, probe binary, script, manifest,
prior report, and RDKit version. Reproduce it with:

```sh
cargo build --locked -p chematic-3d --example etkdg_quality_seed_budget_probe
PYTHONPATH=/path/to/rdkit-2026.09.1 python3 scripts/benchmarks/etkdg_quality_seed_budget_probe.py \
  --binary target/debug/examples/etkdg_quality_seed_budget_probe \
  --output validation/results/etkdg-quality-third-seed-batch-2026-10-22.json \
  --prior-report validation/results/etkdg-quality-lowest-unfinished-2026-10-22.json \
  --profile etkdg_quality_third_seed_batch_development_v6 \
  --rdkit-confs 60
```

Elapsed time, memory, conformer diversity, experimental geometry, and a fresh
holdout have not been measured. The additional search and retries add work, so
these results do not establish overall 3D, force-field, or conformer
superiority over RDKit.
