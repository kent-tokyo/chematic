# ETKDG quality profile: retry an unfinished minimum

Date: 2026-10-22. This is a source comparison on an already used 12-molecule
Tier B development slice. It is not a held-out or released-package result.

## Change

The quality profile embeds 40 distinct-seed candidates in two batches,
optimizes each with MMFF94 for at most 200 iterations, and returns at most 20.
It now continues the lowest finite-energy candidate for at most 400 additional
iterations whenever that candidate has not converged, even if another
candidate has converged. Previously the retry ran only when **none** of the
returned candidates had converged. It replaces the candidate only if the retry
converges, does not raise its energy, has finite coordinates, and preserves
the input stereochemistry. The compatibility API is unchanged.

## Development comparison

RDKit 2026.09.1 used ETKDGv3, seed 42, 40 candidates and up to 200 MMFF94
iterations per candidate. CheMatic used the same initial candidate and
iteration limits plus the conditional retry. This is a comparison of a
quality-oriented policy with **unequal optimization budgets**, not an
equal-compute ranking.

| Metric across 12 molecules | Previous profile | New profile |
|---|---:|---:|
| Returned candidates | 240 | 240 |
| Converged returned candidates | 188 | 190 |
| Molecules with a lower raw minimum | — | 2 |
| Raw minimum-energy regressions | — | 0 |
| Molecules with more converged returned candidates | — | 2 |
| Returned candidate stereo violations | 0 | 0 |

On `chembl_tier_b_0051`, the best converged energy changes from
−43.184048 to −46.609310 kcal/mol; RDKit's 40-candidate best converged
energy is −43.184048. On `chembl_tier_b_0095`, it changes from 102.590755
to 102.077675 kcal/mol; RDKit's is 102.590756. These are two strict
converged-energy wins within this development slice, purchased with an extra
retry. The source profile has a converged candidate on `chembl_tier_b_0154`,
where RDKit's 40 candidates have none, as in the preceding report.

The raw minimum is at or below RDKit's on 11/12 molecules. The remaining
converged-energy gap on `chembl_tier_b_0146` is 0.308065 kcal/mol. All 240
returned conformers pass the input CIP and double-bond stereochemistry checks.
Independent RDKit MMFF94 scoring of the same atom-indexed coordinates differs
by at most 9.95e-14 kcal/mol after graph and bond-order checks. The
[machine record](../validation/results/etkdg-quality-unfinished-minimum-retry-2026-10-22.json)
contains row-level results and hashes for the source, binary, script,
manifest, and previous report.

## Reproduce and limits

```sh
cargo build --locked -p chematic-3d --example etkdg_quality_seed_budget_probe
/path/to/rdkit-2026.09.1/bin/python scripts/benchmarks/etkdg_quality_seed_budget_probe.py \
  --binary target/debug/examples/etkdg_quality_seed_budget_probe \
  --output validation/results/etkdg-quality-unfinished-minimum-retry-2026-10-22.json \
  --prior-report validation/results/etkdg-quality-convergence-retry-2026-10-22.json \
  --profile etkdg_quality_unfinished_minimum_retry_development_v4
```

This does not measure elapsed time, memory, conformer diversity, experimental
geometry, populations, cross-platform numerical agreement, a fresh held-out
corpus, or published artifacts. The Python binding was not rebuilt in this
measurement. The minimum force-field energy is one quality metric and cannot
establish overall 3D or conformer superiority.
