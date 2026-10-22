# ETKDG quality profile: retry the lowest unfinished candidate

Date: 2026-10-22. This is a source comparison against RDKit 2026.09.1 on an already used 12-molecule Tier B development slice. It is not a fresh holdout or a published-package result.

The quality profile embeds 40 distinct-seed candidates, optimizes each with MMFF94 for at most 200 iterations, and returns at most 20. It now retries the lowest finite-energy **unfinished** returned candidate for at most 400 more iterations, even when the overall lowest-energy candidate has converged. The earlier policy retried only when the overall minimum was unfinished. A retry replaces its candidate only if it converges, does not raise energy, has finite coordinates, and preserves stereochemistry. The compatibility embedding API is unchanged.

| Metric on 12 molecules | Previous profile | Candidate |
| --- | ---: | ---: |
| Returned conformers | 240 | 240 |
| Converged returned conformers | 190 | 191 |
| Molecules with more converged conformers | — | 1 |
| Molecules with fewer converged conformers | — | 0 |
| Raw minimum-energy regressions | — | 0 |
| Returned candidate CIP / double-bond stereo violations | 0 / 0 | 0 / 0 |

The additional converged candidate is on `chembl_tier_b_0191`; its energy changes from 2.726637671956862 to 2.726637541568934 kcal/mol. Across the slice, CheMatic's raw minimum is at or below RDKit's on 11/12 molecules. Two molecules have a strictly lower converged minimum than RDKit's 40-candidate result, under the unequal retry budget. The remaining gap on `chembl_tier_b_0146` is 0.308065 kcal/mol. RDKit MMFF94 rescoring of all 240 returned CheMatic coordinate sets differs by at most 9.95e-14 kcal/mol after atom-indexed graph and bond-order checks.

The [machine record](../validation/results/etkdg-quality-lowest-unfinished-2026-10-22.json) pins source, probe binary, manifest, script, prior report and RDKit version. Reproduce it with:

```sh
cargo build --locked -p chematic-3d --example etkdg_quality_seed_budget_probe
/path/to/rdkit-2026.09.1/bin/python scripts/benchmarks/etkdg_quality_seed_budget_probe.py \
  --binary target/debug/examples/etkdg_quality_seed_budget_probe \
  --output validation/results/etkdg-quality-lowest-unfinished-2026-10-22.json \
  --prior-report validation/results/etkdg-quality-unfinished-minimum-retry-2026-10-22.json \
  --profile etkdg_quality_lowest_unfinished_development_v5
```

The retry adds optimization work. This comparison does not establish an equal-compute win or measure elapsed time, memory, conformer diversity, experimental geometry, or populations. The Python binding and published wheels were not rebuilt for this packet.

Local verification: `cargo test --locked -p chematic-3d`, `cargo clippy --locked -p chematic-3d --all-targets -- -D warnings`, and `cargo fmt --all -- --check` passed.
