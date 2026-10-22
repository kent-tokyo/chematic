# ETKDG quality profile: pooled returned conformers

Date: 2026-10-22. This is a source candidate comparison on an already used
development slice, not a held-out test or a released-package result.

## Change

The quality profile embeds and MMFF94-optimizes two batches of up to 20
conformers. Previously it returned all 20 from whichever batch had the lower
finite minimum energy. The second batch was discarded even when it contained
more converged conformers. The profile now retains the global lowest finite
energy candidate, then fills the remaining 19 places with converged candidates
from both batches, in energy order. Non-converged candidates fill any remaining
places. Equal scores retain deterministic batch and seed order. It still
returns at most 20 conformers and uses the same 40-candidate search budget.

This favors a usable returned ensemble while preserving the previous raw
minimum-energy guarantee. A low-energy candidate that did not converge can
still occupy the first place; its energy is not evidence of a converged
minimum. The public compatibility embedding API is unchanged.

## Development comparison

The paired script selects the same 12 single-fragment, non-macrocyclic Tier B
structures as the preceding [seed-budget report](2026-10-22-etkdg-quality-seed-budget.md).
RDKit 2026.09.1 uses ETKDGv3, seed 42, 40 candidates, MMFF94, and 200
optimization iterations per candidate. CheMatic uses the same settings in
two 20-candidate batches. This corpus has already influenced development.

| Metric | Previous quality profile | Pooled profile |
|---|---:|---:|
| Returned conformers | 20 per structure | 20 per structure |
| Converged returned conformers, 12 structures total | 149 | 187 |
| Structures with more converged returned conformers | — | 7/12 |
| Structures with fewer converged returned conformers | — | 0/12 |
| Raw minimum-energy regressions above 1e-6 kcal/mol | — | 0/12 |
| All returned conformers: CIP / double-bond stereo violations | — | 0 / 0 |

Against RDKit's 40 candidates, CheMatic's raw minimum energy was lower or
within 1e-6 kcal/mol on 9/12 structures. For the 11 structures where both
engines had a converged candidate, converged minimum energies agreed within
1e-6 kcal/mol on 10/11. CheMatic was higher by 0.308065 kcal/mol on the
remaining structure and did not have a strict converged-energy win. The twelfth
structure had no converged candidate in either engine. The maximum
independent RDKit MMFF94 re-score discrepancy across all 240 returned
CheMatic coordinates was 9.95e-14 kcal/mol. Atom-indexed graphs and bond
orders were checked before scoring.

The row-level data, exact source and binary hashes, manifest hash, and prior
report hash are in
[`validation/results/etkdg-quality-pooled-2026-10-22.json`](../validation/results/etkdg-quality-pooled-2026-10-22.json).

## Reproduce and limits

```sh
cargo build -p chematic-3d --example etkdg_quality_seed_budget_probe --locked
/path/to/rdkit-2026.09.1/bin/python scripts/benchmarks/etkdg_quality_seed_budget_probe.py \
  --binary target/debug/examples/etkdg_quality_seed_budget_probe \
  --output validation/results/etkdg-quality-pooled-2026-10-22.json
```

Rust 3D library tests: 709 passed, 4 ignored. Clippy and formatting passed.
The comparison used a local source build; the Python binding was not rebuilt
in this measurement. The two engines returned different ensemble sizes
(20 versus 40), so the table does not establish a toolkit-wide quality win.
Experimental geometry, diversity, elapsed time, memory, a fresh held-out
corpus, and published artifacts remain unmeasured here.
