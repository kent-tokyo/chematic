# ETKDG quality profile: distinct candidate seeds

Date: 2026-10-22. Source candidate at `10346020`; local wheel only. This is
development evidence for one conformer-search metric, not a broad RDKit quality
claim.

## Change

The quality profile formerly embedded the same `num_confs` seeds twice when
there were no fifth- or sixth-order experimental torsion coefficients. Its
compatibility and corrected gradient paths are identical in that case. It now
uses seeds `1..N` and `N+1..2N`, optimizes both ensembles with MMFF94, and
returns the selected batch of up to `N` conformers with the lower finite best
energy.
For a molecule with fifth- or sixth-order coefficients, it still compares the
two gradient paths on the same seeds. The public compatibility and corrected
APIs retain their existing behavior. For the identical-gradient case, a zero
seed or seed-range overflow falls back to the initial ensemble.

This profile performs about twice the embedding and optimization work of the
20-conformer compatibility call. It retains only 20 conformers even when it
examines 40 candidates; downstream conformer diversity was not measured.

## Development comparison

The fixed Tier B manifest was used in earlier development, so this is **not a
held-out set**. Before scoring, the script selected single-fragment,
non-macrocyclic molecules with at most 30 heavy atoms, 3–10 RDKit rotatable
bonds, MMFF parameters, and no fifth- or sixth-order ETKDG torsion coefficient.
It sorted eligible rows by SHA-256 of SMILES and took the first 12. Both engines
used ETKDGv3, seed 42, and MMFF94 with up to 200 iterations. CheMatic examined
20 + 20 distinct seeds; RDKit 2026.09.1 examined 40. Results are minimum final
MMFF94 energy in kcal/mol, including candidates that did not converge in 200
iterations. The tolerance for numerical agreement is 1e-6 kcal/mol.

| Tier B row | CheMatic quality − former profile | CheMatic quality − RDKit 40 | CheMatic / RDKit converged candidates |
|---|---:|---:|---:|
| 0192 | −0.988815 | < 1e-6 | 20 / 40 |
| 0177 | 0 | < 1e-6 | 13 / 30 |
| 0051 | 0 | +0.006473 | 0 / 1 |
| 0159 | 0 | < 1e-6 | 17 / 30 |
| 0174 | < 1e-6 | < 1e-6 | 20 / 40 |
| 0146 | −1.754894 | +0.308065 | 13 / 26 |
| 0126 | −0.804823 | < 1e-6 | 14 / 29 |
| 0000 | < 1e-6 | < 1e-6 | 20 / 39 |
| 0095 | 0 | +0.010105 | 11 / 27 |
| 0169 | 0 | < 1e-6 | 20 / 38 |
| 0191 | < 1e-6 | < 1e-6 | 1 / 7 |
| 0154 | −0.528759 | −0.527364 | 0 / 0 |

The new profile lowered the best energy relative to the former profile on
4/12 rows and raised it on none. Compared with RDKit's 40 candidates, it was
strictly lower on 1/12, within numerical tolerance on 8/12, and higher on
3/12. **The strict win is row 0154, where neither engine converged any
candidate**. Thus these results do not demonstrate a converged-energy win or
general conformer accuracy superiority. The largest remaining energy deficit
was +0.308065 kcal/mol. CheMatic's retained ensemble sometimes has fewer
converged candidates than RDKit's 40-candidate ensemble.

The script checked identical atom-indexed graphs and bond orders, then
independently re-scored both CheMatic best geometries with RDKit MMFF94. The
maximum same-coordinate energy discrepancy was 9.60e-14 kcal/mol. The best
geometries had zero input-versus-3D CIP violations and zero input-versus-3D
double-bond stereo violations on this set. The full row-level data and pinned
source, manifest, and executable hashes are in
[`validation/results/etkdg-quality-seed-budget-2026-10-22.json`](../validation/results/etkdg-quality-seed-budget-2026-10-22.json).

The earlier six molecules with a nonzero sixth-order term were rerun after the
change: two meaningful energy wins and four numerical ties versus RDKit 40,
matching the previous development comparison. The corrected-gradient branch
therefore retained its measured behavior. RDKit has already merged its own
sixth-order gradient correction ([issue #8715](https://github.com/rdkit/rdkit/issues/8715),
[PR #9699](https://github.com/rdkit/rdkit/pull/9699)), so that separate
comparison is specific to the installed 2026.09.1 release.

## Reproduce and limits

Check out the pinned `10346020` source candidate before running this v1
command. Later source uses the same probe entry point for the pooled v2 report.

```sh
cargo build --release -p chematic-3d --example etkdg_quality_seed_budget_probe --locked
PYTHONPATH=/path/to/rdkit-2026.09.1 python3 scripts/benchmarks/etkdg_quality_seed_budget_probe.py \
  --binary target/release/examples/etkdg_quality_seed_budget_probe \
  --output validation/results/etkdg-quality-seed-budget-2026-10-22.json
```

Rust library tests: 708 passed, 4 ignored; Clippy and formatting passed.
Python embedding tests: 13 passed against a locally built CPython 3.13 arm64
wheel (SHA-256 `215140e73794d9b7427c0c980c62bdffe4badc23c1d65bb2dec66c0f535265e7`).
The wheel is not a published artifact.

This comparison does not measure experimental geometry, conformer coverage or
diversity, elapsed time, memory, a fresh held-out corpus, or convergence beyond
200 iterations. Those measurements are required before claiming overall 3D,
force-field, or conformer superiority over RDKit.
