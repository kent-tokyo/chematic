# Symmetry-aware conformer clustering

Date: 2026-10-22. Source implementation and Python 3.13 development wheel on
macOS arm64. Comparator: RDKit 2026.09.1 from the local evaluation environment.
This is a focused correctness check, not a force-field or conformer-search
quality ranking.

## Change

`ConformerEnsemble::cluster_conformers_by_rms_symmetric` applies the existing
automorphism-aware aligned RMSD to greedy leader clustering. Atom mappings are
computed once per ensemble. The fixed-index clustering method is retained for
callers that explicitly want that metric. Python exposes the new method as
`Mol.cluster_conformers_by_rms_symmetric(conformers, rmsd_threshold=0.5)` and
rejects wrong atom counts, non-finite coordinates and non-finite thresholds.

## Focused fixtures

| Input and coordinate edit | Fixed-index leaders | Symmetric leaders | Independent RDKit `GetBestRMS` for the swapped pair |
|---|---|---|---:|
| `FC(F)(F)C`: swap two F coordinates; third conformer moves one F by 4 Å | `[0, 1, 2]` | `[0, 2]` | 0 Å |
| `CC(=O)[O-]`: swap terminal O coordinates | `[0, 1]` | `[0]` | 0 Å |
| `C[N+](=O)[O-]`: swap terminal O coordinates | `[0, 1]` | `[0]` | 0 Å |

The independent RDKit checks use the same atom-indexed coordinate arrays as
the Rust/Python tests. With `symmetrizeConjugatedTerminalGroups=False`, RDKit
returns 0.2562223321058943 Å for each terminal-oxygen fixture. The original
`rmsd_symmetric` tests also check that isotope labels are respected.

Verification: all 43 `conformer::tests` passed; the locally built Python wheel
passed all 6 tests in `test_rdkit_conformer_rms.py`; `cargo check -p chematic-py`,
targeted Clippy with warnings denied, formatting, and the refactor budget gate
passed. These fixtures do not measure performance, geometry against experiment,
minimum MMFF energy, conformer populations or the published package. The
separate [ETKDG quality record](2026-10-22-etkdg-unfinished-minimum-retry.md)
still has a 0.308065 kcal/mol gap on its already used development slice.
