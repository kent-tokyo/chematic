# Conformer symmetry RMSD: terminal groups and isotopes

This source change improves the symmetry-aware RMSD used by conformer
deduplication. It does not change 3D embedding, force-field energies,
minimization, or the A6 conformer-quality counts.

## Reproduction

- Source base: `7d773b72` (`origin/main` at branch creation).
- Oracle: RDKit `2025.09.3` on macOS arm64, CPython 3.13.
- Command: `cargo run -p chematic-3d --example rmsd_symmetric_oracle_dump`
  followed by `python3 scripts/rmsd_symmetric_oracle_check.py` on its JSONL.
- Result: 8/8 cases agreed with RDKit `GetBestRMS` within 0.001 Å.
  Acetate and nitro now agree under the default terminal-group
  symmetrization. The isotope-labelled isobutane probe remains distinct
  (`0.164695 Å` in both implementations).
- On the same acetate input, the previous main returned `0.356340 Å`;
  the candidate and RDKit return approximately zero.
- Rust library tests: 706 passed, 4 ignored, 0 failed.

The fixed coordinates and molecular graphs are exposed development cases,
not an independent holdout. This result establishes a narrow RMSD and
deduplication correction; it does not establish overall superiority to RDKit
for 3D geometry, force-field accuracy, or conformer quality.
