# Integration with concurrent main update

While PR #827 waited for required CI, PR #826 entered main as
c74391a2 (full revision recorded in post-main-merge.json).
It included the original 0bae2e0b history and new checked IO functionality.

Merge that main into this candidate, retaining all new main IO features.
The MMFF94/3D files match candidate 8c773ba7 except two documentation lines
describing out-of-plane coverage centers. Preserve the new optional diagnostic
argument for compatibility; per-term energies remain always available and
duplicate evaluation/output is removed. Add a regression for both argument values.

Post-integration local checks: 959 Rust unit/integration tests passed, four ignored;
47 Python diagnostic tests passed; scoped formatting passed. The Rust rerun
downloaded missing cached dependencies and did not use --offline. The Python
diagnostic uses the previously built candidate extension: its FF numerical
implementation is unchanged, and this is not a newly rebuilt IO binding check.

The 95-test Python run, wheel, shared-coordinate energy and 265-input geometry
measurements elsewhere in this directory refer to candidate 8c773ba7 before the
IO merge. Preserve them and their manifest hashes unchanged.
GitHub PR #827 records CI for the merged candidate separately.
