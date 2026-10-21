# ETKDG sixth-order gradient: exposed conformer development

Date: 2026-10-21. Source candidate only; not a published-package or broad
RDKit superiority claim.

## Design

RDKit 2026.09.1's ETKDG sixth Fourier torsion term has an energy/gradient
inconsistency ([issue #8715](https://github.com/rdkit/rdkit/issues/8715)).
RDKit merged its own correction to master on 2026-10-16
([PR #9699](https://github.com/rdkit/rdkit/pull/9699)); this comparison is
specifically against the installed 2026.09.1 release, not current RDKit master.
CheMatic already offered the exact derivative for
single-molecule embedding. This change exposes that derivative for an ensemble,
then offers a quality profile that runs both the compatibility and corrected
20-candidate ensembles and returns whichever complete ensemble has the lower
finite best MMFF94 energy. The compatibility path is unchanged.

Six exposed development molecules were selected by structure before scoring;
`rdDistGeom.GetExperimentalTorsions` reports a nonzero sixth-order coefficient
of -0.9 for each. No frozen holdout was inspected for this change. All engines
used seed 42, ETKDGv3, MMFF94, and 200 MMFF iterations. RDKit used 40
compatibility candidates; CheMatic's quality profile used 20 compatibility
plus 20 corrected candidates. The work budget is equal in candidate count,
although setup, ETKDG and MMFF runtimes were not measured. Best energy is in
kcal/mol. A difference smaller than 1e-6 is treated as numerical agreement.

| SMILES | RDKit 40 best | CheMatic quality best | CheMatic minus RDKit |
|---|---:|---:|---:|
| `CS(=O)(=O)NCC` | -72.327666 | -72.327666 | < 1e-6 |
| `CCS(=O)(=O)NCC` | -75.820009 | -75.820009 | < 1e-6 |
| `CS(=O)(=O)NCCO` | -57.565141 | -57.565141 | < 1e-6 |
| `CS(=O)(=O)NCCc1ccccc1` | -42.663092 | -42.745568 | -0.082477 |
| `CS(=O)(=O)NCCC` | -72.936518 | -72.936518 | < 1e-6 |
| `CS(=O)(=O)NCCF` | -57.340336 | -57.913871 | -0.573535 |

All 240 RDKit candidates and all 240 CheMatic candidates converged. The
RDKit MMFF94 implementation independently re-scored all 240 CheMatic final
coordinate sets; the maximum same-coordinate energy difference was
5.69e-14 kcal/mol. The
corrected-only 20-candidate ensemble was worse than the compatibility-only
ensemble by 1.515 kcal/mol for `CS(=O)(=O)NCCO`; selecting the better complete
ensemble prevents that regression on the chosen MMFF metric. The quality
profile is thus better than RDKit's 40-candidate run for 2/6 selected rows and
equivalent within 1e-6 for 4/6. It retains 20 conformers, not all 40.

This is a deliberately targeted development slice. It does not test broad
chemical coverage, independent geometry or stereo quality, conformer diversity,
timing, memory, or a fresh held-out set. A broad 3D claim still requires those
gates and comparison against the exact published artifacts.

## Reproduce

```sh
cargo run --release -p chematic-3d --example etkdg_sixfold_ensemble_probe > /tmp/etkdg-sixfold.jsonl
python scripts/etkdg_sixfold_ensemble_probe.py /tmp/etkdg-sixfold.jsonl
```

The Python script requires RDKit 2026.09.1. The Rust example emits every
candidate's status, energy, and coordinates, so the comparison can be
recalculated without
re-running CheMatic. The library regression test pins the corrected seed-42
coordinates to the existing single-conformer API and checks the quality
profile's best-energy selection. The built CPython 3.13 arm64 candidate wheel
had SHA-256
`4ffdc9e04ccf27f4835e2744ea35d86dc03bb696dd81c770048c02998736361`;
the Python embedding test file passed 13/13 with that wheel. The wheel is a
local candidate, not a release artifact.
