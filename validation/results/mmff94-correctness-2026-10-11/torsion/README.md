# Torsion correctness, 2026-10-11

The local CheMatic 1.0.20 source candidate agrees with RDKit 2026.03.6
torsion energies within 7.2e-14 kcal/mol on all 262 comparable frozen inputs.
This step follows the [angle/stretch-bend correction](../angle-stretch-bend/README.md)
and uses that freshly measured candidate as its baseline.

## Changes

1. Torsion enumeration now requires four distinct atoms. A three-membered
   ring's closed path `i-j-k-i` previously contributed a spurious constant:
   0.708 kcal/mol for cyclopropane, 0.670 for oxirane, and 0.830 for aziridine
   in the explicit-H regression fixtures. Genuine torsions involving ring
   atoms and hydrogens remain present.
2. The additive `mmff94_torsion_energy_with_bond` API tries the existing
   parameter table, then the ordered Halgren empirical rules using central
   atom metadata and the actual central bond. This resolves the missing
   iminium-N/aromatic-C barrier without adding molecule-specific constants
   to the table. Carbon, nitrogen, oxygen, silicon, phosphorus, and sulfur
   use the defined empirical constants. Unknown types/elements, terminal
   centers, and unsupported central bond kinds remain unresolved by this
   empirical path. Linear centers retain explicit no-term-by-design
   accounting; other empirical zero barriers count as resolved.
3. The numerical and prepared torsion energy paths share enumeration and
   the energy kernel. The 3D force-field coverage report uses the same
   four-atom constraint and parameter resolver, with regression tests for
   triangle/square/chain counts and iminium coverage.

Reference implementations were inspected at RDKit's fixed release:
[torsion enumeration](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/ForceFieldHelpers/MMFF/Builder.cpp)
and [empirical parameter selection](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp).
The atom metadata registry retains its existing Release_2026_03_3 provenance.
The table-only lookup API is unchanged. This step does not establish complete
torsion parameter-lookup parity outside the measured corpus.

## Same-coordinate measurements

All 265 input statuses are unchanged: 262 comparable, two RDKit embedding
failures, and one declared unsupported. All paired coordinate hashes and
RDKit total energies match the baseline exactly; frozen manifest hashes
are identical. No input was removed.

Absolute total-energy differences, kcal/mol:

| Measure | Before this step | After this step |
| --- | ---: | ---: |
| Median | 0.004270 | 0.003907 |
| P90, nearest rank | 0.042537 | 0.036463 |
| Maximum | 0.723842 | 0.125920 |
| Within 1 kcal/mol | 262/262 | 262/262 |

Maximum absolute per-term differences after the change:

| Term | kcal/mol |
| --- | ---: |
| Torsion | 7.11e-14 |
| Bond | 4.97e-14 |
| Angle | 7.11e-14 |
| Stretch-bend | 1.78e-14 |
| Electrostatic | 1.14e-13 |
| Out-of-plane | 0.172758 |
| vdW | 0.079247 |

Torsion energies changed in 17 inputs: eleven contained closed triangle
paths and six had missing iminium barriers. Twelve absolute total errors
improved, five worsened, and 245 were unchanged within 1e-8 kcal/mol.
The five increases expose errors previously cancelled by missing torsions.
Input 99 now has the largest total difference, 0.125920 kcal/mol; its
out-of-plane difference is -0.172758 and vdW difference is +0.046838.
Both remaining term classes need further work.

Analytic versus central-difference total gradients were checked on inputs
85, 99, 137, 142, and 254, using delta 1e-5 angstrom. Maximum component
error was 8.93e-8 kcal/mol/angstrom. This checks internal energy/gradient
consistency, not equality with the full RDKit gradient.

## Regression checks

- Rust: 242 force-field tests and 610 3D tests passed. One force-field
  doc test, ten 3D unit tests, and four 3D doc tests remain ignored.
- Python: 53 diagnostic, pipeline-v2, and ensemble-v2 tests passed against
  the rebuilt extension.
- `cargo clippy -p chematic-ff --all-targets --offline -j 2 -- -D warnings`,
  scoped formatting, and `git diff --check` passed.
- Before the fix, the new triangle enumeration test failed, and six of
  eight live-oracle torsion fixtures failed; the square-ring and butane
  controls passed. The iminium fixtures are twisted to 70 degrees to
  expose the missing barrier away from its nearly planar minimum.
- The strict-MMFF94 stereo-safe geometry run succeeded on all 265 inputs.
  The common scorer marked all 265 geometrically sound, with no gross
  clashes, stereo violations, or unevaluable declared stereo. Configuration:
  seed 20260801, eight attempts, 200 minimizer iterations, 20-second timeout.
  The scorer uses CheMatic geometry/stereo routines, not an independent QM
  or experimental reference. Its historical RDKit rows only populate the
  second scoring lane; the reported 265 results are fresh CheMatic rows.

## Evidence and reproduction

[candidate-summary.json](candidate-summary.json) records the complete energy
and build commands, versions, corpus hashes, and artifact hashes.
The extension SHA-256 is
`747487b6474ee9ba452f50f802c3c981eeba10cddfb774281f0e710eb3d57e77`;
it remained unchanged through energy measurement, Python tests, and the
geometry run. The installed RDKit artifact hash pins `rdBase.so`, not every
RDKit binary or an archived wheel.

[comparison.json](comparison.json) contains the paired comparisons, changed
inputs, remaining outliers, relevant source hashes, and geometry artifacts.
[implementation.patch](implementation.patch) is a partial snapshot relative
to HEAD, including preceding MMFF corrections. Other concurrent workspace
changes are excluded; HEAD plus this patch is not a complete build snapshot.
The build and verification logs are retained alongside these results.

This is an uncommitted dirty-tree source candidate, not a published package.
The energy command returned exit 2 as intended: every recorded gate except
`source_tree_clean` passed. Speed, conformer-ensemble coverage, and independent
structure accuracy were not measured.

To repeat the geometry check in the rebuilt source-extension environment:

```sh
python validation/results/mmff94-correctness-2026-10-11/run_quality.py \
  --output /tmp/mmff94-torsion-quality-repeat.jsonl
target/debug/examples/pipeline_v2_vs_rdkit_common_scorer \
  --chematic-rows /tmp/mmff94-torsion-quality-repeat.jsonl \
  --rdkit-rows validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.jsonl \
  --output /tmp/mmff94-torsion-quality-repeat-scored.jsonl
```
