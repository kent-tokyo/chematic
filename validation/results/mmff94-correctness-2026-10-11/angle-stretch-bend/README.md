# Angle and stretch-bend correctness, 2026-10-11

The local CheMatic 1.0.20 source candidate now agrees with RDKit 2026.03.6
for angle and stretch-bend energies to below 1e-13 kcal/mol on all 262
comparable frozen inputs. This follows the [O-H/vdW correction](../README.md)
and uses that freshly measured candidate as the before record.

## Changes and causes

1. The angle resolver maintained an incomplete duplicate equivalence table.
   Type 58 was missing, so benzyl-pyridinium angles in inputs 89 and 96 used
   `ka=1.179, theta0=106.327` instead of `1.160, 107.963`. The resolver now
   reads the existing generated numeric-type registry. A missing angle
   classification goes to the empirical rule rather than borrowing type 0.
2. The stretch-bend resolver borrowed type-0 constants after a classified
   table miss. Biphenyl-like connections therefore used -0.411 instead of
   the periodic-row default +0.3. Lookup now preserves the requested
   classification and reaches the periodic-row default on a miss.
3. Equal outer atom types with different bond types require a separate
   orientation tie-break. Butadiene's constants are now correctly oriented
   as `(0.219, 0.250)` or `(0.250, 0.219)` when the triple is reversed.
4. Linear centers use `143.9325 * ka * (1 + cos(theta))`, with its finite
   cosine gradient, and have no stretch-bend term. Nonlinear angles use
   cubic coefficient -0.006981317; angle/stretch-bend unit conversions are
   calculated from 143.9325 and pi/180 instead of rounded decimal factors.
5. The numerical and prepared angle/stretch-bend energy paths share kernels.
   The registry documentation was updated through its generator; generated
   parameter data is unchanged.

Rules were checked against the fixed RDKit release's
[angle kernel](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/ForceField/MMFF/AngleBend.cpp),
[stretch-bend kernel](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/ForceField/MMFF/StretchBend.cpp),
[parameter lookup](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/ForceField/MMFF/Params.h),
and [molecule parameter selection](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp).
The registry retains its existing frozen Release_2026_03_3 provenance.
Live-oracle regression fixtures cover biphenyl, butadiene, a
benzyl-pyridinium motif, acetaminophen, and a deliberately bent nitrile.

The table-only public angle and stretch-bend lookup functions can now return
`None` where they previously borrowed a different classification's value.
Production resolved paths apply the empirical/periodic-row fallback instead.
No API signature or package version was changed.

## Measurements

All 262 paired rows have identical coordinate hashes and RDKit total energies.
The full 265-row accounting is unchanged: 262 comparable, two RDKit embedding
failures, one declared unsupported. No inputs were removed.

Absolute total-energy differences, kcal/mol:

| Measure | Before this step | After this step |
| --- | ---: | ---: |
| Median | 0.005825 | 0.004270 |
| P90, nearest rank | 0.710687 | 0.042537 |
| Maximum | 2.104696 | 0.723842 |
| Within 1 kcal/mol | 242/262 | 262/262 |

Maximum absolute per-term differences after the change:

| Term | kcal/mol |
| --- | ---: |
| Angle | 7.11e-14 |
| Stretch-bend | 1.78e-14 |
| Bond | 4.97e-14 |
| Electrostatic | 1.14e-13 |
| Torsion | 0.708000 |
| Out-of-plane | 0.172758 |
| vdW | 0.079247 |

144 absolute total differences improved and 118 worsened, with a 1e-8
kcal/mol tolerance. Removing angle/stretch-bend errors can expose other
terms previously cancelled in the total. Per-term deltas are retained in
[candidate.jsonl](candidate.jsonl), so the total is not used as proof that
every molecule improved in every respect.

The largest remaining total error is input 85, 0.723842 kcal/mol, including
torsion +0.708000 and vdW +0.017069. Several other large residuals share the
same +0.708000 torsion offset. Torsion is the next numerical-correctness
target, followed by out-of-plane and vdW residuals.

Analytic versus central-difference gradients were checked on inputs
89, 96, 166, 168, and 231 with delta 1e-5 angstrom. Maximum component error
was 1.33e-7 kcal/mol/angstrom. A Rust test also checks bent and exactly
straight nitrile gradients, including finite values at 180 degrees.
These checks establish internal gradient consistency, not whole-force-field
RDKit gradient equality.

## Geometry regression and tests

The stereo-safe strict-MMFF94 single-conformer run succeeded on all 265
inputs. The common scorer marked all 265 geometrically sound with no gross
clashes, stereo violations, or unevaluable declared stereo. The configuration
matches the preceding run: seed 20260801, eight attempts, 200 force-field
iterations, 20-second timeout. Scoring uses CheMatic geometry/stereo routines;
it is not an independent QM or experimental reference.

- `cargo test -p chematic-ff --offline -j 2 --quiet`: 238 passed,
  one ignored doc test.
- `cargo test -p chematic-3d --offline -j 2 --quiet`: 608 passed,
  ten ignored unit tests and four ignored doc tests.
- `cargo clippy -p chematic-ff --all-targets --offline -j 2 -- -D warnings`: passed.
- Diagnostic, pipeline-v2, and ensemble-v2 Python tests: 45 passed against
  the rebuilt extension. The two new parameter fixtures and the live
  energy fixture failed before the fixes.
- `cargo fmt -p chematic-ff --check` and `git diff --check`: passed.

## Evidence and reproduction

[candidate-summary.json](candidate-summary.json) records the complete
energy command, build command, extension/oracle hashes, versions, and frozen
corpus hashes. The new extension SHA-256 is
`4f7eaec2fa2cfc5305a9c2ccefa82db8a66ec83fbe4bd03f359db751d3615d7d`;
the preceding candidate's hash is retained in its own summary.
[comparison.json](comparison.json) records paired counts, remaining outliers,
relevant file hashes, and geometry artifact hashes.
[implementation.patch](implementation.patch) is relative to HEAD and includes
the earlier changes to the shared Python test file. The preceding O-H/vdW
implementation patch remains in the parent directory.

This is an uncommitted dirty-tree source candidate, not a published package.
The energy command correctly returned exit 2: all recorded gates passed
except `source_tree_clean`. Relevant source files form a partial snapshot;
other concurrent working-tree changes are excluded, so HEAD alone does not
reproduce the build. The installed RDKit artifact hash pins `rdBase.so`, not
every RDKit binary or a complete wheel. Speed and ensemble coverage were not
measured.

To repeat the geometry check in the rebuilt source-extension environment,
run from the repository root:

```sh
python validation/results/mmff94-correctness-2026-10-11/run_quality.py \
  --output /tmp/mmff94-angle-quality-repeat.jsonl
target/debug/examples/pipeline_v2_vs_rdkit_common_scorer \
  --chematic-rows /tmp/mmff94-angle-quality-repeat.jsonl \
  --rdkit-rows validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.jsonl \
  --output /tmp/mmff94-angle-quality-repeat-scored.jsonl
```

Historical RDKit geometry rows are used only for the scorer's second lane;
the totals reported here refer to the fresh CheMatic candidate rows. There
is no new speed or ensemble-quality comparison against RDKit.
