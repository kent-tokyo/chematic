# MMFF94 correctness development, 2026-10-11

Local CheMatic 1.0.20 source candidate at HEAD
`7f84ae4a98399086ce9f389c099f3fb77b66066e`, compared with installed
RDKit 2026.03.6 on macOS arm64 / Python 3.13.6. Changes are uncommitted.
This is a numerical-correctness and single-conformer regression diagnostic.
It does not establish superiority in conformer coverage, accuracy, or speed.

## Implementation

- Assign explicit O-bound hydrogen types from their chemical environment:
  ordinary hydroxyl/N-O-H 21, carboxylic/phospho acid 24,
  phenol/enol/HO-C=N 29, sulfur-bound hydroxyl 33. Water H remains 31.
- Suppress the vdW radius asymmetry correction when either partner is a
  donor. Apply the existing donor/acceptor scaling afterward.
- Add optional seven-term RDKit energy diagnostics with
  `--rdkit-term-breakdown`. Verify the indexed bond graph, as well as element
  order, before transferring coordinates.
- Add regression fixtures for H types/charges, symmetric vdW combining,
  live RDKit vdW energy, isolated-term totals, and invalid coordinate mapping.

The rules were checked against RDKit's pinned
[hydrogen typer](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp)
and [nonbonded kernel](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/ForceField/MMFF/Nonbonded.cpp),
and the installed oracle. Implicit-H force-field typing and the heavy-atom
typing of isolated water were not addressed in this change.

## Same-coordinate energies

Baseline: the historical 2026-09-23 v1.0.19 diagnostic, not a newly built
prechange baseline. Candidate: the source extension built during this run.
All 262 comparable rows have identical coordinate hashes and identical RDKit
total energies between the two records. All 265 terminal statuses agree:
262 comparable, two RDKit embedding failures, one declared unsupported.

Absolute total-energy differences, in kcal/mol:

| Measure | Historical baseline | Current candidate |
| --- | ---: | ---: |
| Median | 0.217545 | 0.005825 |
| P90, nearest rank | 1.144679 | 0.710687 |
| Maximum | 9.874153 | 2.104696 |
| Within 1 kcal/mol | 232/262 | 242/262 |
| Within 5 kcal/mol | 260/262 | 262/262 |

The two previous >5 kcal/mol outliers improved from 9.874153 to 0.501320
(input 166), and from 6.505362 to 0.047858 (input 231).
Across paired rows, 163 absolute total differences improved, 70 were
unchanged, and 29 worsened (1e-8 kcal/mol comparison tolerance).
Total-energy cancellation means aggregate improvement is not proof of
per-term equality for every molecule.

The largest remaining total difference is input 89, 2.104696 kcal/mol:
angle +1.173394 and stretch-bend +0.911031 dominate. Maximum absolute
per-term differences across the corpus are angle 1.289958,
stretch-bend 1.409731, torsion 0.708000, oop 0.172758, vdW 0.079247;
bond and electrostatic differences are below 2e-13 kcal/mol on these inputs.
These are priorities for the next numerical-correctness step.

Analytic versus central-difference gradients were checked on inputs 166
and 231, with delta 1e-5 angstrom. Maximum component error was
1.30e-7 kcal/mol/angstrom. This checks internal energy/gradient consistency,
not RDKit gradient agreement.

## Geometry and stereo regression

The fixed A+B 265-molecule corpus was run with
`PipelineV2Config.stereo_safe`, force field `mmff94_bond_angle_strict`,
seed 20260801, eight embedding attempts, 200 force-field iterations, and a
20-second per-molecule timeout. The other configuration fields come from
`scripts/public_package_3d_chematic.py::pipeline_config`.

All 265 inputs succeeded. The separate common scorer found all 265
geometrically sound, without gross clashes, stereo violations, or
unevaluable declared stereo. The scorer uses CheMatic geometry/stereo
routines, so this is a regression check rather than an independent
experimental or QM reference. It evaluates a single conformer per input;
ensemble coverage was not measured. Timing was collected during concurrent
work and is not used for a speed claim.

## Validation

- `cargo test -p chematic-ff --offline -j 2`: 234 passed; one ignored doc test.
- `cargo test -p chematic-3d --offline -j 2 --quiet`: 608 passed, ten ignored.
  Doc tests were rerun separately after a temporary disk-space failure:
  `cargo test -p chematic-3d --doc --offline -j 2 --quiet`, four ignored.
- `cargo clippy -p chematic-ff --all-targets --offline -j 2 -- -D warnings`: passed.
- Python diagnostic tests: six passed against the rebuilt extension.
- Python pipeline-v2 and ensemble-v2 tests: 34 passed.
- `cargo fmt -p chematic-ff --check` and `git diff --check`: passed.

## Artifacts and reproduction

[candidate-summary.json](candidate-summary.json) records the exact energy
command, build command, versions, manifest hashes, compiled-extension hash,
and installed RDKit `rdBase.so` hash. The latter identifies an installed
shared library, not an archived RDKit wheel or every RDKit binary.
[candidate.jsonl](candidate.jsonl) contains the same-coordinate results and
per-term deltas. [comparison.json](comparison.json) pins the historical
baseline, relevant source files, and geometry artifacts.
[implementation.patch](implementation.patch) captures only this change's
implementation and tests; other concurrent changes are excluded.

The source tree was dirty. The diagnostic therefore correctly exits 2 with
`gate_passed: false` solely because `source_tree_clean` is false. No release
gate or published package is claimed. Relevant source hashes and the patch
are a partial snapshot; HEAD alone cannot regenerate this dirty build.

The successful local build used:

```sh
VIRTUAL_ENV=/private/tmp/chematic-mmff-development-20261011 \
  .venv-v1gate/bin/maturin develop --release --offline -j 2 \
  -m crates/chematic-py/Cargo.toml -- -C strip=none -C debuginfo=1
```

In a Python environment with that source extension and the pinned RDKit
oracle installed, repeat the geometry run and common scoring from the
repository root:

```sh
python validation/results/mmff94-correctness-2026-10-11/run_quality.py \
  --output /tmp/mmff94-quality-repeat.jsonl
target/debug/examples/pipeline_v2_vs_rdkit_common_scorer \
  --chematic-rows /tmp/mmff94-quality-repeat.jsonl \
  --rdkit-rows validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.jsonl \
  --output /tmp/mmff94-quality-repeat-scored.jsonl
```

The RDKit rows above only supply the scorer's second engine lane. This
report makes no new timing or quality comparison against that historical
lane. The CheMatic geometry totals come from the fresh candidate rows.

## Subsequent development

The latest [optimizer recovery and coordinate JSON report](optimizer-recovery/README.md)
records the next source candidate: bounded descending recovery steps,
final-state convergence flags and exact scorer coordinate decoding.
Its fixed-coordinate force-field results retain the original 265-input
protocol; pipeline convergence and stereo measurements are reported separately.
