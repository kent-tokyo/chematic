# Water and conjugated nitrogen cations, 2026-10-11

This local CheMatic 1.0.20 source candidate corrects MMFF94 typing and
charge assignment for water, amidinium and guanidinium, closing the
[two previously recorded typing gaps](../oop-vdw/known-holdout-typing-gaps.json).
The separate
20-input regression set includes resonance spellings, substituted and
connected cationic groups, nitro-substituted amidinium, neutral guanidine,
and alcohol/ether/hydroxide/oxonium/oxenium controls. These development
regressions are outside the existing frozen 265-input diagnostic; they
are not an unseen test set or a change to that diagnostic's denominator.

## Implementation

- Water oxygen uses OH2 (70), with explicit or stored implicit hydrogens.
  The ordered degree/valence tests also distinguish oxonium O+ (49) and
  oxenium O=+ (51). Their explicit H types are 31, 50 and 52 respectively;
  hydroxide H remains ordinary hydroxyl H (21).
- Amidinium/guanidinium carbon uses CNN+/CGD+ (57). The participating
  nitrogen types are NCN+ (55) and NGD+ (56), independent of which nitrogen
  bears the input positive charge. Nitro/nitrate N retains priority over
  the conjugated N rule. Neutral guanidine retains types 3/9/40.
- Nitrogen-bound H now follows its assigned heavy atom type, including
  HNN+/HGD+ (36), rather than duplicating the heavy-atom structural tests.
- Formal charge is averaged over the complete connected component of
  nitrogen types 55/56/81 reached through carbon types 57/80. This handles
  a shared nitrogen between two cationic groups without double-counting
  its charge. It does not complete the separate, pre-existing aromatic
  imidazolium carbon typing gap.

The rules were checked against RDKit's fixed-release
[AtomTyper.cpp](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/ForceFieldHelpers/MMFF/AtomTyper.cpp):
heavy C/N/O assignment, `setMMFFHydrogenType`, and `computeMMFFCharges`.
Parameter tables and public API signatures were unchanged. Rust regression
fixtures preserve independent RDKit 2026.03.6 type/charge values and check
both atom-order and bond-insertion-order reversal.

## Validation

See [comparison.json](comparison.json) for the final measurements and
artifact/source SHA-256 values. The frozen corpus remains 265 inputs:
262 comparable, two RDKit embedding failures and one declared unsupported.
The manifest, protocol, coordinate hashes, reference energies and candidate
energies are unchanged from the [preceding candidate](../oop-vdw/README.md).
Maximum total-energy error remains 9.095e-13 kcal/mol. The direct gradient
comparison covers all 37,407 Cartesian components in the same 262 inputs.

[Separate regression measurements](candidate-holdout.json) compare all
20 inputs at shared explicit-H coordinates, with a deterministic coordinate
perturbation to expose nonzero energy/force terms. Maximum absolute total
energy error is 2.843e-13 kcal/mol, maximum partial-charge error 2.221e-16 e,
and maximum analytic-gradient component error 3.127e-13 kcal/mol/angstrom.
Every per-atom numeric type matches. Scaled internal analytic/central-
difference gradient error is at most 6.444e-8 with a 1e-5 angstrom step.
The script also records all seven energy terms and coordinate hashes.

Fresh 3D generation on the frozen 265-input set produces the same coordinates
as the preceding candidate. All 265 pass the project's common geometry/stereo
scorer, with zero gross clashes, violated or unevaluable declared stereo.
The optimizer reports convergence for 88/265 at 200 iterations, unchanged.
Geometry soundness and optimizer convergence are separate results. This
pipeline profile generally operates on implicit-H molecular graphs; it is
separate from the explicit-H fixed-coordinate force-field comparison.

Rust: 860 passed, 15 ignored across `chematic-ff` and `chematic-3d`.
Python: 92 passed across the force-field diagnostic, pipeline v2, and both
ensemble suites. Strict FF Clippy, scoped formatting and diff checks pass.
The two new Rust test functions cover 20 explicit-H fixtures in four
relabeling/insertion configurations and ten implicit-H typing controls.

## Reproduction and evidence boundary

From the repository root, with the recorded development environment:

```sh
VIRTUAL_ENV=/private/tmp/chematic-mmff-development-20261011 .venv-v1gate/bin/maturin develop --release --offline -j 2 -m crates/chematic-py/Cargo.toml -- -C strip=none -C debuginfo=1
cargo test -p chematic-ff -p chematic-3d --offline -j 2 --quiet
cargo clippy -p chematic-ff --all-targets --offline -j 2 -- -D warnings
/private/tmp/chematic-mmff-development-20261011/bin/python -m pytest scripts/tests/test_mmff94_same_explicit_h_energy.py crates/chematic-py/tests/test_pipeline_v2.py crates/chematic-py/tests/test_conformer_ensemble.py crates/chematic-py/tests/test_conformer_ensemble_v2.py -q
PYTHONPATH=. /private/tmp/chematic-mmff-development-20261011/bin/python validation/results/mmff94-correctness-2026-10-11/cation-water/measure_holdout.py candidate
/private/tmp/chematic-mmff-development-20261011/bin/python validation/results/mmff94-correctness-2026-10-11/cation-water/compare_gradients.py
/private/tmp/chematic-mmff-development-20261011/bin/python validation/results/mmff94-correctness-2026-10-11/cation-water/compare.py
```

The exact frozen-energy command is in [candidate-summary.json](candidate-summary.json).
Its exit status 2 records the intentionally dirty source tree; all other
provenance/accounting gates pass. The extension and `rdBase.so` are hashed
before and after measurements. The latter pins that installed binary only,
not the entire RDKit distribution.

[implementation.patch](implementation.patch) is cumulative relative to the
recorded HEAD and includes earlier force-field work. It excludes unrelated
concurrent work and is not a complete build snapshot. This is uncommitted
source-candidate evidence; no package was released. Agreement with RDKit's
MMFF94 implementation does not establish superior physical accuracy,
conformer quality, search coverage or runtime speed. Optimizer convergence,
remaining typing/charge gaps and conformer search remain development work.
