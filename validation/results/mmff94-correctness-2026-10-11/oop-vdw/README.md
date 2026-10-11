# Out-of-plane and long-range vdW correctness, 2026-10-11

On the frozen 265-input diagnostic, the local CheMatic 1.0.20 source
candidate now matches RDKit 2026.03.6 total MMFF94 energies within
9.1e-13 kcal/mol on all 262 comparable inputs. This follows the
[torsion correction](../torsion/README.md), which supplies the fresh baseline.
The numerical result applies to these inputs and shared coordinates;
known typing gaps outside this corpus and optimizer convergence remain.

## Changes and causes

1. Each trigonal center requires three Wilson out-of-plane angles, one
   per substituent relative to the plane of the other two. Only one was
   evaluated before, making the result depend on neighbor insertion order.
   The builder now enumerates all three, and numerical/prepared paths share
   that enumeration and kernel. Multiplying a single angle by three would
   also be incorrect for asymmetric geometries.
2. OOP parameter lookup now searches the generated MMFFDEF outer-type
   equivalence levels before the final wildcard default. The previous
   exact/wildcard resolver selected wrong parameters at 92 centers across
   six distinct type combinations in the comparable corpus. A
   [source-derived lookup audit](oop-parameter-audit.json) checked 3,556
   degree-three centers against the live oracle: 92 mismatches for the
   previous lookup model, zero for the equivalence model. This audit models
   lookup in Python; energy measurements and regression tests exercise Rust.
3. The eligible-center list now matches every central type in the OOP
   table, including N-oxide type 69 and cationic nitrogen type 55. Zero
   parameter rows are included. The coverage report still counts parameter
   centers, each supplying three angles. Unknown outer types remain
   unresolved rather than being hidden by a wildcard parameter.
4. OOP energy and gradient use the exact 143.9325 conversion with angles
   in radians, replacing the rounded 0.043844 degree conversion.
5. The ordinary vdW evaluator and its analytic gradient now evaluate all
   topology-prepared pairs. A hard 10-angstrom cutoff had been applied even
   in the documented no-cutoff path, dropping distant attractive terms.
   The explicitly selected cutoff energy/gradient and neighbor-list path
   retain their 10-angstrom boundary. The numerical minimizer's vdW path
   also retains long-range pairs.

The rules were checked against RDKit's fixed-release
[OOP enumeration](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/GraphMol/ForceFieldHelpers/MMFF/Builder.cpp),
[OOP kernel](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/ForceField/MMFF/OopBend.cpp),
and [parameter lookup](https://github.com/rdkit/rdkit/blob/Release_2026_03_6/Code/ForceField/MMFF/Params.h).
The registry retains its existing Release_2026_03_3 provenance. No package
version or API signature was changed; default vdW behavior and the public
eligible-center inventory were corrected.

## Fixed-coordinate energy measurements

All 265 statuses are unchanged: 262 comparable, two RDKit embedding
failures, and one declared unsupported. Every comparable coordinate hash
and RDKit total energy is identical to the preceding candidate. Manifest
hashes are unchanged and no input was removed.

Absolute total-energy differences, kcal/mol:

| Measure | Before this step | After this step |
| --- | ---: | ---: |
| Median | 0.003907 | 8.53e-14 |
| P90, nearest rank | 0.036463 | 2.91e-13 |
| Maximum | 0.125920 | 9.09e-13 |
| Within 1 kcal/mol | 262/262 | 262/262 |

Maximum absolute per-term differences after the change:

| Term | kcal/mol |
| --- | ---: |
| Bond | 4.97e-14 |
| Angle | 7.11e-14 |
| Stretch-bend | 1.78e-14 |
| Torsion | 7.11e-14 |
| Out-of-plane | 2.78e-16 |
| vdW | 2.84e-13 |
| Electrostatic | 1.14e-13 |

At a 1e-8 kcal/mol change tolerance, 205 total errors improved, zero
worsened, and 57 were unchanged. OOP changed in 173 inputs and vdW in 197.
These are numerical agreement results on RDKit-generated coordinates,
not an independent assessment of the physical accuracy of MMFF94.

## Gradients, structures, and tests

[Direct gradient comparison](gradient-comparison.json) re-created and
verified the exact coordinate hashes before comparing all 37,407 Cartesian
components in 262 inputs with RDKit `CalcGrad`. Maximum absolute difference
was 3.8111e-7 kcal/mol/angstrom; maximum scaled difference was 7.6091e-8,
using `abs(actual-reference)/(1+abs(reference))`. Both passed the preset
1e-6 absolute and 1e-7 scaled thresholds. The three noncomparable inputs
are retained as explicit not-compared rows. The reproducible script is
[compare_gradients.py](compare_gradients.py).

Internal central-difference gradients were also checked on inputs 91, 94,
99, 192, 222, and 231, with delta 1e-5 angstrom. Maximum component error
was 1.4112e-7 kcal/mol/angstrom. Five deliberately nonplanar regression
motifs compare complete analytic gradients directly with RDKit as well.

- Rust: 247 force-field tests and 611 3D tests passed. One force-field doc
  test, ten 3D unit tests, and four 3D doc tests remain ignored.
- Python: 59 diagnostic, pipeline-v2, and ensemble-v2 tests passed against
  the rebuilt extension.
- Strict all-target force-field Clippy, scoped formatting, and
  `git diff --check` passed.
- Before the implementation fix, two new force-field unit tests, two OOP
  parameter tests, the N-oxide coverage test, and all six live-oracle
  OOP/vdW fixtures failed. Their logs are retained.
- The unchanged stereo-safe strict-MMFF94 profile generated all 265
  structures, and the common scorer accepted all 265 as geometrically
  sound with no gross clashes, stereo violations, or unevaluable declared
  stereo. All reported the requested force field without fallback.

Geometry acceptance does not imply optimizer convergence: the reported
converged count was 88/265, versus 90/265 before this step, at the fixed
200-iteration budget. Improving that remains work. The geometry scorer
uses CheMatic routines rather than independent QM/experimental references.
Its historical RDKit rows only fill the second scoring lane. There is no
new speed, RDKit convergence, or ensemble-quality comparison.

## Scope, evidence, and next work

The fixed corpus does not cover all MMFF94 typing chemistry. Two freshly
checked external cases are retained in
[known-holdout-typing-gaps.json](known-holdout-typing-gaps.json): water oxygen
is assigned type 6 instead of 70, and guanidinium
`NC(=[NH2+])N` does not receive the oracle's 56/57 family and attached-H
types. They were never part of this fixed 265-input score. Correcting
those families and improving convergence are the next development targets.

[candidate-summary.json](candidate-summary.json) records complete energy
and build commands, versions, frozen manifest hashes, and extension/oracle
hashes. The extension SHA-256 is
`cf6f1e7f50122880dfc7160e52946fd8da0d55caee7f9daad767b44ce4a216d8`;
it remained unchanged through energy measurements, tests, geometry, and
direct gradient comparison. The RDKit artifact hash pins `rdBase.so`, not
every RDKit binary or an archived wheel.

[comparison.json](comparison.json) records paired energy counts, source
and artifact hashes, geometry results including convergence, and logs.
[implementation.patch](implementation.patch) is a partial snapshot relative
to HEAD, including preceding MMFF corrections. Other concurrent dirty-tree
changes are excluded; HEAD plus this patch is not a complete build snapshot.
This is an uncommitted source candidate, not a published package. The energy
command intentionally returned exit 2 because `source_tree_clean` is false;
all its other recorded gates passed.

Reproduce the additional parameter/gradient checks from the repository root
in the rebuilt source diagnostic Python environment:

```sh
python validation/results/mmff94-correctness-2026-10-11/oop-vdw/audit_oop_parameters.py
python validation/results/mmff94-correctness-2026-10-11/oop-vdw/compare_gradients.py
python validation/results/mmff94-correctness-2026-10-11/run_quality.py \
  --output /tmp/mmff94-oop-vdw-quality-repeat.jsonl
target/debug/examples/pipeline_v2_vs_rdkit_common_scorer \
  --chematic-rows /tmp/mmff94-oop-vdw-quality-repeat.jsonl \
  --rdkit-rows validation/results/public-package-3d-rdkit-2026.3.6-2026-09-22.jsonl \
  --output /tmp/mmff94-oop-vdw-quality-repeat-scored.jsonl
```
