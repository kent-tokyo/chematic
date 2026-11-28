# Near-planar carbon with an implicit hydrogen

The opt-in tetrahedral geometry screen previously used ligand-tetrahedron
weights and 45–165° angle limits. An exposed synthetic example showed those
checks could accept `F[C@H](Cl)Br` coordinates for which RDKit 2026.09.2 did
not assign a tetrahedral 3D tag. The regression test retains that exact example.
The new additional screen computes the absolute scalar triple product of the
three normalized real C–neighbor bond directions. It is zero for coplanar
bond directions. For carbon with exactly three graph neighbors and one implicit
hydrogen, the screen requires a value of at least 0.15. This threshold is an
experimental margin, not an RDKit-derived rule or a physical-validity test.
It also applies when selecting alternative stereo repairs. The ordinary
non-screened stereo verifier is unchanged.

## Seeded synthetic comparison

`synthetic.py` generates three independent random unit directions and bond
lengths uniformly in 0.8–2.0 Å for `F[C@H](Cl)Br`. The first seed (20261128,
30,000 geometries) was used while choosing the threshold; the second seed
(920261127, 50,000 geometries) was used afterwards as an independent synthetic
check. Both use the same single molecule and coordinate distribution. This is
not an independent chemical corpus or the reserved molecular holdout.

The script compares the previous screen and new screen with RDKit's 3D chiral
tag on the same coordinates. Counts below refer to cases the respective screen
would report as satisfied or violated; all other cases were unevaluable.

| Synthetic seed | RDKit unassigned, reported before | RDKit unassigned, reported after | RDKit-agreeing, reported before | RDKit-agreeing, reported after |
|---|---:|---:|---:|---:|
| Threshold selection, 30,000 | 12 | 0 | 13,469 | 13,347 |
| Independent seed, 50,000 | 23 | 0 | 22,608 | 22,410 |

No raw-vs-RDKit disagreements occurred among RDKit-assigned geometries in this
single-molecule synthetic sample. Zero retained RDKit-unassigned cases here does
not establish that 0.15 catches all uncertain structures. The screen remains
conservative, withholds many other decisions, and does not prove stereochemical
or physical correctness. The seed-labelled JSON outputs retain the complete
contingency counts. Reproduce once with each wheel installed:

```sh
python validation/results/2026-11-27-implicit-h-planarity/synthetic.py > output.json
```

## Exposed development corpus

The same `../2026-11-27-tetrahedral-geometry-screen/run.py` protocol was rerun
on its exposed 1,345-molecule development input. Of its 2,055 centers where
both the raw verifier and RDKit returned comparable tetrahedral decisions,
20 previously retained RDKit-agreeing centers became unevaluable. The other
raw and screened statuses were unchanged, including all 190 raw disagreements
already screened out. Thus 955 decisions were retained and agreed, while 1,100
were withheld, of which 910 raw decisions agreed and 190 disagreed with RDKit.
The remaining 39 RDKit/chematic non-comparable centers are also recorded.
`development-after.json` retains the full output; the previous output is in the
preceding directory. This is development evidence, not a holdout evaluation.

The 576-case `../2026-11-27-stereo-repair-alternatives/run.py` evaluation was
also rerun: the acceptance counts stayed at 206 unchanged baseline successes,
2 recovered cases and 368 failures. This compares acceptance, not every
coordinate or energy; details are in `alternatives-after.json`.

Baseline is commit `aa451f9f`, local debug wheel SHA-256
`58ff82ddb167bdec2036d8b275099d7a046d58c5ca08fc4b2f1ec9ebbe5163c9`.
The candidate wheel uses that commit plus the accompanying code changes,
SHA-256 `dd1fc9a0527dd3ccbe7add431e646c44d9ab9fed143f6dc27b6af55a4c16ec50`.
Both were built locally on macOS arm64 for CPython 3.13.6; neither is a
published wheel. Synthetic tests used NumPy 2.5.3 and RDKit 2026.09.2.
