# Square-planar constrained embedding: development probe

This local probe compares the existing `rdkit_embed` with the new opt-in
`rdkit_embed_square_planar` at the same seeds. It uses six synthetic platinum
examples spanning SP1–SP3 with four ligands or three ligands and a vacancy,
plus the previously exposed iodine and selenium development cases from the
[wheel pilot](2026-11-26-3d-wheel-pilot.md). The eight molecules and seeds 1
and 42 were selected for a small development check, not as an independent
holdout or population sample. The two exposed cases have already influenced
this development.

Each input is parsed and hydrogens are added. Both methods run one attempt per
seed, using the same build and input. A returned structure passes the selected
angle policy only when every declared trans angle is at least 150° and every
measured cis angle is in 60°–120°. These thresholds are the opt-in method's
own acceptance policy, so its successful results pass by construction.
This is not an independent stereochemical or geometry-quality assessment.

| Method | Attempts | Returned | Returned and passed angle policy | Errors |
| --- | ---: | ---: | ---: | ---: |
| Existing embedding | 16 | 10 | 0 | 6 |
| Opt-in square-planar embedding | 16 | 5 | 5 | 11 |

All six existing-method errors were the four-ligand Pt inputs. All 11 opt-in
errors were final angle-policy failures after a constrained embedding; rejected
coordinates are not returned. Neither iodine attempt passed; one selenium
attempt did. Passing this narrow angular screen does not establish bond-length
quality, energy, other declared stereo, general square-planar stereochemistry,
or agreement with crystal structures. These results do not establish RDKit
parity. The separate square-planar angle diagnostic tests compare the measured
angles on RDKit 2026.09.2 coordinates for a small fixed sample.

The host was macOS 27.0.1 arm64 with Python 3.13.6 and `RAYON_NUM_THREADS=1`.
The source-built CPython wheel used release settings with LTO disabled and 16
codegen units. The wheel was not published. Wheel SHA-256:
`ab939282aa2989807e2b7c7a1b285c7acd52f2af36df2f365358d276ec280e05`;
extension SHA-256:
`72fb1cc50141a06318ec5b7bcb0f5268ee7d1f4ada444f883e8f42fefc684440`.
The [probe](../scripts/validation/square_planar_embedding_development.py)
SHA-256 is `09d48a538dc8ceb5aacaae2ec5047436b11eeec3cf2803b827da0534f1ac4e8b`.
Input SMILES, seeds, returned coordinates and individual errors are recorded in
`validation/results/2026-11-27-square-planar-embedding/development.json`.
Run the probe from the repository root with the built wheel installed and
`--wheel PATH --output PATH` to reproduce the format.
