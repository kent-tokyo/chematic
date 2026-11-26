# COD local geometry development comparison

The separate development partition of the reserved COD sample contained 24
records with neutral three-coordinate boron and seven with a center having at
least five heavy-atom neighbors. All 31 were selected without looking at
coordinates or outcomes. The holdout partition was not scored. This is a
development comparison, not a representative or independent estimate.

`cod_local_geometry.py` compares the published chematic 1.2.5 wheel's
`generate_3d()` with RDKit 2026.09.2 ETKDGv3 seed 42 followed by MMFF94 when
available. Each engine had three seconds including process startup and
cleanup on macOS arm64, Python 3.13.6, with one thread. The experimental SDF
and selection split hashes are in the per-record JSON. The script maps the
first substructure match to the reference, then sorts all bond lengths and
all neighbor-pair angles independently at each center. This comparison is
invariant to neighbor ordering but does not validate individual ligand
correspondence, cis/trans assignment, global geometry, or all stereo labels.

| Quantity | Published chematic | Candidate | RDKit |
| --- | ---: | ---: | ---: |
| Molecules with coordinates | 31/31 | 31/31 | 31/31 |
| Boron centers | 36 | 36 | 36 |
| B angle absolute error median / p95 | 10.16° / 49.72° | 2.81° / 18.11° | 3.03° / 9.59° |
| B bond-length absolute error median | 0.149 Å | 0.115 Å | 0.039 Å |
| High-coordinate centers | 8 | 8 | 8 |
| High-coordinate angle absolute error median / p95 | 9.05° / 32.43° | 2.42° / 13.47° | 4.57° / 26.47° |
| High-coordinate bond-length absolute error median | 0.231 Å | 0.249 Å | 0.052 Å |
| Process elapsed median / p95 | 0.636 / 1.581 s | 0.444 / 1.331 s | not compared |

Angle percentiles pool the individual angle errors, so they are not
confidence intervals or molecule-level statistics. Of 36 paired boron
centers, mean angle error improved by more than 0.01° at 35 and worsened at
one. Of eight high-coordinate centers, five improved and three worsened.
The larger boron errors and bond-length gap remain. RDKit reported MMFF
setup unavailable on 26 molecules, unfinished on three, and converged on
two; its unavailable cases are embedded rather than MMFF-optimized
coordinates. The two chematic builds also differ in release build settings,
and these single-run elapsed times include process overhead.

The candidate uses the existing neutral-trigonal-boron and expanded-octet
angle heuristics in the DREIDING path and prepares invariant topology once
instead of rebuilding it during finite-difference energy evaluations. When
two substituents at a trigonal boron coincide in the initial placement, it
projects topology constraints before minimization if declared stereo passes.
It is
not a complete DREIDING parameterization for unsupported elements. The
candidate library SHA-256 is
`0f39f49e3754da4a2a5cb71126a2b21bdd61c46c1c510bff3b32d5cbde697f47`.
An earlier candidate without precomputed nonbonded exclusions timed out on one
input; its successor completed all 31. Separately, 20 selected development
molecules without boron or a five-neighbor center returned exactly the same
coordinates as the published wheel before the coincident-boron repair. This
checks a bounded regression sample rather than all other inputs.

Fixing angle projection globally initially changed the native ETKDG output
for 19 of those 20 molecules. On the existing RDKit 3D stereo checker, three
previously passing records then had new atom-stereo violations and one
previously failing record passed. This candidate was not retained. The native
ETKDG entry point now keeps its earlier angle-projection behavior while the
explicit constraints API and coincident-boron repair use the corrected
projection. On this 20-record sample, its coordinates again matched the
published wheel exactly. The three newly broken examples are protected by
RDKit 3D stereo regression tests. The compatibility choice is provisional:
it preserves known behavior while a stereo-safe ETKDG projection is developed.
After that adjustment, the 31 local-geometry records again produced identical
sorted center angles and lengths to the candidate in the table. The later
candidate library SHA-256 was
`af8b3a615fc163232d6d277f11bf5829d7334b94eef5adf986dfad5f37af3893`.

The same coarse UFF-based heavy-atom screen was applied to generated
structures and their experimental references. Of 31 references, 22 passed.
On those 22 records, 14 published-wheel structures, 17 candidate structures
and 19 RDKit structures passed. Three candidate records improved from failure
to pass and none changed from pass to failure. This screen checks supported
bond lengths, angles and heavy-atom clashes. All bonds and angles were
supported in this sample. It does not establish complete stereochemical or
conformational validity.

The script snapshots, input lists and per-record results are in
`validation/results/2026-11-26-3d-wheel-pilot/`. The source SDF is reconstructed
from `validation/manifests/cod-reserved-2026-11-26-source.json` as described
in the reserved-holdout protocol.
