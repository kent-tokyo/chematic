# Platinum Diverse development baseline

The 2017_01 Platinum Diverse SDF was downloaded from the University of Vienna.
Its ZIP SHA-256 is `e65e7d6a8ca55c55d1d1dda9e8805251714f45859c10998dfabe4752831bc212`.
The SDF SHA-256 and per-record split are in
`validation/results/platinum-diverse-scaffold-split-2026-11-21.json`. The
split excludes scaffolds in the repository's Tier A and Tier B manifests.
It contains 480 development, 1,932 holdout and 447 excluded records. This
checks two known development corpora, so other previous exposure remains
unassessed. The holdout was not scored in this run.

The macOS arm64 run uses the published chematic 1.2.2 wheel and RDKit
2026.09.2, seed 42, up to 20 candidates and at most 200 MMFF94 iterations.
The generated graphs are checked atom by atom and bond by bond. Success means
at least one returned structure has heavy-atom symmetry-aware RMSD at most
1 Å to the crystal reference and retains the specified tetrahedral and
double-bond stereochemistry. The full development denominator includes
failures. Percentile intervals use paired molecule bootstrap samples.

| Stage | chematic | RDKit | Difference | Paired 95% interval |
| --- | ---: | ---: | ---: | ---: |
| ETKDG | 299/480 | 298/480 | +0.21 pp | −2.08 to +2.50 pp |
| ETKDG and MMFF94 | 300/480 | 304/480 | −0.83 pp | −3.13 to +1.46 pp |
| Lowest MMFF94 energy only | 151/480 | 148/480 | +0.63 pp | −1.46 to +2.71 pp |

Both implementations returned 20 candidates for all 480 inputs. No declared
stereochemistry violations were detected among the 9,600 raw and 9,600
minimized candidates per engine. This evaluates stereochemistry assigned in
the input on the recorded host.

The source wheel build and other diagnostics ran concurrently on the host.
Times in the row data are exploratory and do not support a speed comparison.
These development data were inspected to guide changes. A separate equal-time
benchmark, independent holdout, geometry checks and Linux published-wheel
reproduction remain necessary.

The force-field pilot checks identical graph and explicit-H coordinates for
21 other development molecules. All three fields agree on total energy and
gradient at those coordinates, and six selected gradient components per
field per molecule agree with central differences. Its result and conditions
are recorded in `validation/results/force-field-coordinate-parity-pilot-2026-11-21.json`.
RDKit MMFF property construction can alter aromaticity on its input, so each
force field is built from a fresh molecule copy.

RDKit 2026.09.2 on this host raises `RuntimeError: bad any cast` for
`MurckoScaffoldSmiles(smiles="Cc1ccccc1")` even in a fresh environment.
The split therefore uses RDKit 2026.09.1. No scaffold operation occurs in
the 2026.09.2 conformer measurements.
