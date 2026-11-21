# Platinum Diverse development baseline

The 2017_01 Platinum Diverse SDF was downloaded from the University of Vienna.
Its ZIP SHA-256 is `e65e7d6a8ca55c55d1d1dda9e8805251714f45859c10998dfabe4752831bc212`.
The SDF SHA-256 and per-record split are in
`validation/results/platinum-diverse-scaffold-split-2026-11-21.json`. The
split excludes scaffolds in the repository's Tier A and Tier B manifests.
It contains 480 development, 1,932 holdout and 447 excluded records. This
checks two known development corpora, so other previous exposure remains
unassessed. The holdout was not scored in this run.

A stricter split additionally excludes the 4,999-SMILES ChEMBL accuracy
corpus, which includes inputs used in earlier conformer experiments. It has
444 development, 1,780 holdout and 635 excluded records. The 480-row
measurements below remain development results under the earlier split.
A broader split also excludes the NCI 5k and browser 10k scaffolds, leaving
404 development, 1,610 holdout and 845 excluded records. Eight exclusion
SMILES fail RDKit sanitization and are listed in the split file; their
scaffolds could not be compared. Other previously used datasets also require
an exposure audit before the holdout is treated as independent.

The macOS arm64 run uses the published chematic 1.2.2 wheel and RDKit
2026.09.2, seed 42, up to 20 candidates and at most 200 MMFF94 iterations.
Their wheel SHA-256 values are respectively
`682c96b2c74b9f43c21efb11a8d8979e9760f7856d911af7c455238560464266` and
`3860ba53505f27cb24d05383027f78e3f82aefd6caf193c217a94d876d91b023`.
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
These development data were inspected to guide changes. Repeated equal-time
measurements, independent holdout and Linux published-wheel reproduction remain
necessary.

The commands below use the SDF from the supplied ZIP and the recorded split.
The baseline script is available at commit `7b912c7`.

```sh
python scripts/benchmarks/platinum_development_probe.py \
  --sdf /path/to/platinum_diverse_dataset_2017_01.sdf \
  --split validation/results/platinum-diverse-scaffold-split-2026-11-21.json \
  --limit 480 --output /tmp/platinum-development.json
python scripts/benchmarks/summarize_platinum_probe.py \
  /tmp/platinum-development.json --output /tmp/platinum-summary.json
```

The force-field pilot checks identical graph and explicit-H coordinates for
21 other development molecules. All three fields agree on total energy and
gradient at those coordinates, and six selected gradient components per
field per molecule agree with central differences. Its result and conditions
are recorded in `validation/results/force-field-coordinate-parity-pilot-2026-11-21.json`.
RDKit MMFF property construction can alter aromaticity on its input, so each
force field is built from a fresh molecule copy.

The candidate wheel exposes the existing corrected sixth-order ETKDG gradient
through Python for one requested ensemble. Its SHA-256 is
`8e84189a3100aa6d7ca711d90eec2ce34f7977a2567defc91b4e498dcf34b9ef`.
On the same 480 development inputs and 20-candidate budget, it recovered
301/480 experimental structures within 1 Å versus RDKit's 304/480. The paired
95% bootstrap interval for the difference is −2.50 to +1.25 percentage points.
The minimum-energy candidate recovered 146/480 versus 148/480. All 9,600
returned candidate conformers retained declared stereochemistry. This candidate
run used the same host and is development evidence.

An equal-time development probe used the 100 lowest SHA-256 record names from
that development set, 0.2 seconds per engine and at most 20 returned
conformers. Each attempt embeds and minimizes one candidate, and only attempts
completed by the deadline count. Engine order alternates between molecules.
Stereochemistry and RMSD scoring occur after timing. The candidate recovered
52/100 structures and RDKit 50/100, with a paired 95% bootstrap difference
interval of −2 to +6 percentage points. Median completed candidates were 7.5
and 4 respectively, and no specified-stereo violations were detected. This
small development sample supports further investigation under the same budget.

The stricter 444-row development split was measured again with a geometry
screen. It uses UFF reference bond lengths with 50% relative tolerance, UFF
reference angles with 35° tolerance, 60°/90° reference angles for three- and
four-membered rings, and a heavy-atom nonbonded cutoff of 0.65 times the sum
of van der Waals radii for pairs separated by at least three graph bonds.
These are screening thresholds. All reference bonds and angles were available
on these inputs. Among 8,880 MMFF-refined conformers per engine, the screen
flagged 37 angle cases in each engine and no bond-length or heavy-atom clash
cases. With the screen and the 1 Å RMSD and stereo conditions applied together,
the candidate recovered 271/444 and RDKit 273/444, with a paired 95% interval
of −2.48 to +1.58 percentage points.

The equal-time rerun uses the stricter split and the wheel with the final
gradient check, SHA-256
`2f0f74d80e95ae9708241825b44470222de360898439a76d5a9e507351af14d1`.
At 0.2 seconds per engine, the 444 inputs gave 227 recoveries for chematic and
202 for RDKit. The paired 95% interval for the difference was +3.15 to
+8.33 percentage points. Median completed candidates were 9 and 5. The run
uses individual embed/minimize calls and alternates engine order per input.
This is development evidence for that specific time budget and host.
Recovery accepts conformers meeting the RMSD, stereo and geometry screens,
including those with a nonzero optimizer status. The row data report final
gradient-checked optimizer statuses separately for both engines.
Requiring the matching conformer to report convergence under the common
gradient threshold reduced recoveries to 174/444 and 158/444 respectively,
with a paired 95% interval for the difference of +1.58 to +5.63 points.

On one input (`07U_3TXO_A`), the prior candidate returned three converged
statuses at nearly linear C–N–C geometries. At those returned coordinates,
independent RDKit evaluation also produced extremely large MMFF gradients.
The checked API now reports a nonzero status when the final gradient is
nonfinite or exceeds 0.01 kcal/mol/Å in any component. The term-level
diagnostic records finite differences near the singular geometry for further
investigation. These checks report the optimizer outcome more conservatively.
For conformer 9, RDKit's torsion gradient at component 16 is approximately
−1.20 million kcal/mol/Å while central differences at steps `1e-8` and
`1e-10` Å are approximately +42,735 and +42,769. This is a near-singular
example and needs a dedicated derivative fix and further numerical analysis.
The checked API leaves the returned coordinates in place and marks the
conformer unfinished when its final gradient exceeds the threshold.

An opt-in retry moves a near-linear three-coordinate nitrogen by 0.02 Å
perpendicular to the near-linear bond and runs a bounded additional MMFF94
minimization. The retry is accepted only if it converges, does not increase
energy, passes the 0.01 kcal/mol/Å final-gradient check, and retains declared
stereochemistry. With 400 additional iterations allowed on `07U_3TXO_A`,
the three checked unfinished conformers reached status 0 and energies of
36.402, 38.199 and 54.036 kcal/mol. Their earlier energies were 172.855,
162.082 and 163.252 kcal/mol. The repair wheel SHA-256 is
`aa69dafd9975d6bd05d263e694c907a863319b09280b9d36c729fdd2569b3053`.
The near-singular MMFF derivative remains an open numerical issue.

RDKit 2026.09.2 on this host raises `RuntimeError: bad any cast` for
`MurckoScaffoldSmiles(smiles="Cc1ccccc1")` even in a fresh environment.
The split therefore uses RDKit 2026.09.1. No scaffold operation occurs in
the 2026.09.2 conformer measurements.
