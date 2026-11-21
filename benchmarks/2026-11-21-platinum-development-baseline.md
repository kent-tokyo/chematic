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
An additional audit against the four MMFF94 post-freeze manifests and both
3D trust-holdout files leaves 404 development, 1,606 holdout and 849 excluded
records. The eight invalid exclusion inputs are unchanged.

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
Molecules with a shared scaffold may be correlated; the reported intervals
do not model scaffold-level clustering.

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

The 444 development inputs above were also checked at one shared ETKDG
coordinate set per molecule. This set was already filtered for RDKit MMFF
support during the preceding generation run. All 444 were prepared by both
engines for MMFF94, MMFF94s and UFF. The 95th percentile of absolute energy
differences was at most `1.99e-13 kcal/mol`, and the 95th percentile over all
gradient-component differences was at most `1.53e-13 kcal/mol/Å`. Maximum
gradient-component differences were `8.85e-6`, `8.87e-6` and `4.30e-4`
respectively. Six sampled coordinates per field per molecule were also
checked by central differences at `1e-5 Å`. Maximum differences were 0.011,
0.055 and `4.31e-7 kcal/mol/Å`, respectively. Both larger MMFF residuals
occurred on `07U_3TXO_A` and warrant further numerical analysis.

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
An opt-in torsion-gradient evaluator now uses the same unrounded plane norms
as the torsion energy, and its torsion-only gradient was checked against
finite differences for near-linear geometries. A second opt-in evaluator also
uses unrounded normals for out-of-plane terms, with term-specific finite
difference checks. Other MMFF terms retain their RDKit-compatible gradients,
and the optimizer has not yet been migrated to the alternate derivatives.
Full near-singular gradient and minimizer validation remain open.

The 444-row development rerun with the retry enabled changed four conformers
in two molecules. Three were the near-linear amine conformers above, and one
was in `NNL_3P0Q_A`. The independent angle screen flagged 33/8,880 returned
chematic conformers, compared with 37/8,880 before the retry. The number of
inputs with a screened conformer within 1 Å remained 271/444 for chematic and
273/444 for RDKit. The same fixed-candidate comparison had concurrent local
workload, so its recorded elapsed times are not used for performance claims.

With the optional 400-iteration retry included in the 0.2-second per-engine
budget, the same development protocol recovered 222/444 for chematic and
192/444 for RDKit. The paired 95% bootstrap interval for the difference was
+4.05 to +9.68 percentage points. Requiring a converged matching conformer
gave 172 and 152 recoveries, with an interval of +2.48 to +6.76 points.
The results changed across the two local runs for both engines, so repeated
measurements on additional hosts remain necessary.
Repeating the repair-enabled protocol produced 224/444 and 199/444 recoveries
(paired 95% interval +3.15 to +8.11 points). Requiring convergence gave
172/444 and 157/444 (interval +1.35 to +5.41 points).

The [hosted Linux source-wheel run](https://github.com/kent-tokyo/chematic/actions/runs/44912969456)
built the pull-request merge commit `e64c00a710c859627aeb3e425158dc4a1d80ce43`
from head `16c19ef4c9b625cd397e1c184a9174941931e541` and base
`0be29688d88c13baebc323418a5f4143cd31836d`. Its wheel SHA-256 was
`d986e7c59fe810023b30002b2fccb132feac0881f030c0067109c152d1b28ad8`.
On 444 development molecules under the geometry and stereo screens, each
engine recovered 272 within 1 Å (paired 95% interval −1.80 to +1.80 points).
The 0.2-second, 100-molecule development sample recovered 46 for chematic
and 37 for RDKit (interval +4 to +15 points). Requiring convergence gave
37 and 33 (interval 0 to +9 points). The hosted source wheel and the smaller
time-budget sample are separate from published-artifact and holdout evidence.

RDKit 2026.09.2 on this host raises `RuntimeError: bad any cast` for
`MurckoScaffoldSmiles(smiles="Cc1ccccc1")` even in a fresh environment.
The split therefore uses RDKit 2026.09.1. No scaffold operation occurs in
the 2026.09.2 conformer measurements.

A conservative graph audit of the eight unsanitizable exclusions is recorded
in `validation/results/platinum-invalid-exclusion-graph-audit-2026-11-22.json`.
Each eligible scaffold was searched as an atom-number subgraph of the
unsanitized exclusion structures with every bond order relaxed. No candidate
scaffold matched. This treats any potential match as an exclusion and does not
assert that the invalid inputs have chemically meaningful sanitized structures.
The audit covers the recorded exclusion datasets; it cannot establish that no
other historical data or methodological choices exposed a holdout molecule.
