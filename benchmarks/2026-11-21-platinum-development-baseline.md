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
These development data were inspected to guide changes. A separate equal-time
benchmark, independent holdout, geometry checks and Linux published-wheel
reproduction remain necessary.

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

RDKit 2026.09.2 on this host raises `RuntimeError: bad any cast` for
`MurckoScaffoldSmiles(smiles="Cc1ccccc1")` even in a fresh environment.
The split therefore uses RDKit 2026.09.1. No scaffold operation occurs in
the 2026.09.2 conformer measurements.
