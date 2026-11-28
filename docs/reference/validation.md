# Validation report

## v1.2.2 published release

The release tag points to `d7962f9501da9ad2f2cdc5eaddf6d665d79d60bc`.
The [channel verification](https://github.com/kent-tokyo/chematic/blob/main/validation/results/release-channel-verification-v1.2.2.json)
records the GitHub Release, 20 crates, npm tarball, 25 PyPI wheels plus one
sdist, docs.rs and project site.

All 26 PyPI hashes match the original build. The
[published-wheel recovery run](https://github.com/kent-tokyo/chematic/actions/runs/44887442455)
passed three runtime and eleven chemistry jobs, with chemistry comparisons
on Linux, Windows, macOS arm64 and macOS x86-64. The gate records a
native-parsed bridgehead amine CIP difference against RDKit 2026.03.6 in the
10,000-row lane.

The [published Python benchmark](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-21-v122-published-comparison.md)
uses exact wheels on one macOS arm64 host. It reports descriptor, fingerprint
and bounded MCS results with their output checks.

## v1.2.3 published release

The tag points to `f63733dec7ae8a63262e9ec0c29cd916b2e7954b`.
The [release workflow](https://github.com/kent-tokyo/chematic/actions/runs/45008595311)
completed successfully on its second attempt after PyPI index propagation.
It verified the hashes of 25 wheels and one source distribution, three runtime
lanes and chemistry comparisons on Linux, Windows and both macOS architectures.
The GitHub Release, npm and 20 crates were also published. The checked-in
[channel record](https://github.com/kent-tokyo/chematic/blob/main/validation/results/release-channel-verification-v1.2.3.json)
still records the pre-publication state.

## v1.2.4 published release

The tag points to `c66a7294aeb82f20c30a9365884e8e16a40bb01a`. The
[release workflow](https://github.com/kent-tokyo/chematic/actions/runs/45158039310)
passed on its second attempt after PyPI index propagation. It verified build
hashes for 25 wheels and one source distribution, runtime checks, and chemistry
comparisons on Linux, Windows and both macOS architectures. The GitHub Release,
npm package and 20 crates were published. The checked-in
[channel record](https://github.com/kent-tokyo/chematic/blob/main/validation/results/release-channel-verification-v1.2.4.json)
retains its pre-publication state and is not evidence of those later checks.

The [checked ensemble guide](checked-conformer-ensembles.md) describes the
opt-in Python methods. The [holdout protocol](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-22-platinum-holdout-protocol.md)
and [macrocycle evaluation](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-11-22-macrocycle-pdb-evaluation.md)
report source-wheel measurements for the identified candidate artifacts;
they do not establish those results for the published wheel.

## v1.2.5 published release

The tag points to `b8c3c56fc83a580e7b18b4d4391c3f6b719ae46b`. The
[release workflow](https://github.com/kent-tokyo/chematic/actions/runs/45724794309)
passed on its second attempt. It ran CI, coverage, package smoke tests,
published-wheel chemistry gates and the npm and crates.io publication jobs.
The checked-in
[channel record](https://github.com/kent-tokyo/chematic/blob/main/validation/results/release-channel-verification-v1.2.5.json)
retains its pre-publication state and is not evidence of those later checks.

The release includes reaction SVG placement and 3D geometry fixes. Their
regression tests cover selected examples; broader 3D quality measurements
remain to be done.

## v1.2.6 release candidate

The source adds Python 3D stereo and square-planar diagnostics, timed and
checked conformer generation, and SMILES and CDXML reader fixes. Their
regression tests cover selected examples. Release artifact qualification and
the published-wheel chemistry gates remain to be run on the tagged source.

## Retained package and source evidence

| Area | Record |
|---|---|
| Published chemistry | [v1.2.0 runtime](https://github.com/kent-tokyo/chematic/actions/runs/39511354658) and [chemistry](https://github.com/kent-tokyo/chematic/actions/runs/39511354905); [v1.1.3 four-wheel row packet](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-14-v113-published-artifact-packet.md) |
| Published timings | [v1.1.8 Labute ASA](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-19-v118-published-labute-asa-cosmolkit.md), [v1.0.42 Linux and macOS](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-11-v1042-published-rdkit-packet.md), and [v1.0.30 browser Morgan](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-03-v1030-published-linux-three-browser-paired20.md) |
| Reactions and CTfile | [Reaction and workflow fixtures](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-12-rdkit-workflows-expanded-source.md) and [CTfile residuals](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-18-trust-gates-and-residuals.md) |
| 3D | [A6 topology bounds](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-12-a6-rdkit-2026-09-bounds-followup.md) and [torsion search](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-24-etkdg-torsion-search-rdkit.md) |
| Independent evaluation | [v1.1.8 holdout score](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-10-19-trust-holdout-v1-score.md) and [row data](https://github.com/kent-tokyo/chematic/blob/main/validation/results/2026-10-19-trust-holdout-v1-v118.json) |
| Directory layout | [Validation directory guide](https://github.com/kent-tokyo/chematic/blob/main/validation/README.md): where inputs, expected values and recorded results live |
| Full inventory | [Benchmark index](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md), [historical inventory](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/archive-index.md), and [versioned results](https://github.com/kent-tokyo/chematic/tree/main/validation/results/) |

## How to interpret results

Use the artifact or source revision, comparator, corpus, options, host and
failure policy of each measurement. Refusals and abstentions stay in their
denominators. A speed comparison requires equivalent work and outputs.
Rust CI coverage applies to its declared scope. Chemistry accuracy is
measured by the operation-specific tests. See the
[benchmark method](benchmark.md) and [compatibility scope](compatibility-scope.md).

## Local consistency checks

```bash
python3 scripts/checks/check_release_docs_consistency.py
python3 scripts/checks/check_benchmark_index.py
python3 scripts/checks/check_roadmap_disposition.py
python3 scripts/checks/check_product_names.py
```

See [compatibility scope](compatibility-scope.md), [benchmark method](benchmark.md),
and the [accuracy plan](../plans/rdkit-accuracy-plan.md) for the support boundary.
