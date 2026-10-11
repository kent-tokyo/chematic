# v1.1.1 eager-preparation regression follow-up

This record compares the published CheMatic v1.1.1 Python wheel with a
source-built candidate that restores lazy aromatic and ring preparation.
It answers a narrow before/after question on one host; it is not evidence for
a future published package or for other platforms.

## Provenance

- Baseline: published `chematic==1.1.1` macOS arm64 CPython 3.13 wheel
- Baseline wheel SHA-256:
  `1589f9f7a88b6271c8fb5631b4e8b92583c890f66b03f6a6793d40aa8998779b`
- Candidate: branch `codex/v1.1.1-evaluation-fixes`, built before the final
  documentation-only changes in this record
- Candidate wheel SHA-256:
  `71d1072c52e8847e0e5ecfef82c00b4d5e626179eb07571a23b61d6a7affe9d7`
- Host: macOS 27.0.1, Apple arm64, Python 3.13.6
- Corpus: first 1,000 rows of `chembl_accuracy_corpus_4999.smi`
- Corpus SHA-256:
  `1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`
- Protocol: 20 alternating AB/BA blocks, fresh process per arm, one untimed
  warm-up call, bootstrap median interval with fixed seed 0

The ratio below is baseline time divided by candidate time. Values above one
favor the candidate. Every lane produced the same output digest in both arms.

## Results

| Operation boundary | Baseline/candidate median | Paired 95% interval | Interpretation |
|---|---:|---:|---|
| SMILES parse only | 2.144x | 2.058-2.229 | Candidate faster |
| Parse plus TPSA | 1.393x | 1.257-1.506 | Candidate faster |
| Parse plus Morgan R2/2048 | 1.186x | 1.159-1.211 | Candidate faster |
| TPSA after precomputation | 1.003x | 0.981-1.117 | Inconclusive |
| Morgan after precomputation | 0.980x | 0.929-0.997 | Candidate about 2% slower |

The change removes eager descriptor/ring preparation from `from_smiles` while
retaining mutation-invalidated caches when callers actually request those
properties. It therefore improves parse-only and one-property pipelines. It
does not improve the already-precomputed Morgan operation, and this record
does not claim otherwise.

## Published v1.1.1 artifact gate

Hosted run
[`38184048499`](https://github.com/kent-tokyo/chematic/actions/runs/38184048499)
reproduced the pinned public-wheel chemistry counts on Windows amd64, macOS
arm64 and macOS x86-64. The Linux x86-64 sdist and the Linux, macOS and Windows
xsmarts-autoconf lanes also passed. The gate covers the 83 reaction fixtures,
310,000 SMARTS cells, CIP, hybridization, MMFF94 typing, MOL stereo round-trip,
the 4,072 unusual-SMILES descriptor rows, the pinned BioTransformer corpus and
the 265-row A6 scorecard. These are v1.1.1 artifact results; the HBD/QED and
lazy-parse fixes above remain source-only until a later release.

## Raw records

- [`parse`](2026-10-12-v111-public-vs-lazy-fix-parse-paired20.json)
- [`TPSA pipeline`](2026-10-12-v111-public-vs-lazy-fix-tpsa-pipeline-paired20.json)
- [`Morgan pipeline`](2026-10-12-v111-public-vs-lazy-fix-morgan-pipeline-paired20.json)
- [`TPSA precomputed`](2026-10-12-v111-public-vs-lazy-fix-tpsa-precomputed-paired20.json)
- [`Morgan precomputed`](2026-10-12-v111-public-vs-lazy-fix-morgan-precomputed-paired20.json)
