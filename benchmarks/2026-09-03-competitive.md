# Competitive benchmark record — 2026-09-03

Status: complete for schematic and RDKit; Open Babel was not installed and is
recorded as `not_installed`, not as a zero or a pass.

This is an operation-specific measurement, not a claim of universal
superiority. The run used schematic 1.0.1 from commit `bfaf1d52f9094702eda011451d6479ff0265d3cf`,
RDKit 2025.09.3, Python 3.13.6, macOS 26.5.2 arm64, and the repository's
current benchmark commands. The stale 0.89.0 run in the older result directory
was rejected by the version preflight and is not included here.

## Results

| Operation | schematic | RDKit | Relative result |
|---|---:|---:|---:|
| Import | 13.3 ms | 117.5 ms | schematic 8.8× faster |
| SMILES parse, 5,000 | 2.17 µs/mol | 48.65 µs/mol | schematic 22.3× faster |
| Canonical SMILES, 5,000 | 101.61 µs/mol | 24.97 µs/mol | schematic 4.1× slower |
| ECFP4, 5,000 diverse | 29.22 µs/mol | 88.77 µs/mol | schematic 3.0× faster |
| SDF read, 365 | 588.38 µs/mol | 99.15 µs/mol | schematic 5.9× slower |
| SDF write, 365 | 334.5 µs/mol | 39.15 µs/mol | schematic 8.5× slower |

The descriptor-agreement lane also completed on 4,999 molecules. Its detailed
JSON is in the result directory; array-family agreement remains intentionally
reported as compatibility data rather than a single accuracy score.

## Reproduction

```bash
python3 scripts/validate_competitive_benchmark_manifest.py
python3 scripts/run_competitive_benchmark.py \
  --result-dir validation/results/competitive-benchmark-v1.0.1-2026-09-03
```

The run is resumable with `--resume`. Raw state and per-operation logs are
preserved at
`validation/results/competitive-benchmark-v1.0.1-2026-09-03/`. Corpus hashes:

- `scripts/descriptor_census_corpus.smi`:
  `d6f2ba3f128296f935007f0b0813aa97b6ebcc2457e014ddca2213ddd655276c`
- `scripts/chembl_accuracy_corpus_4999.smi`:
  `1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`

