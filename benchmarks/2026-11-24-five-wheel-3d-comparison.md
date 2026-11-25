# Five published wheels: 100-molecule 3D development comparison

This is an exploratory comparison of five published macOS arm64 CPython 3.13 wheels, on 100 name-hash-selected records from the already-scored v3 COD cohort. It is not an independent holdout, does not establish statistical noninferiority, and cannot directly measure improvement against the corrected candidate path. The three chematic versions used their common `rdkit_embed_multiple_confs_mmff` compatibility API, not the experimental corrected 1.2.4 API. All installed package files were checked against their wheel archives before the run. Every worker reported the expected installed version; chematic and RDKit worker extension hashes were recorded per row. Scorers, wheel hashes, input and selection hashes, and candidate seed schedule are recorded in the raw JSONL and the freeze manifest.

Each engine prepared the same canonical isomeric SMILES and explicit-H graph, in a separate spawned worker. Atom and indexed bond order were checked before using candidate coordinates. Engines were run serially, rotating the starting order by input. Each engine made up to 20 independent single-conformer attempts at seeds 42, 84, ... 840, with 200 MMFF iterations and a 10-second generation deadline. Startup, parsing and H addition were excluded from generation timing; observed elapsed includes worker termination/cleanup. This repeated-single-conformer protocol is not equivalent to native multiple-conformer embedding. RDKit and COS can return unoptimized candidates: only status 0 is counted as converged. A timeout or an engine error is retained in the denominator; fixed-count recovery requires all 20 attempts completed. RMSD <= 1 Å to the reference plus declared-stereo and supported bond/angle/clash checks defines recovery. Stereo violations and unevaluable outcomes fail the check. A separately reported 0.2-second prefix accepts candidates delivered by that time during the same 10-second run; it does **not** terminate the worker at 0.2 seconds.

| Engine | Fixed recovered / 100 | Fixed converged recovered | 0.2 s prefix recovered | Status: ok / error / timeout | Median elapsed s | Median peak worker RSS MiB |
|---|---:|---:|---:|---|---:|---:|
| chematic 1.2.2 | 84 | 51 | 66 | 94 / 5 / 1 | 0.451 | 88.5 |
| chematic 1.2.3 | 84 | 51 | 66 | 94 / 5 / 1 | 0.453 | 88.4 |
| chematic 1.2.4 | 84 | 51 | 67 | 94 / 5 / 1 | 0.449 | 88.5 |
| RDKit 2026.09.2 | 84 | 51 | 62 | 99 / 0 / 1 | 0.759 | 67.7 |
| COSMolKit 0.6.3 | 91 | 55 | 62 | 99 / 1 / 0 | 0.590 | 90.0 |

Each chematic release had five boron MMFF type-assignment errors, on the same five structures. COS had one error when an MMFF force field was requested without a generated conformer. There were no input or scoring errors, missing RMSD scores, or missing stereo checks. All recorded candidates had zero strict stereo violations and unevaluable stereo in this sample. A failed or timed-out worker may lack a measured RSS, and the measured high-water mark is sampled after completed attempts; startup/import memory is included. Its count and maxima are in the summary. Wall-clock medians include incomplete and failed runs and are not pure throughput estimates. Runs used a shared single machine, and 100 exposed molecules are too few to infer broad quality or speed rankings.

A first experimental run was interrupted after one completed row when a missing geometry-support metadata field was discovered. The code and freeze were corrected before this complete run; all affected rows are treated as exposed development data. An earlier 20-molecule probe used a modified Python initializer for the 1.2.4 wheel; this 100-molecule run used file-verified unmodified packages.

The empirical 95th percentile (nearest-rank, all 100 rows with elapsed time) and maximum, in seconds, are:

| Engine | p95 | Maximum |
|---|---:|---:|
| chematic 1.2.2 | 2.145 | 10.004 |
| chematic 1.2.3 | 2.141 | 10.004 |
| chematic 1.2.4 | 2.138 | 10.004 |
| RDKit 2026.09.2 | 2.696 | 10.003 |
| COSMolKit 0.6.3 | 2.871 | 8.334 |

These are elapsed times until completion, error or timeout, including cleanup and excluding worker startup. As in the median column, they include failed rows; a timeout censors the underlying runtime. They should not be treated as p95 time to produce 20 valid conformers.
