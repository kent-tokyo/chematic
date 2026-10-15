# Python unambiguous organometallic cleanup (2026-10-14)

This candidate is not adopted: the complete matrix found no strict baseline-relative speed win. Its isolated source, outputs and timing records are retained as diagnostic evidence.

RDKit-compatible model construction ranks all atoms to select a metal bond during organometallic cleanup. When every affected nonmetal donor has one eligible single bond to a metal, that selection is unique. This change skips ranking in that case. It preserves pre-cleanup ring perception and its errors, cache state, hydrogen handling, and the ranked fallback for ambiguous or parallel bond selection. Nothing moves into parsing.

Both 1.1.4 wheels are fresh local source builds from the frozen attached worktree. The baseline includes the prior unmerged relevant-cycle and V3000 writer candidates; the candidate adds only this cleanup change, its tests, and a changelog entry. The last verified base is `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared Git metadata became unresponsive; this run does not establish current main, commit, release, or registry-artifact status. The source hashes and self-contained source archive identify what was tested.

## Correctness

- Frozen rank-based cleanup versus the new path: 64 metal-neighbour layouts with shared centres, existing dative bonds, and both ring-cache states; no partial rewrite before ranked fallback; explicit hydrogen/isotope/charge/stereo and invalid inputs; internally constructed parallel bonds; all 10,000 ordinary inputs. Full graph, properties, caches, errors, and subsequent sanitization are compared.
- Rust SMILES library: 391 passed, 1 ignored. Python suite: 1,389 passed, 1 skipped. Clippy for the changed crate and its targets, workspace formatting, and the refactor budget pass.
- All 14,072 inputs (10,000 ordinary plus 4,072 stress): first/repeated descriptor and four Morgan output streams are byte-identical to baseline, including retained errors. Equality does not convert an error into a valid chemistry result.
- Test-data omissions initially stopped Python collection and the benchmark before any completed blocks. Data was copied unchanged and complete reruns passed. The intermediate pre-parallel-guard timing was stopped, retained, and excluded from final counts.

## Protocol and scope

CPython 3.13.6, Rust 1.97.0, maturin 1.13.3, macOS arm64, release defaults. RDKit 2026.09.1 and COSMolKit 0.5.0rc22. Unfiltered 1,000-input slices at offsets 0 and 5,000 of `rdkit-js-browser-10k-v1.smi`; SHA-256 `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`. These are exposed development inputs, not a sealed holdout. Existing corpus provenance and redistribution terms remain unchanged.

Each slice has 20 rotating four-arm blocks in fresh processes. First use times an operation on a fresh parsed object; repeated use follows one untimed call; parse-inclusive times actual parse plus operation. Native results are retained and normalization happens outside timing. All rows and blocks are retained; every worker checks corpus, extension hashes, and deterministic outputs. Whole-process peak RSS includes imports and all operations. The two final slices run sequentially after builds/tests, without concurrent build or compression.

Tables give comparator time / candidate time, paired-median 95% bootstrap interval (10,000 deterministic resamples), and faster blocks. Above 1 means candidate faster. A strict win requires the interval lower bound above 1 and 20/20 faster blocks. COSMolKit comparisons are unranked if its outputs or candidate outputs disagree with RDKit; no mismatching row is removed. One host and seven operations cannot establish whole-toolkit superiority.

## Initial descriptor calculation

| Offset | Operation / mode | Fresh baseline / candidate | RDKit / candidate | COSMolKit / candidate |
|---:|---|---|---|---|
| 0 | tpsa/first_use | 1.000 [0.990, 1.021]; 10/20 | 1.470 [1.443, 1.505]; 20/20 | 0.512 [0.505, 0.525]; 0/20 |
| 0 | tpsa/parse_inclusive | 0.999 [0.988, 1.018]; 10/20 | 13.149 [13.004, 13.302]; 20/20 | 12.258 [12.133, 12.318]; 20/20 |
| 0 | labute_asa/first_use | 0.974 [0.952, 1.004]; 6/20 | 0.690 [0.672, 0.718]; 0/20 | 0.341 [0.334, 0.364]; 0/20 |
| 0 | labute_asa/parse_inclusive | 1.013 [1.003, 1.032]; 14/20 | 12.345 [12.240, 12.413]; 20/20 | 11.687 [11.514, 11.795]; 20/20 |
| 5000 | tpsa/first_use | 1.005 [0.999, 1.019]; 14/20 | 0.344 [0.338, 0.349]; 0/20 | 0.097 [0.093, 0.102]; 0/20 |
| 5000 | tpsa/parse_inclusive | 1.001 [0.987, 1.017]; 10/20 | 4.272 [4.210, 4.319]; 20/20 | 4.404 [4.349, 4.464]; 20/20 |
| 5000 | labute_asa/first_use | 1.003 [0.993, 1.035]; 11/20 | 0.166 [0.161, 0.170]; 0/20 | 0.067 [0.065, 0.074]; 0/20 |
| 5000 | labute_asa/parse_inclusive | 0.998 [0.990, 1.011]; 9/20 | 3.998 [3.936, 4.039]; 20/20 | 4.301 [4.192, 4.375]; 20/20 |

## Strict comparator wins

| Mode | Demonstrated / eligible |
|---|---:|
| first_use | 13/27 |
| repeated_use | 22/27 |
| parse_inclusive | 27/27 |

## Complete matrix

| Offset | Operation / mode | Fresh baseline / candidate | RDKit / candidate | COSMolKit / candidate |
|---:|---|---|---|---|
| 0 | tpsa/first_use | 1.000 [0.990, 1.021]; 10/20 | 1.470 [1.443, 1.505]; 20/20 | 0.512 [0.505, 0.525]; 0/20 |
| 0 | tpsa/repeated_use | 1.009 [1.002, 1.048]; 15/20 | 9.850 [9.595, 10.119]; 20/20 | 1.724 [1.671, 1.784]; 19/20 |
| 0 | tpsa/parse_inclusive | 0.999 [0.988, 1.018]; 10/20 | 13.149 [13.004, 13.302]; 20/20 | 12.258 [12.133, 12.318]; 20/20 |
| 0 | labute_asa/first_use | 0.974 [0.952, 1.004]; 6/20 | 0.690 [0.672, 0.718]; 0/20 | 0.341 [0.334, 0.364]; 0/20 |
| 0 | labute_asa/repeated_use | 0.986 [0.974, 0.999]; 6/20 | 2.732 [2.562, 2.836]; 20/20 | 1.474 [1.321, 1.544]; 20/20 |
| 0 | labute_asa/parse_inclusive | 1.013 [1.003, 1.032]; 14/20 | 12.345 [12.240, 12.413]; 20/20 | 11.687 [11.514, 11.795]; 20/20 |
| 0 | num_rings/first_use | 1.009 [0.991, 1.027]; 11/20 | 0.035 [0.033, 0.036]; 0/20 | unranked |
| 0 | num_rings/repeated_use | 1.006 [0.978, 1.035]; 10/20 | 1.737 [1.644, 1.800]; 19/20 | unranked |
| 0 | num_rings/parse_inclusive | 1.012 [0.982, 1.026]; 12/20 | 11.232 [10.971, 11.400]; 20/20 | unranked |
| 0 | morgan2_2048/first_use | 0.999 [0.986, 1.015]; 10/20 | 3.089 [3.051, 3.120]; 20/20 | 2.022 [1.983, 2.046]; 20/20 |
| 0 | morgan2_2048/repeated_use | 1.005 [0.998, 1.019]; 14/20 | 5.201 [5.121, 5.249]; 20/20 | 3.253 [3.246, 3.293]; 20/20 |
| 0 | morgan2_2048/parse_inclusive | 1.006 [1.000, 1.018]; 14/20 | 8.573 [8.499, 8.705]; 20/20 | 7.584 [7.473, 7.658]; 20/20 |
| 0 | morgan2_chiral/first_use | 0.997 [0.988, 1.021]; 9/20 | 1.239 [1.234, 1.272]; 20/20 | 0.830 [0.814, 0.837]; 0/20 |
| 0 | morgan2_chiral/repeated_use | 1.008 [0.999, 1.016]; 13/20 | 5.055 [5.011, 5.088]; 20/20 | 3.169 [3.116, 3.238]; 20/20 |
| 0 | morgan2_chiral/parse_inclusive | 1.002 [0.992, 1.011]; 10/20 | 4.446 [4.417, 4.490]; 20/20 | 3.905 [3.872, 3.941]; 20/20 |
| 0 | morgan3_2048/first_use | 1.004 [0.993, 1.013]; 11/20 | 3.370 [3.333, 3.392]; 20/20 | 2.251 [2.235, 2.279]; 20/20 |
| 0 | morgan3_2048/repeated_use | 1.007 [0.995, 1.011]; 13/20 | 4.910 [4.871, 4.934]; 20/20 | 3.176 [3.150, 3.217]; 20/20 |
| 0 | morgan3_2048/parse_inclusive | 1.006 [0.997, 1.009]; 13/20 | 8.025 [7.951, 8.075]; 20/20 | 7.031 [6.959, 7.063]; 20/20 |
| 0 | maccs/first_use | 1.003 [1.001, 1.006]; 15/20 | 1.886 [1.882, 1.892]; 20/20 | 2.171 [2.161, 2.174]; 20/20 |
| 0 | maccs/repeated_use | 1.004 [0.999, 1.010]; 13/20 | 1.954 [1.943, 1.964]; 20/20 | 2.250 [2.244, 2.256]; 20/20 |
| 0 | maccs/parse_inclusive | 1.003 [1.001, 1.008]; 15/20 | 2.114 [2.108, 2.123]; 20/20 | 2.387 [2.381, 2.395]; 20/20 |
| 5000 | tpsa/first_use | 1.005 [0.999, 1.019]; 14/20 | 0.344 [0.338, 0.349]; 0/20 | 0.097 [0.093, 0.102]; 0/20 |
| 5000 | tpsa/repeated_use | 1.017 [0.981, 1.030]; 12/20 | 9.870 [9.243, 10.248]; 20/20 | 1.296 [1.241, 1.317]; 19/20 |
| 5000 | tpsa/parse_inclusive | 1.001 [0.987, 1.017]; 10/20 | 4.272 [4.210, 4.319]; 20/20 | 4.404 [4.349, 4.464]; 20/20 |
| 5000 | labute_asa/first_use | 1.003 [0.993, 1.035]; 11/20 | 0.166 [0.161, 0.170]; 0/20 | 0.067 [0.065, 0.074]; 0/20 |
| 5000 | labute_asa/repeated_use | 1.012 [0.985, 1.025]; 12/20 | 2.642 [2.574, 2.717]; 20/20 | 1.199 [1.127, 1.325]; 18/20 |
| 5000 | labute_asa/parse_inclusive | 0.998 [0.990, 1.011]; 9/20 | 3.998 [3.936, 4.039]; 20/20 | 4.301 [4.192, 4.375]; 20/20 |
| 5000 | num_rings/first_use | 0.996 [0.955, 1.015]; 10/20 | 0.092 [0.089, 0.095]; 0/20 | 0.074 [0.068, 0.086]; 0/20 |
| 5000 | num_rings/repeated_use | 0.998 [0.983, 1.017]; 10/20 | 1.662 [1.596, 1.726]; 20/20 | 1.039 [1.017, 1.075]; 16/20 |
| 5000 | num_rings/parse_inclusive | 1.003 [0.988, 1.035]; 12/20 | 11.531 [11.340, 11.634]; 20/20 | 12.611 [12.426, 12.851]; 20/20 |
| 5000 | morgan2_2048/first_use | 0.993 [0.987, 1.005]; 7/20 | 1.389 [1.365, 1.403]; 20/20 | 0.892 [0.873, 0.907]; 0/20 |
| 5000 | morgan2_2048/repeated_use | 0.992 [0.977, 1.000]; 7/20 | 4.674 [4.594, 4.689]; 20/20 | 2.873 [2.809, 2.931]; 20/20 |
| 5000 | morgan2_2048/parse_inclusive | 1.001 [0.988, 1.027]; 11/20 | 4.235 [4.202, 4.276]; 20/20 | 4.099 [4.041, 4.124]; 20/20 |
| 5000 | morgan2_chiral/first_use | 1.001 [0.993, 1.017]; 10/20 | 1.414 [1.393, 1.438]; 20/20 | 0.909 [0.900, 0.920]; 1/20 |
| 5000 | morgan2_chiral/repeated_use | 1.000 [0.992, 1.021]; 10/20 | 4.721 [4.630, 4.819]; 20/20 | 2.907 [2.868, 2.933]; 20/20 |
| 5000 | morgan2_chiral/parse_inclusive | 1.003 [0.985, 1.012]; 11/20 | 4.261 [4.202, 4.293]; 20/20 | 4.101 [4.052, 4.124]; 20/20 |
| 5000 | morgan3_2048/first_use | 1.005 [0.989, 1.018]; 12/20 | 1.609 [1.594, 1.627]; 20/20 | 1.062 [1.052, 1.077]; 19/20 |
| 5000 | morgan3_2048/repeated_use | 0.997 [0.988, 1.001]; 8/20 | 4.542 [4.510, 4.625]; 20/20 | 2.915 [2.869, 2.954]; 20/20 |
| 5000 | morgan3_2048/parse_inclusive | 1.009 [0.987, 1.016]; 11/20 | 4.197 [4.107, 4.247]; 20/20 | 3.992 [3.949, 4.025]; 20/20 |
| 5000 | maccs/first_use | 1.002 [1.000, 1.006]; 14/20 | 1.624 [1.621, 1.629]; 20/20 | 1.742 [1.737, 1.746]; 20/20 |
| 5000 | maccs/repeated_use | 1.004 [0.998, 1.011]; 13/20 | 1.728 [1.720, 1.739]; 20/20 | 1.856 [1.847, 1.868]; 20/20 |
| 5000 | maccs/parse_inclusive | 1.005 [0.996, 1.007]; 13/20 | 1.799 [1.786, 1.806]; 20/20 | 1.928 [1.913, 1.942]; 20/20 |

## Regressions and remaining work

All baseline intervals wholly below 1:
- Offset 0, `labute_asa/repeated_use`: 0.986 [0.974, 0.999]; 6/20

The preceding [root-pruning report](2026-10-14-python-cycle-root-pruning.md) retained a chiral-Morgan parse-inclusive regression against its own baseline. This run uses that candidate inside its fresh baseline and therefore cannot establish that the earlier regression against main is resolved. Eager/lazy ring perception and first-use descriptor/fingerprint gaps remain visible above. Existing bounded ring/fallback limitations are unchanged.

## Evidence and reproduction

- [Candidate patch](2026-10-14-python-unambiguous-metal-model-candidate.patch), [source hashes](2026-10-14-python-unambiguous-metal-model-source.json), and [verification index](2026-10-14-python-unambiguous-metal-model-verification.json).
- [Slice 0](2026-10-14-python-unambiguous-metal-model-slice0-paired20.json) and [slice 5000](2026-10-14-python-unambiguous-metal-model-slice5000-paired20.json): every final block, load averages, versions, agreements, timings, RSS, and output hashes.
- [Self-contained tested source](../validation/results/python-unambiguous-metal-model-source-2026-10-14.tar.gz). Extract `candidate/`; to recover baseline, replace the two files from `baseline-overrides/`, remove `crates/chematic-smiles/src/rdkit/sanitize_metal_tests.rs`, and verify every baseline manifest hash. The prerequisite candidates are already included.
- Exact baseline/candidate descriptor and Morgan dumps, build/check logs, setup failure, and intermediate raw blocks are retained under `validation/results/python-unambiguous-metal-model-*`. Measured wheel paths and hashes are recorded; wheels are local build artifacts.

Build each source with `maturin build --release --manifest-path crates/chematic-py/Cargo.toml --interpreter python3.13`, install into separate target directories, then use the pinned comparator environment:

```sh
python scripts/bench_python_speed_frontier.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --offset 0 --limit 1000 --blocks 20 \
  --base-revision cac480fb23cd8adeec6e5e9bed1f18677f87b392-plus-unmerged-cycle-and-writer \
  --source-manifest benchmarks/2026-10-14-python-unambiguous-metal-model-source.json \
  --baseline-wheel <baseline.whl> --baseline-path <baseline-install> \
  --candidate-wheel <candidate.whl> --candidate-path <candidate-install> \
  --output <slice0.json>
```

Repeat with `--offset 5000`. Recheck all source and wheel hashes before and after timing.

## Separate adaptive ASA follow-up

The initial repeated-use ASA interval below 1 prompted one separate 20-block four-arm run on slice 0. The original matrix is retained, and results are not pooled.

| Mode | Baseline / candidate |
|---|---|
| labute_asa/first_use | 0.999 [0.986, 1.007]; 9/20 |
| labute_asa/repeated_use | 0.987 [0.983, 1.018]; 8/20 |
| labute_asa/parse_inclusive | 0.994 [0.981, 1.002]; 6/20 |
