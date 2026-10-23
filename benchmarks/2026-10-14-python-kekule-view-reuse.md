# Python Kekulé and descriptor view reuse (2026-10-14)

Initial descriptor calculations spend time rebuilding aromatic representations. This candidate avoids a no-op Kekulé intermediate copy when the input has no aromatic atom flags or aromatic bonds and already has bond-ordered adjacency. Kekulization on such an input always succeeds with an empty result; the old intermediate was an exact clone. Aromaticity is assigned on the same original graph and uses the same verdict and output builder. Other inputs retain the original path. No computation is moved into parsing.

Descriptors also borrow the shared, immutable successful aromatic view when it already carries the no-further-dative marker and needs no ring-triple normalization. Externally seeded views without that marker retain the copy path. Inputs with aromatic endpoint triple bonds retain the previous normalization copy. Errors retain the original cleanup/general-model fallback. The original model gate, sanitizer, stereo cleanup, aromaticity verdict, ring limits and error boundaries are unchanged. The unsuccessful metal-cleanup experiment is retained separately and is not part of this candidate.


The complete matrix observes median baseline-relative first-use speedups of 1.183 for TPSA and 1.187 for ASA on slice 5,000, each faster in 19/20 blocks. Those do not meet the stricter 20/20 demonstrated-win criterion. The separate follow-up reproduces the slice-5,000 ring parse-inclusive slowdown (0.979 [0.962, 0.991], 4/20 faster); Morgan radius-3 repeat use remains slightly slower (0.997 [0.993, 0.999], 4/20). Slice-0 MACCS repeat use also has an interval below 1 in the follow-up. This is an isolated candidate with unresolved regressions, not a recommendation to merge or a claim of complete speed superiority.

## Source and correctness

Both 1.1.4 wheels are fresh local source builds. Baseline includes the preceding unmerged relevant-cycle and V3000 writer improvements; candidate changes only the two aromatic-view consumers, their tests, and the changelog. Last verified base: `cac480fb23cd8adeec6e5e9bed1f18677f87b392`. Shared Git metadata is unresponsive, so this run does not establish current-main, commit, release or registry-artifact status. The frozen source hashes and archive identify the tested source.

- Frozen pre-change aromatic view versus the new path on all 10,000 ordinary inputs, plus annotated/stereo/tag/R-group/group-bond and edited-bond cases. Atoms, indexed bonds, adjacency, hydrogen counts, metadata and success/error boundaries are checked. Frozen descriptor copy-path graphs are compared on the same corpus, including repeat calls and ring-triple cases.
- Native tests: chematic-chem 951 passed / 1 ignored; chematic-perception 257 passed / 1 ignored. Python 1,389 passed / 1 skipped. Release Clippy for both crates and all targets, workspace formatting and refactor budget pass. Missing copied test-data files initially stopped checks; complete reruns passed after copying the unchanged fixtures.
- All 14,072 ordinary/stress inputs retain byte-identical descriptor and four Morgan output streams against baseline, including errors. This verifies serialized outcomes; errors are not counted as valid chemistry.

## Measurement

CPython 3.13.6; Rust 1.97.0; maturin 1.13.3; macOS arm64; release defaults. RDKit 2026.09.1 and COSMolKit 0.5.0rc22. Fixed unfiltered 1,000-input slices at offsets 0 and 5,000; corpus SHA-256 `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`. These are exposed development inputs, not sealed holdout data. Existing provenance and licenses remain unchanged.

Each slice uses 20 four-arm blocks with rotating order in fresh processes. First use times the operation on fresh parsed objects; repeated use follows one untimed call; parse-inclusive times actual parsing plus the operation in one loop. Native results are retained; normalization happens outside timing. All rows, blocks, mismatches, hashes, versions and whole-process peak RSS are retained. RSS includes imports and all operations, not just operation allocation. Builds/tests finish before measurement; slices run sequentially without concurrent build/compression.

Tables show comparator time / candidate time, paired-median 95% bootstrap interval (10,000 deterministic resamples), and faster blocks. Above 1 means candidate faster. Strict wins require the lower interval endpoint above 1 and 20/20 faster blocks. Candidate/RDKit or COSMolKit/RDKit disagreement suppresses the corresponding rank; no mismatching rows are removed. One host and seven operations do not establish whole-toolkit superiority.

## Descriptor results

| Offset | Operation / mode | Fresh baseline / candidate | RDKit / candidate | COSMolKit / candidate |
|---:|---|---|---|---|
| 0 | tpsa/first_use | 0.989 [0.982, 0.999]; 5/20 | 1.474 [1.456, 1.490]; 20/20 | 0.511 [0.503, 0.516]; 0/20 |
| 0 | tpsa/parse_inclusive | 0.994 [0.975, 1.004]; 6/20 | 12.924 [12.765, 13.104]; 20/20 | 12.048 [11.944, 12.139]; 20/20 |
| 0 | labute_asa/first_use | 0.993 [0.990, 0.999]; 5/20 | 0.699 [0.681, 0.716]; 0/20 | 0.336 [0.331, 0.343]; 0/20 |
| 0 | labute_asa/parse_inclusive | 0.993 [0.986, 1.000]; 6/20 | 12.260 [12.144, 12.305]; 20/20 | 11.590 [11.550, 11.634]; 20/20 |
| 5000 | tpsa/first_use | 1.183 [1.172, 1.191]; 19/20 | 0.409 [0.404, 0.415]; 0/20 | 0.110 [0.108, 0.114]; 0/20 |
| 5000 | tpsa/parse_inclusive | 1.145 [1.122, 1.162]; 19/20 | 4.894 [4.842, 4.925]; 20/20 | 5.069 [4.996, 5.115]; 20/20 |
| 5000 | labute_asa/first_use | 1.187 [1.168, 1.206]; 19/20 | 0.197 [0.189, 0.203]; 0/20 | 0.080 [0.079, 0.083]; 0/20 |
| 5000 | labute_asa/parse_inclusive | 1.141 [1.129, 1.148]; 19/20 | 4.669 [4.579, 4.674]; 20/20 | 4.984 [4.911, 5.019]; 20/20 |

## Strict competitor wins

| Mode | Demonstrated / eligible |
|---|---:|
| first_use | 13/27 |
| repeated_use | 24/27 |
| parse_inclusive | 27/27 |

## Complete matrix

| Offset | Operation / mode | Fresh baseline / candidate | RDKit / candidate | COSMolKit / candidate |
|---:|---|---|---|---|
| 0 | tpsa/first_use | 0.989 [0.982, 0.999]; 5/20 | 1.474 [1.456, 1.490]; 20/20 | 0.511 [0.503, 0.516]; 0/20 |
| 0 | tpsa/repeated_use | 1.001 [0.991, 1.027]; 10/20 | 10.020 [9.799, 10.113]; 20/20 | 1.746 [1.621, 1.835]; 20/20 |
| 0 | tpsa/parse_inclusive | 0.994 [0.975, 1.004]; 6/20 | 12.924 [12.765, 13.104]; 20/20 | 12.048 [11.944, 12.139]; 20/20 |
| 0 | labute_asa/first_use | 0.993 [0.990, 0.999]; 5/20 | 0.699 [0.681, 0.716]; 0/20 | 0.336 [0.331, 0.343]; 0/20 |
| 0 | labute_asa/repeated_use | 0.991 [0.983, 1.017]; 8/20 | 2.709 [2.636, 2.753]; 20/20 | 1.383 [1.355, 1.432]; 20/20 |
| 0 | labute_asa/parse_inclusive | 0.993 [0.986, 1.000]; 6/20 | 12.260 [12.144, 12.305]; 20/20 | 11.590 [11.550, 11.634]; 20/20 |
| 0 | num_rings/first_use | 0.998 [0.993, 1.000]; 6/20 | 0.034 [0.033, 0.035]; 0/20 | unranked |
| 0 | num_rings/repeated_use | 1.001 [0.986, 1.013]; 10/20 | 1.746 [1.709, 1.778]; 20/20 | unranked |
| 0 | num_rings/parse_inclusive | 0.993 [0.979, 1.008]; 9/20 | 11.207 [11.154, 11.296]; 20/20 | unranked |
| 0 | morgan2_2048/first_use | 0.990 [0.986, 0.995]; 5/20 | 3.103 [3.049, 3.118]; 20/20 | 2.006 [1.988, 2.015]; 20/20 |
| 0 | morgan2_2048/repeated_use | 0.993 [0.985, 0.997]; 5/20 | 5.202 [5.147, 5.234]; 20/20 | 3.186 [3.168, 3.220]; 20/20 |
| 0 | morgan2_2048/parse_inclusive | 0.998 [0.989, 1.006]; 9/20 | 8.589 [8.496, 8.628]; 20/20 | 7.540 [7.469, 7.579]; 20/20 |
| 0 | morgan2_chiral/first_use | 1.000 [0.989, 1.006]; 10/20 | 1.262 [1.256, 1.268]; 20/20 | 0.820 [0.812, 0.826]; 0/20 |
| 0 | morgan2_chiral/repeated_use | 0.995 [0.989, 1.002]; 7/20 | 5.052 [4.974, 5.087]; 20/20 | 3.099 [3.069, 3.114]; 20/20 |
| 0 | morgan2_chiral/parse_inclusive | 1.004 [0.996, 1.007]; 13/20 | 4.450 [4.421, 4.480]; 20/20 | 3.917 [3.884, 3.956]; 20/20 |
| 0 | morgan3_2048/first_use | 0.997 [0.984, 0.999]; 5/20 | 3.357 [3.341, 3.389]; 20/20 | 2.234 [2.221, 2.246]; 20/20 |
| 0 | morgan3_2048/repeated_use | 0.992 [0.986, 0.998]; 4/20 | 4.850 [4.800, 4.946]; 20/20 | 3.123 [3.086, 3.154]; 20/20 |
| 0 | morgan3_2048/parse_inclusive | 0.994 [0.988, 1.000]; 6/20 | 7.974 [7.927, 7.999]; 20/20 | 6.945 [6.907, 6.977]; 20/20 |
| 0 | maccs/first_use | 0.997 [0.996, 0.999]; 3/20 | 1.877 [1.869, 1.886]; 20/20 | 2.164 [2.157, 2.168]; 20/20 |
| 0 | maccs/repeated_use | 0.997 [0.993, 1.001]; 7/20 | 1.946 [1.938, 1.953]; 20/20 | 2.240 [2.237, 2.246]; 20/20 |
| 0 | maccs/parse_inclusive | 0.999 [0.996, 1.004]; 8/20 | 2.104 [2.096, 2.113]; 20/20 | 2.376 [2.371, 2.382]; 20/20 |
| 5000 | tpsa/first_use | 1.183 [1.172, 1.191]; 19/20 | 0.409 [0.404, 0.415]; 0/20 | 0.110 [0.108, 0.114]; 0/20 |
| 5000 | tpsa/repeated_use | 0.994 [0.935, 1.044]; 9/20 | 9.709 [8.724, 10.242]; 20/20 | 1.319 [1.264, 1.368]; 19/20 |
| 5000 | tpsa/parse_inclusive | 1.145 [1.122, 1.162]; 19/20 | 4.894 [4.842, 4.925]; 20/20 | 5.069 [4.996, 5.115]; 20/20 |
| 5000 | labute_asa/first_use | 1.187 [1.168, 1.206]; 19/20 | 0.197 [0.189, 0.203]; 0/20 | 0.080 [0.079, 0.083]; 0/20 |
| 5000 | labute_asa/repeated_use | 1.012 [0.997, 1.059]; 14/20 | 2.638 [2.539, 2.679]; 20/20 | 1.214 [1.145, 1.324]; 18/20 |
| 5000 | labute_asa/parse_inclusive | 1.141 [1.129, 1.148]; 19/20 | 4.669 [4.579, 4.674]; 20/20 | 4.984 [4.911, 5.019]; 20/20 |
| 5000 | num_rings/first_use | 0.976 [0.932, 1.004]; 6/20 | 0.090 [0.089, 0.096]; 0/20 | 0.072 [0.067, 0.075]; 0/20 |
| 5000 | num_rings/repeated_use | 1.015 [0.973, 1.041]; 12/20 | 1.680 [1.626, 1.764]; 20/20 | 1.040 [0.995, 1.066]; 14/20 |
| 5000 | num_rings/parse_inclusive | 0.980 [0.965, 0.998]; 6/20 | 11.476 [11.268, 11.548]; 20/20 | 12.413 [12.095, 12.554]; 20/20 |
| 5000 | morgan2_2048/first_use | 1.056 [1.042, 1.067]; 18/20 | 1.482 [1.465, 1.495]; 20/20 | 0.942 [0.932, 0.956]; 0/20 |
| 5000 | morgan2_2048/repeated_use | 0.997 [0.987, 1.015]; 10/20 | 4.704 [4.646, 4.747]; 20/20 | 2.860 [2.815, 2.895]; 20/20 |
| 5000 | morgan2_2048/parse_inclusive | 1.051 [1.041, 1.064]; 18/20 | 4.507 [4.452, 4.522]; 20/20 | 4.307 [4.287, 4.331]; 20/20 |
| 5000 | morgan2_chiral/first_use | 1.063 [1.049, 1.072]; 18/20 | 1.522 [1.497, 1.548]; 20/20 | 0.974 [0.965, 0.984]; 3/20 |
| 5000 | morgan2_chiral/repeated_use | 0.998 [0.984, 1.009]; 10/20 | 4.758 [4.669, 4.807]; 20/20 | 2.880 [2.841, 2.902]; 20/20 |
| 5000 | morgan2_chiral/parse_inclusive | 1.051 [1.032, 1.060]; 18/20 | 4.474 [4.364, 4.500]; 20/20 | 4.314 [4.263, 4.347]; 20/20 |
| 5000 | morgan3_2048/first_use | 1.050 [1.031, 1.069]; 18/20 | 1.712 [1.644, 1.727]; 20/20 | 1.113 [1.094, 1.136]; 17/20 |
| 5000 | morgan3_2048/repeated_use | 0.992 [0.984, 0.998]; 5/20 | 4.524 [4.341, 4.574]; 20/20 | 2.868 [2.829, 2.896]; 20/20 |
| 5000 | morgan3_2048/parse_inclusive | 1.045 [1.034, 1.060]; 16/20 | 4.393 [4.325, 4.429]; 20/20 | 4.177 [4.135, 4.187]; 20/20 |
| 5000 | maccs/first_use | 0.999 [0.996, 1.006]; 9/20 | 1.623 [1.616, 1.629]; 20/20 | 1.741 [1.735, 1.748]; 20/20 |
| 5000 | maccs/repeated_use | 1.002 [0.992, 1.007]; 11/20 | 1.728 [1.714, 1.744]; 20/20 | 1.847 [1.842, 1.859]; 20/20 |
| 5000 | maccs/parse_inclusive | 1.002 [0.995, 1.007]; 11/20 | 1.797 [1.785, 1.807]; 20/20 | 1.921 [1.914, 1.932]; 20/20 |

## Regressions and limits

Baseline-relative intervals wholly below 1:
- Offset 0, `tpsa/first_use`: 0.989 [0.982, 0.999]; 5/20
- Offset 0, `labute_asa/first_use`: 0.993 [0.990, 0.999]; 5/20
- Offset 0, `morgan2_2048/first_use`: 0.990 [0.986, 0.995]; 5/20
- Offset 0, `morgan2_2048/repeated_use`: 0.993 [0.985, 0.997]; 5/20
- Offset 0, `morgan3_2048/first_use`: 0.997 [0.984, 0.999]; 5/20
- Offset 0, `morgan3_2048/repeated_use`: 0.992 [0.986, 0.998]; 4/20
- Offset 0, `maccs/first_use`: 0.997 [0.996, 0.999]; 3/20
- Offset 5000, `num_rings/parse_inclusive`: 0.980 [0.965, 0.998]; 6/20
- Offset 5000, `morgan3_2048/repeated_use`: 0.992 [0.984, 0.998]; 5/20

One separate adaptive 20-block follow-up was selected from the complete matrix for operations with baseline-relative intervals wholly below 1. No pooling, discarded blocks, or recursive retry is used. All its modes follow:

| Offset | Operation / mode | Baseline / candidate | RDKit / candidate | COSMolKit / candidate |
|---:|---|---|---|---|
| 0 | labute_asa/first_use | 0.994 [0.985, 1.000]; 6/20 | 0.665 [0.653, 0.685]; 0/20 | 0.311 [0.308, 0.316]; 0/20 |
| 0 | labute_asa/repeated_use | 1.006 [0.988, 1.027]; 10/20 | 2.704 [2.656, 2.784]; 20/20 | 1.348 [1.311, 1.397]; 19/20 |
| 0 | labute_asa/parse_inclusive | 0.991 [0.981, 1.002]; 8/20 | 12.058 [11.997, 12.181]; 20/20 | 11.449 [11.340, 11.528]; 20/20 |
| 0 | maccs/first_use | 0.999 [0.997, 1.001]; 7/20 | 1.883 [1.878, 1.885]; 20/20 | 2.167 [2.160, 2.170]; 20/20 |
| 0 | maccs/repeated_use | 0.995 [0.992, 0.998]; 4/20 | 1.945 [1.937, 1.950]; 20/20 | 2.237 [2.232, 2.243]; 20/20 |
| 0 | maccs/parse_inclusive | 1.002 [0.996, 1.003]; 11/20 | 2.106 [2.097, 2.112]; 20/20 | 2.374 [2.366, 2.383]; 20/20 |
| 0 | morgan2_2048/first_use | 0.994 [0.987, 1.004]; 8/20 | 3.044 [3.015, 3.064]; 20/20 | 1.986 [1.977, 2.005]; 20/20 |
| 0 | morgan2_2048/repeated_use | 0.998 [0.994, 1.004]; 9/20 | 5.201 [5.151, 5.232]; 20/20 | 3.197 [3.166, 3.224]; 20/20 |
| 0 | morgan2_2048/parse_inclusive | 0.999 [0.991, 1.006]; 10/20 | 8.595 [8.558, 8.660]; 20/20 | 7.576 [7.505, 7.623]; 20/20 |
| 0 | morgan3_2048/first_use | 0.996 [0.982, 1.007]; 8/20 | 3.343 [3.336, 3.387]; 20/20 | 2.249 [2.234, 2.265]; 20/20 |
| 0 | morgan3_2048/repeated_use | 0.996 [0.992, 1.003]; 7/20 | 4.855 [4.800, 4.900]; 20/20 | 3.110 [3.100, 3.130]; 20/20 |
| 0 | morgan3_2048/parse_inclusive | 0.996 [0.993, 1.008]; 8/20 | 7.979 [7.924, 8.007]; 20/20 | 6.953 [6.880, 6.972]; 20/20 |
| 0 | tpsa/first_use | 0.994 [0.982, 1.005]; 7/20 | 1.552 [1.529, 1.582]; 20/20 | 0.529 [0.524, 0.536]; 0/20 |
| 0 | tpsa/repeated_use | 1.006 [0.993, 1.035]; 12/20 | 9.949 [9.847, 10.030]; 20/20 | 1.721 [1.635, 1.773]; 20/20 |
| 0 | tpsa/parse_inclusive | 0.998 [0.977, 1.009]; 9/20 | 13.047 [12.951, 13.201]; 20/20 | 12.182 [12.092, 12.251]; 20/20 |
| 5000 | morgan3_2048/first_use | 1.062 [1.057, 1.064]; 20/20 | 1.598 [1.591, 1.616]; 20/20 | 1.065 [1.060, 1.073]; 20/20 |
| 5000 | morgan3_2048/repeated_use | 0.997 [0.993, 0.999]; 4/20 | 4.577 [4.531, 4.617]; 20/20 | 2.865 [2.844, 2.874]; 20/20 |
| 5000 | morgan3_2048/parse_inclusive | 1.039 [1.033, 1.047]; 20/20 | 4.371 [4.344, 4.398]; 20/20 | 4.142 [4.121, 4.206]; 20/20 |
| 5000 | num_rings/first_use | 0.992 [0.982, 1.002]; 8/20 | 0.098 [0.095, 0.100]; 0/20 | 0.074 [0.071, 0.081]; 0/20 |
| 5000 | num_rings/repeated_use | 1.012 [0.977, 1.033]; 11/20 | 1.739 [1.670, 1.802]; 20/20 | 1.080 [1.028, 1.106]; 17/20 |
| 5000 | num_rings/parse_inclusive | 0.979 [0.962, 0.991]; 4/20 | 11.223 [11.106, 11.317]; 20/20 | 12.289 [12.125, 12.405]; 20/20 |

First-use descriptor and ring-count gaps remain visible above. Prior [root-pruning](2026-10-14-python-cycle-root-pruning.md) measurements retained a chiral-Morgan parse-inclusive regression against that run's baseline. This candidate includes root pruning in both arms and therefore cannot establish that its earlier regression against main has been resolved. Existing bounded ring/fallback chemistry limitations are unchanged.

## Evidence and reproduction

- [Candidate patch](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-14-python-kekule-view-reuse-candidate.patch), [source manifest](2026-10-14-python-kekule-view-reuse-source.json), [verification/hashes](2026-10-14-python-kekule-view-reuse-verification.json), [slice 0](2026-10-14-python-kekule-view-reuse-slice0-paired20.json), [slice 5000](2026-10-14-python-kekule-view-reuse-slice5000-paired20.json).
- [Self-contained source archive](../validation/results/python-kekule-view-reuse-source-2026-10-14.tar.gz): extract `candidate/`; for baseline, restore all three `baseline-overrides/` files, remove both new test files listed in the manifest, and verify the baseline hashes before building. Prior candidate prerequisites are included.
- Exact baseline/candidate output dumps and all build/check logs: `validation/results/python-kekule-view-reuse-*`. Measured local wheel hashes and paths are recorded; these are not published wheels.
- [Unadopted metal-cleanup experiment](2026-10-14-python-unambiguous-metal-model.md): its source, complete full matrices, separate adaptive repeat, outcomes, setup failures and intermediate blocks remain available. Phase diagnostics motivated the change in direction; their unpaired timings are not speedup evidence.

Build each source with `maturin build --release --manifest-path crates/chematic-py/Cargo.toml --interpreter python3.13`; install into separate target directories. In the pinned comparator environment:

```sh
python scripts/bench_python_speed_frontier.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --offset 0 --limit 1000 --blocks 20 \
  --base-revision cac480fb23cd8adeec6e5e9bed1f18677f87b392-plus-unmerged-cycle-and-writer \
  --source-manifest benchmarks/2026-10-14-python-kekule-view-reuse-source.json \
  --baseline-wheel <baseline.whl> --baseline-path <baseline-install> \
  --candidate-wheel <candidate.whl> --candidate-path <candidate-install> \
  --output <slice0.json>
```

Repeat with `--offset 5000`. Recheck source, wheel and corpus hashes before and after timing.
