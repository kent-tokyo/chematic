# Python speed frontier: Morgan optimization — 2026-10-14

The three measured Morgan first-use lanes improve by 1.127–1.144x over the source baseline. Repeated-use improvements are 1.149–1.230x. All 56,288 Morgan output checks retain exact RDKit agreement. The complete speed target remains open: chiral Morgan, TPSA and Labute ASA still lose to COSMolKit in prepared calls, and first-use rings lose to RDKit.

## Scope and provenance

- Baseline: clean `ffe682c32cb4c76de0d309649427202a7d8e956d` source, version 1.1.3.
- Candidate: that commit plus the [archived six-file patch](https://github.com/kent-tokyo/chematic/blob/0bc130f24a7d163535e0f4ad59eb5a155847cf07/benchmarks/2026-10-14-python-speed-frontier-candidate.patch). The [source manifest](2026-10-14-python-speed-frontier-source.json) pins every changed source file and build settings.
- Both wheels were built with release optimization, LTO, one codegen unit and `native-inchi`, in a fresh target directory dedicated to this comparison. These are source-candidate measurements.
- Comparators: PyPI RDKit 2026.09.1 and COSMolKit 0.5.0rc22; installed extension hashes and runtime versions are retained in the [raw timing record](2026-10-14-python-speed-frontier-paired20.json).
- Host: macOS arm64, CPython 3.13.6, 10 logical CPUs. Load averages are recorded. There is no additional-host result.
- Timing corpus: first 1,000 rows of `validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi`, SHA-256 `f48777e96f74738336f47f682fbff5af6775a3474fe25daa28b595340feee95f`.
- Baseline wheel SHA-256: `33aa0d0d6e4cef0e0be7cdbe0b8f2482a8af8b1e9ab793ddf6bd4cc0ac28355e`.
- Candidate wheel SHA-256: `7e77ff905a393315fa4036063e349e44e9691a3596486404fad62ae356890709`.

## What changed

1. Compact 128-bit Morgan bond environments now preserve the general expansion's tetrahedral stereo hash and suppression lifecycle. Inputs with more than 128 bonds retain the general path.
2. Python fingerprint results are packed from little-endian 64-bit words instead of testing every bit. Variable widths preserve padding and bit numbering.
3. The narrow Fe(II) coordination check inspects neighbor charges only after confirming the relevant element, charge and degree.
4. The COSMolKit rc22 MACCS timing adapter uses its current `fingerprint_maccs()` method. Missing methods fail the timing run.

## Protocol

Twenty sequential, rotating four-arm blocks use fresh interpreters. Every arm occupies each position five times. Each operation retains native return objects; normalization occurs outside timing. Each block verifies every output, and the loaded source extensions must match the supplied wheels. An aggregation error was repaired after the timing blocks; a fresh validation pass reconfirmed all saved output and extension hashes before aggregation.

- First use: fresh parsed molecules, parsing excluded. RDKit/COSMolKit perform perception during parsing; chematic also has lazy perception costs here.
- Repeated use: fresh molecules receive one untimed operation each, then the operation is timed again.
- Parse inclusive: parsing and the operation run inside the same timed loop.

Ratios are comparator time / candidate time; above 1 favors chematic. Intervals are paired bootstrap 95% intervals. A strict demonstrated win additionally requires all 20 blocks to win and all compared outputs to agree. No failure, unsupported call or mismatching row is counted as a speed win. Numeric descriptor tolerance is absolute `1e-9`; fingerprints and ring counts require exact equality.

## First-use results

| Operation | Baseline / candidate | RDKit / candidate | COSMolKit / candidate |
|---|---:|---:|---:|
| TPSA | 0.994x [0.982, 1.012] | 1.352x [1.311, 1.378] | 0.460x [0.441, 0.474] |
| Labute ASA | 1.000x [0.994, 1.017] | 0.735x [0.686, 0.746] | 0.335x [0.314, 0.343] |
| Ring count | 0.992x [0.967, 1.001] | 0.019x [0.018, 0.020] | not ranked (6 output differences) |
| Morgan R2, plain | 1.144x [1.128, 1.155] | 3.132x [3.092, 3.195] | 2.022x [1.994, 2.049] |
| Morgan R2, chiral | 1.127x [1.124, 1.142] | 1.274x [1.264, 1.285] | 0.828x [0.814, 0.835] |
| Morgan R3, plain | 1.128x [1.116, 1.138] | 3.379x [3.340, 3.448] | 2.234x [2.196, 2.258] |
| MACCS | 1.003x [0.999, 1.008] | 1.888x [1.885, 1.893] | 2.180x [2.171, 2.208] |

Plain Morgan R2/R3 and MACCS pass the strict speed gate against both comparators. Chiral Morgan improves against baseline and beats RDKit in all blocks, but its COSMolKit ratio remains 0.828x. TPSA has a favorable RDKit interval but wins only 19/20 blocks. Labute ASA loses to both comparators. First-use ring count remains expensive; COSMolKit ring counts differ from the RDKit 2026.09.1 reference on 6/1,000 rows, so its timings are not ranked.

## Parse-inclusive results

| Operation | RDKit / candidate | COSMolKit / candidate |
|---|---:|---:|
| TPSA | 12.651x [12.358, 12.795] | 11.845x [11.643, 11.975] |
| Labute ASA | 12.321x [12.059, 12.601] | 11.612x [11.531, 11.805] |
| Ring count | 8.560x [8.469, 8.656] | not ranked (6 output differences) |
| Morgan R2, plain | 8.560x [8.247, 8.684] | 7.544x [7.462, 7.599] |
| Morgan R2, chiral | 4.461x [4.417, 4.490] | 3.942x [3.890, 3.957] |
| Morgan R3, plain | 7.904x [7.802, 8.024] | 6.894x [6.854, 6.968] |
| MACCS | 2.114x [2.104, 2.126] | 2.389x [2.377, 2.398] |

All seven pipelines pass the strict gate against RDKit. All six output-equivalent COSMolKit pipelines pass; ring count is not ranked. Baseline/candidate Morgan pipeline intervals favor the candidate, but only 18–19/20 blocks win, so they do not pass this packet's stricter every-block gate.

## Repeated use and remaining work

Plain Morgan R2 is 5.175x RDKit / 3.193x COSMolKit; plain R3 is 4.837x / 3.083x. Both pass the strict gate. Chiral R2 is 1.484x RDKit but 0.915x COSMolKit. Repeated ring count beats RDKit by 1.760x, while TPSA and Labute ASA remain slower than both comparators.

The next measured deficits are shared stereo preparation for chiral Morgan, descriptor preparation/result reuse for TPSA and Labute ASA, and initial ring perception. Reaction, 3D, batch/multithreaded, WASM, Linux/Windows and published-package measurements remain outside this packet. The exposed corpus also needs an independent holdout before a broader performance claim.

## Correctness and tests

- [Morgan parity record](2026-10-14-python-speed-frontier-parity.json): normal 10,000 plus unusual 4,072 inputs, radius 2/3, chirality off/on, 2,048 bits. Every one of the 56,288 cells matches both the baseline and RDKit 2026.09.1 exactly.
- Every timed operation matches the baseline on 1,000/1,000 rows. Candidate/RDKit agreement passes every timed lane; COSMolKit ring count is the six-row exception described above.
- `chematic-fp`: 339 unit tests and 6 integration tests passed on the isolated final source.
- Python fingerprints, bulk operations, similarity index, pickle, RDKit compatibility and fingerprint fixtures: 320 passed.
- Benchmark gate tests: 8 passed. Clippy passed with warnings denied for the final fingerprint/Python source and for all fingerprint targets; Rust formatting and diff checks passed.
- The shared working tree contained unrelated in-progress 3D/compatibility changes at the start. The isolated build keeps those changes outside the measured source.

## Reproduction

The archived patch reconstructs the six candidate files from the pinned baseline; their hashes were checked against the manifest. Build the baseline wheel first, apply the patch, then build the candidate with identical settings. Install each into a separate `pip --target` directory in a CPython 3.13 environment containing the pinned comparators.

```sh
python scripts/bench_python_speed_frontier.py \
  --corpus validation/benchmark_corpora/rdkit-js-browser-10k-v1.smi \
  --limit 1000 --blocks 20 \
  --baseline-wheel "$FRONTIER_BASELINE_WHEEL" --baseline-path "$FRONTIER_BASELINE_PATH" \
  --candidate-wheel "$FRONTIER_CANDIDATE_WHEEL" --candidate-path "$FRONTIER_CANDIDATE_PATH" \
  --base-revision ffe682c32cb4c76de0d309649427202a7d8e956d \
  --source-manifest benchmarks/2026-10-14-python-speed-frontier-source.json \
  --output frontier.json
```

[`verify_python_morgan_frontier.py`](../scripts/verify_python_morgan_frontier.py) regenerates the two-corpus R2/R3 byte dumps for each wheel and RDKit. The compressed [baseline/candidate dump](../validation/results/python-morgan-frontier-candidate-2026-10-14.jsonl.gz) is byte-identical for both source wheels; the separate [RDKit dump](../validation/results/python-morgan-frontier-rdkit-2026-10-14.jsonl.gz) carries its reference-version headers. Hashes are in the parity record.
