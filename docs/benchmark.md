# Benchmarks

Benchmark claims in chematic are operation-, corpus-, version-, and machine-
specific. Unsupported operations, failures, and non-equivalent APIs are never
counted as wins. Dated raw records are indexed in
[benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks).

## Current status

The current release line is **v1.0.21**. The newest completed public-package performance
record remains the 2026-09-23 v1.0.20 browser/3D comparison against RDKit 2026.03.6;
v1.0.21 has no new package-performance remeasurement. Older
similarity, streaming, and operation timing records remain pinned to their
recorded source/release versions.

The current records cover four separate evidence types:

1. performance: operation timing, scaling, and source A/B comparisons;
2. compatibility: record accounting, failure behavior, and cross-engine contracts;
3. artifacts: WASM size, wheel installation, and other release-adjacent evidence;
4. validation: chemistry parity and correctness, which is not interchangeable with speed.

## Current Parse + compatible Morgan gate

Merge commit `7d98dcd3` is faster than pinned `@rdkit/rdkit@2026.3.6` for the
declared browser operation: SMILES parse + radius-2/2048-bit compatible Morgan
fingerprint + the same consumed 256-byte packed output. In GitHub-hosted
Chromium, Firefox, and WebKit, the paired 95% lower speedup bounds are **3.11x,
2.36x, and 3.67x**, respectively. The fixed corpus is exact on all 9,999
supported rows; one Fe(II) coordination structure remains an explicit typed
refusal. A separate ChEMBL 5k Chromium gate is 5,000/5,000 exact.

The registry-installed v1.0.20 npm package passes both fixed-corpus lanes:
1.398x parse-inclusive (95% lower bound 1.363x) and 3.511x prepared (lower
bound 3.407x), with 9,999/9,999 supported rows bit-exact. The same release's
PyPI wheel reaches 265/265 independently sound, stereo-clean and clash-free
MMFF94 stereo-safe outputs, but that quality-equivalent lane is only 0.944x
RDKit speed (lower bound 0.861x). See the complete
[v1.0.20 public-package record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-23-public-package-fingerprint-3d-v1.0.20.md).

## Published-source timing summary

| Operation | chematic | RDKit | Scope |
|---|---:|---:|---|
| Canonical SMILES | 24.95 µs/mol | 25.58 µs/mol | v1.0.2 code, 5,000-molecule corpus, macOS arm64 |
| Canonical SMILES | 18.27 µs/mol | 26.82 µs/mol | Independent 5,000-entry ChEMBL-derived corpus |
| SDF graph/property read | 9.48 µs/mol | 99.96 µs/mol | 365-record `egfr.sdf`, graph-only supplier |
| SDF serialization-only write | 7.62 µs/mol | 79.54 µs/mol | Same corpus; automatic 2D layout excluded |

These rows show a lead only for the named operation and environment. Canonical
strings need not match RDKit's spelling, and a writer-only measurement is not a
full depiction or export benchmark. See the dated canonical and SDF records.

## Source A/B timing summary

The v1.0.7 hot-path record reports alternating paired medians:

| Lane | Speedup |
|---|---:|
| Canonical SMILES | 1.176x |
| File-backed SDF read | 1.180x |
| Reused-buffer SDF serialization | 1.419x |
| SMILES parse | 1.034x; below the 1.10x target |

The historical v1.0.5 follow-up reports separate canonical, SDF-read, and
V2000-write A/B results. Rejected experiments, raw pairs, exact-output checks,
and load caveats remain in the linked records; do not merge their numbers with
the current release summary.

## Similarity search

The [2026-09-11 similarity-search record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-11-similarity-search-v1.0.12.md)
is the current three-lane gate using a fixed 4,500-entry library and 500-query
split. Native/native, RDKit-compatible/RDKit, and cross-profile overlap are
reported separately; failures are excluded from the valid-input scope and
counted independently. The result is not a claim that one engine is universally
faster or more accurate.

## Accuracy and compatibility

Performance does not establish chemistry parity. The 4,999-molecule
ChEMBL-derived descriptor suite reports metric-specific agreement and known
residuals in [`validation.md`](validation.md). Canonical identity,
aromaticity/CIP modes, fingerprint definitions, and Experimental 3D/MMFF94 have
separate contracts in [`compatibility-scope.md`](compatibility-scope.md).

Streaming records primarily measure record accounting, failure behavior, and
boundary semantics. The ten-format Rust and cross-engine matrices, including
plain/gzip stages, are indexed under
[Streaming and cross-engine contracts](https://github.com/kent-tokyo/chematic/tree/main/benchmarks).
Use [`scripts/validate_streaming_cross_engine_matrix.py`](https://github.com/kent-tokyo/chematic/blob/main/scripts/validate_streaming_cross_engine_matrix.py)
for the fail-closed matrix check.

## Artifact size

The published v1.0.15 npm WASM asset is **4,005,280 bytes raw / 1,460,499
bytes gzip**. Official RDKit.js 2026.03.6 is **7,333,095 / 2,379,975 bytes**
under the same local file-compression method. This is an artifact measurement,
not an internet-transfer or full application-size claim. Exact hashes and
commands are in the [isolated browser comparison](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md).

## Reproduction entry points

| Purpose | Entry point |
|---|---|
| Accuracy vs RDKit | `pip install chematic rdkit`; `python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi --json /tmp/bench5k.json`; `python scripts/gen_validation_report.py /tmp/bench5k.json` |
| Similarity search | [v1.0.12 record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-11-similarity-search-v1.0.12.md) |
| Parse + compatible Morgan vs official RDKit.js | [2026-09-20 record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-20-parse-morgan-rdkitjs.md) |
| Published-package fingerprint and 3D comparison | [v1.0.20 record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-23-public-package-fingerprint-3d-v1.0.20.md) |
| A3 similarity search rerun | [v1.0.13 record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-12-similarity-search-a3-v1.0.13.md) |
| Hot-path A/B | [2026-09-05 record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-05-hotpath-110.md) |
| File streaming | [2026-09-04 record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-streaming-formats.md) |
| Cross-engine contracts | [benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks) |
| WASM artifact and official RDKit.js gate | [v1.0.12 record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-11-official-rdkit-js-v1.0.12.md) |
| Official RDKit.js browser gate | [v1.0.15 record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md) |
| Full inventory | [benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks) |

## Hardware and interpretation

Hardware varies by snapshot; read each record's header. The 2026-06 records
used Apple M2 / macOS 14, while later records moved to Apple M4 / macOS 26.
Rust `cargo bench` numbers are comparable to Python numbers only when the
record says the machine, corpus, and operation boundary match.

Results vary by CPU, load, corpus, and fixture choice. Treat them as scoped
measurements, not SLA guarantees. The 2026-07 snapshot demonstrated that a
headline throughput ratio can be fixture-sensitive.

## Rules for new results

Every new result must record:

- source revision and package versions;
- corpus identity and hash;
- hardware, OS, language/runtime, and build profile;
- exact operation boundary and configuration;
- warm-up, repetitions, aggregation, failure policy, and raw output location;
- correctness, ranking, or byte-equivalence checks relevant to the operation.

Do not relabel source-level A/B data as a published artifact result, compare a
streaming API with a materializing API without saying so, or generalize one
corpus to all chemistry workloads.
### Streaming parser failure taxonomy

`streaming_benchmark` reports a deterministic `failure_kinds` object for all
ten supported formats. For parser errors, the keys are the Rust parser error
variant names and the values are counts accumulated across repetitions;
resource-limit failures are recorded as `ResourceLimit`. This field is
diagnostic coverage evidence, not a claim of exhaustive parser-state coverage
or cross-engine error taxonomy parity.
The 120-case bounded taxonomy result is stored in
[`validation/results/streaming-failure-taxonomy-v1.0.10.json`](https://github.com/kent-tokyo/chematic/blob/main/validation/results/streaming-failure-taxonomy-v1.0.10.json)
and can be regenerated with `python3 scripts/check_streaming_failure_taxonomy.py`
after building the streaming example.
