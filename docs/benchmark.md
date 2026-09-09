# Benchmarks

Benchmark claims in chematic are operation-, corpus-, version-, and machine-
specific. Unsupported operations, failures, and non-equivalent APIs are never
counted as wins. Dated raw records are indexed in
[`benchmarks/README.md`](../benchmarks/README.md).

## Current status

The current published release is **v1.0.11**. The newest checked-in performance
records are historical source-level or release-line measurements; they
must not be presented as v1.0.11 measurements unless explicitly rerun.

The current records cover four separate evidence types:

1. performance: operation timing, scaling, and source A/B comparisons;
2. compatibility: record accounting, failure behavior, and cross-engine contracts;
3. artifacts: WASM size, wheel installation, and other release-adjacent evidence;
4. validation: chemistry parity and correctness, which is not interchangeable with speed.

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

The [2026-09-09 similarity-search record](../benchmarks/2026-09-09-similarity-search-v1.0.9.md)
is a v1.0.9 historical source measurement using a fixed 4,500-entry library and
500-query split against RDKit. Search latency and top-k ranking overlap are
reported separately because native chematic ECFP4 and RDKit Morgan fingerprints
are not bit-identical. The result is not a claim that one engine is universally
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
[`Streaming and cross-engine contracts`](../benchmarks/README.md#streaming-and-cross-engine-contracts).
Use [`scripts/validate_streaming_cross_engine_matrix.py`](../scripts/validate_streaming_cross_engine_matrix.py)
for the fail-closed matrix check.

## Artifact size

The optimized v1.0.10 release-line WASM artifact was measured at **3.73 MB raw /
1.36 MB gzip** with `wasm-pack 0.13.1` and `wasm-opt 130`. This is a dated
artifact measurement, not a permanent release bundle-size guarantee. Exact
hashes and commands are in the
[WASM artifact record](../benchmarks/2026-09-09-wasm-size-v1.0.10.md).

## Reproduction entry points

| Purpose | Entry point |
|---|---|
| Accuracy vs RDKit | `pip install chematic rdkit`; `python scripts/bench5k.py scripts/chembl_accuracy_corpus_4999.smi --json /tmp/bench5k.json`; `python scripts/gen_validation_report.py /tmp/bench5k.json` |
| Similarity search | [`2026-09-09-similarity-search-v1.0.9.md`](../benchmarks/2026-09-09-similarity-search-v1.0.9.md) |
| Hot-path A/B | [`2026-09-05-hotpath-110.md`](../benchmarks/2026-09-05-hotpath-110.md) |
| File streaming | [`2026-09-04-streaming-formats.md`](../benchmarks/2026-09-04-streaming-formats.md) |
| Cross-engine contracts | [`benchmarks/README.md`](../benchmarks/README.md#streaming-and-cross-engine-contracts) |
| WASM artifact | [`2026-09-09-wasm-size-v1.0.10.md`](../benchmarks/2026-09-09-wasm-size-v1.0.10.md) |
| Full inventory | [`benchmarks/README.md`](../benchmarks/README.md) |

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
[`validation/results/streaming-failure-taxonomy-v1.0.10.json`](../validation/results/streaming-failure-taxonomy-v1.0.10.json)
and can be regenerated with `python3 scripts/check_streaming_failure_taxonomy.py`
after building the streaming example.
