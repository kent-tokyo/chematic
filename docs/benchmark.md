# Benchmarks

Benchmark claims in chematic are operation-, corpus-, version-, and machine-
specific. Unsupported operations, failures, and non-equivalent APIs are never
counted as wins. Dated raw summaries and reproduction commands live in
the repository's [`benchmarks/`](https://github.com/kent-tokyo/chematic/tree/main/benchmarks) directory.

## Current status

The published release is **v1.0.5**. The newest hot-path measurements were
recorded before the release artifact was built and remain source-level evidence;
they must not be generalized beyond their named corpus and configuration.

### Published-source comparisons

| Operation | chematic | RDKit | Scope |
|---|---:|---:|---|
| Canonical SMILES | 24.95 µs/mol | 25.58 µs/mol | v1.0.2 code, 5,000-molecule descriptor corpus, medians, macOS arm64 |
| Canonical SMILES | 18.27 µs/mol | 26.82 µs/mol | v1.0.2 code, independent 5,000-entry ChEMBL-derived corpus |
| SDF graph/property read | 9.48 µs/mol | 99.96 µs/mol | v1.0.2 code, 365-record `egfr.sdf`, graph-only supplier |
| SDF serialization-only write | 7.62 µs/mol | 79.54 µs/mol | v1.0.2 code, same corpus, layout disabled |

These measurements show a lead only for the stated operations and environment.
Canonical strings need not match RDKit's spelling, and the SDF writer row does
not include automatic 2D layout. See the
[canonical](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-canonical-fast-path.md) and
[SDF](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-sdf-fast-path.md) records.

### Unreleased source-level follow-up

| Lane | v1.0.4 baseline | v1.0.5 source | Speedup |
|---|---:|---:|---:|
| Canonical SMILES, 5,000 molecules | 120.74 ms | 108.28 ms | **1.115x** |
| File-backed SDF read, 365 records x 50 | 145.949 ms | 129.201 ms | **1.130x** |
| V2000 SDF serialization, returned `String` | 3.669 µs/record | 1.206 µs/record | **3.042x** |
| V2000 SDF serialization, reusable buffer | 3.669 µs/record | 1.169 µs/record | **3.139x** |

The [follow-up record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-05-hot-path-follow-up.md) contains
the exact A/B protocol, rejected experiments, output checks, and commands.

## Accuracy and compatibility

Performance does not establish chemistry parity. The 4,999-molecule
ChEMBL-derived descriptor suite reports metric-specific agreement and known
residuals in [`validation.md`](validation.md). Canonical identity,
aromaticity/CIP modes, fingerprint definitions, and Experimental 3D/MMFF94
have separate contracts in [`compatibility-scope.md`](compatibility-scope.md).

## Artifact size

The optimized v1.0.2-candidate WASM artifact was measured at **3.30 MB raw /
1.21 MB gzip** with `wasm-pack 0.13.1` and `wasm-opt 130 -O3`. This is a dated
artifact measurement, not a permanent bundle-size guarantee. Exact hashes and
commands are in the [WASM artifact record](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-wasm-size.md).

## Reproduction index

| Record | Purpose |
|---|---|
| [2026-09-05 hot paths](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-05-hot-path-follow-up.md) | Fixed-version source A/B for canonical and SDF |
| [2026-09-04 RDKit/Open Babel](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-rdkit-openbabel.md) | In-process and CLI comparisons kept separate |
| [2026-09-04 canonical](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-canonical-fast-path.md) | Two 5,000-molecule canonical runs |
| [2026-09-04 SDF](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-sdf-fast-path.md) | Graph/property read and serialization-only write |
| [2026-09-04 streaming formats](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-streaming-formats.md) | File-backed SDF/MOL/XYZ runner |
| [2026-09-04 MMFF94/3D](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-mmff94-3d.md) | Experimental local microbenchmarks |
| [2026-09-04 WASM](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-04-wasm-size.md) | Artifact size and digest |
| [2026-09-03 competitive](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/2026-09-03-competitive.md) | Resumable six-operation v1.0.1 run |
| [Older snapshots](https://github.com/kent-tokyo/chematic/blob/main/benchmarks/README.md) | Earlier versioned measurements and hardware notes |

## Rules for new results

Every new result must record:

- source revision and package versions;
- corpus identity and hash;
- hardware, OS, language/runtime, and build profile;
- exact operation boundary and configuration;
- warm-up, repetitions, aggregation, failure policy, and raw output location;
- correctness or byte-equivalence checks relevant to the optimization.

Do not relabel source-level A/B data as a published artifact result, compare a
streaming API with a materializing API without saying so, or generalize one
corpus to all chemistry workloads.
