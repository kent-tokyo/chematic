# v1.0.4-based hot-path follow-up

Measured 2026-09-05 on Apple arm64, macOS 26.5.2. The optimized side is an
unreleased local working tree based on chematic v1.0.4; the package version was
kept at `1.0.4`. These numbers must not be presented as measurements of a
published artifact.

## Result

All three targeted lanes exceeded the additional 1.10x source-level speed
gate.

| Lane | v1.0.4 baseline median | Optimized median | Speedup |
| --- | ---: | ---: | ---: |
| Canonical SMILES, 5,000 molecules | 120.74 ms | 108.28 ms | **1.115x** |
| File-backed SDF graph/property read, 365 records x 50 | 145.949 ms | 129.201 ms | **1.130x** |
| V2000 SDF serialization, returned `String` | 3.669 us/record | 1.206 us/record | **3.042x** |
| V2000 SDF serialization, reusable buffer | 3.669 us/record | 1.169 us/record | **3.139x** |

The SDF reader's optimized measurements were taken with 100 repeats and divided
by two to match the 50-repeat baseline work count. The raw optimized median was
258.402 ms for 36,500 records. Writer rows compare the pre-change v1.0.4
serializer median with the retained implementation; the reusable-buffer row is
the path used by Python `SDWriter`.

## Retained changes

- Canonical E/Z preparation now skips double bonds that cannot contribute an
  ambiguous two-substituent carrier. This avoids unnecessary ring perception
  and marker-resolution work without changing the candidate set used by the
  writer.
- The file-backed SDF reader accumulates bytes and validates UTF-8 once per
  record, instead of once per physical line. Fixed-width count and common zero-Z
  fields use narrow fast paths with the previous generic parsers retained as
  fallbacks.
- The graph-only V2000 parser returns after the declared atom and bond blocks;
  the diagnostic parser retains the full `M  END` scan.
- V2000 counts, zero-coordinate atom tails, and bond rows use allocation-free
  fixed-width integer emission. Streaming writers reuse one caller-owned record
  buffer.

Two experiments were rejected: direct growth into `String` in the SDF reader
regressed throughput by about 2%, and an alternate canonical final-rank
normalization was neutral (120.99 ms versus 121.08 ms). Neither is retained.

## Measurement protocol

Canonical results are medians of 11 alternating baseline/current process pairs
over `scripts/descriptor_census_corpus.smi`. Parsing is outside the Tier C
canonical timing. The baseline executable was built from v1.0.4 before the
working-tree change; current and baseline processes were alternated to reduce
load-order bias.

SDF read results are medians of independent release-mode processes using
RDKit's 365-record `egfr.sdf`. SDF write results are medians of nine rounds over
the same corpus, 50 repetitions per round. The checked-in writer benchmark
compares its internal implementations in alternating order and asserts
byte-for-byte equality before timing.

```bash
# Canonical throughput; prints a deterministic output digest as an additional
# accidental-output-change guard.
cargo run --release -p chematic-smiles --example canonical_throughput -- \
  scripts/descriptor_census_corpus.smi 11

# File-backed SDF read.
cargo run --release -p chematic-mol --example streaming_benchmark -- \
  --format sdf --path /path/to/egfr.sdf --repeats 100

# SDF serializer comparison and reusable-buffer lane.
cargo run --release -p chematic-mol --example sdf_writer_ab -- \
  /path/to/egfr.sdf
```

## Correctness and static gates

- `cargo test --workspace --all-targets --locked`: pass; no failures. Existing
  explicitly ignored long-running or known-gap tests remain ignored.
- Python focused binding suite: 284 passed (`test_io.py`,
  `test_rdkit_compat.py`, `test_canonical_diff.py`, and
  `test_cross_binding_contract.py`).
- `cargo fmt --all -- --check`: pass.
- `cargo clippy --workspace --all-targets --all-features --locked -- -D warnings`:
  pass.
- New regressions verify fixed-width output against Rust formatting, reusable
  SDF output byte identity and buffer clearing, and invalid-UTF-8 rejection in
  the byte-buffered file reader.

The public Python benchmark was not used for the source-level speed claim in
this record because macOS `StorageManagementService` consumed substantial CPU
during the final session. Core release executables and alternating A/B runs
provide the reported evidence; a clean-idle Python rerun remains useful before
the next artifact release.
