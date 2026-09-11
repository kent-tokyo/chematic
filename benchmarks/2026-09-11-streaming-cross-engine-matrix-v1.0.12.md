# Streaming cross-engine matrix: v1.0.12 candidate

This is the 2026-09-11 refresh of the ten-format same-input record/failure
contract matrix. The machine-readable result is
[`cross-engine-matrix-v1.0.12.json`](../validation/results/cross-engine-matrix-v1.0.12.json).

The run used workspace version `1.0.12`, 20 repetitions per fixture, RDKit
`2025.09.3`, and Open Babel `3.2.1`. All available lanes returned the expected
record count with zero failures. The fixture SHA-256 values and input byte
counts are retained in the JSON result.

This is contract evidence, not a speed ranking. schematic uses a Rust
file-backed or materialized runner, RDKit uses Python block parsing where
available, and Open Babel starts a CLI process per repetition. A like-for-like
throughput gate and full binding-level semantic parity remain open.

Reproduce with:

```sh
cargo build -p chematic-mol --example streaming_benchmark --locked --offline
python3 scripts/check_streaming_cross_engine_matrix.py \
  --binary target/debug/examples/streaming_benchmark \
  --repeats 20 \
  --output validation/results/cross-engine-matrix-v1.0.12.json
```
