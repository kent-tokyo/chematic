# Streaming cross-engine matrix: v1.0.10 candidate

This is the 2026-09-09 refresh of the ten-format same-input contract matrix.
The complete machine-readable result is
[`2026-09-09-streaming-cross-engine-matrix-v1.0.10.json`](2026-09-09-streaming-cross-engine-matrix-v1.0.10.json).

The run used the current workspace version `1.0.10`, 20 repetitions per
fixture, RDKit `2025.09.3`, and Open Babel `3.2.1`. Every available lane
returned the expected record count with zero failures. Extended XYZ, CML,
CDXML, mmCIF, and PDB do not have an RDKit lane in this matrix.

The JSON includes `records_per_second` for operational context, but these
values are not a speed ranking: schematic uses a Rust file-backed or
materialized runner, RDKit uses Python block parsing, and Open Babel starts a
CLI process per repetition. A like-for-like throughput gate remains open.

| Format | schematic records/s | RDKit records/s | Open Babel records/s |
|---|---:|---:|---:|
| SDF | 11,450.52 | 37.70 | 5.29 |
| V2000 MOL | 26,680.01 | 36.35 | 5.86 |
| XYZ | 34,240.48 | 39.97 | 5.73 |
| Extended XYZ | 25,060.69 | — | 6.06 |
| V3000 MOL | 13,605.05 | 24.32 | 2.95 |
| MOL2 | 11,107.25 | 20.69 | 2.96 |
| CML | 9,543.31 | — | 2.92 |
| CDXML | 9,437.30 | — | 2.90 |
| mmCIF | 6,236.36 | — | 2.81 |
| PDB | 22,520.42 | — | 2.97 |

Reproduce with:

```sh
cargo build -p chematic-mol --example streaming_benchmark
python3 scripts/check_streaming_cross_engine_matrix.py \
  --binary target/debug/examples/streaming_benchmark \
  --repeats 20 \
  --output /tmp/chematic-streaming-cross-engine-matrix.json
```
