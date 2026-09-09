# Same-input Extended XYZ streaming contract — 2026-09-08

This v1.0.9 development-candidate contract check reads the checked-in
two-frame Extended XYZ fixture 20 times. It requires chematic and RDKit to
report the same valid-frame count, zero failures, and identical input-byte
accounting.

Fixture: `benchmarks/fixtures/streaming.extxyz`, 268 bytes,
SHA-256 `0094ab92ea84aa6acd826004e1849f2f2cb750c81d02e061c4156dc25c2830c7`.

| Engine | Boundary | Frames | Failures | Input bytes | Seconds | Frames/s |
|---|---|---:|---:|---:|---:|---:|
| chematic | Rust file-backed `BufRead` Extended XYZ reader | 40 | 0 | 5,360 | 0.001249 | 32,036.32 |
| RDKit 2025.09.3 | Python `MolFromXYZBlock` after frame splitting | 40 | 0 | 5,360 | 0.598595 | 66.82 |

The parser boundaries differ, so these numbers do not establish same-process
parity or a general throughput advantage. The contract only establishes that
the identical fixture produces the same frame/failure and byte accounting.

Reproduce:

```text
python3 scripts/check_streaming_cross_engine.py \
  --format extxyz --repeats 20 \
  --output /tmp/chematic-streaming-cross-engine-extxyz.json
```
