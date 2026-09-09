# Same-input XYZ streaming contract — 2026-09-08

This is a v1.0.9 development-candidate contract check, not a speed ranking.
It reads the checked-in two-frame XYZ fixture 20 times and requires chematic
and RDKit to report the same valid-record count, zero failures, and identical
input-byte accounting.

Fixture: `benchmarks/fixtures/streaming.xyz`, 161 bytes,
SHA-256 `3c717e4520a3e1b7b90ace86b72e73cce4dda38de7ae29aa811733d2ec227221`.

| Engine | Boundary | Records | Failures | Input bytes | Seconds | Records/s |
|---|---|---:|---:|---:|---:|---:|
| chematic | Rust file-backed `BufRead` XYZ reader | 40 | 0 | 3,220 | 0.000870 | 45,968.19 |
| RDKit 2025.09.3 | Python `MolFromXYZBlock` after frame splitting | 40 | 0 | 3,220 | 0.351125 | 113.92 |

The parser boundaries differ, so these numbers do not establish same-process
parity or a general throughput advantage. The contract only establishes that
the identical fixture produces the same record/failure accounting.

Reproduce:

```text
python3 scripts/check_streaming_cross_engine.py --format xyz --repeats 20
```
