# Same-input V3000 streaming contract — 2026-09-08

This v1.0.9 development-candidate check reads the checked-in one-record V3000
fixture 20 times. It requires chematic and RDKit to report the same valid
record count, zero failures, and identical source-byte accounting.

Fixture: `benchmarks/fixtures/ethanol.v3000`, 298 bytes,
SHA-256 `64ac11b3de1c186deae555d54d6a47e7f3c72607ef0a90b0d93f467496bbd3ad`.

| Engine | Boundary | Records | Failures | Input bytes | Seconds | Records/s |
|---|---|---:|---:|---:|---:|---:|
| chematic | Rust V3000 materialized one-shot parser | 20 | 0 | 5,960 | 0.002563 | 7,802.21 |
| RDKit 2025.09.3 | Python `MolFromMolBlock` block parser | 20 | 0 | 5,960 | 0.350324 | 57.09 |

The parser boundaries differ, and chematic's current V3000 runner path is
materialized rather than a true streaming reader. These measurements do not
establish same-process parity or a throughput ranking; the contract covers
record, failure, and source-byte accounting only.

Reproduce:

```text
python3 scripts/check_streaming_cross_engine.py --format v3000 --repeats 20
```
