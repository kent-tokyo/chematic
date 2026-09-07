# Same-input V2000 MOL streaming contract — 2026-09-08

This v1.0.9 development-candidate check reads the two V2000 MOL records in
the checked-in SDF fixture 20 times. It requires chematic and RDKit to report
the same valid-record count, zero failures, and identical source-byte
accounting.

Fixture: `benchmarks/fixtures/streaming.sdf`, 633 bytes,
SHA-256 `627e7814694f7d00a9212712ec697f2174d30327b459d983d2a10f717cae0465`.

| Engine | Boundary | Records | Failures | Input bytes | Seconds | Records/s |
|---|---|---:|---:|---:|---:|---:|
| chematic | Rust file-backed `BufRead` V2000 MOL reader | 40 | 0 | 12,660 | 0.002981 | 13,417.57 |
| RDKit 2025.09.3 | Python `MolFromMolBlock` after block splitting | 40 | 0 | 12,660 | 0.348670 | 114.72 |

RDKit receives the same two MOL blocks extracted from the source file, while
chematic reads the file-backed stream. The figures therefore do not establish
same-process parity or a throughput ranking; the contract covers record,
failure, and source-byte accounting only.

Reproduce:

```text
python3 scripts/check_streaming_cross_engine.py --format mol --repeats 20
```
