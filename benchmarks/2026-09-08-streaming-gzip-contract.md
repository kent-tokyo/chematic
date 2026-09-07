# Gzip SDF streaming contract — 2026-09-08

This v1.0.9 development-candidate check uses a deterministic gzip encoding of
the checked-in two-record SDF fixture and runs it 20 times. Both engines report
40 records and zero failures.

The fixture is 633 decompressed bytes and 134 compressed bytes per repetition;
the decompressed payload SHA-256 is
`627e7814694f7d00a9212712ec697f2174d30327b459d983d2a10f717cae0465`.

| Engine | Boundary | Records | Failures | Reported input bytes |
|---|---|---:|---:|---:|
| chematic | Rust gzip file-backed `BufRead` reader | 40 | 0 | 2,680 compressed |
| RDKit 2025.09.3 | Python block parser over decompressed fixture | 40 | 0 | 12,660 decompressed |

The byte counters intentionally describe different stages: compressed bytes
for chematic's gzip stream and decompressed bytes for RDKit's block input.
This is a record/failure contract, not same-process parity or a throughput
comparison.

Reproduce:

```text
python3 scripts/check_streaming_gzip_contract.py --repeats 20
```
