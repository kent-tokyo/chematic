# Gzip Extended XYZ streaming contract — 2026-09-08

This v1.0.9 development-candidate check uses a deterministic gzip encoding of
the checked-in two-frame Extended XYZ fixture and runs it 20 times. Both
engines report 40 records and zero failures.

The fixture is 268 decompressed bytes and 113 compressed bytes per repetition;
the decompressed payload SHA-256 is
`0094ab92ea84aa6acd826004e1849f2f2cb750c81d02e061c4156dc25c2830c7`.

| Engine | Boundary | Records | Failures | Reported input bytes |
|---|---|---:|---:|---:|
| chematic | Rust gzip file-backed `BufRead` Extended XYZ reader | 40 | 0 | 2,260 compressed |
| RDKit 2025.09.3 | Python XYZ block parser over decompressed Extended XYZ frames | 40 | 0 | 5,360 decompressed |

The byte counters intentionally describe different stages. This is a
record/failure contract, not same-process parity or a throughput comparison.

Reproduce:

```text
python3 scripts/check_streaming_gzip_contract.py --format extxyz --repeats 20
```
