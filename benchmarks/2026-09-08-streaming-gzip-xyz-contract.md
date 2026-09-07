# Gzip XYZ streaming contract — 2026-09-08

This v1.0.9 development-candidate check uses a deterministic gzip encoding of
the checked-in two-frame XYZ fixture and runs it 20 times. Both engines report
40 records and zero failures.

The fixture is 161 decompressed bytes and 65 compressed bytes per repetition;
the decompressed payload SHA-256 is
`3c717e4520a3e1b7b90ace86b72e73cce4dda38de7ae29aa811733d2ec227221`.

| Engine | Boundary | Records | Failures | Reported input bytes |
|---|---|---:|---:|---:|
| chematic | Rust gzip file-backed `BufRead` XYZ reader | 40 | 0 | 1,300 compressed |
| RDKit 2025.09.3 | Python XYZ block parser over decompressed frames | 40 | 0 | 3,220 decompressed |

The byte counters intentionally describe different stages. This is a
record/failure contract, not same-process parity or a throughput comparison.

Reproduce:

```text
python3 scripts/check_streaming_gzip_contract.py --format xyz --repeats 20
```
