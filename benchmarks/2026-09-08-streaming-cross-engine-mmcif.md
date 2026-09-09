# Same-input mmCIF streaming contract — 2026-09-08

This v1.0.9 development-candidate check reads the checked-in one-record mmCIF
fixture 20 times. It requires chematic and Open Babel to report the same valid
record count, zero failures, and identical source-byte accounting.

Fixture: `benchmarks/fixtures/minimal.mmcif`, 698 bytes,
SHA-256 `d583d8047758d0e804c98729796d76eea3fa87f3b5b8ed85962717e784e98998`.

| Engine | Boundary | Records | Failures | Input bytes |
|---|---|---:|---:|---:|
| chematic | Rust mmCIF materialized one-shot parser | 20 | 0 | 13,960 |
| Open Babel 3.2.1 | CLI conversion per repetition | 20 | 0 | 13,960 |

The process and parser boundaries differ, so this is not a same-process parity
or throughput claim; it covers record, failure, and source-byte accounting
only.

Reproduce:

```text
python3 scripts/check_streaming_cml_openbabel.py --format mmcif --repeats 20
```
