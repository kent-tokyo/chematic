# Same-input CML streaming contract — 2026-09-08

This v1.0.9 development-candidate check reads the checked-in one-record CML
fixture 20 times. It requires chematic and Open Babel to report the same valid
record count, zero failures, and identical source-byte accounting.

Fixture: `benchmarks/fixtures/ethanol.cml`, 383 bytes,
SHA-256 `0d89c3e1250890147debd2ae097b27f8128f73b8d44bb60cbe3279242053e06c`.

| Engine | Boundary | Records | Failures | Input bytes |
|---|---|---:|---:|---:|
| chematic | Rust CML materialized one-shot parser | 20 | 0 | 7,660 |
| Open Babel 3.2.1 | CLI conversion per repetition | 20 | 0 | 7,660 |

CML has no equivalent RDKit block lane in this benchmark. The process and
parser boundaries differ, so this is not a same-process parity or throughput
claim; it covers record, failure, and source-byte accounting only.

Reproduce:

```text
python3 scripts/check_streaming_cml_openbabel.py --repeats 20
```
