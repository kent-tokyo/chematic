# Same-input CDXML streaming contract — 2026-09-08

This v1.0.9 development-candidate check reads the checked-in one-record CDXML
fixture 20 times. It requires chematic and Open Babel to report the same valid
record count, zero failures, and identical source-byte accounting.

Fixture: `benchmarks/fixtures/ethanol.cdxml`, 298 bytes,
SHA-256 `788130936f9a9cf46bb2c64d54e05bb9da39a707e04199f8ec4ae97416436857`.

| Engine | Boundary | Records | Failures | Input bytes |
|---|---|---:|---:|---:|
| chematic | Rust CDXML materialized one-shot parser | 20 | 0 | 5,960 |
| Open Babel 3.2.1 | CLI conversion per repetition | 20 | 0 | 5,960 |

The process and parser boundaries differ, so this is not a same-process parity
or throughput claim; it covers record, failure, and source-byte accounting
only.

Reproduce:

```text
python3 scripts/check_streaming_cml_openbabel.py --format cdxml --repeats 20
```
