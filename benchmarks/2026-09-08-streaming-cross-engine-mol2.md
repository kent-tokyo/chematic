# Same-input MOL2 streaming contract — 2026-09-08

This v1.0.9 development-candidate check reads the checked-in one-record MOL2
fixture 20 times. It requires chematic and RDKit to report the same valid
record count, zero failures, and identical source-byte accounting.

Fixture: `benchmarks/fixtures/ethanol.mol2`, 361 bytes,
SHA-256 `80144aad97332e2f58eb42bb2da1472587a9dfc2466f9d2fe081eb9834c2224d`.

| Engine | Boundary | Records | Failures | Input bytes | Seconds | Records/s |
|---|---|---:|---:|---:|---:|---:|
| chematic | Rust MOL2 materialized one-shot parser | 20 | 0 | 7,220 | 0.003190 | 6,270.16 |
| RDKit 2025.09.3 | Python `MolFromMol2Block` block parser | 20 | 0 | 7,220 | 0.399155 | 50.11 |

The parser boundaries differ, and chematic's current MOL2 path is materialized
rather than a true streaming reader. These measurements do not establish
same-process parity or a throughput ranking; the contract covers record,
failure, and source-byte accounting only.

Reproduce:

```text
python3 scripts/check_streaming_cross_engine.py --format mol2 --repeats 20
```
