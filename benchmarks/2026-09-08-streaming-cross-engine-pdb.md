# Same-input PDB streaming contract — 2026-09-08

This v1.0.9 development-candidate check reads the checked-in one-record PDB
fixture 20 times. It requires chematic and Open Babel to report the same valid
record count, zero failures, and identical source-byte accounting.

Fixture: `benchmarks/fixtures/minimal.pdb`, 286 bytes,
SHA-256 `69d1e5a2d6cbdb39377325cbbf077face011ef0893718e28b42d0563901ad15b`.

| Engine | Boundary | Records | Failures | Input bytes |
|---|---|---:|---:|---:|
| chematic | Rust PDB materialized one-shot parser | 20 | 0 | 5,720 |
| Open Babel 3.2.1 | CLI conversion per repetition | 20 | 0 | 5,720 |

The PDB reader retains its explicit lenient/line-limit boundary. The process
and parser boundaries differ, so this is not a same-process parity or
throughput claim; it covers record, failure, and source-byte accounting only.

Reproduce:

```text
python3 scripts/check_streaming_cml_openbabel.py --format pdb --repeats 20
```
