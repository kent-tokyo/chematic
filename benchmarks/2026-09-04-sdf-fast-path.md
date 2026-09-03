# SDF fast-path benchmark — 2026-09-04

The Python `SDMolSupplier` now uses chematic's lightweight graph/property
reader. The writer benchmark uses `SDWriter(..., compute2d=False)` and removes
RDKit conformers, so this is a serialization-only comparison. It must not be
read as a comparison of automatic 2D depiction.

Environment: chematic 1.0.1, RDKit 2025.09.3, Python 3.13.6, macOS 26.5.2
arm64; corpus is RDKit's bundled `egfr.sdf` (365 records).

| Operation | chematic | RDKit | Result |
|---|---:|---:|---:|
| SDMolSupplier read | 16.17 µs/mol | 110.8 µs/mol | chematic 6.8× faster |
| serialization-only write | 10.34 µs/mol | 83.01 µs/mol | chematic 8.0× faster |

The full layout-enabled `SDWriter` remains a separate path and is not claimed
to beat RDKit. Reproduce with:

```bash
TMPDIR=/private/tmp/chematic-bench python3 scripts/bench_sdf.py --rdkit --json
```
