# SDF fast-path benchmark — 2026-09-04

The Python `SDMolSupplier` now uses chematic's lightweight graph/property
reader. The writer benchmark uses `SDWriter(..., compute2d=False)` and removes
RDKit conformers, so this is a serialization-only comparison. It must not be
read as a comparison of automatic 2D depiction.

Environment: chematic 1.0.1, RDKit 2025.09.3, Python 3.13.6, macOS 26.5.2
arm64; corpus is RDKit's bundled `egfr.sdf` (365 records).

| Operation | chematic | RDKit | Result |
|---|---:|---:|---:|
| SDMolSupplier read | 323.73 µs/mol | 103.45 µs/mol | RDKit 3.1× faster |
| serialization-only write | 9.91 µs/mol | 84.39 µs/mol | chematic 8.6× faster |

The full layout-enabled `SDWriter` remains a separate path and is not claimed
to beat RDKit. Reproduce with:

```bash
TMPDIR=/private/tmp/chematic-bench python3 scripts/bench_sdf.py --rdkit --json
```

