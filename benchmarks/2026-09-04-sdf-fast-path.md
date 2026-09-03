# SDF fast-path benchmark — 2026-09-04

The Python `SDMolSupplier` now uses chematic's lightweight graph/property
reader. The writer benchmark uses `SDWriter(..., compute2d=False)` and removes
RDKit conformers, so this is a serialization-only comparison. It must not be
read as a comparison of automatic 2D depiction.

Environment: code included in chematic 1.0.2, measured before the version-only
metadata bump when the local wheel still reported 1.0.1; source optimization
through `41acb3b4`, RDKit 2025.09.3, Python 3.13.6, macOS 26.5.2 arm64; corpus
is RDKit's bundled `egfr.sdf` (365 records).

## Follow-up optimization

Seven fresh processes were measured before and after the follow-up. Medians
exclude neither outliers nor warm-up processes. Both implementations parsed or
wrote all 365 records in every run.

| Operation | Before | After | Improvement | RDKit after |
|---|---:|---:|---:|---:|
| SDMolSupplier read | 11.95 µs/mol | **9.48 µs/mol** | **1.26×** | 99.96 µs/mol |
| serialization-only write | 10.16 µs/mol | **7.62 µs/mol** | **1.33×** | 79.54 µs/mol |

The after medians are 10.5× and 10.4× faster than RDKit respectively on this
machine and corpus. The read path retains strict malformed/non-finite Z
validation; only X/Y conversion and coordinate storage that the graph-only
supplier cannot expose are skipped. The writer result applies only when 2D
layout is disabled.

## Initial fast-path snapshot

| Operation | chematic | RDKit | Result |
|---|---:|---:|---:|
| SDMolSupplier read | 16.17 µs/mol | 110.8 µs/mol | chematic 6.8× faster |
| serialization-only write | 10.34 µs/mol | 83.01 µs/mol | chematic 8.0× faster |

The full layout-enabled `SDWriter` remains a separate path and is not claimed
to beat RDKit. Reproduce with:

```bash
for i in 1 2 3 4 5 6 7; do
  TMPDIR=/private/tmp/chematic-bench \
    python3 scripts/bench_sdf.py --rdkit --json
done
```
