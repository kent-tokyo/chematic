# Orthorhombic periodic-neighbor cell-list parity — 2026-09-09

The periodic neighbor search now has a bounded cell-list path for
orthorhombic lattices. It indexes replicated site images in Cartesian cells,
queries a bounded 125-cell neighborhood with two-cell floating-point boundary
padding, and applies the exact `distance <= cutoff`
predicate, excludes only the trivial self image, and sorts by the public
`(center_index, neighbor_index, image)` key.

The regression fixture covers three sites, negative image shifts, cross-site
images, and a non-cubic orthorhombic cell, and compares the complete result set
with the existing reciprocal-vector bounded search. Result keys match exactly;
floating-point displacement and distance fields match within `1e-12`.

```text
CARGO_TARGET_DIR=/private/tmp/chematic-target-crystal-cell \
  cargo test -p chematic-crystal --offline --lib --tests
```

The command passed 95 unit tests, 4 integration tests, and 9 additional test
targets. This is parity and bounded-allocation evidence, not a throughput
claim. Triclinic cells and oversized replicated indexes deliberately retain
the existing validated fallback.
