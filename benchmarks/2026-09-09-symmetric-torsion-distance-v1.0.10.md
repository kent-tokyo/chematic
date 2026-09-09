# Symmetry-aware torsion distance: v1.0.10 candidate

The 3D torsion-motif module now exposes a bounded diagnostic that minimizes
mean circular dihedral difference over topology-preserving self-mappings. This
prevents a terminal-atom relabeling from appearing as a conformational change.
The result is normalized to `[0, 1]`.

The focused local contract and a ten-molecule change-sensitivity corpus passed:

```text
running 2 tests
test torsion_motif::tests::symmetric_torsion_distance_is_zero_for_terminal_swap ... ok
test torsion_motif::tests::symmetric_torsion_distance_is_bounded_and_detects_change ... ok
test torsion_motif::tests::symmetric_torsion_distance_corpus_is_bounded_and_change_sensitive ... ok
```

Reproduce with:

```sh
cargo test -p chematic-3d --offline --lib torsion_motif::tests::symmetric_torsion_distance
```

The broader 3D crate regression also passed:

```text
cargo test -p chematic-3d --offline
594 passed; 0 failed; 10 ignored
```

The machine-readable record includes this full-crate regression separately
from the focused three-test diagnostic gate.

The corpus contains ten flexible molecules; all ten produced a valid measured
distance and exercised at least four rotatable cases. This is a local API and
invariant slice. It is intentionally not presented as RDKit TFD parity, a full
conformer-quality corpus measurement, or completion of the remaining
symmetry-aware TFD roadmap item.
