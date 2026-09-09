# Reaction and 3D focused gate — v1.0.10

The local Rust regression gate passed:

- `chematic-rxn`: 209/209 tests
- `chematic-3d`: 596/596 executed tests, with 10 explicitly ignored experimental long-run tests

Reproduce with:

```sh
CARGO_TARGET_DIR=/private/tmp/chematic-target-roadmap-focus cargo test -p chematic-rxn --offline --lib
CARGO_TARGET_DIR=/private/tmp/chematic-target-roadmap-focus cargo test -p chematic-3d --offline --lib
```

The reaction lane covers typed documents, parsing, query matching, bounded
transformation, stereo preservation, retro templates, and invalid-product
handling. The 3D lane covers geometry, stereo constraints, force-field
diagnostics, ensemble policies, symmetry-aware metrics, and format contracts.

This is a local regression gate only. It does not replace curated precision /
recall reports, independent force-field oracles, full RDKit TFD parity, or
external and browser evidence.
