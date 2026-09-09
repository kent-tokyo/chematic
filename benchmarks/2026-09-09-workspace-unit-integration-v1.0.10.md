# Workspace unit/integration gate — v1.0.10

The current working tree passed both workspace Rust test scopes:

```sh
CARGO_TARGET_DIR=/private/tmp/chematic-target-workspace-check cargo test --workspace --lib --quiet
CARGO_TARGET_DIR=/private/tmp/chematic-target-workspace-check cargo test --workspace --tests --quiet
```

All executed library and integration suites passed. Tests marked ignored remained ignored and are not promoted to executed coverage by this record.

This is local Rust evidence only; Python packaging, browser/WASM runtime parity, external review, and publication remain separate gates.
