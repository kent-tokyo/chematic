# chematic fuzz targets

These targets use `cargo-fuzz` and keep parser resource limits explicit so a
fuzz run exercises the same bounded contracts used by the public bindings.

```sh
cargo fuzz run smiles
cargo fuzz run xyz
cargo fuzz run moljson
cargo fuzz run formats
```

For the v1.0 gate, use bounded campaigns and minimize each corpus after the
run:

```sh
cargo fuzz run smiles -- -runs=100000
cargo fuzz run xyz -- -runs=100000
cargo fuzz run moljson -- -runs=100000
cargo fuzz run formats -- -runs=100000
cargo fuzz cmin smiles
cargo fuzz cmin xyz
cargo fuzz cmin moljson
cargo fuzz cmin formats
```

Any file created under `fuzz/artifacts/` is a failure candidate and must be
triaged before release. Keep minimized, reproducible regressions in the
corresponding `fuzz/corpus/` directory.

The targets are intentionally independent of the main workspace. Their
`Cargo.lock` is generated locally by `cargo fuzz` and is not committed.
