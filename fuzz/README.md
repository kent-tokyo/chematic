# chematic fuzz targets

These targets use `cargo-fuzz` and keep parser resource limits explicit so a
fuzz run exercises the same bounded contracts used by the public bindings.

```sh
cargo fuzz run smiles
cargo fuzz run xyz
cargo fuzz run moljson
```

The targets are intentionally independent of the main workspace. Their
`Cargo.lock` is generated locally by `cargo fuzz` and is not committed.
