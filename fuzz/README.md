# chematic fuzz targets

These targets use `cargo-fuzz` to exercise the public parsers. The XYZ,
MolJSON, and formats targets cap input size; the SMILES target also exercises
canonicalization.

```sh
cargo +nightly fuzz run smiles
cargo +nightly fuzz run xyz
cargo +nightly fuzz run moljson
cargo +nightly fuzz run formats
```

For a bounded campaign, use a fixed run count and minimize each corpus after
the run:

```sh
cargo +nightly fuzz run smiles -- -runs=100000
cargo +nightly fuzz run xyz -- -runs=100000
cargo +nightly fuzz run moljson -- -runs=100000
cargo +nightly fuzz run formats -- -runs=100000
cargo +nightly fuzz cmin smiles
cargo +nightly fuzz cmin xyz
cargo +nightly fuzz cmin moljson
cargo +nightly fuzz cmin formats
```

Any file created under `fuzz/artifacts/` is a failure candidate and must be
triaged before release. Keep minimized, reproducible regressions in the
corresponding `fuzz/corpus/` directory.

The six historical `formats` crash inputs from v1.0.7 were replayed against
the current target on 2026-10-20 without a failure and moved to its corpus.
The historical `smiles` crash input was already present in its corpus. The
remaining `smiles` timeout input exits successfully, but still takes several
seconds on a local arm64 build; it is retained for performance triage.

The fuzz targets have their own workspace and a committed `fuzz/Cargo.lock`
to pin their dependencies. `fuzz/target/` is generated build output.
