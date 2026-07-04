# Contributing to chematic

Thank you for your interest in contributing! chematic is a pure-Rust cheminformatics library — contributions of all sizes are welcome.

## Quick start

```bash
git clone https://github.com/kent-tokyo/chematic.git
cd chematic
cargo test --workspace --lib --quiet   # 211 tests, all should pass
cargo clippy --workspace -- -D warnings
```

No C compiler, no Python, no external tools required for the core crates.

## Where to contribute

| Area | Crate | Difficulty |
|------|-------|-----------|
| Bug fixes anywhere | any | Easy |
| New molecular descriptors | `chematic-chem` | Easy — add a `pub fn` + tests |
| SMARTS pattern fixes | `chematic-smarts` | Medium |
| New file formats | `chematic-mol` | Medium |
| Force field parameters | `chematic-ff` | Medium–Hard |
| Kekulization / aromaticity | `chematic-core` | Hard |

**Easiest entry point**: `chematic-core` has zero external dependencies (no `[dependencies]` in `Cargo.toml`) and is pure safe Rust. Adding a helper to `molecule.rs` or a new descriptor to `chematic-chem/src/descriptors.rs` is a great first contribution.

## Branch naming

`feat/*` new feature · `fix/*` bug fix · `docs/*` docs only · `release/*` version bump + CHANGELOG + tag

## Main branch protection (recommended settings)

GitHub → Settings → Branches → Add rule for `main`:

- [x] Require a pull request before merging
- [x] Require status checks: `Test`, `Clippy`, `Format Check`, `Build Check`
- [x] Do not allow bypassing the above settings

## Development workflow

```bash
# Run tests for a single crate
cargo test -p chematic-chem --lib --quiet

# Run a specific test by name
cargo test -p chematic-chem --lib -- hba_count

# Check native InChI feature (requires C compiler)
cargo test -p chematic-inchi --features native-inchi --lib --quiet
```

## Code style

- **No unsafe code** — `#![forbid(unsafe_code)]` is enforced in all crates, except the `native-inchi` feature's C FFI boundary (`crates/chematic-inchi/src/native/`).
- **No external C/C++ deps** by default (the `native-inchi` feature is the only exception).
- **WASM-compatible** — core crates must not use `std::fs`, threads, or OS-specific APIs.
- **Justify unsafe/unwrap/expect/panic!** outside `#[cfg(test)]` with a `// safe: ...` comment explaining why it cannot fail (e.g. `// safe: index bounds-checked above`).
- Format: `cargo fmt --workspace` before committing.
- Lints: `cargo clippy --workspace -- -D warnings` must pass (zero warnings).

## Adding a new descriptor

1. Add a `pub fn my_descriptor(mol: &Molecule) -> f64` in `crates/chematic-chem/src/descriptors.rs`
2. Export it in `crates/chematic-chem/src/lib.rs`
3. Add at least one `#[test]` with a known value (e.g. benzene, aspirin)
4. If it belongs in the Python API, add it to `Mol` in `crates/chematic-py/src/lib.rs`

## Pull requests

- Keep PRs focused — one feature or fix per PR.
- Include tests. New functions without tests will be asked to add them.
- Update `CHANGELOG.md` under `## [Unreleased]` with a brief description.
- CI must pass (`cargo test --workspace --lib --quiet` + `cargo clippy --workspace -- -D warnings`).

## Reporting issues

Use the [Bug report](.github/ISSUE_TEMPLATE/bug_report.md) template. A minimal SMILES reproduction is the most helpful thing you can provide.

## Questions

Open a [Discussion](https://github.com/kent-tokyo/chematic/discussions) — no need to file an issue for questions.
