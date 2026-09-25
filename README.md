# chematic

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)

Pure-Rust cheminformatics for Python, Rust, and the browser: local-first,
bounded, and explicit about unsupported chemistry.

[日本語](README_ja.md) · [中文](README_zh.md) · [Documentation](https://kent-tokyo.github.io/chematic/) · [Playground](https://kent-tokyo.github.io/chematic/playground/)

## Install

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

Python wheels need no C/C++ compiler. All bindings share the Rust core.

### v1.0.25 release boundary

v1.0.25 adds SMILES atom-output order plus fragment and Rust-reaction source
atom provenance, corrects the named RDKit-compatible Python HBA profile, and
fixes plain-writer ring-closure spelling. These are scoped compatibility
changes, not a claim of complete RDKit parity. See [validation](docs/validation.md)
and the [CHANGELOG](CHANGELOG.md).

## Use it

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
print(mol.mw, mol.logp, mol.tpsa)
print(mol.has_substructure("[OH]"))
```

```rust
use chematic::smiles::parse;

let mol = parse("c1ccccc1").expect("valid SMILES");
println!("{}", mol.atom_count());
```

```js
import init, { parse_smiles, get_descriptors_json } from "@kent-tokyo/chematic";

await init();
const mol = parse_smiles("CCO");
console.log(JSON.parse(get_descriptors_json(mol)));
```

The optional `rdkit_compat` Python namespace covers a documented subset;
unsupported options fail explicitly rather than silently approximating RDKit.

## Scope

Stable selected paths include SMILES/SMARTS, descriptors, fingerprints,
similarity and substructure search, common MOL/SDF I/O, and documented
Rust/Python/Node/WASM bindings. 3D, pKa/ADMET, IUPAC, Markush/polymer,
rich CDXML editing, and RDKit-style APIs remain experimental or bounded.

`canonical_smiles()` is not a universal identity key. Use the fail-closed
`canonical_smiles_stable_key()` only within its documented domain.

## Find the right guide

- [Getting started and cookbook](https://kent-tokyo.github.io/chematic/)
- [Compatibility scope](docs/compatibility-scope.md) and [RDKit migration](docs/rdkit-migration.md)
- [Validation](docs/validation.md) and [benchmark methodology](docs/benchmark.md)
- [Formats and bindings](docs/format-capabilities.md)
- [MCP server](crates/chematic-mcp/README.md)

## Development

```bash
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

## License

Licensed under Apache-2.0 or MIT. See [NOTICE](NOTICE) for attribution.
