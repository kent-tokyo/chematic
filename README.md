# chematic

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg?branch=main&event=push)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/github/deployments/kent-tokyo/chematic/github-pages?label=docs)](https://kent-tokyo.github.io/chematic/)
[![Coverage](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml/badge.svg?branch=main&event=push)](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/kent-tokyo/chematic/branch/main/graph/badge.svg)](https://app.codecov.io/github/kent-tokyo/chematic)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)
[![Ask DeepWiki](docs/assets/ask-deepwiki.svg)](https://deepwiki.com/kent-tokyo/chematic)

chematic is a pure-Rust cheminformatics library for searching, analyzing, and drawing molecular structures.
Python, Rust, Node.js, and browser applications share the same Rust core.

[日本語](README_ja.md) · [中文](README_zh.md) · [Documentation](https://kent-tokyo.github.io/chematic/) · [Playground](https://kent-tokyo.github.io/chematic/playground/)

## Install

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

The Python wheels install without a C/C++ compiler. In the browser, WebAssembly processes molecules on the user's device.

### v1.2.5 release boundary

The current release is **v1.2.5**. See [release notes](https://github.com/kent-tokyo/chematic/releases/tag/v1.2.5) and [artifact verification](docs/reference/validation.md). Upcoming changes are listed under [Unreleased](CHANGELOG.md#unreleased).

## Molecular search and analysis

- Read molecules from SMILES or SDF and calculate molecular weight, LogP, and topological polar surface area (TPSA).
- Search for a substructure with SMARTS, or find similar molecules using Morgan fingerprints and Tanimoto similarity.
- Find the maximum common substructure (MCS) of a set of molecules.
- Draw molecular structures as SVG for web pages and reports.

In a comparison of the published Python packages on the same 1,000 molecules,
chematic 1.2.2 and RDKit 2026.9.2 agreed on ring counts and Morgan fingerprints
(radius 2, 2,048 bits) for every molecule. All TPSA values agreed within the
specified tolerance. [Comparison methods and results](benchmarks/2026-11-21-v122-published-comparison.md)

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

- [When to choose chematic](docs/reference/choosing-chematic.md)
- [Getting started and cookbook](https://kent-tokyo.github.io/chematic/)
- [Compatibility scope](docs/reference/compatibility-scope.md) and [RDKit migration](docs/reference/rdkit-migration.md)
- [Validation](docs/reference/validation.md) and [benchmark methodology](docs/reference/benchmark.md)
- [Formats and bindings](docs/reference/format-capabilities.md)
- [MCP server](crates/chematic-mcp/README.md)

## Development

```bash
cargo test --workspace --all-targets --locked --exclude chematic-py
cargo test -p chematic-py --lib --no-default-features --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

## License

Licensed under Apache-2.0 or MIT. See [NOTICE](NOTICE) for attribution.
