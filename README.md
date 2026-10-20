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

Pure-Rust cheminformatics for Python, Rust, and the browser when molecular data
must stay local and unsupported chemistry must fail explicitly.

[日本語](README_ja.md) · [中文](README_zh.md) · [Documentation](https://kent-tokyo.github.io/chematic/) · [Playground](https://kent-tokyo.github.io/chematic/playground/)

## Install

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

Python wheels need no C/C++ compiler. All bindings share the Rust core.

### v1.2.0 release boundary

v1.2.0 adds checked XBS and TURBOMOLE/POSCAR input paths, typed Pattern
fingerprints, bounded 3D geometry diagnostics, and stronger fail-closed
canonical-SMILES, reaction, R-group, and standardization handling. Internal
chemistry, I/O, CLI, and evidence modules were split without intentionally
breaking the public 1.x API. Exact performance and compatibility claims remain
tied to their dated artifacts; source improvements are not relabelled as
published-package evidence. See
[validation](docs/validation.md), [benchmarks](benchmarks/README.md), and the
[CHANGELOG](CHANGELOG.md).

## Choose chematic when

- chemistry must run locally in a browser or inside a Rust application;
- Python, Rust, Node, and WebAssembly should share one implementation;
- rejected inputs, resource limits, and lossy conversions must stay visible;
- adoption claims must point to versioned validation and benchmark records.

Choose RDKit for maximum ecosystem and chemistry breadth, and Open Babel for
broad command-line format conversion. Compare pure-Rust alternatives on your
own corpus and exact operation. See [When to choose chematic](docs/choosing-chematic.md)
for the decision guide and current limits.

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

- [When to choose chematic](docs/choosing-chematic.md)
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
