# chematic

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![License](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)

Cheminformatics for Python, Rust, and the browser. chematic is built in pure
Rust, with bounded input handling, typed errors, and optional native InChI.

[日本語](README_ja.md) · [中文](README_zh.md) · [Documentation](https://kent-tokyo.github.io/chematic/) · [Live demo](https://kent-tokyo.github.io/chematic/playground/)

## Install

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

Python needs no C/C++ compiler. Rust and WebAssembly builds use the same core.

## Python

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
print(mol.mw, mol.logp, mol.tpsa)
print(mol.has_substructure("[OH]"))

report = chematic.report([mol], names=["aspirin"])
report.save("report.html")
```

For a small RDKit-compatible subset:

```python
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import Descriptors

mol = Chem.MolFromSmiles("CCO")
print(Descriptors.MolWt(mol))
```

This is not a full RDKit clone. Unsupported options fail explicitly; see the
[migration guide](docs/rdkit-migration.md).

## Rust

```rust
use chematic::smiles::parse;

let mol = parse("c1ccccc1")?;
println!("{}", mol.atom_count());
# Ok::<(), Box<dyn std::error::Error>>(())
```

The workspace also contains focused crates for SMILES/SMARTS, descriptors,
fingerprints, reactions, SDF/MOL/CDXML, 2D/3D, crystal formats, MCP, and
bindings. See the [format matrix](docs/format-capabilities.md).

## JavaScript / WebAssembly

```js
import init, { parse_smiles, get_descriptors_json } from "@kent-tokyo/chematic";

await init();
const mol = parse_smiles("CCO");
console.log(mol.atom_count());
console.log(JSON.parse(get_descriptors_json(mol)));
```

The WASM package supports molecule parsing, descriptors, fingerprints,
reactions, selected 2D/3D operations, and format conversion. The current
export surface is documented in [the WASM README](crates/chematic-wasm/README.md).

## What is stable

- SMILES and SMARTS parsing/writing
- Descriptors, fingerprints, Tanimoto similarity, and substructure search
- SDF/MOL V2000/V3000 and selected chemical formats
- Python, Node/WASM, and Rust bindings
- Bounded parsing, typed failures, and deterministic batch APIs

Experimental or intentionally bounded areas include 3D generation, pKa/ADMET
screening, IUPAC names, Markush/polymer expansion, CDXML editing, and the
RDKit-compatible subset. Canonical SMILES is not always a safe deduplication
key; use the fail-closed `canonical_smiles_stable_key()` API where required.

See [compatibility scope](docs/compatibility-scope.md), [validation](docs/validation.md),
and [error and resource limits](docs/error-and-limits.md) for exact guarantees.

## For research software users

Use the [researcher guide](docs/researchers.md) for a short path from
installation to a reproducible result. It explains which APIs are stable,
which comparisons are version- and corpus-pinned, how unsupported cases are
reported, and how to cite a specific release. The [benchmark index](docs/benchmark.md)
keeps performance claims separate from correctness and compatibility evidence.

## MCP server

`chematic-mcp` provides local chemistry tools over stdio for MCP-compatible
agents. It performs no network access except for the optional name lookup
tool. See [its README](crates/chematic-mcp/README.md).

## Development

```bash
cargo build --workspace
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

Benchmark methodology and dated results are kept in [docs/benchmark.md](docs/benchmark.md)
and [benchmarks/](benchmarks/). Release history is in [CHANGELOG.md](CHANGELOG.md).

## License

Licensed under either Apache License 2.0 or MIT, at your option. See
[NOTICE](NOTICE) for attribution details.
