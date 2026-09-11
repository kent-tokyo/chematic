# chematic

[English](README.md) | [中文](README_zh.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![ライセンス](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)

Python、Rust、ブラウザで使えるケモインフォマティクスライブラリです。
コアはPure Rustで実装し、入力上限・型付きエラー・決定的なバッチAPIを備えています。

[公式サイト](https://chematic.io/) · [ドキュメント](https://kent-tokyo.github.io/chematic/) · [ライブデモ](https://kent-tokyo.github.io/chematic/playground/)

## インストール

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

Python版はC/C++コンパイラなしで導入できます。

## Python

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")
print(mol.mw, mol.logp, mol.tpsa)
print(mol.has_substructure("[OH]"))

report = chematic.report([mol], names=["aspirin"])
report.save("report.html")
```

RDKit互換subsetも利用できます。

```python
from chematic import rdkit_compat as Chem
from chematic.rdkit_compat import Descriptors

mol = Chem.MolFromSmiles("CCO")
print(Descriptors.MolWt(mol))
```

完全なRDKit互換ではありません。非対応オプションは明示的にエラーになります。
詳細は[RDKit移行ガイド](docs/rdkit-migration.md)を参照してください。

## Rust

```rust
use chematic::smiles::parse;

let mol = parse("c1ccccc1")?;
println!("{}", mol.atom_count());
# Ok::<(), Box<dyn std::error::Error>>(())
```

SMILES/SMARTS、記述子、フィンガープリント、反応、SDF/MOL/CDXML、2D/3D、
結晶形式などの機能を個別crateとして提供しています。
[形式対応表](docs/format-capabilities.md)も参照してください。

## JavaScript / WebAssembly

```js
import init, { parse_smiles, get_descriptors_json } from "@kent-tokyo/chematic";

await init();
const mol = parse_smiles("CCO");
console.log(mol.atom_count());
console.log(JSON.parse(get_descriptors_json(mol)));
```

分子解析、記述子、フィンガープリント、反応、2D/3D処理、形式変換を利用できます。
詳細は[WASM README](crates/chematic-wasm/README.md)を参照してください。

## 安定性と制限

安定した機能はSMILES/SMARTS、記述子、フィンガープリント、類似度・部分構造検索、
SDF/MOL入出力、Rust/Python/Node/WASMバインディングです。

3D、pKa/ADMET、IUPAC名、Markush/polymer展開、CDXML編集、RDKit互換subsetは
実験的またはbounded subsetです。`canonical_smiles()`は常にdedup/cache keyに
使えるとは限らないため、必要な場合はfail-closedな
`canonical_smiles_stable_key()`を使用してください。

正確な契約は[互換性範囲](docs/compatibility-scope.md)、[検証](docs/validation.md)、
[エラーとリソース制限](docs/error-and-limits.md)を参照してください。

## MCPサーバー

`chematic-mcp`はMCP対応エージェント向けのローカルstdioサーバーです。
詳細は[chematic-mcp README](crates/chematic-mcp/README.md)を参照してください。

## 開発

```bash
cargo build --workspace
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

ベンチマークは[benchmark guide](docs/benchmark.md)、変更履歴は[CHANGELOG](CHANGELOG.md)にあります。

## ライセンス

Apache License 2.0 または MIT License のいずれかで利用できます。帰属表示は
[NOTICE](NOTICE)を参照してください。
