# chematic

[English](README.md) | [中文](README_zh.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![ライセンス](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)

Python、Rust、ブラウザ向けのPure Rustケモインフォマティクスです。ローカルで動作し、
入力上限と型付きエラーにより、未対応の化学表現を曖昧に処理しません。

[公式サイト](https://chematic.io/) · [ドキュメント](https://kent-tokyo.github.io/chematic/) · [Playground](https://kent-tokyo.github.io/chematic/playground/)

## インストール

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

Python wheelはC/C++コンパイラを必要とせず、各バインディングは同じRustコアを使います。

### v1.0.25 の対応範囲

v1.0.25はSMILESの原子出力順、断片・Rust反応生成物の由来原子APIを追加し、
RDKit互換Python HBA profileとplain writerの環閉鎖表記を修正します。これは
範囲を明示した互換性改善であり、完全なRDKit互換の主張ではありません。
[検証報告](docs/validation.md)と[CHANGELOG](CHANGELOG.md)を参照してください。

## 使い方

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
console.log(JSON.parse(get_descriptors_json(parse_smiles("CCO"))));
```

`rdkit_compat`は文書化されたsubsetです。未対応オプションは近似せず明示的に失敗します。

## 範囲

SMILES/SMARTS、記述子、フィンガープリント、類似度・部分構造検索、主要なMOL/SDF I/O、
文書化済みのRust/Python/Node/WASMバインディングが安定対象です。3D、pKa/ADMET、
IUPAC、Markush/polymer、CDXML編集、RDKit型APIは実験的またはbounded subsetです。

`canonical_smiles()`は汎用の同一性キーではありません。必要な場合は、文書化された範囲で
fail-closedに動作する`canonical_smiles_stable_key()`を使用してください。

## ガイド

- [導入とCookbook](https://kent-tokyo.github.io/chematic/)
- [互換性範囲](docs/compatibility-scope.md)・[RDKit移行](docs/rdkit-migration.md)
- [検証](docs/validation.md)・[ベンチマーク方法](docs/benchmark.md)
- [形式とバインディング](docs/format-capabilities.md)
- [MCPサーバー](crates/chematic-mcp/README.md)

## 開発

```bash
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

## ライセンス

Apache-2.0またはMIT。帰属表示は[NOTICE](NOTICE)を参照してください。
