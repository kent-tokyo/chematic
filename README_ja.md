# chematic

[English](README.md) | [中文](README_zh.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg?branch=main&event=push)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![Docs](https://img.shields.io/github/deployments/kent-tokyo/chematic/github-pages?label=docs)](https://kent-tokyo.github.io/chematic/)
[![Coverage](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml/badge.svg?branch=main&event=push)](https://github.com/kent-tokyo/chematic/actions/workflows/coverage.yml)
[![codecov](https://codecov.io/gh/kent-tokyo/chematic/branch/main/graph/badge.svg)](https://app.codecov.io/github/kent-tokyo/chematic)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![ライセンス](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)
[![Ask DeepWiki](docs/assets/ask-deepwiki.svg)](https://deepwiki.com/kent-tokyo/chematic)

分子データを外部へ送らず、未対応の化学表現を曖昧に処理したくない場面に向けた、
Rust製のケモインフォマティクスライブラリです。Python、Rust、Node、ブラウザから
同じコアを利用できます。

[公式サイト](https://chematic.io/) · [ドキュメント](https://kent-tokyo.github.io/chematic/) · [Playground](https://kent-tokyo.github.io/chematic/playground/)

## インストール

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

Python wheelはC/C++コンパイラを必要とせず、各バインディングは同じRustコアを使います。

### v1.2.0 の対応範囲

現在のリリースは **v1.2.0** です。変更内容は[CHANGELOG](CHANGELOG.md)を参照して
ください。精度・性能の主張は、[検証報告](docs/reference/validation.md)と
[ベンチマーク索引](benchmarks/README.md)に記載した成果物、コーパス、実行環境の
範囲に限ります。

## chematicが向いている場面

- ブラウザ内またはRustアプリ内で化学処理を完結させたい
- Python・Rust・Node・WebAssemblyで同じ実装を使いたい
- 入力拒否、リソース上限、変換時の情報損失を結果として扱いたい
- バージョン付きの検証結果を確認してから採用したい

機能範囲とエコシステムの広さを優先するならRDKit、幅広いコマンドライン形式変換を
優先するならOpen Babelが第一候補です。選定条件と未対応範囲は
[選定ガイド](docs/reference/choosing-chematic.md)にまとめています。

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

`rdkit_compat`の互換範囲は文書化されています。未対応オプションは近似せず、
明示的に失敗します。

## 範囲

SMILES/SMARTS、記述子、フィンガープリント、類似度・部分構造検索、主要なMOL/SDF I/O、
文書化済みのRust/Python/Node/WASMバインディングが安定対象です。3D、pKa/ADMET、
IUPAC、Markush/polymer、CDXML編集、RDKit型APIは実験的または制限付きです。

`canonical_smiles()`は汎用の同一性キーではありません。必要な場合は、文書化された範囲で
安全側に停止する`canonical_smiles_stable_key()`を使用してください。

## ガイド

- [選定ガイド](docs/reference/choosing-chematic.md)
- [導入とCookbook](https://kent-tokyo.github.io/chematic/)
- [互換性範囲](docs/reference/compatibility-scope.md)・[RDKit移行](docs/reference/rdkit-migration.md)
- [検証](docs/reference/validation.md)・[ベンチマーク方法](docs/reference/benchmark.md)
- [形式とバインディング](docs/reference/format-capabilities.md)
- [MCPサーバー](crates/chematic-mcp/README.md)

## 開発

```bash
cargo test --workspace --all-targets --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

## ライセンス

Apache-2.0またはMIT。帰属表示は[NOTICE](NOTICE)を参照してください。
