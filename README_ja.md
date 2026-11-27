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

chematicは、分子構造の検索・解析・描画に使えるRust製のケモインフォマティクスライブラリです。
Python、Rust、Node.js、ブラウザから同じRustコアを利用できます。

[公式サイト](https://chematic.io/) · [ドキュメント](https://kent-tokyo.github.io/chematic/) · [Playground](https://kent-tokyo.github.io/chematic/playground/)

## インストール

```bash
pip install chematic
cargo add chematic --features "smiles,perception,chem,3d,fp"
npm install @kent-tokyo/chematic
```

Python wheelはC/C++コンパイラを必要とせず、ブラウザではWebAssemblyを使って端末内で分子を処理できます。

### v1.2.5 の対応範囲

現在のリリースは **v1.2.5** です。[変更履歴](https://github.com/kent-tokyo/chematic/releases/tag/v1.2.5)と[検証結果](docs/reference/validation.md)を確認できます。次版の変更は[未公開の項目](CHANGELOG.md#unreleased)にまとめています。

## 分子の検索と解析

- SMILESやSDFから分子を読み込み、分子量、LogP、トポロジカル極性表面積（TPSA）を計算する。
- SMARTSで指定した部分構造を探す。構造が似た分子を探すには、MorganフィンガープリントとTanimoto係数を使う。
- 複数の分子から最大共通部分構造（MCS）を求める。
- 分子構造をSVGで描画し、ブラウザやレポートに表示する。

chematic 1.2.2とRDKit 2026.9.2の公開Pythonパッケージで、同じ1,000分子を比較しました。
全件で環数とMorganフィンガープリント（半径2、2,048ビット）が一致し、TPSAもすべて所定の許容誤差内でした。[比較条件と測定結果](benchmarks/2026-11-21-v122-published-comparison.md)

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
cargo test --workspace --all-targets --locked --exclude chematic-py
cargo test -p chematic-py --lib --no-default-features --locked
cargo clippy --workspace --all-targets --locked -- -D warnings
```

## ライセンス

Apache-2.0またはMIT。帰属表示は[NOTICE](NOTICE)を参照してください。
