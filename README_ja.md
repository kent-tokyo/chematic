# chematic

[English](README.md) | [中文](README_zh.md)

Pure Rust 製のケモインフォマティクスライブラリ。RDKit の代替を目指す、C/C++ FFI ゼロの Rust 実装。

---

## ライブデモ

**[https://kent-tokyo.github.io/chematic/](https://kent-tokyo.github.io/chematic/)** — 記述子計算、薬らしさルール、類似度比較をブラウザ上の WebAssembly で実行できるインタラクティブデモ。

---

## 設計目標

**Pure Rust、FFI ゼロ**
rdkit-sys も openbabel バインディングも使用しない。すべてのアルゴリズムを安全な Rust で実装する。

**WASM 対応、軽量**
コアクレートは `wasm32-unknown-unknown` に無修正でコンパイルできる。バイナリサイズは数百 KB 程度（C++ FFI ラッパーの数十 MB とは対照的）。

**化学ドメイン固有の実装**
汎用グラフライブラリのラッパーではなく、ケクレ化・Hückel 芳香族性・CIP 立体化学・SSSR など、化学ドメイン固有のアルゴリズムを Rust でスクラッチ実装する。

**再現性と決定性**
フィンガープリントは固定不変量順序の FNV-1a ハッシュを使用する。同じ SMILES 入力から常に同じビット列が得られる。乱数なし、プラットフォーム依存なし。

---

## 現在のステータス

全フェーズ完了。544 テスト、全パス。

| クレート               | 説明                                                                                         | テスト数 |
|------------------------|----------------------------------------------------------------------------------------------|---------|
| `chematic-core`        | Atom, Bond, Molecule, Element, ケクレ化（依存ゼロ）                                         | 30      |
| `chematic-smiles`      | OpenSMILES パーサー、ライター、正規 SMILES                                                  | 52      |
| `chematic-perception`  | SSSR (Balducci-Pearlman)、Huckel 芳香族性認識                                               | 14      |
| `chematic-mol`         | MOL/SDF V2000+V3000 パーサーとライター                                                      | 37      |
| `chematic-depict`      | 2D SVG 描画（CPK カラー・アトム/ボンドハイライト）                                          | 15      |
| `chematic-chem`        | 記述子、BRICS フラグメント化、QED、標準化、Murcko スキャフォルド、CIP 立体化学             | 216     |
| `chematic-fp`          | ECFP4/6、MACCS 166-bit、位相的パス、AtomPair、Torsion FP、Tanimoto/Dice                    | 44      |
| `chematic-smarts`      | SMARTS（再帰・原子価・ハイブリッド化対応）、VF2 部分構造一致、MCS                          | 76      |
| `chematic-3d`          | 3D 座標生成、PDB/XYZ ファイル形式                                                           | 25      |
| `chematic-rxn`         | 反応 SMILES パーサーとライター                                                               | 15      |
| `chematic-wasm`        | WebAssembly バインディング — npm: `@kent-tokyo/chematic`                                    | 18      |
| `chematic`             | フィーチャーフラグ付きアンブレラクレート（全サブクレート）                                  | 1       |

```
cargo test --workspace   # 544 テスト、全パス
```

---

## クイックスタート

### アンブレラクレートを使う場合

```toml
# Cargo.toml
[dependencies]
chematic = { git = "https://github.com/kent-tokyo/chematic", features = ["smiles", "fp"] }
```

```rust
use chematic::smiles::{parse, canonical_smiles};
use chematic::fp::ecfp4;
```

### 個別クレートを使う場合

```toml
# Cargo.toml
[dependencies]
chematic-smiles     = { git = "https://github.com/kent-tokyo/chematic" }
chematic-perception = { git = "https://github.com/kent-tokyo/chematic" }
chematic-fp         = { git = "https://github.com/kent-tokyo/chematic" }
```

```rust
use chematic_smiles::{parse, canonical_smiles};
use chematic_perception::{find_sssr, assign_aromaticity};
use chematic_fp::{ecfp4, tanimoto_ecfp4};

fn main() {
    let benzene = parse("c1ccccc1").unwrap();
    let toluene = parse("Cc1ccccc1").unwrap();

    // 環認識と芳香族性
    let rings = find_sssr(&benzene);
    println!("環数: {}", rings.ring_count()); // 1

    // フィンガープリント類似度
    let sim = tanimoto_ecfp4(&benzene, &toluene);
    println!("Tanimoto(ベンゼン, トルエン): {sim:.3}"); // ~0.5

    // 正規 SMILES
    println!("{}", canonical_smiles(&benzene)); // c1ccccc1
}
```

---

## SMARTS 部分構造検索

```rust
use chematic_smiles::parse;
use chematic_smarts::{parse_smarts, find_matches};

let mol = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap(); // アスピリン
let query = parse_smarts("[$(C(=O)O)]").unwrap();   // カルボン酸/エステル C
let matches = find_matches(&query, &mol);
println!("C(=O)O 基の数: {}", matches.len()); // 2
```

---

## 分子記述子計算

```rust
use chematic_smiles::parse;
use chematic_chem::{molecular_weight, tpsa, logp_crippen, fsp3, qed, lipinski_passes};

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
println!("分子量:   {:.2}", molecular_weight(&aspirin)); // ~180.16
println!("TPSA:     {:.2}", tpsa(&aspirin));             // ~63.6
println!("LogP:     {:.2}", logp_crippen(&aspirin));     // ~1.2
println!("Fsp3:     {:.3}", fsp3(&aspirin));             // ~0.111
println!("QED:      {:.3}", qed(&aspirin));              // ドラッグライクネス
println!("Lipinski: {}", lipinski_passes(&aspirin));     // true
```

---

## BRICS フラグメント化

```rust
use chematic_smiles::parse;
use chematic_chem::brics_fragments;

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
let frags = brics_fragments(&aspirin);
println!("フラグメント数: {}", frags.len()); // ≥ 2
```

---

## フィンガープリント

```rust
use chematic_smiles::parse;
use chematic_fp::{ecfp4, atom_pair_fp, torsion_fp};

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
let caffeine = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();

let sim_ecfp4    = ecfp4(&aspirin).tanimoto(&ecfp4(&caffeine));
let sim_atompair = atom_pair_fp(&aspirin).tanimoto(&atom_pair_fp(&caffeine));
let sim_torsion  = torsion_fp(&aspirin).tanimoto(&torsion_fp(&caffeine));
```

---

## 2D 構造式 SVG 出力

```rust
use chematic_smiles::parse;
use chematic_depict::depict_svg;

let caffeine = parse("Cn1cnc2c1c(=O)n(c(=O)n2C)C").unwrap();
let svg = depict_svg(&caffeine);
std::fs::write("caffeine.svg", svg).unwrap();
```

### ハイライト付き描画

```rust
use std::collections::HashSet;
use chematic_smiles::parse;
use chematic_depict::depict_svg_highlighted;

let mol = parse("c1ccncc1").unwrap(); // ピリジン
let n_idx = mol.atoms().find(|(_, a)| a.element.atomic_number() == 7)
               .map(|(i, _)| i).unwrap();
let svg = depict_svg_highlighted(&mol, &HashSet::from([n_idx]), &HashSet::new());
// → N 原子が黄色の丸でハイライト、N ラベルが青（CPK カラー）
```

---

## JavaScript / TypeScript（WebAssembly）

```sh
npm install @kent-tokyo/chematic
```

```js
import init, { parse_smiles, tanimoto_ecfp4, tanimoto_atom_pair, brics_fragment_count } from '@kent-tokyo/chematic';

await init();

const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // アスピリン
console.log(mol.molecular_weight()); // ~180.16
console.log(mol.logp_crippen());     // ~1.2
console.log(mol.qed());              // ドラッグライクネス [0,1]
console.log(mol.fsp3());             // sp3 炭素割合
console.log(brics_fragment_count(mol)); // BRICS フラグメント数

const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));    // ECFP4 類似度
console.log(tanimoto_atom_pair(mol, caffeine)); // AtomPair 類似度
```

---

## 他のケモインフォマティクスライブラリとの比較

| 機能                               | chematic                      | RDKit (rdkit-sys)  | OpenBabel FFI  | chemcore / purr   |
|------------------------------------|-------------------------------|--------------------|----------------|-------------------|
| 実装言語                           | Pure Rust                     | Rust + C++ FFI     | Rust + C++ FFI | Pure Rust         |
| WASM ターゲット                    | 対応                          | 非対応             | 非対応         | 部分対応          |
| バイナリサイズ（コア）             | 約 500 KB                     | 約 50 MB           | 約 20 MB       | 約 200 KB         |
| OpenSMILES パーサー                | 完全実装                      | 完全実装           | 完全実装       | 部分実装          |
| SMILES ライター / 正規 SMILES      | 対応                          | 対応               | 対応           | 非対応            |
| ケクレ化                           | 対応                          | 対応               | 対応           | 非対応            |
| 芳香族性認識                       | 対応 (Huckel 則)              | 対応               | 対応           | 部分対応          |
| 環認識 (SSSR)                      | 対応                          | 対応               | 対応           | 非対応            |
| SDF/MOL V2000+V3000                | 対応                          | 対応               | 対応           | 非対応            |
| 2D 描画 (SVG、CPK カラー)          | 対応                          | 対応               | 対応           | 非対応            |
| ECFP フィンガープリント            | 対応 (ECFP4/6)                | 対応               | 対応           | 非対応            |
| AtomPair / Torsion FP              | 対応                          | 対応               | 対応           | 非対応            |
| MACCS フィンガープリント           | 対応 (166-bit 構造キー)       | 対応               | 対応           | 非対応            |
| SMARTS / 部分構造検索              | 対応 (VF2 + 再帰 SMARTS)      | 対応               | 対応           | 非対応            |
| 分子記述子計算                     | 対応 (MW/LogP/TPSA/Fsp3/QED/…) | 対応             | 対応           | 非対応            |
| BRICS フラグメント化               | 対応                          | 対応               | 非対応         | 非対応            |
| 3D 座標生成                        | 対応（ルールベース）          | 対応 (ETKDG)       | 対応           | 非対応            |
| PDB/XYZ ファイル形式               | 対応                          | 対応               | 対応           | 非対応            |
| CIP 立体化学 (R/S、E/Z)            | 対応                          | 対応               | 対応           | 非対応            |
| 力場エネルギー最小化               | 対応（ルールベース）          | 対応 (UFF/MMFF)    | 対応           | 非対応            |
| 反応 SMILES/SMIRKS                 | 対応                          | 対応               | 対応           | 非対応            |
| unsafe Rust                        | なし                          | 多数               | 多数           | なし              |
| メンテナンス状況 (2026)            | 活発                          | 活発               | 最小限         | アーカイブ済み    |

注:
- バイナリサイズは有効化する機能により異なる概算値。
- chemcore と purr はアーカイブ済み。chematic はそのスコープを包括する。

---

## ロードマップ

### Phase 1 — 基盤（完成）
コア型定義、OpenSMILES パース/ライター、ケクレ化、正規 SMILES。

### Phase 2 — 分子認識（完成）
SSSR (Balducci-Pearlman + GF(2))、Huckel 芳香族性認識、SDF/MOL V2000+V3000、2D SVG 描画。

### Phase 3 — 化学インテリジェンス（完成）
分子記述子（MW、LogP、TPSA、Fsp3）、QED、BRICS フラグメント化、
ECFP4/6 フィンガープリント、SMARTS + VF2（再帰 SMARTS・原子価・ハイブリッド化対応）、
分子標準化（塩除去・電荷中和）、Murcko スキャフォルド、CIP R/S および E/Z 立体化学。

### Phase 4 — 類似性と検索（完成）
MACCS 166 ビット構造キー ✓、位相的パス FP ✓、AtomPair FP ✓、Topological Torsion FP ✓、MCS ✓、互変異性体正規化 ✓。

### Phase 5 — 3D 化学（完成）
ルールベース 3D 座標生成、PDB/XYZ ファイル形式、UFF ライク最小化 ✓。

### Phase 6 — RDKit パリティ（完成）
反応 SMILES/SMIRKS ✓、フィーチャーフラグ付きアンブレラクレート ✓、
WASM npm パッケージ `@kent-tokyo/chematic` ✓、CPK 彩色 + ハイライト描画 ✓、
ChEMBL 37 全量検証（2,897,819 分子 / 100.000%）✓。

---

## リポジトリ構成

```
chematic/
├── Cargo.toml               ワークスペースルート
├── CHANGELOG.md             バージョン履歴
├── crates/
│   ├── chematic-core/       Atom, Bond, Molecule, Element、ケクレ化
│   ├── chematic-smiles/     OpenSMILES パーサー、ライター、正規 SMILES
│   ├── chematic-perception/ SSSR 環認識、Huckel 芳香族性認識
│   ├── chematic-mol/        MOL/SDF V2000+V3000 パーサーとライター
│   ├── chematic-depict/     2D SVG 描画エンジン（CPK カラー、ハイライト）
│   ├── chematic-chem/       分子記述子、BRICS、QED、標準化、スキャフォルド
│   ├── chematic-fp/         ECFP4/6、MACCS、パス、AtomPair、Torsion FP
│   ├── chematic-smarts/     SMARTS パーサー + VF2 部分構造一致（再帰 SMARTS）
│   ├── chematic-3d/         3D 座標生成、PDB/XYZ ファイル形式
│   ├── chematic-rxn/        反応 SMILES パーサーとライター
│   └── chematic/            フィーチャーフラグ付きアンブレラクレート
└── tasks/
    ├── todo.md              全フェーズロードマップチェックリスト（日本語）
    └── lessons.md           開発の教訓
```

---

## 開発コマンド

```bash
cargo build --workspace      # 全クレートのビルド
cargo test --workspace       # 全テストの実行（544 件）
cargo check --workspace      # ビルドなしの型チェック
cargo clippy --workspace     # リント
```

---

## ライセンス

Apache License 2.0 または MIT License のいずれかで利用可能。
