# chematic

[English](README.md)

Pure Rust 製のケモインフォマティクスライブラリ。RDKit の代替を目指す、C/C++ FFI ゼロの Rust 実装。

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

Phase 1〜3 および Phase 4（MACCS・パスフィンガープリント）、Phase 5（座標生成 + ファイル形式）が完了。
263 テスト、全パス。

| クレート               | 説明                                                                    | テスト数 |
|------------------------|-------------------------------------------------------------------------|---------|
| `chematic-core`        | Atom, Bond, Molecule, Element, ケクレ化（依存ゼロ）                     | 30      |
| `chematic-smiles`      | OpenSMILES パーサー、ライター、正規 SMILES                              | 50      |
| `chematic-perception`  | SSSR (Balducci-Pearlman)、Huckel 芳香族性認識                           | 14      |
| `chematic-mol`         | MOL/SDF V2000+V3000 パーサーとライター                                  | 34      |
| `chematic-depict`      | 2D SVG 描画（環・鎖テンプレート）                                       | 14      |
| `chematic-chem`        | 記述子、標準化（塩除去・電荷中和）、Murcko スキャフォルド               | 38      |
| `chematic-fp`          | ECFP4/6、MACCS 166-bit 構造キー、位相的パス FP、Tanimoto/Dice 類似度   | 31      |
| `chematic-smarts`      | SMARTS パーサー、VF2 部分構造一致                                       | 34      |
| `chematic-3d`          | 3D 座標生成、PDB/XYZ ファイル形式                                       | 15      |

```
cargo test --workspace   # 263 テスト、全パス
```

---

## クイックスタート

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
    let arom = assign_aromaticity(&benzene);
    println!("芳香族原子数: {}", arom.aromatic_atom_count()); // 6

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
let query = parse_smarts("C=O").unwrap();
let matches = find_matches(&query, &mol);
println!("C=O 基の数: {}", matches.len()); // 2
```

---

## 分子記述子計算

```rust
use chematic_smiles::parse;
use chematic_chem::{molecular_weight, tpsa, lipinski_passes};

let aspirin = parse("CC(=O)Oc1ccccc1C(=O)O").unwrap();
println!("分子量:    {:.2}", molecular_weight(&aspirin)); // ~180.16
println!("TPSA:      {:.2}", tpsa(&aspirin));             // ~63.6
println!("Lipinski:  {}", lipinski_passes(&aspirin));     // true
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

---

## 他のケモインフォマティクスライブラリとの比較

| 機能                           | chematic               | RDKit (rdkit-sys)  | OpenBabel FFI  | chemcore / purr   |
|--------------------------------|------------------------|--------------------|----------------|-------------------|
| 実装言語                       | Pure Rust              | Rust + C++ FFI     | Rust + C++ FFI | Pure Rust         |
| WASM ターゲット                | 対応                   | 非対応             | 非対応         | 部分対応          |
| バイナリサイズ（コア）         | 約 500 KB              | 約 50 MB           | 約 20 MB       | 約 200 KB         |
| OpenSMILES パーサー            | 完全実装               | 完全実装           | 完全実装       | 部分実装          |
| SMILES ライター                | 対応                   | 対応               | 対応           | 非対応            |
| 正規 SMILES                    | 対応                   | 対応               | 対応           | 非対応            |
| ケクレ化                       | 対応                   | 対応               | 対応           | 非対応            |
| 芳香族性認識                   | 対応 (Huckel 則)       | 対応               | 対応           | 部分対応          |
| 環認識 (SSSR)                  | 対応                   | 対応               | 対応           | 非対応            |
| SDF/MOL V2000                  | 対応                   | 対応               | 対応           | 非対応            |
| SDF/MOL V3000                  | 対応                   | 対応               | 対応           | 非対応            |
| 2D 描画 (SVG)                  | 対応                   | 対応               | 対応           | 非対応            |
| ECFP フィンガープリント        | 対応 (ECFP4/6)         | 対応               | 対応           | 非対応            |
| SMARTS / 部分構造検索          | 対応 (VF2)             | 対応               | 対応           | 非対応            |
| 分子記述子計算                 | 対応 (MW/LogP/TPSA/…)  | 対応               | 対応           | 非対応            |
| 3D 座標生成                    | 対応（ルールベース）   | 対応 (ETKDG)       | 対応           | 非対応            |
| PDB/XYZ ファイル形式           | 対応                   | 対応               | 対応           | 非対応            |
| CIP 立体化学 (R/S)             | 予定                   | 対応               | 対応           | 非対応            |
| MACCS フィンガープリント       | 対応 (166-bit 構造キー) | 対応               | 対応           | 非対応            |
| 力場エネルギー最小化           | 予定                   | 対応 (UFF/MMFF)    | 対応           | 非対応            |
| 反応 SMILES/SMIRKS             | 予定                   | 対応               | 対応           | 非対応            |
| unsafe Rust                    | なし                   | 多数               | 多数           | なし              |
| メンテナンス状況 (2026)        | 活発                   | 活発               | 最小限         | アーカイブ済み    |

注:
- "chematic" 列は現在の実装に加え、全フェーズ完了後の最終予定状態を示す。
- バイナリサイズは有効化する機能により異なる概算値。
- chemcore と purr はアーカイブ済み。chematic はそのスコープを包括する。

---

## ロードマップ

### Phase 1 — 基盤（完成）
コア型定義、OpenSMILES パース/ライター、ケクレ化、正規 SMILES。80 テスト。

### Phase 2 — 分子認識（完成）
SSSR (Balducci-Pearlman + GF(2))、Huckel 芳香族性認識、SDF/MOL V2000+V3000、2D SVG 描画。63 テスト追加。

### Phase 3 — 化学インテリジェンス（完成）
分子記述子（MW、LogP、TPSA、Lipinski）、ECFP4/6 フィンガープリント、SMARTS + VF2 部分構造検索、
分子標準化（塩除去・電荷中和）、Murcko スキャフォルド。残り: CIP R/S 立体化学割り当て。

### Phase 4 — 類似性と検索（一部完成）
MACCS 166 ビット構造キー ✓、位相的パスフィンガープリント ✓。
残り: 最大共通部分構造 (MCS)、互変異性体正規化。

### Phase 5 — 3D 化学（一部完成）
ルールベース 3D 座標生成、PDB/XYZ ファイル形式。
残り: UFF 力場エネルギー最小化。

### Phase 6 — RDKit パリティ
WASM パッケージ (npm: chematic)、反応 SMILES/SMIRKS、フィーチャーフラグ付きアンブレラクレート、ChEMBL スケール検証。

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
│   ├── chematic-depict/     2D SVG 描画エンジン
│   ├── chematic-chem/       分子記述子、標準化、スキャフォルド
│   ├── chematic-fp/         ECFP4/6、MACCS、位相的パス FP、類似度計算
│   ├── chematic-smarts/     SMARTS パーサー + VF2 部分構造一致
│   └── chematic-3d/         3D 座標生成、PDB/XYZ ファイル形式
└── tasks/
    ├── todo.md              全フェーズロードマップチェックリスト（日本語）
    └── lessons.md           開発の教訓
```

---

## 開発コマンド

```bash
cargo build --workspace      # 全クレートのビルド
cargo test --workspace       # 全テストの実行（260+ 件）
cargo check --workspace      # ビルドなしの型チェック
cargo clippy --workspace     # リント
```

---

## ライセンス

Apache License 2.0 または MIT License のいずれかで利用可能。
