# chematic

[English](README.md) | [中文](README_zh.md)

Pure Rust 製のケモインフォマティクスライブラリ。RDKit の代替を目指す、**C/C++ FFI ゼロ**の Rust 実装。

> **なぜ C/C++ ゼロが重要か？**
> RDKit.js、Indigo WASM、OpenBabel はいずれも C++ コードを Emscripten でコンパイルしています。
> そのため **30〜50 MB の WASM バイナリ**、複雑なビルドツールチェーン、プラットフォーム固有のビルドエラーが生じます。
> chematic は `wasm-pack build` 一発で **〜550 KB の WASM バンドル**を生成します。
> `cmake`・`clang`・`-sys` クレート・`build.rs` での C コンパイルは依存ツリー全体にわたって一切使用しません。

---

## ライブデモ

**[https://kent-tokyo.github.io/chematic/](https://kent-tokyo.github.io/chematic/)** — 記述子計算、薬らしさルール、フィンガープリント類似度、3D ビューア、反応スキーム、SAR 解析をブラウザ上の WebAssembly で実行できるインタラクティブデモ。

---

## 設計目標

**Pure Rust、C/C++ FFI ゼロ — 保証済み**
`rdkit-sys`・`openbabel-sys`・`cc` ビルド依存・`bindgen` なし。SSSR 環認識から ECFP フィンガープリント、力場最小化まで、すべてのアルゴリズムを 100% 安全な Rust で実装。依存ツリー全体を FFI フリーで検証済み。

**WASM 対応、軽量**
全クレートが `wasm32-unknown-unknown` に無修正でコンパイルできる。npm パッケージ `@kent-tokyo/chematic` は **〜550 KB**（C++ FFI 代替の 30〜50 MB とは対照的）。`cmake`・`emcc`・Emscripten ツールチェーン不要。

**100+ WebAssembly API**
WASM レイヤーは記述子・フィンガープリント・スキャフォルド解析・立体異性体列挙・3D ジオメトリ・多様性選択・MMP 分析・R-group 分解・分子編集などを網羅した 100 以上の関数を公開。TypeScript 型定義付き。

**化学ドメイン固有の実装**
汎用グラフライブラリのラッパーではなく、ケクレ化・Hückel 芳香族性・CIP 立体化学・SSSR 環認識・Gasteiger 電荷・MaxMin/Butina 多様性ピッキングを Rust でスクラッチ実装。

**再現性と決定性**
フィンガープリントは固定不変量順序の FNV-1a ハッシュを使用。同じ SMILES 入力から常に同じビット列が得られる。乱数なし、プラットフォーム依存なし。

---

## 現在のステータス

全フェーズ完了。**933 テスト、全パス。C/C++ 依存ゼロ。**

| クレート               | 説明                                                                                                                                      | テスト数 |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `chematic-core`        | Atom, Bond, Molecule, Element, ケクレ化（依存ゼロ）；ミュータブル API・`fragments`・`validate_valence`・`formula_with_isotopes` | 48      |
| `chematic-smiles`      | OpenSMILES パーサー、ライター、正規 SMILES                                                                                               | 57      |
| `chematic-perception`  | SSSR、Hückel 芳香族性、`apply_aromaticity`・`aromatize`・`kekulize_inplace`・`assign_stereo_from_2d`                                     | 18      |
| `chematic-mol`         | MOL/SDF V2000+V3000（R/W）、CML（R/W）、CDXML（R）；`SdfRecord`（coords+props）、MDL RXN V2000 読み書き                                  | 61      |
| `chematic-depict`      | 2D SVG 描画（CPK カラー・ハイライト・グリッド）、DepictData、`suggest_bond_direction`・反応 SVG                                          | 39      |
| `chematic-chem`        | 40+ 記述子（`xlogp3` 含む）、BRICS（`BricsConfig`）、QED、標準化、CIP、IFG、Gasteiger、`expand_abbreviation`                             | 226     |
| `chematic-fp`          | ECFP2/4/6、FCFP4/6、MACCS 166-bit、TopoPF、AtomPair、Torsion FP — bitvec + Tanimoto/Dice                                               | 50      |
| `chematic-smarts`      | SMARTS（再帰・原子価）、VF2（`MatchConfig`）、MCS（`AtomCompare`/`BondCompare`/ring-awareness）                                          | 84      |
| `chematic-3d`          | 3D 座標生成、力場最小化、形状記述子、ConformerEnsemble、PDB/XYZ 形式                                                                    | 68      |
| `chematic-rxn`         | 反応 SMILES/SMIRKS — `run_reactants`（生成物原子価バリデーション付き）                                                                   | 28      |
| `chematic-wasm`        | **100+ WASM エクスポート** — npm: `@kent-tokyo/chematic`                                                                                 | 162     |
| `chematic-iupac`       | ローカル IUPAC 命名（Pure Rust・オフライン）— アルカン、シクロアルカン、アルコール、アミン、ハロアルカン                                | 8       |
| `chematic`             | フィーチャーフラグ付きアンブレラクレート（`iupac` フィーチャー追加）                                                                    | 1       |

```
cargo test --workspace   # 933 テスト、全パス
```

---

## クイックスタート

### アンブレラクレートを使う場合

```toml
# Cargo.toml
[dependencies]
chematic = { version = "0.1.23", features = ["smiles", "fp", "chem", "mol", "depict"] }
```

### 個別クレートを使う場合

```toml
# Cargo.toml
[dependencies]
chematic-smiles     = "0.1.23"
chematic-perception = "0.1.23"
chematic-fp         = "0.1.23"
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

## JavaScript / TypeScript（WebAssembly）

> **〜550 KB、C/C++ 依存ゼロ。** ブラウザ・Node.js に対応。
> RDKit.js は Emscripten ビルドで〜30 MB。

```sh
npm install @kent-tokyo/chematic
```

```js
import init, {
  parse_smiles, canonical_tautomer, murcko_scaffold,
  tanimoto_ecfp4, tanimoto_ecfp6, tanimoto_maccs,
  brics_fragments_json, mcs_smiles_json,
  get_descriptors_json, enumerate_stereo_isomers_json,
  sdf_to_records_json, sdf_from_records_json,
  mmp_pairs_json, rgroup_decompose_json,
  mol_with_atom_added, mol_with_atom_charge, mol_with_atom_element,
  depict_data_json, cpk_color,
} from '@kent-tokyo/chematic';

await init();

// ── パース・記述子 ───────────────────────────────────────────
const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // アスピリン
console.log(mol.molecular_weight()); // ~180.16
console.log(mol.qed());              // ドラッグライクネス [0,1]

// 全記述子を一括取得（JSON オブジェクト）
const desc = JSON.parse(get_descriptors_json(mol));
console.log(desc.mw, desc.tpsa, desc.logP);

// ── 立体異性体列挙 ───────────────────────────────────────────
const isomers = JSON.parse(enumerate_stereo_isomers_json(parse_smiles('C(F)(Cl)Br')));
// ["[C@@H](F)(Cl)Br","[C@H](F)(Cl)Br"]

// ── SAR 解析 ─────────────────────────────────────────────────
const smiles_json = '["CCc1ccccc1","CCCc1ccccc1","CCCCc1ccccc1"]';
const pairs = JSON.parse(mmp_pairs_json(smiles_json));
// [{mol_a:"...", mol_b:"...", core:"...", fragment_a:"...", fragment_b:"..."}]

const rgroups = JSON.parse(rgroup_decompose_json(smiles_json, 'c1ccc(*)cc1'));
// [{matched:true, r1:"C(C)[*]"}, ...]

// ── 分子編集 API ─────────────────────────────────────────────
const mol2 = mol_with_atom_added(mol, 'N');
const mol3 = mol_with_atom_charge(mol, 0, 1); // 原子 0 を +1 に
const mol4 = mol_with_atom_element(mol, 0, 'O'); // 原子 0 を O に変更
```

---

## 他のケモインフォマティクスライブラリとの比較

| 観点                               | **chematic**             | RDKit.js (WASM)  | OCL.js | Indigo WASM |
|------------------------------------|--------------------------|------------------|--------|-------------|
| **C/C++ 依存**                     | **ゼロ — Pure Rust**     | C++（Emscripten）| △      | C++（Emscripten）|
| **WASM バイナリサイズ**            | **〜550 KB**             | 〜30 MB          | 〜5 MB | 〜10 MB     |
| **ビルド要件**                     | `cargo build` のみ       | cmake + clang    | —      | Emscripten SDK |
| 記述子の豊富さ                     | **◎ 40+**                | ○ 〜30           | △      | △           |
| FP 種類・設定自由度                | **◎ 7 種 bitvec + 類似度**| ◎               | ○      | △           |
| 立体化学（CIP + 列挙）             | **◎**                    | ○                | ○      | △           |
| 3D + コンフォーマー管理            | **◎**                    | ○                | △      | △           |
| 多様性選択（MaxMin/Butina）        | **◎**                    | ○                | ✗      | ✗           |
| MMP 分析                           | ✓                        | ✓                | ✗      | ✗           |
| R-group 分解                       | ✓                        | ✓                | ✗      | ✗           |
| 分子編集 API                       | **◎ with_atom_* 系**     | ○                | ○      | ○           |
| CML 読み書き                       | ✓                        | ✓                | ✓      | ✓           |
| CDXML 読み込み（複数フラグメント・立体化学）| ✓               | ✓                | ✓      | ✓           |
| InChI / InChIKey                   | ✗（C ライブラリ依存）    | ✓                | ✓      | ✓           |
| unsafe Rust                        | **なし**                 | —                | —      | —           |

---

## ロードマップ

### Phase 1〜6（完成）
基盤・分子認識・化学インテリジェンス・類似性・3D・エコシステム。

### Phase 7（完成）
拡張記述子・多様性・SA スコア・EState・IFG・Gasteiger・VSA。

### Phase 8（v0.1.20〜v0.1.22、完成）
100+ WASM エクスポート・CML/CDXML・Mutable Molecule API・DepictData・MMP・R-group・ConformerEnsemble・SDF/V3000 write・MCS ring-awareness 制約。

### Phase 15（v0.1.29〜32、完成）
ミュータブル `Molecule`（`add/remove_atom/bond`・`fragments`・`is_connected`）、
`assign_stereo_from_2d`（ウェッジ結合→R/S）、`aromatize`/`kekulize_inplace`、
`depict_reaction_svg`、`SdfRecord`（coords+properties 統合）、MDL RXN V2000 読み書き、
`expand_abbreviation`（30 略号）、`formula_with_isotopes`。

### Phase 14（v0.1.28、完成）
`xlogp3()` (Cheng 2007 原子型)、`chematic-iupac`（純 Rust オフライン IUPAC 命名）、
`BricsConfig { min_fragment_size }`、`MatchConfig { max_matches }`、
`McsConfig { atom_compare: AtomCompare, bond_compare: BondCompare }` でヘテロ環 scaffold hopping 対応。

### Phase 13（v0.1.27、完成）
`MolMetadata::default().with_name("アスピリン").with_comment("...")` — MOL/SDF メタデータ用 fluent builder。

### Phase 12（v0.1.26、完成）
`atom_color_rgb(atomic_number: u8) -> [u8; 3]` — hex 解析なしで CPK カラーを RGB バイトトリプルとして取得。

### Phase 11（v0.1.25、完成）
`suggest_bond_direction(mol, atom, layout) -> f64`（ラジアン）: sp2/sp3 角度オフセット + 最大最小分離角選択による化学的に自然な新規結合方向提案。`BOND_LEN` 定数を公開。

### Phase 10（v0.1.24、完成）
`validate_valence(mol) -> Vec<ValenceError>` 公開 API（chematic-core + chematic-perception 経由で参照可能）、`run_reactants` が過原子価の生成物セットを自動除外。

### Phase 9（v0.1.23、完成）
`Element::vdw_radius()` / `covalent_radius()`（Bondi/Alvarez テーブル 118 元素）、
`Molecule::implicit_hydrogen_count()` / `total_formula()`（暗黙的 H を含む Hill 式）、
`apply_aromaticity()`（ケクレ化分子 → 芳香族フラグ適用 Molecule）、
`with_atom_aromatic()` / `with_bond_order()` immutable update API 拡張、
`minimize_uff()` エイリアス（UFF 力場最小化の発見性向上）。

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
│   ├── chematic-mol/        MOL/SDF V2000+V3000、CML、CDXML
│   ├── chematic-depict/     2D SVG 描画エンジン（CPK カラー、DepictData）
│   ├── chematic-chem/       記述子、BRICS、QED、MMP、標準化、CIP
│   ├── chematic-fp/         ECFP4/6、FCFP4/6、MACCS、AtomPair、Torsion FP
│   ├── chematic-smarts/     SMARTS パーサー + VF2 部分構造一致（再帰 SMARTS）
│   ├── chematic-3d/         3D 座標生成、ConformerEnsemble、PDB/XYZ 形式
│   ├── chematic-rxn/        反応 SMILES/SMIRKS
│   └── chematic/            フィーチャーフラグ付きアンブレラクレート
└── tasks/
    ├── todo.md              全フェーズロードマップチェックリスト（日本語）
    └── lessons.md           開発の教訓
```

---

## ライセンス

Apache License 2.0 または MIT License のいずれかで利用可能。
