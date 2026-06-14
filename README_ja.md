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

全フェーズ完了 + v0.2.7–v0.2.9（MMFF94 完全スタック: 電荷 → エネルギーパラメータ → 幾何最適化器）+ **v0.2.10（L-BFGS 幾何最適化器 + WASM 5関数追加: MMFF94 minimize/L-BFGS/エネルギー内訳/トーションスキャン/電荷 + Demo 充実: MMFF94 UI・FF比較・トーションスキャン・電荷マップ）。** 1,951 テスト、全パス。C/C++ 依存ゼロ。

最新リリース: **v0.2.10**（2026-06-14）— L-BFGS 最小化器 + トーションスキャン + MMFF94 WASM API + Demo 充実

| クレート               | 説明                                                                                                                                      | テスト数 |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `chematic-core`        | Atom, Bond, Molecule, Element, ケクレ化（依存ゼロ）；ミュータブル API・`fragments`・`validate_valence`・`formula_with_isotopes`・`StereoGroup`/`StereoGroupKind` | 48      |
| `chematic-smiles`      | OpenSMILES パーサー、ライター、正規 SMILES、**CXSMILES メタデータ対応**                                                                  | 57      |
| `chematic-perception`  | SSSR、Hückel 芳香族性 + 反芳香族性（4n+2 則）、`apply_aromaticity`・`aromatize`・`kekulize_inplace`・`assign_stereo_from_2d`・`assign_ez_from_2d`・`cip_ez_descriptor` | 34      |
| `chematic-mol`         | MOL/SDF V2000+V3000（R/W、2D 座標付き）、CML（R/W）、CDXML（R）；`SdfRecord`（coords+props）、MDL RXN V2000 読み書き；V3000 ステレオグループ COLLECTION R/W | 63      |
| `chematic-depict`      | 2D SVG（CPK カラー・ハイライト・グリッド）、`detect_crossings`・`render_svg_with_metadata`・反応 SVG；Y座標系ドキュメント整備  | 43      |
| `chematic-chem`        | 70+ 記述子、タウトマー スコアリング、スキャフォルド ネットワーク、BRICS、QED、標準化、分子ハッシング、立体化学、`parse_condensed`、CIP、IFG、Gasteiger、VSA（EState+Labute）、`isotope_distribution`、`num_amide_bonds`、`num_ester_bonds` | 375     |
| `chematic-fp`          | ECFP2/4/6、FCFP4/6、MACCS 166-bit、TopoPF、AtomPair、Torsion、Layered、Pattern、Pharmacophore、Reaction FP — bitvec + Tanimoto/Dice；バルク類似度 | 50      |
| `chematic-ff`          | **MMFF94 完全スタック**: 数値原子型 (1–99)、電荷 (PBCI 99 + CHG 498, Halgren eq.15)、全エネルギーパラメータ (Bond 493 / Angle 2,245 / Torsion 926 / VdW 95)、steepest descent + **L-BFGS** 最適化器、トーションスキャン、エネルギー内訳；DREIDING 原子型付け | 98      |
| `chematic-smarts`      | SMARTS（再帰・原子価）、VF2（`MatchConfig`）、MCS（`match_chiral_tag` によるキラリティマッチング対応）                                    | 87      |
| `chematic-3d`          | 3D 座標生成、力場最小化、形状記述子、ConformerEnsemble（RMSD 剪定付き）、PDB/XYZ 形式                                                  | 147     |
| `chematic-rxn`         | 反応 SMILES/SMIRKS、`find_reaction_center` — `run_reactants`（生成物原子価バリデーション付き）                                        | 30      |
| `chematic-inchi`       | InChI/InChIKey 生成；formula/connectivity/hydrogen/stereo/charge/isotope レイヤー；環クロージャー対応                                | 28      |
| `chematic-wasm`        | **115+ WASM エクスポート** — npm: `@kent-tokyo/chematic` v0.2.10；MMFF94 minimize/L-BFGS/エネルギー内訳/トーションスキャン/電荷 | 175     |
| `chematic-iupac`       | ローカル IUPAC 命名（Pure Rust・オフライン）— アルカン/シクロアルカン、アルケン/アルキン、アルコール、アミン、ハロアルカン、アルデヒド、ケトン、カルボン酸、エステル、アミド、ベンゼン、芳香族ヘテロ環 | 14      |
| `chematic`             | フィーチャーフラグ付きアンブレラクレート（`iupac`, `inchi` フィーチャー追加）                                                         | 1       |

```
cargo test --workspace --lib   # 1,691 ライブラリテスト、全パス
```

---

## クイックスタート

### アンブレラクレートを使う場合

```toml
# Cargo.toml
[dependencies]
chematic = { version = "0.1.94", features = ["smiles", "fp", "chem", "mol", "depict"] }
```

### 個別クレートを使う場合

```toml
# Cargo.toml
[dependencies]
chematic-smiles     = "0.1.94"
chematic-perception = "0.1.94"
chematic-fp         = "0.1.94"
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

## CXSMILES でメタデータを保持

```rust
use chematic_smiles::parse_cxsmiles;

let cx = parse_cxsmiles("CCO |$ethanol$,atomProp:1.role.acceptor,^2:0|").unwrap();
// cx.atom_labels: ["ethanol"]
// cx.atom_props: [(atom: 1, key: "role", value: "acceptor")]
// cx.atom_radicals: [None, 2, None]

// CX 情報を保持したまま正規化
let canonical = chematic_smiles::write_cxsmiles(&cx);
println!("{}", canonical); // CCO |$ethanol$,atomProp:1.role.acceptor,^2:0|
```

---

## 標準化パイプラインと監査レポート

```rust
use chematic_chem::{StandardizationPipeline, StandardizeOptions};

let opts = StandardizeOptions {
    largest_fragment_only: true,
    neutralize_charges: true,
    remove_explicit_h: false,
    canonical_tautomer: false,
};

let pipeline = StandardizationPipeline::new(opts);
let (standardized, report) = pipeline.run(&mol);

// ステータスを確認（Unchanged / Modified / CompletedWithWarnings）
println!("Status: {:?}", report.status);

// 各ステップの変更を追跡
for step in &report.steps {
    if step.changed {
        println!("  {}: {} atoms → {} atoms",
            step.step.as_str(),
            step.before.atoms,
            step.after.atoms
        );
    }
}

// 警告（金属結合、原子価エラーなど）を確認
for warning in &report.warnings {
    println!("⚠️  {}: {}", warning.code, warning.message);
}
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

## 最近の開発（v0.1.89–v0.1.95）

**v0.1.95**: 真フィンガープリントアルゴリズム + IUPAC 命名拡張
- **MHFP 正規化ハッシュ**: Morgan 式循環フラグメントハッシュでアトムインデックス依存を解消。
- **ERG 薬理ノードタイプ**: DONOR/ACCEPTOR/POSITIVE/NEGATIVE/HYDROPHOBIC を正確に付与。ピリジン N と ピロール N を区別。
- **Reaction FP**: XOR 差分エンコーディングがデフォルトであることを確認。
- **IUPAC 命名**: 位置番号付きケトン、カルボン酸、エステル、アミド、ベンゼン、芳香族ヘテロ環に対応。
- 1,657 ライブラリテスト、全パス。

**v0.1.91–v0.1.94: RDKit ギャップ分析クロージャー（A1–A5, B3）**
- **v0.1.91**: True MHFP（構造フラグメント ハッシング）、True ERG（Ertl 2017 機能基）
- **v0.1.92**: パス FP に結合型を含む、InChI 立体層パース（`/t`, `/b`）
- **v0.1.93**: 完全な多球 CIP 優先度ルール（chematic-perception に移動、循環依存回避）
- **v0.1.94**: SA Score コーパス拡張（145 → 188 FDA 分子、1034 → 1415 ユニークフラグメント）

**v0.1.88–v0.1.90**: InChI 立体層、Brenk SMARTS、再イオン化、グループ標準化

**v0.1.69–v0.1.87**: RDKit ギャップ分析初期 — SSSR、ケクレ化、CIP、3D 幾何、WASM API 成熟度

詳細な歴史的ロードマップ（Phase 1–16, v0.1.14–v0.1.33）は `tasks/todo.md` を参照。

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
