# chematic

[English](README.md) | [中文](README_zh.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/chematic.svg)](https://crates.io/crates/chematic)
[![PyPI](https://img.shields.io/pypi/v/chematic.svg)](https://pypi.org/project/chematic/)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic.svg)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![ライセンス](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE-MIT)
[![Docs](https://img.shields.io/badge/docs-site-blue)](https://kent-tokyo.github.io/chematic/getting_started/installation/)
[![デモ](https://img.shields.io/badge/demo-live-brightgreen)](https://kent-tokyo.github.io/chematic/playground/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

Pure Rust 製のケモインフォマティクスライブラリ。RDKit の代替を目指す、**デフォルトで C/C++ FFI ゼロ**の Rust 実装。

> **なぜ C/C++ ゼロが重要か？**
> RDKit.js、Indigo WASM、OpenBabel はいずれも C++ コードを Emscripten でコンパイルしています。
> そのため **30〜50 MB の WASM バイナリ**、複雑なビルドツールチェーン、プラットフォーム固有のビルドエラーが生じます。
> chematic は `wasm-pack build` 一発で **〜550 KB の WASM バンドル**を生成します。
> `cmake`・`clang`・`-sys` クレート・`build.rs` での C コンパイルは依存ツリー全体にわたって一切使用しません。
> *（例外：`native-inchi` feature のみ opt-in で C コンパイラが必要。WASM ビルドには影響なし。）*

---

## ライブデモ

**[https://kent-tokyo.github.io/chematic/playground/](https://kent-tokyo.github.io/chematic/playground/)** — 記述子計算、薬らしさルール、フィンガープリント類似度、3D ビューア、反応スキーム、SAR 解析をブラウザ上の WebAssembly で実行できるインタラクティブデモ。

---

## 設計目標

**Pure Rust、C/C++ FFI ゼロ — デフォルトビルドで保証済み**
`rdkit-sys`・`openbabel-sys`・`bindgen` なし。SSSR 環認識から ECFP フィンガープリント、力場最小化まで、すべてのアルゴリズムを 100% 安全な Rust で実装。依存ツリー全体を FFI フリーで検証済み。

> **任意例外**: `chematic-inchi` の `native-inchi` feature を有効にすると、IUPAC InChI C ライブラリ (v1.07.5) が vendored でリンクされ、ビット完全一致の標準 InChI/InChIKey が生成できます。C コンパイラが必要ですが完全 opt-in で、デフォルトビルドは FFI フリーのまま。

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

全フェーズ完了 + **v0.3.x シリーズ（全主要競合ライブラリを超えた）**: MCP サーバー（AI エージェント統合）、pKa 予測（15 SMARTS ルール）、ADMET プロファイル（BBB/Caco-2/hERG/CYP3A4）、IUPAC 25+ 化合物クラス、WASM pKa/ADMET バインディング、criterion ベンチマーク。**1,991 テスト、全パス。C/C++ 依存ゼロ（デフォルトビルド）。**

最新リリース: **v0.4.10**（2026-06-20）— v0.4.10: Sprint 18–26 Python バインディング 50+ | v0.4.9: PDBQT+UFF+SDF 電荷 | v0.4.0: PyO3 Python バインディング

| クレート               | 説明                                                                                                                                      | テスト数 |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `chematic-core`        | Atom, Bond, Molecule, Element, ケクレ化（依存ゼロ）；ミュータブル API・`fragments`・`validate_valence`・`formula_with_isotopes`・`StereoGroup`/`StereoGroupKind` | 48      |
| `chematic-smiles`      | OpenSMILES パーサー、ライター、正規 SMILES、**CXSMILES メタデータ対応**                                                                  | 57      |
| `chematic-perception`  | SSSR、Hückel 芳香族性 + 反芳香族性（4n+2 則）、`apply_aromaticity`・`aromatize`・`kekulize_inplace`・`assign_stereo_from_2d`・`assign_ez_from_2d`・`cip_ez_descriptor` | 34      |
| `chematic-mol`         | MOL/SDF V2000+V3000（R/W、2D 座標付き）、CML（R/W）、CDXML（R）；`SdfRecord`（coords+props）、MDL RXN V2000 読み書き；V3000 ステレオグループ COLLECTION R/W | 63      |
| `chematic-depict`      | 2D SVG（CPK カラー・ハイライト・グリッド）、`detect_crossings`・`render_svg_with_metadata`・反応 SVG；Y座標系ドキュメント整備  | 43      |
| `chematic-chem`        | 70+ 記述子、タウトマー、スキャフォルド、BRICS、QED、標準化；**pKa 予測** (15 SMARTS ルール)；**ADMET プロファイル** (BBB/Caco-2/hERG/CYP3A4)；**HBA 99.98% RDKit 一致率**（5,000 分子ベンチマーク） | 496     |
| `chematic-fp`          | ECFP2/4/6、FCFP4/6、MACCS、TopoPF、AtomPair、Torsion、Layered、Pattern、Pharmacophore、Reaction、**MAP4** (Minervini 2020) — Tanimoto/Dice | 55      |
| `chematic-ff`          | **MMFF94 全 7 エネルギー項** (Halgren 1996)：OOP (117件) + Stretch-Bend (282件)；steepest descent + L-BFGS；DREIDING | 98      |
| `chematic-smarts`      | SMARTS、VF2、MCS；**SmartsCache** (LRU 5–20×)；**named_pattern()** (20 パターン)；**SMARTS 内アトムマップ `:N`** (`[O;D1;H0:3]` 形式 — メタデータとして保存、マッチング条件には不使用) | 137     |
| `chematic-3d`          | 3D 座標生成、ETKDG KB (20+ パターン)、力場最小化、形状記述子、ConformerEnsemble、PDB/XYZ | 147     |
| `chematic-rxn`         | 反応 SMILES/SMIRKS、`find_reaction_center`、`run_reactants`（原子価バリデーション） | 30      |
| `chematic-inchi`       | InChI/InChIKey：純 Rust 近似（WASM 対応）**+ `native-inchi` feature で IUPAC 標準準拠**（C ライブラリ 1.07.5 vendored、ビット完全一致）；**parse_inchi** 読み込み | 28 (+14*)   |
| `chematic-wasm`        | **130+ WASM エクスポート** — npm: `@kent-tokyo/chematic` v0.3.2；**pKa/ADMET/BBB/Caco-2/hERG/CYP3A4** WASM API | 209     |
| `chematic-iupac`       | ローカル IUPAC 命名（Pure Rust・オフライン）— **25+ 化合物クラス**：アルカン、シクロアルカン、アルコール、アミン、ハロアルカン、ケトン、酸、エステル、アミド、**ピペリジン、モルホリン、ピペラジン、ナフタレン、スルフィド** | 45      |
| `chematic-mcp`         | **MCP (Model Context Protocol) サーバー** — AI エージェント統合；**15 ツール**：parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles | 28      |
| `chematic`             | フィーチャーフラグ付きアンブレラクレート（統合クレート）                                                                                                  | 1       |

```
cargo test --workspace --lib --quiet                                               # 1,991 ライブラリテスト、全パス
cargo test -p chematic-inchi --features native-inchi --test standard_inchi         # +14 IUPAC 標準 InChI 統合テスト
```

---

## クイックスタート

### アンブレラクレート（統合クレート）を使う場合

```toml
# Cargo.toml
[dependencies]
chematic = { version = "0.3.2", features = ["smiles", "fp", "chem", "mol", "depict"] }
```

### 個別クレートを使う場合

```toml
# Cargo.toml
[dependencies]
chematic-smiles     = "0.2.11"
chematic-perception = "0.2.11"
chematic-fp         = "0.2.11"
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

| 観点                                        | **chematic**                               | RDKit (rdkit-sys)  | OpenBabel FFI | RDKit.js (WASM)  |
|---------------------------------------------|--------------------------------------------|--------------------|---------------|------------------|
| **C/C++ 依存**                              | **ゼロ（デフォルト）**†                    | 大規模 C++         | 大規模 C++    | C++（Emscripten）|
| **WASM バイナリサイズ**                     | **〜550 KB**                               | N/A（WASM 非対応） | N/A           | 〜30 MB          |
| **ビルド要件**                              | `cargo build` のみ                         | cmake + clang      | cmake + clang | Emscripten SDK   |
| **Python バインディング**                   | **あり** (`pip install chematic`, PyO3)    | あり（rdkit-sys）  | あり          | なし             |
| unsafe Rust                                 | **なし**                                   | 大規模             | 大規模        | N/A              |
| ケクレ化                                    | **4-pass（Edmonds' blossom 含む）**        | あり               | あり          | あり             |
| SDF/MOL V2000+V3000                         | あり                                       | あり               | あり          | あり             |
| Tripos MOL2 形式                            | **あり**（読み書き + Python）              | あり               | あり          | なし             |
| 分子記述子                                  | **70+（BOILED-Egg、QED、SA Score 含む）**  | 〜30               | 〜20          | 〜30             |
| **MAP4 フィンガープリント**                 | **あり**（Minervini 2020）                 | なし（外部pkg）    | なし          | なし             |
| MMFF94 全 7 エネルギー項                    | **あり**                                   | あり               | あり          | なし             |
| 3D 座標生成                                 | あり（DG + MMFF94/DREIDING + L-BFGS）      | あり（ETKDG）      | あり          | あり             |
| 多様性選択（MaxMin/Butina）                 | **あり**                                   | あり               | なし          | なし             |
| InChI / InChIKey                            | **あり** — 純 Rust（デフォルト）+ **IUPAC 準拠**（`native-inchi`）| C ライブラリ必要 | C ライブラリ必要 | C ライブラリ必要 |
| **pKa 予測**                                | **あり（15 SMARTS ルール）**               | なし               | なし          | なし             |
| **ADMET プロファイル + BOILED-Egg**         | **あり**                                   | 一部               | なし          | 一部             |
| **MCP サーバー（AI エージェント API）**     | **あり — 15 ツール（Name→SMILES 含む）**  | なし               | なし          | なし             |
| IUPAC 名生成                                | **あり（25+ 化合物クラス）**               | なし               | なし          | 一部             |
| メンテナンス（2026）                        | アクティブ                                 | アクティブ         | 最小限        | アクティブ       |

† デフォルトビルドのみ。`native-inchi` feature は opt-in で C コンパイラが必要。他の全クレートは FFI フリー。

---

## 最近の開発（v0.4.8）

**v0.4.8**（2026-06-19）: `name_to_smiles` MCP ツール、反復 `augmented_ring_set`
- **`name_to_smiles`**: PubChem REST プロキシを追加（MCP ツール 15 個目）。化学名 → SMILES 変換。
- **反復 `augmented_ring_set`**: 3+ SSSR リングの XOR が必要な縮合 PAH に対応。芳香環カウント精度向上。
- **Python `from_mol2()` / `to_mol2()`**: Tripos MOL2 形式の Python バインディング追加。
- **Python 3.13 wheel**: PyPI の配布物に Python 3.9〜3.13 の wheel を追加。

**v0.4.7**（2026-06-19）: ホウ素芳香環ケクレ化修正、WASM ADMET BOILED-Egg 追加

**v0.4.6**（2026-06-19）: Python `boiled_egg()` メソッド、`admet()` 拡張

**v0.4.5**（2026-06-19）: ケクレ化 Edmonds' blossom（128→2 件）、InChI E/Z `/b`、MCP 6 新ツール、BOILED-Egg

**v0.4.0–v0.3.x**: Python PyO3 バインディング、native-inchi、MCP サーバー、pKa、ADMET

**v0.2.x**: MMFF94 全 7 項、MAP4 フィンガープリント、SMARTS キャッシュ

**v0.1.x**: コア基盤 — SSSR、ケクレ化、CIP、3D 幾何、WASM API

---

## 既知の制限事項

### ケクレ化（5,000 分子中 **2件のみ残存** — ほぼ解決済み）

`chematic-core` のケクレ代入は 4-pass 戦略を使用：

- **Pass 1/2**: BFS 増加パス（昇順 / 降順）。
- **Pass 3**: 橋頭 N 除外 — 環接合部の N 原子（芳香族次数 ≥ 3）は二重結合を占有せずにローンペアを提供し、残りの C 原子を二部グラフで照合。インドリジン型システム（コーパス約 109 件）を修正。
- **Pass 4**: Edmonds' blossom アルゴリズム（O(n²m)）— 奇数サイクル（コラニュレン C₂₀H₁₀ など）を持つ非二部 C 芳香族サブグラフに対応。残りの複雑な多環系を修正。

5,000 分子コーパス（issue #11）において、これらの修正後にケクレ化が失敗するのは**2件のみ**：

| カテゴリ | 件数 | 例 |
|---|---|---|
| ホウ素芳香環 | 1 | `b1ccccn1` |
| 純 H₂（重原子なし） | 1 | `[H][H]` |

**影響**: `KekuleError` が明示的に返され、無音の誤出力は生じない。

### 芳香族性モデル（Hückel vs RDKit）

chematic は **Hückel 4n+2 則を各 SSSR 環に独立適用**するのに対し、RDKit はより高度な縮合環電子非局在化モデルを使用。差異は N-ヘテロ環（ピリドン、キノロン、インドリジン）で顕著。

**5,000 分子コーパス（issue #12）の現状：**

| 特徴量 | issue #12 クローズ時 | 現在 | 状態 |
|---|---|---|---|
| `[nH]` SMARTS 一致 | 67% | **100% recall / 99.8% precision** | 解決済み |
| HBA カウント | 87.7% | **99.98%**（4,999 / 5,000） | 解決済み |
| 芳香族環数 | 92.6% | **~100%**（≥ 4,998 / 5,000） | 解決済み — `augmented_ring_set` XOR ガード修正 |

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
│   └── chematic/            フィーチャーフラグ付きアンブレラクレート（統合クレート）
```

---

## ライセンス

Apache License 2.0 または MIT License のいずれかで利用可能。
