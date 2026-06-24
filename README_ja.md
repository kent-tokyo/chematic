# chematic

[English](README.md) | [中文](README_zh.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![crates.io](https://img.shields.io/crates/v/chematic.svg)](https://crates.io/crates/chematic)
[![PyPI](https://img.shields.io/pypi/v/chematic.svg)](https://pypi.org/project/chematic/)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic.svg)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![ライセンス](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue.svg)](LICENSE-MIT)
[![Docs](https://img.shields.io/badge/docs-site-blue)](https://kent-tokyo.github.io/chematic/)
[![デモ](https://img.shields.io/badge/demo-live-brightgreen)](https://kent-tokyo.github.io/chematic/playground/)
[![Open in Colab](https://colab.research.google.com/assets/colab-badge.svg)](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

Python・Rust・ブラウザ向けケモインフォマティクスライブラリ。

**デフォルトで速く、設計で安全なケモインフォマティクス。**  
Pure Rust · C/C++ ゼロ · Python · WebAssembly · [ライブデモ](https://kent-tokyo.github.io/chematic/playground/)

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **導入方法** | `pip install chematic` | conda / cmake が必要 | Python バインディングなし |
| **ブラウザ向けバンドル** | **504 KB** | 提供なし | ~30 MB（60× 大きい） |
| **バッチ FP 速度** | **3.6 µs/mol**（5–14× 高速） | 20–50 µs/mol | — |
| **メモリ安全性** | コンパイラが保証（Rust） | C++ | C++ |
| **ソースビルド** | `cargo build` のみ | cmake + clang + Boost | Emscripten SDK |

すべての数値は再現可能です — [ベンチマーク詳細](https://kent-tokyo.github.io/chematic/benchmark/)を参照。  
WASM サイズ: chematic **504 KB** · RDKit.js ~30 MB · Indigo WASM ~40 MB

---

## chematic を使うべき場面

**chematic が適している場合：**

- ブラウザで化学計算を動かしたい（WASM、504 KB、サーバー不要）
- C++ ツールチェーンなしの Pure Rust スタックが必要
- `pip install rdkit` が困難な環境（Cloudflare Workers、Lambda、組み込み）にデプロイする
- AI エージェントを構築し、ネイティブな MCP ツール統合が必要
- バッチ処理で高スループットが必要（ECFP4: RDKit の 5〜14 倍高速）
- `pip install chematic` がどこでも動くシンプルさを求めている

**RDKit が適している場合：**

- 20 年以上の実績と最大のエコシステム互換性が必要
- ML 補助のコンフォーマー生成が必要（RDKit の ETKDGv3 の方が 3D 品質が高い）
- `native-inchi` feature を有効にせずビット完全な標準 InChI が必要
- RDKit Python API 向けのコミュニティプラグインに依存している

---

## クイックスタート

### インストール

```bash
pip install chematic  # C/C++ コンパイラ不要
```

```python
import chematic

mol = chematic.from_smiles("CC(=O)Oc1ccccc1C(=O)O")  # アスピリン

# Jupyter では mol とセルに書くだけで 2D 構造が自動表示される
mol

# 70+ 記述子をプロパティで取得
print(mol.mw, mol.logp, mol.tpsa)            # 180.16  1.31  63.6
print(mol.lipinski_passes, mol.pains_passes)  # True   True

# サブ構造検索
mol.has_substructure("[OH]")    # True
mol.find_matches("[CX3](=O)O")  # → [[1, 2, 3], [7, 8, 9]]

# バッチ処理（並列、numpy 対応）
fps = chematic.bulk.ecfp4(["CCO", "c1ccccc1"])  # (2, 2048) uint8

# ワンライナーで DataFrame
df = chematic.descriptors_df(["CCO", "c1ccccc1", "CC(=O)O"])
df[["mw", "logp", "tpsa", "qed"]]
```

Rust・JavaScript の詳細な使用例は [ドキュメント](https://kent-tokyo.github.io/chematic/) を参照してください。

---

## AI / LLM 開発者向け

chematic は **MCP（Model Context Protocol）サーバー**をネイティブに搭載した、初のケモインフォマティクスライブラリです。

```json
// Claude Desktop (~/.config/claude/claude_desktop_config.json)
{
  "mcpServers": {
    "chematic": { "command": "chematic-mcp" }
  }
}
```

MCP 対応エージェントから呼び出せる 15 の化学ツール：

| ツール | 機能 |
|---|---|
| `name_to_smiles` | 化合物名（"アスピリン"、"カフェイン"…）を PubChem 経由で SMILES に変換 |
| `calc_properties` | MW、LogP、TPSA、HBA/HBD、QED、SA Score、pKa、ADMET |
| `smarts_match` | 部分構造検索 |
| `pains_check` / `brenk_check` | アッセイ干渉・反応性フラグ付け |
| `generate_3d` | 3D 座標生成（ETKDG + MMFF94） |
| `find_mcs` | 最大共通部分構造 |
| その他 9 ツール | `ecfp4`、`tanimoto`、`canonical_smiles`、`admet_profile`、`boiled_egg`、`sa_score`、`lipinski_check`… |

---

## Pure Rust の理由

### 速い

Rust のゼロコスト抽象化と所有権モデルはオーバーヘッドをソースレベルで排除します。
chematic の ECFP4 フィンガープリントバッチは **3.6 µs/mol** — 同じハードウェアで
RDKit Python API の 5〜14× 高速。GIL なし、インタープリタオーバーヘッドなし、
`_sys` クレート内の FFI 呼び出しコストなし。

### 安全

デフォルトの依存ツリー全体で、15,000 行以上の Rust コードに **`unsafe` ブロックは約 6 個**のみ。
C++ のヒープ破壊なし。不正な SMILES 入力によるセグメンテーション違反なし。
`-sys` クレートによるプラットフォーム固有のビルド失敗なし。
コンパイラがすべての呼び出し箇所でメモリ安全性を保証します。

> `native-inchi` feature は唯一の opt-in 例外 — ビット完全一致の標準 InChI 用に
> IUPAC InChI C ライブラリ (v1.07.5) を vendored でリンクします。他の全クレートは FFI フリーのまま。

### どこでも動く

Pure Rust は Emscripten・`cmake`・`clang` なしで `wasm32-unknown-unknown` にネイティブでコンパイルされます。
npm パッケージ `@kent-tokyo/chematic` は **504 KB gzip** — RDKit.js の 60 分の 1。
1 つのコードベースが Linux・macOS・Windows・あらゆるブラウザで動作します。

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

<details>
<summary>全機能比較（30 以上の機能）</summary>

| 観点 | **chematic** | RDKit (rdkit-sys) | OpenBabel FFI | RDKit.js (WASM) |
|---|---|---|---|---|
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


</details>
---

## JavaScript / TypeScript（WebAssembly）

**504 KB gzip — RDKit.js の 60 分の 1。** Emscripten・cmake 不要。ブラウザ・Node.js どちらでも動作。

```sh
npm install @kent-tokyo/chematic
```

```js
import init, { parse_smiles, get_descriptors_json, tanimoto_ecfp4,
               generate_3d_minimized_pdb, enumerate_stereo_isomers_json,
               maxmin_picks_ecfp4_json } from '@kent-tokyo/chematic';

await init();

const mol = parse_smiles('CC(=O)Oc1ccccc1C(=O)O'); // アスピリン
console.log(mol.molecular_weight(), mol.qed(), mol.lipinski_passes());

const desc = JSON.parse(get_descriptors_json(mol));  // 全記述子を一括取得
const caffeine = parse_smiles('Cn1cnc2c1c(=O)n(c(=O)n2C)C');
console.log(tanimoto_ecfp4(mol, caffeine));           // ECFP4 類似度: 0.26

const pdb = generate_3d_minimized_pdb(mol);           // 3D 座標生成
const isomers = JSON.parse(enumerate_stereo_isomers_json(parse_smiles('C(F)(Cl)Br')));
const picks = JSON.parse(maxmin_picks_ecfp4_json('["CC","c1ccccc1","CCO","CCCC"]', 2));
```

130+ のエクスポート関数が記述子・フィンガープリント・3D・反応・多様性選択・SDF を網羅。
全エクスポートは [WASM API リファレンス](https://kent-tokyo.github.io/chematic/) を参照。
---

## クレート一覧

| クレート               | 説明                                                                                                                                      | テスト数 |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `chematic-core`        | Atom, Bond, Molecule, Element, ケクレ化（依存ゼロ）；ミュータブル API・`fragments`・`validate_valence`・`formula_with_isotopes`・`StereoGroup`/`StereoGroupKind` | 48      |
| `chematic-smiles`      | OpenSMILES パーサー、ライター、正規 SMILES、**CXSMILES メタデータ対応**                                                                  | 57      |
| `chematic-perception`  | SSSR、Hückel 芳香族性 + 反芳香族性（4n+2 則）、`apply_aromaticity`・`aromatize`・`kekulize_inplace`・`assign_stereo_from_2d`・`assign_ez_from_2d`・`cip_ez_descriptor` | 34      |
| `chematic-mol`         | MOL/SDF V2000+V3000（R/W、2D 座標付き）、CML（R/W）、CDXML（R）；`SdfRecord`（coords+props）、MDL RXN V2000 読み書き；V3000 ステレオグループ COLLECTION R/W | 63      |
| `chematic-depict`      | 2D SVG（CPK カラー・ハイライト・グリッド）、`detect_crossings`・`render_svg_with_metadata`・反応 SVG；Y座標系ドキュメント整備  | 43      |
| `chematic-chem`        | 70+ 記述子、タウトマー、スキャフォルド、BRICS、QED、標準化；**pKa 予測** (15 SMARTS ルール)；**ADMET プロファイル** (BBB/Caco-2/hERG/CYP3A4)；**HBA 99.98% RDKit 一致率**（5,000 分子ベンチマーク）；**TPSA ±1.0 Å² / LogP ±0.3 / HBD 100%** RDKit 一致（175 分子バルク回帰） | 496     |
| `chematic-fp`          | ECFP2/4/6、FCFP4/6、MACCS、TopoPF、AtomPair、Torsion、Layered、Pattern、Pharmacophore、Reaction、**MAP4** (Minervini 2020) — Tanimoto/Dice | 55      |
| `chematic-ff`          | **MMFF94 全 7 エネルギー項** (Halgren 1996)：OOP (117件) + Stretch-Bend (282件)；steepest descent + L-BFGS；DREIDING | 98      |
| `chematic-smarts`      | SMARTS、VF2、MCS；**SmartsCache** (LRU 5–20×)；**named_pattern()** (20 パターン)；**SMARTS 内アトムマップ `:N`** (`[O;D1;H0:3]` 形式 — メタデータとして保存、マッチング条件には不使用) | 137     |
| `chematic-3d`          | 3D 座標生成、ETKDG KB (40 パターン、adaptive noise)、力場最小化、形状記述子、ConformerEnsemble、PDB/XYZ | 147     |
| `chematic-rxn`         | 反応 SMILES/SMIRKS、`run_reactants`/`run_reactants_strict`；**`retro_disconnect()`** — 60 retro-SMIRKS テンプレート (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ランク付き | 30      |
| `chematic-inchi`       | InChI/InChIKey：純 Rust 近似（WASM 対応）**+ `native-inchi` feature で IUPAC 標準準拠**（C ライブラリ 1.07.5 vendored、ビット完全一致）；**parse_inchi** 読み込み | 28 (+14*)   |
| `chematic-wasm`        | **130+ WASM エクスポート** — npm: `@kent-tokyo/chematic` v0.4.13；**pKa/ADMET/BBB/Caco-2/hERG/CYP3A4** WASM API | 209     |
| `chematic-iupac`       | ローカル IUPAC 命名（Pure Rust・オフライン）— **25+ 化合物クラス**：アルカン、シクロアルカン、アルコール、アミン、ハロアルカン、ケトン、酸、エステル、アミド、**ピペリジン、モルホリン、ピペラジン、ナフタレン、スルフィド** | 45      |
| `chematic-mcp`         | **MCP (Model Context Protocol) サーバー** — AI エージェント統合；**15 ツール**：parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles | 28      |
| `chematic`             | フィーチャーフラグ付きアンブレラクレート（統合クレート）                                                                                                  | 1       |

```
cargo test --workspace --lib --quiet                                               # 211 ライブラリテスト、全パス
cargo test -p chematic-inchi --features native-inchi --test standard_inchi         # +14 IUPAC 標準 InChI 統合テスト
```

---

## 最近の開発（v0.4.x）

**v0.4.19**（2026-06-23）: **PDF/EPS 出力、ChemicalJSON、新記述子、WASM −38.5%**
- `chematic-depict`: `depict_pdf()` / `depict_eps()` — PDF・EPS 出力（Pure Rust、外部ツール不要）
- `chematic-mol`: **ChemicalJSON** — `parse_cjson()` / `write_cjson()` で Avogadro2 / MolSSI 相互運用
- `chematic-chem`: 新記述子 4 件 — `schultz_mti()`, `gutman_mti()`, `vabc()`（Bondi vdW 体積）, `gravitational_index()`
- `chematic-3d`: **Spectrophores** 3D フィンガープリント（ファーマコフォアシェルエンコーディング）
- `chematic-py`: `mol.to_pdf()`, `mol.to_eps()`, `mol.to_cjson()`, `from_cjson()`；`bulk.substructure_match(smarts, mols)` 並列 VF2；`estate_all()` / `ring_bundle` を bulk に追加
- **WASM バンドル: 819 → 504 KB gzip（−38.5%）** — `tiny_skia` オプション化、インライン SHA-256、`opt-level="z" lto=true codegen-units=1`

**v0.4.18**（2026-06-23）: **Python API 拡充 + ベンチマーク公開**
- `chematic-py`: **Jupyter 自動表示** — セルに `mol` と書くだけで 2D 構造が表示（`_repr_svg_()` フック）；`mol.has_substructure(smarts)`, `mol.find_matches(smarts)`；`from_smiles_list()`, `descriptors_df()`
- `chematic-chem`: `chi_all()` — Hall-Kier 連結指数 10 本を 1 パスで計算；`cns_mpo_from_parts()`；`pains_passes_and_matches()` / `brenk_passes_and_matches()` — 1 回のスキャンでフラグと名前を返す
- ドキュメント: ベンチマークページ追加（ECFP4 5–14× 高速、5 000 分子コーパス 100% 一致）

**v0.4.16–v0.4.17**（2026-06-22–23）: **SSSR 共有パフォーマンス改善**
- `chematic-smarts`: `find_matches_with_rings()` — バッチ全パターンで `RingSet` を 1 回だけ計算して共有
- `chematic-chem`: Crippen 117 SSSR → 1/呼び出し；PAINS ~480 → 1；QED 113 → 1；pKa 42 → 1；新 API `logp_and_mr()`, `logd_from_logp()`, `pka_both()`
- `chematic-fp`: MHFP 増分 BFS — 分子あたり 3N → N BFS 操作（radius=2 時）

**v0.4.15**（2026-06-21）: **TPSA 校正 + 反応 E/Z ステレオ**
- `chematic-chem`: TPSA ±0.1 Å² 校正 — **HBA 100%、HBD 100%、芳香環数 100%**（5 000 分子コーパス）；TPSA 86.7% → 93.3%（5 000 分子）、175 分子薬様セットで 100%
- `chematic-rxn`: `run_reactants` に E/Z 幾何フィルタ追加 — SMIRKS `/`/`\` に基づく `smirks_ez_stereo_ok()` / `ez_stereo_outward()`

**v0.4.14**（2026-06-21）: **トポロジー記述子 + ステレオ正確性**
- `chematic-chem`: トポロジー記述子 8 件 — `petitjean_index()`, `graph_eccentricities()`, `graph_diameter()`, `graph_radius()`, `eccentric_connectivity_index()`, `hosoya_index()`, `moran_autocorr()`, `geary_autocorr()`
- `chematic-3d`: GETAWAY HATS-matrix 19 次元；`whim_getaway_combined()` が 29 次元に
- `chematic-smiles`: アレン累積二重結合ステレオ `C=C=C` `@`/`@@` — ラウンドトリップ安定
- `chematic-smarts`: `[kN]` 環サイズプリミティブ；VF2 クエリ原子数 > 対象時の早期終了
- `chematic-rxn`: パリティ対応 SMIRKS キラリティマッチ；product bracket クリーンアップ

**v0.4.13**（2026-06-21）: **テンプレート逆合成 + 記述子修正**
- `chematic-rxn`: `retro_disconnect()` — 60 retro-SMIRKS テンプレート（AmideBond / Ester / Ether / CNBond / CCBond / CSBond）、SA Score ランク付き；Python `mol.retro_disconnect(reaction_class=...)`
- `chematic-3d`: ETKDG トーション KB 28 → 40 パターン；adaptive noise
- `chematic-chem`: `hbd_count()` に S-H（チオール）追加；TPSA nitro-N / 芳香族オキシドブリッジ / Kekulé-N 修正

**v0.4.9–v0.4.12**（2026-06-19–21）: **AutoDock、UFF、SMARTS アトムマップ、環認識**
- `chematic-mol`: AutoDock PDBQT 読み書き；`write_sdf_with_charges`（部分電荷）
- `chematic-ff`: 金属・有機金属向け UFF 力場（Zn、Fe、Cu…）
- `chematic-smarts`: SMARTS アトムマップ `:N`（`[O;D1;H0:3]` 形式、メタデータとして保存）
- `chematic-perception`: 多縮合芳香環の `augmented_ring_set` 反復更新（bench5k 222 件全修正）
- MCP: 15 番目のツール `name_to_smiles`（PubChem REST プロキシ）

**v0.4.5–v0.4.7**（2026-06-19）: **ケクレ化 blossom + BOILED-Egg + InChI E/Z**
- Edmonds' blossom アルゴリズム導入（128 → 2 失敗）；InChI `/b` E/Z レイヤー；BOILED-Egg + Python/WASM バインディング

**v0.4.0–v0.4.4**（2026-06-17–18）: **PyO3 Python バインディング + native-inchi**
- `chematic-py`: PyO3/maturin バインディング — `from_smiles()`, `Mol.aromatic_ring_count`, `Mol.descriptors()`
- `native-inchi` feature: IUPAC 標準 InChI（vendored C lib v1.07.5）
- HBA 書き直し: RDKit と 99.98% 一致（5,000 分子ベンチマーク）

---

## 既知の制限事項

- **ケクレ化**: 5,000 分子中 2 件のみ失敗 — ホウ素芳香環（`b1ccccn1`）と `[H][H]`。`KekuleError` を明示的に返し、無音の誤出力は生じない。
- **芳香族性モデル**: Hückel 4n+2 則を各 SSSR 環に独立適用（RDKit は縮合環電子非局在化モデルを使用）。N-ヘテロ環で差異あり。5,000 分子コーパスの現状: HBA/HBD/芳香環数 **100%**、TPSA **93.3%**（±0.1 Å²）。

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

---

chematic が役に立ったら、[GitHub スター](https://github.com/kent-tokyo/chematic)をいただけると他の方への発見につながります。
