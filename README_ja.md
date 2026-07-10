# chematic

[English](README.md) | [中文](README_zh.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![docs.rs](https://docs.rs/chematic/badge.svg)](https://docs.rs/chematic)

![Pure Rust](https://img.shields.io/badge/Pure%20Rust-zero%20C%2B%2B-orange?logo=rust)
![WASM](https://img.shields.io/badge/WASM-504%20KB-blueviolet?logo=webassembly)
![MCP](https://img.shields.io/badge/MCP-agent%20ready-purple)
[![ライセンス](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)
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

**機能の成熟度（早見表）：**

| 機能 | ステータス |
|---|---|
| SMILES / SMARTS / フィンガープリント / 記述子 | 安定 |
| 3D 配座生成（DG + MMFF94） | 実験的 |
| pKa / ADMET | ルールベーススクリーニング（臨床用途不可） |
| IUPAC 名生成 | 部分実装（25+ クラス） |
| Pure-Rust InChI | 近似値（完全精度には `native-inchi` feature を有効化） |

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
- ML 補助のトーション補正による出版品質の 3D 構造が必要（RDKit の ETKDGv3）
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

# 190+ 記述子をプロパティで取得
print(mol.mw, mol.logp, mol.tpsa)            # 180.16  1.31  63.6
print(mol.lipinski_passes, mol.pains_passes)  # True   True

# サブ構造検索
mol.has_substructure("[OH]")    # True
mol.find_matches("[CX3](=O)O")  # → [[1, 2, 3], [7, 8, 9]]

# 自然言語サマリー（LLM / MCP エージェント向け）
print(mol.describe())
# → "Molecular weight 180.2 Da, formula C9H8O4. LogP 1.31 (mildly lipophilic)..."

# 2分子の構造差分
ibuprofen = chematic.from_smiles("CC(C)Cc1ccc(CC(C)C(=O)O)cc1")
d = mol.diff(ibuprofen)  # {"summary": "+C7, -O2. ΔLogP +2.75 ...", "delta_mw": 66.1, ...}

# バッチ処理（並列、numpy 対応）
fps = chematic.bulk.ecfp4(["CCO", "c1ccccc1"])  # (2, 2048) uint8

# ワンライナーで DataFrame
df = chematic.descriptors_df(["CCO", "c1ccccc1", "CC(=O)O"])
df[["mw", "logp", "tpsa", "qed"]]
```

Rust・JavaScript の詳細な使用例は [ドキュメント](https://kent-tokyo.github.io/chematic/) を参照してください。

### 動作確認

```python
import chematic
chematic.doctor()
# chematic v0.4.29
# Python 3.12.x  |  darwin arm64
#
# Descriptor accuracy (benchmark 2026-06, v0.4.29 vs RDKit 2026.03.3):
#   MW / HBA / HBD / ARC  100%   (4,999-mol ChEMBL subset)
#   TPSA                  98.1%
#   LogP (Crippen)        ~99%
# ...
```

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
| 分子記述子                                  | **190+（MQN×42、BCUT2D、ESOL、LogD、XLogP3、BOILED-Egg 含む）**  | 〜30               | 〜20          | 〜30             |
| **MAP4 フィンガープリント**                 | **あり**（Minervini 2020）                 | なし（外部pkg）    | なし          | なし             |
| MMFF94 全 7 エネルギー項                    | **あり**                                   | あり               | あり          | なし             |
| 3D 座標生成                                 | あり（DG + MMFF94/DREIDING + L-BFGS）      | あり（ETKDG）      | あり          | あり             |
| 多様性選択（MaxMin/Butina）                 | **あり**                                   | あり               | なし          | なし             |
| InChI / InChIKey                            | **あり** — 純 Rust（デフォルト）+ **IUPAC 準拠**（`native-inchi`）| C ライブラリ必要 | C ライブラリ必要 | C ライブラリ必要 |
| **pKa 予測**                                | **あり（15 SMARTS ルール）**               | なし               | なし          | なし             |
| **ADMET プロファイル + BOILED-Egg**         | **あり**                                   | 一部               | なし          | 一部             |
| **MCP サーバー（AI エージェント API）**     | **あり — 20 ツール（Name→SMILES 含む）**  | なし               | なし          | なし             |
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
| `chematic-core`        | Atom, Bond, Molecule, Element, ケクレ化（依存ゼロ）；ミュータブル API・`fragments`・`validate_valence`・`formula_with_isotopes`・`StereoGroup`/`StereoGroupKind` | 71      |
| `chematic-smiles`      | OpenSMILES パーサー、ライター、正規 SMILES、**CXSMILES メタデータ対応**                                                                  | 109      |
| `chematic-perception`  | SSSR、Hückel 芳香族性 + 反芳香族性（4n+2 則）、`apply_aromaticity`・`aromatize`・`kekulize_inplace`・`assign_stereo_from_2d`・`assign_ez_from_2d`・`cip_ez_descriptor` | 101      |
| `chematic-mol`         | MOL/SDF V2000+V3000（R/W、2D 座標付き）、CML（R/W）、CDXML（R）；`SdfRecord`（coords+props）、MDL RXN V2000 読み書き；V3000 ステレオグループ COLLECTION R/W | 130      |
| `chematic-depict`      | 2D SVG（CPK カラー・ハイライト・グリッド）、`detect_crossings`・`render_svg_with_metadata`・反応 SVG；Y座標系ドキュメント整備  | 64      |
| `chematic-chem`        | 190+ 記述子値（71 関数）、タウトマー、スキャフォルド、BRICS、QED、標準化；**pKa 予測** (15 SMARTS ルール)；**ADMET プロファイル** (BBB/Caco-2/hERG/CYP3A4)；**HBA 99.98% RDKit 一致率**（4,999 分子 ChEMBL ベンチマーク）；**TPSA ±0.1 Å² 98.1% / LogP ±0.01 96.5% / HBD 100%** RDKit 一致 | 662     |
| `chematic-fp`          | ECFP2/4/6、FCFP4/6、MACCS、TopoPF、AtomPair、Torsion、Layered、Pattern、Pharmacophore、Reaction、**MAP4** (Minervini 2020) — Tanimoto/Dice | 185      |
| `chematic-ff`          | **MMFF94 全 7 エネルギー項** (Halgren 1996)：OOP (117件) + Stretch-Bend (282件)；steepest descent + L-BFGS；DREIDING | 98      |
| `chematic-smarts`      | SMARTS、VF2、MCS；**SmartsCache** (LRU 5–20×)；**named_pattern()** (20 パターン)；**SMARTS 内アトムマップ `:N`** (`[O;D1;H0:3]` 形式 — メタデータとして保存、マッチング条件には不使用) | 142     |
| `chematic-3d`          | 3D 座標生成、ETKDG KB (40 パターン、adaptive noise)、力場最小化、形状記述子、ConformerEnsemble、PDB/XYZ | 265     |
| `chematic-rxn`         | 反応 SMILES/SMIRKS、`run_reactants`/`run_reactants_strict`；**`retro_disconnect()`** — 60 retro-SMIRKS テンプレート (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ランク付き | 137      |
| `chematic-inchi`       | InChI/InChIKey：純 Rust 近似（WASM 対応）**+ `native-inchi` feature で IUPAC 標準準拠**（C ライブラリ 1.07.5 vendored、ビット完全一致）；**parse_inchi** 読み込み | 96 (+16*)   |
| `chematic-wasm`        | **130+ WASM エクスポート** — npm: `@kent-tokyo/chematic` v0.4.19；**pKa/ADMET/BBB/Caco-2/hERG/CYP3A4** WASM API | 211     |
| `chematic-iupac`       | ローカル IUPAC 命名（Pure Rust・オフライン）— **25+ 化合物クラス**：アルカン、シクロアルカン、アルコール、アミン、ハロアルカン、ケトン、酸、エステル、アミド、**ピペリジン、モルホリン、ピペラジン、ナフタレン、スルフィド** | 47      |
| `chematic-mcp`         | **MCP (Model Context Protocol) サーバー** — AI エージェント統合；**20 ツール**：parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles, retrosynthesis, smiles_to_moljson, moljson_to_smiles, representation_router, molecule_context_pack | 31      |
| `chematic`             | フィーチャーフラグ付きアンブレラクレート（統合クレート）                                                                                                  | 1       |

```
cargo test --workspace --lib --quiet                                               # 2,366 ライブラリテスト、全パス
cargo test -p chematic-inchi --features native-inchi --test standard_inchi         # +16 IUPAC 標準 InChI 統合テスト
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
- ドキュメント: ベンチマークページ追加（ECFP4 5–14× 高速、4,999 分子 ChEMBL コーパス 100% 一致）

**v0.4.16–v0.4.17**（2026-06-22–23）: **SSSR 共有パフォーマンス改善**
- `chematic-smarts`: `find_matches_with_rings()` — バッチ全パターンで `RingSet` を 1 回だけ計算して共有
- `chematic-chem`: Crippen 117 SSSR → 1/呼び出し；PAINS ~480 → 1；QED 113 → 1；pKa 42 → 1；新 API `logp_and_mr()`, `logd_from_logp()`, `pka_both()`
- `chematic-fp`: MHFP 増分 BFS — 分子あたり 3N → N BFS 操作（radius=2 時）

**v0.4.15**（2026-06-21）: **TPSA 校正 + 反応 E/Z ステレオ**
- `chematic-chem`: TPSA ±0.1 Å² 校正 — **HBA 100%、HBD 100%、芳香環数 100%**（4,999 分子 ChEMBL コーパス）；TPSA 86.7% → 93.3%（4,999 分子）、175 分子薬様セットで 100%
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
- HBA 書き直し: RDKit と 99.98% 一致（4,999 分子 ChEMBL ベンチマーク）

---

## 既知の制限事項

- **カノニカルSMILESは完全には自己安定でない（最優先の既知課題）**: `canonical_smiles(parse(x))` は `x` の入力走査順によらず不動点であるべきです — これはカノニカルSMILESが存在する理由そのものです（重複排除・ハッシュ・キャッシュ/DBキー）。5,000分子ChEMBLサブセット・worst-of-10走査・単一の制御された実行で測定: 生の `.smiles` 自己安定性は**86.02%**（13.98%不安定、699/5000）。同一サンプルで原因を切り分け: ステレオのみを除去（`canonical_smiles_mode("nostereo")`）すると自己安定性は**99.68%**（0.32%不安定）まで上昇 — 検証した生の不安定例は全て `@`/`@@`/`/`/`\` の位置のみが異なり骨格の違いはなく、13.98ポイントのうち約13.6ポイントはステレオ記述子の非カノニカル化（ライターが入力の原子順序に依存せず `@`/`@@`・結合方向を再計算していない）が原因です。残る0.8%はより深いカノニカルランクのタイブレーク問題で、根本原因は未特定です。この問題はSSSR修正とは独立した既存のバグであり、「将来課題」として先送りしません — `canonical_smiles()` の出力を重複排除やキャッシュキーとして使うあらゆるコードに影響します。下記のMurckoスキャフォールドの残存不安定性も同じ根本原因です。
- **環知覚（SSSR）は非決定的・非最小だったが修正済み**: 旧 `find_sssr` は単一の全域木から非木辺ごとに基本閉路を1つずつ生成していたため、木の形によっては最小でない環（例: ナフタレン `c1ccc2ccccc2c1` が決定的に `[6,10]` を返す。正しくは `[6,6]`）を返していました。現在の `find_sssr` は Horton アルゴリズム（全頂点×全辺の最短路木から候補閉路を生成、O(V·E) 候補、決定性のためのカノニカルランクによるタイブレーク）を用いており、真に最小重み・決定的な基底を返します。5,000分子ChEMBLサブセット・worst-of-10走査で測定: 自己安定性 **100%**（旧: 50.6%）、単一パースでのRDKitとの環サイズ一致率 **98.9%**（旧: 72.4%）— 残る約1.1%の差はRDKit自身の `GetSymmSSSR` が対称縮合環系（例: キュバン、μ=5に対しRDKitは6環を返す）でトポロジー的最小数より多い環を正当に返すことによるもので、chematic側のバグではありません（完全な対称化＝Vismara relevant cyclesは将来課題）。下流への効果（同コーパス）: 環サイズSMARTS `[r5]`/`[r6]` **0%不安定**（旧: 29〜55%）、`NumAromaticRings` **0%**（旧: 約4%）、`RingCount`/MW/TPSA/HBA/HBD/LogP/MRは元々0%のまま。旧SSSRのバグが偶然、別の未解決な芳香族性バグを覆い隠していた2件については下記の芳香族性モデルの項を参照してください。詳細: `scripts/ringinfo_parity.py`。
- **Murckoスキャフォールド: 環トポロジーは安定化、SMILES文字列出力は上記カノニカルSMILESの問題を引き継ぐ**: 以前報告した「100%不安定」自体が測定ハーネスのバグでした（`Mol` オブジェクトをPythonの同一性で比較しており、実際の結果に関わらず常に「不安定」と判定していた）。このスクリプトバグは修正済みです（`scripts/ring_collateral_damage.py`）。正しく測定すると、5,000分子・worst-of-10で `scaffold().smiles` の生文字列比較は依然として**約45%**不安定で、上記と同じステレオ/結合方向のカノニカル化の不備が主因です。正規化（`apply_aromaticity().canonical_smiles_mode("nostereo")`）後は**0.8%**の残差があります — 直接の切り分け（不一致となった変異体間で `scaffold().heavy_atoms`＝原子数を比較し、検証した40件全てで原子数が一致）により、Murckoの環選択バグではなく上記と**同じ一般的なカノニカルSMILESライターのランキング問題**であることを確認済みです: 異なる文字列は同一の原子集合が異なる形（環closureの番号付け違い、ケージ構造の分岐順序違い）でシリアライズされたものです。`scaffold()` は正しい環系を確実に抽出します。走査順が異なる入力間で文字列として比較したい場合は、生の `.smiles` ではなく `mol.apply_aromaticity().canonical_smiles_mode("nostereo")` を使ってください。
- **芳香族性モデル**: Hückel 4n+2 則を各 SSSR 環に独立適用（RDKit は縮合環電子非局在化モデルを使用）。N-ヘテロ環で差異あり。4,999 分子 ChEMBL コーパスの現状: HBA/HBD/芳香環数 **100%**、TPSA **98.1%**（±0.1 Å²）。Kekulized入力でのworst-of-10芳香族性パリティ: **96.3%**（`scripts/aromaticity_atom_parity.py`）— 上記SSSR修正後も数値はビット単位で不変（SSSRのバグと芳香族性のギャップが独立していたことを確認）。芳香族性ギャップの根本原因は別途 `aromatic_context` バイパス機構にあり、未修正です。SSSR修正により2分子（アズレン、プリン）が明確に退行したことが判明しています — 旧来の壊れたSSSRが、これら非交互環系・架橋頭部を多く含む構造に対して `aromatic_context` のバグを偶然覆い隠していたためです。この2分子は5,000分子の測定コーパスに**そもそも含まれていません**（直接検索で確認済み — ChEMBL由来の薬様コーパスには裸のアズレン・プリンは含まれない）。したがって96.3%という数値が不変なのはこの2分子を「見ていない」からであり、「影響がゼロ」だからではありません。両分子は根本原因をコード内に記載した上で `chematic-perception` のテストスイートに `#[ignore]` 付き既知回帰として記録済みで、`aromatic_context` 修正を待っています。
- **TPSA 残差**: 1.9% はアジド基・マクロライドラクトン・ホスファゼン等の特殊化学構造に集中。

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
