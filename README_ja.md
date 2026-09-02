# chematic

[English](README.md) | [中文](README_zh.md)

[![CI](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml/badge.svg)](https://github.com/kent-tokyo/chematic/actions/workflows/ci.yml)
[![PyPI](https://img.shields.io/pypi/v/chematic?logo=pypi)](https://pypi.org/project/chematic/)
[![crates.io](https://img.shields.io/crates/v/chematic?logo=rust)](https://crates.io/crates/chematic)
[![npm](https://img.shields.io/npm/v/@kent-tokyo/chematic?logo=npm)](https://www.npmjs.com/package/@kent-tokyo/chematic)
[![docs.rs](https://docs.rs/chematic/badge.svg)](https://docs.rs/chematic)
[![ライセンス](https://img.shields.io/badge/license-MIT%20OR%20Apache--2.0-blue)](LICENSE-MIT)

[Open in Colab](https://colab.research.google.com/github/kent-tokyo/chematic/blob/main/notebooks/quickstart.ipynb)

Python・Rust・ブラウザ向けケモインフォマティクスライブラリ。

**デフォルトで速く、設計で安全なケモインフォマティクス。**  
Pure Rust · C/C++ ゼロ · Python · WebAssembly · [公式サイト](https://chematic.io/) · [ライブデモ](https://kent-tokyo.github.io/chematic/playground/)

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **導入方法** | `pip install chematic` | `pip install rdkit`（公式prebuiltホイール）または conda | `npm install @rdkit/rdkit`、Python バインディングなし |
| **ブラウザ向けバンドル** | **2.94 MB raw / 1.10 MB gzip** | 該当なし（Python/C++ライブラリ） | 6.91 MB raw* |
| **バッチ FP 速度** | **~78 µs/mol**（2–3× 高速） | ~160–235 µs/mol | — |
| **メモリ安全性** | コンパイラが保証（Rust） | C++ | C++ |
| **ソースビルド** | `cargo build` のみ | cmake + clang + Boost | Emscripten SDK |

\* RDKit.js の gzip転送時サイズは未計測のため、rawサイズ同士で比較している。RDKit.js は
現在メンテナ移行中(詳細は同リポジトリを参照)。

すべての数値は再現可能です — [ベンチマーク詳細](https://kent-tokyo.github.io/chematic/benchmark/)を参照。  
WASM サイズ(raw、2026-08-21計測、`wasm-pack build --target web --release` + `wasm-opt -O3`
のクリーンビルド、commit `ef7dc25`): chematic **2.94 MB**(**1.10 MB gzip**) · RDKit.js **6.91 MB**
(`@rdkit/rdkit@2025.3.4-1.0.0`の`RDKit_minimal.wasm`、unpkg.com で確認) · Indigo(Ketcher向けビルド)
**11.24 MB**(`indigo-ketcher@1.45.1`のメイン`.wasm`、jsDelivr で確認) — chematic の raw WASM
バイナリは現在、RDKit.js よりおよそ2.3倍、Indigo の Ketcher向けビルドよりおよそ3.8倍小さい
(raw同士の比較)。

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

- ブラウザで化学計算を動かしたい（WASM、1.10 MB gzip、サーバー不要）
- C++ ツールチェーンなしの Pure Rust スタックが必要
- RDKit の導入が困難・非対応な環境（Cloudflare Workers、Lambda、組み込み）にデプロイする
  (RDKit 自体も公式`pip install rdkit`ホイールを提供しているが、通常のCPython環境を前提とする)
- AI エージェントを構築し、ネイティブな MCP ツール統合が必要
- バッチ処理で高スループットが必要（ECFP4: RDKit の 2〜3 倍高速、Rayon 並列）
- `pip install chematic` がどこでも動くシンプルさを求めている

**RDKit が適している場合：**

- 20 年以上の実績と最大のエコシステム互換性が必要
- ML 補助のトーション補正による出版品質の 3D 構造が必要（RDKit の ETKDGv3）
- `native-inchi` feature を有効にせずビット完全な標準 InChI が必要
- RDKit Python API 向けのコミュニティプラグインに依存している

---

## インターフェースを選ぶ

- [Rust](#クイックスタート)
- [Python](#クイックスタート)
- [WebAssembly / Node.js](README.md#javascript--typescript-webassembly)
- [材料・シミュレーション用フォーマット](docs/format-capabilities.md) — mmCIF, PQR, QCSchema, ORCA, Gaussian Cube, OpenDX, LAMMPS
- [RDKit からの移行](docs/rdkit-migration.md) — 機能ごとの Supported / Partial / Not-supported 対応表

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
# chematic v0.23.0
# Python 3.12.x  |  darwin arm64
#
# Descriptor accuracy (benchmark 2026-06, v0.4.22 vs RDKit 2026.03.3 --
# descriptor calculation paths unchanged through v0.8.0, not re-measured since):
#   MW / HBA / HBD / ARC  100%   (4,999-mol ChEMBL subset)
#   TPSA                  100%
#   LogP (Crippen)        100%*  (max Δ = 1.1×10⁻¹³)
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

MCP 対応エージェントから呼び出せる 20 の化学ツール（全リストは [`chematic-mcp` README](crates/chematic-mcp/README.md) 参照）：

| ツール | 機能 |
|---|---|
| `name_to_smiles` | 化合物名（"アスピリン"、"カフェイン"…）を PubChem 経由で SMILES に変換(外部通信を行う唯一のツール) |
| `calc_properties` | MW、exact mass、Crippen LogP、TPSA、HBD、HBA、rotatable bonds、QED |
| `smarts_match` | 部分構造検索 |
| `pains_check` / `brenk_check` | アッセイ干渉・反応性フラグ付け |
| `generate_3d` | rule-based配置 + DREIDING力場最小化による3D座標生成 |
| `find_mcs` | 最大共通部分構造 |
| その他 13 ツール | `ecfp4`、`tanimoto`、`canonical_smiles`、`admet_profile`、`boiled_egg`、`sa_score`、`lipinski_check`、`retrosynthesis`、`smiles_to_moljson`、`moljson_to_smiles`、`representation_router`、`molecule_context_pack`、`parse_smiles` |

**Transport**: stdio（標準入出力経由の JSON-RPC 2.0）のみ。ローカルプロセスとして動作し、公開された Remote MCP エンドポイント・認証・公開サービス SLA は存在しない。remote 対応のリファクタは検討中だが未実装。

**Protocol**: 同一の stdio コネクション上でレガシー（`2024-11-05` 形式の `initialize` ハンドシェイク）と、MCP `2026-07-28` のステートレス方言（`server/discover`、リクエストごとの `_meta`、キャッシュ可能な `tools/list`、`structuredContent`）の両方に対応。詳細は [`chematic-mcp` README](crates/chematic-mcp/README.md#protocol-eras) を参照。Remote HTTP・OAuth・Tasks 拡張・MCP Apps は引き続き未対応。

---

## Pure Rust の理由

### 速い

Rust のゼロコスト抽象化と所有権モデルはオーバーヘッドをソースレベルで排除します。
chematic の ECFP4 フィンガープリントバッチは多様な分子コーパスで **~78 µs/mol** — 同じハードウェアで
RDKit Python API の 2〜3× 高速（全 CPU コアで Rayon 並列化）。GIL なし、インタープリタオーバーヘッドなし、
`_sys` クレート内の FFI 呼び出しコストなし。

### 安全

chematic 自身の約 180,700 行の Rust コード(tokei計測、全20クレート、2026-08-21時点)には
`unsafe` ブロックが1ファイルの外では**ゼロ**個 — `unsafe {}` 9個 + `unsafe extern "C"` 宣言1個のみで、
すべて opt-in の `native-inchi` FFI 層に限定されています(下記参照)。
C++ のヒープ破壊なし。不正な SMILES 入力によるセグメンテーション違反なし。
`-sys` クレートによるプラットフォーム固有のビルド失敗なし。
コンパイラが chematic 自身が書いたすべての呼び出し箇所でメモリ安全性を保証します。

> `native-inchi` feature は唯一の opt-in 例外 — ビット完全一致の標準 InChI 用に
> IUPAC InChI C ライブラリ (v1.07.5) を vendored でリンクします。他の全クレートは
> FFI フリー・unsafe フリーのまま。この数値は chematic 自身のソースのみを指し、
> 依存ツリー全体ではありません — `depict` feature(SVG/PDF/EPS 描画)はフォント・
> 画像レンダリングスタック(resvg/usvg/rustybuzz/tiny-skia/zune-jpeg)を引き込み、
> これらは unsafe フリーでは**ありません**。実測値は下記の比較表脚注を参照。

### どこでも動く

Pure Rust は Emscripten・`cmake`・`clang` なしで `wasm32-unknown-unknown` にネイティブでコンパイルされます。
npm パッケージ `@kent-tokyo/chematic` は **1.10 MB gzip**(raw 2.94 MB)— RDKit.js の
`RDKit_minimal.wasm`(raw 6.91 MB)と raw同士で比較しておよそ2.3分の1。
1 つのコードベースが Linux・macOS・Windows・あらゆるブラウザで動作します。

---

## 他のケモインフォマティクスライブラリとの比較

| 観点                                        | **chematic**                               | RDKit (rdkit-sys)  | OpenBabel FFI | RDKit.js (WASM)  |
|---------------------------------------------|--------------------------------------------|--------------------|---------------|------------------|
| **C/C++ 依存**                              | **ゼロ（デフォルト）**†                    | 大規模 C++         | 大規模 C++    | C++（Emscripten）|
| **WASM バイナリサイズ**                     | **raw 2.94 MB(gzip 1.10 MB)**              | N/A（WASM 非対応） | N/A           | raw 6.91 MB      |
| **ビルド要件**                              | `cargo build` のみ                         | cmake + clang      | cmake + clang | Emscripten SDK   |
| **Python バインディング**                   | **あり** (`pip install chematic`, PyO3)    | あり（rdkit-sys）  | あり          | なし             |
| unsafe Rust                                 | **自クレートはなし**‡                      | 大規模             | 大規模        | N/A              |
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
| **pKa 予測**                                | **あり（23 SMARTS ルール）**               | なし               | なし          | なし             |
| **ADMET プロファイル + BOILED-Egg**         | **あり**                                   | 一部               | なし          | 一部             |
| **MCP サーバー（AI エージェント API）**     | **あり — 20 ツール（Name→SMILES 含む、stdio のみ）**  | なし               | なし          | なし             |
| IUPAC 名生成                                | **あり（25+ 化合物クラス）**               | なし               | なし          | 一部             |
| メンテナンス（2026）                        | アクティブ                                 | アクティブ         | 最小限        | アクティブ       |

† デフォルトビルドのみ。`native-inchi` feature は opt-in で C コンパイラが必要（IUPAC InChI C ライブラリ v1.07.5 の vendoring）。これは C/C++ FFI 固有の話 — 下記の `depict` feature は純 Rust の描画クレートを引き込むため、unsafe フリーではなくても C コンパイラ依存は追加しません（‡参照）。

‡ chematic 自身の約 180,700 行の Rust コード(tokei計測、2026-08-21): `native-inchi` の9個の FFI ブロックを除き unsafe フリー（上記「安全」参照）— chematic 自身が書いたコードについての実測済みの主張であり、コンパイラによるチェックが一切効かない RDKit/OpenBabel の C++ FFI unsafe とは、たとえ個数が同程度でも種類が根本的に異なります。**依存ツリー全体については成り立ちません**: opt-in の `depict` feature（SVG/PDF/EPS 描画）は resvg/usvg/rustybuzz/tiny-skia/zune-jpeg を引き込み、これらは純 Rust ですが unsafe フリーではありません — 実測（`unsafe fn`/`impl`/`trait`/`{` の出現数）: tiny-skia 151、zune-jpeg 79、rustybuzz 14、image 8、fontdb 3、tiny-skia-path 3(この範囲だけで合計 258)。`chematic-py`（`pip install chematic`）と npm パッケージはどちらも `chematic-depict` に直接依存するため、これは実際の2つのインストール経路の両方に当てはまります。


</details>
---

## JavaScript / TypeScript（WebAssembly）

**1.10 MB gzip — RDKit.js の raw WASM と比べておよそ2.3分の1。** Emscripten・cmake 不要。ブラウザ・Node.js どちらでも動作。

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

218+ のエクスポート関数(+ `MolHandle`/`DepictOptions` のクラスメソッド、2026-08-21計測)が記述子・フィンガープリント・3D・反応・多様性選択・SDF を網羅。
全エクスポートは [WASM API リファレンス](https://kent-tokyo.github.io/chematic/) を参照。
---

## クレート一覧

| クレート               | 説明                                                                                                                                      | テスト数 |
|------------------------|-------------------------------------------------------------------------------------------------------------------------------------------|---------|
| `chematic-core`        | Atom, Bond, Molecule, Element, ケクレ化（依存ゼロ）；ミュータブル API・`fragments`・`validate_valence`・`formula_with_isotopes`・`StereoGroup`/`StereoGroupKind` | 132      |
| `chematic-smiles`      | OpenSMILES パーサー、ライター、正規 SMILES、**CXSMILES メタデータ対応**                                                                  | 202      |
| `chematic-perception`  | SSSR、Hückel 芳香族性 + 反芳香族性（4n+2 則）、`apply_aromaticity`・`aromatize`・`kekulize_inplace`・`assign_stereo_from_2d`・`assign_ez_from_2d`・`cip_ez_descriptor` | 194      |
| `chematic-mol`         | MOL/SDF V2000+V3000（R/W、2D 座標付き）、CML（R/W）、CDXML（R）；`SdfRecord`（coords+props）、MDL RXN V2000 読み書き；V3000 ステレオグループ COLLECTION R/W；**2Dウェッジ/ハッシュのtetrahedral parity + E/Z二重結合方向を読み込み時に自動認識**（`read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics`、型付きopt-in診断）；**PDBx/mmCIF**（R/W、chain/altloc/model/occupancy/B-factor保持、Open Babel本体はread-only）；**PQR**（R/W）；**QCSchema JSON**（`Molecule`/`AtomicInput`/`AtomicResult`、MolSSIスキーマ、Bohr↔Å変換）；**ORCA**（input R/W・未知blockを損失なく保持、output R：final geometry/trajectory/energy/frequencies/termination/convergenceを型付きで取得）；新設の共有`VolumetricGrid`型 + **Gaussian Cube**（R/W、大規模グリッド向けstreaming-*input* `CubeFileReader` — パース後のvoxel配列自体は引き続き全てメモリ上、非直交axes対応、Bohr/Ångström単位を明示タグ化）+ **OpenDX/APBS scalar field**（R/W）— single-dataset限定、multi-dataset Cubeは黙って切り詰めず型付きでreject | 476     |
| `chematic-depict`      | 2D SVG（CPK カラー・ハイライト・グリッド）、`detect_crossings`・`render_svg_with_metadata`・反応 SVG；Y座標系ドキュメント整備  | 75      |
| `chematic-chem`        | 190+ 記述子値（71 関数）、タウトマー、スキャフォルド、BRICS、QED、標準化；**pKa 予測** (23 SMARTS ルール)；**ADMET プロファイル** (BBB/Caco-2/hERG/CYP3A4)；**HBA 100% RDKit 一致率**（4,999 分子 ChEMBL ベンチマーク）；**TPSA ±0.1 Å² 100% / LogP 100%\* / HBD 100%** RDKit 一致 | 724     |
| `chematic-fp`          | ECFP2/4/6、FCFP4/6、MACCS、TopoPF、AtomPair、Torsion、Layered、Pattern、Pharmacophore、Reaction、**MAP4** (Minervini 2020) — Tanimoto/Dice | 266      |
| `chematic-ff`          | **MMFF94 全 7 エネルギー項** (Halgren 1996)：OOP (117件) + Stretch-Bend (282件)；steepest descent + L-BFGS；DREIDING | 198      |
| `chematic-smarts`      | SMARTS、VF2、MCS；**SmartsCache** (LRU 5–20×)；**named_pattern()** (20 パターン)；**SMARTS 内アトムマップ `:N`** (`[O;D1;H0:3]` 形式 — メタデータとして保存、マッチング条件には不使用) | 169     |
| `chematic-3d`          | 3D 座標生成、ETKDG KB (40 パターン、adaptive noise)、力場最小化、形状記述子、ConformerEnsemble、PDB/XYZ | 540     |
| `chematic-rxn`         | 反応 SMILES/SMIRKS、`run_reactants`/`run_reactants_strict`；**`retro_disconnect()`** — 60 retro-SMIRKS テンプレート (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ランク付き | 180      |
| `chematic-inchi`       | InChI/InChIKey：純 Rust 近似（WASM 対応）**+ `native-inchi` feature で IUPAC 標準準拠**（C ライブラリ 1.07.5 vendored、ビット完全一致）；**parse_inchi** 読み込み；**検証付きcanonical SMILES重複排除**（`dedup`モジュール、legacy CIPで未解決の指定済みtetrahedral stereoに対してfail-closed）；**accurate-CIP dedup preflight**（issue #161、legacy CIPで未解決のstereocentreに対して検証capabilityを回復）；**indexed graph relation API**（`compare_indexed_graph_relation`、直交する`GraphStrictness`/`AtomMapPolicy`軸） | 108 (+16*)   |
| `chematic-cip`         | opt-inの高精度CIPエンジン（`assign_cip_accurate_experimental`、階層的digraph、Rules 1a/1b/2/4b/5、RDKit互換MANCUDE分数原子番号）— デフォルトの`assign_cip()`/`CipMode::LegacyFast`は変更なし | —       |
| `chematic-wasm`        | **218+ WASM エクスポート**（+クラスメソッド、2026-08-21計測） — npm: `@kent-tokyo/chematic`（crates.io/PyPIと同期して公開）；**pKa/ADMET/BBB/Caco-2/hERG/CYP3A4** WASM API | 276     |
| `chematic-iupac`       | ローカル IUPAC 命名（Pure Rust・オフライン）— **25+ 化合物クラス**：アルカン、シクロアルカン、アルコール、アミン、ハロアルカン、ケトン、酸、エステル、アミド、**ピペリジン、モルホリン、ピペラジン、ナフタレン、スルフィド** | 56      |
| `chematic-mcp`         | **MCP (Model Context Protocol) サーバー** — AI エージェント統合；**20 ツール**：parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles, retrosynthesis, smiles_to_moljson, moljson_to_smiles, representation_router, molecule_context_pack；dual-eraプロトコル（レガシー`2024-11-05` + モダン`2026-07-28`ステートレス方言）、全20ツールに`structuredContent`/`outputSchema` | 82      |
| `chematic-py`          | PyO3 Python バインディング（`pip install chematic`）；**`PeriodicStructure.from_cif()`/`.from_poscar()`, `Lattice`, `Site`**（周期／結晶構造 — `chematic-crystal`の最初のホスト言語バインディング）；**`from_cif(text, expand_symmetry=True)`** はデフォルトでCIF自身が明示するsymmetry operationsを展開してfull unit cellを生成（`expand_symmetry=False`で非対称単位のみ — space-groupデータベースは無く、名前/番号からのoperation生成も無し） | 300+    |
| `chematic`             | フィーチャーフラグ付きアンブレラクレート（統合クレート）                                                                                                  | 1       |

```
cargo test --workspace --lib --quiet                                               # 3,912 ライブラリテスト、全パス（2026-08-21 時点）
cargo test -p chematic-inchi --features native-inchi --test standard_inchi         # +16 IUPAC 標準 InChI 統合テスト
```

---

## 最近の開発

**v0.23.0**（2026-08-30）: **MCS精度修正（挙動変更）、RDKit互換fingerprintを2種追加、MCS bindingフル対応**
- `chematic-smarts`：**挙動変更** — `find_mcs`のデフォルト`AtomCompare::Elements`が、芳香族性の一致を要求しなくなった。RDKitの同名`rdFMCS.AtomCompare.CompareElements`と完全に一致する挙動（ライブオラクルで確認済み——RDKitは芳香族性を原子側の制約として一切エンコードせず、結合タイプのクエリのみで表現する）。ライブRDKitオラクルとの一致率は、3つの既存corpusで74.6%/68.2%/70.4%から88.4%/88.5%/97.0%へ向上。従来の厳密な元素+芳香族一致を復元する`AtomCompare`モードは現時点で存在しない
- `chematic-py`/`chematic-wasm`/`chematic-mcp`/`chematic-chem`：`find_mcs`の結果再構築が4つのbinding全てでヘテロ原子や芳香族性を静かに失っていた問題を修正（`QueryMolecule`→具体的な`Molecule`への変換が複合atomクエリを正しく展開していなかった）——上記の修正の測定中に発見
- `chematic-fp`/`chematic-py`：`rdkit_rdk_fp`/`rdkit_layered_fp` — RDKit互換の`Chem.RDKFingerprint`/`Chem.LayeredFingerprint`移植版。6種のfingerprint parityシリーズを完了（ライブRDKitオラクルに対し3corpusで100%/100%/99.44%、100%/100%/99.46%のbit-exact一致）
- `chematic-py`/`chematic-wasm`：`McsConfig`/`McsOutcome`の全フィールドを`find_mcs` bindingに公開（`match_charge`/`match_isotope`/`atom_compare`/`bond_compare`/`timeout_ms`等、従来Rust限定だったもの）
- 詳細は`CHANGELOG.md`の`[0.23.0]`section参照

**v0.22.0**（2026-08-29）: **WASM ensemble binding新規追加、canonicalizationハング修正、3員環embedding修正**
- `chematic-wasm`：`chematic_3d::embed_ensemble_v2`向けの新規`embed_ensemble_v2_json` bindingを追加。Python binding（`Mol.conformer_ensemble_v2()`）を`pipeline_v2.rs`既存の規約（camelCase JSONキー、`schemaVersion: 1` envelope）に沿って踏襲——純粋に追加のみ
- `chematic-smiles`：複数の対称領域が同時に未解決な分子で`canonical_smiles`/`canonical_atom_order`がハングしうる問題（issue #421）——automorphism backtracking探索に内部ステップ上限がなかった。常時有効なステップ上限を追加し、超過時は安全側（automorphismと証明できない扱い）にフォールバックするよう修正
- `chematic-3d`：3員環（シクロプロパン/エポキシド/アジリジン/チイラン）がdistance-geometry embedding段階でfail closedしていた問題——リング閉環ペアが同時に「1-3」かつ直接結合している場合に、汎用角度制約が正しい（より厳しい）結合長制約を上書きしていた。既に結合済みの近傍ペアには角度制約をスキップするよう修正。strict-MMFF94 3Dコーパス 252/265 → 263/265
- 詳細は`CHANGELOG.md`の`[0.22.0]`section参照

**v0.21.0**（2026-08-27）: **`McsConfig`の電荷/同位体マッチングと型付きtimeout結果、正確性修正4件**
- `chematic-smarts`：`McsConfig`に`match_charge`/`match_isotope`フィールドを追加（既存の`match_chiral_tag`と同様、デフォルト`false`）。また新規`McsOutcome` enum（`Exhaustive`/`TimedOut`）を`find_mcs_with_config_checked`経由で提供し、timeoutで探索が打ち切られたかを明示的に報告（最適とは限らない結果をexhaustiveなものと区別できずに返す状態を解消）——純粋に追加のみ、`find_mcs`/`find_mcs_with_config`は無変更
- `chematic-smarts`：`find_mcs`のbranch-and-bound探索が不完全だった問題——`grow()`が各探索ノードでfrontierの最初の候補しか試さず、除外して別候補を試す手段がなかったため、真に大きい共通部分構造を見逃すケースがあった（最小反例：`OC(N)N` vs `NC(N)`は真の3原子ではなく2原子を返していた）。標準的なinclude/exclude branch-and-boundで修正
- `chematic-chem`：`disconnect_metals`が金属結合切断後にdative結合由来の形式電荷を中和せず放置していた問題（issue #403）——RDKit同梱のNCI Diversity Setホールドアウトで4,999件中34件が非冪等だった。切断直後に価電子推論でH数を再計算するよう修正。NCIホールドアウト34件→0件、新規11件の金属錯体ホールドアウトも追加
- `chematic-chem`：`normalize_zwitterion`が、移動可能なプロトンがどちらの原子にもない永続的な電荷分離構造（例：diazo-N,N'-dioxide）に対し、負電荷側原子にプロトンを無から生成し分子式を静かに変えていた問題（issue #407）——両原子が実際にプロトンをやり取りできる場合のみ移動するよう修正。開発用corpusの残差は4→1（残る1件は無関係、issue #402/#415参照）
- `chematic-chem`：`canonical_tautomer`が、融合・架橋環系で化学的に不正な過剰原子価窒素を生成しうる問題（issue #415）——2つの芳香族水素移動機構双方に、受理前のkekulization検証を追加
- 詳細は`CHANGELOG.md`の`[0.21.0]`section参照

**v0.20.1**（2026-08-26）: **canonical SMILES／standardizeの正確性修正3件（パッチリリース、破壊的変更なし）**
- `chematic-smiles`：coupled E/Zのcanonicalizationが再canonicalize時に幾何をsilentlyに変えてしまう問題（issue #390）——canonical writerのE/Zマーカー機構における独立した2つのバグ、再現には両方の修正が必要。修正後、実際の290化合物コーパスで検証（idempotency 290/290、独立したRDKit InChIKeyとの一致290/290、修正前の289/290から改善）
- `chematic-chem`：`standardize()`がいくつかの再構築経路でstereoテーブルをsilentlyに失う問題（issue #399）——`standardize.rs`の8関数が素の`MoleculeBuilder`で再構築する際に`stereo_neighbor_order`/`bond_directions`/`stereo_groups`を引き継いでおらず、ring-closing/ring-openingの役割次第で`@`/`@@`がフリップしていた。開発用corpusのstandardize経由idempotencyは615/519→68/60（#392以前の水準に正確に一致）、NCI holdout（未使用の実分子4,999件、一度だけ実行）でステレオ関連の失敗は0件
- `chematic-smiles`：canonical writerのring-closure bond markerがclosure相手側の芳香族性を見ていなかった問題（issue #395）——個別には芳香族な2原子をつなぐ非芳香族の"fusion"結合（例：`c1-2`）が、再parse時にsilentlyに芳香族結合へ変わっていた。開発用corpusのbare-parse idempotencyは73/57→0/0（完全な修正）、独立したRDKit InChI oracleで全10,000件corpus一致（不一致0件）
- 上記2件のstandardize経由修正を組み合わせた結果、そのcorpusの残差は0/4まで改善——残る4件は新規に発見・ファイルした2つの未修正issue（#407、#402系統）に起因し、本リリースには含めていない
- 詳細は`CHANGELOG.md`の`[0.20.1]`section参照

**v0.20.0**（2026-08-25）: **立体化学を保った3D構造生成、connectivity-ordered座標エンジン、`remove_hydrogens`の同一性正確性修正**
- `chematic-3d`/`chematic-py`/`chematic-wasm`：`PipelineV2Config::stereo_safe(force_field_policy)`——環に組み込まれた宣言済み立体中心（テストステロン、コレステロール等）に対する実際のギャップを一括解決する設定。従来`repair_tetrahedral_center`は暗黙のHを反映すべき座標を持たなかった。29分子×5 seedコーパスで測定：144/145（99.3%）がcorrect_and_ok、silently wrongは0件、テストステロン・コレステロールとも全宣言立体中心・全seedで5/5成功
- `chematic-3d`：`generate_coords_connectivity_ordered`——新しい公開の代替3D配置エンジン（issue #256/#255）。「まず全ての環、その後で鎖」ではなく、真の連結順序で環・鎖原子を配置。測定：raw geometryのsoundnessが差分コーパスで10/33→33/33（退行ゼロ）、post-UFFの結合長違反率は0.0000まで到達（旧エンジン自身の基準値より良好）。`generate_coords`自体は完全に無変更——デフォルト動作の切り替えではなく、選択可能な代替として提供。既存の呼び出し元はいずれも新エンジンへルーティングされていない
- `chematic-3d`：`rescue_with_distance_geometry_v2`（UFF破局的破綻からの救済ブリッジ）が再試行時に宣言済みキラリティを強制するようになった（issue #210、部分修正）。58分子コーパスで退行ゼロ、既知の残存5分子のうち1件が新たに成功。このブリッジ単体では残り4分子は未解決（上記の`stereo_safe`では解決済み）
- `chematic-chem`：`remove_hydrogens`が同位体標識水素（`[2H]`、`[3H]`）を破壊したり、呼び出しのたびに宣言済み立体・E/Z情報を静かに失ったりしなくなった——いずれも下流利用者の947万化合物規模の実世界コーパススキャンで発見された実際の出荷済みバグ。影響を受けない分子は変化しないことを確認済み、元の調査で見つかった同一性不一致290件中289件を解消
- `chematic-py`：`Mol.conformer_ensemble_v2(config)`が`embed_ensemble_v2`（決定的な複数コンフォーマー生成、力場ごとにスコープされたエネルギーランキング、全試行のprovenance）を公開——既存の`conformer_ensemble()`を置き換えるのではなく併設。新しいbest-of-10ベンチマークにより大規模でも堅牢に動作することを確認（約250/265分子、RDKit比で中央値RMSD 2.147Å／TFD 0.344）——RDKitのコンフォーマー*選択*との一致は別の未確立の主張であることに注意
- 既知の制約：`generate_coords`はまだ新エンジンへルーティングされていない、issue #210の残り4分子はこのブリッジ単体では未解決、issue #390（無関係の単一分子E/Z正確性残課題）は未解決、本リリース専用の新規全コーパス再測定は実施していない
- 詳細は`CHANGELOG.md`の`[0.20.0]`section参照

**v0.18.0**（2026-08-20）: **v0.17.0の7新形式へのPython/WASMバインディング、MMFF94 atom-typing修正、言語間一貫性の仕上げ**
- `chematic-ff`：issue #337のアリールイソチオシアネート累積二重結合CSP炭素の誤タイプを修正（`getTotalDegree() == 2`、RDKitの実際のルールに合わせた厳密な上位集合）。残る6/8分子はRDKit自身のKekulization/MMFF芳香族性認識に起因する真正のアーティファクトと再診断し（negative-control fragmentsで直接検証）、誠実な残課題として開示
- `chematic-py`：v0.17.0の7形式（mmCIF・PQR・ORCA・QCSchema・Gaussian Cube・OpenDX・LAMMPS data/dump）全てにPythonバインディングを追加——これまでRust限定だった機能。`VolumetricGrid`/`LammpsDumpFrame`のpyclassはNumPy配列プロパティを持ち、`to_opendx`/`to_opendx_lossy`のfail-closed分離も忠実に維持。`py.typed`マーカーが実際にビルド済みwheelへ含まれることを検証済み（ソースツリーではなく新規venvへのwheel installに対して`mypy --strict`が通過）
- `chematic-wasm`：同じ7形式へのWASM（wasm-bindgen）バインディング、加えて既存のJSON文字列APIに追加する形で`js_sys::Float64Array`/`Uint32Array`を返す5つの新規関数（大きな数値グリッド/行データをJSON往復なしで取得）——このクレート初のtyped-array対応
- 言語間の一貫性：同一の4つの小規模fixture（Cube・OpenDX・mmCIF・LAMMPS triclinic dump）がRust・Python・WASMの3言語entry pointで同一結果を返すことを独立検証、不一致なし
- 詳細は`CHANGELOG.md`の`[0.18.0]`section参照

**v0.17.0**（2026-08-17）: **フォーマット/Python/材料科学の相互運用性拡充、およびMMFF94電荷・結合次数の精度修正2件**
- `chematic-mol`：square-planar（`@SP1`/`@SP2`/`@SP3`相当）立体化学のMOL/SDF read/write（3D座標からの再認識）、PDBx/mmCIF・PQR・QCSchema JSON・ORCA入出力、CIF明示的symmetry操作の展開（Rust/Python両方）、共有`VolumetricGrid`型 + Gaussian Cube/OpenDX I/O、LAMMPS data/dump（trajectory）I/O（`Molecule`とは独立した専用document型）
- `chematic-py`：`chematic-crystal`の`Lattice`/`PeriodicStructure`/`Site`向け新規Pythonバインディング——開発中に既存の実バグも発見・修正（`to_cif()`が未展開symmetryのCIFを誤ってP1と再宣言していた問題）
- `chematic-ff`：MMFF94結合次数分類の修正（`torsions_missing` 257→0）とMMFF94 BCI部分電荷の修正（RDKitの原子タイプ由来形式電荷が未計算だった根本原因を特定）。本番`pipeline_v2_mmff94_strict`：240/265→241/265
- リリース整備：宣言MSRVの実態不一致を修正（1.88へ引き上げ、専用CIジョブで継続検証）
- 詳細は`CHANGELOG.md`の`[0.17.0]`section参照

**v0.16.0**（2026-08-15）: **周期構造の相互運用性（CIF/POSCAR/FPS）と一般化立体化学基盤**
- `chematic-mol`：新設のoptional `crystal` featureが既存のCIF reader/writerを`chematic_crystal::PeriodicStructure`へ橋渡し（`parse_cif_periodic_structure`/`write_cif_periodic_structure`）——セルパラメータを`Lattice`へ、`_atom_site_occupancy`を`Occupancy`へ、disorderを共有する複数の`_atom_site_*`行を1つの`PeriodicSite`の複数species listへ統合。新設の`CifSymmetryStatus`列挙型により、真にP1な CIFと、symmetryを宣言しているがこのparserが未展開のCIFを区別（後者を暗黙にP1扱いしない）。`chematic-crystal`自体は`chematic-mol`/`Molecule`から独立したまま（依存方向は`chematic-mol`→`chematic-crystal`のoptional依存のみ）
- `chematic-crystal`：POSCAR/CONTCAR（VASP構造ファイル形式）のnative read/write——`parse_poscar`/`parse_contcar`/`write_poscar`、VASP 5のみ対応、2種類のscale-factor記法、Direct/Cartesian座標、selective dynamics、ion velocities、CONTCARのpredictor-corrector（MD再開用セクション）はそのまま保存（VASP公式ドキュメントも数値レイアウトを規定していないため）
- `chematic-fp`：新規`fps`モジュール——chemfp/OpenBabelで普及しているテキストベースのフィンガープリント交換形式「FPS（Fingerprint file format）」のstreaming read/write。16進ビット順序はchemfp仕様と照合済み、bit-vector表現は既存の`BitVec2048`/`BitVecN`をそのまま流用
- `chematic-core`：新規`stereo_geometry`モジュール——立体配置を、配位幾何（`Tetrahedral`/`SquarePlanar`、将来のTBP/octahedral拡張に備え`#[non_exhaustive]`）と、その幾何の回転対称群下でのリガンド順列の同値類として表現。tetrahedralはA4（位数12）、square-planarは素朴な位数4の面内回転のみの群ではなく、trans-pair分割のS4安定化群（位数8）。`chematic-smiles`内の2つの独立した手書きremappingアルゴリズムを置き換え、`@`/`@@`/`@SP1`/`@SP2`/`@SP3`の既存意味は完全維持（88件のfixtureでbyte-identicalなcanonical SMILES回帰確認済み）。開発中に見つかった実害バグも修正——square-planar中心が`chematic-3d`で誤ってtetrahedral判定に流れ込み、浮動小数点ノイズでSatisfied/Violatedが決まっていた問題。開発過程で一時的に発生したallene端立体中心のparity回帰も発見・修正し、golden-value testで固定
- `pipeline_v2` vs RDKit 2026.03.4ベンチマークのrelease-grade再測定（2026-08-06時点の古い数値を置き換え）：`mmff94_strict`が149/265→239/265。新たな発見として、torsion parameter不足が現在のMMFF94の主要な残存ギャップであること（`complete_bonded_term_gated`失敗の71%がtorsion欠如起因、OOP・bond起因は0%）——次のMMFF94ロードマップ項目への直接的な根拠
- 詳細は`CHANGELOG.md`の`[0.16.0]`section参照

**v0.15.0**（2026-08-14）: **`chematic-crystal`——周期（結晶）構造の基盤crate新設、MMFF94 Bond/Angle empirical rule対応（issue #227）**
- 新規crate `chematic-crystal`：周期（結晶）構造の表現とジオメトリ計算——`Lattice`（三斜晶系対応、行列/逆行列/逆格子ベクトルをvalidation済み）、`FractionalCoord`/`CartesianCoord`、`PeriodicSite`/`SiteSpecies`/`Occupancy`（disorder対応可能な複数species設計）、`PeriodicStructure`は近似（`round()`）ではなく厳密な周期最小像距離計算（等距離の周期像が複数ある場合は辞書順最小のimageへ決定論的に解決）、cutoff近傍探索、diagonal supercellを提供。`chematic_core::Molecule`（結合グラフ）の拡張ではなく意図的に独立した型。optionalな`serde` feature、facadeの`crystal` featureは`full`に含まれる（`default`は空のまま変更なし）。symmetry・CIF parser変更・Python/WASM/MCPバインディングは今回未対応
- `chematic-ff`：HalgrenのMMFF.V eq. 18-20 empirical Bond-stretch/Angle-bend ruleを移植（`mmff94_bond_energy_resolved`/`mmff94_angle_energy_resolved`——既存の`mmff94_bond_energy`/`mmff94_angle_energy`のシグネチャは変更なし）。既存の完全一致テーブル/`eqLevel`ラダー探索の後にのみ試行するため、実データテーブルのヒットを上書きすることはない。この過程で、RDKit実データのAngleテーブルにある97行（中心原子タイプのみで決まる汎用`theta0`デフォルト値）がchematicの移植版に欠落していたことを発見・復元。1タプルのみ意図的に未解決（fail-closed）のまま——外側原子タイプに等価クラステーブルのエントリがなく、RDKit実装のC++コード自体が未定義動作（nullチェックなし参照）を起こす箇所であり、live oracleの返り値を明確な解決メカニズムに帰属させられなかったため。既存の5件のMMFF94原子タイプ判定ギャップも修正し、RDKitの`eqLevel`原子タイプ等価ラダーをAngleテーブル探索へ移植（いずれも今回のempirical rule対応の前提作業）。265分子Wave 1コーパスでの実測（本番の最小化パス経由、いずれもper-molecule完全突合で検証済み・regressionゼロ）：v0.14.1からのリリース全体での変化は178/265ではなく**158/265→248/265（失敗107→17）**；empirical rule対応そのものの寄与分（同一リリース内で先にmergeされた原子タイプ判定/`eqLevel`修正を除いた差分）は**178/265→248/265（失敗87→17）**。最終状態でも`MinimizationFailed`のまま残る3分子は、いずれもv0.14.1時点で既にnon-`Ok`（MissingParameters）だった分子——原子タイプ判定修正で初めて実パラメータが揃ったことで、それまで隠れていた既存の幾何問題が可視化されたものであり、既存の成功例からの後退ではない
- 詳細は`CHANGELOG.md`の`[0.15.0]`セクション参照

**v0.14.1**（2026-08-12）: **抗がん白金配位化学の互換性修正、Extended XYZ（extxyz）読み書き対応**
- `chematic-core`: `valence_inferred_hcount`が`BondOrder::Dative`結合のdonor側を通常の共有結合と同じに扱っていたため、implicit水素数計算が誤っていた——`N->[Pt]Cl`のような括弧なしdative donorが`NH3`ではなく`NH2`と計算されていた。donor側dative結合はvalence合計に0を寄与するよう修正。白金配位化学ベンチマークで発見したが、Fe/Co/Pd/Ruのacceptorでも検証済みの一般的な修正（白金固有ではない）
- `chematic-mol`: MDL bond type 9（dative/coordinate結合——RDKit実装が`Bond::BondType.DATIVE`をV3000で書く際に使う規約）が、V2000・V3000両readerで`BondOrder::Single`へ暗黙に丸められ、配位結合の意味情報が読み込み時に静かに失われていた。両readerともcode 9を`BondOrder::Dative`として解釈するよう修正、V3000のwriterもcode 9を出力するよう修正
- `chematic-chem`: `avg_mass`/`mono_mass`が軽い主族元素約24種のみをカバーし、それ以外の全元素で`atomic_number as f64`へ静かにfallbackしていた——遷移金属・ランタノイド・アクチノイド・重い後周期元素は全てエラーなしで大きく誤った質量を返していた（白金：原子番号78、実際の質量約195Daのところ「78.0 Da」を返していた）。`Element`が持つ全118元素へ拡張、RDKitの周期表データを出典として使用。既存約24元素の値はそのまま維持（セレンなど、本プロジェクトの値が現行IUPAC標準でRDKit側が2013年以前の旧値を採用しているケースは既存値を優先）
- `chematic-mol`: Extended XYZ（extxyz）形式の新規対応——`parse_extxyz`/`write_extxyz`、`ExtxyzReader`/`ExtxyzWriter`、`parse_extxyz_all`。既存のmulti-frame `XyzFrame`型の拡張として実装（ASEの`Lattice=`セル行列、型付きper-atom `Properties=`列、任意の`key=value`フレームメタデータ）。プレーンXYZファイルはextxyz readerを通しても無変更で往復する。Python: `from_extxyz`/`from_extxyz_all`/`to_extxyz`。WASM: `mol_from_extxyz`/`extxyz_frame_json`/`to_extxyz_json`。**Breaking（Rust APIのみ）**: `XyzFrame`に公開フィールド3つ追加、`XyzError`にvariant 7つ追加、`write_extxyz`の戻り値が`Result<String, XyzError>`に変更——これは既にcrates.io公開済みのv0.14.0 Rust APIに対する実際の破壊的変更であり、未リリースAPIへの変更ではない点に注意
- 白金配位化学の立体化学（square-planar cis/trans、例：cisplatinとtransplatinの区別）は依然として表現不可能——今回のリリースでは測定のみ行い、意図的に未修正（`validation/platinum/FEASIBILITY.md`参照）
- 詳細は`CHANGELOG.md`の`[0.14.1]`セクション参照

**v0.14.0**（2026-08-11）: **立体化学を考慮したdistance geometry——宣言済みE/Zをbound-matrix制約として強制、`enforce_chirality`とpost-minimization stereo verificationの組み合わせ対応、Python/WASM公開**
- `chematic-3d`: v0.13.0のissue #285 release-gate waiverの根本原因を特定・修正——`apply_vdw_bounds`の汎用non-bonded Van der Waals下限が、宣言されたE/Z立体化学に関わらず宣言済みE/Zアルケンの1-4置換基ペアに適用され、正しいcis配座が構造的にサンプリング対象から除外されていた。新規`apply_declared_ez_bounds`（`enforce_chirality`時のみ）が、汎用Van der Waals下限の適用より前に、same-side/opposite-sideの1-4距離の解析的境界をbond matrixへ交差させることで、post-hocな修復・retry・reflectionではなく構造的に正しい配座へ到達可能にした。四面体キラリティ（distance matrixでは原理的に表現不可能——分子とその鏡像はpairwise distanceが同一）と異なり、宣言済みE/Zはcis/transが異なるスカラー距離であるため、genuinely distance-representable。265分子コーパスの宣言済みE/Zサブセット（39分子）で測定：stereo-satisfied 22→42、violated 23→3、pipeline成功率・健全性は無変化
- `chematic-3d`: `embed_pipeline_v2`のconfig検証ゲートが従来`stereo_policy`が`Ignore`以外なら`enforce_chirality: true`を拒否していたが、コーパス測定によりこれが誤りと判明——`enforce_chirality`はembedding時の正しさのみを保証し、force-field minimization（宣言済みstereoの概念を持たない）が正しくembedされたE/Z結合を事後的に境界の反対側へ動かしうる（実分子2件で確認、force field無しの再実行で検証）。`enforce_chirality: true`は`StereoPolicy::VerifyOnly`とも併用可能になり、既存のpost-minimizationゲートがこの失敗モードをsilentな誤stereo「成功」ではなく型付きエラーとして検出する
- `chematic-py`、`chematic-wasm`: `enforce_chirality`（デフォルト`false`）が`PipelineV2Config`/`PipelineV2Config.safe()`（Python）および`enforceChirality` JSONフィールド（WASM）で実際に設定可能なパラメータ/フィールドに——どちらのbindingもこれまでこのフィールドを一切引き渡していなかったため、上記の修正はPython・WASM双方から到達不可能だった
- `chematic-rxn`: `suzuki_biaryl`のretro-template修正（issue #294）——`[c:1][c:2]`は実際のbiaryl結合に一切マッチせず、環内芳香族結合のみにマッチしていた（このクレートのSMILES規約では、隣接する2つの芳香族原子間に明示的な結合トークンがない場合デフォルトで芳香族結合になるため）。`[c:1]-[c:2]`に修正。副次的発見：`DEFAULT_TEMPLATES`59件中14件が一切パースされていなかった——issue #296として記録、本修正では対応せず
- opt-in限定——`enforce_chirality: false`が引き続き全箇所でデフォルト。デフォルトのconformer経路（`generate_coords_etkdg`/`Mol.conformer_ensemble()`）は無変更
- 詳細は`CHANGELOG.md`の`[0.14.0]`セクション参照

**v0.13.0**（2026-08-10）: **MMFF94 stretch-bend／torsionパラメータ選択パリティ（両方breaking）、per-atomステレオセンターAPI、E/Z完全性判定、macrocycle検出、notation非依存atropisomer検出・割り当て、XYZ入出力**
- `chematic-ff`: `mmff94_stbn`/`mmff94_stbn_type_only`が`MMFF94_STBN`テーブルのlookup keyとして、これまで代用していた粗いangle type（0-8）ではなく、RDKit実装の本来の細かいstretch-bend type（`getMMFFStretchBendType`、0-11）を使うように（issue #227）— 265分子コーパスで427件のstretch-bend routing候補中220件が、RDKitの汎用Dfsb周期表デフォルト値から正しい専用パラメータへ移行。`angle_type_for`のring-offset式もRDKit実装の`getMMFFAngleType`に合わせて修正。**Breaking**: `mmff94_stbn`/`mmff94_stbn_type_only`の先頭`u8`引数は`stretch_bend_type`（`angle_type`ではない）— 新規`pub stretch_bend_type_for`で計算
- `chematic-ff`: `torsion_type_for`が、atom-type membershipのみに頼る旧分類ではなく、j-k結合の実際のMMFF bond type（`bond_type_for`を再利用）とRDKit実装のlocal bond-adjacencyベースring 4/5-override判定を使うように — 1,107件の欠落torsion候補のうち76.9%が分類修正のみで解決、コーパス全体13,530件のtorsion中1,792件の「値はあるが誤ったパラメータ」だったsilent wrong-parameter populationも是正（99.1%をRDKit oracleで直接検証、新規消失0件）。**Breaking**: `torsion_type_for`のシグネチャが`(rings, i, j, k, l, tj, tk)`から`(mol, i, j, k, l, ti, tj, tk, tl)`へ変更
- `chematic-mol`: XYZ／multi-frame XYZ読み書き（`parse_xyz`/`write_xyz`、`XyzReader`/`XyzWriter`）— 明示的水素は実atomとして保持、結合次数推定なし、原子数不一致・非有限座標はfail closed
- `chematic-perception`: `stereo_centers(&Molecule) -> Vec<(AtomIdx, bool)>`がper-atomのtetrahedral stereocenter分類を公開（issue #263、従来は集計カウントのみ）— 追加時に2件のバグを発見・修正: 負電荷atomでのMorgan-rankヘルパーの`u64` overflow（issue #267）、implicit hydrogenのrank-0センチネルが実atomの正規化rank 0と衝突し正しく指定されたstereocenterを見落としていた問題
- `chematic-chem`: `ez_completeness(&Molecule) -> EzCompleteness`（issue #264）が宣言済みE/Z二重結合の specified/unspecified/total を、RDKit実装のstereo-bond適格性ルール（末端・対称結合を除外、8員未満のring結合をBFS最短閉路で除外——SSSRのみでは検出できないbridged-bicyclic系（norbornene等）も正しく扱う）に沿って報告
- `chematic-chem`: `detect_atropisomers`/`assign_atropisomer_chirality`が完全にnotation非依存に（issue #262、#276）— 検出はSSSRベース（別々の環に属する2つの芳香族炭素、両環ともortho置換あり）で、SMILESが環間結合を明示的に書くか暗黙にするかに依存しない。chirality割り当て側の冗長なbond-order判定も、検出側の分類と一致するよう修正
- `chematic-perception`: `is_macrocycle(ring: &[AtomIdx]) -> bool`（issue #266）— `chematic-3d`側に重複していたハードコードされた閾値を単一の共有述語に統一
- **v0.13.0リリースゲート注記**: 265分子Wave 1コーパス中2分子（`chembl_tier_b_0126`/`0168`）で立体中心充足数の1件退行を確認、根本原因はv0.12.0から存在するdistance-geometry埋め込み段階の既存バグ（今回修正したtorsion分類バグに偶然マスクされていた）と特定 — RDKit自身のMMFF94も同一の出発座標を与えると同じ挙動を示すことを確認済み。明示的waiverの下でリリース（issue #285）— 今回のMMFF94修正が新たに生んだ不具合ではない
- 詳細は`CHANGELOG.md`の`[0.13.0]`セクション参照

**v0.12.0**（2026-08-09）: **MMFF94 stretch-bend本番修正（breaking）、fused/multi-ring分子の3D初期構造修正**
- `chematic-ff`: `mmff94_stbn`がRDKit実装の29行周期表行stretch-bendデフォルトへfallbackするように — 265分子コーパスでstretch-bend欠落2,107→0。**Breaking**: `mmff94_stbn`に`atomic_num_{i,j,k}: u8`が必須引数として追加（旧動作は`mmff94_stbn_type_only`として維持）、Python生`PipelineV2Config(...)`コンストラクタに`gate_mmff94_stretch_bend`が新規必須引数として追加（`.safe(...)`は無変更）
- `chematic-3d`: `dg::generate_coords`がfused/multi-ring分子で原子座標衝突・結合stretchを起こさなくなった（issue #185/#252）— UFF minimizer自体のバグではなく、生成された初期構造自体の問題と判明。265分子コーパスで`MinimizationFailed` 28件全てが解消、regression 0件。fused-ring seam方向（issue #255）とchain-bridged ring island配置（issue #256）は既知の未修正課題として別issue化
- 詳細は`CHANGELOG.md`の`[0.12.0]`セクション参照

**v0.11.0**（2026-08-04）: **MMFF94 O2CM typing coverage、SMIRKS/CDXMLステレオ正当性、2D/3Dレイアウト修正**
- `chematic-ff`: O2CM末端酸素typing gapを解消（issue #227 Priority 1A-3）— 265分子Wave 1コーパスでatom-type一致率98.82%→99.37%、strict-gate最小化成功123→130/265
- `chematic-rxn`: SMIRKS product chirality割り当てをparity-aware化 — reorderされたmapped neighbor順序に応じて`@`/`@@`flagを正しく反転/検証
- `chematic-mol`: CDXMLリーダーがdirectional wedgeからtetrahedral stereoを認識するように（RDKit issue #9359）
- `chematic-depict`: 独立した(非fused)環系が2D上で同一座標に衝突しなくなった
- `chematic-3d`: ETKDG macrocyclic amideの1-4距離boundを真のcis/trans役割で分割
- 詳細は`CHANGELOG.md`の`[0.11.0]`セクション参照

**v0.10.1**（2026-08-02）: **MMFF94数値atom typing修正（正当性ホットフィックス）**
- `chematic-ff`: MMFF94が誤ったelementのparameter行へ衝突し、その結果の物理的に誤ったenergyを成功として返し得るバグ（issue #227「furan collision」）を修正 — 芳香族atom typerがRDKit実装の5員環/6員環alpha/beta-heteroatom分類を実装していなかったのが根本原因。pinしたRDKitソースから移植し、provenance付きnumeric-typeレジストリと、このバグの再発を構造的に不可能にするconstruction-time semantic-compatibility invariant（不整合はfail closed = `NumericTypeError`）を追加。このinvariantが同種のバグをさらに2件検出: protonated amine Nとanionic Oが互いのelementのparameter行として誤typeされていた。265分子コーパス（本番API）で測定: MMFF94最小化成功 44 → 102、pin済みRDKitオラクル比較で6693原子中cross-element型不一致 0件（91.83%完全一致）
- coverage完成リリースではなく正当性ホットフィックス — issue #227は未完了のまま（`MissingParameters` 140件、`MinimizationFailed` 22件が残存、stretch-bendは未gating、full-corpus energy/gradient parityも未実施）。MMFF94結果をキャッシュしている場合は`CHANGELOG.md`の`[0.10.1]` Migration notesを参照
- 詳細は`CHANGELOG.md`の`[0.10.1]`セクション参照

**v0.10.0**（2026-08-01）: **match-level SMIRKS reaction適用、MRVの2Dステレオ認識、shared E/Z carrier bond修正**
- `chematic-rxn`: `find_reaction_matches`/`apply_reaction_match`（issue #225）— SMIRKSのマッチ列挙と個々の適用を分離する公開API
- `chematic-mol`: MRVリーダーが2Dウェッジ/ハッシュtetrahedralとE/Zステレオを自動認識するように（issue #202）
- `chematic-smiles`: shared E/Z carrier bondをjoint component solverで解決（issue #149）— 18件中10件が完全にpermutation-invariantに
- 詳細は`CHANGELOG.md`の`[0.10.0]`セクション参照

**v0.9.0**（2026-08-01）: **Python/WASM向けopt-in 3D embedding pipeline v2、WASM対応monotonic clock**
- `chematic-py`/`chematic-wasm`: `embed_pipeline_v2` — torsion知識付き距離幾何 + stereo検証/repair + policy-gated force fieldを1本化、全stage分のevidenceを返す。既定の3D APIは無変更
- `chematic-3d`/`chematic-smarts`: `wasm32-unknown-unknown`実環境で`Instant::now()`が無条件panicするバグを修正（issue #219, #221）
- 詳細は`CHANGELOG.md`の`[0.9.0]`セクション参照

それ以前の全バージョン履歴（v0.1〜v0.8.1を含む、各リリースのコーパス単位before/after数値・根本原因・migration notesまで）は
[CHANGELOG.md](CHANGELOG.md)を参照。

<details>
<summary>v0.7.1以前の開発履歴（旧形式、参考用）</summary>

**v0.7.1**（2026-07-27）: **`run_reactants`/`canonical_smiles()` パフォーマンス修正**
- `chematic-smiles`: `canonical_smiles()` が毎回 `CanonicalWriter::write_all()` を無駄にもう一度呼んでいたバグを修正（individualize-refineのタイブレーク解決時に勝者は既に一度書き込まれていたが、それを捨てて再度書き直していた）。`be5dbb1`（0.4.26）の正しさ修正で生じた本物の性能回帰で、対称性の高い分子（単純な環、ケージ構造、`CF3`/`tBu`系置換基）で最も顕著（外部利用者の`run_reactants`/`apply_retro`回帰報告経由で発覚、chematic 0.4.30 vs 0.4.25で単体`canonical_smiles()`が45-48倍遅い）。純粋なリファクタリングで出力はバイト単位で同一（`be5dbb1`自身のgolden-stringテスト、issue #50のE/Z回帰スイートを含む既存テスト全てで検証済み）。あわせて`chematic-rxn`に`perf-instrumentation`機能と`reaction_transform_perf_report`ベンチマーク例を追加。既知の残課題として真の自己同型オービット枝刈りは未着手
- 詳細は`CHANGELOG.md`の`[Unreleased]`セクション参照

**v0.7.0**（2026-07-26）: **MOL/SDFから2Dウェッジ/ハッシュ＋E/Z立体化学を自動認識、検証付きcanonical SMILES重複排除、CIP Rule-5リン修正、native InChI明示的水素/同位体修正**
- `chematic-mol`/`chematic-perception`: MOL V2000/V3000/SDFリーダーが読み込み時にtetrahedralウェッジ/ハッシュparity（PR #154）とE/Z二重結合方向（PR #162）を自動認識するようになった（CIP非依存）。型付きopt-in診断（`StereoDiagnostic`/`EzDirectionDiagnostic`）は不正・曖昧な入力に対して推測しない。広域corpus検証（4,999分子、RDKit 2026.03.3比較）：E/Z — RDKit解決済み622二重結合、semantic inversion 0件、false positive 0件。構築中に見つかったV2000 MDLコード4バグ・V3000 `CFG`バグ2件・V2000 writerバグも修正。PR #162自身のcorpus突合により**新規・未修正のギャップ**を発見: `canonical_smiles()`が一部の分子で正しく`write()`されたE/Zマーカーを失うケースがある（aromaticity非対応のcarrierグルーピングが原因）— まだissue化していない
- `chematic-inchi`: 検証付きcanonical SMILES重複排除（`dedup`モジュール）を新規実装 — 高速なcanonical SMILES候補バケット化をnative InChIによる検証済み同一性と突合。指定済み立体中心のlegacy CIPランク付けが未解決の場合は誤マージのリスクを避けてfail-closed（`VerificationUnavailable`）にする(実際の5,000分子corpus検証で見つかった本物のfalse `VerifiedDuplicate`を解消)。follow-upとして[#161](https://github.com/kent-tokyo/chematic/issues/161)（accurate CIP preflightで保守的なケースの大半を回復できる可能性）を追跡
- `chematic-cip`: CIP Rule-5擬似不斉r/sのリン修正
- `chematic-inchi`: native InChI（`native-inchi` feature）の明示的水素/同位体変換修正
- 詳細は `CHANGELOG.md` の `[Unreleased]` セクションを参照

**v0.6.0**（2026-07-25）: **RDKit bit-exact ECFP4のクロス言語stable API化、canonical SMILESのE/Zマーカー一貫性、opt-in芳香族flag authoritative降格**
- `chematic-fp`/Python/WASM: RDKit bit-exact Morgan/ECFP4パスをクロス言語opt-in APIとして公開（Python `Mol.rdkit_ecfp4()`、WASM `rdkit_ecfp4_bitvec()`）。radius×fpSizeの20セル全てを個別にlive RDKit oracleで再検証（closed enumで未対応値は構築不可）。Rust/Python/WASMを実際にビルド・実行してbyte-identicalを確認済み。`ecfp4()`の既定動作は無変更
- `chematic-smiles`: 入力のatom順序によってE/Z方向マーカーが異なる側鎖bondへ載る問題を修正 — permutation invariance **93.0% → 98.1%**（264/282件が収束、残り18件はshared-carrier-bondの破損リスクを避けるため意図的に未解決、issue [#149](https://github.com/kent-tokyo/chematic/issues/149)で追跡）。また、bridged-bicyclic canonicalizationの既存バグ疑惑を調査し**誤診断と判明**（RDKit InChIで検証した結果、問題とされた2つのSMILESは実際には異なる分子だった）
- `chematic-perception`: 新規opt-in `apply_aromaticity_authoritative_experimental` — 芳香族flagの昇格/降格をHückelモデルの計算結果に対して双方向で忠実にする（既定の`apply_aromaticity`/`apply_aromaticity_ex`は無変更、byte-identical確認済み）。fused diazine（quinazoline/quinoxaline/purine型）のring-fusion bond誤分類も修正、既存の32件のfalse-positive regression pinも副次的に解消
- 詳細な項目別数値と既知の制限事項は `CHANGELOG.md` を参照

**v0.5.0**（2026-07-23）: **CIP不要の2Dウェッジ由来local parity、charge対応kekulization（従来失敗していた6分子クラス）、PAINS/Brenkのbudget付きmatching**
- `chematic-perception`: `local_parity_from_wedges`/`apply_local_parity_from_wedges` を新規追加 — CIPランキングを一切使わず、wedge/hash結合と2D座標から直接 `Atom.chirality`/`stereo_neighbor_order` を計算。符号規約はRDKitの生のchiral tagに対して実測で決定（類推ではない）。まだどのreaderのデフォルトparseからも呼ばれない
- `chematic-core`: `kekulize()` の原子マッチング規則がcharge非対応でTelluriumも欠落していた問題を修正 — tropylium、imidazolium、pyridinium、pyrylium、tellurophene、phospholeがRDKitとbond単位で完全一致するKekulé構造でkekulize成功するようになった。`Element::normal_valences()` にTelluriumの実証済みエントリを追加し、ECFP4の芳香族性不一致も解消
- `chematic-perception`: charge対応Hückel π電子計算 — tropylium、imidazolium、pyridinium、pyryliumがRDKitの芳香族atom/bond flagと完全一致（tellurophene/phosphole、および芳香族flagのauthoritative降格修正は別途対応中）
- `chematic-smiles`: writerの2件のバグを修正 — bracket強制原子（isotope/charge/atom-map）が暗黙水素を出力し忘れていた問題（`[NH4+]` → `[N+]`）、および隣接する二重結合を持たないwedge結合が意味のないSMILES `/`・`\` トークンとして出力されていた問題
- `chematic-smarts`/`chematic-chem`: 対称性の高いターゲットでPAINS/Brenk部分構造探索が数分間ハングする可能性があった問題を修正 — VF2に明示的な探索budgetを導入し、`Found`/`NotFound`/`BudgetExhausted` の三値を返すことで、探索打ち切りを黙ってfalse negativeへ畳み込まない設計に
- 詳細な項目別数値と既知の制限事項は `CHANGELOG.md` を参照

**v0.4.29**（2026-07-10）: **Kabsch回転バグ修正 + SDF V3000/CDXML書き込み、Avalon FP、O3A**
- `chematic-3d`: `align_coords` のKabsch回転が逆方向に計算されるバグを修正（純粋な並進以外のアライメントでRMSDが大きく水増しされていた。v0.4.28からcrates.io/PyPI/npmで公開されていた）；O3A原子対応のための `correspondence_search`
- `chematic-mol`: SDF V3000書き込み配線；CDXML書き込み
- `chematic-fp`: Avalonフィンガープリント

**v0.4.28**（2026-07-09）: **SMARTS性能改善、レジストリ再同期**
- `chematic-smarts`: 存在チェックの早期終了 — `bulk.substructure_search` がRDKitより2.2倍高速に
- v0.4.23〜v0.4.27はgitタグが未pushでcrates.ioのみ最新（PyPI/npm/GitHub Releasesが遅れていた）だったため、この版で3レジストリを再同期

**v0.4.27**（2026-07-04）: **記述子修正、RWMol/FCFP、veridict CIゲート**
- `chematic-chem`: `kappa1-3`、`balaban_j`、`labute_asa`、`bcut2d`、`hall_kier_alpha` 記述子修正
- `chematic-fp`: `useFeatures=True` FCFP
- `chematic-mol`: RWMol インプレース編集
- CI: veridictベースの性能/Criterion/精度ドリフト回帰ゲート；統合テストのCIカバレッジギャップ修正

**v0.4.26**（2026-06-29）: **反応でのE/Zステレオ転写 + 検証Sprint 6/7**
- `chematic-rxn`: `run_reactants()` で反応物の `/`/`\` 二重結合幾何が生成物に保持されるように（従来は変換で失われていた）
- 検証: RDKitに対するカノニカルSMILES差分検証（Sprint 6）；SMARTS/芳香族性差分テスト + I/O互換性（rdkit_compat Sprint 7）；残存するRDKit差異の根本原因をMorganランクではなく芳香族性ラウンドトリップと特定

**v0.4.25**（2026-06-29）: **`chematic.rdkit_compat` レイヤー**
- `chematic-py`: RDKit API互換レイヤー（Sprint 1〜5）— Morgan `bitInfo`、Fingerprint/Mol/Atom/Bond/RingInfo互換性、RDKitとの差分テスト；ストリーミング `SDMolSupplier`/`SDWriter`/`Mol.GetProp`
- `chematic-perception`: `AromaticityAlgorithm::RdkitLike` — Se/TeカルコゲンをRDKitのモデルに合わせて処理

**v0.4.24**（2026-06-29）: **CIP Rule 5、架橋頭部/回転可能結合/TPSA/MRを100%に、HDFフィンガープリント**
- `chematic-chem`: CIP Rule 5立体タイブレーク（ステレオセンター一致率 99.8% → 99.98%）；架橋頭部検出 98.5% → 100%；回転可能結合 99.1% → 100%；TPSA 100%；モル屈折率 97.5% → 100%（3環XOR拡張）— いずれも5,000分子ChEMBLコーパス
- `chematic-py`: `bulk.descriptors_array()` 列指向numpy出力；真のストリーミングSDF（`SdfFileReader`/`iter_sdf_batched`）；`screen()` 化合物フィルタワークフロー
- LLM/RAG: 表現ルーター（`to_llm_text`, `best_representation`）、分子コンテキストパック、**Hyper-Dimensional Fingerprints（HDF）** — 学習不要の密な分子ベクトル

**v0.4.23**（2026-06-26）: **LogP 96.5% → 99.7%**
- `chematic-chem`: `crippen_anchor_sets` を `uniquify: false` に修正し、対称な三重結合（内部アルキン）がVF2マッチで両方向とも得られるように（従来は片方が汎用 `[#6]` 値にフォールバックしていた）

**v0.4.22**（2026-06-26）: **CITATION.cff + `chematic.doctor()`**
- `chematic-py`: `doctor()` 自己診断機能；README に信頼性マトリクスを追加

**v0.4.21**（2026-06-25）: **LLM/Jupyter向けHTML/Markdownレポート**
- `chematic-py`: `chematic.report()` 自己完結型HTML化合物グリッド、`chematic.compare()`、`mol.review()` Markdown解析
- ドキュメント: `benchmarks/`/`validation/` 再現可能な精度履歴

**v0.4.20**（2026-06-25）: **ETKDGトーションKB 44 → 80パターン、`mol.describe()`/`diff()`**
- `chematic-3d`: 6員環/5員環脂肪族環の椅子形/封筒形コンフォメーション；SMARTSベースのトーションルールを高精度事前チェック層として追加
- `chematic-py`: LLM/MCPエージェント向け `mol.describe()`/`mol.diff(other)`；`bulk.generate_3d`/`tanimoto_matrix`/`standardize`

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

それ以前の v0.4.x の開発（テンプレート逆合成、AutoDock/UFF、ケクレ化 blossom
アルゴリズム、PyO3 バインディング、native-inchi）と v0.1〜v0.3 の全履歴:
[CHANGELOG.md](CHANGELOG.md)

</details>

---

## 既知の制限事項

- **`canonical_smiles()` はE/Zステレオ化学に対して部分的に正規化されました — それでも重複排除やキャッシュキーとしてはまだ安全ではありません**: 孤立した単純なE/Z二重結合には `/N=N/` と `\N=N\` のような2通りの等価な正しい記法があり、ライターは従来どちらを出力するか正規化していませんでした。一般ケースを修正済み: ある二重結合と、それに幾何学的に連動する全ての方向性結合(共役鎖全体を含む)からなる「連結E/Z系」ごとに、カノニカル書き込み順で最初に現れる方向性結合が常に `/` になるよう正規化するようにしました(入力の綴りに依存しません)。5,000分子ChEMBLコーパス・worst-of-10で測定: E/Z限定の自己不安定性(四面体ステレオを除去)は**9.76%→5.50%**に改善(275/5000が依然不安定)。この変更による構造的正しさへの影響はなし(ChEMBL **0/5000**、非環式ポリエンコーパス **0/33** を再検証済み)。残る275件は**100%が見た目上の問題のみ**(RDKit的にはどの変異体も同一分子で、破損はゼロ)と確認済みですが、原因は**根本解明しきれていない混在プール**です: 約半数は特定のモチーフ(2つ以上の環外二重結合を持つ小員環、例えば交差共役した環状ジイミン)に一致し、そこでは「1つの系」とみなすべき物理結合の集合が入力の綴りに対して不変になっていません。残り半数は未特定です。これが完全に解消するまでは、ステレオを持つ分子の約18分子に1つ(旧: 約10分子に1つ)が、同一分子に対して2通りの、それぞれ単独では正しい `canonical_smiles()` 文字列を生成しうるため、現時点では重複排除やキャッシュキーとして使わないでください。この点が重要な用途では、当面 `apply_aromaticity()` で正規化した文字列を独自の重複排除キーとして使ってください。
- **カノニカルSMILESの構造的破損 — 修正済み**: 修正前は `canonical_smiles(parse(x))` が `x` の入力走査順によっては(単に綴りが異なるだけの等価な文字列ではなく)別の立体異性体を静かに出力することがありました。5,000分子ChEMBLサブセット・worst-of-10走査・RDKit検証済みの構造正しさで測定: **4.28%（214/5000）** の分子で、少なくとも1つの変異体が誤った分子として読み戻されていました。根本原因は独立した2件のパーサーバグで(当初疑われていた「共役二重結合のマーカーは結合をまたいで幾何学的に連動している」という診断は誤りであったことが判明— 下記参照)、いずれも実在する分子と最小構成の回帰テストで確認済みです: (1) 環closure(`/`/`\`) の方向マーカーを閉環側の出現位置で読み取る際、開環→閉環方向へのフリップをせず生の値のまま保存しており、共役E/Z鎖の連結結合がたまたま環closureを経由する場合に破損していた。(2) 自身の分岐の中で環closureの相手が閉じるステレオ中心は、再利用可能な環番号ではなく出現ごとに一意なIDで隣接原子順序の解決を行うべきところ、番号ベースで解決していたため、SMILES文字列の後方で同じ番号が無関係な環に再利用された際に、ステレオ中心の隣接原子順序が静かに乗っ取られ破損していました。**両修正後: ChEMBLコーパスで構造正しさは100%（0/5000）** — これは、rankingという第3の無関係な修正の有無を含む独立した3通りの再構成手順で3重に確認済みです。両方の根本原因が環closure固有のものである一方、レチノイド・カロテノイド・プロスタグランジン・ロイコトリエン・ポリエン系マクロライドは長い共役系を**非環式**の鎖として持っており、ChEMBLのランダムサンプリングにはほぼ含まれないため、これら5クラスの実在化合物33種からなる専用コーパス(トレチノイン、β-カロテン、リコペン、アンフォテリシンB、ロイコトリエンB4など、`scripts/polyene_corpus.csv`)でも独立に再検証しました: **worst-of-30で0/33（0.00%）**。同一コーパスの未修正コードでは12/33（36.36%）の破損を陽性対照として確認済みです(12件全てが環closureを多く含む構造で、完全非環式のリコペンを含め純粋な非環式の例は未修正コードでも一度も破損しませんでした)。これにより当初の「共役系全般が原因」という診断が誤りであったことが直接確定し、残存する破損クラスは見つかっていません。骨格限定・四面体限定の自己**安定性**もいずれも0%に到達(旧: 0.16%、4.36%)。全ステレオを含めた生の自己安定性は**86.02%→90.28%**（不安定率13.98%→9.72%）に改善 — 残差は全て上記の非破損な方向正規化ギャップによるもので、破損の残存ではありません。往復不変性（`canonical(parse(canonical(m))) == canonical(m)`）はわずかに改善（**98.26%→98.32%**）— この指標はそもそも破損クラスを直接測定していなかったためです。
- **環知覚（SSSR）は非決定的・非最小だったが修正済み**: 旧 `find_sssr` は単一の全域木から非木辺ごとに基本閉路を1つずつ生成していたため、木の形によっては最小でない環（例: ナフタレン `c1ccc2ccccc2c1` が決定的に `[6,10]` を返す。正しくは `[6,6]`）を返していました。現在の `find_sssr` は Horton アルゴリズム（全頂点×全辺の最短路木から候補閉路を生成、O(V·E) 候補、決定性のためのカノニカルランクによるタイブレーク）を用いており、真に最小重み・決定的な基底を返します。5,000分子ChEMBLサブセット・worst-of-10走査で測定: 自己安定性 **100%**（旧: 50.6%）、単一パースでのRDKitとの環サイズ一致率 **98.9%**（旧: 72.4%）— 残る約1.1%の差はRDKit自身の `GetSymmSSSR` が対称縮合環系（例: キュバン、μ=5に対しRDKitは6環を返す）でトポロジー的最小数より多い環を正当に返すことによるもので、chematic側のバグではありません（完全な対称化＝Vismara relevant cyclesは将来課題）。下流への効果（同コーパス）: 環サイズSMARTS `[r5]`/`[r6]` **0%不安定**（旧: 29〜55%）、`NumAromaticRings` **0%**（旧: 約4%）、`RingCount`/MW/TPSA/HBA/HBD/LogP/MRは元々0%のまま。旧SSSRのバグが偶然、別の未解決な芳香族性バグを覆い隠していた2件については下記の芳香族性モデルの項を参照してください。詳細: `scripts/ringinfo_parity.py`。
- **Murckoスキャフォールド: 環トポロジーと正規化文字列出力はいずれも完全に安定**: 以前報告した「100%不安定」自体が測定ハーネスのバグでした（`Mol` オブジェクトをPythonの同一性で比較しており、実際の結果に関わらず常に「不安定」と判定していた）。このスクリプトバグは修正済みです（`scripts/ring_collateral_damage.py`）。上記のカノニカルSMILES破損修正後に5,000分子・worst-of-10で再測定した結果、正規化後（`apply_aromaticity().canonical_smiles_mode("nostereo")`）の自己安定性は**100%（0/5000が不安定）**に到達しました（旧: 0.8%残差）— この残差が上記と同じカノニカルSMILESの構造的破損であったことが確認され、今回で完全に解消しました。正規化なしの生の同位体的 `scaffold().smiles` 文字列比較は**79.30%**安定（不安定率20.70%、旧: 約45%。上記のE/Z部分正規化修正による変化はほぼなし — スキャフォールドはその修正が効く側鎖モチーフの大半を除去してしまうため）— 残りはMurcko固有の問題ではなく、上記の部分的に未解決な `/`/`\` 方向正規化ギャップ（別問題）によるものです。`scaffold()` は正しい環系を確実に抽出します。走査順が異なる入力間で文字列として比較したい場合は、生の `.smiles` ではなく `mol.apply_aromaticity().canonical_smiles_mode("nostereo")` を使ってください。
- **芳香族性モデル**: Hückel 4n+2 則を各 SSSR 環に独立適用（RDKit は縮合環電子非局在化モデルを使用）。N-ヘテロ環で差異あり。4,999 分子 ChEMBL コーパスの現状: HBA/HBD/芳香環数 **100%**、TPSA **100%**（±0.1 Å²）。Kekulized入力でのworst-of-10芳香族性パリティ: **96.3%**（`scripts/aromaticity_atom_parity.py`）— 上記SSSR修正後も数値はビット単位で不変（SSSRのバグと芳香族性のギャップが独立していたことを確認）。芳香族性ギャップの根本原因は別途 `aromatic_context` バイパス機構にあり、未修正です。SSSR修正により2分子（アズレン、プリン）が明確に退行したことが判明しています — 旧来の壊れたSSSRが、これら非交互環系・架橋頭部を多く含む構造に対して `aromatic_context` のバグを偶然覆い隠していたためです。この2分子は5,000分子の測定コーパスに**そもそも含まれていません**（直接検索で確認済み — ChEMBL由来の薬様コーパスには裸のアズレン・プリンは含まれない）。したがって96.3%という数値が不変なのはこの2分子を「見ていない」からであり、「影響がゼロ」だからではありません。両分子は根本原因をコード内に記載した上で `chematic-perception` のテストスイートに `#[ignore]` 付き既知回帰として記録済みで、`aromatic_context` 修正を待っています。
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
│   ├── chematic-crystal/    周期結晶構造: 格子、PBC、近傍探索、supercell（Molecule非依存）
│   └── chematic/            フィーチャーフラグ付きアンブレラクレート（統合クレート）
```

---

## ライセンス

Apache License 2.0 または MIT License のいずれかで利用可能。

---

chematic が役に立ったら、[GitHub スター](https://github.com/kent-tokyo/chematic)をいただけると他の方への発見につながります。
