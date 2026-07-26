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
Pure Rust · C/C++ ゼロ · Python · WebAssembly · [ライブデモ](https://kent-tokyo.github.io/chematic/playground/)

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **導入方法** | `pip install chematic` | conda / cmake が必要 | Python バインディングなし |
| **ブラウザ向けバンドル** | **719 KB** | 提供なし | ~30 MB（~42× 大きい） |
| **バッチ FP 速度** | **~78 µs/mol**（2–3× 高速） | ~160–235 µs/mol | — |
| **メモリ安全性** | コンパイラが保証（Rust） | C++ | C++ |
| **ソースビルド** | `cargo build` のみ | cmake + clang + Boost | Emscripten SDK |

すべての数値は再現可能です — [ベンチマーク詳細](https://kent-tokyo.github.io/chematic/benchmark/)を参照。  
WASM サイズ: chematic **719 KB** · RDKit.js ~30 MB · Indigo WASM ~40 MB

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

- ブラウザで化学計算を動かしたい（WASM、719 KB、サーバー不要）
- C++ ツールチェーンなしの Pure Rust スタックが必要
- `pip install rdkit` が困難な環境（Cloudflare Workers、Lambda、組み込み）にデプロイする
- AI エージェントを構築し、ネイティブな MCP ツール統合が必要
- バッチ処理で高スループットが必要（ECFP4: RDKit の 2〜3 倍高速、Rayon 並列）
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
# chematic v0.7.0
# Python 3.12.x  |  darwin arm64
#
# Descriptor accuracy (benchmark 2026-06, v0.7.0 vs RDKit 2026.03.3):
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

---

## Pure Rust の理由

### 速い

Rust のゼロコスト抽象化と所有権モデルはオーバーヘッドをソースレベルで排除します。
chematic の ECFP4 フィンガープリントバッチは多様な分子コーパスで **~78 µs/mol** — 同じハードウェアで
RDKit Python API の 2〜3× 高速（全 CPU コアで Rayon 並列化）。GIL なし、インタープリタオーバーヘッドなし、
`_sys` クレート内の FFI 呼び出しコストなし。

### 安全

chematic 自身の約 15,000 行の Rust コードには **`unsafe` ブロックが約 6 個**のみで、
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
npm パッケージ `@kent-tokyo/chematic` は **719 KB gzip** — RDKit.js の約 42 分の 1。
1 つのコードベースが Linux・macOS・Windows・あらゆるブラウザで動作します。

---

## 他のケモインフォマティクスライブラリとの比較

| 観点                                        | **chematic**                               | RDKit (rdkit-sys)  | OpenBabel FFI | RDKit.js (WASM)  |
|---------------------------------------------|--------------------------------------------|--------------------|---------------|------------------|
| **C/C++ 依存**                              | **ゼロ（デフォルト）**†                    | 大規模 C++         | 大規模 C++    | C++（Emscripten）|
| **WASM バイナリサイズ**                     | **〜550 KB**                               | N/A（WASM 非対応） | N/A           | 〜30 MB          |
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
| **pKa 予測**                                | **あり（15 SMARTS ルール）**               | なし               | なし          | なし             |
| **ADMET プロファイル + BOILED-Egg**         | **あり**                                   | 一部               | なし          | 一部             |
| **MCP サーバー（AI エージェント API）**     | **あり — 20 ツール（Name→SMILES 含む、stdio のみ）**  | なし               | なし          | なし             |
| IUPAC 名生成                                | **あり（25+ 化合物クラス）**               | なし               | なし          | 一部             |
| メンテナンス（2026）                        | アクティブ                                 | アクティブ         | 最小限        | アクティブ       |

† デフォルトビルドのみ。`native-inchi` feature は opt-in で C コンパイラが必要（IUPAC InChI C ライブラリ v1.07.5 の vendoring）。これは C/C++ FFI 固有の話 — 下記の `depict` feature は純 Rust の描画クレートを引き込むため、unsafe フリーではなくても C コンパイラ依存は追加しません（‡参照）。

‡ chematic 自身の約 15,000 行の Rust コード: `native-inchi` の約6個の FFI ブロックを除き unsafe フリー（上記「安全」参照）— chematic 自身が書いたコードについての実測済みの主張であり、コンパイラによるチェックが一切効かない RDKit/OpenBabel の C++ FFI unsafe とは、たとえ個数が同程度でも種類が根本的に異なります。**依存ツリー全体については成り立ちません**: opt-in の `depict` feature（SVG/PDF/EPS 描画）は resvg/usvg/rustybuzz/tiny-skia/zune-jpeg を引き込み、これらは純 Rust ですが unsafe フリーではありません — 実測（`unsafe fn`/`impl`/`trait`/`{` の出現数）: tiny-skia 151、zune-jpeg 79、rustybuzz 14、image 8、fontdb 3、tiny-skia-path 3(この範囲だけで合計 258)。`chematic-py`（`pip install chematic`）と npm パッケージはどちらも `chematic-depict` に直接依存するため、これは実際の2つのインストール経路の両方に当てはまります。


</details>
---

## JavaScript / TypeScript（WebAssembly）

**719 KB gzip — RDKit.js の約 42 分の 1。** Emscripten・cmake 不要。ブラウザ・Node.js どちらでも動作。

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
| `chematic-mol`         | MOL/SDF V2000+V3000（R/W、2D 座標付き）、CML（R/W）、CDXML（R）；`SdfRecord`（coords+props）、MDL RXN V2000 読み書き；V3000 ステレオグループ COLLECTION R/W；**2Dウェッジ/ハッシュのtetrahedral parity + E/Z二重結合方向を読み込み時に自動認識**（`read_mol_with_diagnostics`/`read_mol_v3000_with_diagnostics`、型付きopt-in診断） | 130+     |
| `chematic-depict`      | 2D SVG（CPK カラー・ハイライト・グリッド）、`detect_crossings`・`render_svg_with_metadata`・反応 SVG；Y座標系ドキュメント整備  | 64      |
| `chematic-chem`        | 190+ 記述子値（71 関数）、タウトマー、スキャフォルド、BRICS、QED、標準化；**pKa 予測** (15 SMARTS ルール)；**ADMET プロファイル** (BBB/Caco-2/hERG/CYP3A4)；**HBA 99.98% RDKit 一致率**（4,999 分子 ChEMBL ベンチマーク）；**TPSA ±0.1 Å² 98.1% / LogP ±0.01 96.5% / HBD 100%** RDKit 一致 | 662     |
| `chematic-fp`          | ECFP2/4/6、FCFP4/6、MACCS、TopoPF、AtomPair、Torsion、Layered、Pattern、Pharmacophore、Reaction、**MAP4** (Minervini 2020) — Tanimoto/Dice | 185      |
| `chematic-ff`          | **MMFF94 全 7 エネルギー項** (Halgren 1996)：OOP (117件) + Stretch-Bend (282件)；steepest descent + L-BFGS；DREIDING | 98      |
| `chematic-smarts`      | SMARTS、VF2、MCS；**SmartsCache** (LRU 5–20×)；**named_pattern()** (20 パターン)；**SMARTS 内アトムマップ `:N`** (`[O;D1;H0:3]` 形式 — メタデータとして保存、マッチング条件には不使用) | 142     |
| `chematic-3d`          | 3D 座標生成、ETKDG KB (40 パターン、adaptive noise)、力場最小化、形状記述子、ConformerEnsemble、PDB/XYZ | 265     |
| `chematic-rxn`         | 反応 SMILES/SMIRKS、`run_reactants`/`run_reactants_strict`；**`retro_disconnect()`** — 60 retro-SMIRKS テンプレート (AmideBond/Ester/Ether/CNBond/CCBond/CSBond) + SA Score ランク付き | 137      |
| `chematic-inchi`       | InChI/InChIKey：純 Rust 近似（WASM 対応）**+ `native-inchi` feature で IUPAC 標準準拠**（C ライブラリ 1.07.5 vendored、ビット完全一致）；**parse_inchi** 読み込み；**検証付きcanonical SMILES重複排除**（`dedup`モジュール、legacy CIPで未解決の指定済みtetrahedral stereoに対してfail-closed） | 96 (+16*)   |
| `chematic-cip`         | opt-inの高精度CIPエンジン（`assign_cip_accurate_experimental`、階層的digraph、Rules 1a/1b/2/4b/5、RDKit互換MANCUDE分数原子番号）— デフォルトの`assign_cip()`/`CipMode::LegacyFast`は変更なし | —       |
| `chematic-wasm`        | **130+ WASM エクスポート** — npm: `@kent-tokyo/chematic`（公開版は`0.5.0`；crates.io/PyPIは`0.6.0`まで進んでおり、npm公開が遅れている既知のギャップ）；**pKa/ADMET/BBB/Caco-2/hERG/CYP3A4** WASM API | 211     |
| `chematic-iupac`       | ローカル IUPAC 命名（Pure Rust・オフライン）— **25+ 化合物クラス**：アルカン、シクロアルカン、アルコール、アミン、ハロアルカン、ケトン、酸、エステル、アミド、**ピペリジン、モルホリン、ピペラジン、ナフタレン、スルフィド** | 47      |
| `chematic-mcp`         | **MCP (Model Context Protocol) サーバー** — AI エージェント統合；**20 ツール**：parse_smiles, calc_properties, ecfp4, tanimoto, smarts_match, canonical_smiles, find_mcs, generate_3d, pains_check, brenk_check, sa_score, admet_profile, boiled_egg, lipinski_check, name_to_smiles, retrosynthesis, smiles_to_moljson, moljson_to_smiles, representation_router, molecule_context_pack | 31      |
| `chematic`             | フィーチャーフラグ付きアンブレラクレート（統合クレート）                                                                                                  | 1       |

```
cargo test --workspace --lib --quiet                                               # 2,746 ライブラリテスト、全パス（2026-07-26 時点）
cargo test -p chematic-inchi --features native-inchi --test standard_inchi         # +16 IUPAC 標準 InChI 統合テスト
```

---

## 最近の開発

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

---

## 既知の制限事項

- **`canonical_smiles()` はE/Zステレオ化学に対して部分的に正規化されました — それでも重複排除やキャッシュキーとしてはまだ安全ではありません**: 孤立した単純なE/Z二重結合には `/N=N/` と `\N=N\` のような2通りの等価な正しい記法があり、ライターは従来どちらを出力するか正規化していませんでした。一般ケースを修正済み: ある二重結合と、それに幾何学的に連動する全ての方向性結合(共役鎖全体を含む)からなる「連結E/Z系」ごとに、カノニカル書き込み順で最初に現れる方向性結合が常に `/` になるよう正規化するようにしました(入力の綴りに依存しません)。5,000分子ChEMBLコーパス・worst-of-10で測定: E/Z限定の自己不安定性(四面体ステレオを除去)は**9.76%→5.50%**に改善(275/5000が依然不安定)。この変更による構造的正しさへの影響はなし(ChEMBL **0/5000**、非環式ポリエンコーパス **0/33** を再検証済み)。残る275件は**100%が見た目上の問題のみ**(RDKit的にはどの変異体も同一分子で、破損はゼロ)と確認済みですが、原因は**根本解明しきれていない混在プール**です: 約半数は特定のモチーフ(2つ以上の環外二重結合を持つ小員環、例えば交差共役した環状ジイミン)に一致し、そこでは「1つの系」とみなすべき物理結合の集合が入力の綴りに対して不変になっていません。残り半数は未特定です。これが完全に解消するまでは、ステレオを持つ分子の約18分子に1つ(旧: 約10分子に1つ)が、同一分子に対して2通りの、それぞれ単独では正しい `canonical_smiles()` 文字列を生成しうるため、現時点では重複排除やキャッシュキーとして使わないでください。この点が重要な用途では、当面 `apply_aromaticity()` で正規化した文字列を独自の重複排除キーとして使ってください。
- **カノニカルSMILESの構造的破損 — 修正済み**: 修正前は `canonical_smiles(parse(x))` が `x` の入力走査順によっては(単に綴りが異なるだけの等価な文字列ではなく)別の立体異性体を静かに出力することがありました。5,000分子ChEMBLサブセット・worst-of-10走査・RDKit検証済みの構造正しさで測定: **4.28%（214/5000）** の分子で、少なくとも1つの変異体が誤った分子として読み戻されていました。根本原因は独立した2件のパーサーバグで(当初疑われていた「共役二重結合のマーカーは結合をまたいで幾何学的に連動している」という診断は誤りであったことが判明— 下記参照)、いずれも実在する分子と最小構成の回帰テストで確認済みです: (1) 環closure(`/`/`\`) の方向マーカーを閉環側の出現位置で読み取る際、開環→閉環方向へのフリップをせず生の値のまま保存しており、共役E/Z鎖の連結結合がたまたま環closureを経由する場合に破損していた。(2) 自身の分岐の中で環closureの相手が閉じるステレオ中心は、再利用可能な環番号ではなく出現ごとに一意なIDで隣接原子順序の解決を行うべきところ、番号ベースで解決していたため、SMILES文字列の後方で同じ番号が無関係な環に再利用された際に、ステレオ中心の隣接原子順序が静かに乗っ取られ破損していました。**両修正後: ChEMBLコーパスで構造正しさは100%（0/5000）** — これは、rankingという第3の無関係な修正の有無を含む独立した3通りの再構成手順で3重に確認済みです。両方の根本原因が環closure固有のものである一方、レチノイド・カロテノイド・プロスタグランジン・ロイコトリエン・ポリエン系マクロライドは長い共役系を**非環式**の鎖として持っており、ChEMBLのランダムサンプリングにはほぼ含まれないため、これら5クラスの実在化合物33種からなる専用コーパス(トレチノイン、β-カロテン、リコペン、アンフォテリシンB、ロイコトリエンB4など、`scripts/polyene_corpus.csv`)でも独立に再検証しました: **worst-of-30で0/33（0.00%）**。同一コーパスの未修正コードでは12/33（36.36%）の破損を陽性対照として確認済みです(12件全てが環closureを多く含む構造で、完全非環式のリコペンを含め純粋な非環式の例は未修正コードでも一度も破損しませんでした)。これにより当初の「共役系全般が原因」という診断が誤りであったことが直接確定し、残存する破損クラスは見つかっていません。骨格限定・四面体限定の自己**安定性**もいずれも0%に到達(旧: 0.16%、4.36%)。全ステレオを含めた生の自己安定性は**86.02%→90.28%**（不安定率13.98%→9.72%）に改善 — 残差は全て上記の非破損な方向正規化ギャップによるもので、破損の残存ではありません。往復不変性（`canonical(parse(canonical(m))) == canonical(m)`）はわずかに改善（**98.26%→98.32%**）— この指標はそもそも破損クラスを直接測定していなかったためです。
- **環知覚（SSSR）は非決定的・非最小だったが修正済み**: 旧 `find_sssr` は単一の全域木から非木辺ごとに基本閉路を1つずつ生成していたため、木の形によっては最小でない環（例: ナフタレン `c1ccc2ccccc2c1` が決定的に `[6,10]` を返す。正しくは `[6,6]`）を返していました。現在の `find_sssr` は Horton アルゴリズム（全頂点×全辺の最短路木から候補閉路を生成、O(V·E) 候補、決定性のためのカノニカルランクによるタイブレーク）を用いており、真に最小重み・決定的な基底を返します。5,000分子ChEMBLサブセット・worst-of-10走査で測定: 自己安定性 **100%**（旧: 50.6%）、単一パースでのRDKitとの環サイズ一致率 **98.9%**（旧: 72.4%）— 残る約1.1%の差はRDKit自身の `GetSymmSSSR` が対称縮合環系（例: キュバン、μ=5に対しRDKitは6環を返す）でトポロジー的最小数より多い環を正当に返すことによるもので、chematic側のバグではありません（完全な対称化＝Vismara relevant cyclesは将来課題）。下流への効果（同コーパス）: 環サイズSMARTS `[r5]`/`[r6]` **0%不安定**（旧: 29〜55%）、`NumAromaticRings` **0%**（旧: 約4%）、`RingCount`/MW/TPSA/HBA/HBD/LogP/MRは元々0%のまま。旧SSSRのバグが偶然、別の未解決な芳香族性バグを覆い隠していた2件については下記の芳香族性モデルの項を参照してください。詳細: `scripts/ringinfo_parity.py`。
- **Murckoスキャフォールド: 環トポロジーと正規化文字列出力はいずれも完全に安定**: 以前報告した「100%不安定」自体が測定ハーネスのバグでした（`Mol` オブジェクトをPythonの同一性で比較しており、実際の結果に関わらず常に「不安定」と判定していた）。このスクリプトバグは修正済みです（`scripts/ring_collateral_damage.py`）。上記のカノニカルSMILES破損修正後に5,000分子・worst-of-10で再測定した結果、正規化後（`apply_aromaticity().canonical_smiles_mode("nostereo")`）の自己安定性は**100%（0/5000が不安定）**に到達しました（旧: 0.8%残差）— この残差が上記と同じカノニカルSMILESの構造的破損であったことが確認され、今回で完全に解消しました。正規化なしの生の同位体的 `scaffold().smiles` 文字列比較は**79.30%**安定（不安定率20.70%、旧: 約45%。上記のE/Z部分正規化修正による変化はほぼなし — スキャフォールドはその修正が効く側鎖モチーフの大半を除去してしまうため）— 残りはMurcko固有の問題ではなく、上記の部分的に未解決な `/`/`\` 方向正規化ギャップ（別問題）によるものです。`scaffold()` は正しい環系を確実に抽出します。走査順が異なる入力間で文字列として比較したい場合は、生の `.smiles` ではなく `mol.apply_aromaticity().canonical_smiles_mode("nostereo")` を使ってください。
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
