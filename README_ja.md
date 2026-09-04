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
デフォルトはPure Rust · optional native InChI C FFI · Python · WebAssembly · [公式サイト](https://chematic.io/) · [ライブデモ](https://kent-tokyo.github.io/chematic/playground/)

### v1.0.5 の対応範囲

v1.0.5 は、v1.0.0 の bounded な互換性契約を維持しつつ、typed reaction
document、document-level CDXML 編集、明示的で bounded な Markush/polymer
展開、結晶組成集計、安全性を高めた UFF rescue、canonical/SDF
ホットパス改善を追加します。完全な任意構造 CDXML 編集、
複雑な Markush/polymer 展開、完全な RDKit `RWMol` 互換、3D の完全な
ETKDG/MMFF94 互換は対象外です。再現可能なローカル gate は
[v1.0 local release gate](docs/v1.0-local-release-gate.md) を参照してください。
アルゴリズムと第三者由来コードの境界は
[実装 provenance](docs/implementation-provenance.md) に記録しています。
Spectrophores は patent/FTO 状態が独立に確認されるまで Rust/Python API
から意図的に撤去しています。測定条件と限定範囲は[ベンチマーク文書](docs/benchmark.md)
に記録しています。

| | chematic | RDKit (Python) | RDKit.js (WASM) |
|---|---|---|---|
| **導入方法** | `pip install chematic` | `pip install rdkit`（公式prebuiltホイール）または conda | `npm install @rdkit/rdkit`、Python バインディングなし |
| **ブラウザ向けバンドル** | **3.30 MB raw / 1.21 MB gzip** | 該当なし（Python/C++ライブラリ） | 6.91 MB raw* |
| **ECFP4バッチ** | **54.7 µs/mol** | 94.3 µs/mol | — |
| **Canonical SMILES** | **24.95 / 18.27 µs/mol** | 25.58 / 26.82 µs/mol | — |
| **SDF graph read / serialization-only write** | **9.48 / 7.62 µs/mol** | 99.96 / 79.54 µs/mol | — |
| **メモリ安全性** | コンパイラが保証（Rust） | C++ | C++ |
| **ソースビルド** | `cargo build` のみ | cmake + clang + Boost | Emscripten SDK |

\* RDKit.js の gzip転送時サイズは未計測のため、rawサイズ同士で比較している。RDKit.js は
現在メンテナ移行中(詳細は同リポジトリを参照)。

canonical/SDF の行は 2026-09-04 macOS arm64 の中央値であり、記録した corpus と
処理境界に限定されます。[ベンチマーク詳細](https://kent-tokyo.github.io/chematic/benchmark/)を参照。
chematic の WASM サイズは v1.0.2 リリース候補を `wasm-pack 0.13.1` +
`wasm-opt 130 -O3` でビルドして 2026-09-04 に計測: **3.30 MB raw**(**1.21 MB gzip**)。
比較対象は履歴として固定した RDKit.js **6.91 MB**
(`@rdkit/rdkit@2025.3.4-1.0.0`の`RDKit_minimal.wasm`、unpkg.com で確認) · Indigo(Ketcher向けビルド)
**11.24 MB**(`indigo-ketcher@1.45.1`のメイン`.wasm`、jsDelivr で確認) — chematic の raw WASM
バイナリは現在、RDKit.js よりおよそ2.1倍、Indigo の Ketcher向けビルドよりおよそ3.8倍小さい
(raw同士の比較)。詳細は[artifact 記録](benchmarks/2026-09-04-wasm-size.md)を参照。

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

- ブラウザで化学計算を動かしたい（WASM、1.21 MB gzip、サーバー不要）
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
# chematic v1.0.5
# Python 3.12.x  |  darwin arm64
#
# Descriptor accuracy (2026-08-23, v0.18.0 vs RDKit 2026.03.4):
#   MW                    99.82% within ±0.01 Da
#   HBA / HBD / ARC       100%   (4,999-mol ChEMBL subset)
#   TPSA                  100%
#   LogP (Crippen)        100%*  (max Δ = 1.1×10⁻¹³)
#   CIP R/S/E/Z           99.74% (opt-in accurate engine)
# ...
```

---

## AI / LLM 開発者向け

chematic はローカルAI agent連携向けの **MCP（Model Context Protocol）サーバー**を提供します。

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

common chemistry coreはsafe Rustで、公開されたuntrusted-input pathには有限のdefaultとtyped errorがあります。`native-inchi`はIUPAC InChI C libraryを使う明示的なopt-in FFI例外です。依存crate自身のunsafe codeは別境界として扱います。詳細は[SECURITY](SECURITY.md)を参照してください。

記録済みv0.18.0 ECFP4バッチ中央値は、同じ5,000分子corpusとApple M4環境で54.7 µs/mol、RDKitは94.3 µs/molでした。WASM artifactはv1.0.2候補で3.30 MB raw / 1.21 MB gzipです。どちらも日付と条件を固定した測定であり、一般性能保証ではありません。

## 比較の考え方

| 観点 | chematic | RDKit |
|---|---|---|
| 導入 | pure-Rust core、Python wheel、WASM/Node | 最大規模の機能・ecosystem |
| Browser | native WASM | RDKit.jsは別distribution |
| 3D/MMFF94 | Experimental | 成熟したETKDG/force-field workflow |
| 互換性 | 名前付きsubsetとtyped failure | 基準実装 |

機能差と非対応範囲は[RDKit比較](docs/rdkit-comparison.md)と[互換性範囲](docs/compatibility-scope.md)、速度は[benchmark](docs/benchmark.md)を参照してください。

---

## JavaScript / TypeScript（WebAssembly）

**1.21 MB gzip — RDKit.js の raw WASM と比べておよそ2.1分の1。** Emscripten・cmake 不要。ブラウザ・Node.js どちらでも動作。

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

| 分野 | crate |
|---|---|
| 分子グラフと同一性 | `chematic-core`, `chematic-smiles`, `chematic-perception`, `chematic-cip` |
| 検索・記述子・fingerprint | `chematic-smarts`, `chematic-chem`, `chematic-fp` |
| ファイル・反応 | `chematic-mol`, `chematic-rxn`, `chematic-inchi`, `chematic-iupac` |
| 2D・3D・材料 | `chematic-depict`, `chematic-3d`, `chematic-ff`, `chematic-crystal`, `chematic-ewald` |
| 利用者向けinterface | `chematic`, `chematic-py`, `chematic-wasm`, `chematic-cli`, `chematic-mcp` |

詳細は[形式対応表](docs/format-capabilities.md)、[言語binding](docs/language-bindings.md)、各crateのREADMEを参照してください。

---

## 最近の開発

**未リリース:** #210 のlegacy UFF stereo rescue残差を解消し、canonical SMILESとSDFのhot pathを改善しました。固定条件の測定は[benchmarks](benchmarks/)に記録しています。

**v1.0.5（2026-09-05）:** #210 のlegacy UFF stereo-rescue残差を解消し、canonical SMILESとSDFのhot pathを改善しました。v1.0.4の機能も現行リリースに含まれます。Spectrophoresは独立したpatent/FTO確認まで公開APIから除外しています。

公開リリースの要約は[CHANGELOG](CHANGELOG.md)を参照してください。詳細な開発記録は、
そこからリンクされたarchiveに保持しています。

---

## 既知の制限事項

- `canonical_smiles()` は表現であり、dedup/cache keyではありません。`canonical_smiles_stable_key()` の `None` を処理してください。
- Aromaticity/CIPはdefaultとopt-in modelを明示的に分離し、universal RDKit parityを主張しません。
- 3D生成とMMFF94はExperimentalです。
- Python `RWMol`、CDXML編集、Markush/polymer expansionは意図的なbounded subsetです。
- pure-Rust InChIは近似です。標準IUPAC InChIには`native-inchi`を使います。

正確な契約は[互換性範囲](docs/compatibility-scope.md)、[検証](docs/validation.md)、[errorとresource limit](docs/error-and-limits.md)を参照してください。

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
著作権表示: Kentaro Tanabe (kent-tokyo)。再配布時の帰属表示は [`NOTICE`](NOTICE) を
参照してください。

---

chematic が役に立ったら、[GitHub スター](https://github.com/kent-tokyo/chematic)をいただけると他の方への発見につながります。
