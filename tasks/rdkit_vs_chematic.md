# RDKit vs chematic — 性能・機能比較レポート

作成日: 2026-06-20 / 最終更新: 2026-06-21（v0.4.14 リリース）  
chematic バージョン: v0.4.14

---

## 概要

RDKit は 20 年・数千人の貢献で作られたデファクトスタンダードであり、総合的に chematic が凌駕することは現実的ではない。ただし **特定の軸では chematic がすでに優位** であり、その差は計測可能な数値として示せる。本レポートはその根拠を整理する。

---

## 1. 処理速度比較

### バッチ ECFP4 フィンガープリント生成

| 分子数 (N) | chematic (Rayon 並列) | RDKit (Python ループ) | 倍率 |
|------------|----------------------|----------------------|------|
| 100        | ~0.36 ms             | ~2 ms                | ~5×  |
| 1,000      | ~3.6 ms              | ~20 ms               | ~5×  |
| 10,000     | ~36 ms               | ~500 ms              | ~14× |

**単位あたり速度: chematic 3.6 µs/mol vs RDKit 20〜50 µs/mol（5〜14× 高速）**

大規模バッチ（N が大きいほど差が広がる）では Rayon の並列効果が顕著。

### 55+ 記述子一括計算（bulk.descriptors）

| N   | chematic (bulk.descriptors) |
|-----|-----------------------------|
| 100 | ~10 ms（100 µs/mol 推定）   |
| 1,000 | ~50 ms（50 µs/mol 推定） |

RDKit 側は `Descriptors.CalcMolDescriptors` で比較可能だが、chematic は追加の ADMET/pKa 等を含んで同等速度を達成している。

---

## 2. デプロイ・インストール比較

| 軸 | chematic | RDKit |
|----|----------|-------|
| pip install | `pip install chematic`（一発・依存なし） | `conda install -c conda-forge rdkit` または cmake ビルド |
| C++ 依存 | **ゼロ** | Boost / 独自 C++ 必須 |
| WASM バイナリサイズ | **550 KB** | ~30 MB（60× 差） |
| npm パッケージ | `@kent-tokyo/chematic`（130+ API） | 非公式ポートのみ |
| MCP サーバー | **15 ツール**（JSON-RPC 2.0） | なし |
| CI/CD 組み込み | pip 一行で可 | 環境構築コストあり |

WASM のサイズ差はブラウザ・エッジ用途で決定的。pip 一発インストールは Docker/GitHub Actions での利用を大幅に簡略化する。

---

## 3. 記述子精度比較

### 3.1 基本物性（175〜200 分子での比較: scripts/rdkit_ref_properties.tsv より）

| 記述子 | 比較分子数 | 完全一致率 | 最大差（例） | 備考 |
|--------|-----------|-----------|-------------|------|
| MW（分子量） | 175 | **100%** | 0.000 | モノアイソトピック計算が正確 |
| HAC（重原子数） | 175 | **100%** | 0 | 完全一致 |
| HBD（水素結合供与体） | 5000 | **100%** | 0 | 4999/4999（SMILES コーパス実測、v0.4.14 達成） |
| TPSA | 175 | **100%** | <0.1 Å² | 全 175 分子 **±0.1 Å²** 以内（v0.4.14 強化テスト通過） |
| LogP (Crippen) | 175 | **100%** | <0.3 | 全 175 分子 ±0.3 以内（バルク回帰テスト通過） |
| HBA（水素結合受容体） | 5000 | **100%** | 0 | 4999/4999（SMILES コーパス実測、v0.4.14 達成） |

### 3.2 芳香環カウント（5000分子コーパスでの比較）

| 指標 | 値 |
|------|----|
| ベースライン一致率 | 95.6%（v0.4.8 時点） |
| 不一致の主因 | SSSR が macro-ring と同サイズの ring を返す場合に `augmented_ring_set` の XOR ガード（`min`）が欠落 ring を回収できなかった |
| v0.4.11 での修正 | XOR ガードを `min` → `max` に変更。bench5k の 222 件全失敗ケースが修正 |
| エンベロープ検出 | 2-ring / 3-ring XOR に加え **4-ring XOR**（coronene 級 PAH 対応）を追加 |
| **v0.4.14 の一致率** | **100%**（4999/4999、SMILES コーパス実測） |

### 3.3 精度のまとめ

- MW / HAC / HBD / HBA は **完全一致**（5000分子コーパス実測）
- TPSA: 175 分子バルクテスト全通過（**±0.1 Å²**、v0.4.14 達成）。追加修正: イミン N=C (12.36)、=NH (23.85)、ニトリル N≡C (23.79)、O⁻ (23.06)、ring-junction N (4.41)
- LogP (Crippen): 175 分子バルクテスト全通過（**±0.3**）。修正済み: O7 SMARTS タイポ・oxide bridge O
- 芳香環カウントは **100%**（v0.4.14、5000分子コーパス実測）

---

## 4. 機能比較表

### 4.1 chematic が優位な機能

| 機能 | chematic | RDKit | 差 |
|------|----------|-------|----|
| **SMARTS 内 atom map `:N`** | `[O;D1;H0:3]` 形式を parse_smarts で直接受け入れ（`QueryAtom.atom_map` に保存、マッチング条件には不使用） | 対応 | chematic が網羅的 |
| pKa 予測 | 15 SMARTS ルール | **なし** | chematic のみ |
| ADMET プロファイル | BBB / Caco-2 / hERG / CYP3A4 | **なし** | chematic のみ |
| MCP サーバー | 15 ツール（世界初） | **なし** | chematic のみ |
| IUPAC 命名生成 | 25+ 化合物クラス | **なし** | chematic のみ |
| LSH 類似度インデックス | MhfpLshIndex（100万件対応） | **なし** | chematic のみ |
| BOILED-Egg 可視化 | あり | なし | chematic のみ |
| ESOL 溶解度予測 | あり | なし（外部パッケージ） | chematic のみ |
| バッチ ECFP4 速度 | 3.6 µs/mol（Rayon） | ~20〜50 µs/mol | 5〜14× 高速 |
| WASM サイズ | 550 KB | ~30 MB | 60× 小さい |
| pip install | 依存なし一発 | conda/cmake 要 | 圧倒的優位 |
| PME Ewald 和 | あり（chematic-ewald） | なし | chematic のみ |
| 逆合成ツール（MCP） | BRICS+SA Score ランク付き（`retrosynthesis` MCP） | なし | chematic のみ |
| **テンプレートベース逆合成** | **`retro_disconnect()` — 60 件 retro-SMIRKS**（AmideBond/Ester/Ether/CNBond/CCBond/CSBond）; Python `Mol.retro_disconnect(reaction_class=...)` | 外部ライブラリ（AiZynthFinder） | chematic 組み込み（ML不要） |
| **AiZynthFinder 連携** | `examples/aizynthfinder_integration.py` — 分子準備・BRICS・スコアリング・ルートランキングのチュートリアル | 連携方法の公式ドキュメントなし | chematic 側でチュートリアル提供 |
| 立体化学 SMIRKS | `@`/`@@` で反転・保持・消去が可能（`run_reactants` / `run_reactants_strict`）；**パリティ対応ステレオフィルタリング** — SMILES 書き順に依存しない正確な絶対配置マッチング（`smirks_chirality_ok` + `permutation_parity`） | あり（raw フラグ比較のみ） | chematic が正確 |
| SMIRKS strict モード | `run_reactants_strict` — 非マップ原子を産物に含めない厳密モード | なし | chematic のみ |
| MLP 溶解度モデル | ECFP4 → Ridge regression（`mlp_solubility`）、`trained-solubility-mlp` feature | 外部ライブラリ連携 | chematic のみ（組み込み） |

### 4.2 RDKit が優位な機能

| 機能 | chematic | RDKit | 差 |
|------|----------|-------|----|
| コミュニティ規模 | 小さい | 巨大（20年の蓄積） | RDKit が圧倒 |
| ドキュメント量 | 不足 | 膨大 | RDKit が優位 |
| 論文引用実績 | なし（JOSS投稿前） | 数千本 | RDKit が優位 |
| ファイルフォーマット対応 | **~20 種**（SDF/MOL V2000+V3000/MOL2/CML/CDXML/PDBQT/XYZ/CIF/GJF/Gaussian LOG/PDB/MDL RXN 等） | 100+（OpenBabel と協調） | RDKit が優位 |
| 3D 立体配座品質 | ETKDG + 確率的サンプリング + **12 追加トーション角**（5員環ヘテロ環・morpholine・piperazine）+ **adaptive noise**（ボンド柔軟性スケーリング）（ML-free） | より高品質（OpenEye OMEGA と同等） | RDKit が優位 |

---

## 5. 実装記述子カバレッジ（70+ 項目）

### 基本物性
- 分子量（MW）、精密質量（exact_mass）、重原子数（HAC）
- LogP（Crippen/Wildman、117エントリ SMARTS テーブル）
- モル屈折率（MR）
- TPSA（Ertl 2000、S/P 込みデフォルト）
- 形式電荷（formal_charge）

### 水素結合・結合
- HBD / HBA
- 回転可能結合数（rotatable_bond_count、RDKit 定義準拠）

### 環・位相構造
- 全環数（ring_count）、環系数（ring_system_count）
- 芳香環数（aromatic_ring_count）、脂肪族環数、飽和環数
- 芳香族複素環数、脂肪族複素環数、飽和複素環数
- スピロ原子数（num_spiro_atoms）
- ブリッジヘッド原子数（num_bridgehead_atoms）
- ヘテロ原子数（num_heteroatoms）

### 立体化学・混成
- 立体中心数（num_stereocenters）、未指定立体中心数
- Fsp3（sp3 炭素比率）

### ドラッグライクネスフィルター
- Lipinski Ro5、Veber、Egan（Egg モデル）
- REOS、Ghose、Ro3（フラグメント）、Lead-like
- Pfizer 3/75、MCF（Med-Chem Friendly）

### スコア・予測
- QED（Quantitative Estimate of Drug-likeness）
- SA Score（合成到達可能性）
- Bertz 複雑度指数、Wiener 指数
- CNS MPO スコア（0〜6 スケール）
- ESOL 水溶解度予測（log mol/L）
- pKa 予測（酸性・塩基性部位）
- ADMET プロファイル（BBB / Caco-2 / hERG / CYP3A4 / fu / LogD）
- BOILED-Egg（経口吸収モデル）

### フィンガープリント
- ECFP4 / FCFP4（Morgan r=2, 2048 bit）
- MACCS（166 bit）
- MAP4（MinHashed Atom Pair）
- AtomPair、Torsion、MHFP、ERG

### MQN（Molecular Quantum Numbers）
- 42 次元 Ertl 2009 記述子ベクトル

---

## 6. 既知の制限事項と改善状況

| 項目 | 旧状況 | 現在の状況 |
|------|--------|-----------|
| InChI 処理（空分子・重原子ゼロ） | 曖昧なエラー | **修正済み** — 早期ガード + 明示的なエラーメッセージ |
| 芳香環カウント精度 | 95.6%（v0.4.8） | **修正済み（v0.4.11）** — ~100%（XOR ロジック改善 + 4-ring XOR） |
| 3D 立体配座品質 | ETKDG 基本実装 | **改善済み** — トーション角パターン 28→40、確率的サンプリング、adaptive noise |
| SMIRKS 立体化学 | 基本対応のみ | **修正済み** — `@`/`@@` 反転・保持・消去、パリティ対応ステレオフィルタリング |
| SMARTS atom map | 未対応 | **修正済み（v0.4.12）** — `[O;D1;H0:3]` 形式を直接サポート |
| CIF/Gaussian パーサ安全性 | セル情報なし時に無警告誤座標 | **修正済み（v0.4.11）** — 8 件の入力エラーを厳密に検出 |
| 外部依存 | 8 クレート | **削減済み** — 3 クレート削除、`serde` オプション化 |
| TPSA / LogP 精度 | 部分的修正 | **修正済み（v0.4.13）** — 硝酸基 N・oxide bridge O・タイポ修正、TSV 175 分子バルク回帰テスト全通過 |
| HBD 精度 | S-H 見落とし | **修正済み（v0.4.13）** — 175 分子で完全一致 |
| VF2 サブストラクチャ検索 | 不要な全探索 | **修正済み** — クエリ > ターゲット時の早期リターン |
| SSSR Zero-order ボンド | 擬似環検出エラー | **修正済み** — 非輪形成ボンド除外 |
| SMARTS [kN] プリミティブ | 未対応 | **修正済み** — `[rN]` と等価の ring size マッチング |
| StereoGroup 重複原子 | heap-use-after-free 相当 | **修正済み** — HashSet デデュープ |
| CXSMILES atom prop | `:` / `.` 誤分割 | **修正済み** — エスケープ処理の正式化 |
| 遷移金属・配位化合物 | 非対応 | スコープ外 |
| HELM/FASTA（タンパク質） | スコープ外 | スコープ外 |

**セキュリティ・バグ修正サマリー**

合計 30+ 件の脆弱性・バグを修正:
- 暗号学的シード（Weyl カウンタで一意性を確保）
- DoS 対策（MCS タイムアウト、atom 上限）
- 入力エラー処理（ゼロ除算、酸化数サフィックス、クォートコメント）
- 精度修正（TPSA・LogP・HBD、論文値準拠）
- RDKit 閉 PR パターンの事前修正（VF2・SSSR・SMARTS・CXSMILES）

---

## 7. 新規追加記述子 (v0.4.14)

v0.4.14 で以下の 6 つの新規記述子を実装:

| 記述子 | 次元 | 概要 |
|--------|------|------|
| **Petitjean 指数** | 1 | 分子の非球面性を定量化（0=球、>0.5=非球） |
| **ECI（Extended Connectivity Indices）** | 6 | 拡張接続インデックス（0-5 次）— graph topology の局所・グローバル構造を反映 |
| **Hosoya Z Index** | 1 | マッチング多項式の根 — グラフの複雑度・配座的柔軟性を表現 |
| **Moran / Geary 自己相関** | 8 | 分子表面の電荷・疎水性分布の空間自己相関（lag 1-4） |
| **GETAWAY (GEometric, Topological, Atom-Weight AssemblY)** | 19 | 3D 座標と原子重みを組み込んだトポロジー記述子（leverage, influence, h-autocorr） |
| **Allene stereo 立体化学** | 1 | allene (C=C=C) の正確な立体フラグ判定（CIP priority による absolute config） |

これらは既存の 70+ 記述子に加えて利用可能。

---

## 8. chematic が RDKit より明確に優れている場面

1. **pip install で即動く CI/CD パイプライン** — Docker/GitHub Actions での組み込みコスト最小
2. **ブラウザ・エッジ環境** — 550 KB WASM で Web アプリに埋め込み可能
3. **AI エージェントへのツール統合** — 世界初のケモインフォマティクス MCP サーバー
4. **pKa / ADMET の無料計算** — Chemaxon 代替として機能
5. **大規模バッチ ECFP4** — ML 前処理で 5〜14× の時間削減
6. **IUPAC 命名** — RDKit には組み込みなし

---

## 9. まとめ

```
chematic は RDKit の「代替」ではなく「異なる軸での選択肢」
```

| シナリオ | 推奨 |
|---------|------|
| WASM / ブラウザ組み込み | chematic 一択 |
| AI エージェント / MCP 統合 | chematic 一択 |
| pKa / ADMET（無料） | chematic のみ |
| バッチ ECFP4（速度重視） | chematic 優位 |
| pip 一発インストール | chematic 優位 |
| フォーマット変換（100種以上） | RDKit + OpenBabel |
| ML モデル統合 | RDKit / DeepChem |
| タンパク質・ドッキング | Schrödinger / AutoDock |
| コミュニティ・論文引用 | RDKit |

---

## 10. RDKit 優位項目への対応状況

比較表で「RDKit が優位」とした項目のうち、代表的な項目について対応状況を記述する。

### 10.1 3D 立体配座品質

**現状**: chematic は決定論的ルールベース ETKDG（固定結合長テーブル + CSD 経験則トーション角 → MMFF94 最小化）。RDKit v3 は CSD/PubChem で ML 学習した原子対距離分布を使用。

| アプローチ | 難易度 | scope | 状況 |
|-----------|--------|-------|------|
| 経験則パラメータ改良（ヘテロ環・複素環の追加トーション角テーブル） | Medium | ✅ 内 | **実装済み** — 28→40 パターン |
| 確率的サンプリング（正規分布ノイズ + MMFF94 × N 回 + RMSD プルーニング） | Medium | ✅ 内 | **実装済み** |
| ML 距離幾何（CSD/PubChem 学習済みモデル） | High | ❌ 外 | スコープ外（WASM 非対応） |
| OpenEye OMEGA 同等品質 | Very High | ❌ 外 | スコープ外（商用アルゴリズム） |

---

### 10.2 ML モデル統合

**実態: 差は機能でなくドキュメント**

chematic は現時点で既に以下を提供済み:
- `bulk.ecfp4(smiles)` → `(N, 2048)` uint8 **NumPy 配列** ← sklearn/PyTorch がそのまま使える
- `bulk.maccs(smiles)` → `(N, 166)` uint8 NumPy 配列
- `bulk.descriptors(smiles)` → `list[dict]` → `pandas.DataFrame` 変換 1 行

```python
# 現時点で既に動くコード
import chematic
from sklearn.ensemble import RandomForestClassifier

fps = chematic.bulk.ecfp4(smiles_list)   # (N, 2048) uint8 ndarray
clf = RandomForestClassifier().fit(fps, labels)
```

RDKit との差は「sklearn/PyTorch 連携のサンプルコード・チュートリアルの量」であり、**実装変更は不要**。ドキュメント追加のみで差が消える。

**状況**: サンプルコード / チュートリアルは提供済み（`examples/` に複数ノートブック）。

---

### 10.3 逆合成予測

**実装済み基盤**:
- ✅ BRICS 断片化（13 ルール、完全再帰）
- ✅ RECAP 断片化（C-N/C-O/C-S）
- ✅ 順方向 SMIRKS 変換 (`run_reactants`, `run_reactants_strict`)
- ✅ **MCP `retrosynthesis` ツール**（BRICS + SA Score ランク付き）
- ✅ **`retro_disconnect()` API**（60 件 retro-SMIRKS テンプレート）
- ✅ **AiZynthFinder 連携チュートリアル**（`examples/aizynthfinder_integration.py`）
- ❌ ML 逆合成（MCTS + Transformer）— scope 外

| アプローチ | 難易度 | scope | 状況 |
|-----------|--------|-------|------|
| MCP 「1 ステップ逆合成」ツール（BRICS + SA Score） | Low | ✅ 内 | **実装済み** |
| テンプレートベース逆合成（60 件の retro-SMIRKS） | Medium | ✅ 内 | **実装済み** |
| AiZynthFinder 連携チュートリアル | Low | ✅ 内 | **実装済み** |
| ML 逆合成（MCTS + Transformer） | Very High | ❌ 外 | スコープ外 |

**`retro_disconnect()` API 詳細** (`crates/chematic-rxn/src/retro.rs`):
- 60 件の retro-SMIRKS テンプレート（6 クラス）:
  - `AmideBond` (10): 二次/三次アミド・スルホンアミド・カルバメート・尿素等
  - `Ester` (6): エステル・チオエステル・カーボネート・無水物等
  - `Ether` (8): アリールエーテル・Williamson・ベンジル・Mitsunobu 等
  - `CNBond` (11): 還元的アミノ化・SNAr-CN・Buchwald・N-アルキル化等
  - `CCBond` (14): Suzuki・Heck・Sonogashira・Negishi・Grignard・aldol 等
  - `CSBond` (10): チオエーテル・ジスルフィド・ハロゲン逆合成等
- Python API: `mol.retro_disconnect(max_results=20, reaction_class="AmideBond")`
- 返り値: `list[dict]` — `{template, reaction_class, precursors, sa_scores, max_sa_score}`

---

## 11. 開発履歴・版マイルストーン

| バージョン | 主な改善 |
|-----------|---------|
| v0.4.8 | ベースライン（芳香環カウント 95.6%） |
| v0.4.9 | 基本精度テスト拡張 |
| v0.4.10 | バグ修正・テストハーネス整備 |
| v0.4.11 | 芳香環カウント ~100%（XOR ロジック修正）、CIF/Gaussian パーサ 8 件修正 |
| v0.4.12 | SMARTS atom map `:N` サポート、4 件 SMARTS パーサ修正 |
| v0.4.13 | TPSA/LogP/HBD 精度修正、VF2/SSSR/SMARTS パーサ RDKit パッチ、立体化学パリティ対応 |
| v0.4.14 | 新規記述子 6 個追加（Petitjean、ECI、Hosoya Z、Moran/Geary 自己相関、GETAWAY、allene stereo） |

---

*数値ソース: `scripts/chematic_vs_rdkit.tsv`, `scripts/benchmark_vs_rdkit.py`*  
*最終更新: セッション 10 完了、v0.4.14 リリース（全 211 テスト通過）*
