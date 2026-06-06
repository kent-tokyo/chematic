# chematic — 全フェーズロードマップ（RDKit 代替を目指す Pure Rust 実装）

目標: C/C++ FFI なし、WASM ネイティブで動く、RDKit の全機能をカバーする単一 Rust クレートエコシステム。

制約: 外部 C/C++ FFI ゼロ、コアクレートは WASM 互換、petgraph 不使用、Python バインディングなし。

---

## Phase 1 — 基盤（完了）

- [x] Cargo ワークスペース構築
- [x] chematic-core: Element（118元素）、Atom、Bond、Molecule、MoleculeBuilder
- [x] chematic-smiles: OpenSMILES パーサー（有機サブセット / ブラケット / 芳香族 / 環 / 分岐 / 立体 / 非連結）
- [x] chematic-smiles: DFS SMILES ライター（環クロージャー番号の正確な割り当て）
- [x] テスト 50件 全パス（aspirin / caffeine / glucose / NaCl 等のラウンドトリップ）
- [x] ワイルドカード原子 `*`: Atom::wildcard() 追加、パーサーで [*] を正しく処理
- [x] 暗黙的 H 数計算: chematic-core/src/valence.rs に implicit_hcount(mol, idx) -> u8
      （結合次数の総和、最小適合価数を選択、電荷による調整）
- [x] ケクレ化: chematic-core/src/kekulization.rs
      （芳香族部分グラフの最大マッチングで二重結合を割り当て）
- [x] 正規 SMILES: chematic-smiles/src/canonical.rs
      （Morgan ランク反復 -> 正規 DFS 順序）

---

## Phase 2 — 分子認識（完了）

新規クレート: chematic-perception、chematic-mol、chematic-depict

- [x] SSSR（最小環集合）: Balducci-Pearlman + GF(2) Gaussian elimination
      find_sssr(mol) -> RingSet  [chematic-perception/src/sssr.rs]
- [x] 芳香族性認識: Hückel 4n+2 π 電子モデル
      assign_aromaticity(mol) -> AromaticityModel
      対応: ベンゼン、ピリジン、ピロール、フラン、ナフタレン  [chematic-perception/src/aromaticity.rs]
- [x] SDF/MOL ファイル形式: V2000 パーサー+ライター、SDF マルチ分子イテレーター
      parse_mol / write_mol / SdfReader  [chematic-mol/]
- [x] SDF V3000 パーサー（拡張ブロック）
      M  V30 BEGIN/END CTAB, ATOM, BOND ブロック対応
      行継続（末尾 `-`）、CHG= / MASS= / HCOUNT= / aamap 対応
      [chematic-mol/src/mol3000.rs]
- [x] 2D 描画エンジン（SVG）: chematic-depict クレート
      - 鎖テンプレート（ジグザグ、BOND_LEN=40px、±30° 交互）
      - 環テンプレートライブラリ（3〜8員環: r = BOND_LEN / (2 sin(π/n))）
      - 縮合環の貪欲配置（重心を基準に外方向を選択）
      - 非 C 原子のラベル表示（白背景 rect + テキスト）
      - ウェッジ/ダッシュ立体結合、二重/三重/芳香族結合の SVG 描画
      [chematic-depict/]
- [x] ChEMBL サンプルセットとのラウンドトリップ検証（大規模テスト）
      ChEMBL から 1000分子以上の SDF を取得し、parse -> write -> parse で一致確認

---

## Phase 3 — 化学インテリジェンス（完了）

新規クレート: chematic-chem、chematic-fp、chematic-smarts

- [x] 分子記述子（chematic-chem）:
      - molecular_weight（平均同位体質量）、exact_mass（モノアイソトピック質量）
      - heavy_atom_count（重原子数）
      - hbd_count（水素結合ドナー: NH, OH）
      - hba_count（水素結合アクセプター: N, O）
      - rotatable_bond_count（非環式単結合 + 非末端重原子間 + アミド除外）
      - tpsa（位相的極性表面積: Ertl 2000 原子タイプ別ルックアップテーブル）
      - logp_crippen（簡略化 Crippen-Wildman 原子寄与テーブル）
      - lipinski_passes（MW<=500, HBD<=5, HBA<=10, LogP<=5）
      - fsp3（sp3 炭素の割合）
      - aromatic_ring_count（芳香族環数：SSSR 経由）
      [chematic-chem/src/descriptors.rs]
- [x] QED スコア（chematic-chem/src/qed.rs）:
      - Bickerton et al. 2012 (Nature Chemistry) の 8 指標幾何平均
      - 7-parameter ADS（Asymmetric Double Sigmoidal）関数 — RDKit と同一パラメータ
      - 113 Brenk 2008 構造アラート SMARTS（第 8 指標）
      - qed(mol) -> f64 in [0, 1]
- [x] Molar Refractivity（chematic-chem/src/descriptors.rs）:
      - Wildman-Crippen 加成モデル（LogP と同一原子タイプフレームワーク）
      - molar_refractivity(mol) -> f64
- [x] 薬物様性フィルター（chematic-chem/src/descriptors.rs）:
      - Veber: TPSA ≤ 140 Å²、回転可能結合数 ≤ 10
      - Egan: TPSA ≤ 131.6 Å²、LogP ≤ 5.88
      - REOS: MW / LogP / HBD / HBA / 電荷 / 重原子数の 6 基準
      - Ghose: MW 160–480, LogP −0.4–5.6, 重原子 20–70, MR 40–130
- [x] 互変異性体ルール拡張（chematic-chem/src/tautomer.rs）:
      - 5 ルール → 15 ルール（チオアミド、チオ-イミノール、チオ-ケト-エノール、6 種クロスヘテロ原子 1,3 プロトン移動）
- [x] BRICS フラグメント化（chematic-chem/src/brics.rs）:
      - Dien et al. 2008 の 16 環境ルールに基づく結合切断
      - brics_bonds(mol) -> Vec<(AtomIdx, AtomIdx)>
      - brics_fragments(mol) -> Vec<Molecule>（[*] アタッチメントポイント付き）
- [x] ECFP / Morgan フィンガープリント（chematic-fp）:
      - 設定可能半径（ECFP4 = r2、ECFP6 = r3）
      - 原子不変量: 原子番号、電荷、次数、H 数、環内フラグ、芳香族フラグ
      - ハッシュ: FNV-1a 64bit（再現性・決定性、to_le_bytes で決定論的）
      - 出力: BitVec2048（[u64; 32]）、fold で 1024/512/256 ビットに畳み込み可
      - Tanimoto 係数、Dice 係数
      [chematic-fp/]
- [x] SMARTS パーサー + VF2 評価器（chematic-smarts）:
      - 原子プリミティブ: [#6] [!C] [a] [A] [r5] [D3] [H2] [R]
      - 結合プリミティブ: ~ @ - = # :
      - 論理演算子: & , ; !（優先順位: NOT > 高優先 AND > OR > 低優先 AND）
      - 再帰 SMARTS `$(...)`: ネスト対応、VF2 アンカー付きマッチング
      - 拡張プリミティブ: [vN] 原子価、[xN] 環結合数、[^N] ハイブリッド化
      - [XN] 全結合数（重原子次数 + 暗黙的 H 数）、[RN] 環帰属数
      - 明示的中性電荷 [+0]/[-0] 対応
      - QueryMolecule 型
      - find_matches(query, mol) -> Vec<HashMap<usize, AtomIdx>>
- [x] 分子標準化 + Murcko スキャフォルド（chematic-chem）:
      - 塩除去: largest_fragment(mol) -> Molecule（最大フラグメント選択）
      - 電荷中和: neutralize_charges(mol) -> Molecule（カルボキシレート/アンモニウム/プロトン化エーテル対応）
      - Murcko スキャフォルド: murcko_scaffold(mol) -> Molecule（固定点リンカー展開）
      - ジェネリック Murcko: generic_murcko_scaffold(mol) -> Molecule
      [chematic-chem/src/standardize.rs, scaffold.rs]
- [x] CIP 立体化学（chematic-chem）:
      - 四面体中心の R/S 割り当て（CIP 優先順位規則）
      - 二重結合の E/Z 割り当て
      - CIPCode 列挙型を Atom に格納

---

## Phase 4 — 類似性・検索・標準化（完了）

- [x] MACCS 166 ビット構造キー（chematic-fp/src/maccs.rs）
      - 標準 MACCS 166-bit SMARTS ベース構造キー
      - `maccs(mol) -> BitVec2048`
      - key 164 = `[!#6;!#1]`（任意ヘテロ原子）などの標準パターン
- [x] 位相的パスフィンガープリント（chematic-fp/src/topo_path.rs）
      - DFS パス列挙、max_len=7（デフォルト）
      - `topo_path(mol, &config) -> BitVec2048`
      - FNV-1a ハッシュ、正規化（前後方向の小さい方を選択）
- [x] AtomPair フィンガープリント（chematic-fp/src/atom_pair.rs）
      - Carhart et al. 1985、原子ペア + BFS 位相距離エンコード
      - `atom_pair_fp(mol) -> BitVec2048`
- [x] Topological Torsion フィンガープリント（chematic-fp/src/atom_pair.rs）
      - Nilakantan et al. 1987、4原子パスエンコード
      - `torsion_fp(mol) -> BitVec2048`
- [x] 最大共通部分構造（MCS）: McGregor または FMCS アルゴリズム
      find_mcs(mols: &[Molecule]) -> QueryMolecule
- [x] 互変異性体正規化（ルールベース）
- [x] 2D ウェッジ結合からの立体認識

---

## Phase 5 — 3D 化学（完了）

新規クレート: chematic-3d（外部依存ゼロ）

- [x] 3D 座標生成（ルールベース DFS 配置）:
      - 理想結合長テーブル（元素ペア + 結合次数別）
      - sp3/sp2/sp 混成別結合角（109.5° / 120° / 180°）
      - 環テンプレートを XY 平面に配置 + 連鎖を分岐として伸長
      - 非連結フラグメントを X 方向にオフセット
      [chematic-3d/src/dg.rs]
- [x] 3D ファイル形式:
      - PDB パーサー/ライター（ATOM/HETATM レコード、距離ベース結合推定）
      - XYZ パーサー/ライター
      [chematic-3d/src/pdb.rs, xyz.rs]
- [x] UFF 力場エネルギー最小化（Pure Rust）:
      - 結合伸縮、角度変形、二面角、VDW、静電相互作用
      - 勾配降下 / LBFGS 最小化

---

## Phase 6 — エコシステムと RDKit パリティ（完了）

新規クレート: chematic-wasm、chematic-rxn、chematic（アンブレラ）

- [x] WASM パッケージ（chematic-wasm）:
      - wasm-bindgen バインディング: パース、ライター、フィンガープリント、記述子計算
      - npm パッケージ: @kent-tokyo/chematic
        - v0.1.3: tpsa/mw/hba/hbd/lipinski/ecfp4
        - v0.1.4: logp/fsp3/qed/exact_mass/rotbonds/aromatic_ring_count/
                  tanimoto_atom_pair/tanimoto_torsion/brics_fragment_count
        - (unreleased): molar_refractivity/formal_charge_sum/veber_passes/
                        egan_passes/reos_passes/ghose_passes
        （"chematic" はnpmが chromatic と類似として拒否 → スコープ付きで公開）
- [x] 2D 描画強化（chematic-depict）:
      - CPK カラーリング（N=青, O=赤, S=黄, Cl=緑, F=黄緑, Br=茶, I=紫, P=橙）
      - render_svg_highlighted / depict_svg_highlighted（黄色ハイライト + 橙ボンド）
- [x] 反応 SMILES / SMIRKS（chematic-rxn）:
      - 反応 SMILES パーサー（>> 区切り）
      - アトムアトムマッピング
      - 反応テンプレート適用
- [x] アンブレラクレート（chematic）:
      - フィーチャーフラグで全サブクレートを再エクスポート
      - feature "full": 全機能有効
      - feature "wasm": 3D と重い依存を無効化
- [x] 検証とベンチマーク:
      - [x] 175 分子データセットで物性値精度を RDKit と定量比較 (docs/rdkit_comparison.md)
            MW: MAE=0.0002 Da, r=1.0000 | HAC: r=1.0000 | HBD: r=0.9974
            TPSA: MAE=0.081 Å², r=0.9999 | HBA: MAE=0.137, r=0.9750
            LogP: MAE=0.134, r=0.9847 (改善: v0.1.0 MAE=1.346 → v0.1.3 MAE=0.298 → Sprint C MAE=0.141 → Sprint D MAE=0.134)
            ECFP4 Tanimoto: Spearman r=0.917 (50×50 ペア)
      - [x] Sprint C — RDKit 品質改善（LogP MAE 0.298→0.141, TPSA MAE 0.759→0.324 Å²）:
            - junction C アトム型修正（縮合芳香族環：naphthalene/indole 等）
            - vinyl C アトム型符号修正（C=C の Crippen 寄与 +0.2274）
            - 硝酸基 N の TPSA 修正（[N+](=O)[O-]: N=41.44 Å², O-=0 Å²）
            - イミン N の TPSA 修正（C=N の非環式 N: 12.89 Å²）
      - [x] Sprint D — RDKit 品質改善 + 不足機能追加（LogP MAE 0.141→0.134, TPSA MAE 0.324→0.081 Å²）:
            - イミン N-H の TPSA 修正（C=N-H: 23.79 Å²、metformin/arginine 誤差解消）
            - リン酸 P の TPSA 修正（P=O あり: 26.88 Å² vs P=O なし: 34.14 Å²）
            - リン酸 P の LogP 修正（P=O あり: +0.7933 vs P=O なし: -0.3451）
            - 5 種リング記述子追加: num_aromatic_heterocycles, num_aliphatic_heterocycles,
              num_saturated_heterocycles, num_spiro_atoms, num_bridgehead_atoms
            - 互変異性体ルール 15 → 20（ルール 16〜20: O→N, O→O, N→C, C→O, C→N）
      - [x] Sprint E — グアニジニウム N LogP 修正 + タウトマー 1,2-shift（LogP MAE 0.134→0.117）:
            - グアニジニウム/アミジン N の LogP 修正（Wildman-Crippen N14 型: -0.335）:
              イミン =N（直接 C=N 二重結合）と隣接グアニジニウム N（C=N 隣の N）を検出
              metformin 誤差 2.07 → ~0.00、arginine 改善
            - 互変異性体 1,2-shift 追加（pyrazole N1H↔N2H 等）:
              find_direct_aromatic_matches + transfer_hydrogen_aromatic（結合次数変更なし）
              enumerate_tautomers: H-assignment フィンガープリントで位置異性体を識別
              canonical_tautomer: 最小 H-assignment で N1H/N2H を同一正規形に収束
      - [x] マルチエージェントセキュリティ/バグ/リファクタリング審査:
            - [Security] 再帰 SMARTS $(...) に深さ上限 8 を追加（SmartsError::RecursionDepthExceeded）
            - [Security] リングクロージャー unwrap() → expect()（不変条件を明文化）
            - [Bug] clone_mol / transfer_hydrogen_aromatic の .ok() → .expect()（サイレントボンド欠落防止）
            - [Refactor] mol_fingerprint の FNV-1a マジック数 → 名前付き定数
            - [Refactor] TPSA 硝酸基検出: 2 回の neighbors スキャン → 1 回の fold
      - [x] criterion による全ホットパスのベンチマーク
      - [x] ChEMBL 37 全量バリデーション: **2,897,819 分子 / 100.000% 成功**（parse + roundtrip）
              curl chembl_37_chemreps.txt.gz | gzip -d | awk | validate_smiles でストリーム検証

---

## 現在のテスト数

| クレート               | テスト数 | 状態     |
|------------------------|---------|---------|
| chematic-core          | 48      | 完了     |
| chematic-smiles        | 57      | 完了     |
| chematic-perception    | 18      | 完了     |
| chematic-mol           | 61      | 完了     |
| chematic-depict        | 39      | 完了     |
| chematic-chem          | 226     | 完了     |
| chematic-fp            | 50      | 完了     |
| chematic-smarts        | 84      | 完了     |
| chematic-3d            | 68      | 完了     |
| chematic-rxn           | 28      | 完了     |
| chematic-wasm          | 162     | 完了     |
| chematic               | 1       | 完了     |
| chematic-iupac         | 8       | 完了     |
| **合計**               | **933** | —        |

---

## 最終クレート構成

    chematic/
    crates/
      chematic-core/        Phase 1  — Atom, Bond, Molecule, Element
      chematic-smiles/      Phase 1  — SMILES パース/ライター/正規化
      chematic-perception/  Phase 2  — SSSR、芳香族性認識
      chematic-mol/         Phase 2  — SDF/MOL V2000+V3000 ファイル形式
      chematic-depict/      Phase 2  — 2D SVG 描画（CPK カラー、ハイライト）
      chematic-chem/        Phase 3  — 分子記述子、BRICS、QED、標準化、CIP 立体化学
      chematic-fp/          Phase 3  — ECFP、MACCS、パス、AtomPair、Torsion FP
      chematic-smarts/      Phase 3  — SMARTS + VF2 部分構造検索（再帰 SMARTS 対応）
      chematic-3d/          Phase 5  — ルールベース 3D 座標、PDB/XYZ 形式
      chematic-rxn/         Phase 6  — 反応 SMILES/SMIRKS
    chematic/               Phase 6  — フィーチャーフラグ付きアンブレラクレート

---

## フェーズ間の依存関係

    Phase 1（コア + SMILES）
      -> Phase 1 ケクレ化
        -> Phase 2（SSSR、芳香族性、SDF、2D 描画）
          -> Phase 3（記述子、ECFP、SMARTS）
            -> Phase 3 SMARTS -> Phase 4（MACCS、MCS、標準化）
              -> Phase 5（3D、力場）
                -> Phase 6（WASM、反応、検証）

---

## Phase 7 — RDKit 完全対等（未着手）

RDKit と比較して未実装の主要機能を優先度順に列挙する。
制約: FFI ゼロ・WASM 互換は変更しない。

### Tier 1 — 高優先度（製薬/化学情報処理ユーザーが最も必要とする機能）

#### 7-1. 反応 SMIRKS 適用（RunReactants）✅ Sprint J 完了
  - [x] RDKit: `rxn.RunReactants(reactants)` → 生成物 SMILES の列挙
  - 実装場所: chematic-rxn/src/transform.rs（実装済み）
  - `run_reactants(smirks, reactants) -> Result<Vec<Vec<Molecule>>, TransformError>`
  - VF2 サブグラフ同形 + BFS 置換基引き継ぎ + カルテシアン積列挙

#### 7-2. トポロジカル記述子 ✅ Sprint G 完了
  - [x] Wiener index（全原子ペア距離の総和）
  - [x] Hall–Kier Kappa 指標 κ1 / κ2 / κ3
  - [x] 分子接続性指標 Chi χ0v / χ1v / χ2v / χ3v / χ4v（Kier–Hall）
  - [x] Bertz 複雑度（BertzCT）
  - [x] Labute 近似表面積（LabuteASA）
  - 実装場所: chematic-chem/src/topo_descriptors.rs（実装済み）
  - 難易度: 中（距離行列が基盤 → 一度実装すれば残りは派生）

#### 7-3. 明示的 H 管理 ✅ Sprint G 完了
  - [x] `add_hydrogens(mol) -> Molecule` — 全暗黙的 H を明示的原子に変換
  - [x] `remove_hydrogens(mol) -> Molecule` — 明示的 H 原子を暗黙的に戻す
  - 現状: implicit_hcount() による暗黙的 H 計算のみ
  - 実装場所: chematic-chem/src/hydrogen.rs（実装済み）
  - 難易度: 低〜中

#### 7-4. SVG グリッド描画 ✅ Sprint G 完了
  - [x] `depict_svg_grid(mols, cols) -> String` — 複数分子を格子状に並べた SVG
  - RDKit: `Draw.MolsToGridImage`
  - 実装場所: chematic-depict/src/grid.rs（実装済み）
  - 難易度: 低（既存 depict_svg を組み合わせるだけ）

---

### Tier 2 — 中優先度（QSAR・3D ワークフロー）

#### 7-5. 形状記述子（3D 座標が必要） ✅ Sprint H 完了
  - [x] 慣性主軸モーメント PMI1 / PMI2 / PMI3
  - [x] 正規化主軸比 NPR1 / NPR2
  - [x] 回転半径（Radius of Gyration）
  - [x] 球面性（Asphericity）・偏心率（Eccentricity）
  - [x] 最良平面比（PBF: Plane of Best Fit）
  - RDKit: `rdMolDescriptors.CalcPMI`, `CalcNPR1/2`, `CalcRadiusOfGyration` 等
  - 実装場所: chematic-3d/src/shape_descriptors.rs（実装済み、3×3 Jacobi eigensolver 手実装）
  - 難易度: 中（固有値分解が必要、nalgebra または手実装）

#### 7-6. コンフォーマー管理 ✅ Sprint I 完了
  - [x] Molecule に複数コンフォーマー（座標セット）を保持する構造
  - [x] `add_conformer()` / `get_conformer()` / `get_conformer_mut()` / `remove_conformer()`
  - [x] コンフォーマー間 RMSD 計算（`conformer_rmsd_no_align` / `conformer_rmsd`）
  - 設計: chematic-core 変更なし。外部コンテナ `ConformerEnsemble` として chematic-3d に実装
  - 実装場所: chematic-3d/src/conformer.rs（実装済み）
  - Kabsch アライメント付き RMSD は既存 jacobi3 を再利用

#### 7-7. UFF パラメータ改善 ✅ Sprint K 完了
  - [x] 元素ペア別理想結合長テーブル（C-C/C-N/C-O/C-S/C-F/C-Cl/C-Br/C-H 等 30+ ペア）
  - [x] 混成軌道判定（SP/SP2/SP3）に基づく理想結合角（O:104.5°, N:107°, S:99° 等）
  - [x] 元素別 UFF/Bondi VDW 半径による VDW 反発エネルギー
  - 実装場所: chematic-3d/src/minimize.rs（改善済み）
  - テスト: 68（旧 58）、+10 新テスト（結合長精度・混成軌道・対称性）
  - MMFF94 フル実装（パラメータテーブル 8 種、原子タイプ 95 種）は工数過大で保留

#### 7-8. 3D からの立体化学割り当て ✅ Sprint H 完了
  - [x] 3D 座標から R/S・E/Z を自動計算（AssignStereochemistryFrom3D）
  - 現状: SMILES の wedge/dash から CIP 割り当てのみ → 3D 符号付き体積＋二面角で独立計算可能
  - 実装場所: chematic-3d/src/stereo3d.rs（新規、1-sphere CIP 優先度で中核を担当）
  - 難易度: 中

---

### Tier 3 — 低優先度（ニッチ・高難度）

#### 7-9. 確率的 3D 埋め込み（ETKDG 相当）
  - [ ] 距離ジオメトリ法（Distance Geometry）による初期座標生成
  - [ ] 実験的ねじれ角分布による改良（ET-DG の "ET" 部分）
  - RDKit: `AllChem.EmbedMolecule` / `EmbedMultipleConfs`
  - 難易度: 非常に高（距離行列の固有値分解 + 実験的ライブラリ必要）

#### 7-10. ハッシュベース FP の密なカウント形式
  - [x] Morgan FP のカウントベクター形式（ビットでなく整数カウント）
  - [x] `GetMorganFingerprint(mol, radius)` → `{hash: count}` 形式
  - 実装場所: chematic-fp/src/ecfp.rs の拡張
  - 難易度: 低（既存 ECFP の出力形式を変えるだけ）

#### 7-11. InChI / InChIKey
  - [ ] 標準 InChI 文字列の生成
  - [ ] InChIKey（27文字ハッシュ）の生成
  - **制約**: IUPAC 公式実装は C ライブラリのみ。Pure Rust では未完成実装のみ存在。
  - FFI ゼロ方針と相反するため、pure Rust 実装が成熟するまで保留
  - 難易度: 非常に高 or FFI 許容が必要

---

### スコープ外（FFI ゼロ方針と相反、または工数が過大）

- ETKDG の完全再現（確率的サンプリング + DG）
- InChI（C ライブラリが唯一の正式実装）
- ML ベース予測モデル（LogP, solubility 等）
- HELM / FASTA 記法（ペプチド/タンパク質）
- 遷移金属・錯体化合物への対応（配位化学）

---

## 実装推奨順序（Sprint G〜）

```
Sprint G: ✅ 7-2（トポロジカル記述子）+ 7-3（明示的 H 管理）+ 7-4（SVG グリッド）
          → コード追加のみ、破壊的変更なし、テスト +38（582→620）
Sprint H: ✅ 7-5（形状記述子）+ 7-8（3D から立体化学）
          → chematic-3d/src/shape_descriptors.rs + stereo3d.rs、テスト +15（620→635）
Sprint I: ✅ 7-6（コンフォーマー管理）
          → chematic-3d/src/conformer.rs、Kabsch RMSD、テスト +14（623→637）
Sprint J: ✅ 7-1（RunReactants）
          → chematic-rxn/src/transform.rs、VF2 + BFS 置換基引き継ぎ、テスト +11（612→623）
Sprint K: ✅ 7-7（UFF パラメータ改善）
          → 元素別結合長・混成軌道角・VDW 半径、テスト +10（637→646）
Sprint L: ✅ Sprint L audit — セキュリティ/バグ/リファクタリング審査（0.1.5 → 0.1.6）
Sprint M: ✅ SMARTS ハイライト表示 + クリックハイライト + 反応スキーム（demo 0.1.11）
Sprint N: ✅ タブ UI + 3D インタラクティブビューア（demo 0.1.12）
Sprint P: ✅ SDF/MOL WASM バインディング + EState インデックス + パスフィンガープリント WASM（v0.1.14）
Sprint Q: ✅ IFG + SA Score + Gasteiger 電荷 + VSA 記述子 + MaxMin/Butina（v0.1.15）
          → テスト: 697 → 736（+39）
Sprint R: ✅ E/Z 二重結合立体化学 SMILES 出力（v0.1.16）
Sprint S: ✅ SA スコア フラグメントテーブル実装（v0.1.17）
          → テスト: 742 → 743（+1）
Sprint T: ✅ per-atom カラーハイライト + 名前付き官能基検出 + 原子情報 API（demo v0.1.18）
Sprint U: ✅ インタラクティブ記事向け WASM 利便性 API（v0.1.19）

## Phase 8 — WASM 機能拡充・ファイル形式・編集 API（v0.1.20〜v0.1.21）

Sprint V–AA: ✅ WASM エクスポート 84 → 103 に拡張（v0.1.20）
  - Murcko / 互変異性体 / 標準化 / MACCS / 一括記述子 / MOL 2D座標修正
  - PAINS/CIP 詳細 / ECFP6 / Dice / 3D 形状記述子 / MaxMin・Butina / MCS
  - V3000 読み込み / 3D 最小化 / SDF プロパティ読み書き / SMARTS ハイライトグリッド
  - XYZ/PDB I/O / per-atom 記述子 / SSSR / カスタム ECFP / 立体異性体列挙
  - BRICS SMILES / AtomPair・Torsion bitvec / FCFP6 / SDF 書き込み
  - FCFP4/6 bitvec / Dice ECFP6 / write_smiles / 反応 SMILES 正規化
  - ConformerEnsemble WASM / R-group 分解 / MMP 分析
  - CML read/write / CDXML read / Mutable API / DepictData / SDF・V3000 write / CPK
  - テスト: 743 → 863（+120）

Sprint v0.1.21: ✅ Mutable API 拡張・SDF/CDXML 機能強化（v0.1.21）
  - chematic-core: with_atom_charge, with_atom_element, with_bond_added → (Mol, BondIdx)
  - chematic-mol: parse_mol_with_coords, parse_sdf_with_coords, parse_cdxml_all, CDXML 立体化学
  - chematic-depict: depict_data_with_coords
  - WASM: mol_with_atom_charge, mol_with_atom_element, cdxml_to_smiles_json, mol_block_coords_json, depict_data_with_coords_json
  - テスト: 863 → 869（+6）

Sprint v0.1.22: ✅ MCS ring-awareness constraints（Issue #1）
  - ring_matches_ring_only: McGregor 探索フェーズで SSSR を使いリング↔非リングのクロスマッチをブロック
  - complete_rings_only: 探索後の反復後処理で mol[0] の部分リングを除去
  - テスト: 869 → 877（+8）

Sprint v0.1.23: ✅ Element 半径 API・implicit H 補完・芳香族性適用 API（v0.1.23）
  - chematic-core: Element::vdw_radius() / covalent_radius()（Bondi 1964 + Alvarez 2013/2008 テーブル 118 元素、不明は 1.70/0.77 フォールバック）
  - chematic-core: Molecule::implicit_hydrogen_count(idx)（valence::implicit_hcount のラッパー）
  - chematic-core: Molecule::total_formula()（暗黙的 H を含む Hill 式 — CH4, C2H6O 等）
  - chematic-core: Molecule::with_atom_aromatic() / with_bond_order()（immutable update API 拡張）
  - chematic-perception: apply_aromaticity(mol) → ケクレ化分子に芳香族フラグと BondOrder::Aromatic を適用した新 Molecule を返す
  - chematic-3d: minimize_uff() エイリアス（既存 minimize() の UFF 最小化を名前で発見しやすくする）
  - テスト: 877 → 886（+9）

Sprint v0.1.24: ✅ validate_valence 公開 API + run_reactants 生成物フィルタリング（v0.1.24）
  - chematic-core: ValenceError 構造体 + validate_valence(mol) -> Vec<ValenceError>（元素別 normal_valences + 形式電荷調整）
  - chematic-core/lib.rs: ValenceError / validate_valence を re-export
  - chematic-perception/lib.rs: chematic_core から re-export（chematic::perception::validate_valence で参照可能）
  - chematic-rxn: run_reactants で生成物に validate_valence を適用し過原子価の product set を除外
  - テスト: 886 → 893（+7）

Sprint v0.1.25: ✅ suggest_bond_direction 公開 API（v0.1.25）
  - chematic-depict/layout.rs: suggest_bond_direction(mol, atom, layout) -> f64（ラジアン）
    - 既存結合角度を収集 → 30° グリッド + 化学的オフセット（sp2±120°、sp3ジグザグ±150°）の候補から最大最小分離角を選択
  - chematic-depict/lib.rs: BOND_LEN・suggest_bond_direction を re-export
  - draw 側の 30° 総当たり独自実装の置き換えが可能に
  - テスト: 893 → 897（+4）

Sprint v0.1.26: ✅ atom_color_rgb 公開 API（v0.1.26）
  - chematic-depict/svg.rs: atom_color_rgb(atomic_number: u8) -> [u8; 3]（atom_color と同一 CPK 値、hex 解析なし）
  - chematic-depict/lib.rs: re-export
  - draw 側の hex パーサー独自実装の置き換えが可能に（egui::Color32::from_rgb 直接利用）
  - テスト: 897 → 900（+3）

Sprint v0.1.27: ✅ MolMetadata builder API（v0.1.27）
  - chematic-mol/mol2000.rs: MolMetadata::with_name(name) -> Self、with_comment(comment) -> Self
  - SDF エクスポート時に MolMetadata::default().with_name("...").with_comment("...") で名前・コメントを設定可能に
  - テスト: 900 → 902（+2）

Sprint v0.1.27-ext: ✅ E/Z 二重結合立体化学（2D 座標から）+ 拡張 StereoGroup + 同位体分布（v0.1.27）
  - chematic-perception: assign_ez_from_2d(mol, coords) — 2D 座標の外積から E/Z を割り当て
                          cip_ez_descriptor(mol, bond_idx, coords) -> Option<CipCode> — 特定結合の E/Z 返却
                          [crates/chematic-perception/src/stereo2d.rs; lib.rs に再エクスポート]
                          アルゴリズム: 二重結合ベクトル vs 置換基位置ベクトルの 2D 外積 + 1-sphere CIP 優先度
  - chematic-core: StereoGroupKind 列挙型（Absolute / Or(u32) / And(u32)）、StereoGroup 構造体
                   Molecule.stereo_groups フィールド + stereo_groups() / set_stereo_groups() / add_stereo_group() メソッド
                   MoleculeBuilder に add_stereo_group() メソッド + from_molecule() が stereo_groups をコピー
                   [crates/chematic-core/src/stereo_group.rs; lib.rs に再エクスポート]
  - chematic-mol: V3000 パーサーが BEGIN COLLECTION / MDLV30/STEABS / MDLV30/STEOR<n> / MDLV30/STEAND<n> を解析
                  V3000 ライターが stereo_groups 存在時に COLLECTION ブロックを出力
                  ラウンドトリップテスト追加
                  [crates/chematic-mol/src/mol3000.rs]
  - chematic-chem: isotope_distribution(mol, resolution) -> Vec<(f64, f64)>
                   (m/z, 相対強度) ペアを返す、基準ピーク=1.0 で正規化
                   resolution パラメータで指定 Da 以内のピークをマージ
                   H/C/N/O/F/Si/P/S/Cl/Br/I/Se/Na/K/As の同位体対応
                   明示的同位体ラベル（atom.isotope）を優先使用
                   [crates/chematic-chem/src/isotope_distribution.rs; lib.rs に再エクスポート]
  - chematic（アンブレラクレート）: lib.rs の //! モジュールドキュメント全面改訂（機能表・クイックスタート例・フィーチャーフラグ表）
                   Cargo.toml: description 更新、categories に parser-implementations/rendering 追加
                   [package.metadata.docs.rs] セクション追加: features=["full"], rustdoc-args=["--cfg","docsrs"]

Sprint v0.1.28: ✅ 全残タスク実装（v0.1.28）
  - Issue C: BricsConfig { min_fragment_size } + brics_fragments_with_config（chematic-chem）
  - Issue E: MatchConfig { max_matches } + find_matches_with_config（chematic-smarts）
  - Issue A: AtomCompare / BondCompare enum + McsConfig フィールド追加（chematic-smarts）
            AnyHeavyAtom モードで異種ヘテロ環間の scaffold hopping MCS が動作
  - ⑩ xlogp3.rs: xlogp3() / xlogp3_per_atom() — Cheng 2007 原子型貢献テーブル（chematic-chem）
  - ⑪ chematic-iupac: 新クレート（Pure Rust、ネットワーク不要）
            直鎖アルカン/アルケン/アルキン、シクロアルカン、アルコール/アミン/ハロアルカン
            IupacError::NotSupported で未対応構造を明示
  - テスト: 902 → 915（+13）

Sprint v0.1.29: ✅ Mutable Molecule API + Fragments + MoleculeBuilder::from_molecule（v0.1.29）
  - Molecule::add_atom / remove_atom / add_bond / remove_bond / set_charge / set_element / set_cip_code
  - Molecule::is_connected() / fragments() → 連結成分分割
  - MoleculeBuilder::from_molecule(mol)
  - テスト: 915 → 924（+9）

Sprint v0.1.30: ✅ 2D 立体化学 + Aromatize/Kekulize in-place（v0.1.30）
  - chematic-perception: stereo2d.rs 新規（assign_stereo_from_2d / apply_stereo_from_2d）
  - chematic-perception: aromatize(mol: &mut Molecule) / kekulize_inplace(mol: &mut Molecule)
  - テスト: 924 → 926（+2）

Sprint v0.1.31: ✅ SdfRecord 拡張 + 反応 SVG + 化学略号（v0.1.31）
  - SdfRecord: coords: Vec<(f64,f64)> + meta: MolMetadata + properties: HashMap<String,String> 追加
  - chematic-depict: depict_reaction_svg / depict_reaction_svg_opts（反応物→矢印→生成物 SVG）
  - chematic-chem: expand_abbreviation / abbreviations（30 略号テーブル）
  - テスト: 926 → 929（+3）

Sprint v0.1.32: ✅ MDL RXN ファイル + formula_with_isotopes（v0.1.32）
  - chematic-mol: parse_rxn_file / write_rxn_file（MDL RXN V2000 フォーマット）
  - chematic-core: Molecule::formula_with_isotopes()（²H・¹³C 等の同位体ラベル付き分子式）
  - テスト: 929 → 933（+4）

Sprint v0.1.25: ✅ P2 機能完成・リリース（2026-06-06）
  - detect_crossings: 2D レイアウト品質評価（結合交差検出）
  - invert_stereocenter: R/S キラリティ反転（ウェッジ結合反転）
  - enumerate_stereoisomers: 立体異性体全列挙（2^n、最大 64）
  - render_svg_with_metadata: SVG メタデータ埋め込み（SMILES）
  - find_reaction_center: 反応中心分析（broken/formed bonds + changed atoms）
  - テスト: 865 → 935（+70）
  - cargo & npm publish 完了
  - CHANGELOG / README 全言語更新済み

Sprint v0.1.26: 🔨 Issue D + P3 Features（実装中）

### Completed ✅
  - [x] Issue D (matchChiralTag): `McsConfig.match_chiral_tag` 実装済み
        - R/S 鏡像体マッチング制御（default: false）
        - 3 つの新規テスト追加（enantiomer blocking/allowing）
        - 実装場所: crates/chematic-smarts/src/mcs.rs
        - テスト: 87 tests all passing
  
  - [x] parse_condensed(): "CH3COOH" → structure parsing 実装済み
        - condensed formula 字句解析 + 官能基置換
        - 実装場所: crates/chematic-chem/src/condensed.rs（新規）
        - テスト: 10 件実装（基本的なケースカバー）
        - Note: H-count digits (CH3) 処理は将来の改善対象
        - 実装: parse_condensed(input) → Result<Molecule, CondensedError>

### In Progress 🔨
  - [ ] WASM bindings
        - find_reaction_center_json(reaction_smiles) → JSON
        - standardize_smiles(mol, opts) → SMILES
  - [ ] Demo updates
        - "Stereo" タブ追加（立体異性体列挙）
        - "Reaction" タブ拡張（broken/formed bonds ハイライト）

### Deferred (High Difficulty):
  - [ ] PolymerSMILES (SRU notation): [CH2CH2]n support
        - 推奨: v0.2.0 リメジャー実装

## テスト現況
- chematic-smarts: 87 tests (v0.1.25 84 + matchChiralTag 3)
- chematic-chem: 248 tests (v0.1.25 238 + parse_condensed 10)
- 全体: 948 tests passing (v0.1.25 935 + 13)

## 次のステップ
- Step 5: WASM binding 追加完了
- Step 6: Demo HTML 更新完了
- Step 7: Version 0.1.26 bump + CHANGELOG + 全テスト確認
```

---

## Issue 候補（Issue #1 類似パターン — 今後発生しうる問題）

Issue #1 のパターン: **アルゴリズムが位相的には正しい結果を返すが、化学的に無意味な結果になる**（RDKit にある制約オプションが chematic に未実装で、移行時にサイレントに誤った結果が生成される）。

以下は同パターンで将来 Issue 化する可能性が高い項目。

### ✅ Issue 候補 A (🔴 高): MCS — `atomCompare` / `bondCompare` レベル（Sprint v0.1.28 で解決）
  - **状態**: ✅ Sprint v0.1.28 で実装済み（`AtomCompare::Elements/AnyHeavyAtom/Any`, `BondCompare::OrderOrAromatic/Any`）
  - **実装済み場所**: `crates/chematic-smarts/src/mcs.rs` (McsConfig struct lines 48-69)
  - `find_mcs_with_config` で `McsConfig { atom_compare, bond_compare, ... }` を使用
  - キラリティ比較は別 Issue D を参照

### ✅ Issue 候補 B (🔴 高): `run_reactants` — 生成物の原子価バリデーションなし（Sprint v0.1.24 で解決）
  - **現状**: SMIRKS 適用後の生成物 Molecule に valence チェックなし
  - **症状**: 四級窒素へのアルキル化で `[N](C)(C)(C)(C)` (valence 5) が無音で生成される
  - **RDKit 対応**: デフォルトで `sanitizeMols=True`（原子価違反で生成物を除外）
  - **対象**: `crates/chematic-rxn/src/transform.rs`
  - **実装**: 生成物ごとに `valence::bond_order_sum > max_valence` を検査し除外（or `TransformError::InvalidProduct`）
  - **推奨 Sprint**: v0.1.24（正確性問題のため優先）

### ✅ Issue 候補 C (🟡 中): BRICS — `minFragmentSize` オプションなし（Sprint v0.1.28 で解決）
  - **状態**: ✅ Sprint v0.1.28 で実装済み（`BricsConfig { min_fragment_size }` + `brics_fragments_with_config`）
  - **実装済み場所**: `crates/chematic-chem/src/brics.rs` (BricsConfig lines 69-76, brics_fragments_with_config)
  - min_fragment_size で 1-2 原子の無意味なフラグメントを除外可能

### 🔨 Issue 候補 D (🟡 中): MCS — `matchChiralTag` オプションなし（Sprint v0.1.26 で実装中）
  - **状態**: 実装予定（v0.1.26）— キラル SAR 解析向けの重要な機能
  - **症状**: R/S 鏡像体間の MCS が「全原子一致」になる（化学的には別化合物）
  - **RDKit 対応**: `matchChiralTag=True`
  - **対象**: `crates/chematic-smarts/src/mcs.rs`
  - **実装**: `McsConfig { match_chiral_tag: bool }` (default: false) を追加、`atoms_compatible` で chirality チェック
  - **テスト**: R-Ala vs S-Ala で動作確認

### ✅ Issue 候補 E (🟡 中): `find_matches` — マッチ数上限なし（Sprint v0.1.28 で解決）
  - **状態**: ✅ Sprint v0.1.28 で実装済み（`MatchConfig { max_matches }` + `find_matches_with_config`）
  - **実装済み場所**: `crates/chematic-smarts/src/match_vf2.rs` (lines 73-78, match_recursive lines 102-104)
  - max_matches でメモリ爆発を防止可能

### Issue 候補 F (🟡 中): VF2 部分構造検索 — キラリティ考慮なし
  - **状態**: 未着手
  - **症状**: `[C@@H](F)(Cl)Br` の SMARTS でキラルな部分構造を検索すると、R 体・S 体の両方がヒットする
            （`matchChiralTag` 相当の VF2 サイドの制約が未実装）
  - **RDKit 対応**: `RWMol.HasSubstructMatch(query, useChirality=True)`
  - **対象**: `crates/chematic-smarts/src/match_vf2.rs`
  - **実装方針**: `atoms_compatible()` で SMARTS 側の chirality `@`/`@@` フラグと分子側の `CipCode::R/S` を照合
  - **優先度**: 中（MCS matchChiralTag の VF2 拡張として自然な対）

### Issue 候補 G (🟢 低): ECFP フィンガープリント — キラリティ不変量なし
  - **状態**: 未着手
  - **症状**: R/S 鏡像体（例: R-Ala vs S-Ala）の ECFP4 Tanimoto が 1.0 になる
            （chirality は原子不変量に含まれておらず、鏡像体を区別しない）
  - **RDKit 対応**: `GetMorganFingerprintAsBitVect(mol, radius, useChirality=True)`
  - **対象**: `crates/chematic-fp/src/ecfp.rs`
  - **実装方針**: `EcfpConfig { use_chirality: bool }` を追加し、初期不変量の計算で `atom.cip_code` を組み込む
  - **優先度**: 低（ECFP は通常キラリティなしで利用される。StereoGroup 対応後に検討）
