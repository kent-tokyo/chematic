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
      - [x] criterion による全ホットパスのベンチマーク
      - [x] ChEMBL 37 全量バリデーション: **2,897,819 分子 / 100.000% 成功**（parse + roundtrip）
              curl chembl_37_chemreps.txt.gz | gzip -d | awk | validate_smiles でストリーム検証

---

## 現在のテスト数

| クレート               | テスト数 | 状態     |
|------------------------|---------|---------|
| chematic-core          | 30      | 完了     |
| chematic-smiles        | 52      | 完了     |
| chematic-perception    | 14      | 完了     |
| chematic-mol           | 37      | 完了     |
| chematic-depict        | 15      | 完了     |
| chematic-chem          | 216     | 完了     |
| chematic-fp            | 44      | 完了     |
| chematic-smarts        | 75      | 完了     |
| chematic-3d            | 25      | 完了     |
| chematic-rxn           | 15      | 完了     |
| chematic-wasm          | 18      | 完了     |
| chematic               | 1       | 完了     |
| **合計**               | **542** | —        |

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
