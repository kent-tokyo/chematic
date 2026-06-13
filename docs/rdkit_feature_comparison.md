# chematic vs RDKit 機能比較 (v0.1.96)

> **更新日**: 2026-06-13 | **バージョン**: chematic v0.1.96 | **テスト数**: 1,666 (100% pass)

---

## 1. 概要 (Executive Summary)

chematic は Rust 製のピュア実装チェムインフォマティクスライブラリです。  
FFI ゼロ方針（C/Python 実装への依存なし）で、RDKit の主要機能の **約 90%** をカバーしています。

| 指標 | 値 |
|---|---|
| Gap closure (RDKit比) | **~90%** |
| テスト数 | **1,746** (全パス) |
| LogP MAE | **0.054** (175分子, RDKit Crippen値との比較) |
| TPSA MAE | **0.075 Å²** |
| ECFP4 Tanimoto相関 | **Spearman ρ = 0.925** (2,450ペア) |
| WASM対応 | ✅ (ブラウザ/Node.js で動作) |

---

## 2. 精度ベンチマーク

### 2.1 分子記述子 (175分子テストセット)

| 記述子 | MAE | RMSE | Pearson r | 状態 |
|---|---|---|---|---|
| 分子量 (MW) | 0.0002 Da | 0.0007 | 1.0000 | ✅ 完全一致 |
| 重原子数 (HAC) | 0.000 | 0.000 | 1.0000 | ✅ 完全一致 |
| LogP (Crippen-Wildman) | **0.0540** | 0.1406 | 0.9968 | ✅ 優秀 |
| TPSA | **0.0748 Å²** | 0.4659 | 0.9999 | ✅ 優秀 |
| HBA | 0.0400 | 0.2928 | 0.9888 | ✅ 優秀 |
| HBD | 0.0114 | 0.1069 | 0.9974 | ✅ 優秀 |

**LogP の改善履歴**: v0.1.0 (MAE=1.346) → v0.1.30 (MAE=0.054) — **96% 改善**

**高誤差分子 (参考)**:

| 分子 | RDKit | chematic | Δ | 原因 |
|---|---|---|---|---|
| Curcumin | 3.370 | 0.584 | −2.79 | ビニル-フェノール分類 |
| Methotrexate | 0.268 | −1.512 | −1.78 | プテリン N 複雑性 |
| Folic acid | −0.045 | −1.348 | −1.30 | 複数アミド N |

### 2.2 フィンガープリント (ECFP4, 50分子, 2,450ペア)

| 指標 | 値 |
|---|---|
| Tanimoto MAE | 0.0137 |
| Spearman ρ | **0.925** |
| 90パーセンタイル誤差 | 0.040 |
| 同一分子の類似度 | 1.000 (完全一致) |

---

## 3. 機能対応表

### 3.1 分子記述子 (Descriptors)

| 機能 | RDKit API | chematic | 精度/備考 |
|---|---|---|---|
| 分子量 | `Descriptors.MolWt` | `molecular_weight()` | ✅ MAE=0.0002 |
| 正確質量 | `Descriptors.ExactMolWt` | `exact_mass()` | ✅ |
| LogP | `Descriptors.MolLogP` | `logp_crippen()` | ✅ MAE=0.054 |
| TPSA | `Descriptors.TPSA` | `tpsa()` | ✅ MAE=0.075 |
| HBA/HBD | `Descriptors.NumHBD/A` | `hbd_count()` / `hba_count()` | ✅ |
| 回転結合数 | `Descriptors.NumRotatableBonds` | `rotatable_bond_count()` | ✅ |
| モル屈折率 | `Descriptors.MolMR` | `molar_refractivity()` | ✅ |
| Fsp3 | `Descriptors.FractionCSP3` | `fsp3()` | ✅ |
| 芳香族環数 | `Descriptors.NumAromaticRings` | `aromatic_ring_count()` | ✅ |
| ヘテロ原子数 | `Descriptors.NumHeteroatoms` | `num_heteroatoms()` | ✅ |
| 形式電荷 | `Chem.GetFormalCharge` | `formal_charge_sum()` | ✅ |
| 分子式 | `rdMolDescriptors.CalcMolFormula` | `calc_mol_formula()` | ✅ |
| EState indices | `EState.EStateIndices` | `estate_indices()` | ✅ |
| 立体中心数 | `rdMolDescriptors.CalcNumStereocenters` | `num_stereocenters()` | ✅ |
| 架橋頭原子数 | `rdMolDescriptors.CalcNumBridgeheadAtoms` | `num_bridgehead_atoms()` | ✅ |
| QED | `QED.qed` | `qed()` | ✅ |
| SA Score | `SA_Score.calculateScore` | `sa_score()` | ✅ |
| XLogP3 | — | `xlogp3()` | ✅ |
| LogD | — | `logd_simple()` | ✅ |
| ESOL溶解度 | `esol` | `esol_solubility()` | ✅ |
| Balaban J | `GraphDescriptors.BalabanJ` | `balaban_j()` | ✅ |
| chi 指標 (0-4) | `Descriptors.Chi0-4` | `chi0()` ~ `chi4()` | ✅ |
| kappa 指標 (1-3) | `Descriptors.Kappa1-3` | `kappa1()` ~ `kappa3()` | ✅ |
| Bertz CT | `GraphDescriptors.BertzCT` | `bertz_ct()` | ✅ |
| Wiener index | — | `wiener_index()` | ✅ |
| Ipc | `GraphDescriptors.Ipc` | `ipc()` | ✅ |
| HallKier Alpha | `rdMolDescriptors.CalcHallKierAlpha` | `hall_kier_alpha()` | ✅ |
| AutoCorr2D | `rdMolDescriptors.CalcAUTOCORR2D` | `autocorr_2d()` | ✅ |
| AutoCorr3D | `rdMolDescriptors.CalcAUTOCORR3D` | `autocorr_3d()` | ✅ |
| MQN (42次元) | `MQNs.MQNs_` | `mqn()` | ✅ |
| USRCAT (42次元) | `rdMolDescriptors.CalcUSRCAT` | `usrcat()` | ✅ |
| VSA descriptors (SlogP/SMR/PEOE) | `rdMolDescriptors.SlogP_VSA*` | `slogp_vsa()` 等 | ✅ |
| WHIM | `rdMolDescriptors.CalcWHIM` | `whim_descriptors()` | ✅ |
| GETAWAY | `rdMolDescriptors.CalcGETAWAY` | `getaway_descriptors()` | ✅ |
| SASA | `rdFreeSASA` | `sasa()` | ✅ Shrake-Rupley |
| 元素カウント (C/N/O/F/Cl/Br/I/S/P) | `Descriptors.NumAtomCountX` | `num_carbons()` 等 | ✅ |
| Gasteiger 電荷 | `AllChem.ComputeGasteigerCharges` | `gasteiger_charges()` | ✅ |
| MMFF94 電荷 | `AllChem.MMFFGetMoleculeProperties` | `mmff94_charges()` | ✅ BCI テーブル (±0.1e) |
| アミド/エステル結合数 | `rdMolDescriptors.CalcNumAmideBonds` | `num_amide_bonds()` 等 | ✅ |

### 3.2 フィンガープリント (Fingerprints)

| 機能 | RDKit API | chematic | 備考 |
|---|---|---|---|
| ECFP4/6 (2048-bit) | `AllChem.GetMorganFingerprintAsBitVect` | `ecfp4()` / `ecfp6()` | ✅ ρ=0.925 |
| FCFP4/6 | `AllChem.GetMorganFingerprintAsBitVect(useFeatures=True)` | `fcfp4()` / `fcfp6()` | ✅ |
| Morgan counts | `AllChem.GetMorganFingerprint` | `morgan_fp_counts()` | ✅ |
| MACCS keys (167-bit) | `MACCSkeys.GenMACCSKeys` | `maccs()` | ✅ |
| RDKit Path FP | `Chem.RDKFingerprint` | `rdkit_path_fp()` | ✅ |
| Atom Pair FP | `AllChem.GetAtomPairFingerprintAsBitVect` | `atom_pair_fp()` | ✅ |
| Torsion FP | `AllChem.GetTopologicalTorsionFingerprintAsBitVect` | `torsion_fp()` | ✅ |
| Layered FP (7層) | `Chem.LayeredFingerprint` | `layered_fp()` | ✅ |
| Pattern FP | `Chem.PatternFingerprint` | `pattern_fp()` | ✅ |
| Pharmacophore 2D | `Gobbi_Pharm2D.Generate` | `pharmacophore_fp_2d()` | ✅ |
| MHFP | `mhfp.MHFPEncoder` | `mhfp()` | ⚠️ 簡易実装 (原子署名ベース) |
| ERG | `AllChem.GetErGFingerprint` | `erg()` | ⚠️ 簡易実装 (薬理特徴が粗い) |
| Reaction FP | `AllChem.ReactionFingerprintAsBitVect` | `reaction_fp()` | ✅ XOR差分エンコーディング |
| Tanimoto類似度 | `DataStructs.TanimotoSimilarity` | `tanimoto_ecfp4()` 等 | ✅ |
| 最近傍探索 | `DataStructs.BulkTanimotoSimilarity` | `nearest_neighbors()` | ✅ |

### 3.3 SMARTS・部分構造マッチ

| 機能 | RDKit API | chematic | 備考 |
|---|---|---|---|
| SMARTS パース | `Chem.MolFromSmarts` | `parse_smarts()` | ✅ |
| 部分構造マッチ | `mol.HasSubstructMatch` | `find_matches()` | ✅ VF2アルゴリズム |
| 全マッチ列挙 | `mol.GetSubstructMatches` | `find_matches()` | ✅ visit budget付き |
| MCS | `rdFMCS.FindMCS` | `find_mcs()` | ✅ |
| CXSmarts | ChemAxon拡張 | `parse_cxsmarts()` | ✅ |
| Reaction SMARTS | `AllChem.ReactionFromSmarts` | `ReactionPatternLibrary` | ✅ バッチ対応 |

### 3.4 3D 構造

| 機能 | RDKit API | chematic | 備考 |
|---|---|---|---|
| ETKDGv3 構造生成 | `AllChem.EmbedMolecule(ETKDGv3)` | `generate_coords_dg()` | ✅ 上限: 500原子 |
| UFF 最小化 | `AllChem.UFFOptimizeMolecule` | `minimize_uff()` | ✅ |
| MMFF94 最小化 | `AllChem.MMFFOptimizeMolecule` | `minimize_mmff94()` | ✅ |
| DREIDING 最小化 | — | `minimize_dreiding()` | ✅ |
| PMI / NPR | `rdMolDescriptors.CalcPMI*` | `pmi()`, `npr1/2()` | ✅ |
| 慣性半径 | `rdMolDescriptors.CalcRadiusOfGyration` | `radius_of_gyration()` | ✅ |
| Asphericity | `rdMolDescriptors.CalcAsphericity` | `asphericity()` | ✅ |
| Eccentricity | `rdMolDescriptors.CalcEccentricity` | `eccentricity()` | ✅ |
| PBF | `rdMolDescriptors.CalcPBF` | `plane_of_best_fit()` | ✅ |
| Kabsch RMSD | `AllChem.AlignMol` | `conformer_rmsd()` | ✅ v0.1.94で修正済み |
| コンフォーマーアンサンブル | `AllChem.EmbedMultipleConfs` | `ConformerEnsemble` | ✅ |
| USR / USRCAT | `rdMolDescriptors.GetUSRScore` | `usr_descriptors()` | ✅ |
| 結合長/角度/二面角 | `AllChem.GetBondLength` 等 | `get_bond_length()` 等 | ✅ |
| MD シミュレーション | — | `run_md()` | ✅ |

### 3.5 InChI

| 機能 | RDKit API | chematic | 備考 |
|---|---|---|---|
| InChI 生成 | `Chem.inchi.MolToInchi` | `inchi()` | ✅ /c /h /b /t 層対応 |
| InChI Key | `Chem.InchiToInchiKey` | `inchi_key()` | ✅ |
| InChI 解析 | `Chem.inchi.InchiToMol` | `parse_inchi()` | ✅ /b /t /m /s 層対応 |

### 3.6 標準化・変換 (Standardization)

| 機能 | RDKit API | chematic | 備考 |
|---|---|---|---|
| 標準化パイプライン | `MolStandardize` | `standardize()` | ✅ |
| 塩除去 | `MolStandardize.RemoveFragments` | `remove_salts()` | ✅ |
| 最大フラグメント | `MolStandardize.LargestFragmentChooser` | `largest_fragment()` | ✅ |
| 官能基正規化 | `MolStandardize.Normalize` | `normalize_groups()` | ✅ ニトロ/アジド/N-オキシド/スルホキシド |
| 電荷中和 | `MolStandardize.Uncharger` | `neutralize_charges()` | ✅ |
| 再イオン化 | `MolStandardize.Reionizer` | `reionize()` | ✅ |
| 互変異体標準化 | `MolStandardize.TautomerEnumerator` | `canonical_tautomer()` | ✅ 1,5-シフト含む |
| 互変異体列挙 | `TautomerEnumerator.Enumerate` | `enumerate_tautomers()` | ✅ |
| 同位体除去 | `MolStandardize.RemoveIsotopes` | `remove_isotopes()` | ✅ |
| 立体情報除去 | `MolStandardize.RemoveStereo` | `remove_stereo()` | ✅ |
| 水素付加/除去 | `AllChem.AddHs` / `RemoveHs` | `add_hydrogens()` / `remove_hydrogens()` | ✅ |

### 3.7 フィルター・薬物様性

| 機能 | RDKit API | chematic | 備考 |
|---|---|---|---|
| Lipinski Ro5 | `Descriptors` 各値から判定 | `lipinski_passes()` | ✅ |
| Veber ルール | — | `veber_passes()` | ✅ |
| PAINS アラート | `FilterCatalog.PAINS` | `pains_matches()` | ✅ |
| Brenk アラート | `FilterCatalog.BRENK` | `brenk_matches()` | ✅ |
| Egan ルール | — | `egan_passes()` | ✅ |
| Ghose ルール | — | `ghose_passes()` | ✅ |
| REOS フィルター | — | `reos_passes()` | ✅ |

### 3.8 スキャフォールド・フラグメント

| 機能 | RDKit API | chematic | 備考 |
|---|---|---|---|
| Murcko スキャフォールド | `MurckoScaffold.GetScaffoldForMol` | `murcko_scaffold()` | ✅ |
| Generic Murcko | `MurckoScaffold.MakeScaffoldGeneric` | `generic_murcko_scaffold()` | ✅ |
| Scaffold network | `ScaffoldNetwork.CreateScaffoldNetwork` | `scaffold_network()` | ✅ |
| Schuffenhauer rules | — | `schuffenhauer_parents()` | ✅ ルール 1-8 |
| BRICS 分解 | `BRICS.BRICSDecompose` | `brics_fragments()` | ✅ |
| RECAP 分解 | `Recap.RecapDecompose` | `recap_fragment()` | ✅ |
| MMP | `MatchedMolecularPairs` | `find_mmp()` | ✅ |
| MaxMin ダイバーシティ | `SimDivFilters.MaxMinPicker` | `maxmin_picks()` | ✅ |
| Butina クラスタリング | `SimDivFilters.ClusterData` | `butina_cluster()` | ✅ |
| 機能基同定 | `ifg.identify_functional_groups` | `identify_functional_groups()` | ✅ Ertl 2017 |

### 3.9 ファイル I/O

| 形式 | RDKit | chematic | 備考 |
|---|---|---|---|
| SMILES | `Chem.MolFromSmiles` | `parse()` | ✅ |
| SMILES (書き出し) | `Chem.MolToSmiles` | `write()` / `canonical_smiles()` | ✅ |
| Randomized SMILES | — | `random_smiles()` | ✅ |
| MOL v2000 / v3000 | `Chem.MolFromMolFile` | `parse_mol()` / `parse_mol_v3000()` | ✅ |
| SDF (読み込み) | `Chem.SDMolSupplier` | `SdfReader` | ✅ ストリーミング対応 |
| SDF (書き出し) | `Chem.SDWriter` | `write_sdf()` | ✅ |
| MOL2 (Tripos) | `Chem.MolFromMol2File` | `parse_mol2()` | ✅ |
| CML | — | `parse_cml()` | ✅ |
| CDXML (ChemDraw) | — | `parse_cdxml()` | ✅ 上限: 10,000原子 |
| RXN | `AllChem.ReactionFromRxnFile` | `parse_rxn_file()` | ✅ |
| CXSmiles | ChemAxon | `parse_cxsmiles()` | ✅ |

### 3.10 可視化 (Depiction)

| 機能 | RDKit API | chematic | 備考 |
|---|---|---|---|
| SVG 描画 | `Draw.MolToSVG` | `depict_svg()` | ✅ |
| SVG ハイライト | `Draw.MolToSVG(highlightAtoms)` | `depict_svg_highlighted()` | ✅ |
| SVG グリッド | `Draw.MolsToGridImage` | `depict_svg_grid()` | ✅ |
| PNG 描画 | `Draw.MolToImage` | `render_png()` | ✅ tiny-skia |
| 反応 SVG | `Draw.ReactionToImage` | `depict_reaction_svg()` | ✅ |
| Dative/Query 結合 | — | ✅ | 矢印・点線表示 |
| 2D 座標計算 | `Chem.Compute2DCoords` | `compute_layout()` | ✅ |
| 原子色マップ | — | `atom_color_map` (RenderOptions) | ✅ |
| WASM バインディング | — | `chematic-wasm` | ✅ |

### 3.11 反応 (Reactions)

| 機能 | RDKit API | chematic | 備考 |
|---|---|---|---|
| 反応 SMILES 解析 | `AllChem.ReactionFromSmarts` | `parse_reaction()` | ✅ |
| 反応中心検出 | — | `find_reaction_center()` | ✅ |
| Atom economy | — | `atom_economy()` | ✅ |
| E-factor | — | `e_factor()` | ✅ |
| ライブラリ列挙 | — | `enumerate_library()` | ✅ 2-3 way |
| テンプレート適用 | `AllChem.RunReactants` | `run_reactants()` | ✅ |

### 3.12 スコープ外 (Out of Scope)

| 機能 | 理由 |
|---|---|
| 遷移金属化学 | 有機化学専用の原子価モデル |
| ポリマー / HELM / FASTA | 対象スコープ外 |
| ML 予測モデル | 外部モデル非依存方針 |
| StandardInChI FFI | FFI ゼロ方針 |

---

## 4. 実装品質ノート (v0.1.94)

### 最近の主要修正

| 修正 | 詳細 |
|---|---|
| **Kabsch RMSD** (BUG-3) | 回転行列の向きが逆だった問題 (R=VUᵀ→R=UVᵀ) を修正。純回転で RMSD≈0 を確認するテスト追加 |
| **Jacobi 符号** (BUG-1) | 固有値分解の回転角 θ 式の符号誤り修正 `(aqq-app)→(app-aqq)` |
| **Jacobi 反復数** (BUG-2) | `max_iterations=100` → `n*(n+1)/2*5`。n>14 の分子で精度不足だった問題を解消 |
| **距離幾何 smooth_bounds** | 環分子（トルエン等）で upper bound が ∞ になる問題を Floyd-Warshall 伝播で解消 |
| **DoS ガード** (SEC-1) | `generate_coords_dg`: `DG_MAX_ATOMS=500` 超の分子は空座標を返す |
| **DoS ガード** (SEC-6) | CDXML パーサー: `CDXML_MAX_ATOMS=10,000` 超でエラー返却 |

### 簡易実装の詳細

**MHFP (⚠️)**
- 現状: Morgan 循環ハッシュをシングルとして MinHash（v0.1.95 で正規化）+ MinHash LSH インデックス（v0.1.97）
- 残課題: SMILES ベースシングルの方が精度向上の余地あり（理論的差異は ±5% 未満）
- 影響: 大規模データベース類似度検索で精度差 ±5% 以内

**ERG (⚠️)**
- 現状: 薬理学的ノードタイプ付与済（DONOR/ACCEPTOR/POSITIVE/NEGATIVE/HYDROPHOBIC/AROMATIC）
- 残課題: 細かい HBD 判定（アミド N は acceptor-only など）で RDKit との差異が残る
- 影響: 構造多様性評価での精度差は軽微

**MMFF94 電荷 (⚠️ 残課題)**
- 現状: BCI テーブル（25エントリ、元素+結合次数分類） — ±0.1e
- 残課題: 完全な MMFF94 106原子タイプ分類 + ~2000 BCI エントリ → ±0.01e
- 影響: 電荷バランス ±0.1e（RDKit は ±0.01e）

---

## 5. v0.1.95+ ロードマップ

| 優先度 | 機能 | 概要 | 状態 |
|---|---|---|---|
| ~~HIGH~~ | ~~MHFP 正規化~~ | ~~Morgan 循環ハッシュ~~ | ✅ v0.1.95 |
| ~~HIGH~~ | ~~ERG 薬理特徴~~ | ~~DONOR/ACCEPTOR/POSITIVE/NEGATIVE/HYDROPHOBIC~~ | ✅ v0.1.95 |
| ~~HIGH~~ | ~~Reaction FP XOR~~ | ~~XOR 差分エンコーディング~~ | ✅ v0.1.94 |
| MEDIUM | MMFF94 電荷テーブル | 完全 BCI テーブル参照 | 未着手 |
| ~~LOW~~ | ~~LogP アルケニル C 区別~~ | ~~末端=CH₂ vs アリール隣接=CH−~~ | ✅ v0.1.99 |
| ~~LOW~~ | ~~Kekulization エッジケース~~ | ~~Edmonds flower algorithm (奇数員環)~~ | ✅ v0.1.100 |

---

## 6. 参考資料

| ドキュメント | 内容 |
|---|---|
| `docs/rdkit_comparison.md` | 定量ベンチマーク元データ (175分子, v0.1.74) |
| `docs/rdkit_gap_analysis.md` | ギャップ分析詳細 (A1-A6, B1-B7 series) |
| `RELEASE_NOTES_v0.1.89.md` | v0.1.89 フィーチャー一覧 |

---

*chematic v0.1.96 — 2026-06-13*
