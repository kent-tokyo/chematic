# Changelog（日本語）

All notable changes to chematic will be documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

v0.1.8 以前の変更履歴は [CHANGELOG.md](CHANGELOG.md) を参照。

---

## [Unreleased]

---

## [0.2.9] — 2026-06-14

### Added — MMFF94 幾何最適化器（Halgren 1996 完全力場）

**`crates/chematic-ff/src/mmff94_minimizer.rs`** (新規, ~590 行):

- `mmff94_total_energy(mol, coords)` — 5 項全エネルギー評価
- `minimize_mmff94_full(mol, coords, max_iter)` — steepest descent 幾何最適化
- `MinimizeResult { energy, rmsd, converged, iterations }`

**実装エネルギー項:**

| 項 | 式 |
|----|----|
| 結合伸縮 | `(143.9325×kb/2)×ΔR²×(1−cs×ΔR+(7/12)cs²ΔR²)` (cubic補正, cs=2.0) |
| 角度変形 | `(0.043844×ka/2)×Δθ²×(1−0.007×Δθ)` (Δθ in degrees) |
| **二面角 (初実装)** | `(v1/2)(1+cosφ)+(v2/2)(1−cos2φ)+(v3/2)(1+cos3φ)` |
| vdW | buffered 14-7: `ε×t⁷×(t⁷−2)`, Slater-Kirkwood 結合則 |
| 静電 | `332.0716×qi×qj/(r+0.05)`, 1-4 スケーリング 0.75 |

- 有限差分勾配 (δ=1e-4 Å)、収束閾値 1e-4、1-2/1-3 除外、1-4 静電 0.75
- テスト: +6（gauche/anti torsion エネルギー差、vdW 近距離反発、二面角、最適化でエネルギー低下）
- **合計: 1,947 テスト、全通過**

---

## [0.2.8] — 2026-06-14

### Added — MMFF94 全エネルギーパラメータ（Halgren 1996 Tables IV–VII）

**`crates/chematic-ff/src/mmff94_energy.rs`** (新規, ~4,000 行):

RDKit `Code/ForceField/MMFF/Params.cpp` (BSD ライセンス) から `gh api` でダウンロードして逐語抽出。Copyright © Merck and Co., Inc., 1994–1996。

| テーブル | エントリ数 | インデックス |
|---------|-----------|-------------|
| `MMFF94_BOND_ENERGY` (Table IV) | 493 | (bond_type, type_i, type_j) |
| `MMFF94_ANGLE_ENERGY` (Table V) | 2,245 | (angle_type, type_i, type_j, type_k) |
| `MMFF94_TORSION_ENERGY` (Table VI) | 926 | (tors_type, type_i, type_j, type_k, type_l) |
| `MMFF94_VDW_ENERGY` (Table VII) | 95 | type_i |

- ソート済みスライスに対する O(log n) バイナリサーチ
- Torsion wildcard フォールバック階層: exact → reversed → wildcard ends → both wildcards → generic
- **RDKit Python API クロス検証**: C-C-C-C torsion v1=0.103/v2=0.681/v3=0.332 ✅、C-C bond kb=4.258/r0=1.508 ✅
- **PBCI 99件全値一致** vs Params.cpp ✅
- テスト: +11
- **合計: 1,941 テスト、全通過**

---

## [0.2.7] — 2026-06-14

### Added — Canonical SMILES ステレオパリティ補正 + MMFF94 忠実電荷

#### Canonical SMILES ステレオパリティ補正（RDKit issue #8775 先行解決）

**`crates/chematic-smiles/src/parser.rs`**:
- `StereoEntry` enum で parse 時の隣接原子順序を記録（`Atom`・`ImplicitH`・`PendingRing` の 3 種）
- リングクロージャー完了時に `PendingRing` を解決し `set_stereo_neighbor_order()` で保存

**`crates/chematic-smiles/src/canonical.rs`**:
- `corrected_chirality()` — parse 時順序 vs canonical 書き出し順序を比較、奇置換を検出して `@`/`@@` を自動反転

**`crates/chematic-core/src/molecule.rs`**:
- `stereo_neighbor_order: HashMap<AtomIdx, Vec<u32>>` フィールド追加
- `STEREO_H_SENTINEL: u32 = u32::MAX` 定数

検証: L-アラニン（N優先 vs C優先書き出し）、アミノシクロペンタン、フルオロシクロヘキサン。全テスト hard-assert。

#### MMFF94 忠実電荷（Halgren 1996 equation 15）

**`crates/chematic-ff/src/mmff94_numeric.rs`** (新規, ~1,300 行):
- `assign_mmff94_numeric_types()` — 環サイズ考慮芳香族型付け (5 員環 C→37/38, 6 員環 C→63, 5 員環 N→40/58, 6 員環 N→67)
- `MMFF94_PBCI` 99 件、`MMFF94_CHG` 498 件（RDKit Params.cpp 由来）
- `mmff94_charges_numeric()` — equation 15 完全実装
- グリシン参照ログ照合: C-O 結合 O→−0.28 ✅、N-H 結合 N→−0.36 ✅

テスト: +15
- **合計: 1,930 テスト、全通過**

---

## [0.1.102] — 2026-06-13

### Added — IUPAC 命名拡充 Round 2: チオール・アルコール位置ロカント・二置換ベンゼン・メチルシクロアルカン

- **チオール** (`name_thiol()`): `CS` → "methanethiol"、`CCS` → "ethanethiol"
- **アルコール位置ロカント**: `CCCO` → "propan-1-ol"（旧 "propanol"）、`CC(O)C` → "propan-2-ol"
- **二置換ベンゼン** (`name_disubstituted_benzene()`): ortho(2)/meta(3)/para(4) プレフィックスを環 BFS 距離で決定
- **メチルシクロアルカン**: `CC1CCCCC1` → "methylcyclohexane"
- テスト数: 1,691 → 1,695 (+4)、全通過。

---

## [0.1.101] — 2026-06-13

### Added — IUPAC 命名拡充 Round 1

- **分枝アルカン**: `name_branched_alkane()` + `find_longest_c_chain()` + `format_substituents()` を新設。IUPAC 最小ロカント規則で主鎖選択、メチル/エチル置換基に対応。例: `CC(C)C` → "2-methylpropane"、`CC(C)(C)C` → "2,2-dimethylpropane"。
- **一置換ベンゼン慣用名**: `name_monosubstituted_benzene()` を新設。OH → phenol、CH3 → toluene、NH2 → aniline、Cl → chlorobenzene、Br → bromobenzene、CHO → benzaldehyde、COOH → benzoic acid、CN → benzonitrile。
- **ニトリル**: `is_nitrile()` + `name_nitrile()` で R-C≡N を検出し "ethanenitrile" 形式で命名。
- テスト数: 1,686 → 1,691 (+5 テスト関数 / 12+ アサーション)、全通過。

---

## [0.1.100] — 2026-06-13

### Improved — Kekulization: エッジケーステスト + 降順フォールバック

- **エッジケーステスト 5 件追加** (`chematic-core/src/kekulization.rs`):
  - biphenylene（4 員環シクロブタジエンブリッジ）
  - anthracene（3 直線縮合 6 員環）
  - 大型 4 環 PAH（pyrene 類似）
  - biphenylene 二重結合数検証
  - 決定論チェック（同一分子 → 同一結果）
- **降順フォールバック**: 主パス（昇順）が未マッチ原子を残した場合、降順で自動再試行。コスト O(2·V·E)、グリーディ順序依存の行き詰まりを解消。
- `run_matching_pass()` ヘルパー関数を抽出して 2 パスで共用。
- 実際に 5 分子全て Pass 1 で成功→実用的有機化合物では Edmonds は不要と判明。
- テスト数: 1,681 → 1,686 (+5)、全通過。

---

## [0.1.99] — 2026-06-13

### Improved — LogP Crippen: エノン系 vinyl C 分類

- α,β不飽和カルボニル系（C=C-C=O）の内部 vinyl C を精細化。
- 既存の `neighbor_has_carbonyl()` ヘルパーを活用して「Ar隣接」「末端」「エノン」「その他」の 4 層分類を実現。
  - Ar隣接 =CH−: 0.2640（継続）
  - 末端 =CH₂: 0.1551（継続）
  - エノン =CH− (C=C-C=O): **0.1302（新規）**
  - その他内部: 0.2274（継続）
- 根拠: C=O の電子求引効果により β-vinyl 炭素の疎水性が低下。
- テスト +4（MVK、カルコン、クロトン酸、エノン vs アルケン比較）。
- テスト数: 1,677 → 1,681 (+4)、全通過。

---

## [0.1.98] — 2026-06-13

### Added — WASM API: MMFF94 電荷・MHFP フィンガープリント・MinHash LSH インデックス

`chematic-wasm` に 3 つの新しいエクスポートを追加（テスト +5）：

- `mmff94_charges_json(mol)` → `[q0, q1, ..., qN]` — MMFF94 BCI 部分電荷の配列（±0.1e）。
- `mhfp_hashes_json(mol)` → `{"num_hashes":128,"hashes":[...]}` — MinHash FP を JSON 形式で返す。
- `tanimoto_mhfp_smiles(smi1, smi2)` → `f64` — 2 つの SMILES 間の MHFP Tanimoto 類似度。
- `MhfpLshHandle` — ステートフルな JS クラス。`add_smiles(smiles)` で分子を登録し、`query_json(smiles, threshold)` で閾値以上の近傍を高速検索。
- テスト数: 1,672 → 1,677 (+5)、全通過。

---

## [0.1.97] — 2026-06-13

### Added — MinHash LSH インデックス (`chematic-fp`)

- 新型 `MhfpLshIndex` を `crates/chematic-fp/src/lsh.rs` に追加。MHFP フィンガープリントへのバンド分解型 LSH で、大規模ライブラリの近似類似度検索を準線形時間で実行。
- デフォルト: 16 バンド × 8 行（128ハッシュ対応）、`with_bands(bands, rows)` でカスタマイズ可能。
- `add(fp) -> usize` でフィンガープリント登録、`query(fp, threshold) -> Vec<(usize, f64)>` で閾値以上の全エントリをソート済みで返す。
- 再現率: P(s=0.8) ≈ 94%、P(s=0.7) ≈ 90%（確率的な近似検索）。
- テスト +6（空クエリ、自己類似度、閾値フィルタリング、類似分子発見、降順ソート、カスタムバンド）。
- テスト数: 1,666 → 1,672 (+6)、全通過。

---

## [0.1.96] — 2026-06-13

### Improved — MMFF94 BCI 部分電荷

- `mmff94_charges()` を電気陰性度近似（±0.5e）から Bond Charge Increment (BCI) テーブルベース（±0.1e）に置換。
- 新モジュール `chematic-chem/src/mmff94_bci.rs`: C=O (+0.47)、C–O (+0.04)、C–N (+0.10)、C–N(アミド) (+0.31)、N–H (-0.16)、O–H (-0.30)、C–F (+0.22) など 25 エントリ（Halgren 1996）。
- 電荷保存を保証: `sum(q) == sum(形式電荷)`。
- テスト +9（エタノール、アセトン、メチルアミン、酢酸、アンモニウム、アセタート、クロロメタン、イミダゾール、アミド vs アミン BCI 比較）。
- テスト数: 1,657 → 1,666 (+9)、全通過。

---

## [0.1.95] — 2026-06-13

### Added — `chematic-iupac` ローカル命名範囲の拡張

- 位置番号付きケトンに対応: `propan-2-one`、`butan-2-one`、`pentan-3-one`。
- カルボン酸、エステル、第一/第二アミドに対応: `ethanoic acid`、`methyl ethanoate`、`ethanamide`。
- 非置換ベンゼンと主要な芳香族ヘテロ環に対応: ピリジン、フラン、チオフェン、ピロール、イミダゾール、ピリミジン。
- IUPAC ユニットテストを 8 → 14 に拡張。

### Fixed — CI Clippy 互換性

- `cargo clippy --workspace -- -D warnings` が現在の stable Clippy で通過。
- deprecated な `total_hcount` は crate root で互換維持しつつ、互換 re-export のみ警告を抑制。
- 新しい Clippy lint に合わせて ECFP 反復ループ、condensed formula のガード、DG 座標生成器の doc コメントを更新。

### Improved — 真フィンガープリントアルゴリズム (A4/A5/A6)

- **MHFP 正規化ハッシュ**: アトムインデックス依存の byte 署名を Morgan 式循環フラグメントハッシュに置換。同一分子を異なる SMILES 順で入力しても同一 FP を生成。テスト +3。
- **ERG 薬理ノードタイプ**: `assign_pharmacophore_features()` を新設。DONOR (N-H/O-H)・ACCEPTOR (無 H のピリジン型 N、O、F)・POSITIVE/NEGATIVE (形式電荷)・HYDROPHOBIC (純 C/H) を正確に付与。ピリジン N とピロール N を区別。テスト +5。
- **Reaction FP**: `use_xor: true` (XOR 差分) がデフォルトであることを確認、比較表を ⚠️ → ✅ に更新。
- テスト数: 1,649 → 1,657 (+8)、全通過。

---

## [0.1.37] — 2026-06-08

### Added — mol_transforms API + ランダム SMILES 生成

#### `chematic-3d` — 分子幾何操作

**NEW**: 結合長・角度・二面角の測定・変換の公開 API:
- `get_bond_length(coords, a, b) -> f64` — 結合長（Å）
- `get_bond_angle(coords, a, center, b) -> f64` — 角度（ラジアン）
- `get_bond_angle_deg(...)` — 角度（度）
- `get_dihedral(coords, a, b, c, d) -> Option<f64>` — 二面角（ラジアン）
- `get_dihedral_deg(...)` — 二面角（度）
- `set_dihedral(coords, mol, a, b, c, d, angle_rad) -> Coords3D` — D 側部分樹を回転
- `compute_centroid(coords) -> [f64; 3]` — 全原子の重心
- `center_on_origin(coords) -> Coords3D` — 原点に移動
- `transform_conformer(coords, 4x4_matrix) -> Coords3D` — 4×4 同次行列を適用

**内部関数の公開**:
- `dihedral()`、`compute_angle()`、`rotate_around_axis()` （従来は private）
- 新規ファイル: `crates/chematic-3d/src/mol_transforms.rs`

#### `chematic-smiles` — SMILES 多様性生成

**NEW**: ML データ拡張向けのランダム SMILES 生成:
- `random_smiles(mol, seed) -> String` — xorshift64 RNG で原子順を permute
- `random_smiles_vect(mol, count, seed) -> Vec<String>` — N 件のユニークな変異体生成
- アルゴリズム: Fisher-Yates shuffle
- ユースケース: データ拡張（同一分子の多様な SMILES 表現）

#### `chematic-wasm` — mol_transforms + ランダム SMILES の WASM バインディング

**NEW**: 4 つの WASM エクスポート関数:
- `get_bond_length_json(smiles, a, b) -> f64` — SMILES から結合長
- `get_dihedral_json(smiles, a, b, c, d) -> JsValue` — SMILES から二面角（度）
- `set_dihedral_json(smiles, a, b, c, d, angle_deg) -> Result<String>` — PDB ブロック返却
- `random_smiles_json(smiles, count, seed) -> Result<String>` — SMILES の JSON 配列

### テスト カバレッジ

- **chematic-3d**: +5 mol_transforms テスト（結合長、角度、重心、行列変換）
- **chematic-smiles**: +6 random_smiles テスト（決定性、一意性、ラウンドトリップ）
- **合計**: 1,120 → 1,151 (+31 新テスト)
- 全テスト合格 ✅

---

## [0.1.36] — 2026-06-08

### Fixed — Issue #1 監査: 位相的に正しいが化学的に無意味な結果

#### `chematic-smarts` — VF2・MCS 正確性修正

**BUG-2: SMARTS `[h]` プリミティブ（暗黙的 H 数）**
- **修正**: パーサーが `[H]`/`[H2]`（合計 H 数）と `[h]`/`[h2]`（暗黙的 H のみ）を正しく区別するようになった
- `AtomPrimitive::ImplicitHCount(u8)` バリアント を query.rs に追加
- parser.rs を修正して、要素フォールスルーの前に小文字 `h` を処理
- `eval_atom_primitive()` に `implicit_hcount()` を使用したマッチアームを追加
- **影響**: 芳香族 H が明示的原子に誤ってマッチするサイレント誤マッチを防止

**BUG-3: MCS `maximize_bonds` タイブレイク**
- **修正**: mcs.rs の `grow()` 関数を修正して、原子数が等しい場合に結合数をタイブレイク条件として使用
- 条件追加: `|| (maximize_bonds && mapping.size == best.size && mapping.bond_count > best.bond_count)`
- デフォルト `maximize_bonds=true` で RDKit 動作と同期
- **影響**: 複数の同等サイズマッチが存在する場合、MCS が一貫した結果を返すようになった

**BUG-4: SMARTS `/\` 幾何立体結合（E/Z）**
- **修正**: `BondPrimitive` enum に `Up` と `Down` バリアントを追加して幾何立体化学対応
- parser.rs を修正: `is_bond_token()` が `/` と `\` 文字を認識
- `consume_bond_prim()` を修正して `/` を `Up`、`\` を `Down` として解析
- `eval_bond_primitive()` に `BondPrimitive::Up` と `Down` マッチアームを追加
- **影響**: `/C=C\` のような SMARTS クエリが E/Z に設定された二重結合を正しくマッチできるようになった

### テスト カバレッジ

- **chematic-smarts**: 124 テスト全パス（BUG-2/3/4 検証含む）
- **ワークスペース**: 1,120+ テスト全パス
- clippy 警告なし

### 注記

- **Issue #1 パターン**: 監査で、アルゴリズムが位相的には正しいが化学的に無意味な結果を返すバグ 3 件を発見
  - 根本原因: RDKit には制約オプションがあるが、chematic では実装されていなかった
  - 例: VF2 `uniquify`（v0.1.33 で修正）、MCS 環認識（v0.1.22 で修正）
  - 本 Sprint で 3 つの欠落していた正確性制約を追加

---

## [0.1.32] — 2026-06-07

### Added — 3D ジオメトリ・座標処理・WASM 安定性の強化

#### `chematic-3d` — 距離幾何制約充足

- **NEW**: `build_constraints()` & `satisfy_constraints()` による反復制約射影法
  - 理想結合距離（±0.05 Å 許容度）と原子価角（±5° 許容度）を強制
  - 典型分子で 5～10 回の反復で収束
  - 反復あたり O(n²)；中規模分子（< 1000 原子）に適用可能
- **NEW**: `generate_and_minimize_constrained()` 高水準 API：DG → 制約充足 → DREIDING
  - ひずみのある・問題のある分子の幾何品質を改善
- **NEW**: BondConstraint & AngleConstraint 構造体と強制メソッド

#### `chematic-mol` — V3000 座標復元機能

- **NEW**: `parse_mol_v3000_with_coords()` 関数が (Molecule, MolMetadata, Vec<(f64, f64)>) を返す
  - MOL V3000 原子ブロックから 2D 座標を復元（従来は破棄していた）
  - V2000 の `parse_mol_with_coords()` API パターンに統一
  - ラウンドトリップで 2D 座標を保持できるように

#### `chematic-perception` — 芳香族性モデルの厳密化

- **NEW**: `RingAromaticity` enum が Aromatic/Antiaromatic/NonAromatic を区別
- **NEW**: Hückel 4n+2 則と反芳香族性（4n）検出
  - C、N（H-依存）、O、S 原子の π 電子数計算
  - `ring_classifications()`、`antiaromatic_rings()`、`has_antiaromaticity()` メソッド
  - 特異系検出：シクロブタジエン、シクロオクタテトラエン、アンヌレン

### Fixed — WASM ビルド信頼性・座標系の明確化

#### `chematic-3d` — WASM RNG シード設定（MD 必須）

- **CRITICAL WASM**: `fastrand = "2.4"` が wasm32-unknown-unknown ターゲット向けに `features = ["js"]` を指定
  - MD 速度初期化が暗号学的ランダム性を使用するように修正（従来は固定シード）
  - WASM ビルドが非決定的（物理的に意味のある）軌跡を生成
  - ネイティブビルドは影響なし；feature は WASM ターゲットのみ有効化

#### `chematic-depict`、`chematic-mol` — 座標系ドキュメント

- **明確化**: `compute_layout()` は SVG Y-down ピクセル座標を生成（化学的 Y-up ではない）
- **明確化**: `parse_cml()` は化学的 Y-up 規約を返す；SVG レンダリング時は Y を反転すること
- **明確化**: `parse_cdxml()` は ChemDraw Y-down（SVG 互換、変換不要）を返す
- 座標系バグを防ぐため詳細なドキュメント文字列を追加

### Changed — エラー処理の完全性

#### 13 のエラー型が `std::error::Error` トレイトを実装

- **高優先度**: `SmartsError`、`ValenceError`、`StereoError` — Display + Error を追加
- **その他**: `CmlError`、`CdxmlError`、`Mol2Error`、`RxnParseError`、`MolError`、`IupacError`、`ConformerError`、`RxnError`、`TransformError` — Error トレイトを追加
- 標準的なエラー処理パターン（`.source()`、`Box<dyn Error>` など）が使用可能に

### Test Coverage

- **+44 個の新テスト**: 制約充足（12）、芳香族性（16）、V3000 座標（2）、エラー型（14）
- **合計**: 全クレートで 171 テスト合格
- WASM ビルド検証済み（js feature での退行なし）

---

## [0.1.30] — 2026-06-07

### Fixed — 重大な物理計算・静電計算バグ修正

#### `chematic-ewald` — PME メッシュインデックス計算バグ

- **CRITICAL**: `interpolate_charges_to_mesh()` の 3D→1D メッシュインデックス計算が 97% の電荷データを喪失していた
- 間違った `isqrt()` 近似を適切な 3D インデックス公式に置き換え: `linear_idx = ix + iy*M0 + iz*M0*M1`
- 非立方形メッシュでも完全に電荷を分配できるように修正

#### `chematic-3d` — 分子動力学・力場

- **CRITICAL**: Maxwell-Boltzmann 速度初期化のユニット換算係数を修正（0.01038 kcal/mol → amu·Ų/fs²）
  - 従来コードは速度を 347 倍小さく生成していた（初期運動エネルギーが極度に低い）
  - 目標温度に対応した物理的に正しい熱速度を生成するように修正
  
- **CRITICAL**: VDW エネルギー計算をハードコード値から DREIDING パラメータを使用するように修正
  - `r_eq = 2.0 Å`（全ペア共通）を元素別 DREIDING VDW 半径に置き換え
  - Lennard-Jones 12-6 ポテンシャルの完全な実装（魅力的分散項も追加）
  - Lorentz-Berthelot 混合則を実装
  - 1-2 および 1-3 結合ペアの除外を実装

### Changed — crates.io SEO・発見性改善

#### Cargo.toml メタデータの改善

- 最上位クレートに `wasm` カテゴリを追加（crates.io ブラウズトラフィック向け）
- キーワード `rdkit` → `drug-discovery` に変更（ターゲットユーザー層との適合性向上）
- 全 15 クレートの説明文を改善（DREIDING・Velocity Verlet・SPME など具体的なアルゴリズム名を明示）
- 最上位クレートを `all-features = true` + WASM ターゲット指定で docs.rs 設定を強化
- 全 14 サブクレートに `[package.metadata.docs.rs]` セクションを追加（従来は 14/15 欠落）
- `homepage` をワークスペースレベルに昇格（重複排除）
- `chematic-iupac` に 5 番目のキーワード `systematic-name` を追加

### Test Coverage

- **+122 個の新規テスト** DREIDING・MD・SPME モジュール追加
- 合計: 1,073 テスト成功（失敗 0）
- Phase 2–3 全実装を検証済み

---

## [0.1.26] — 2026-06-06

### Added — Sprint v0.1.26–v0.1.28: 立体化学・スキャフォールド・3D類似度など大規模追加

#### `chematic-core` — 強化立体グループ

- `StereoGroup` / `StereoGroupKind`（`Absolute` / `Or(u32)` / `And(u32)`）— ChemDraw V3000 互換の AND/OR/ABSOLUTE ステレオグループ
- `Molecule::stereo_groups()`, `set_stereo_groups()`, `add_stereo_group()`
- `MoleculeBuilder::add_stereo_group()`、`from_molecule()` が stereo_groups をコピー

#### `chematic-perception` — 立体化学・環知覚 API 拡充

- `assign_ez_from_2d(mol, coords)` — 2D 座標の幾何学的外積から E/Z を割り当て
- `cip_ez_descriptor(mol, bond_idx, coords) -> Option<CipCode>`
- `ring_membership(mol) -> Vec<Vec<usize>>` — 原子ごとの所属環インデックス
- `ring_sizes_for_atom(mol, atom_idx) -> Vec<usize>`
- `is_fused_ring_system(mol) -> bool` — 縮合環システムの検出
- `validate_stereo(mol) -> Vec<StereoError>` — `ImpossibleCenter` / `ConflictingWedges` / `RedundantStereo` 検出
- `stereo_completeness(mol) -> StereoCompleteness` — 指定済み/未指定立体中心数

#### `chematic-mol` — 新フォーマット・V3000 立体グループ

- MOL V3000 `BEGIN COLLECTION` ブロック: `MDLV30/STEABS`, `MDLV30/STEOR<n>`, `MDLV30/STEAND<n>` の解析・書き出し
- Tripos MOL2: `parse_mol2()` / `write_mol2()` — `@<TRIPOS>MOLECULE`, `ATOM`, `BOND` 対応

#### `chematic-chem` — 新アルゴリズム・記述子

- `isotope_distribution(mol, resolution) -> Vec<(f64,f64)>` — 多項式畳み込みによる同位体分布（H/C/N/O/F/S/Cl/Br/I 等対応）
- `assign_cip()` でアレン軸不斉を検出・割り当て（`>C=C=C<` パターン）
- `TautomerConfig` + `canonical_tautomer_with_config()` / `enumerate_tautomers_with_config()` — ルールセット・反復上限を設定可能に
- `scaffold_network(mol) -> Vec<Molecule>` — Schuffenhauer 2007 階層スキャフォールド分解
- `schuffenhauer_parents(mol) -> Vec<Molecule>` — 直接の親スキャフォールド
- `esol_solubility(mol) -> f64` — Delaney 2004 ESOL 水溶性予測（log mol/L）
- `logd_simple(mol, ph) -> f64` / `logd_profile()` — Henderson-Hasselbalch LogD
- `randic_index()`, `zagreb_index_m1()`, `topological_distance_matrix()` — トポロジー指数

#### `chematic-fp` — キラリティ対応 FP・一括類似度検索

- `EcfpConfig::use_chirality: bool` — R/S 感受性 Morgan FP（デフォルト false、後方互換性あり）
- `nearest_neighbors(query, db, k, FpType) -> Vec<(usize, f64)>` — 線形 Tanimoto 検索
- `nearest_neighbors_from_fp(query_fp, db_fps, k)` — 事前計算 FP からの検索
- `FpType` 列挙型: `Ecfp4`, `Ecfp6`, `Ecfp4Chiral`, `Fcfp4`, `Maccs`, `TopoPath`

#### `chematic-smarts` — 同位体・キラリティプリミティブ

- `AtomPrimitive::Isotope(u16)` — `[13C]` を同位体制約としてパース
- `AtomPrimitive::Chirality(u8)` — `[@]` / `[@@]` をキラリティ制約としてパース
- `MatchConfig::use_chirality: bool` — `[@]`/`[@@]` を対象分子に照合（デフォルト false）
- `MatchConfig::use_isotopes: bool` — `[13C]` を対象分子に照合（デフォルト false）

#### `chematic-smiles` — 新ユーティリティ

- `canonical_atom_order(mol) -> Vec<usize>` — Morgan ランク DFS 順序
- `equivalent_atom_classes(mol) -> Vec<usize>` — 対称クラス番号
- `are_atoms_equivalent(mol, a, b) -> bool`
- `parse_smi_file(s)` / `write_smi_file(records)` — `.smi` タブ/スペース区切りフォーマット

#### `chematic-rxn` — 反応メトリクス・バランス確認

- `balance_check(rxn) -> BalanceResult` — 元素バランス確認と差分レポート
- `atom_economy(rxn) -> f64` — Trost 原子経済性 %
- `e_factor(waste, product) -> f64` — Sheldon E ファクター
- `pmi_rxn(all_masses, product) -> f64` — Process Mass Intensity
- `reaction_mass_efficiency(reactants, product) -> f64`
- `find_reaction_center` をクレートルートに再エクスポート

#### `chematic-3d` — アライメント・形状認識

- `align_coords(reference, mobile) -> AlignResult` — Kabsch 最適重ね合わせ
- `apply_alignment(mobile, result) -> Vec<[f64;3]>` — 変換後座標を生成
- `rmsd_no_align(a, b) -> f64` — 回転なし RMSD
- `usr_descriptors(coords) -> [f64;12]` — Ballester-Richards USR 12 モーメント記述子
- `usr_similarity(a, b) -> f64` — Soergel 距離類似度 ∈ [0, 1]

#### `chematic`（アンブレラ）— docs.rs・WASM

- `//!` モジュールドキュメント全面改訂（機能表・クイックスタート・フィーチャーフラグ表）
- `[package.metadata.docs.rs]` で `features = ["full"]` ビルドを設定
- WASM バインディング追加: `find_reaction_center_json`, `standardize_smiles`, `balance_check_json`, `nearest_neighbors_json`, `mol2_to_smiles`, `smiles_to_mol2`

### テスト

- 全クレートで約 200 件追加（前回: ~933件 → 今回: ~1133件）

---

## [0.1.25] — 2026-06-06

### Added — P2 features: 2D layout quality + stereochemistry manipulation + reaction analysis

#### `chematic-depict` — 2D Layout Quality & Metadata

- `detect_crossings(layout, mol) -> Vec<(BondIdx, BondIdx)>` — identify bond crossing pairs for layout quality assessment
- `render_svg_with_metadata(mol, layout, opts, smiles) -> String` — embed SMILES in SVG `<metadata>` tags for image-based structure recovery

#### `chematic-chem` — Stereochemistry Manipulation

- `invert_stereocenter(mol, idx) -> Molecule` — flip R↔S configuration by inverting wedge bonds (Up↔Down)
- `enumerate_stereoisomers(mol) -> Vec<Molecule>` — generate all 2^n stereoisomers from unspecified stereocenters (max 2^6 = 64)

#### `chematic-rxn` — Reaction Center Analysis

- `ReactionCenter { broken_bonds, formed_bonds, changed_atoms }` structure
- `find_reaction_center(rxn) -> ReactionCenter` — identify bonds broken/formed and atoms changed using atom_map matching

### Tests

- 865 + 70 = 935 tests, all pass

---

## [0.1.24] — 2026-06-06

### Added — P1 features: atom label generation + standardization + molecular hashing

#### `chematic-depict` — Atom Display Labels

- `HPosition` enum (Left, Right, Up, Down) for H position hints
- `AtomLabel` struct with symbol, h_count, h_position
- `atom_display_label(mol, idx) -> String` — condensed notation ("CH₃", "NH₂", "OH")
- `atom_label_with_h(mol, idx) -> AtomLabel` — structured label data with H positioning

#### `chematic-chem` — Molecule Standardization

- `StandardizeOptions { canonical_tautomer, neutralize_charges, remove_explicit_h, largest_fragment_only }`
- `standardize(mol, opts) -> Molecule` — chain transformations: largest_fragment → neutralize → remove_h → tautomer

#### `chematic-chem` — Molecular Hashing

- `mol_hash(mol) -> u64` — FNV-1a hash of canonical SMILES
- `are_identical(a, b) -> bool` — compare canonical SMILES

### Tests

- 865 + 70 = 935 tests, all pass

---

## [0.1.23] — 2026-06-06

### Added — Element API expansion + implicit H computation + aromaticity application

#### `chematic-core` — Element Radius & Implicit Hydrogen

- `Element::vdw_radius() -> f64` — Van der Waals radius (Bondi 1964 + Alvarez 2008/2013)
- `Element::covalent_radius() -> f64` — covalent radius
- `Molecule::implicit_hydrogen_count(idx) -> u8` — implicit H count via valence rules
- `Molecule::total_formula() -> String` — Hill notation including implicit H

#### `chematic-core` — Immutable Update API

- `Molecule::with_atom_aromatic(idx, aromatic) -> Molecule`
- `Molecule::with_bond_order(idx, order) -> Molecule`

#### `chematic-perception` — Aromaticity Application

- `apply_aromaticity(mol) -> Molecule` — apply aromatic flags and BondOrder::Aromatic to Kekulized structure

#### `chematic-rxn` — Alias

- `minimize_uff()` — alias for `minimize()` for discoverability

### Tests

- 877 + 9 = 886 tests, all pass

---

## [0.1.21] — 2026-06-06

### Added — Mutable Molecule API 拡張・SDF/CDXML 機能強化・DepictData with user coords

#### `chematic-core` — Mutable Molecule API 拡張

- `Molecule::with_atom_charge(idx: AtomIdx, charge: i8) -> Molecule` — 指定原子の形式電荷を変更した新 Molecule を返す
- `Molecule::with_atom_element(idx: AtomIdx, el: Element) -> Molecule` — 指定原子の元素を変更した新 Molecule を返す（chirality・hydrogen_count・aromatic フラグはリセット）

### Changed

#### `chematic-core` — 破壊的変更 (注意)

- `Molecule::with_bond_added` の戻り値を `Result<Molecule, MolError>` から `Result<(Molecule, BondIdx), MolError>` に変更。新しく追加された結合のインデックスも同時に返すようになった。

#### `chematic-mol` — SDF/MOL V2000 座標取得

- `parse_mol_with_coords(input)` を新規追加し、`parse_mol` はそのラッパーに変更。V2000 atom block の x/y 座標（bytes 0–19）を `Vec<(f64, f64)>` として返す。
- `parse_sdf_with_coords(input) -> Result<Vec<(Molecule, MolMetadata, Vec<(f64, f64)>)>, MolParseError>` を追加。

#### `chematic-mol` — CDXML 複数フラグメント対応

- `parse_cdxml_all(input) -> Result<Vec<(Molecule, Vec<(f64, f64)>)>, CdxmlError>` を追加。`<fragment>` 要素ごとに独立した Molecule を返す。
- `parse_cdxml()` は `parse_cdxml_all` の wrapper に変更（最初の fragment のみ返す互換 API を維持）。

#### `chematic-mol` — CDXML 立体化学読み取り

- `<b>` 要素の `Display` 属性を読み取り、くさび結合を BondOrder に変換:
  - `"WedgeBegin"` / `"WedgedHashBegin"` → `BondOrder::Up`
  - `"Hash"` / `"Dash"` / `"WedgeEnd"` / `"WedgedHashEnd"` → `BondOrder::Down`

#### `chematic-depict` — DepictData with user coordinates

- `depict_data_with_coords(mol: &Molecule, coords: &[(f64, f64)]) -> DepictData` を追加。ユーザーが用意した 2D 座標から DepictData を生成する（`compute_layout` を呼ばない）。
- `compute_depict_data` を内部的に `depict_data_from_layout` ヘルパーを通じて実装するよう整理。

#### WASM 新規エクスポート

- `mol_with_atom_charge(mol, idx, charge)` → `MolHandle`
- `mol_with_atom_element(mol, idx, element_symbol)` → `MolHandle`
- `cdxml_to_smiles_json(cdxml)` → 全フラグメントの canonical SMILES の JSON 配列
- `mol_block_coords_json(mol_block)` → V2000 MOL の 2D座標 JSON `[[x,y],...]`
- `depict_data_with_coords_json(mol, coords_json)` → ユーザー指定座標で DepictData JSON を生成

### Tests

- 869 tests、全パス（前版 863 から +6）

---

## [0.1.20] — 2026-06-06

### Added — Sprint V〜CC: WASM 機能拡充・ファイル形式・編集 API

#### WASM API（84 → 103 エクスポート）

**Sprint V — Scaffold / Tautomer / 標準化 / MACCS / 一括記述子 / MOL 2D座標**
- `murcko_scaffold`, `generic_murcko_scaffold`, `canonical_tautomer`, `enumerate_tautomers_json`
- `largest_fragment`, `neutralize_charges`
- `maccs_bitvec`, `tanimoto_maccs`, `get_descriptors_json`（40+ 記述子を JSON で一括返却）
- `to_mol_block` 2D座標修正（`compute_layout` + スケーリングで実座標を出力）

**Sprint W — PAINS / CIP / ECFP6 / Dice / 3D形状記述子 / MaxMin・Butina / MCS**
- `pains_matches_json`, `cip_assignments_json`
- `ecfp6_bitvec`, `tanimoto_ecfp6`, `dice_ecfp4`, `dice_maccs`
- `shape_descriptors_json`（PMI, NPR, asphericity, eccentricity, radiusOfGyration）
- `maxmin_picks_ecfp4_json`, `butina_cluster_ecfp4_json`
- `mcs_smiles_json`

**Sprint X — V3000読み込み / 3D最小化 / SDF プロパティ / SMARTS ハイライトグリッド**
- `mol_from_v3000_block`, `generate_3d_minimized_pdb`
- `sdf_to_records_json`（name + properties の JSON 配列）
- `depict_svg_grid_highlighted`（SMARTS マッチ原子を黄色ハイライト）

**Sprint Y — XYZ/PDB I/O / per-atom 記述子 / SSSR / カスタム ECFP / 立体異性体列挙**
- `mol_from_xyz`, `to_xyz`, `mol_from_pdb`
- `logp_per_atom_json`, `mr_per_atom_json`, `labute_asa_per_atom_json`
- `sssr_rings_json`（原子インデックス配列の JSON 配列）
- `ecfp_bitvec_custom(mol, radius, nbits)`
- `enumerate_stereo_isomers_json`（未指定立体中心の全異性体、上限 64 組）

**Sprint Z — BRICS SMILES / FP bitvec / FCFP6 / SDF 書き込み**
- `brics_fragments_json`（SMILES 配列）
- `atom_pair_bitvec`, `torsion_bitvec`（各 256 bytes）
- `tanimoto_fcfp6`
- `sdf_from_records_json`（プロパティ付き SDF 書き出し）

**Sprint AA — FCFP4/6 bitvec / Dice ECFP6 / write_smiles / 反応正規化**
- `fcfp4_bitvec`, `fcfp6_bitvec`
- `dice_ecfp6`
- `write_smiles`（非正規化 SMILES）
- `normalize_reaction_smiles`

**Sprint BB — ConformerEnsemble / R-group 分解**
- `ConformerHandle` クラス: `add_generated_conformer`, `add_minimized_conformer`, `get_conformer_pdb`, `conformer_rmsd`
- `rgroup_decompose_json(smiles_json, core_smarts)` → `[{"matched":true,"r1":"..."}]`

**Sprint CC — MMP 分析**
- `mmp_pairs_json(smiles_json)` → `[{"mol_a":"...","mol_b":"...","core":"...","fragment_a":"...","fragment_b":"..."}]`

**CML / CDXML ファイル形式**（ゼロ外部依存の手書き XML パーサー）
- `mol_from_cml`, `to_cml`（CML 読み書き）
- `mol_from_cdxml`（ChemDraw XML 読み込みのみ、書き込みは仕様非公開のため未実装）

**Mutable Molecule API**
- `mol_with_atom_added(mol, element_symbol)` → MolHandle
- `mol_with_bond_added(mol, a, b, order)` → MolHandle
- `mol_with_atom_removed(mol, idx)` → MolHandle
- `mol_with_bond_removed(mol, idx)` → MolHandle
- `mol_next_atom_idx(mol)` → u32

**SDF / V3000 書き込み**
- `smiles_array_to_sdf(smiles_json)` — 2D座標付き SDF 生成
- `to_mol_v3000_block(mol)` — MOL V3000 形式文字列

**DepictData**
- `depict_data_json(mol)` → `{"atoms":[{"idx","element","x","y","label","color"}],"bonds":[{"idx","atom1","atom2","kind"}]}`
  egui / HTML5 Canvas などカスタムレンダラー向け構造化描画データ

**CPK カラー**
- `cpk_color(element_symbol)` → CSS hex 文字列

#### Rust ライブラリ拡張

**`chematic-core`**
- `Molecule::with_atom_added(&self, atom)` → `(Molecule, AtomIdx)`
- `Molecule::with_bond_added(&self, a, b, order)` → `Result<Molecule, MolError>`
- `Molecule::with_atom_removed(&self, idx)` → `(Molecule, Vec<Option<AtomIdx>>)`
- `Molecule::with_bond_removed(&self, idx)` → `Molecule`

**`chematic-mol`**
- `cml` モジュール新規: `parse_cml`, `write_cml`, `CmlError`
- `cdxml` モジュール新規: `parse_cdxml`, `CdxmlError`
- `write_sdf(records)` — 複数分子 + メタデータの SDF 書き出し
- `write_mol_v3000(mol, meta, coords)` — MOL V3000 ライター

**`chematic-depict`**
- `DepictData`, `DepictAtom`, `DepictBond`, `DepictBondKind` 構造体新規
- `compute_depict_data(mol) -> DepictData`
- `RenderOptions::with_cpk_colors_for(mol)` — CPK カラー一括設定
- `atom_color(atomic_number) -> &'static str` を `pub` に昇格

**`chematic-chem`**
- `mmp` モジュール新規: `find_mmp(mols) -> Vec<MmpPair>`
- `chematic-smiles` を `dev-dependencies` から `dependencies` に昇格（MMP の canonical_smiles 使用のため）

### Changed

- `criterion` dev-dependency: 0.5 → 0.8（`chematic-fp`, `chematic-smiles`）
- README: ゼロ C/C++ 依存の訴求を強化、WASM バイナリサイズ比較（~550 KB vs RDKit.js ~30 MB）を明記

### Tests

- 863 tests、全パス（前版 736 から +127）

---

## [0.1.19] — 2026-06-02

### Added — Sprint U: WASM 利便性 API

**SMILES-string-in 系フリー関数** (`crates/chematic-wasm/src/lib.rs`):
- `smiles_to_svg_highlighted(smiles, atoms, bonds, color)` — SMILES 文字列から直接ハイライト SVG を 1 コール生成（JS: `Uint32Array` で原子・結合インデックスを渡す）
- `match_smarts_smiles(smiles, smarts)` — SMILES + SMARTS 文字列のみで SMARTS マッチング（`parse_smiles` + `smarts_match_atoms` の 1-call wrapper）
- `tanimoto_smiles(smiles1, smiles2)` — SMILES 文字列のみで Tanimoto 類似度計算（ECFP4）
- `mol_block_from_smiles(smiles)` — SMILES から直接 MOL V2000 ブロック生成

**結合情報 API**:
- `get_bond_info(mol, bond_idx)` → `{"bondOrder":1.5,"isAromatic":true,"isInRing":true,"atomFrom":0,"atomTo":1}`
- `get_bond_between(mol, atom1, atom2)` → 同 JSON + `bondIdx` フィールド。原子インデックスペアから結合を検索（SMARTS マッチ結果からの自然なフロー）

**`get_atom_info` 拡張**:
- `totalHydrogens` フィールドを追加（明示的 H + 暗黙的 H の合計）

**InChIKey は未実装**（C ライブラリ級の複雑さのため Phase 3 以降）

---

## [0.1.18] — 2026-06-02

### Added — Sprint T: API

**per-atom カラーハイライト** (`crates/chematic-depict/src/svg.rs`, `crates/chematic-wasm/src/lib.rs`):
- `RenderOptions.atom_color_map: HashMap<AtomIdx, String>` — 原子ごとに異なる色で円ハイライト
- `DepictOptions.set_atom_color(idx, color)` — WASM API。`set_highlight_atoms` と共存可能（per-atom 色が優先）

**名前付き官能基検出** (`crates/chematic-chem/src/named_groups.rs` 新規):
- `detect_named_functional_groups(mol) -> Vec<NamedGroup>` — 20グループの SMARTS パターンテーブル
- 返却: `{"name":"hydroxyl","atoms":[3]}` 形式の JSON 配列（WASM: `detect_functional_groups(mol)`）
- カルボン酸 → carboxyl + hydroxyl + carbonyl のように重複グループを全列挙。JS 側でユニーク化可能

**原子情報取得** (`crates/chematic-wasm/src/lib.rs`):
- `get_atom_info(mol, idx) -> String` — `{"element":"C","hybridization":"sp2","charge":0,"isAromatic":false}`
- 混成軌道 (sp/sp2/sp3) を結合次数から計算。範囲外 idx → `"null"`

**MOL V2000 出力 WASM バインド** (`crates/chematic-wasm/src/lib.rs`):
- `to_mol_block(mol) -> String` — MOL V2000 形式文字列。座標はすべて 0.0

---

## [0.1.17] — 2026-06-01

### Changed (`chematic-chem`) — Sprint S: SA スコア フラグメントテーブル実装

**SA スコアのフラグメント頻度テーブルを実データに置き換え** (`crates/chematic-chem/src/sa_score.rs`):
- 従来: 10 件のダミーエントリ（任意の u32 ハッシュ、意味のないスコア）
- 新規: 145 分子の検証済みコーパスから生成した 1034 件のリアルエントリ（u64 FNV-1a ハッシュ、i16 対数頻度スコア）
- ハッシュ互換性修正: 旧実装の非公開 32-bit FNV-1a を廃止し、`chematic_fp::morgan_fp_counts` を直接使用（ECFP と同一スキーム）
- スコアエンコード: `i16 = (log10(freq_in_corpus) × 1000.0) as i16`; デフォルト -5000（テーブル未登録断片）
- 検索: ソート済みスライスへの `partition_point` バイナリサーチ（O(log 1034)）

### Added (`tools/gen_sa_table`) — コーパスからテーブルを再生成するオフラインツール

**新規ツール** (`tools/gen_sa_table/`):
- 145 件の検証済み SMILES（chematic テストスイート + デモプリセット + 既知医薬品）を内蔵
- `morgan_fp_counts(mol, 2)` を呼び出し、分子横断の断片頻度を計算
- ソート済み `static FRAGMENT_SCORES: &[(u64, i16)]` を標準出力に出力
- ファイル引数で任意の SMILES コーパスにも対応（ChEMBL など）

### Tests (`chematic-chem`)
- `taxol_harder_than_aspirin` — Taxol (SA スコア高) > Aspirin (SA スコア低) の順序確認

---

## [0.1.16] — 2026-06-01

### Fixed (`chematic-smiles`) — Sprint R: E/Z 二重結合立体化学 SMILES 出力

**正規 SMILES ライターの E/Z 方向バグを修正** (`crates/chematic-smiles/src/canonical.rs`):
- `write_chain()` の子ボンド方向修正: DFS トラバーサル方向が保存方向と逆の場合（`bond.atom1 == nb`）に Up/Down を反転するように修正。修正前は `F/C=C/Cl`（E）の正規形が Z として解釈される可能性があった
- `dfs_mark()` の環クロージャ方向修正: open atom（`neighbor`）では正しい方向の Up/Down を記録し、close atom（`atom`）では Single を記録してコンフリクトを回避

### Tests (`chematic-smiles`, `chematic-chem`)
- `test_ez_e_stable` — `C/C=C/C` の正規化が安定
- `test_ez_z_stable` — `C/C=C\C` の正規化が安定
- `test_ez_fluoro_e_stable` — `F/C=C/Cl` の正規化が安定
- `test_ez_fluoro_z_stable` — `F/C=C\Cl` の正規化が安定
- `test_ez_e_ne_z` — E と Z の正規 SMILES が異なる文字列
- `test_canonical_preserves_ez` (`cip.rs`) — 正規化後も `assign_cip` が正しい E/Z コードを返す

---

## [0.1.15] — 2026-05-31

### Added (`chematic-chem`) — Sprint Q: 官能基識別 + SA スコア + Gasteiger 電荷 + VSA 記述子

**官能基識別** (`chematic-chem/src/ifg.rs`、新規):
- `identify_functional_groups(mol) -> Vec<FunctionalGroup>` — Ertl (2017) アルゴリズム: ヘテロ原子 + 隣接 C をマーク → BFS 接続成分 = 官能基
- `FunctionalGroup { atom_indices: Vec<usize>, atom_types: String }` — 原子インデックスと元素記号文字列
- 7 テスト: ヘキサン（官能基なし）、酢酸（O あり）、ピリジン（N を含む 1 基）、アスピリン（複数）、アニリン、クロロベンゼン（Cl）

**Gasteiger-Marsili PEOE 部分電荷** (`chematic-chem/src/gasteiger.rs`、新規):
- `gasteiger_charges(mol) -> Vec<f64>` — 12 反復、ダンピング 0.5^(iter+1)
- 電気陰性度パラメータ: χ(q) = a + b·q + c·q²（C/N/O/S/F/Cl/Br/I/P/H 対応）
- 暗黙的 H を明示的 H に展開してから PEOE を実行; 重原子分の電荷のみ返す
- 5 テスト: メタノール O < C、水 O が負、電荷の合計≈0

**VSA 記述子** (`chematic-chem/src/vsa.rs`、新規):
- `slogp_vsa(mol) -> Vec<f64>` — 12 ビン (RDKit SlogP_VSA1–12)
- `smr_vsa(mol) -> Vec<f64>` — 10 ビン (RDKit SMR_VSA1–10)
- `peoe_vsa(mol) -> Vec<f64>` — 14 ビン (RDKit PEOE_VSA1–14)
- 各ビンに Labute ASA 寄与を集計; ビン境界は RDKit MolSurf.py と同一
- `logp_crippen_per_atom`、`mr_per_atom`、`labute_asa_per_atom` — 総和関数が移譲する per-atom 変形を追加

**SA スコア** (`chematic-chem/src/sa_score.rs`、新規):
- `sa_score(mol) -> f64` — [1, 10] 範囲; 1 = 合成容易、10 = 困難
- 複雑度成分: スピロ原子 × 0.25 + 架橋頭炭素 × 0.35 + マクロ環 × 0.30 + 不斉中心 × 0.10 + (環数−1)×0.05 + 環結合比 × 0.50 + サイズペナルティ
- **注**: フラグメントスコア成分（Ertl 2009 の断片頻度テーブル）は未実装; 現在の実装は複雑度ベースの近似

**多様性ピッキング + クラスタリング** (`chematic-chem/src/diversity.rs`、新規):
- `maxmin_picks(mols, n, sim_fn) -> Vec<usize>` — MaxMin 多様性ピッキング（最大-最小距離を繰り返し選択）
- `butina_cluster(mols, cutoff, sim_fn) -> Vec<Vec<usize>>` — Butina クラスタリング（類似度閾値ベース）
- `sim_fn: Fn(&Molecule, &Molecule) -> f64` — フィンガープリントに依存しない汎用インターフェース

Tests: 697 → 736 (+39 new tests)

### Added (`chematic-wasm`) — Sprint Q WASM バインディング

6 新規関数:
- `identify_functional_groups(mol) -> String` — JSON 配列 `[{"atoms":[0,1],"types":"CN"},…]`
- `gasteiger_charges_json(mol) -> String` — JSON 配列 `[q0, q1, …]`（重原子のみ）
- `sa_score(mol) -> f64` — 合成アクセシビリティスコア [1, 10]
- `slogp_vsa_json(mol) -> String` — JSON 配列（12 要素）
- `smr_vsa_json(mol) -> String` — JSON 配列（10 要素）
- `peoe_vsa_json(mol) -> String` — JSON 配列（14 要素）

### Added (`demo/index.html`) — Sprint Q UI 更新

- IFG（官能基識別）パネル: 分子ロードで即時更新
- 記述子テーブルに SA Score + Labute ASA を追加
- バージョンバッジ: v0.1.14 → v0.1.15

---

## [0.1.14] — 2026-05-31

### Added (`chematic-chem`) — EState インデックス

**EState インデックス** (`chematic-chem/src/estate.rs`、新規):
- `estate_indices(mol) -> Vec<f64>` — Hall & Kier (1991) 電子状態インデックス; 全重原子に対して per-atom 値を返す
- `max_estate(mol) -> f64`, `min_estate(mol) -> f64`, `sum_estate(mol) -> f64` — 集計記述子
- intrinsic state I_i = ((2/n)² · δᵛ + 1) / δ; 扰动 S_i = I_i + Σ (I_i − I_j) / r²_{ij} (BFS 距離)

### Added (`chematic-fp`) — パスフィンガープリント

**パス FP** (`chematic-fp/src/path_fp.rs`、新規):
- `path_fp(mol) -> BitVec2048` — 長さ 1〜7 の単純パスを DFS 列挙し FNV-1a ハッシュ; 2048 ビット
- `tanimoto_topo_path(a, b) -> f64` — パス FP の Tanimoto 係数

### Added (`chematic-wasm`) — Sprint P WASM バインディング

- `mol_from_sdf_block(block) -> MolHandle` — SDF/MOL V2000 ブロックから分子を生成
- `sdf_to_smiles_json(sdf) -> String` — SDF 文字列から SMILES JSON 配列
- `estate_indices_json(mol) -> String` — EState インデックスの JSON 配列
- `tanimoto_path(a, b) -> f64` — パス FP Tanimoto
- `MolHandle` に `sum_estate`, `max_estate`, `min_estate` メソッドを追加

---

## [0.1.13] — 2026-05-31

### Added (`chematic-wasm`) — panic hook + 反応 SVG

- `wasm_bindgen(start)` で `console_error_panic_hook` を設定; WASM パニックがブラウザコンソールに詳細を出力するように
- 反応 SVG に矢印（`→`）と試薬ラベルを追加

---

## [0.1.12] — 2026-05-31

### Added (`demo/index.html`) — タブ UI + 3D ビューア

- タブ切り替え UI: 2D 描画・3D ビューア・類似度・反応スキーム・薬らしさ
- 3D インタラクティブビューア: WebGL ベース（マウスドラッグで回転、ホイールズームを追加）

---

## [0.1.11] — 2026-05-31

### Added (`demo/index.html`) — SMARTS ハイライト + クリックハイライト + 反応スキーム

- SMARTS 検索結果の原子をクリックでハイライト
- SMIRKS 反応スキーム UI: 反応物→生成物の SVG 表示
- クリックで原子インデックスと元素情報を表示

---

## [0.1.10] — 2026-05-31

### Added (`chematic-wasm`) — 原子データ属性 + Kekulé 表示

- SVG 原子ラベルに `data-atom-idx` 属性を追加（JavaScript クリックハンドラ用）
- Kekulé 表示モード: 芳香族ボンドを単結合/二重結合で交互に表示
- npm bundler ターゲットビルドに修正（ES モジュール形式）

---

## [0.1.9] — 2026-05-31

### Fixed (`chematic-depict`) — 単原子 SMILES の描画

単原子分子（`"O"`, `"C"`, `"N"` 等）が空白または誤表記の SVG を返していた問題を修正。

- **`"C"` (メタン)**: 骨格式ルール（炭素はラベル不要）が孤立炭素にも適用され SVG が空になっていた → `CH4` を表示するよう修正。
- **`"O"` (水)**: ラベルが `OH2` と表示されていた → 分子式スタイル `H2O` に修正。
- 一般に atom_count == 1 の分子は Hill 記法の分子式（`H2O`、`CH4`、`NH3` 等）でラベルを表示。

### Added (`chematic-depict`) — `RenderOptions` + `render_svg_opts`

```rust
let opts = RenderOptions {
    width: Some(240), height: Some(240),
    background: "transparent".into(),
    dark: true,
    ..Default::default()
};
depict_svg_opts(&mol, &opts)
```

- `width` / `height`: SVG の `width=` / `height=` 属性を上書き（`None` = 自動）。
- `padding`: 分子外周の余白（デフォルト 20.0）。
- `background`: 背景色。`"transparent"` で背景 rect + ラベル背景を省略。
- `dark`: `true` のとき結合線を白、炭素ラベルを白に変更（ダークモード対応）。
- `highlight_atoms` / `highlight_bonds` / `highlight_color`: ハイライト機能を既存の `render_svg_highlighted` と統一。

### Added (`chematic-wasm`) — `is_valid_smiles` + `DepictOptions` + `depict_svg_opts`

**`is_valid_smiles(smiles: string): boolean`**
```js
is_valid_smiles("CCO")      // true
is_valid_smiles("")         // false
is_valid_smiles("[INVALID]") // false
```

**`DepictOptions` クラス**
```js
const opts = new DepictOptions();
opts.set_background("transparent");
opts.set_dark(true);
opts.set_width(240);
opts.set_height(240);
opts.set_highlight_atoms([0, 1]);
opts.set_highlight_color("#FF6B6B");
mol.depict_svg_opts(opts);
```

---

[Unreleased]: https://github.com/kent-tokyo/chematic/compare/v0.1.21...HEAD
[0.1.21]: https://github.com/kent-tokyo/chematic/compare/v0.1.20...v0.1.21
[0.1.20]: https://github.com/kent-tokyo/chematic/compare/v0.1.19...v0.1.20
[0.1.19]: https://github.com/kent-tokyo/chematic/compare/v0.1.18...v0.1.19
[0.1.18]: https://github.com/kent-tokyo/chematic/compare/v0.1.17...v0.1.18
[0.1.17]: https://github.com/kent-tokyo/chematic/compare/v0.1.16...v0.1.17
[0.1.16]: https://github.com/kent-tokyo/chematic/compare/v0.1.15...v0.1.16
[0.1.15]: https://github.com/kent-tokyo/chematic/compare/v0.1.14...v0.1.15
[0.1.14]: https://github.com/kent-tokyo/chematic/compare/v0.1.13...v0.1.14
[0.1.13]: https://github.com/kent-tokyo/chematic/compare/v0.1.12...v0.1.13
[0.1.12]: https://github.com/kent-tokyo/chematic/compare/v0.1.11...v0.1.12
[0.1.11]: https://github.com/kent-tokyo/chematic/compare/v0.1.10...v0.1.11
[0.1.10]: https://github.com/kent-tokyo/chematic/compare/v0.1.9...v0.1.10
[0.1.9]: https://github.com/kent-tokyo/chematic/compare/v0.1.8...v0.1.9
