# chematic を COSMolKit より強い化学ツールキットにする開発指示書

> **来歴**: 2026-08-26、ユーザーから「そのまま開発AIに渡してください」という指示とともに
> 提供された標準指示書（指示書）。以後この track に取り組む AI セッションは、本ドキュメント
> を出発点として読むこと。要約と現在のステータスは `ROADMAP.md` の
> `## COSMolKit Competitive Advantage (post-100-point ambition layer)` セクションを参照。
> 誤字 (schematic → chematic) のみ機械的に修正済み。内容そのものはユーザー提供のまま。

> **2026-08-26 時点のステータス（本ドキュメント第5節「最優先」に対応）**:
> - **issue #399 は修正済み**（PR #404、draft、未マージ）。`standardize()` パイプライン内
>   の8関数が `stereo_neighbor_order`/`bond_directions`/`stereo_groups` を再構築時に
>   silently drop していた問題。開発用corpus 1,134失敗→128失敗（#392以前の水準まで復帰）。
>   NCI holdout（4,999件、blind, run once）でステレオ関連の失敗は0件。詳細は issue #399 の
>   コメント履歴（診断コメント→修正→最終結果コメント）を参照。
> - **issue #395 は未解決（OPEN）のまま**。canonical SMILES のring-closure
>   explicit-bond-order記法に起因する非冪等性（10,000件corpusの約130件）。#399とは別の
>   root cause。
> - 本指示書の第5節が明記する通り、**#395が解決するまで、Phase A（COSMolKit競合スナップ
>   ショット）以降のCOSMolKit固有の作業には着手しない**。次の具体的な一手は #395 である。

---

## 1. あなたの役割

あなたはpure Rust化学情報処理ライブラリchematicの主任開発者です。

目標はCOSMolKitの機能を表面的にコピーすることではありません。

次の3条件を同時に満たし、検証可能な総合品質でCOSMolKitを上回ることが目標です。

- COSMolKitが強いPython・batch・3D・protein・ML workflowを同等以上にする
- chematicが既に強いWASM・browser・材料・反応・説明可能性をさらに伸ばす
- 化学的identity、stereo、失敗時の安全性では明確に上回る

「APIが存在する」「テストが通る」だけで勝ったと判定してはいけません。実分子corpus、独立oracle、未使用holdout、性能測定で評価してください。

## 2. 競合の基準を固定する

作業開始時にCOSMolKitの比較対象を固定してください。

最低限、以下を記録します。

- crates.io最新版
- PyPI最新版
- GitHub mainのHEAD SHA
- Rust crate version
- Python package version
- README
- dev/parity_scope.md
- docs.rs公開API
- CI状態
- 公開benchmarkやparity fixture
- ライセンス

比較レポートには必ず、

```
COSMolKit comparison snapshot:
version:
commit:
checked_at:
```

を記載してください。

競合が更新された場合も、作業途中で基準を入れ替えてはいけません。新しい版との比較は別roundで行います。

## 3. 勝敗の定義

以下の5軸で評価します。

| 軸 | 重み |
|---|---|
| 化学的正確性・identity | 30 |
| 3D・force field・conformer | 20 |
| Python・batch・ML workflow | 20 |
| Rust・WASM・deployment | 15 |
| 機能範囲・資料・運用品質 | 15 |

総合点だけでごまかさず、各軸の勝敗を個別に出してください。

### 勝利条件

以下をすべて満たした場合のみ、「総合的にCOSMolKitを上回った」と報告してください。

- 化学的identityで非劣性ではなく優位
- 主要descriptor／fingerprint／SMARTSで同等以上
- 3D成功率とstereo安全性で同等以上
- Python batch usabilityで同等以上
- Rust API品質で同等以上
- WASM／browserでは明確に優位
- COSMolKitにあってchematicにない主要機能について、実装済みまたは明確な代替価値がある
- public releaseと第三者再現可能なbenchmarkが存在する

## 4. 現在の競合認識

COSMolKitには、少なくとも次の強みがあります。

**COSMolKitの強み**

- Rust core＋Python package
- copy-on-write／value-style molecule operations
- 明示的なin-place API
- 構造化されたunsupported error
- Python向けarray access
- NumPy／PyTorch／model-building志向
- ordered batch processing
- record単位の失敗保持
- 並列batch処理
- progress reporting
- SDF indexed datasetとchunk読み込み
- SMILES／SDF／MOL2／XYZ
- PDB／mmCIF
- Protein high-level API
- ETDG／ETKDG系preset
- multi-conformer generation
- UFF／MMFF
- Morgan／MACCS
- SMARTS／substructure
- InChI
- RDKit parityを前面に出した検証

COSMolKitは、公開変換のvalue semantics、batch-native処理、array-oriented data accessを製品アイデンティティにしています。

**chematicの現在の勝ち筋**

- WASM／browser
- Rust／Python／WASMで同一coreを利用
- typed failureとfail-closed設計
- detailed provenance
- canonical／standardization／parent identity
- reaction chemistry
- crystal／periodic／Ewald／材料系
- CDXML／CMLなどを含む広い化学I/O
- stereo-safe 3D
- connectivity-ordered初期座標
- 説明可能な変換
- 軽量deployment
- chemical stack全体との連携

この優位を捨てて、COSMolKitのPython cloneになってはいけません。

## 5. 最優先：正しさの負債を解消する

COSMolKitとの機能競争へ進む前に、現在openのcanonical identity問題を解決してください。

優先順位：

1. issue #399 — **完了（2026-08-26、PR #404 draft）**
2. issue #395 — **未解決、次の一手**
3. standardize／canonicalの未使用holdout
4. InChIKey／CIP／graph identityの独立確認

必須条件：

- bare canonical idempotency：既存corpusで100%
- standardized canonical idempotency：既存corpusで100%
- NCI first_5K.smi 4,999件holdoutの結果を報告
- mirror imageがcollapseしない
- atom permutation invariant
- isotope保存
- E/Z保存
- aromaticity／bond order保存
- failure ceilingの数字だけを変更してgreenにすることは禁止

この段階を完了するまで、新しいprotein機能やfingerprint種類の追加を開始しないでください。

## 6. Phase A：COSMolKitとの差分表を作る

コードを書く前に、競合比較表を作成してください。

ファイル：

```
docs/competitive/cosmolkit_gap_matrix.md
```

最低限、以下を比較します。

**Core chemistry**
molecule model / immutability／value semantics / explicit mutation API / sanitization / valence / aromaticity / Kekulization / ring perception / symmetrized SSSR / tetrahedral stereo / double-bond stereo / atropisomer／non-tetrahedral stereo / isotope / atom maps / substance groups

**Notation and search**
SMILES parse / canonical SMILES / random SMILES / CXSMILES / SMARTS / recursive SMARTS / substructure matching / uniquify / chirality-aware matching / query serialization

**Fingerprints and descriptors**
Morgan bit / Morgan count / sparse fingerprint / MACCS / topological fingerprint / Avalon / provenance output / molecular formula / MW／exact mass / HBD／HBA / TPSA / LogP／MR / QED / ring descriptors / fragment descriptors

**3D**
DG / ETDG / ETKDG / ETKDGv3 / multi-conformer / macrocycle / stereo enforcement / UFF / MMFF94 / MMFF94s / optimization failure reporting / conformer alignment / RMSD／symmetry-aware RMSD / TFD / shape similarity

**I/O**
MOL V2000／V3000 / SDF / indexed SDF dataset / MOL2 / XYZ / PDB / mmCIF / CDXML / CML / InChI / binary serialization

**Batch and ML**
ordered batch / per-record errors / skip／keep／raise / configurable parallelism / progress / streaming / chunking / NumPy arrays / PyTorch adapters / graph export / zero-copy possibility / stable schema / dataset fingerprinting

**Protein**
PDB parsing / mmCIF parsing / chain iteration / residue iteration / atom selection / ligand extraction / neighborhood selection / alternate locations / insertion codes / biological assembly / nucleic acids / protein–ligand mixed structure

**Deployment**
Rust / Python / WASM / JavaScript / browser / Web Worker / serverless / edge / bundle size / startup time / memory ceiling

各項目には次を記載します。

```
chematic status:
COSMolKit status:
evidence:
quality level:
known limitation:
winner:
priority:
```

READMEの宣伝文だけを根拠にしてはいけません。ソース、テスト、公開API、実測結果を確認してください。

## 7. Phase B：Python・batch workflowで追いつき、超える

COSMolKitの最も明確な優位の一つはbatch-native Python APIです。

chematicにも統一的なbatch APIを追加してください。

**API案**

```python
batch = MoleculeBatch.from_smiles(
    smiles,
    errors="keep",
    n_jobs=8,
)

prepared = (
    batch
    .standardize()
    .add_hydrogens()
    .compute_descriptors()
    .fingerprint_morgan()
)
```

**必須仕様**

- 入力順を保持
- record IDを保持
- raise／keep／skip
- record単位のtyped error
- stage単位のprovenance
- 決定的な出力順
- 並列数を明示可能
- cancellation
- timeout
- memory budget
- progress callback
- streaming iterator
- chunk処理
- CSV／SDF／JSONL出力
- partial failureを成功扱いしない
- Python例外とRust errorの対応表
- 同じinput＋config＋versionで同じ結果

**COSMolKitを超えるための追加要件**

COSMolKitのbatch transformをそのまま模倣するだけでは不十分です。

chematicでは各recordについて、

- どのstageが実行されたか
- 何が変更されたか
- 何が除去されたか
- どのfallbackを使ったか
- どのwarningが出たか
- なぜabstainしたか
- 入力identityと出力identity
- config hash
- library version

を取得可能にしてください。

**受け入れ条件**

- 100万SMILESのstreaming test
- peak memoryを測定
- 1／2／4／8 worker比較
- invalid record混在
- 順序再現性
- cancellation後の状態
- Python GILを不必要に保持しない
- COSMolKitと同一条件でthroughput比較
- 正しさを落として速度を稼がない

## 8. Phase C：SDF dataset workflow

COSMolKitにはSDFのbyte-range indexを作り、全体をメモリへ読み込まず個別recordやchunkを読む設計があります。

chematicにも以下を追加してください。

```rust
let dataset = SdfDataset::open(path)?;
let record = dataset.get(100)?;
for batch in dataset.batches(1024) {
    // ...
}
```

**必須機能**

- record offset index
- lazy parse
- random access
- batch iterator
- parallel parse
- malformed record recovery
- property access without full molecule parse
- V2000／V3000
- compressed input方針
- index cache
- file modification detection
- deterministic record numbering

**独自優位**

- chemical validation status
- parse warning provenance
- standardization pipelineの直接接続
- Arrow／Parquet export
- browser File API対応の設計
- WASMでのchunk parse
- reproducible dataset manifest

## 9. Phase D：ML-ready data export

COSMolKitはNumPy／PyTorch志向を明確に掲げています。

chematicはAPI名だけで対抗せず、stable schemaを設計してください。

**必須export**

- atom feature matrix
- bond index
- bond feature matrix
- coordinates
- conformer batch
- graph-level descriptors
- atom maps
- formal charges
- aromatic flags
- isotope
- chiral tags
- CIP labels
- masks
- failure metadata

**対象**

- NumPy
- PyTorch
- PyArrow
- JSON schema
- Rust ndarray互換
- WASM TypedArray

**設計条件**

- feature orderをversioned enumで定義
- dtypeを固定
- missing valueを明示
- atom orderingを記録
- canonical orderとinput orderを選択可能
- conformerとatom indexの整合性保証
- schema hash
- migration guide

**勝利条件**

同じデータをRust、Python、WASMでexportしたとき、型と値が一致すること。

## 10. Phase E：3D比較

COSMolKitはETKDG v3 preset、multi-conformer、UFF/MMFFを公開APIとして提供しています。

名称が同じだから同等だと判断してはいけません。

**比較対象**

single conformer coverage / best-of-10 coverage / stereo preservation / macrocycle coverage / UFF parameter coverage / MMFF parameter coverage / force-field convergence / catastrophic geometry / runtime / deterministic seed / failure reporting

**corpus**

現在の265分子corpus / 既存33分子connectivity corpus / NCI holdout / macrocycle subset / steroid subset / charged molecules / fused／bridged／spiro / isotope / metal-containing molecules /独立experimental conformer subset

**比較arm**

chematic legacy / chematic stereo-safe / chematic connectivity-ordered / chematic ensemble v2 / COSMolKit default / COSMolKit ETKDGv3 / COSMolKit multi-conformer / RDKit ETKDGv3

**指標**

success coverage / independently sound / silently wrong / stereo violations / worst bond ratio / gross clashes / RMSD / symmetry-aware RMSD / TFD / energy ranking / runtime / peak memory

**禁止事項**

- COSMolKitの出力を真値としない
- RDKitの出力も真値としない
- 競合に合わせるためだけにweight調整しない
- attempt数が違う比較をしない
- timeout条件が違う比較をしない
- failed moleculeを分母から除外しない

**勝利条件**

固定corpusと独立holdoutの両方で、

- coverageがCOSMolKit以上
- silently wrongが0
- stereo failure率がCOSMolKit以下
- bond／clash品質が非劣性
- experimental RMSD／TFDが非劣性
- failure reasonの情報量では優位

## 11. Phase F：Protein／structural biology

COSMolKitはProtein.from_pdb()、Protein.from_mmcif()、chain／residue／atom iterationを公開しています。

この領域を無視したまま「総合で勝った」とは言えません。

ただし、巨大なprotein toolkitを一気に作らないでください。

**Phase F1：読み込みとimmutable model**
PDB / mmCIF / model / chain / residue / atom / altloc / occupancy / B-factor / insertion code / hetero atom / water / ligand / nucleic acid

**Phase F2：selection**

```rust
chain("A")
residues(10..50)
within(5.0, ligand)
protein_only()
ligands()
waters()
atoms_by_element()
```

**Phase F3：chemical bridge**
ligandをchematic Moleculeへ変換 / bond-order assignmentのconfidence／abstention / protein–ligand contact / H-bond候補 / salt bridge候補 / metal coordination候補 / residue annotation

**chematic独自優位**

Protein APIを孤立させず、

- ligand standardization
- conformer validation
- reaction／adduct workflows
- WASM visualization
- structured provenance

へ接続してください。

**受け入れ条件**

- PDB/mmCIF reference corpus
- chain/residue counts
- altloc
- insertion code
- multiple models
- malformed records
- COSMolKit比較
- Biopython／Gemmiによる独立cross-check

## 12. Phase G：WASMとbrowserで圧倒する

COSMolKitはWASM、JavaScript bindings、browser-native workflowsを未完成項目として挙げています。ここはchematicの最重要な差別化領域です。

**必須demo**

ブラウザのみで次を行えるアプリを作成してください。

- SMILES parse
- canonicalize
- standardize
- descriptor
- fingerprint
- SMARTS search
- 2D depiction
- 3D generation
- conformer選択
- SDF drag-and-drop
- batch processing
- result export

サーバー送信なしで動作すること。

**非機能条件**

Web Worker / progress / cancellation / memory ceiling / chunk処理 / CSP対応 / no eval / deterministic seed / mobile browserの最低限検証 / bundle size計測 / cold-start計測

**勝利条件**

COSMolKitがPython環境を必要とするworkflowを、chematicではbrowserだけで完結できること。

## 13. Phase H：材料・結晶・反応で差を広げる

COSMolKitとの比較を小分子cheminformaticsだけに限定しないでください。

chematicには次の独自スタックがあります。

- chematic-crystal
- periodic structure
- Ewald
- reaction chemistry
- materials stackとの連携
- mikiwame／gugen／kizashi等への接続可能性
- 統一データモデル

次の関係を明示してください。

```
Molecule
Reaction
CrystalStructure
PeriodicSystem
ProteinStructure
ConformerEnsemble
```

完全に同じ型へ押し込まず、共通するprovenance、element、coordinates、properties、serializationを共有します。

**公開demo**

最低限、次のend-to-end例を作ります。

- SMILES → standardize → 3D → descriptor
- SDF dataset → batch fingerprint → Arrow
- CIF → symmetry／periodic analysis
- reaction SMILES → mapping／fingerprint
- PDB → ligand抽出 → molecule analysis
- browser-only local compound explorer

これにより、「COSMolKitより小分子機能が一つ多い」ではなく、分子から材料・反応・構造生物まで接続できるRust chemical platformとして勝ちます。

## 14. API品質比較

COSMolKitのvalue-style APIは強い設計です。

chematicのpublic APIを監査してください。

**確認事項**

- mutationかnew valueかが名前から分かるか
- hidden mutationがないか
- error時に入力が破壊されないか
- partial mutationがないか
- panicがpublic APIへ漏れないか
- defaultが安全か
- opt-in best pathが発見しやすいか
- Rust／Python／WASMで意味が一致するか
- serialization後もidentityが保たれるか

**改善方針**

既存APIを一括破壊変更してはいけません。

additive API → deprecation → migration guide → compatibility test → versioned schema

の順で進めてください。

## 15. 検証規律

各機能には最低限、以下を用意してください。

positive fixture / negative fixture / known-broken fixture / holdout / atom-order permutation / determinism / mirror distinctness / malformed input / timeout / cancellation / independent oracle / competitor comparison

**比較結果の表現**

次を区別してください。

exact parity / chemical equivalence / non-inferior / superior / unsupported / unmeasured

測定していない項目を「対応」と書いてはいけません。

## 16. ソース利用とライセンス

COSMolKitやRDKitの公開ソースを参照する場合：

- ライセンスを確認
- 参照箇所を記録
- 著作権表示を維持
- コピーか独立実装かを明示
- algorithm sourceをコメントと文書に記録
- 非公開／商用コードを推測しない
- benchmark結果からコードを逆推定しない

機能名やAPI名の類似だけで互換を主張しないでください。

## 17. PR運用

一つのPRには一つの主目的だけを入れてください。

推奨分割：

diagnosis / fixture／measurement / core implementation / Python binding / WASM binding / benchmark / docs / default routing

default変更とアルゴリズム実装を同じPRに入れないでください。

merge、issue close、tag、publishは明示的な指示があるまで行わないでください。

tasks/*.mdや個人用作業ログを公開repositoryへ追加しないでください。

## 18. 最初の作業

最初からprotein APIやbatch APIを実装してはいけません。

以下の順で開始してください。

**Step 1：競合スナップショット**
COSMolKitのversion、commit、公開API、docs、test、benchmarkを記録する。

**Step 2：gap matrix**
`docs/competitive/cosmolkit_gap_matrix.md`を作る。

**Step 3：再現可能な比較harness**
新規に以下を作る。

```
validation/cosmolkit_comparison/
```

内容：

- pinned environment
- input corpus manifest
- chematic runner
- COSMolKit runner
- result schema
- scorer
- report generator
- raw result保存方針

Python package、Rust crateのどちらを比較しているか明記する。

**Step 4：現在点を採点**
各領域を100点満点で評価し、根拠を記載する。

**Step 5：最大3件だけ選ぶ**
次の開発候補から、点数への寄与とリスクを算定し、最大3件を選ぶ。

- #399／#395解消
- batch API
- SDF indexed dataset
- ML graph export
- protein reader
- 3D comparison
- WASM compound explorer

ただし、#399／#395が未解決なら、原則としてその修正を最初に行う。

## 19. 初回報告形式

最後に日本語で次を報告してください。

1. 比較したCOSMolKitのversion／commit
2. COSMolKitの主要な強み
3. chematicが既に勝っている領域
4. chematicが負けている領域
5. 未測定領域
6. 100点満点の分野別比較
7. 総合点
8. 最も点数へ効く上位3施策
9. 最初に着手する1施策
10. 作成／変更ファイル
11. 実行したテスト
12. Git状態
13. PR状態
14. 未解決リスク

推測と測定結果を明確に分離してください。

## 20. 最終方針

COSMolKitに勝つための中心戦略は次の通りです。

COSMolKitのPython・batch・protein・ML-ready workflowへ追いつく。

そのうえで、

chematicのWASM・browser・材料・反応・provenance・fail-closed化学処理を組み合わせ、競合が簡単に追随できない総合chemical platformにする。

単なるRDKit互換Rustライブラリ同士の機能数競争にしないでください。

目標は、

COSMolKitより機能が一つ多い

ではなく、

小分子、反応、結晶、protein、browser、batchを、検証可能な一つのRust coreで扱える

状態です。

（提供者注記）現時点のchematicは化学機能の広さ・WASM・材料・説明可能性では勝っていますが、Python batch、indexed SDF、protein API、ML向けデータ導線ではCOSMolKitに負けています。最短で効く順番は、#399/#395 → batch API → SDF indexed dataset → ML export → proteinです。
