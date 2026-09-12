# RDKit同等精度から、独立検証での優位性へ

更新日: 2026-09-12。対象: chematic 1.0.13からの開発系列。
状態: **一部の既存コーパスは合格、A0–A6の全出口は未達**。
[ROADMAP](../ROADMAP.md) の実行計画。P0–P6は製品領域、A0–A6は精度の作業単位として維持する。
今回の更新はA1.Rの芳香族環数退行の回復、最新候補の再測定、計画と証拠範囲の同期を含む。
リリース判定や独立評価の完了を意味しない。

## 1. 目標と優先順

最初に、指定した2D操作・化学クラスでRDKit互換性を完成させる。
次に、CIPと構造保持を中心に、独立した正解に対して同等または優位かを判定する。
3Dは別プロファイルで完成させる。追加記述子はcore完了後に並行できるが、
CIP・構造同一性の解決を全記述子の完成待ちにはしない。

**直近はA0（P0の証拠契約・未使用評価）。A1.Rのopt-in候補受入れは影響範囲のbinding/workflow回帰を通過済み。その後はA1残項目 → A2 → A3 → A4 → A5最終判定 → A6。**
A5のデータ設計・評価器はA0と並行し、新しいsilent corruptionは領域を問わず最優先とする。

| 評価軸 | 判定基準 | 最終的に主張できる範囲 |
|---|---|---|
| RDKit互換性 | 固定版・同じ設定・同じ入力に対して、宣言範囲の不一致0、対応率100% | 指定操作・プロファイル・評価セットでの同等性 |
| 化学的正しさ | 出典・規則・独立レビューで正解を決め、両エンジンを同じ分母で評価 | A5の条件を満たした領域での同等／非劣性／優位 |
| 検索の有用性 | 独立した活性・関連性ラベルへのPR-AUCやenrichment | 別途評価した用途での検索品質。Morgan一致率とは別 |
| 3D品質 | 同一座標の力場一致と、立体保持・配座品質を分離して評価 | A6で合格した化学クラス・操作のみ |

100%は固定評価セットでの観測目標であり、未知の全分子への保証ではない。
芳香族性・環表現は採用モデルも比較条件に含める。
[RDKit Book: Aromaticity](https://www.rdkit.org/docs/RDKit_Book.html#aromaticity)

## 2. 現在地と証拠の限界

以下は保存済み成果物と検証コードの監査であり、今回再測定した結果ではない。
確認時HEADは `0d951486493f6baab8f3a645a0665987af335a79`、manifestは1.0.13、
作業ツリーには既存の未コミット変更がある。結果の版名だけでは現在のソースとの一致を証明できない。

| 領域 | 保存済み結果 | まだ必要な検証 |
|---|---|---|
| A0 | 8記述子の厳密診断、native固定4件、12件MW holdoutが合格。ECFP4 5,000行のRust/source-built Python/Node-WASM binding gateも合格し、Python extension/WASM hashを記録。互換API導入コミット`5aae7b23`と現行candidateの別wheel比較が成立し、schema-v2 raw rowsの再集計もcandidateで合格。予約済み4件unused holdoutも別測定した。既存2コーパスの再現可能なdedup/split準備は2,000件development＋7,737件exposed holdoutを生成するが、sealedではない | 全項目holdout、独立取得した8,000件sealed未使用評価（既存7,737件に263件を追加しても未使用にはならない）、candidate commit/tag freeze、raw/baseline scorecardの自動生成 |
| A1 core | Source-built Python + RDKit 2025.09.3との8項目は各5,000/5,000 strict、未対応0。MW最大差は約4.0e-9 Da。Fsp3のゼロ価数同位体炭素も修正済み。採用ゲートは `adopted_opt_in` | RDKit版を固定した再評価、全bindingでの採用確認、未使用評価 |
| A1拡張 | Strict rotatable bondsは5,000/5,000。branch provenance fallback後のpotential stereocentersは5,000/5,000。さらに7,737-row exposed holdoutでpotential centers 7,737/7,737、原子別FP=0/FN=0、未解決oracle行0。P=S、N+–O−–N、S(=O)(=S)、芳香族酸素橋の環境境界を追加し、fixes6では8記述子のうち7項目が7,737/7,737 strict、芳香族環数が7,736/7,737。混在芳香族性の広域修正fixes4は29残差へ退行したため採用せず、bounded chordless aromatic cycle候補を固定RDKit 2025.09.3で再測定した結果、8項目すべて7,737/7,737 strict、不一致原因0。promotion gateは`adopted_opt_in`、native defaultは維持。potential centersはfixes3時点の測定 | 未使用/sealed評価、共有認識・全bindingゲートの再検証、他のdescriptor family、独立goldでの確認 |
| A2 | CIP履歴は4,171/4,186、P系15件は未確定。canonical意味比較200/200に加え、RDKit InChIを独立identity oracleとする5,000分子×4 randomized valid spellingの構造同一性は5,000/5,000、失敗0・oracle-invalid除外0。最新ソースの#503 K=1,024は4/28 divergent components、cross-correspondence failure 0。代表残差の内部診断では、共有carrierを選ぶと相手系の唯一の方向情報を失い、代替carrierはcanonical DFSのring-close側で出力できないため、現行solverは安全にabstainしている | 4残差を安定して解くcarrier/DFS表現、CIP独立ラベル、未使用評価、異性体の誤統合検査 |
| A3 | Rust/source-built Python/Node-WASMのk=1/10/100が500クエリ×4,500ライブラリで丸め前スコアを含め0不一致。独立RDKit exhaustive oracleも全kで一致。RDKit oracle fixtureは34成功+1 unsupported-bond errorを全bindingで検証。`score >= threshold`の閾値APIと、RDKit由来の直上/同値/直下を含む96ケースのcross-bindingゲートが合格。さらにRDKit-parity芳香族性レーンで、チェックイン済み5,000行のraw Morgan residual 59件をすべて解消し、RDKit 2025.09.3と5,000/5,000 strict一致（前処理エラー0）となった | raw/provenance packetの保存、全binding・検索ゲートの再実行、実際に未使用の評価入力 |
| A4 | 固定した過去のSMARTS診断では比較可能155,633セル中155,618一致（99.9904%）、全155,651セルでは155,618一致（99.9788%）。現行ソースの直接比較をmatch順序正規化後に再集計すると、16クエリ×5,000分子の155,651セル中、145,579一致、RDKit側unsupported 10,042、chematic明示拒否18、残差12。従来の80,000/80,000表記はこの全分母を表さないため採用しない。共有symmetrized-SSSR＋cage fallbackの実験的opt-in hybrid selectorは全SMARTS回帰と代表cageケースに合格したが、全コーパス再測定は未完了。反応presence 6/6、限定V3000往復、SMIRKS product set 7/7一致 | hybrid laneの全コーパス測定、12残差と18拒否の分類、supported/unsupported SMARTSポリシー固定、原子対応match集合、反応生成物、標準化の現行再測定、typed metadataの独立比較 |
| A5 | 4件のgold候補、2件のblind欄、manifest構造検査 | 絶対正解、実際の未使用入力、非実装者レビュー、統計評価器。現状は独立評価未実施 |
| A6 | MMFF94型IDは6,681/6,698（99.76%相当、unsupported probe 1件を除く比較対象6,697件中16残差）。現行候補のstrict bond+angleは265/265（tier A 65、tier B 200、失敗0）。BCI電荷は6,665/6,698 exact、unsupported probe除外の比較対象で6,665/6,693（99.58%、残差28原子） | 全エネルギー項・勾配・最適化・配座品質。型・電荷の残差解消と、265/265を全力場合格としない |

証拠の参照先:

- A0比較: [before/after scorecard](../validation/results/descriptor-baseline-candidate-v1.0.13-measured-raw.json)、[baseline raw](../validation/results/descriptor-rdkit-diagnostics-v1.0.13-baseline-5aae7b23-raw.json)、[candidate raw](../validation/results/descriptor-rdkit-diagnostics-v1.0.13-candidate-raw.json)
- 未使用評価: [4件holdout summary](../validation/results/descriptor-rdkit-unused-holdout-v1.0.13.json)、[candidate raw](../validation/results/descriptor-rdkit-unused-holdout-v1.0.13-candidate-raw.json)
- 評価分割準備: [split manifest](../validation/rdkit_accuracy_split_v2/manifest.json)、[builder](../scripts/build_rdkit_accuracy_split.py)。既存公開コーパスのexact-SMILES dedup結果であり、2,000件development＋7,737件exposed holdout、sealedではない。263件追加して件数を8,000件にしても未使用評価にはならない。出自と開発入力との重複を確認した未使用集合を別途確保・凍結する必要があり、このsplitはsealed合格証跡には算入しない
- 記述子: [8項目診断](../validation/results/descriptor-rdkit-diagnostics-v1.0.13.json)、[MW holdout](../validation/results/descriptor-rdkit-holdout-v1.0.13.json)、[potential centers](../validation/results/descriptor-stereocenter-rdkit-parity-v1.0.13.json)、[native固定回帰](../validation/results/native-descriptor-regression-v1.0.13.json)
- A1拡張再評価: [退行前 fixes3](../validation/results/descriptor-rdkit-exposed-holdout-v1.0.13-source-2025.09.3-after-a1-fixes3.json)、[退行実験 fixes4・29残差](../validation/results/descriptor-rdkit-exposed-holdout-v1.0.13-source-2025.09.3-after-a1-fixes4.json)、[fixes6・1残差](../validation/results/descriptor-rdkit-exposed-holdout-v1.0.13-source-2025.09.3-after-a1-fixes6.json)、[chordless候補・RDKit固定再測定](../validation/results/descriptor-rdkit-exposed-holdout-v1.0.13-source-2025.09.3-after-a1-chordless.json)、[promotion gate](../validation/results/descriptor-rdkit-promotion-gate-v1.0.13-after-a1-chordless.json)、[potential centers exposed holdout](../validation/results/descriptor-stereocenter-rdkit-exposed-holdout-v1.0.13-source-2025.09.3-after-a1-fixes3.json)
- 立体・同一性: [CIP履歴](validation.md#cip-rsez-label-agreement)、[並べ替え](../validation/results/cip-order-invariance-v1.0.12.json)、[canonical意味比較](../validation/results/canonical-cross-engine-a2-v1.0.13.json)、[#503現行K=1,024要約](../validation/results/ez_shared_carrier_coupling_mechanism_audit_summary_1024_2026-09-12.json)、[詳細JSONL](../validation/results/ez_shared_carrier_coupling_mechanism_audit_1024_2026-09-12.jsonl)
- 検索: [2026-09-12比較](../benchmarks/2026-09-12-similarity-search-a3-v1.0.13.json)、[binding間検索](../validation/results/rdkit-search-cross-binding-parity-v1.0.13.json)

最新のA3.4採用後再測定では、同一4,500ライブラリ/500クエリに対する
Rust・source-built Python・Node/WASMのk=1/10/100を比較し、1,500/1,500件で
結果集合、順位、丸め前スコアが一致した。これはbinding間整合性の証拠であり、
未使用入力による独立評価やRDKit検索との独立再測定を置き換えない。
- ワークフロー: [SMARTS集計](../validation/results/rdkit_ring_parity_diagnosis_summary.json)、[反応presence](../validation/results/reaction-rdkit-quality-v1.0.13.json)、[標準化履歴](../validation/results/standardization-phase1-holdout-1.0.8.json)
- 3D: [型比較](../validation/results/mmff94_type_parity_227_postfix.json)、[strict bond+angle](../validation/results/mmff94_strict_gate_remeasure_227_v1.0.13.json)、[角度gap](../validation/results/mmff94-angle-gap-diagnosis-v1.0.13.json)

canonicalの200件はRDKit 2025.09.3との再canonical化による意味比較であり、
独立したグラフ正解器の証明ではない。記述子・検索の2026.03.6結果と版を混ぜて集計しない。
native対Morganのtop-10一致は今回の保存値で72.7%だが、これは別定義間の一致率である。
nativeの定義を変えてこの数字を100%にすることは、互換性・化学精度の改善目標に含めない。

### 完了状態を修正する理由

- [A0集約ゲート](../scripts/check_rdkit_accuracy_gate.py) は8項目の診断を検査するが、
  holdoutには正の件数と失敗0しか要求しない。[holdout runner](../scripts/check_descriptor_rdkit_holdout.py)
  が実測するのはMWの±0.01 Daのみ。全8項目の厳密holdout合格ではない。
- [立体中心runner](../scripts/descriptor_stereocenter_rdkit_parity.py) の不一致数は最大50例の
  表示用配列長に依存し、例外も同じ配列に入る。今回の5件は上限未満だが、
  今後の評価では総数とサンプルを分け、空入力・古いbindingも検出する必要がある。
- 旧[検索binding runner](../scripts/rdkit_search_cross_binding_parity.py) はk=10固定で、
  スコアを6桁に切り捨てる。現在のfull-precision runnerはk=1/10/100を検証し、
  [閾値ゲート](../scripts/rdkit_search_threshold_gate.py) はRDKit由来96ケースの
  直上/同値/直下と0/1境界を全bindingで検証する。旧runnerの丸め契約は維持する。
- [A5 manifest](../validation/manifests/rdkit_accuracy_adjudication_v1.json) のgold候補は主に
  不変性・構造保持条件で、独立した絶対ラベルではない。blind欄も公開済み共鳴コーパスを参照し、
  reviewerは未記入。構造検査合格と独立gold合格を分ける。

A0・A3は「既存の限定ゲート合格、全出口は作業中」に戻す。
既存成果を取り消す変更ではなく、当初の全項目・全k・独立評価要件に状態を合わせる。

## 3. 共通の評価契約

### 3.1 対象・比較器・API

基準比較器は保存済み主系列と同じ **RDKit 2026.03.6** とし、実行時の配布物・版・
digestを照合して凍結する。2025.09.3、将来版、`@rdkit/rdkit` は別lane。
基準版の変更はmanifestを改版して両版を再測定し、暗黙に置換しない。

初期2D範囲は一般有機分子、塩、荷電・芳香族・縮合/架橋/大環状構造、明示H、
主要同位体、四面体・二重結合立体とする。具体的な元素/同位体表、電荷・環・分子サイズ境界、
各APIの対応表はA0.1で機械可読化する。金属配位、非四面体立体、Markush/ポリマーは
別の境界集合に残し、後から分母を狭めて成功させない。

`native`、`rdkit_compat`、規則準拠profile、`experimental` は別契約にする。
各操作について次を固定する:

- chematic APIとRDKit API、全引数・既定値、sanitize、H追加/除去、芳香族性、
  電荷/tautomer標準化、原子map、対称matchの重複処理。
- stable-keyの同一性契約とcanonical文字列表現の契約。
  文字列一致をグラフ/立体一致の代用にせず、InChIも単独の正解器にしない。
- FPはMorgan radius=2、2048 bitsを最初に採用。chirality、bond types、
  ring membership、count/bit、空FP規則を明記し、他の設定は個別採用。
- 制限時間、探索上限、最大入力サイズ、実行環境とseed。片側だけ前処理しない。

### 3.2 分母・誤差・出自

入力全体Nと有効な宣言範囲N_validを実行前に固定する。
範囲外・不正入力も境界集合として報告し、受理すべき有効入力の未対応はN_validに残す。
各エンジン・各操作の `ok / unsupported / invalid_input / error / timeout / not_measured`
を独立保存し、Nに対して完全に会計する。正当な拒否は安全性合格に数え、化学的正解には数えない。

| 指標 | 定義・ゲート |
|---|---|
| 対応率 | 正常な化学結果数/N_valid。全入力Nを分母にした値も別表示 |
| 比較可能集合一致率 | 両者成功した入力での一致率。これだけで全体同等としない |
| 全入力有用一致率 | 正常かつ一致した結果数/N_valid。拒否・timeoutで上げられない |
| 8記述子 | HBA/HBD/aromatic ring countは完全一致、MW/TPSA/LogP/MR/Fsp3は絶対差<=1e-6 |
| 中心・ラベル・match・FP | 原子対応後の集合/カテゴリ/整数/bit/ID/countが完全一致 |
| 検索 | 全クエリのIDと順序一致、recall=1.0、Tanimoto絶対差<=1e-12。順位は丸め前 |
| 重大誤り | 誤った確定立体、構造/立体消失、stable-key誤統合、無告知metadata損失は0 |

8記述子には従来の公開許容差での値も併記する。MAE、中央値、p95、最大差、
全不一致IDと構造クラス別集計を出す。bindingごとの欠測も集約で消さない。
native既定値はタグbaselineとcandidateを別成果物で比較し、同じbinaryの二重使用を拒否する。

各結果はrun ID、source commit、dirty差分と未追跡ソースのdigest、lockfile、
Rust binary/Python拡張/JS/WASMのdigest・読込パス、比較器版・配布物digest、
corpus SHA-256、分割、設定、seed、コマンド、日時、OS/CPUを持つ。
行ごとの生結果から集計を再計算し、要約JSONだけの自己申告で通さない。
欠測値は `null` と型付き理由とし、NaN/Inf・空集合・重複ID・欠落行はゲート失敗にする。

### 3.3 コーパス設計と未使用評価

| セット | 初期設計 | 利用方法 |
|---|---|---|
| 既存回帰 | 既存5,000件、12件MW、CIP155行、#149/#503、#337、標準化10件 | 開発で既に見たデータ。回帰用として維持 |
| 新規一般分子 | 10,000件: 開発2,000／封印評価8,000 | ライセンス・取得版・重複監査後に凍結。必要標本数に応じ50,000件まで拡張を検討 |
| 難例 | 300件以上: 開発200／封印評価100以上 | 同位体、ラジカル、P/S、共鳴、架橋/大環、塩、対称性、立体を層別 |
| A5 gold pilot | 200件を目安に規則・絶対期待値を整備 | 統計設計・判定手順の校正用。正式な同等性証明には使い回さない |
| 変形回帰 | 通常16、高リスク256、#503は1,024並べ替え | 分子数に加算せず、変換自体の元素・同位体・電荷・立体保存を先に検査 |

取得・分割・独立gold作成は今後の作業であり、上記件数の新データが既にあるとはしない。
親分子・塩・互変異性体を同じ群に置き、scaffoldで分割する。
scaffoldを持たない分子は親構造/seriesをcluster単位として事前定義する。
重複検出は開発対象のstable-keyだけに依存せず、第三者ツールとグラフ照合を併用する。

PR/nightlyは開発・回帰群のみ。封印評価はcandidate固定後に実行し、開示履歴を保存する。
修正に利用した評価ケースは次回から回帰群に移し、未使用分を補充する。
engine名の盲検化、入力の未使用性、正解の独立レビューはそれぞれ別に確認する。

## 4. 作業パッケージと出口

各項目は `planned → implemented → regression_passed → heldout_passed → adopted`
を区別する。A5にはさらに `independently_reviewed` が必要。
下記チェックは残作業。既存合格は第2節に集約し、狭い合格を全体完了に繰り上げない。

### A0 — 評価ゲートと基準固定

目的: 偽の100%を防ぎ、以後の変更を同じ条件で比較する。

- [x] **A0.1** `validation/manifests/rdkit_accuracy_v2.json`へ改版し、操作対応表・比較器pin・
  期待件数・許容差・split・artifact digestを必須化する。既存v1は履歴として保持。
  `scripts/check_rdkit_accuracy_manifest.py`が実在artifactのSHA-256とsealed evaluation / candidate
  freeze要件をfail-closedで検証する。通常実行は開発manifestを検証し、リリース判定では
  `--require-sealed`を付けて`status=sealed`、sealed split、凍結SHA-256、candidate commit/tagを
  必須化する。未凍結の8,000件評価はstatusで明示し、完了扱いにしない。
- [x] **A0.2** 記述子holdoutを8項目・厳密許容差・全bindingへ拡張。
  native固定4件に加え、変更した認識クラスをタグbaseline/candidateで比較する。
  8項目のRust/Python/Node-WASM binding一致は5,000行および7,737行exposed holdoutで合格済み
  （化学的正解判定とは別ゲート）。
- [x] **A0.3** stereocenter runnerの総不一致数、例外、表示サンプルを分離。
  expected row count、読込binary、versionを実測し、中心数に加えて原子別結果を保存。
- [x] **A0.4** ゲート負例を自動化。HBAのみ誤り、厳密TPSA差、NaN/Inf、
  boolを整数として渡す入力、空/全unsupported、50件超の不一致、行欠落/重複、
  stale binary、異なるcorpus/版、baseline使い回し、summaryとrawの不整合を必ず拒否。
  現在はvalidatorがbool件数、空の成功集合、manifest外フィールド、欠落したchematic
  provenance、重複/不正ID、非有限値を拒否し、focused negative testsで固定している。
  manifest validatorはcorpusの存在、論理行数、SHA-256も検証し、改変・差し替えを拒否する
  回帰テストを追加した。7,737行exposed binding summaryもv2 manifestへ登録した。
  artifactは空ファイルを拒否し、比較環境の失敗で生成された空JSON/JSONLが
  測定済み成果物として受理されない負例を回帰固定した。
  `check_descriptor_raw_accounting.py`はbaseline/candidateのschema-v2 `raw_rows`から
  field別parsed/matches/strict_matches/mismatchesを再計算し、stale summaryとraw hash不一致を拒否する。
  baseline/candidate packet validatorはsource/artifact identityの再利用とcomparator/corpusの
  不一致を拒否する。現行タグv1.0.13は互換8記述子APIを持たないため、native値をbaselineに
  代用せず、互換API導入コミット`5aae7b23`をbaselineとして明示した。baseline/candidateは
  同一corpus・RDKit版で比較し、field別coverage差は未対応を隠さず改善結果として保存する。
  追加した非有限値、raw行欠落、比較器/コーパス不一致、51件超の不一致を含むfocused negative回帰は通過済み。
  stale binary、raw/summary一致、50件超の不一致はA0.5の実測packetへ接続する。
- [x] **A0.5** 互換API導入コミット`5aae7b23`とcandidateを別ビルド・同じ環境で再測定し、
  `baseline / candidate / delta / N / coverage / mismatch / provenance` のscorecardを作る。
  保存済み結果の履歴は上書きしない。タグv1.0.13は互換API導入前なのでbaselineには使わず、
  `validation/results/descriptor-baseline-candidate-v1.0.13-measured.json`にnot-measured理由を
  混同せず、互換API導入時点の実測を保存した。

作業先: `scripts/check_rdkit_accuracy_gate.py`、`check_descriptor_promotion_gate.py`、
`check_descriptor_rdkit_holdout.py`、`descriptor_stereocenter_rdkit_parity.py`、
`check_native_descriptor_regression.py`、`validation/manifests/`。

**出口:** 全negative fixtureが期待どおり失敗し、正例は通る。全行と成果物を追跡でき、
baseline/candidate比較が再現できる。A0は評価器の完成であり、化学的不一致を隠して通すものではない。

### A1 — 共通認識・記述子・立体中心

- [x] **A1.R（受入れ確認、P1/P2）** bounded chordless-aromatic-cycle候補で芳香族環数の退行を解消した。
  fixes3/fixes4/fixes6の全行差分で新規・継続・解消した不一致を分け、反例を最小化する。
  現候補は固定RDKit 2025.09.3に対して8項目すべて7,737/7,737 strictで、
  promotionは`adopted_opt_in`。単体テスト合格だけでは採用しない。
  native既定動作を保ち、共有認識に依存する立体中心・A2/A3/A4/A6の関連fixtureと
  Rust/Python/Node/WASMゲートを最新ビルドで再検証した。perception 204実行+1 ignored、
  fingerprint 309、SMARTS 185の回帰、およびpromotion gate `adopted_opt_in`が通過済み。
  fixes3の中心集合合格は
  fixes4以降の証拠として流用しない。影響範囲のbinding/workflow再検証と、この公開済み集合とは別の未使用評価が必要。
- [ ] **A1.1** MWとexact massを分け、元素/同位体表、D/T/13C/15N/18O、
  明示・暗黙H、価数・電荷・ラジカルを原子別に診断。表にない同位体の質量数代用は禁止。
  HBA/HBD/TPSA、Crippen LogP/MRは原子/fragment寄与で追跡する。既知同位体は
  nuclide massを使い、未知同位体は`NaN`へfail-closedする実装と回帰を追加した。
  5員環の芳香族酸素架橋ではmorphine/codeineのRDKit LogP残差を解消し、6員環の
  単純環状エーテルを誤昇格させない負例を維持した。descriptor/ChEMBL各5,000行の
  再測定では、P-H phosphonate、thioamide N、charged sulfurのTPSA環境修正後、
  8項目すべて5,000/5,000 strict一致した。証拠は
  `validation/results/descriptor-rdkit-a1.1-tpsa-fixed-v1.0.13.json`。
  未使用評価と広い原子別corpusは未完了。
- [x] **A1.2** 残っていたpotential-center 3分子を原子対応で最小化した。
  中性sp3 NのRDKit公式環拘束条件に、全隣接枝の環内次数2・非アミド/非スルホンアミド・
  provenance-preserving branch比較を適用し、5,000行でFP=0/FN=0を確認した。
  結果は `validation/results/descriptor-stereocenter-rdkit-parity-v1.0.13-source-2025.09.3-atom-level-after-ring-n-official-rule.json`。
  中心数の一致だけでFPとFNを相殺せず、存在・指定/未指定・元素/結合環境・環制約を別列で保持する。
- [ ] **A1.3** 有限深度のbranch比較と中心判定を分離し、対称性/環の同値判定を整備。
  現行深度16をさらに増やすだけで一般解とはしない。芳香族S例外を電荷・配位数・結合環境で
  検証し、芳香族N/P・非対象Sの負例も置く。探索予算超過は未解決とし、確定判定を捏造しない。
  現在はcarbon branch provenanceと中性sp3 ring-Nの狭い条件を実装済みだが、巨大縮合環の
  芳香族環基底差による巨大ケージ状分子1件が残るため、A1.3全体は未完了。
- [ ] **A1.4** core8項目、Strict rotatable bonds、potential centersを既存＋新規群で検証。
  7,737-row exposed holdoutでは8項目中7項目とpotential centersがstrict合格し、中心集合は
  FP/FN=0、RDKit側の5件の明示的fallbackを含め未解決0だった（fixes3時点）。混在芳香族性の広域修正はfixes4で芳香族環数を1残差から29残差へ退行させたため戻し、fixes6で1残差へ回復した。他7項目は全件strict合格を維持。A1.Rの残差修正と中心集合の再検証、
  Strict rotatable bondsの同じ未使用群での再測定、Rust/Python/Node/WASM全bindingの拡張評価は残る。
  Rust/Python/Node/WASMに存在しない項目は公開APIを整備するか未測定として残す。
  共有認識変更後の現行wheelによる7,737-row再測定は8項目すべてstrict 7,737/7,737、
  parse失敗0・unsupported0だった。証拠は`validation/results/descriptor-rdkit-a1.4-recheck-v1.0.13.json`。
- [ ] **A1.5** 拡張familyをAPI棚卸し→既定値/定義→寄与診断→独立oracle→bindingの順で採用。
  順番は追加環指標 → VSA/MQN → Kappa/Hall–Kier/Bertz/Balaban/BCUT2D。
  familyごとに整数/浮動小数の許容差と対応範囲を結果を見る前に固定する。

作業先: `chematic-core`、`chematic-perception`、`chematic-chem`、
`scripts/descriptor_core_parity.py`、`descriptor_rdkit_diagnostics.py`、
`descriptor_rotatable_rdkit_parity.py`、`descriptor_stereocenter_rdkit_parity.py`、各binding。
共有認識を修正したらA2/A3/A4/A6の関連fixtureも同時実行する。

**出口:** 宣言した各familyで不一致0・coverage100%・native未承認変更0。
potential centersは原子別FP=FN=0も必要。core合格と全family完了は別状態として記録する。
A1.5はA2以降と並行でき、未完なら「追加記述子は未完」と明示する。

### A2 — CIP・canonical・構造同一性

RDKitもpotential centerの検出と絶対配置の割当を分けているため、
legacy、`FindPotentialStereo`、modern `rdCIPLabeler` を別列で評価する。
[RDKit Book: Stereochemistry](https://www.rdkit.org/docs/RDKit_Book.html#stereochemistry)

- [x] **A2.0** `num_unspecified_stereocenters`のfalse positiveを修正した。
  旧実装は4価炭素だけを見ていたため、普通のCH2/CH3まで未指定立体中心に数えていた。
  A1のRule-5-aware potential-center集合を再利用し、`C`、`CC`、isobutaneは0、
  `C(F)(Cl)Br`は1となる回帰を追加した。
- [ ] **A2.1** 既存CIP全コーパスを基準版で再測定し、P系15件、負電荷共鳴、
  同位体、擬不斉、duplicate-nodeを規則別に分類。R/S、r/s、E/Z、非中心、
  指定なし、未解決を原子/結合mapで区別する。JSONLの残差分類コーパス181件は
  固定RDKit 2025.09.3でpotential-centerの件数不一致0、parse失敗0を再確認したが、
  これはCIPラベル割当の正解証明ではない。atom R/Sとbond E/Zを別比較する
  `scripts/cip_label_parity_audit.py`も追加し、181行ではatom map 155/181、
  E/Z bond map 14/14、未解決理由52件（`oracle_unstable`/`tied`）を記録した。
  atom側の26不一致はすべてP中心かつ`oracle_unstable`で、確定ラベル同士の不一致は0件。
  証拠は `validation/results/cip-residual-stereocenter-v1.0.13.json` と
  `validation/results/cip-label-parity-a2-v1.0.13.json`。
- [ ] **A2.2** 最小反例に正負の近傍構造を追加し、中心検出とラベル割当のどちらが原因か診断。
  規則準拠の修正はA5 goldに接続する。RDKitが返すことだけを正解の根拠にしない。
- [x] **A2.3** #149/#503を最新ソースでK=1,024まで再測定。
  現行ソースで28 coupled components中4件の残差と0件の対応関係破綻を再現し、
  source commit・dirty状態・corpus hashを要約へ保存した。これは残差の再現ゲートであり、
  修正完了やcanonical winnerの証明ではない。
  小分子の全列挙、共有E/Z carrier、aromatic stash、ring closureのfixtureで修正を判定する。
- [ ] **A2.4** parse→write→parseの原子/結合/立体保存、再出力の冪等性、
  atom-renumbering不変性を検証。対掌体、E/Z、同位体、電荷、結合異性体の負例で
  stable-key誤統合を調べる。RDKit InChI identity oracleによる5,000分子×4 randomized
  valid spellingの構造同一性は失敗0で通過した。さらに、現行公開コーパス3種の
  14,999行でstandardize→canonicalizeの冪等性が14,999/14,999通過した
  （`validation/results/canonical-idempotency-a2.4-current-v1.0.13.json`）。これは
  回帰証拠であり、未使用データによる独立正確性やRDKitの文字列一致を意味しない。
  写像はatom-mapを答えとして利用する経路と分離する。
- [ ] **A2.5** 適用範囲を広げた不変性・衝突・立体goldを全bindingで実行。
  誤った確定結果を避ける既存fail-closed境界は、代替実装の合格まで維持する。
  代表的な表現不安定リン酸ケースについては、PythonとNode/WASMの両方で
  accurate assignmentを空にし、`oracle_unstable`/`oracleUnstable`を返す回帰を追加済み。

作業先: `chematic-chem/src/cip.rs`、`chematic-cip`、`chematic-smiles`、
`scripts/cip_accurate_full_corpus_report.py`、`canonical_cross_engine.py`、
`canonical_structural_correctness.py`、`check_ez_residual_evidence.py`。

**出口:** 対応範囲内の中心・ラベル不一致0、coverage100%、並べ替え差0、
構造/立体損失0、stable-key誤統合0。P系の拒否やK=1,024残差が残ればA2未完。
canonicalのcross-engine文字列の綴り一致は必須にしない。

### A3 — Morgan・検索の完全性

- [ ] **A3.1** FPを atom invariants → environment → raw ID → duplicate suppression →
  sparse count → folding → bitInfoの順で独立RDKit oracleと比較。
  既存5,000件に新規群を加え、全入力IDとエンジン別状態を保つ。
- [ ] **A3.2** 同一ライブラリ・クエリに対する総当たり参照を用意。
  bit Tanimotoはintersection/unionの整数比を参照し、交差積で同点と順位を判定する。
  ゼロunionの規則を先に固定する。丸めは表示時のみ。
- [ ] **A3.3** k=1/10/100、threshold直上/同値/直下、chunk境界、近接スコア、
  同点・重複分子・空FP・空DB・k>件数・失敗行・取消をfixture化。
  同点は元入力index昇順。chunk再結合がglobal brute-forceと一致することを検査する。
  閾値APIは実装済みで、`score >= threshold`、有限値かつ `[0,1]` の入力契約を
  Rust/Python/WASMで共有する。RDKit由来96ケースのcross-binding実測は合格済みである。
- [ ] **A3.4** 同じ生成成果物でRust、source-built Python、Node/WASMを実行。
  binding一致と独立RDKit一致を別ゲートにし、検索ベンチマークのrecall>=0.99を
  このexact profileでは全クエリ1.0かつ順位一致へ強化する。
  WASMの `RdkitSearchIndex.search_json` 自体も6桁へ切り捨てているため、
  runnerだけでなくfull-precisionの型付き結果またはJSON経路を整備する。
  既存の丸め契約は明示して保ち、新しいexact経路をopt-inで公開・検証する。
- [ ] **A3.5** 既存4,500/500の除外0を回帰条件に固定し、封印評価用ライブラリ/クエリを
  別途凍結。native対nativeも同じID集合で一致を確認し、native対Morganは診断列に残す。

検索oracle、threshold、cross-bindingの各レポートは、source commit、dirty差分と
未追跡ソースのdigest、Cargo.lock、実行Pythonとプラットフォームを共通provenanceとして
保存する。比較実行中に既存レポートを切り詰めない原子的書き込みも共通契約とする。
`check_rdkit_search_reports.py`は3レポートが同一corpus・同一source/diff provenanceを
共有し、すべての比較階層でmismatch=0であることを横断検証する。

作業先: `chematic-fp`、`scripts/ecfp_rdkit_*_parity.py`、
`bench_similarity_search_vs_rdkit.py`、`rdkit_search_cross_binding_parity.py`、
`rdkit_ecfp4_*_cross_binding_parity.py`、各binding。

**出口:** 宣言設定で全bit/count/由来情報一致、coverage100%、
全クエリの各k・thresholdでID/順位一致、recall=1.0、スコア差<=1e-12。
将来の近似検索は別profileとし、このexact契約を緩めない。

### A4 — 標準化・SMARTS・反応・V3000

- [ ] **A4.1** fragment selection、中和、tautomerを個別API契約で再測定。
  履歴の標準化10件中酢酸塩2件を、現在のfragment_parent/charge_parentの規則に照らして判定。
  200件以上の層別fixture、冪等性、原子/電荷収支、塩の順序不変性を用意する。
- [ ] **A4.2** SMARTSの18拒否・15残差を環モデル/芳香族性/写像に分類。
  既存31 patternの全セルを維持し、100 pattern以上へ段階拡張。
  存在判定から原子対応match集合のFP/FNへ拡張し、recursive/立体/明示H/対称重複、
  `[R]`/`[r]`/`[k]`、探索打切りを検査する。入力順依存の特例で一致させない。
- [ ] **A4.3** 反応presence 6件を回帰として保持し、50 template×各4 reactant条件を
  初期設計とする生成物評価を追加。正例・無反応・競合部位・立体条件を含め、
  product集合、map、価数、原子/電荷収支、立体保持/反転、複数生成物、失敗を評価する。
  同じ規則templateの比較と実験反応予測精度は別主張にする。無反応・荷電変換・分解を含む
  拡張51件（14 template、無反応・荷電変換・分解・複数条件を含む）ではRDKitとchematicの
  atom-map付きcanonical product setが51/51一致した。RDKitのembedding重複は集合比較で除去する。
  立体不一致の安全側拒否は互換性分母から分離し、専用ゲート1/1で確認した。
  証拠は `validation/results/reaction-product-parity-v1.0.13-expanded.json` と
  `validation/results/reaction-product-stereo-safety-v1.0.13.json`。
- [ ] **A4.4** V3000はSGROUP/COLLECTION、同位体、enhanced stereo、配位、
  ENDPTS/ATTACH、query/polymer境界をfixture matrix化。
  chematic→RDKit、RDKit→chematic、Indigoとの両方向で、原子/結合/typed metadataを比較。
  最初に各対応機能の正例・負例・順序変更を最低1組ずつ置く。現行実装では
  SGROUP/COLLECTION、同位体、enhanced stereoの限定範囲に加え、V3000結合行の
  `ENDPTS=`/`ATTACH=`等をopaque属性として保存・再出力し、専用回帰が2/2通過した。
  意味解釈や両方向の外部エンジン一致はまだ主張しない。
- [ ] **A4.5** opaque保存、型付き解釈、編集/展開可能性を別列にし、
  単に再読込SMILESが同じでもmetadata損失を合格させない。
  無制限Markush/ポリマー展開は対象外だが、保存・明示拒否の契約は検証する。

作業先: `chematic-smarts`、`chematic-rxn`、`chematic-mol`、
`scripts/standardization_holdout.py`、`rdkit_ring_parity_diagnosis.py`、
`reaction_rdkit_quality_gate.py`、`crates/chematic-mol/tests/v3000_semantic_interop.rs`。

**出口:** 宣言操作の全有効入力でcoverage100%、正解match/product集合FP=FN=0、
無告知の構造/立体/metadata損失0。未導入Indigoは `local-toolchain`、
利用許諾・外部fixture提供待ちは `external` とし、片側のみで両方向完了としない。

### A5 — 独立正解で同等・優位を判定

最初の対象はCIPと構造保持。規則の根拠は
[IUPAC Blue Book P-9](https://iupac.qmul.ac.uk/BlueBook/P9.html) 等の版・節をケース別に記録する。
RDKit側のアルゴリズム由来は
[rdCIPLabeler documentation](https://www.rdkit.org/docs/source/rdkit.Chem.rdCIPLabeler.html)
で別に追跡し、規則の改訂提案と採用規則を混同しない。

- [ ] **A5.1 local-open:** gold候補を `source / rule_version / input_digest / atom_mapping /
  expected_absolute / rationale / reviewer / unresolved_reason` のケースにする。
  不変性だけのケースはmetamorphic regressionへ分類し、絶対ラベルgoldと区別する。
- [ ] **A5.2 local-open:** 200件pilotで通常例・負例・難例を事前選定し、
  手順、cluster、endpoint、必要標本数を設計する。公開済みblind欄は開発例に戻す。
  正式な封印集合は別に用意し、入力・正解・両engine出力の開示履歴を残す。
- [x] **A5.3 local-open:** paired集計・不確実性の評価器を実装し、
  全一致、両者全失敗、片側拒否、同一scaffold多数、極小標本で負例検証する。
  `scripts/evaluate_a5_paired.py`が全行を保持したカテゴリ集計とcluster単位の
  10,000回paired bootstrap（seedを記録）を生成する。未解決行は分母の監査に残し、
  正解には数えない。独立レビュー済み入力が未取得のため、正式な判定packet生成は
  外部レビュー後に行う。
- [ ] **A5.4 external:** 非実装者が絶対期待値をレビューし、不一致を規則に基づいて裁定。
  第三エンジンは差の発見補助とし、多数決を正解の定義にはしない。
- [ ] **A5.5:** candidate固定後に正式集合で実行し、独立レビュー済み正解に照らして公開する。
  規則準拠profileと版固定互換profileが異なる場合は両方の結果を示す。

主指標は、事前に妥当性と正解を確定した有効入力全体N_validに対する
「正しい結果を返した割合」。拒否・timeout・誤答は成功0として残す。
各入力で両者正解／chematicのみ正解／RDKitのみ正解／両者非正解を集計し、
各engineの拒否等も別軸で保存する。正解自体が未解決のケースは黙って除外せず、
N全体の内訳と最良/最悪の感度分析を出す。結論が変わる場合は未証明とする。

scaffold/series単位のpaired cluster bootstrap（10,000反復）で95%信頼区間を出す。
少数・全一致で幅0となる場合は同等性の証明に使わず、事前指定した保守的手法を
採用するか証拠不足とする。pilotから必要cluster数を設計し、±0.1 ppに不足すれば
正式データを増やす。試行後の任意停止・閾値変更で結論を作らない。

| 結論 | 事前条件（差=chematic−RDKit、pp=percentage point） |
|---|---|
| 同等 | 正解率差の95%区間全体が[-0.1, +0.1] pp内、対応率の非劣性、重大誤り回帰0 |
| 非劣性 | 正解率差の95%下限>-0.1 pp、対応率も同じ非劣性条件、重大誤り回帰0 |
| 優位 | 正解率差の95%下限>0、対応率の点推定は低下せず非劣性条件も満たす、重大誤り回帰0 |
| 限定改善／未証明 | 上記未達、標本不足、未解決gold、または選別したchallenge-setのみ |

主endpointはCIPと構造保持の全条件を満たした入力割合とし、領域別結果は副指標にする。
複数の優位性主張を行う場合はfamily-wiseな多重比較調整を事前登録する。
一般分子の評価と難例評価は別集計とし、難例を多く集めた比率を一般分子の正解率にしない。
RDKitの既知誤りだけを集めた集合の勝率を一般分子での優位性に外挿しない。

**出口:** 凍結protocol・レビュー済みgold・生結果・不確実性・再現手順が揃い、
同等または優位の条件を満たす。評価を実施しても非劣性や証拠不足のみなら
A5の「同等以上」目標は未達として残す。ローカル準備は外部レビュー待ちの間も進められる。

LogP/MRのRDKit一致は実測物性の予測精度ではない。
実測LogP/活性検索に進む場合は測定条件、ライセンス、学習重複、scaffold/target分割を
定める別計画にし、A0–A6の達成条件には追加しない。

### A6 — 3D・力場の別プロファイル

- [ ] **A6.1** #337の環代表集合→MMFF atom typing→BCI/charge→parameter coverageを診断。
  MMFF atom typingは6,681/6,698 exact（宣言済みunsupported probe 1件を除く比較対象では
  6,681/6,697、99.76%）まで改善した。今回の構造境界で16残差を解消し、残る比較対象残差は
  0009/0029/0030の3 macrocycle群。結果は
  `validation/results/mmff94-type-parity-a6-macrocycle-v1.0.13.json`に保存する。
  最後のunsupported 1件を集計から落とさない。MMFF94/94s/UFFとfallbackを別列にする。
- [ ] **A6.2** 同一座標・同一H/電荷・同一parameter variantでbond、angle、stretch-bend、
  out-of-plane、torsion、vdW、electrostaticを項別および総和で比較。
  finite differenceによる独立勾配検査も実行する。現行ソースのbounded strict bond+angle
  再測定は265/265（tier A 65、tier B 200、失敗0）で通過したが、これは全項目・全MMFF94
  パリティの証拠ではない。結果は
  `validation/results/mmff94-strict-gate-current-candidate-v1.0.13.json`。同じ現行候補の
  BCI電荷は6,665/6,698 exact、unsupported probeを除く比較対象では6,665/6,693
  （99.58%）で、28原子・3マクロサイクルが残る。結果は
  `validation/results/mmff94-charge-parity-a6-macrocycle-v1.0.13.json`。現行Rustのtier-A
  pipeline（65件）は有限座標・sound・衝突なしが各65/65だが、force-field収束は42/65で、
  23件が200反復上限に達した。これは収束非劣性をまだ示さない候補証拠であり、
  `validation/results/mmff94-pipeline-tier-a-current-candidate-v1.0.13.json`に保存する。
  行単位の外部ハードタイムアウトでA/B全265件を再測定した結果、215件がsuccess、
  50件が明示的timeoutとなり、success行はすべてfinite/sound・gross clashなしだった。
  200反復での収束は54/215 success（A 42/64、B 12/151）。この全行coverage証拠は
  `validation/results/mmff94-hard-timeout-pipeline-a6-v1.0.13.json`に保存し、
  timeoutを成功から除外していない。なお、同一H・同一座標でのRDKit全項目parityは未達。
  診断用に反復上限500も測定したところ、完了行の収束は57/64まで改善した一方、
  atorvastatin_fragmentがtimeoutとなった。budgetを増やすだけでは出口にならず、
  200回既定値は変更しない。RDKitが生成した同一明示H座標を両エンジンへ渡す追加診断では、
  265件中262件が比較可能、2件が埋め込み失敗、1件がschematic非対応となり、絶対エネルギー差は
  中央値101.39、p90 185.71、最大387.69 kcal/molだった。これは配座生成差ではなく、
  パラメータ／項の適用差を示す候補診断であり、A6の合格証拠ではない。結果は
  `validation/results/mmff94-same-explicit-h-energy-a6-265-v1.0.13.json`、runnerは
  `scripts/mmff94_same_explicit_h_energy.py`。この診断でBuffered 14-7の実装誤り
  （`t^7 * (t^7 - 2)`）を特定し、RDKitの
  `t^7 * (1.12R*^7 / (r^7 + 0.12R*^7) - 2)`へ修正した。修正後は262/265件が比較可能で、
  中央値0.262、p90 9.297、最大32.67 kcal/mol、5 kcal/mol以内200/262となった。
  結果は`validation/results/mmff94-same-explicit-h-energy-a6-265-after-vdw-v1.0.13.json`。
  残るmacrocycle/stress残差、勾配、収束、立体、配座品質は未完了。比較は
  `validation/results/mmff94-convergence-sweep-a6-tier-a-v1.0.13.json`。最大残差勾配を
  追加取得した結果、収束群42件の最大は9.998e-5 kcal/mol/Å、非収束群23件の中央値は
  0.0718、最大は7.04 kcal/mol/Åで、単純な反復上限だけでなく局所的な数値／項別問題を
  切り分ける必要がある。結果は
  `validation/results/mmff94-convergence-residual-gradient-a6-tier-a-v1.0.13.json`。
  項別エネルギーを追加取得し、atorvastatin_fragmentでは最終bond 24.75、vdW 85.41、
  torsion 23.91 kcal/mol、最大残差7.04、cholesterolでは最大残差0.0638を確認した。
  これは同一MMFF94関数内の診断であり、RDKitとの項別一致を意味しない。結果は
  `validation/results/mmff94-energy-breakdown-a6-residual-pair-v1.0.13.json`。
  項別勾配では、cholesterolがbond/angle 19.82、nonbonded 22.33、stretch-bend 11.78、
  atorvastatin_fragmentがbond/angle 38.06、nonbonded 37.99、stretch-bend 11.45
  kcal/mol/Åであり、総勾配との相殺が確認された。項別・総和の両方をゲートし、閾値を
  緩めて隠さない。結果は
  `validation/results/mmff94-term-gradient-a6-residual-pair-v1.0.13.json`。
  同一重原子座標をRDKit 2025.09.3へ入力した診断では、総エネルギー差がcholesterol
  14.77、atorvastatin_fragment 51.77 kcal/molだった。ただしimplicit Hと項の規約を
  固定していないため、最終parity値とは扱わない。fixed-H protocolを先に確立する。
  結果は `validation/results/mmff94-rdkit-same-heavy-coordinates-a6-residual-pair-v1.0.13.json`。
- [ ] **A6.3** 小分子/芳香族/荷電/架橋/大環などの層別回帰で、
  最適化の収束、結合/非結合衝突、立体保持、非有限値、打切りを別判定。
  productionの有限差分経路を置換する前に同じ範囲を解析勾配ゲートで通す。
- [ ] **A6.4** 配座評価は新規300分子を初期設計（100開発/200封印）とし、
  seed群0–9・各20配座を初期候補に計算予算を開発pilotで凍結する。
  ETKDGの版/全設定と同じ分子・seed群・予算を使い、seedを共通にしても同一乱数列とは扱わない。
- [ ] **A6.5** 対称性を考慮したRMSDとTFDを同じ固定参照に対して評価。
  RDKit生成配座への近さ、実験構造への近さ、同じ力場上のエネルギーは別列にする。
  実験構造は相・測定条件・利用権を確認し、なければ実験精度は `not_measured` とする。

数値ゲートの初期案（達成値ではない）:

| 量 | 初期許容差・合格条件 |
|---|---|
| atom type / parameter有無 | 宣言した有効範囲で100%一致、隠れたfallback0 |
| 部分電荷 | 原子ごと絶対差<=1e-6 e、全電荷収支一致 |
| 同一座標エネルギー | 項別/総和で絶対差<=1e-4 kcal/mol + 1e-8×参照絶対値 |
| 解析勾配 | 最大成分差<=1e-4 kcal/(mol Å) + 1e-5×参照成分絶対値。特異点を分離しFD刻み1e-4/1e-5/1e-6 Åで安定性確認 |
| 有用配座成功率 | valid geometryかつ立体保持の全入力成功率差の95%下限>-1 pp。重大な立体破壊0 |
| 配座品質 | 共通参照に対する分子別best RMSD/TFDの悪化量p95の95%上限<=0.10 Å / 0.02を候補とする |

上記はプロジェクトの採用基準案であり、RDKit/IUPACの規格値ではない。
比較器の反復誤差・倍精度数値誤差・化学クラスを開発pilotで評価してから凍結し、
変更理由を記録する。封印評価後の緩和は禁止。必要なprecisionを公開APIが出せない場合は
その数値profileを未採用とする。

作業先: `chematic-ff`、`chematic-3d`、既存MMFF94型/項別/勾配runner。
**出口:** クラス別の全項目ゲート、全入力の有用配座成功率と配座品質の非劣性を満たす。
不足クラスは明示し、A6全出口まではExperimentalを維持する。

## 5. 実装順・候補版・運用

### 5.1 最初の実装単位

| 順番 | 変更単位 | 受入証拠 |
|---|---|---|
| 1 | A0.1–A0.4: schema・件数/出自・全項目holdout・negative tests | 正しい記録だけ通り、全破損fixtureが失敗する |
| 2 | A0.5 + A1.2: baseline再測定・5分子の原子別最小反例 | 同条件before/afterと中心FP/FN、原因と未解決事項 |
| 3 | A1.1–A1.4: 共通認識の修正とcore再検証 | 既存/新規開発群の不一致0、native・関連領域回帰0 |
| 4 | A2.1–A2.5: CIP/identity残差と広い不変性 | ラベル・coverage・損失・誤統合の4軸、K=1,024 |
| 5 | A3.1–A3.5: 全k/threshold/独立oracle | 全クエリID・順位・丸め前値・binding一致 |
| 6 | A4.1–A4.5: 標準化→SMARTS→反応生成物→V3000 | 操作別の全分母、FP/FN、metadata保持 |
| 並行 | A1.5、A5.1–A5.3 | 追加family別ゲート、gold/統計評価器/レビューpacket |
| 後続 | A5.4–A5.5、A6 | 独立判定、別3D候補の証拠 |

共通認識の修正は下流のCIP、Morgan、SMARTS、MMFFに波及する。
小さい変更単位で検証し、1件の改善と同時に別クラスを悪化させた変更は採用しない。
未使用データの取得や外部レビューが待ちになっても、既存fixture・評価器・生結果監査を進める。

初期工数目安はA0 2–4、A1 5–10、A2 8–15、A3 3–6、A4 7–14、
A5ローカル5–10、A6 15–30実働日以上。残差の研究・データ取得・独立レビュー次第で変わる
作業規模の見積りであり、納期や自律実行時間の約束ではない。

### 5.2 候補版の出口

| 区切り | 必須条件 | 主張 |
|---|---|---|
| 次の候補 | A0、A1 core8項目の新規評価、A2/A3/A4既存安全回帰、影響binding、native既定値維持 | 改善が確認できた具体的操作。CIP等が残れば全2D同等とはしない |
| 2D互換候補 | A0、宣言したA1 family、A2/A3/A4が封印評価まで全出口合格 | 指定版・2D profile・評価セットでRDKit互換 |
| 独立精度候補 | A5独立レビューと同等または優位条件 | 実証した領域でのみ「同等」「より正確」 |
| 3D候補 | A6の数値・coverage・配座品質ゲート | 指定した力場/構造クラスの3D品質 |

RC番号・公開日・commit/tag/push/publishはこの計画で実行しない。
A1追加family・A5・A6の未達を、次候補へ入らないという理由で完了にしない。

### 5.3 継続検証と成果物

- **PR:** 最小反例と近傍負例、既存回帰、評価器negative tests、native既定値、
  型付き失敗、影響bindingの実ランタイム。認識変更は関連下流fixtureも必須。
- **Nightly:** 固定した開発/回帰群、難例・並べ替え、全binding、基準RDKit。
  重い3D/全列挙はbounded jobに分割し、timeout・未実行を成功扱いしない。
- **Candidate:** 別ビルドのタグbaseline/candidate、封印holdout、全入力会計、
  artifact digest、群別誤差、未対応、native差分を確認する。
- **公開文書:** `docs/validation.md`、`docs/rdkit-comparison.md`、dashboard、
  CHANGELOGに出典付き結果を集約。READMEは短い要約とリンクを維持する。

A0で次schemaのmanifest、A1–A4/A6で操作別raw JSONLと集計、
A5で判定packet・gold・paired report、最後にbefore/after scorecardを生成する。
これらは予定成果物であり、既存v1 manifestやJSONを全要件合格と読み替えない。
進捗は実装・回帰・封印評価・独立レビューを分け、残差IDと次の作業を常に残す。
