# chematic 1.x Trust Release 実行計画

更新日: 2026-09-26。リリース対象は **v1.0.27** です。次の開発では、配布channelの
独立検証と、残るCIP・SMARTS・A6 gateを優先します。

この文書は実行順と合格条件だけを定義します。機能別の優先順位は
[`ROADMAP.md`](https://github.com/kent-tokyo/chematic/blob/main/ROADMAP.md)、未完了項目と依存関係は
[`roadmap-open-work.md`](roadmap-open-work.md)、実測値は
[`validation.md`](validation.md) と
[benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks) を参照してください。

## 目的

Trust Releaseの目的は、機能数を増やすことではありません。利用者が次を確認できる
状態を作ることです。

1. どのAPIが安定・実験的・非対応か。
2. どの版、コーパス、設定、失敗方針で比較したか。
3. Rust/Python/Node/WASMで入力件数と結果が失われないか。
4. 壊れた入力や不確実な化学を、推測せず型付きで拒否できるか。
5. source候補、公開package、公開channelの証拠を混同していないか。

## 現在地

- A0 core-eight descriptor gateは完了しています。2,000件のdevelopmentと
  一度だけ評価した8,000件のsealed holdoutで、8項目すべてが合格しました。
  使用済み入力は露出済みで、今後のholdoutには再利用しません。
- 公開v1.0.20は、固定ブラウザ条件のParse + compatible Morganで
  RDKit.jsより高速かつ、対応9,999件でbit-exactです。
- 公開v1.0.20のMMFF94 stereo-safe laneは265/265件で品質条件を満たしますが、
  RDKitより高速ではありません。current sourceの高速化候補は別証拠です。
- v1.0.21に含まれる#632修正は、RDKit 2026.03.6比較のSMILES semantic差を
  18/10,000から0/10,000へ減らしました。28 component x 1,024 relabelの
  長時間gateは、v1.0.25ベースのclean source（`13d70a2e`）で再実行し0/28 divergentで完走しました。
- v1.0.23はRDKit-defined atom-pair/torsion、MACCS、QED、TPSA、Murckoの一致を
  改善しました。3本の5,000行source laneと43/44操作のmatrixは共有2 vCPU実行であり、
  公開packageまたはWASMの性能・完全互換性の証拠ではありません。
- v1.0.24はその出力を保持したままperceptionとSMARTS存在判定のhot pathを短縮しました。
  v1.0.23との差分は7コーパス・1,402,080出力行で0件です。計時は共有2 vCPUの
  source実行だけであり、公開package、WASM、全環境の性能主張ではありません。
- v1.0.25はatom-output order、断片とRust反応生成物のsource-atom provenance、
  RDKit互換Python HBA profileを追加・修正しました。6操作・46,736分子の既存出力
  differentialは0件ですが、HBA profileとplain writerの意図した修正は別に扱います。
- v1.0.26 release sourceでは、同じ固定laneでCIPが9,994/10,000一致（旧9,770）、
  SMARTS差分が200/310,000セル（旧14,306）です。残差はすべて原因別に分類済みです
  （#634: P oracle不安定4、三価N非対応1、要裁定1。#635: `[Rn]`/`[kn]`の
  ring数の意味194、フェロセン6）。
- v1.0.26 release source（#637）では、同一座標のMMFF94をRDKitの項別energyと
  比べ、262/262行が1 kcal/mol以内（最大0.32、旧9.87）です。
- 残る主要差分は、CIPの要裁定1中心、SMARTSのring数の意味、A6のheavy-atom
  typing（Kekulé/荷電入力の芳香族性）・timeout・conformer qualityです。

## 実行順

### 1. #632 SMILES E/Zを閉じる

状態: 下記の合格条件はv1.0.25ベースのsourceで満たしました（長時間gate 0/28 divergent）。

合格条件:

- focused parser/writer/canonical regressionが通る。
- RDKit 2026.03.6の固定10,000件でgraph差0、semantic差0を維持する。
- 28 component x 1,024 relabelの長時間gateを完走する。完走しない場合は、
  未実行であることをrelease判断に明記する。
- 測定外のcoupled shapeではstable-keyがfail-closedのままである。

### 2. #634 CIP差分を分類して解く

状態: source候補`11a4ea27`で230→6。比較器のE/Z結合同定を二重結合の端点に修正し、
`CipMode.ACCURATE`のE/Zを階層digraphの順位付けに切り替えました。
残る6件は分類済みです（`validation/results/rdkit-rebaseline-residual-classification-v1.0.25-issue634-635-vs-2026.03.6-2026-09-25.json`）。

差分を次の4種類に分けます。

- atom/bond correspondenceまたは比較器の問題;
- RDKit版・表現依存などoracleの問題;
- chematicが明示的に非対応とする領域;
- chematic実装の誤り。

実装修正は最後の分類に対して行います。未知同位体、選択した中心だけのlabel、
atom-order permutation、full/pseudo atrop、負電荷共鳴系を回帰に含めます。
不確実な中心はラベルを推測せず、typed refusalまたはunresolvedにします。

### 3. #635 SMARTS差分を意味単位で解く

状態: source候補`11a4ea27`で14,306→200セル。Python/WASMのSMARTS APIは、
入力の索引を保つperceived aromatic viewに対して照合します（Rust coreは不変）。
残る200セルは、すべて`[Rn]`/`[kn]`のring数の意味（SSSRか、対称化ring集合か）の差です。

差分を原子primitive、結合、芳香族性、再帰SMARTS、ring、stereo、
logical operatorへ分割します。各修正は小さなtruth tableと、RDKit版・設定を固定した
比較を必要とします。parse成功とmatch互換は別の契約です。

### 4. A6 MMFF94の正しさを閉じる

状態: 同一座標でのenergyと各termのgateは、source候補`13d70a2e`で満たしました
（`benchmarks/2026-09-25-mmff94-per-term-energy.md`）。それ以外のgateは未完了です。

速度だけで完了にしません。次を別gateとして扱います。

- 同一座標でのenergyと各term;
- convergence、iteration、timeout、cancel accounting;
- stereo、bond sanity、gross clash;
- fixed cohortと独立holdoutでのconformer quality;
- source候補と公開packageの速度。

新しいembedding optionや3D機能の追加は、この基礎gateより後です。

### 5. 次のRDKit rebaselineを準備する

RDKit 2026.03.6の記録はhistoricalとして固定します。新しい公式stable artifactが
公開された後にだけ、Python wheelとofficial npm packageをhash付きで固定し、
SMILES、CIP、SMARTS、Morgan、binding overhead、browser laneを再実行します。
upstream PRや未公開版から利用可能性を推測しません。

## 常設ゲート

### Compatibility Contract

各比較profileは次を機械可読に保持します。

- chematicと比較対象のversion/hash;
- corpus identityと重複・露出状態;
- operation、options、support domain;
- success、failed、refused、skippedの件数;
- exact、numeric tolerance、semanticのどの比較か;
- 再現commandと生成artifact。

### Binding Contract

共有操作はRustを参照実装とし、Python/Node/WASMで次を守ります。

`input_count = success + failed + refused + skipped`

各結果はoriginal index、stage、typed reasonを保持し、cancel時は未処理範囲を
明示します。`failed == 0`だけを全件成功の意味にしません。

### Parser Security

固定malformed corpusをprocess isolation下で実行し、panic/crash 0、時間・memory
上限、全入力のterminal accounting、typed refusalを確認します。Rustであること
自体を安全性の証拠にはしません。

### Stereo Torture Suite

atom-order permutation、SMILES spelling、file round-trip、selected-center label、
macrocycle/atrop、競合で確認された再現可能な不具合を収録します。競合の不具合を
取り込む目的は優越宣言ではなく、同種の退行を防ぐことです。

## 共通リリース候補の合格条件

- roadmapで今回対象にしたissueが、回帰・証拠・境界説明付きでclose可能である。
- workspace test、clippy、binding contract、documentation/evidence consistency、
  parser-security、dependency/license gateが宣言した環境で通る。
- README 3言語、CHANGELOG、validation、benchmark index、version metadataが一致する。
- source-onlyの数値を公開packageの数値として書いていない。
- tag、push、registry publish、GitHub release、Pages更新は個別に確認する。

## 中断と再開

長時間gateが1時間を超える場合は、安全な区切りで中断し、次を記録します。

- 実行commandとsource SHA;
- 完了した範囲と未完了範囲;
- 中断が正しさ判定へ与える影響;
- 再開command。

中断したgateを合格として扱いません。ローカル候補が通っても、公開packageや
公開channelの確認を省略しません。

## 記録先

- 現在の優先順位: [`ROADMAP.md`](https://github.com/kent-tokyo/chematic/blob/main/ROADMAP.md)
- 未完了と依存関係: [`roadmap-open-work.md`](roadmap-open-work.md)
- 互換性境界: [`compatibility-scope.md`](compatibility-scope.md)
- 検証概要: [`validation.md`](validation.md)
- benchmark方法: [`benchmark.md`](benchmark.md)
- 日付付き実測: [benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)
- 完了した利用者向け変更: [`CHANGELOG.md`](https://github.com/kent-tokyo/chematic/blob/main/CHANGELOG.md)

過去の長い計画本文はGit historyに残します。完了した実装経緯をこの文書へ追記し
続けず、必要な場合だけ日付付きevidenceへリンクします。
