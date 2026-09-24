# chematic 1.x Trust Release 実行計画

更新日: 2026-09-24。公開版は **v1.0.24**、次の開発候補は
**v1.0.25** です。

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
  18/10,000から0/10,000へ減らしました。長時間permutation gateは中断済みで、
  その完了や公開package再測定の成果としては扱いません。
- v1.0.23はRDKit-defined atom-pair/torsion、MACCS、QED、TPSA、Murckoの一致を
  改善しました。3本の5,000行source laneと43/44操作のmatrixは共有2 vCPU実行であり、
  公開packageまたはWASMの性能・完全互換性の証拠ではありません。
- v1.0.24はその出力を保持したままperceptionとSMARTS存在判定のhot pathを短縮しました。
  v1.0.23との差分は7コーパス・1,402,080出力行で0件です。計時は共有2 vCPUの
  source実行だけであり、公開package、WASM、全環境の性能主張ではありません。
- 残る主要差分はCIP #634、SMARTS #635、A6のenergy/term・timeout・
  conformer qualityです。

## 実行順

### 1. #632 SMILES E/Zを閉じる

合格条件:

- focused parser/writer/canonical regressionが通る。
- RDKit 2026.03.6の固定10,000件でgraph差0、semantic差0を維持する。
- 28 component x 1,024 relabelの長時間gateを完走する。完走しない場合は、
  未実行であることをrelease判断に明記する。
- 測定外のcoupled shapeではstable-keyがfail-closedのままである。

### 2. #634 CIP差分を分類して解く

差分を次の4種類に分けます。

- atom/bond correspondenceまたは比較器の問題;
- RDKit版・表現依存などoracleの問題;
- chematicが明示的に非対応とする領域;
- chematic実装の誤り。

実装修正は最後の分類に対して行います。未知同位体、選択した中心だけのlabel、
atom-order permutation、full/pseudo atrop、負電荷共鳴系を回帰に含めます。
不確実な中心はラベルを推測せず、typed refusalまたはunresolvedにします。

### 3. #635 SMARTS差分を意味単位で解く

差分を原子primitive、結合、芳香族性、再帰SMARTS、ring、stereo、
logical operatorへ分割します。各修正は小さなtruth tableと、RDKit版・設定を固定した
比較を必要とします。parse成功とmatch互換は別の契約です。

### 4. A6 MMFF94の正しさを閉じる

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
