# chematic 1.x Trust Release 実行計画

更新日: 2026-10-03。**v1.0.31は公開済み**です。6つの配布経路は確認しました。
WASMの分子式・反応JSON修正は、公開npmでの出力再検証が残っています。

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

- A0の8記述子は、凍結候補でdevelopment 2,000件と一度限りのsealed 8,000件に
  合格しました。使用済みデータは再びsealedとして扱いません。
- 公開v1.0.30の比較では、SMARTSに200/310,000セルの差、CIPに5件の理由付き
  判定保留が残ります。出力監査は完了しましたが、RDKit完全互換ではありません。
- 公開v1.0.30の速度証拠は、出力一致と区間条件を満たすPython 20操作など、
  記録された環境・操作に限ります。3D/MMFF94はExperimentalです。
- 反応83件はchecked sourceで76件が生成物グラフ・原子由来・テンプレートマップの
  全軸で一致。3件は型付き非対応、1件は理由付き拒否、3件は双方無効です。
  Pythonのローカルsource拡張でも同じ内訳ですが、CI wheelと公開artifactは未確認です。
- SMARTSのopt-in source-wheelは309,982/310,000一致し、残り18件を型付きで
  非対応とします。公開v1.0.30の200件の差分を置き換える測定ではありません。
- 3Dは公開v1.0.31 macOSで265/265の構造・立体・clash判定を通過しましたが、
  収束は100/265です。Linuxでの2件の立体失敗と独立conformer品質が残ります。

版・コーパス・残差は[検証報告](validation.md)と
[benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)に固定します。

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

状態: 固定1万行では9,995件一致、4件の`oracle_unstable`と1件の
`lone_pair_center`を理由付きで保留しています。保留を一致件数に含めません。
[行単位の記録](https://github.com/kent-tokyo/chematic/blob/main/validation/results/rdkit-rebaseline-issue634-v1.0.27-candidate-vs-rdkit-2026.03.6-2026-09-28.json)を参照してください。

差分を次の4種類に分けます。

- atom/bond correspondenceまたは比較器の問題;
- RDKit版・表現依存などoracleの問題;
- chematicが明示的に非対応とする領域;
- chematic実装の誤り。

実装修正は最後の分類に対して行います。未知同位体、選択した中心だけのlabel、
atom-order permutation、full/pseudo atrop、負電荷共鳴系を回帰に含めます。
不確実な中心はラベルを推測せず、typed refusalまたはunresolvedにします。

### 3. #635 SMARTS差分を意味単位で解く

状態: 公開v1.0.30は固定31万セル中200セルのmatch-set差分が残ります。
開発中のopt-in source-wheelは309,982件一致、18件を型付き非対応とし、
wrong-confidentな結果は0件です。これは公開packageでの合格ではありません。
ネイティブSSSRを暗黙に変更せず、互換設定を分けます。

差分を原子primitive、結合、芳香族性、再帰SMARTS、ring、stereo、
logical operatorへ分割します。各修正は小さなtruth tableと、RDKit版・設定を固定した
比較を必要とします。parse成功とmatch互換は別の契約です。

### 4. A6 MMFF94の正しさを閉じる

状態: 公開v1.0.31 macOSの265件ではgeometry/stereo/clashが全件通過し、
同一座標の比較可能な262件は総エネルギーが1 kcal/mol以内です。ただし
収束は100件のみです。Linux/Python 3.9では公開版とsource版の両方で53/246行が
立体エラーとなり、原因調査中です。同一座標の各termの旧source測定だけで
現在の公開版のMMFF94同等性は主張しません。

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
