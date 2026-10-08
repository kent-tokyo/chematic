# chematic 1.x Trust Release 実行計画

更新日: 2026-10-07。**v1.0.37リリース系列**です。v1.0.36の6配布経路と
日付付き公開成果物測定は記録済みで、v1.0.37公開物の確認は別に行います。

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
- 反応83件はPR #755のchecked sourceで80件が生成物グラフ・原子由来・テンプレートマップの
  全軸で一致。残り3件は双方無効です。
  Linux/macOSの公開用設定で作ったsource wheelとWASM NodeテストはCI通過。
  公開レジストリのartifactは未測定です。
- SMARTSのopt-in source-wheelは309,982/310,000一致し、残り18件を型付きで
  非対応とします。公開v1.0.30の200件の差分を置き換える測定ではありません。
- 3Dは公開v1.0.31 macOSで265/265の構造・立体・clash判定を通過しましたが、
  収束は100/265です。Linuxでの2件の立体失敗と独立conformer品質が残ります。

版・コーパス・残差は[検証報告](validation.md)と
[benchmark index](https://github.com/kent-tokyo/chematic/tree/main/benchmarks)に固定します。

## 実行順

優先順位と各数値の詳細は[ROADMAP](https://github.com/kent-tokyo/chematic/blob/main/ROADMAP.md)と
[検証報告](validation.md)に置き、この文書では判定条件だけを保持します。

| 順 | ゲート | 合格条件 |
|---:|---|---|
| 1 | 反応・SMARTS | 83反応のgraph/origin/mapを全行分類し、元の57件を維持する。SMARTSはmatch-set、Boolean、typed refusalを別計上し、公開artifactと独立コーパスで再測定する。 |
| 2 | A6 MMFF94 | [#739](https://github.com/kent-tokyo/chematic/issues/739)のLinux立体失敗を解き、atom type、同一座標の各energy term、収束、stereo/clash、独立conformer品質を別gateにする。品質前に速度や新しいembedding optionを成果扱いしない。 |
| 3 | CIP・E/Z | 既知の5件は理由付き棄権を維持し、原子順・SMILES表記・file往復で誤った確信ラベルを出さない。E/Zの既存source回帰を保つ。 |
| 4 | 新RDKit版 | 公式stable公開後に版・hashを固定して再比較する。RDKit 2026.03.6の証拠はhistoricalとして残す。 |

棄権や無効入力を一致に数えません。対応外と実装誤りを分け、修正時は小さな
回帰テストと版固定の全量比較を要求します。

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
