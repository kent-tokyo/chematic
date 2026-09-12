# chematic 1.x Trust Release — 実行計画

更新日: 2026-09-13。期間: 2026-09-14〜2026-12-06（12週間の作業配分案）。
対象: v1.0.14公開系列と、その前段であるv1.0.13候補。日付は納期の保証ではなくレビュー時点。
状態: **計画策定。実装完了・リリース承認ではない**。

2026-09-13の実装進捗: T0のperception予算伝播・実測候補数・fail-closed回帰、
SMARTS dumpの完了フッター/全行会計、RDKit accuracy profileのdashboard表示を実装した。
P0–P2全体の完了やsealed評価の成立を意味しない。

[ROADMAP](../ROADMAP.md) の短期実行順を定める。
[A0–A6精度計画](rdkit-accuracy-plan.md) の完了条件は維持する。
優先度P0/P1/P2は緊急度、既存Phase P0–P6は製品領域、A0–A6は精度目標、
T0–T6は今回の作業ID。番号は相互に置き換えない。

## 1. 戦略と確認した前提

訴求の中心を「ブラウザ・Python・Rust・AI agentsから使える、
対応範囲と失敗を明示した化学計算カーネル」にする。
差別化はインストールから最初の正しい結果までの再現性、
データのローカル処理、型付きAPI、追跡可能な検証証跡で示す。

- RDKit 2026.03.6にはsynthon-space shape searchとAI発見問題の修正がある。
  [公式リリース](https://github.com/rdkit/rdkit/releases/tag/Release_2026_03_6)。
  比較基準版の候補にするが、既存2025.09.3測定の版ラベルを書き換えない。
- COSMolKitはRust・Python・ブラウザと、固定RDKitに対する対応範囲別の
  検証を掲げる直接競合。ChEMBL 37の約290万record・34 phaseという規模は
  **競合自身の報告であり、こちらの独立再現結果ではない**。
  [README](https://github.com/cosmol-studio/COSMolKit)、
  [検証文書](https://github.com/cosmol-studio/COSMolKit/blob/main/VALIDATION.md)。
  規模に加え、入力の適格性、比較対象の内部状態、除外、再現手順をレビューする。
- Open Babel 3.2系は機能・OSS-Fuzz由来の修正を継続している。
  [公式リリース](https://github.com/openbabel/openbabel/releases)。
  CDKのRGroup/CXSMILES対応も確認できる。
  [CDKリリース](https://github.com/cdk/cdk/releases)。
- IndigoとOpenChemLibは表現・描画・相互運用の監視対象。
  [Indigo](https://github.com/epam/Indigo/releases)、
  [OpenChemLib](https://github.com/Actelion/openchemlib/releases)。
  提供された個別不具合・CVE・PR日付は、それぞれの一次資料を確認してから
  回帰ケースへ採用する。全件を今週の確定事実として転載しない。
- 今回の取得では[公式トップ](https://chematic.io/)はv1.0.14、
  Published 2026-09-11を表示し、[PyPI](https://pypi.org/project/chematic/)とも一致。
  0.20.1が現在もトップにあるという前提は採用しない。
  [validationページ](https://chematic.io/validation/)の取得本文だけでは
  pinned版の整合性を確定できないため、直接HTTP・各言語・配信assetを監査する。
  検索キャッシュと公開ページの現状は区別する。
- PyPIのMorgan互換表にはbit-identicalではない旨が残り、
  ローカルにはnative ECFPとRDKit互換APIが別に存在する。
  API名・設定・リリース版単位の対応表が必要。
- MCPは現在のソースで20 tools。全toolがofflineではなく
  `pubchem_lookup`はネットワークを使う。ブラウザのローカル演算と、
  optional native-InChI/FFI・外部取得toolの境界を明示する。

## 2. 最初に訂正すべき検証上の問題（T0の入力）

以前の「共有環モデルでSMARTS不一致12→3件」は、実行中の部分dumpを
全件と誤認した集計であり撤回する。2026-09-13に保存済みdumpを再読した結果:

| 比較lane | 分子 | query | 一致 | 不一致 | chematic拒否 | RDKit非対応 |
|---|---:|---:|---:|---:|---:|---:|
| 旧一時dump（履歴のみ） | 5,021 | 31 | 145,579 | 12 | 18 | 10,042 |
| 旧共有symmetrized SSSR単独（履歴のみ） | 5,021 | 31 | 145,576 | 15 | 18 | 10,042 |
| 旧hybrid（履歴のみ） | 5,021 | 31 | 145,579 | 12 | 18 | 10,042 |

上の3行は履歴比較であり、現行の採否判断には使わない。現行の正本は
RDKit 2025.09.3、`uniquify=True`、各matchの原子index集合を正規化した
永続コーパスlaneである。全155,651セル、parity一致145,588、残差21、
RDKit SMARTS parse error 10,042、atom alignment failure 0。
atom-map対応や全embedding一致を証明する比較ではなく、ビルド出自を固定した
候補採用packetでもない。

再集計した一時入力のSHA-256:
- 共有単独 `e56d5d72db3af0faa97d23ec73410f43e9df8d0d6a6758bc5e2a3d1a7ca3bc05`
- hybrid `673eb82a4dd8317bbca6400382f525ed2bdbbe4fa4c66e0d9e9c0e20b5e58d81`

footer検証済みの永続コーパスlaneの正本artifactは
`validation/results/rdkit-smarts-direct-chembl-5000-footer-verified-v1.0.13.json`。
入力SHA-256は`1c47371dcbe37f4e0a141bf545b72bf238de2761fa3894fa251a552d84728d3e`、
dump SHA-256は`233731d5c2c1d6af07008523628003e5faa53d4a0679a63387dcfe6c9ca734bc`。
RDKit 2025.09.3（installed Python package）、155,651セル、整列失敗0、
parity一致145,588、残差21、RDKit parse error10,042。
artifactのSHA-256は`55816ac9e71df833416476c9b763938f0a098f5b1d72b0f2ca547596ace19563`。
これは現行laneの完了証拠であり、独立baseline/candidate比較やsealed評価ではない。

`build_shared_symmetrized_ring_model`は利用者の正の`max_candidates`を共有探索へ
渡し、cap超過時の`candidates_examined`を実測値として返す経路へ更新済み。
bounded SSSRとSMARTS parityの予算回帰も追加済みである。残るT0は、独立した
baseline/candidateのビルド出自を固定した採用packet、残差分類、hybridの採否である。
環数の多少だけでモデルを選ぶhybridは正しさの証明にならない。
原因となる候補生成・基底・対称化規則を検証し、原子順依存を測る。

## 3. 優先順と担当領域

| 作業ID | 優先度 | 既存Phase / 精度目標 | 成果物・責任範囲 | 主な依存 |
|---|---|---|---|---|
| T0 証拠・予算契約の修復 | P0 | Phase P0/P2、A0/A4 | SMARTS runner、全件会計、hybrid採否 | なし |
| T1 Compatibility Contract | P0 | Phase P0/P2/P3、A0–A4 | API profile manifest、生成dashboard、候補測定packet | T0 |
| T2 公開情報同期 | P0 | Phase P0/P6、A0 | release metadata、registry実測、サイト取込、移行文書 | T1のschema。棚卸しは即時 |
| T3 WASM・MCP利用品質 | P1 | Phase P3、A3/A4 | install smoke、Worker、Explorer、型・エラー契約 | T1/T2 |
| T4 悪意ある入力の回帰 | P1 | Phase P1/P3、A0/A4 | 出典付きcorpus、隔離runner、資源上限CI | T0。収集は即時 |
| T5 Stereo Torture Suite | P1 | Phase P2/P4、A2/A5 | 不変性・round-trip・誤統合・CIP判定 | T1/T4。既存反例調査は即時 |
| T6 維持運用・独立評価・限定3D | P2 | Phase P5/P6、A5/A6 | release運用、第三者packet、3D既存gap解消 | T1–T5 |

主経路は **T0 → T1 → T2 → T3/T4/T5 → T6**。
公開表示の棚卸し、未使用データ取得、gold候補収集は初日から並行する。
新しいsilent corruption・panic・資源制御欠陥は領域にかかわらず先行する。
外部レビュー待ちでも、ローカルの契約・回帰・文書作業は進める。

## 4. 実装単位と受入条件

### T0 — 比較の信頼性を回復する（目安2–4実働日）

- [ ] T0.1 実行開始時にsource commit/diff、binary hash、比較器、入力/query hash、
  期待ID・行数・全セル数を固定。入力parse失敗も1行として保持し、
  存在しないcorpusへのhand-only fallbackは正式ゲートでは失敗にする。
- [ ] T0.2 子processの終了codeを確認してから集計し、完了footerと全件会計が一致した
  場合だけatomicにfinal artifactへ昇格。途中読込、重複/欠落ID、budget exhaustion、
  match上限到達、stale binary、oracle parse失敗を区別する。
  RDKit側のdefault match上限も固定・検知する。
- [ ] T0.3 共有探索まで同じ予算を伝播するか、意味が異なる上限を別の型とstatusにする。
  訪問数・候補数を実測し、0/1/小予算/通常予算・中断・巨大環の負例を追加。
  成功結果は予算を変えても同じ、打切り結果は成功扱いしない。
- [ ] T0.4 独立buildしたbaseline/candidateで全5,021×31を再測定。
  5残差分子、bicyclo/adamantane、周辺構造、原子/結合順の並べ替えを含める。
  退行なし・対象クラス改善・資源制約合格が揃わなければhybridを昇格しない。
  改善がない複雑化は除去候補とし、元の12不一致・18拒否を未達として残す。

既存作業先: `crates/chematic-smarts/examples/rdkit_parity_dump.rs`、
`rdkit_ring_model.rs`、`scripts/rdkit_ring_parity_diagnosis.py`。
出口: 正しい完了packetだけ通る負例検証と、同条件before/after。
「未完了」「失敗」「未測定」を明示できることも必須。

### T1 — API単位のCompatibility Contract（目安4–7実働日）

- [ ] T1.1 `validation/manifests/rdkit_accuracy_v2.json`と
  `validation/cross_binding_contract.json`を参照する操作profile索引を作る。
  SMILES parse/write、canonical string/semantic identity、芳香族性、CIP、
  SMARTS/substructure、ECFP/Morgan、MOL/SDF V2000/V3000を必須行にする。
- [ ] T1.2 各行へAPI名/各binding、Stable/Experimental/Unsupported、
  native/compat profile、oracle版/設定、exact/numeric/semantic、
  全入力/成功/拒否/失敗、許容差、source/binary/corpus hash、
  measured_at、再現コマンド、根拠artifactを格納する。
  要件充足率と化学的正解率を別列にする。
- [ ] T1.3 `scripts/generate_compatibility_dashboard.py`を拡張し、
  JSONとMarkdownを同じ入力から決定的に生成。
  opaque属性保持とtyped意味理解、match集合と原子写像、graphと文字列の
  一致を別行にし、欠測値を0や100%へ変換しない。
- [ ] T1.4 native ECFP、RDKit-compatible Morgan、legacy `rdkit_compat`ラッパーを
  実際の呼出先で分類。radius、nbits、chirality、aromaticity、count/sparse、
  bitInfoごとに保証を書く。nativeの過去recallを互換APIの現状に流用しない。
- [ ] T1.5 2025.09.3の回帰laneを保存し、2026.03.6を別laneで固定する。
  差分分類と全対象操作の再測定後にだけ主要oracle版を変更する。
  RDKit.jsはnpm版と同梱RDKit版を別々に記録する。
- [ ] T1.6 新規10,000件の取得版・ライセンス・出自を固定し、
  開発2,000/封印評価8,000へ分割。既存の開発・回帰全群との
  parent/scaffold/同一構造重複を監査し、重複は開発側へ移す。
  公開データでも「本開発で未使用」は成立しうるが、出自と露出履歴が必要。
  実行前にprotocol・seed・候補buildを凍結し、結果閲覧後の再調整は
  新しい候補試験として記録する。既存7,737件は封印群に再分類しない。

既存A0 validator・raw-accounting・cross-binding runnersを拡張して使う。
出口: 全操作の状態を生成でき、互換と宣言した範囲は全件coverageと既定許容差を満たす。
拒否によって安全性は維持できても、宣言済み対応範囲の互換達成には数えない。

### T2 — release・文書・移行経路（目安3–5実働日）

- [ ] T2.1 `release-metadata/`、`generate_release_metadata.py`、
  `check_release_metadata.py`、`check_release_docs_consistency.py`を拡張する。
  stable release、candidate、各oracle、測定版を別フィールドで管理。
  ソースtreeが1.0.13でも、その後のdirty実装を1.0.13公開成果として表示しない。
- [ ] T2.2 GitHub tag/release、PyPI wheel/sdist、crates.io、docs.rs build、
  npm tarball/dist-tag、公式サイトの直接HTTPを照合する。
  offlineはmanifest検証、onlineは公開物実測。未到達を成功にしない。
- [ ] T2.3 サイト別repoはversion付きmetadataとdigestを検証して取り込む。
  日本語/英語/中国語、CDN import、canonical URL、cache、構造化データも確認。
  全registryが公開完了するまでcurrentを進めず、失敗はpartial publicationと表示。
  immutable artifactは上書きせず、未完了stepだけを再試行できるようにする。
- [ ] T2.4 短いREADME各言語から生成dashboardへ誘導し、
  `docs/rdkit-migration.md`とPyPI原稿をAPI profileに整合。
  fingerprint再作成、保存indexの互換性、strict/relaxed、stereo未対応、
  molecule解放・Worker利用・エラー処理を実行可能な移行例にする。
  既存PyPI配布物のREADMEは編集できないため修正は次リリースへ含める。
- [ ] T2.5 release候補のtarball/wheelをclean環境へインストールし、
  `tsc --noEmit`、mypy/pyright、Rust公開API、MCP schema smokeを通す。
  clone済みworkspaceでの成功だけではpackage利用品質を合格にしない。

出口: 公開channelごとのversion/digest/時刻と実行可能な例。
historical benchmarkの版は維持し、「全ての数字を最新版にする」同期はしない。

### T3 — WASM・MCPを利用可能な製品にする（目安5–8実働日）

- [ ] T3.1 ESM/npm、CDN、Web Workerの最小例を固定し、
  asset配置・CSP・初期化失敗・型エラー・cleanupを検査する。
  Chromium/Firefox/WebKit、Nodeを対象とし欠測engineを明示する。
- [ ] T3.2 Explorerで1万件をrelease必須、10万件を次段階の容量ゲートにする。
  SDF/CSV/SMIのstreaming、固定batch、backpressure、進捗、cancel、
  部分失敗位置、再試行、exportの順序を検査。
- [ ] T3.3 同一input/seed/profileでscalar/batch/Worker/native結果を比較。
  初期化後にネットワークを遮断してローカル機能が動くことを検証する。
  外部取得toolは明示的なonline機能として別扱い。
- [ ] T3.4 RDKit.jsと同じoperation/inputでraw/gzip、cold init、
  parse/write、ECFP、検索、peak memory、UI応答を測る。
  同一host/browser、cold 20回・warm 5回以上、p50/p95、
  failure/coverage、測定APIとメモリ定義を記録する。
  100万分子換算は実測と別列で、線形外挿の仮定を明記する。
- [ ] T3.5 MCP全toolのruntime schemaからtool数と入出力例を生成。
  型付き化学エラー、oversized input、中断、structured outputを確認する。

1万件の暫定予算: timeoutを含む全件会計、cancel応答p95 ≤250ms、
UI heartbeat gap p95 ≤100ms、batch working set ≤256MiB
（索引/入力保持分は別計測）。対象hostで凍結する前の設計値であり実績ではない。
不合格ならbatch/Worker設計を修正する。10万件の上限はW3で実測から事前固定。
出口: 公開packageで3つの導入例、1万件処理、browser結果一致が再現できる。

### T4 — Parser Security Benchmark（目安4–7実働日、以後継続）

既存の800 malformed入力・10 oversized・20 gzipの検証は維持し、
既存security workflowとcorpusを再利用する。以下100件は出典付きの
競合回帰laneとして追加するもので、既存ゲートを縮小するものではない。

- [ ] T4.1 まずSMILES/SMARTS/MOL V2000/V3000/SDFの5形式を対象に
  少なくとも各20件（合計100件）の出典付き回帰を用意する。
  公開issue/OSS-Fuzz/minimized reproのライセンス、対象版、原因、
  byte hashを記録。CVE番号だけを根拠に形式の異なる入力を流用しない。
- [ ] T4.2 各engineをnetworkなしの子process/コンテナへ隔離し、
  wall time・peak RSS・exit/signal・panic・結果statusを記録。
  まず修正版の固定releaseを比較し、脆弱版実行を標準CIにしない。
- [ ] T4.3 入力≤1MiB、1ケース2秒/256MiB、5形式固定corpusを暫定ゲートとする。
  oversizedは事前検出で構造化エラー。重いstress群は別の上限と分母。
  process killは資源制御の証拠であり、parserの正常完了には数えない。
- [ ] T4.4 PRは固定corpus、nightlyは各target 15分のfuzz、
  releaseは各target 1時間と固定regressionを実行する配分を準備。
  panic/crash/limit violation=0を要求し、timeoutも失敗として残す。
  sanitizer/Miri・FFI/依存層の検証を既存security workflowへ接続する。
- [ ] T4.5 valid inputの負例対照も含め、全入力を拒否して安全率を上げない。
  外部ライブラリで再現しないケースもunsupported/not_applicableとして記録。

出口: 再配布可能なcorpus、固定版runner、全件会計、CI raw artifact。
「Rustだからmemory corruptionが存在しない」は採用しない。
unsafe/FFI/依存層とpanic/DoSを含む検証境界を明示する。

### T5 — Stereo Torture Suite（目安8–15実働日、研究残差は別）

- [ ] T5.1 最低300構造をtetrahedral、E/Zと共有carrier、ring/cage、
  負電荷共鳴、P/S、同位体、enhanced stereo、未対応立体へ事前配分。
  curated regressionと未使用challenge群を分離し、同じscaffoldで独立性を水増ししない。
- [ ] T5.2 atom/bond orderとSMILES spellingを固定seedで各32変換。
  #149/#503の既存K=1,024診断は継続し、32変換へ縮小しない。
  SMILES→MOL/SDF V2000/V3000→再parse、canonical再適用、
  stereoisomer keyの衝突を検査する。
- [ ] T5.3 assignment、abstention、incorrect label、情報損失、FP/FNを分離。
  atom-map/bond-mapで照合し、canonical文字列一致だけで正解にしない。
  非対応表現がV2000に落ちる場合、保存成功ではなく明示的拒否または診断を要求。
- [ ] T5.4 RDKit固定oracleと仕様/独立goldを併用。
  RDKit自身が並べ替えで不安定なケースはA5審査へ送り、
  既存のP系abstentionを「RDKitに一致したから」解除しない。
- [ ] T5.5 原因規則・近傍負例を先に固定してからcarrier/DFS/CIPを修正。
  入力ID・特定コーパス・「環が多い方」のような経験則で正解を選ばない。

出口: 宣言範囲でwrong confident label/情報損失/誤統合=0、
coverage後退なし、不変性成立。安全な拒否は残る精度課題として表示する。
A2全出口・A5独立判定の完了条件は変更しない。

### T6 — 維持運用・独立評価・3D（W5以降）

- [ ] T6.1 日曜レビューではRDKit、COSMolKit、Open Babel、Indigo、
  CDK、OpenChemLibを固定対象にする。release/commit、主張/独立測定、
  問題→自社の影響operation→回帰候補→優先度を記録する。
  COSMolKitの再現可能な小コーパスを確認し、共通operationだけで評価する。
  この計画は監視自動化の作成を意味しない。
- [ ] T6.2 patchは互換性を維持する修正、minorはAPI追加、majorは破壊的変更。
  Experimentalも失敗契約を明示する。通常releaseはゲート通過時に集約し、
  security hotfixは別経路。更新頻度自体を品質指標にしない。
- [ ] T6.3 A5の絶対gold・paired CI・非メンテナreview packetを準備。
  第三者不在は未完了のまま。候補報告にはrawの保存先/digest/再現手順を添える。
- [ ] T6.4 A6は既存macrocycle型/BCI/energy/gradient/timeout/stereo/配座品質に集中。
  新しい3Dアルゴリズムやprotein/Markush機能の拡張はこの期間の必須外。
  既知の誤計算は延期理由にせず修正し、A6全条件の完了まではExperimental。
- [ ] T6.5 security advisory→修正→backport→package公開→サイト同期を演習し、
  raw evidenceはCI/release artifactへ保存。Gitには小さな索引・要約・必要fixture、
  圧縮rawのURL/digest/保持期限を置き、再現に必須のartifact欠落は検証失敗にする。

## 5. 12週間の進め方と候補の境界

| 時期 | 主な作業 | レビューで確認するもの |
|---|---|---|
| W1: 9/14–9/20 | T0、T1 schema、T2棚卸し、未使用データ出自確認 | 誤った比較値の撤回、runner全件会計、profile対応表 |
| W2: 9/21–9/27 | T1 dashboard/封印準備、T2 package/サイト同期、T4 corpus | clean checkoutからの再生成、版ずれ検知、負例ゲート |
| W3: 9/28–10/4 | T3 install/Worker/1万件、T4 CI、T5既存残差 | 公開package経路、bounded input、同条件before/after |
| W4: 10/5–10/11 | T1 sealed実行、T5公開suite、RC監査 | Trust RC条件の充足/未達とblocker一覧 |
| W5–W8: 10/12–11/8 | A1–A4未達、10万件、複数browser/OS、A5依頼packet | 操作別compat達成、coverage・速度・memoryの測定 |
| W9–W12: 11/9–12/6 | A5独立判定、A6既存gap、維持運用演習 | 宣言範囲の同等性判断、3D別profile、次期backlog |

工数は各laneの概算で並行分を含む。単独作業で直列化する場合や化学的残差が難航する場合、
W4は「RC監査」を実施し、合格・発売を約束しない。毎週、残工数と依存を更新する。

Trust RCはT0、T1の宣言操作契約、T2の配布候補/文書同期、
T3の1万件導入ゲート、T4固定corpus、T5既存安全性回帰・公開suiteを必須とする。
**従来のA0封印評価、A1 core8項目の未使用評価、影響binding回帰という
次候補の条件も維持する。** 未達ならRC候補作成を延期し、成果物は開発snapshotとして表示する。
全2D互換はA0–A4、独立精度同等/優位はA5、3D同等はA6の全出口が必要。
公開後のchannel同期はT2のrelease完了条件であり、配布前のRC監査と区別する。

## 6. 次に着手する具体的な変更

最初の変更はT0.1–T0.3: SMARTS runnerの完了状態・全件会計、
hybridの実測予算とcap errorを修正する。T0.4で既存/共有/hybridの
同一入力比較を固定し、採否を決める。
次にT1.1–T1.4とT2.1を実装し、既存の生成物・release metadataを再利用する。
T4の出典確認とT5の既存残差台帳は並行で準備する。
本計画では製品コード変更、commit/push、registry公開、サイト配信、第三者への連絡は実行しない。
