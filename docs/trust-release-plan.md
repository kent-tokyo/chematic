# chematic 1.x Trust Release — 実行計画

更新日: 2026-09-21。期間: 2026-09-14〜2026-12-06（12週間の作業配分案）。
対象: v1.0.17公開後のTrust Release候補。日付は納期の保証ではなくレビュー時点。
状態: **T0の測定packetとT1.1–T1.4、現行版T2.6、公開artifact T3.4の一部を実装済み。Trust RCの全出口は未達**。
v1.0.17は公開済み。各版のrelease境界はCHANGELOGを参照する。
公開完了と、この計画の候補条件充足は別に扱う。

2026-09-13の実装進捗: T0のperception予算伝播・実測候補数・fail-closed回帰、
SMARTS dumpの完了フッター/全行会計、独立baseline/candidate packet、API profile
dashboard、release metadata v2、公開channel実測recordを実装した。candidateは
残差12→21の退行で採用停止。P0–P2全体の完了やsealed評価の成立を意味しない。

[ROADMAP](../ROADMAP.md) の短期実行順を定める。
[A0–A6精度計画](rdkit-accuracy-plan.md) の完了条件は維持する。
優先度P0/P1/P2は緊急度、既存Phase P0–P6は製品領域、A0–A6は精度目標、
T0–T6は今回の作業ID。番号は相互に置き換えない。

9月21日の優先順位更新: **A0の三つの不採用packetを保全し、非sealed分類と次freezeの
前提を閉じる**。T1.6凍結候補は一度ずつ評価済みで、いずれも採用しない。
T3.6の高速化はPR #555 (`7d98dcd3`)で統合され、対応範囲のbit一致と複数browser/hostの
速度優位を記録済み。今後は公開package再測定、当初の統計/資源条件との差分確認、回帰維持へ移す。
[第7節](#7-2026-09-21-trust完了に向けた実行順)が新しい実行packet、ROADMAPが優先順を持つ。
新しい競合回帰・BatchResult変更は9月16日の凍結candidateへ混ぜず、別commitで進める。

2026-09-16に公開npm artifact と official RDKit.js の固定10k browser
scorecardをPR #541で公開した。小さいWASM、local no-store ready、parse/writeの
結果は比較可能だが、parse-inclusive fingerprintではRDKit.jsが速い。local readyは
CDN/network測定ではない。3回のprocess-tree peak RSS診断はあるが、共有pageを
重複計上し得るためunique memory・メモリ予算合格・cross-hostの証拠にはしない。
同日にPR #542で`rustls`を0.23.45へ更新し、hosted Security AuditのCargo Auditと
parser-security corpusは成功した。これらはT1.6のunused-data attestationやTrust RC
全体の代替条件ではない。

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
- [公式トップ](https://chematic.io/)、[PyPI](https://pypi.org/project/chematic/)、
  [validationページ](https://chematic.io/validation/)の版整合性はT2で直接監査する。
  過去の取得値や検索キャッシュを現在の配信状態として扱わず、公開版・測定版・
  各言語ページ・配信assetを分けて確認する。v1.0.14の公開workflow成功だけでは
  全ページの同期完了を証明しない。
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

依存関係の基本形は **T0 → T1 → T2 → T3/T4/T5 → T6**。
完了済み工程を再開する意味ではなく、現在の着手順はROADMAPと第6節を正本にする。
公開表示の棚卸し、未使用データ取得、gold候補収集は初日から並行する。
新しいsilent corruption・panic・資源制御欠陥は領域にかかわらず先行する。
外部レビュー待ちでも、ローカルの契約・回帰・文書作業は進める。

## 4. 実装単位と受入条件

### T0 — 比較の信頼性を回復する（目安2–4実働日）

- [x] T0.1 実行開始時にsource commit/diff、binary hash、比較器、入力/query hash、
  期待ID・行数・全セル数を固定。入力parse失敗も1行として保持し、
  存在しないcorpusへのhand-only fallbackは正式ゲートでは失敗にする。
- [x] T0.2 子processの終了codeを確認してから集計し、完了footerと全件会計が一致した
  場合だけatomicにfinal artifactへ昇格。途中読込、重複/欠落ID、budget exhaustion、
  match上限到達、stale binary、oracle parse失敗を区別する。
  RDKit側のdefault match上限も固定・検知する。
- [x] T0.3 共有探索まで同じ予算を伝播するか、意味が異なる上限を別の型とstatusにする。
  訪問数・候補数を実測し、0/1/小予算/通常予算・中断・巨大環の負例を追加。
  成功結果は予算を変えても同じ、打切り結果は成功扱いしない。
- [x] T0.4 独立buildしたbaseline/candidateで全5,021×31を再測定。
  5残差分子、bicyclo/adamantane、周辺構造、原子/結合順の並べ替えを含める。
  退行なし・対象クラス改善・資源制約合格が揃わなければhybridを昇格しない。
  改善がない複雑化は除去候補とする。現行の完了footer付きlaneは21残差・
  RDKit query parse失敗10,042セルを基準に分類する。旧12不一致・18拒否は
  別laneの履歴であり、この分母へ混ぜない。

2026-09-13の独立測定では、`f2302b66` baselineは残差12、`a259507d`
candidateは残差21（いずれも155,651セル、RDKit parse error10,042、整列失敗0）となった。
比較packetは有効だが、candidateは9残差退行のため**採用不可**。出自とハッシュは
`validation/results/rdkit-smarts-baseline-candidate-v1.0.14.json`に記録する。

既存作業先: `crates/chematic-smarts/examples/rdkit_parity_dump.rs`、
`rdkit_ring_model.rs`、`scripts/rdkit_ring_parity_diagnosis.py`。
出口: 正しい完了packetだけ通る負例検証と、同条件before/after。
「未完了」「失敗」「未測定」を明示できることも必須。

### T1 — API単位のCompatibility Contract（目安4–7実働日）

- [x] T1.1 `validation/manifests/rdkit_accuracy_v2.json`と
  `validation/cross_binding_contract.json`を参照する操作profile索引を作る。
  SMILES parse/write、canonical string/semantic identity、芳香族性、CIP、
  SMARTS/substructure、ECFP/Morgan、MOL/SDF V2000/V3000を必須行にする。
- [x] T1.2 各行へAPI名/各binding、Stable/Experimental/Unsupported、
  native/compat profile、oracle版/設定、exact/numeric/semantic、
  全入力/成功/拒否/失敗、許容差、source/binary/corpus hash、
  measured_at、再現コマンド、根拠artifactを格納する。
  要件充足率と化学的正解率を別列にする。
- [x] T1.3 `scripts/generate_compatibility_dashboard.py`を拡張し、
  JSONとMarkdownを同じ入力から決定的に生成。
  opaque属性保持とtyped意味理解、match集合と原子写像、graphと文字列の
  一致を別行にし、欠測値を0や100%へ変換しない。
- [x] T1.4 native ECFP、RDKit-compatible Morgan、legacy `rdkit_compat`ラッパーを
  実際の呼出先で分類。radius、nbits、chirality、aromaticity、count/sparse、
  bitInfoごとに保証を書く。nativeの過去recallを互換APIの現状に流用しない。

T1.1–T1.4は`validation/compatibility_profiles.json`と
`scripts/check_compatibility_profiles.py`で実装し、生成dashboardへ反映済み。
未測定は`not_measured`、既存結果が限定的な領域は`partial`として残す。
- [ ] T1.5 2025.09.3の回帰laneを保存し、2026.03.6を別laneで固定する。
  差分分類と全対象操作の再測定後にだけ主要oracle版を変更する。
  RDKit.jsはnpm版と同梱RDKit版を別々に記録する。

SMARTSについては同じv1.0.14 build・5,021分子・31 queryで2026.03.6 laneを
2026-09-13に追加した。2025.09.3 laneの10,042 oracle parse-error cellsは2026.03.6で
発生せず、残差は21から22になった。結果は
`validation/results/rdkit-smarts-oracle-lanes-v1.0.14.json`へ分離して保存し、
2025.09.3を既定regression oracleのまま維持する。全対象操作の再測定・RDKit.jsの
同梱版記録は未完了なので、T1.5全体は未完了。
- [x] T1.6 新規10,000件の取得版・ライセンス・出自を固定し、
  開発2,000/封印評価8,000へ分割。既存の開発・回帰全群との
  parent/scaffold/同一構造重複を監査し、重複は開発側へ移す。
  `scripts/generate_rdkit_identity_audit.py` と
  `scripts/prepare_sealed_accuracy_cohort.py` は、RDKit固定の
  canonical/parent/Murcko scaffold キー、source hash、candidate freeze、
  unused-data attestation を必須にする準備器として実装済み。candidate tagは
  実Git commitと同じobjectへ解決する**annotated tag**であること、attestationはattestor/
  timestamp/statement/source hashを持ち対象source hashと一致すること、さらにattestation
  timestampがcandidate tagのtagger timestampより前でないことを検証する。旧preflight
  の既存corpusはsealedと表示してはならない。
  `scripts/fetch_chembl_sealed_candidate.py` は既存の小規模ChEMBL入力と
  重ならないoffsetから、response hash・取得日時・CC BY-SA 3.0出典を保持する
  ローカル候補を取得する。取得だけでは未使用性もsealed statusも主張しない。
  公開データでも「本開発で未使用」は成立しうるが、出自と露出履歴が必要。
  実行前にprotocol・seed・候補buildを凍結し、結果閲覧後の再調整は
  新しい候補試験として記録する。既存7,737件は封印群に再分類しない。
  正式封印時は`validation/templates/unused-data-attestation.template.json`を複製して
  maintainerが候補のannotated tagを先に作成してから記入し、
  `--attest-unused --attestation-file`とcandidate commit/tagを同時に渡す。
  テンプレート自体はattestationではない。
  2026-09-13のlocal preflightでは、ChEMBL APIの14,543 single-fragment候補
  （source SHA-256 `f50156bf…d63e61`）からcanonical 1件・scaffold 3,183件を
  既存scopeとの重複として除外し、11,359件のeligible poolを得た。そこから
  2,000 development / 8,000 holdout候補を決定的に抽出し、監査ハッシュと件数だけを
  `validation/results/sealed-cohort-preflight-v1.0.14.json`へ保存した。candidate
  commit/tagとunused-data attestationをまだ得ていないためstatusは
  `prepared_not_sealed`であり、評価結果はまだ算出していない。
  この旧preflightは履歴として保持する。以下の9月16日の別sourceによる封印記録が
  現行であり、旧sourceを後から未使用データへ読み替えない。
- [x] T1.6を2026-09-16に新しい候補後sourceで再実行した。annotated tag
  `trust-eval-candidate-20260916` は`c2682e3aa75c21566c86ce1ade9cbd052838c694`を
  凍結し、その後に取得した11,689-row ChEMBL source
  (`867c6e3f…394954af`)を対象にした。RDKit 2025.09.3 canonical/parent/scaffold
  auditはdescriptor census、ChEMBL accuracy、browser comparison 10kを参照し、
  10,239 eligible rowsから2,000 development / 8,000 sealed holdoutを決定した。
  maintainer unused-data attestation、source/identity/reference/split hashes、response
  hashesは`validation/results/sealed-cohort-preflight-trust-eval-candidate-20260916.json`
  と`validation/attestations/unused-data-chembl-20260916.json`に保存する。raw inputは
  local-onlyで、scoreはまだ計算していない。
- [x] T1.6 evaluation: raw/split hashes、露出履歴、candidate tag、source-built wheelと
  RDKit 2025.09.3を確認してから、凍結候補の8k評価を一度実行した。8,000/8,000件は
  parseできたが、分子量で6件のunsupported、TPSAで46件のstrict mismatchがあり、
  A0候補は不採用となった。生SMILESを含むrawはlocal-onlyで、commit/tag/wheel/corpus/raw
  hashesと全件集計は`validation/results/sealed-descriptor-evaluation-trust-eval-candidate-20260916-20260920.json`
  に記録する。選別的な再試行はしていない。
- [x] このholdoutは露出済みであり、後続候補のtuningや未使用評価に再利用しない。
  新候補のsealed判定は露出履歴と再凍結手順で決める。A5の第三者gold/reviewは別ゲート。
- [x] T1.7 ordinary-V3000 interchange baselineを固定する。PR #544で
  RDKit `2025.09.3` と Indigo `1.46.0` をversion-pinned readerとして、通常V3000の
  semantic round trip、SGROUPの作成・編集後の外部reader受理、relative stereoと
  COLLECTION順序を検証した。結果は
  `validation/results/v3000-{rdkit,indigo,sgroup-external-reader,rdkit-relative-stereo}-v1.0.15.json`
  に保存する。WASMはSGROUP構文をbounded typed inspectionとして公開する。
  coordination chemistry、haptic bond、polymer expansion、ENDPTS/ATTACHの意味解釈、
  typed SGROUPの化学的編集は明示的に対象外であり、opaque retentionを編集可能性や
  相互運用完全性として扱わない。今後の拡張はこの境界ごとに外部reader gateを追加する。

既存A0 validator・raw-accounting・cross-binding runnersを拡張して使う。
出口: 全操作の状態を生成でき、互換と宣言した範囲は全件coverageと既定許容差を満たす。
拒否によって安全性は維持できても、宣言済み対応範囲の互換達成には数えない。

### T2 — release・文書・移行経路（目安3–5実働日）

- [x] T2.1 `release-metadata/`、`generate_release_metadata.py`、
  `check_release_metadata.py`、`check_release_docs_consistency.py`を拡張する。
  stable release、candidate、各oracle、測定版を別フィールドで管理。
  manifestの版が公開タグと同じでも、その後の実装を当該タグの公開成果として表示しない。

v2 metadataはstable release、candidate、oracle lanes、measurement-version policyを
分離する。v1.0.14ではRDKit 2025.09.3をactive regression lane、2026.03.6を
planned separate laneとして記録し、測定済みとは扱わない。
- [x] T2.2 GitHub tag/release、PyPI wheel/sdist、crates.io、docs.rs build、
  npm tarball/dist-tag、公式サイトの直接HTTPを照合する。
  offlineはmanifest検証、onlineは公開物実測。未到達を成功にしない。

2026-09-13のchannel recordではGitHub Release/tag、npm、PyPI、crates.io、docs.rs、
公式サイトを確認した。PyPIはJSON APIと`pip index`、crates.ioはREST APIの403を
回避してworkspace外から`cargo info --registry crates-io chematic@1.0.14`を実行し、
公開crateの取得・解決まで確認した。
`validation/results/release-channel-verification-v1.0.14.json`は
`release_ready=true`を記録する。最初のPyPI 404観測は伝播途中の一時値であり、
現行の公開状態として残さない。
- [x] T2.3 サイト別repoはversion付きmetadataとdigestを検証して取り込む。
  日本語/英語/中国語、CDN import、canonical URL、cache、構造化データも確認。
  全registryが公開完了するまでcurrentを進めず、失敗はpartial publicationと表示。
  immutable artifactは上書きせず、未完了stepだけを再試行できるようにする。

公式サイトrepoの中央`site-data.json`をv1.0.14、npm tarballのWASM raw 4,012,715 bytes /
gzip 1,461,941 bytes、公開済み`chematic-mcp` 1.0.14へ同期した。site dependencyと
lockfileも`@kent-tokyo/chematic@1.0.14` tarball/integrityへ更新し、site CIと
Cloudflare Pages deployの成功、英語ホームとvalidation URLのv1.0.14表示を直接確認した。
- [x] T2.4 短いREADME各言語から生成dashboardへ誘導し、
  `docs/rdkit-migration.md`とPyPI原稿をAPI profileに整合。
  fingerprint再作成、保存indexの互換性、strict/relaxed、stereo未対応、
  molecule解放・Worker利用・エラー処理を実行可能な移行例にする。
  既存PyPI配布物のREADMEは編集できないため修正は次リリースへ含める。

英日中READMEはdashboardへリンク済みで、RDKit migration guideにはnative/RDKit
fingerprint profile、保存index再構築、不完全入力のfail-closed境界、stereo/canonical
identity、WASM handle解放とbrowser/Worker error例への導線を追加した。PyPI source
READMEも同じdashboard/migration guideへ誘導する。既公開v1.0.14のimmutable PyPI
READMEは変更せず、次の配布物へ反映する。
- [x] T2.5 release候補のtarball/wheelをclean環境へインストールし、
  `tsc --noEmit`、mypy/pyright、Rust公開API、MCP schema smokeを通す。
  clone済みworkspaceでの成功だけではpackage利用品質を合格にしない。

公開npm tarballは空の一時directoryへinstallし、Node 24.5.0でWASM byteを明示
初期化してbenzene (`C6H6`, 6 atoms) を処理できた。再現用の
`scripts/smoke_npm_package.mjs`をweb-target CIへ追加し、実測recordは
`validation/results/npm-clean-install-v1.0.14.json`に保存した。これはnpmの
一経路だけであり、PyPI wheel/sdist、Rust crate、TypeScript typecheck、MCPは未測定のため
T2.5全体は未完了のまま。

PyPI v1.0.14 wheelも空venvでimport、benzene処理、mypy、pyrightを通過し、
公開`chematic-mcp` crateは隔離`cargo install --locked`後にlegacy initializeと
20-tool `tools/list`を返した。npmのTypeScript consumer fixtureも`tsc --noEmit`で
通過した。全結果は`validation/results/package-clean-install-v1.0.14.json`に固定する。
隔離downstream Rust consumerも`chematic = '=1.0.14'`をcrates.ioから解決し、
`chematic::smiles::parse`でbenzeneをcompile/runできた。これでT2.5の指定する
Python、npm TypeScript、Rust公開API、MCP schema smokeを記録環境で満たした。
cross-platform coverageは別のT3/T6出口であり、このclean-install smokeの合格範囲へ
含めない。

- [x] T2.6 v1.0.15の公開channel・サイト・clean installを再確認する。
  `validation/results/release-channel-verification-v1.0.15.json`はGitHub Release、npm、
  PyPI、crates.io、docs.rs、Pagesをverifiedとして`release_ready=true`にした。加えて
  run `35061960770`はbinary-only PyPI wheelでLinux CPython 3.9、macOS CPython 3.13、
  Windows CPython 3.13のimport/version/benzene formulaをartifact化した。Linux 3.9は
  v1.0.15の実配布wheel境界であり、newer Linux Pythonへ一般化しない。次候補では同じ
  channel recordとruntime smokeを新しいartifactとして再実行する。

出口: 公開channelごとのversion/digest/時刻と実行可能な例。
historical benchmarkの版は維持し、「全ての数字を最新版にする」同期はしない。

### T3 — WASM・MCPを利用可能な製品にする（目安5–8実働日）

- [ ] T3.1 ESM/npm、CDN、Web Workerの最小例を固定し、
  asset配置・CSP・初期化失敗・型エラー・cleanupを検査する。
  Chromium/Firefox/WebKit、Nodeを対象とし欠測engineを明示する。
- [ ] T3.2 Explorerで1万件をrelease必須、10万件を次段階の容量ゲートにする。
  SDF/CSV/SMIのstreaming、固定batch、backpressure、進捗、cancel、
  部分失敗位置、再試行、exportの順序を検査。
  Explorerは現在10,000件を保持し、250行だけをvirtual表示するWorker経路へ
  移行済み。`scripts/explorer_worker_10k_smoke.mjs` がChromiumで全10,000件
  の解析・状態・DOM上限を確認する。SDF/CSV streaming、partial failure、retry、
  export順序は未完了。
- [ ] T3.3 同一input/seed/profileでscalar/batch/Worker/native結果を比較。
  初期化後にネットワークを遮断してローカル機能が動くことを検証する。
  外部取得toolは明示的なonline機能として別扱い。
  `scripts/explorer_native_worker_parity.mjs` はethanol、ethylamine、不正ring
  closureをnative scalar CLI・native batch CLI・実module Workerで照合する。
  canonical SMILES、主要記述子（1e-9以内）、拒否結果を比較し、Chromium local
  runは通過した。CIにも独立jobを追加した。全browser・network遮断後の動作は未完了。
- [ ] T3.4 RDKit.jsと同じoperation/inputでraw/gzip、cold init、
  parse/write、ECFP、検索、peak memory、UI応答を測る。
  同一host/browser、cold 20回・warm 5回以上、p50/p95、
  failure/coverage、測定APIとメモリ定義を記録する。
  100万分子換算は実測と別列で、線形外挿の仮定を明記する。
  公開1.0.15の10k・3ブラウザ・20反復と3回のRSS診断は実施済み。
  再作成ではなく、公開1.0.17/次候補の別lane、検索、operation意味論、
  resource会計を補う。npm `2026.3.6` / runtime `2026.03.6`を別fieldにする。
- [ ] T3.5 MCP全toolのruntime schemaからtool数と入出力例を生成。
  型付き化学エラー、oversized input、中断、structured outputを確認する。
  `scripts/check_mcp_runtime_inventory.py` は実際の stdio binary から
  20 tool・input/output schema・代表 structured output を確認する。入力
  上限と中断のwire-level回帰は追加で必要。
- [~] T3.6 [PF0–PF5](parse-morgan-performance-plan.md)のsource最適化と複数browser/hostの
  速度gateは統合済み。公開package再測定、当初20対×3セッション/層別bootstrapと
  実施済み10–20対/t区間の差分、p95/resource条件は別途照合して残す。
  native ECFPと互換Morganは別集計、正しさ/coverageの退行は禁止。
- [~] T3.7の既知長batch契約はPR #559で全bindingに固定済み。第7節の残作業は
  中断・stream・Worker/MCP adapter・resource測定であり、基礎会計を未実装として
  扱わない。

1万件の暫定予算: timeoutを含む全件会計、cancel応答p95 ≤250ms、
UI heartbeat gap p95 ≤100ms、batch working set ≤256MiB
（索引/入力保持分は別計測）。対象hostで凍結する前の設計値であり実績ではない。
不合格ならbatch/Worker設計を修正する。10万件の上限はW3で実測から事前固定。
出口: 公開packageで3つの導入例、1万件処理、browser結果一致が再現できる。

### T4 — Parser Security Benchmark（目安4–7実働日、以後継続）

既存の800 malformed入力・10 oversized・20 gzipの検証は維持し、
既存security workflowとcorpusを再利用する。以下100件は出典付きの
競合回帰laneとして追加するもので、既存ゲートを縮小するものではない。

- [~] T4.1 SMILES/SMARTS/MOL V2000/V3000/SDFの5形式を対象に
  少なくとも各20件（合計100件）の出典付き回帰を用意する。
  公開issue/OSS-Fuzz/minimized reproのライセンス、対象版、原因、
  byte hashを記録。CVE番号だけを根拠に形式の異なる入力を流用しない。
- `validation/parser_security_corpus_v1.json` は各形式1 valid controlと19件の
  CheMatic作成 boundary mutation（合計100件）を固定し、Open Babel、RDKit、
  Indigo の公開reproducer URL、原因、対象範囲、各payload hashを記録する。
  upstream入力のバイト列を再配布したものではないため、`origin` と fixture
  licenseを明示する。upstreamの修正版・影響版を固定する比較laneは未完了。
- [ ] T4.2 各engineをnetworkなしの子process/コンテナへ隔離し、
  wall time・peak RSS・exit/signal・panic・結果statusを記録。
  まず修正版の固定releaseを比較し、脆弱版実行を標準CIにしない。
  `scripts/run_isolated_parser_security.py` と
  `crates/chematic-cli/examples/parser_security_case.rs` はこの実行境界を
  実装済み（wall-time 2秒、Linuxのaddress-space 256 MiB、exit/signal/status、
  Linux runnerの`VmHWM` peak RSS記録）。macOSのdiagnostic runはこの環境で
  address-spaceを強制できないため`not_measured`とし、networkもproxy削除を
  network隔離と呼ばない。Linux PR jobは固定runnerを測定前にbuildし、fresh
  network namespace（GitHub hosted runnerの`sudo unshare --net`）で実行して
  raw JSONをartifactへ残す。namespace境界を作れないrunnerではfail-closedとする。
- [~] T4.3 入力≤1MiB、1ケース2秒/256MiB、5形式固定corpusを暫定ゲートとする。
  oversizedは事前検出で構造化エラー。重いstress群は別の上限と分母。
  process killは資源制御の証拠であり、parserの正常完了には数えない。
  隔離runnerはこの上限を既定値にし、各形式にvalid controlを必須にする。
  `parser_security_case` は≤1MiBを構造化拒否し、runnerは各形式のvalid
  controlを要求する。2026-09-13のmacOS functional diagnosticは100/100
  expectation一致・最遅384.83msだったが、memory/networkを強制していないため
  release evidenceではない。
- [~] T4.4 PRは固定corpus、nightlyは各target 15分のfuzz、
  releaseは各target 1時間と固定regressionを実行する配分を準備。
  panic/crash/limit violation=0を要求し、timeoutも失敗として残す。
  sanitizer/Miri・FFI/依存層の検証を既存security workflowへ接続する。
  `security.yml`のLinux PR jobは固定corpus、256MiB address-space、wall-time、
  per-case RSSを実行してartifactを保存する。nightly/release fuzzの時間割は未完了。
- [x] T4.5 valid inputの負例対照も含め、全入力を拒否して安全率を上げない。
  外部ライブラリで再現しないケースもunsupported/not_applicableとして記録。

出口: 再配布可能なcorpus、固定版runner、全件会計、CI raw artifact。
「Rustだからmemory corruptionが存在しない」は採用しない。
unsafe/FFI/依存層とpanic/DoSを含む検証境界を明示する。

### T5 — Stereo Torture Suite（目安8–15実働日、研究残差は別）

- [~] T5.1 最低300構造をtetrahedral、E/Zと共有carrier、ring/cage、
  負電荷共鳴、P/S、同位体、enhanced stereo、未対応立体へ事前配分。
  curated regressionと未使用challenge群を分離し、同じscaffoldで独立性を水増ししない。
  `stereo_torture_suite_development.jsonl` は既存回帰と固定sampleの
  300 unique structuresを固定し、カテゴリ内訳は生成manifestに記録する。
  source追加時に再生成し、古い内訳を現行値として転載しない。
  `stereo_torture_suite_gate.py`は全300件について
  RDKit semantic identityとchematic canonical再parse安定性を検査し、CIで再実行する。
  これは開発回帰だけであり、未使用challenge、事前配分の全カテゴリ、CIP labelの
  絶対正解は未達のまま残す。
- [~] T5.2 atom/bond orderとSMILES spellingを固定seedで各32変換。
  #149/#503の既存K=1,024診断は継続し、32変換へ縮小しない。
  SMILES→MOL/SDF V2000/V3000→再parse、canonical再適用、
  stereoisomer keyの衝突を検査する。
  固定RDKit 2025.09.3 seedでCIP corpus 155構造のSMILES spellingを各32変換し、
  元表記を含む5,115入力のsemantic identity、canonical spelling、再parse idempotencyを
  `stereo_spelling_invariance_gate.py`でCI検査する。atom/bond order 32変換、
  MOL/SDFの全round-trip、stereoisomer key collisionは未達として残す。
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
| W1: 9/14–9/20 実績 | T4.6 core＋binding回帰、ordinary V3000、開発stereo suite、T3.6 source速度gate、T3.7既知長batch契約 | 実装済み境界と残る受入条件。8k scoreは一度実行され、候補は不採用 |
| W2: 9/21–9/27 | A0不採用記録の保全、T5.7/T1.9残回帰、T1.5実測lane準備 | 新freezeを混ぜない全行会計、原子対応/真理値表、artifactごとのprovenance。外部artifact不足はunavailableとして記録 |
| W3: 9/28–10/4 | A2残差、T3.7のstream/cancel/adapter、T1.5新artifactがあれば比較 | canonical fail-closed境界、binding横断の欠落/重複0、旧/新版の差分。未公開artifactは待機 |
| W4: 10/5–10/11 | T3 runtime、T2同期、T3.6残条件とRC監査 | 10k cancel/offline/資源、配布候補hash、Trust RC合否とblocker一覧 |
| W5–W8: 10/12–11/8 | A1–A4残差、10万件容量、複数browser/OS、A5依頼packet | 操作別compat達成、coverage・速度・memoryの測定 |
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

以下は9月19日から継続する技術packet。9月20日以降の優先順と追加の受入条件は
第7節を参照する。T3.6のsource速度gateは統合済みで、再実装を計画しない。

9月19日の競合レビューを既存T/A/Phaseへ統合する。新Phaseや別の並行ロードマップは
作らない。T0の独立packetは完成し、12→21残差の候補は不採用。残差分類と新候補は
未完了だが、測定済みpacketそのものを「これから実装」として繰り返さない。

### 前提確認と競争方針

- [公式npm metadata](https://registry.npmjs.org/@rdkit%2frdkit/2026.3.6)を確認。
  npm版は`2026.3.6`であり、既存測定のruntime版`2026.03.6`と分ける。
  型定義と依存なしの公開packageは利用可能。既存10k recordのSRIはregistryと一致。
  tarball全体・WASM raw・gzipを混同せず、過去rawの版fieldは書き換えない。
- [RDKit #9572](https://github.com/rdkit/rdkit/pull/9572)は9月10日にmerge済み。
  [Pyodide recipe #538](https://github.com/pyodide/pyodide-recipes/pull/538)は確認時点でopen。
  ビルド成立と配布済み製品を区別し、利用可能になった時点でruntime込みの総転送量・
  startup・利用可能APIを別laneで測る。今はPyodideの機能数を追う開発をしない。
- [COSMolKitのCIP統合](https://github.com/cosmol-studio/COSMolKit/commit/1235890a04286c531b9b4a9b52419a6efcdb0253)
  は作者が全workspace合格ではないと明記。直接競合として扱うが、未検証の性能・
  安定性ランキングは作らない。T6.1に安定版/RC、公開WASM hash/サイズ、再現可能な
  benchmark、利用事例を記録し、downloads/starsを実利用者数に読み替えない。

勝ち筋は**限定した化学操作を、導入しやすい型付きAPIで、ローカルかつ応答性を
保って正しく実行できること**。優位性は操作別の測定で示す。全面RDKit互換、
全面3D、protein、polymer編集の拡大は今回の必須にしない。

### 1. T4.6 — 環番号の自社riskを閉じた（P0、2026-09-19）

`crates/chematic-smiles/src/writer.rs`は100以上を`%100`のように出力する一方、
`parser.rs::parse_ring_num`は`%`の後を2桁だけ読み、open-ring tableも100要素。
これはソース上の往復不整合の懸念だった。[Indigo #3867](https://github.com/epam/Indigo/pull/3867)
を回帰の入口にし、次のbounded implementationで解消した。

- [x] parser/writer双方で`%(n)`を扱う。`%nn`（00–99）は既存どおり、拡張形式は
  100以上のみ。writerが以前出した曖昧な`%100`は出さない。
- [x] `%()`、100未満の`%(n)`、閉じ括弧欠落、非数字、`u32` overflow、旧`%100`を
  strictに拒否。open labelは値で配列確保せず、最大100,000個の現在openなlabelへ制限。
- [x] 99/100/101境界のparser fixtureと、121 closureを持つ12×12 graphの
  write→parse count-preservationを追加。`cargo test -p chematic-smiles --lib`は
  226 passed。小番号と既存stereo回帰も同じsuiteで通過した。
- [x] PR #562で同じvalid/invalid fixtureをRust/Python/Node/WASMの共有契約へ追加し、
  `C%(100)CC%(100)`の受理と、短い括弧形式・閉じ括弧欠落・曖昧な旧`%100`の拒否を
  binding横断で固定した。
- [ ] canonical writerのatom順入替とtimeout/メモリ制限を通す。これは残る
  canonical/resource gateであり、上のparser binding完了とは別。

### 2. T1.6 / T1.5 — 精度とoracleの版を分離する（P0、目安2–4実働日）

T1.6の封印前提は完了済み。次は上記のhash/build/protocol監査と未実行の評価である。
2025.09.3 regression、公開2026.03.6 runtime、上流修正commitの診断を別laneにし、
masterのCIP修正が公開npmにも含まれると推測しない。masterを測る場合はSHAで固定。
T4/T5修正は開発fixtureで行い、旧候補のsealed結果を新候補の合格証明に転用しない。
8kだけでA5の±0.1 percentage point同等性を証明できるとも仮定しない。

### 3. T3.4 — 比較の空欄を埋める（P0、目安3–5実働日＋測定待ち）

- [ ] 公開1.0.17を再buildせず測り、次候補buildは別armにする。既存runnerを拡張し、
  package/runtime/oracle版、SRI/SHA、host/browser、全入力/拒否を保存する。
- [ ] parseのsanitize/perception範囲、writeのsemantic preservation、Morganのradius/
  bit数/chirality/profileを先に照合。prepared-objectとparse-inclusiveを別表示する。
  native ECFPとRDKit-compatible Morganの速度/精度を混ぜない。
- [ ] 部分構造検索は同じquery集合・chirality・uniquify・match上限で、hit判定と
  mapped embeddingを別lane化。索引build時間/容量と検索時間を分離し、全候補/timeout/
  unsupportedを保持。MinimalLibにない操作はN/AとしPython版で穴埋めしない。
- [ ] fresh processでcold 20回、warm 5回以上、交互arm、p50/p95と反復ごとの値を保存。
  二つ目のhostは別集計。3ブラウザの結果を混ぜて優位性を主張しない。
- [ ] JS heap snapshot、WASM linear memory、process-tree peak RSS、unique memoryを
  別metricにする。RSSはidle差分・sampling間隔・共有page重複を記録し、測れない指標は
  unavailable。1万/10万の実測と100万への外挿を別表にする。
- [ ] CDN/remote cold startはnetwork条件/cache/transfer encodingを固定した追加lane。
  未測定でもlocal結果は公開可能。速度負けを含めてdashboard/再現コマンドを公開し、
  正しさ・coverageが退行したarmを「高速化成功」にしない。旧追加1.10x SMILES目標は
  復活させない。新しいParse＋Morgan目標はT3.6の別契約で評価する。

### 4. T5.6 — CIP/atropの上流変更を境界テストにする（P1、目安3–5実働日）

[RDKit #9577](https://github.com/rdkit/rdkit/pull/9577)と
[#9190](https://github.com/rdkit/rdkit/pull/9190)はmerge済み。次のケースを既存T5 suiteへ
追加し、宣言した対応域のwrong confident label・情報損失・false mergeを0にする。

- [ ] 選択中心が未選択中心に依存する分子で、全分子assignmentからの射影と比較。
  chematicに同等の選択APIがない場合は能力差として記録し、API追加を合格条件にしない。
- [ ] 未知同位体の番号fallback、full/pseudo atrop、cleanIt/replaceExistingTags相当の
  preserve/replace、古い2D/3D evidence、8/9/10員環境界と負例を固定する。
- [ ] atom/bond mapを保持した32順列とSMILES spelling、V2000/V3000往復で照合。
  未対応atropはラベルなしで成功させず、型付き拒否または情報を保持した明示診断にする。
  RDKit間で判定が変わるケースはA5へ回し、P系OracleUnstableを安易に解除しない。

### 5. T1.8 — attachment labelとcollapseを区別する（P1、目安2–4実働日）

[RDKit #9454](https://github.com/rdkit/rdkit/pull/9454)を参照する。
現行CXSMILESには一般atom-labelのread/writeがあるが、それだけではattachment identityや
collapseの保証にならない。V3000のENDPTS/ATTACH対応とも別契約である。

- [x] `attachment_point_label_number()`は完全な正の`_AP<n>`だけを`u32`として認識し、
  `_AP0`、符号、空白、文字混在、overflowを拒否する。`CxSmiles::marked_attachment_point()`は
  degree-one wildcardかつ有効labelだけを識別する。通常atom、多重degree wildcard、
  label文字列単体をattachment pointとして扱わない。`chematic-smiles` test suiteで検証済み。
- [x] version-pinned RDKit 2025.09.3 gateで、RDKitが書いた`_AP1`、最大`u32`、`_AP0`、
  負数、文字混在、overflow、non-dummy、2つのdegree-one wildcardを
  RDKit→CheMatic CLI→RDKitで検査する。atom-map keyed label、通常topology、
  attachment identityの保持/拒否を8/8で要求する。これはMDL collapse、V3000
  ENDPTS/ATTACH、またはgeneric CX queryの意味論を主張しない。
- [x] canonical SMILESでもmapped wildcardを`[*:n]`として再出力し、mapを持つ
  wildcardのcanonicalize→reparse→canonicalizeでidentityを固定する。これは
  attachment collapseの可否判定を追加するものではない。
- [ ] attachment identityとcollapse可能性を別判定にする。
  label番号はMDL ATTCHPT位置ではない。label-only collapseの位置1規則、direction/query/
  bond制限を仕様化し、情報を保持できないcollapseは明示的に拒否する。
- [ ] 次RCは保持/拒否契約を必須とし、意味論未確定のcollapse API追加は後続minorへ。
  既存ordinary-V3000 gateは維持し、opaque retentionを意味的編集として宣伝しない。

### 6. T3 / T2 — 機能表より導入から完走までを示す（P1、目安3–5実働日）

T3.1–T3.3の未達部分を埋める。公開packageのESM/TS/Worker導入、初期化後offline、
1万件streaming・partial failure・cancel・再試行・export順序を同じ小さな利用例で示す。
暫定UI/メモリ予算はT3記載値を測定前にhostごとに凍結する。3ブラウザ・Nodeの欠測を
表示し、API capabilityと測定版が一致した結果だけをサイト/MCP説明へ反映する。

上記日数は作業配分の目安。上流ビルド、測定host、第三者reviewは別の待ち時間。
計画更新は実装・新ベンチマーク・公開の証拠ではなく、本更新でsealed inputは開封しない。

## 7. 2026-09-21 Trust完了に向けた実行順

この節は次の1〜3か月の追加・再配置を既存T0–T6/A0–A6に統合する。
新しい製品Phaseは増やさない。各項目は計画であり、競合issueの存在だけでは
CheMaticの不具合や優位性が確定したことにはならない。

### この期間の意思決定ルール

1. **Trust Releaseを先に閉じる。** 新しい競合機能、描画機能、3D embedding
   optionは、A0/A2/T3.7の受入を遅らせない。競合由来の情報は、まず
   version-pinned な開発回帰か typed unsupported 境界へ落とす。
2. **旧sealed入力は最適化データにしない。** 既に不採用となった三候補のraw、
   個別差分、集計外の派生情報を修正の標的にも次候補の検証にも使わない。
   A0の分類は公開済みまたは新規の非sealedデータ、仕様、独立生成fixtureだけで行う。
3. **Browserの価値は制御可能性で測る。** サイズや単発速度は operation/version/
   host を固定した補助指標である。主張の中心は、local-only execution、typed errors、
   cancel、resource limits、row accounting、deterministic output の同時成立とする。
4. **比較は同じ操作だけを順位付けする。** API意味論、拒否条件、prepared state、
   input/output、測定hostのいずれかが異なるlaneは、比較記録には残すが速度・互換性の
   勝敗表には入れない。

### 実装順を固定する追加ルール

この週は新機能の数ではなく、証拠と契約を閉じる。各作業の開始条件と出力を次の
ように固定する。

| 順位 | 実行単位 | 次の成果物 | 開始・停止条件 |
|---:|---|---|---|
| P0-1 | A0 evidence recovery | CI検証済みの開発packetを、影響した非sealed分類とbinding会計に追随させる | PR #606で三つの安全な不採用summary、52-row分類、7,737-row binding影響を固定済み。新source、overlap audit、attestation、annotated tagなしに新holdoutを作らない。 |
| P0-2 | T1.5 rebaseline preparation | Python/native/npmのartifact・backend・型/例外・操作設定を持つ実行可能lane | 次RDKit版は公式release artifactを一次確認してから測る。予定版、issue、PR、stubの量だけではlaneを開始しない。 |
| P1-1 | T5.7/T1.9 upstream regressions | identity-renumber stereo、V3000 E/Z query、必要なら選択atom/bond CIPのsource-pinned fixture | まず報告版と現行pinでreproducerを確認する。CheMaticが同じAPIを持たない場合はlossless round-tripかtyped refusalを検査し、互換達成とは呼ばない。 |
| P1-2 | T3.7 controlled batch | unknown-stream schema、Worker/MCP adapter fixture、cancel/export migration example | 既知長と未知長を一つの成功数へ潰さない。全bindingで入力index、stage、終端理由、処理済みprefix、未読状態を表すまでruntime claimを広げない。 |
| P1-3 | Browser proof maintenance | 同一操作契約のWorker/stream/cancel/resource scorecard | parse/write、compatible Morgan、native ECFP、startup、memoryを別laneとし、package/version/hostを固定する。 |
| P2 | A6 fail-closed 3D | 既存3D gapのtyped failure/quality evidence | 上流のMMFF/embedding機能はwatch対象。基礎のtyping、charge、gradient、convergence、timeout、stereoが未完のまま広い3D APIを追加しない。 |

競合ウォッチは週次で一次情報・reproducer・release artifact・影響APIを台帳化する。
それ以外の競合主張は、計画の根拠ではなく未検証のwatch noteに留める。COSMolKitを
含むPure-Rust競合に対しては機能数を追わず、全入力の終端会計、cross-binding同一性、
local-only制御、再現可能な比較packetを差別化の受入条件にする。

### 一次情報を確認して修正した前提

| 出典（9月20日確認） | 確認できた範囲 | 計画への反映 |
|---|---|---|
| [RDKit #9601](https://github.com/rdkit/rdkit/issues/9601) | nanobindの型・例外・引数・lazy inputを整理するopen tracking issue。[stub PR #9613](https://github.com/rdkit/rdkit/pull/9613)はclosed/unmerged | T1.5は配布artifactの実際のwrapper backend/型/例外を検査。nanobind移行完了や次版収録を前提にしない |
| [RDKit latest release](https://github.com/rdkit/rdkit/releases/latest) | 確認時点はRelease_2026_03_6 | 2026.09.1は予定する再比較先。公開日・Python/native/npm同時提供は仮定しない |
| [RDKit #9629](https://github.com/rdkit/rdkit/issues/9629) | 2025.09.4のbicyclic amineでidentity renumber後にstereocenter countが変わるというopen報告 | T5.7でring/cache/state依存を検査。2026.03.6で再現済みとは記載しない |
| [Indigo #3914](https://github.com/epam/Indigo/issues/3914) | 1.48.0rc1のV3000 queryが反対E/Zにもmatchするというopen報告 | T1.9でquery truth tableを検査。既存1.46.0 ordinary-MOL gateと区別 |
| [COSMolKit #1](https://github.com/cosmol-studio/COSMolKit/issues/1) | 0.3.0の利用者報告と現sourceのskipped accessor不足。issue自身も公開packageでの独立再現は未実施と明記 | T3.7の自己完結した件数/失敗契約に反映。競合packageの検証済み欠陥とは宣伝しない |
| [RDKit #9625](https://github.com/rdkit/rdkit/pull/9625) / [#9626](https://github.com/rdkit/rdkit/pull/9626) | protected-atom tautomer / MMFF initial embeddingはいずれもopen PR | T6で監視。今回のTrust RCへ新機能を追加する根拠にはしない |

### 1. T1.6 / A0 — 不採用評価の保全と次freezeの準備（P0、1–3実働日＋評価時間）

- [x] `trust-eval-candidate-20260916`、`trust-eval-candidate-20260920`、
  `trust-eval-candidate-20260920c` を、
  tag、oracle、source/split/raw hash、build hash、attestation、全行会計とともに
  一度だけ評価した。前者は分子量/TPSA、後者は分子量 7,998/8,000 strict（同位体表の
  残差2件）で不採用。第三候補は7項目が8,000/8,000 strictだったが、TPSAが
  7,977/8,000 strict（23件不一致）で不採用。安全な公開summaryは各候補の
  `validation/results/sealed-candidate-trust-eval-*-summary.json` に固定した。
- [x] 全raw/splitはlocal-onlyかつ露出済みとして明記した。再実行・tuning・次candidateの
  未使用評価へ流用しない。raw分子一覧を調査ログやrepositoryへ出さない。
- [x] **非sealedの開発データだけで** A0の宣言済み範囲を分類する。封印結果から特定入力や
  閾値へ合わせ込まず、独立の既知/生成テストと仕様根拠を使う。修正候補ごとに、影響binding
  回帰、変更影響表、公開主張への影響を記録する。descriptorは官能基・電荷・同位体・
  tautomer/芳香族性の非sealedな層別を先に固定し、全体一致率だけを次freezeの根拠にしない。
  公開52官能基/電荷/同位体/tautomer/芳香族性のgreen probeでは、HBA/HBD/TPSA/LogP/MR/
  Fsp3/芳香族環数、`exact_mass`、`heavy_atoms`、`rotatable_bonds`、明示的な
  `rdkit_molecular_weight`は52/52 strictだった。Kekulé 2-pyridoneは、単一carbonylに
  隣接して芳香族化される環内NHをRDKitのaromatic `[nH]` descriptor typeとして扱う公開
  regressionであり、二つのcarbonylに隣接するphthalimide型imideを同型へ昇格させない。
  native
  `molecular_weight`のSeと同位体3件はnative IUPAC mass tableとRDKit mass profileの
  宣言済み差であり、CLIのnative massをRDKit互換値と偽装しない。結果は
  `validation/results/tpsa-functional-group-probe-current-2026-09-20.json`と
  `validation/results/descriptor-functional-group-classification-current-2026-09-20.json`に
  固定する。これは公開development分類であり、unused採用評価ではない。三つの却下summary、
  52-row分類、7,737-row Rust/Python/Node-WASM binding影響を
  `validation/a0-development-packet.json` と
  `scripts/check_a0_development_packet.py` がCIでまとめて検証する。packetはraw rowを
  埋め込まず、候補の採用や次freezeの評価を主張しない。
- [ ] 次の採用判断が必要になった時点で、別source・別抽出・重複除外・source hash・
  attestation・annotated tagを新規に固定する。新しいcohortを作るだけでは合格とせず、
  candidate buildとfixed oracleで一度だけ実行する。
- [ ] A0 packetの成立とA1–A4各指標の合否を別判定にする。失敗時は結果を固定し、
  A5独立gold/reviewを代替条件にしない。

出口: 両不採用結果の追跡可能な保全、次freeze前の非sealed開発証拠、将来の候補ごとに
独立した未使用data手順。次のRDKit公開をこの準備の待ち条件にしない。

### 2. T1.5 / T3.4 — 次期RDKit rebaseline（P0準備、2–4実働日＋公開待ち）

- [x] 既存oracle 2025.09.3/2026.03.6と次のstableを別laneに固定する
  `validation/rdkit_rebaseline_manifest.json` と検証scriptを追加した。
- [ ] source/package/runtime/wrapper backend、設定、corpus hash、toolchain/OSを
  実測laneごとにmanifestへ記録する。
  Python/native/npmは個別に入手確認。欠測laneはunavailableとして待機し、先に出た
  Python packageの版から公式RDKit.js公開を推測しない。
- [ ] 同一の露出済みcorpusでSMILES/CIP/SMARTS/Morganとbindingを比較する。
  Pythonではscalar/batch、入力型・例外・keywordとcall overheadを同じAPI境界で測る。
  Boost/nanobindは配布物の実体を特定し、別backendの結果を統合しない。
- [ ] browserは同じpacked output契約で3 enginesを比較し、parse単独、prepared FP、
  parse+FP、search、startup、memoryを個別に記録する。API意味論が揃わないlaneは
  速度順位の対象から外す。sealed精度8kはこの回帰・性能集合へ流用しない。
- [ ] RDKitの型stub/nanobind移行は、公開artifactで確認できた場合だけ Python laneへ
  追加する。型注釈の量ではなく、`str | bytes`、path-like、iterable、例外、keyword、
  scalar/batchの実行時契約をCheMaticのtyped binding contractと並べて記録する。
- [ ] 差分を自社退行/oracle変更/契約差/未解決へ分類し、goldが必要ならA5へ送る。
  新oracleへ自動追従して既存出力を変更しない。公開dashboardは測定済みの新旧版を併記し、
  移行後も旧raw/再現コマンドをhistoricalとして保存する。

出口: 利用可能なlaneの全行差分、失敗/拒否、速度区間、API変更の移行表。
再比較先の公開遅延は、現行pinでのTrust RC候補監査を停止させない。

### 3. T5.7 / T1.9 — 競合報告由来の回帰（P1、2–4実働日）

- [x] **T5.7 の狭い再現:** PR #558はatom mapを保ったidentity/reordered atom回帰で
  potential stereocenter集合が不変であることを固定した。これはRDKit #9629の報告クラスを
  開発gateへ輸入したもので、RDKitの同一不具合を再現したという主張ではない。
- [x] **T1.9 の保存/拒否境界:** PR #558はV3000→V3000でopaque query属性を保持し、
  unmodelled V3000 query constraintを通常分子形式へ平坦化する変換を拒否する。
- [x] **T5.7 の順列/状態ゲート:** `C1CCN2CCCC2C1`と通常の第三級アミン負例を、32固定
  seed順列、clone/reparse、ring情報の初期化順、descriptor/CIP呼出し順で比較する。atom mapで
  potential-center集合・CIPを、map endpointでbond骨格を照合し、countやcanonical spellingだけで
  不変としない。入力はstereo未指定のためCIP/E/Z assignmentが空であることも固定する。
- [~] **T1.9 query semantics の最初の観測:** RDKit 2025.09.3が生成した対称2-buteneと
  置換基が異なるchloro/bromo alkeneのE/Z V3000 query/targetを4×4で比較し、CheMatic往復前後でRDKitのisomeric identityに加え、
  `HasSubstructMatch(..., useChirality=True)`の真理値表も保存されることを確認した。
  Indigo 1.46.0の観測表も変わらない。一方Indigoはsource自体で両targetにmatchするため、
  Indigoの4行はquery意味保存の合格ではなく外部semantic-lossの再現記録である。
  `validation/results/v3000-indigo-ez-query-truth-table-current-2026-09-21.json`。
- [ ] **残るT1.9 query semantics:** SMILES入力、V3000保存/再読込、cross-engine往復を
  報告版・現行pin・他engineへ広げる。stereo match設定を固定し、unspecified stereoと
  stereo-insensitive設定は別の対照群にする。原子/結合対応、query predicate、
  SGROUP/COLLECTIONを記録し、reader受理や文字列保存をquery意味保存の代替にしない。
- [ ] RDKit/Indigoの報告版と現行pinを分離して再現を記録する。CheMaticに同じquery APIが
  なければ能力差を明記し、情報を失う経路をtyped unsupportedとして拒否する。
  unsupportedでの安全性合格とquery互換性達成は分ける。

既存`stereo_torture_suite_gate.py`と`v3000_*_semantic_gate.py`を拡張する。
新規ケースは開発回帰として出典・入力hash・期待値の根拠を保存し、sealed群には追加しない。
自社でsilent corruptionを再現した場合はP0へ繰り上げる。

### 4. T3.7 — BatchResultの基礎契約を全bindingで固定した（P1、2026-09-20）

現行`SmilesBatchCanonicalizer`は`input_index`とaccepted/rejectedを保持し、WASMの
`canonicalize_smiles_batch_json`はschema v1の`record_count`/`records`を返す。
rejectionは文字列、envelopeの`complete`は処理完了を表す。これを土台にし、公開v1の
意味を破壊せずversioned adapter/opt-in APIとして拡張する。既存binding fixtureを共用する。

- [x] PR #559で排他的な終端区分`success / failed / refused / skipped`、既知長batchの
  `input_count = success + failed + refused + skipped`、0-based original row index、
  parse error stage、`all_succeeded()`をRust/Python/Node/WASMの共通fixtureで固定した。
  `complete`は処理完了、`all_succeeded()`は完了かつ全入力成功を表す。
- [ ] cancel、未知長stream、retry、chunk境界、Worker/MCP adapter、export stageと
  resource/time limitを同じ会計契約へ接続する。cancel後の未読範囲を架空の`skipped`に
  せず、`complete=false`と未読範囲を明示する。
- [x] PR #580でLocal Compound Explorerの既知長importにこの境界を適用した。record capと
  cancelでは`unprocessed`とterminal reasonを表示し、古いWorker応答が後続importの状態を
  上書きしないgeneration guardを追加した。Chromium/Firefox/WebKitのsmokeで確認済み。
  これはExplorerのUI経路だけであり、下記の全binding/10k Worker出口の代替にはしない。
- [x] PR #584（`5678360b`）はExplorerのauto-detected CSV `File.stream()`へunknown-length
  adapterを統合した。EOF前の件数を推測せず、cancel時は完了済み、観測済み
  未終端、`unread input=unknown`、terminal reasonを別々に示す。10k Chromium smokeは
  normal CSV import/exportのinput order/statusに加え、10k CSV中断で`complete=false`と
  未読範囲を確認する。これはブラウザUI adapterの検証であり、Rust/Python/Node/WASM public
  batch schemaやMCP streamの契約完了を意味しない。
- [~] 本stream contractはRust/Python/WASM/Nodeの共通fixtureを使い、public
  opt-in stream APIでunknown-length schema v1を固定する。`finish`は観測済みのpending rowを
  すべて処理する唯一のEOF経路、`stop`は`cancelled`/`time_limit`/`resource_limit`/
  `producer_error`/`consumer_closed`の列挙理由と、処理済みprefix・観測済み未処理・
  `unread_input=unknown`を返す。未知の理由はbinding errorとし、`failed()==0`を全入力成功の
  根拠にしない。Worker/MCP、retry/export、10kのresource測定は未完了。
- [ ] clean-installで型定義とruntimeを照合し、10k Workerで順序/元index/理由/会計の
  一致、消失0、二重計上0、測定済みの資源上限を確認する。

出口: 公開契約、versioned schema、shared fixture、全bindingの実測record、移行例。
追加APIの配布区分はT6.2の互換性方針で決め、先にpatch番号を予約しない。

### 5. A2 → T3/T2 → A6 — 残りの作業配分

- **A2（P1）:** #503の4 componentを原因・負例・K=1,024 gateごとに分割し、P系CIPは
  A5判定前のabstentionを維持。結合E/Z componentは完全なcarrier割当て以外を選ばず、
  証明できない場合はstable keyをfail-closedのまま残す。T5.6、T1.8、T4.6の残binding/
  表現境界も同じ受入表へ統合する。identity-renumber、SMILES spelling、clone/reparse、
  V3000 round-tripの各変換で atom/bond correspondence と confident/abstained outcome を
  比較し、canonical stringやcenter countだけを合格条件にしない。
- **T3/T2（P1）:** 完了したT3.7契約を使う10k Workerでcancel・backpressure・memory/time limit・
  offline・partial exportを測定。既存のcancel p95≤250ms、heartbeat p95≤100ms、
  working set≤256MiBは測定前にhost条件と固定し、未達を記録する。npm clean-install、
  Python型/runtime、3browser+Nodeの欠測を表示し、dashboardはartifactから生成する。
- **T3.6（維持）:** 既存速度の根拠を保ち、20対×3セッション/層別bootstrap、p95、資源、
  公開tarball再測定の未確認条件を埋める。確認済みbrowser結果をnative/Python全般へ外挿しない。
- **A6/T6（P2、誤計算はP0）:** benzene等の基本系を含むtyping/charge/term/gradientを
  same-coordinateで切り分け、収束/timeout/立体保持/配座品質を別集計。新embedding option
  より既存gapを優先する。MMFF由来の初期化選択やprotected-atom tautomerのような上流新機能は
  監視対象に留め、既存の3D pipelineが失敗時に誤った成功を返さないことを先に証明する。
  A5独立reviewは依頼packetまでローカルで進め、第三者判定待ちを明示。

次Trust RCの必須条件に、T3.7のWorker/stream/cancelを含む宣言batch経路と、
T5.7/T1.9の保存・拒否回帰を含める。
T1.6旧候補の合格だけで新RCを承認せず、候補差分・影響gate・公開主張の対応表を確認する。
この計画更新では化学コード、封印入力、ベンチマーク結果、リリース状態は変更しない。
