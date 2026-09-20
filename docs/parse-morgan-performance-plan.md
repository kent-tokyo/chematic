# Parse + ECFP4/Morgan performance plan

更新: 2026-09-20。状態: **計画のみ。現行版の再測定・最適化は未実施**。
対象: v1.0.17を基準に、RDKitより速いParse＋fingerprint経路を作る。
既存T3.4の測定契約を使う **T3.6** の詳細計画であり、新しい製品Phaseではない。
[ROADMAP](../ROADMAP.md) が優先順、[A3](rdkit-accuracy-plan.md) が互換性の出口を持つ。

## 1. 目標と現状

最初の勝利条件は、**同じSMILESから同じ2048-bit Morgan fingerprintを返す
parse-inclusive操作で、固定RDKit.jsを上回ること**。一次対象は単一Worker・単一threadの
Chromium/WASM。Firefox/WebKit、Node、Python/nativeは別々に検証し、未測定の環境へ
「RDKitより高速」を一般化しない。RustとPythonをまたぐ不公平な比較は行わない。

過去の[公開npm 1.0.15比較](../benchmarks/2026-09-16-official-rdkit-js-isolated-browser-v1.0.15.md)
は、固定exposed 10k・macOS arm64・Chromium・20 fresh processesで次の結果だった。
数値は**プロセスごとの平均ms/molの中央値**であり、個々の分子のp50ではない。

| 操作 | CheMatic npm 1.0.15 | RDKit.js runtime 2026.03.6 | 解釈 |
|---|---:|---:|---|
| Parse API呼出し | 0.003645 ms/mol | 0.108980 ms/mol | sanitize/perception範囲の違いを再監査する。単独で互換parseの勝利としない |
| Parse + radius-2 FP | 0.315347 ms/mol | 0.134228 ms/mol | CheMaticは約2.35倍の所要時間。旧経路の同等到達には約57.4%の時間短縮が必要 |

9,999件の対応域と1件のFe(II) typed refusalを区別する。1.0.15のsource rebuildは
別artifactであり、公開packageの数値と混ぜない。上記は課題の大きさの参考で、
**1.0.17の速度でも、新しい同一出力形式ゲートのbaselineでもない**。

### プロファイルを混ぜない

- **互換Morgan（主目標）:** radius=2、2048 bits、chirality=false、bond types=true、
  redundant environments=false。他の既定値・水素・sanitize・芳香族性・ring情報も
  実際の固定oracleに照合し、manifestへ保存する。成功分子は全bit一致が必須。
- **native ECFP4（別目標）:** nativeの既存bit列・意味論を保存して高速化する。
  同じ入力・radius/幅でもhash/知覚が異なるため、RDKit Morganとの速度比較は
  「異なる定義のfingerprint生成」と明記する。互換Morganの勝利条件の代替にしない。
- raw identifier/count/bitInfo、chirality=true、他のradius/幅、検索は独立回帰対象。
  bit-only経路の高速化をdetail APIや全設定の高速化へ外挿しない。

## 2. 測定契約 — PF0で固定する

1. **三つのartifact:** 公開1.0.17を再buildせず測るarm、v1.0.17 tagのsource build、
   最適化candidate。同一toolchain/flagsでsource baselineとcandidateを作る。
   comparatorは既存npm `@rdkit/rdkit@2026.3.6` / runtime `2026.03.6`を固定し、
   package/runtime/Python oracleの版・SRI/SHA・commit・build flags・host/browserを分離記録。
   comparator更新は別laneとし途中で差し替えない。
2. **共通出力:** parse→互換前処理→FP→同じbit順のpacked 256 bytes→解放を主計時区間にする。
   現runnerはCheMatic byte vectorとRDKit `get_morgan_fp()`文字列を比較しているため、
   固定packageの型/APIを調べ、双方のbyte API、または必要な変換込みで揃える。
   オプション文字列の構築は双方でループ外。出力保持・消費も対称にし、全行の一致と
   digestを別の正しさpassで検査する。FP計算を省くstubが失敗するrunnerテストも追加。
   旧runner結果は保存し、新契約のschemaを分ける。
3. **操作分解:** syntax parse、互換前処理、prepared-molecule FP、binding変換を
   profilerで分ける。prepared測定と主目標を混ぜず、別runの中央値の差分からFP時間を
   推定しない。主比較では前処理を計時外へ逃がさない。cache済みSMILES/FPも禁止。
4. **入力:** 既存固定10k exposed集合を開発用に使いhash/入力順/重複を固定。
   小分子・鎖状・芳香族/縮合環・大環/多環・荷電/塩・同位体・大分子別にも報告する。
   malformed/unsupported/timeoutは別安全性laneだが全入力会計に残す。
   新規の未使用performance確認集合は最低2,000件とし、抽出元・seed・クラス別件数・
   対応域をPF0で事前固定。候補freeze後に一度評価し、由来と重複監査を保存。
   対応域のbit一致とrefusal会計も再検査し、不合格集合をtuningに使ったら確認用から退役させる。
   T1.6の8k sealed精度holdoutはprofilingやtuningに使わない。
5. **統計:** idle host・固定電源/温度条件、交互AB/BAの独立fresh-process対を最低20組、
   別の3セッションで実行。warm-upは両arm同じ固定件数とし計時外。
   一組は同じ全corpusの処理で、10k行を独立した反復数として扱わない。
   全raw値、中央値・p95 of process means、mol/s、paired speedupを保存。
   各対の`speedup = RDKit時間 / CheMatic時間`の幾何平均を主推定量にする。
   セッション層別のpaired bootstrapで95% CIを出し、反復数・bootstrap seed/回数・
   停止条件を事前固定。勝つまで反復・外れ値削除・再試行をしない。
6. **資源:** 一次測定はsingle-thread、候補同士も直列実行。起動/download、
   計時中のallocation、WASM pages、JS heap、process-tree RSSを別指標にする。
   profilerを付けたrunは公式速度値に使わない。100k batch/並列化は別laneで、
   index/search/Tanimotoの速度をParse＋FPの速度に読み替えない。

## 3. 実行順 — T3.6 / PF0–PF5

各段階は小さな変更単位で、正しさ→局所計測→同条件end-to-endの順に検証する。
工数は作業配分の目安で、速度達成やリリース日を保証しない。

| 順番 | 作業・目安 | 成果物 / 次へ進む条件 |
|---|---|---|
| PF0 | 測定契約・1.0.17再baseline（1–2日） | 上記manifest、equivalent-output runner、固定corpus、各armのbit/拒否照合、raw baseline。未検証出力ではprofiling以降へ進めない |
| PF1 | hotspotの分解（1–2日） | Rust/WASM別のprofile、allocation数/bytes、前処理・ring・展開・detail集計・変換の割合。支配コストとAmdahl上限から順序を再確認 |
| PF2 | bit-only経路（2–3日） | provenance集計を必要時のみ生成し、公開bit APIを共通kernelのbit sinkへ接続。detail/raw経路は維持。対応全行bit一致・同じ拒否と実測改善 |
| PF3 | 互換前処理の重複削減（2–4日） | normalize/aromaticity/ringの呼出し回数・cloneを計測してから共有。profile・分子revision・設定付きprepared状態とinvalidationsを検証。必要な知覚処理は省略しない |
| PF4 | 展開kernelとbinding（2–4日） | denseなatom/radius状態、隣接列・scratch再利用、環境集合のallocation/clone低減、直接byte出力を個別に比較。nativeも別laneで回帰。batchはscalar主目標と分離 |
| PF5 | 採用・再現・公表（2–3日＋外部host待ち） | 凍結candidate、独立性能確認集合、3browser/2hostの環境別scorecard、cross-binding/検索回帰、CI gateと再現手順。未達/欠測も公表 |

PF1で根拠が出た場合にPF3/PF4の内部順を変える。改善の小さい仮説に2日以上使う前に
profileと投資対効果を再確認する。実測改善がない変更は採用せず、失敗実験も記録する。
native ECFPはPF1から同時計測するが、別定義の高速nativeへのfallbackで主目標を達成しない。

### 2026-09-20 開発時点の進捗（採用判定前）

- **PF0:** runnerをschema 2へ更新し、両armが256 byte・LSB-firstの2048 bit出力を
  timed operation内で生成・消費する契約にした。固定10,000行でChrome、Firefox、WebKit
  の各armが対応9,999/9,999行で完全一致した。残る1行は既知の高配位Fe(II)を両比較から
  明示除外するtyped refusalであり、coverageを増やしたものではない。独立ChEMBL 5,000行は
  5,000/5,000行で完全一致した。
- **PF1:** 5,000行native診断では、bit-only detail全体の約87%相当を、複雑な環系の
  full SSSR構築が占めると判明した。Morgan connectivity invariantが必要とするのは
  `isRingAtom`だけであり、選ばれた最小環リストではないことを確認した。
- **PF2/PF4（実装済み・採用候補）:** bit-only APIはdetail/provenance mapを生成せず、
  WASMも直接packed bytesを返す。不要なnormalization clone、round状態clone、neighborの
  一時Vecを削減した。さらに互換Morgan経路だけは、eligible non-bridge bondの両端を
  iterative DFSで印付ける線形の環原子フラグを使う。公開SSSR APIとその選択結果は変えず、
  代表的な単環・融合環・ケージ・スピロ・分離環でSSSR membershipとの一致をテストする。
- **速度（同一Mac、fresh browser process、同一packed出力契約）:**
  - Chrome 10,000行/20反復: CheMatic **0.039839**、RDKit **0.105116 ms/mol**、
    paired median **2.63x**。
  - Firefox 10,000行/10反復: **0.314531** vs **0.647615 ms/mol**、**2.06x**。
  - WebKit 10,000行/10反復: **0.037204** vs **0.110011 ms/mol**、**2.96x**。
  - 独立ChEMBL 5,000行/20反復（Chrome）: **0.022650** vs **0.139910 ms/mol**、
    **6.17x**。
- **却下した仮説:** suppressed environmentのwinnerだけを保持する専用fast pathは
  10,000行で **0.144554 ms/mol**となり悪化したため、実装を戻した。profileで支配的だった
  ring-membershipだけを最適化し、SSSRの意味論を近似する変更は採用していない。

この進捗は一台の開発host上のcandidate測定であり、raw artifactは候補commitに記録した。
`scripts/check_parse_morgan_rdkit_speed_gate.py` は各recordでpaired log-speedupの95%下限を
算出し、Chrome **2.61x**、Firefox **2.04x**、WebKit **2.86x**、独立ChEMBL **6.09x** と
いずれも1.0を上回る。GitHub-hosted second hostでも
[`run 35483926941`](https://github.com/kent-tokyo/chematic/actions/runs/35483926941) が三browserで
成功し、artifactから再算出した95%下限はChrome **2.90x**、Firefox **2.37x**、WebKit
**3.65x**だった。候補sourceを隔離PythonとNode-WASMで再ビルドしたcross-binding gateも
全3組で5,000/5,000一致した。これによりcandidate sourceのPF5再現・binding条件は満たした。
公開packageを作る場合だけは、publish後のtarballで同じgateを再実行して初めてrelease固有の
優位主張に更新する。

### ソースから確認できる最初の候補（効果量は未測定）

- `crates/chematic-wasm/src/mol_fingerprints.rs::rdkit_ecfp4_bitvec` は
  `rdkit_morgan_ecfp4_experimental` の全detail結果を作り、fingerprintだけを返す。
  `BitVec2048 → BitVecN → bytes`変換も計測対象。
- `crates/chematic-fp/src/rdkit_morgan_ecfp4.rs` は常にsparse countsと
  raw/folded bitInfoのmapを構築する。PF2では出力sinkを分け、不要な集計だけを省く。
  環境抑制・hash・代表選択を別実装へ複製しない。
- `crates/chematic-fp/src/rdkit_morgan_hash.rs::expand_one_pass_with_chirality`
  にはroundごとのneighborhood clone、atomごとのneighbors収集、集合groupingがある。
  dense化やscratchは候補だが、重複環境抑制・代表atom・overflowの意味論を保存する。
- 互換前処理後に`find_sssr`を行う経路がある。前処理内との重複や費用をprofileで
  確認し、単なるring-membership代替で結果が同じと推測しない。

## 4. 採用ゲート

### 正しさ・安全性（速度より先）

- 開発集合の対応域で**全bit一致、既存成功の新規refusal=0、全入力会計100%**。
  元から未対応のFe(II)も消さず、両armのstatus/分母を残す。対応を広げた場合は
  事前固定の共通域速度と拡張域を別集計する。
- sparse counts/raw/folded bitInfo、環境代表・抑制、atom順/SMILES表記入替
  （atom-index付き出力は対応mapで照合）、
  芳香族/Kekulé、同位体・明示水素・halogen正規化を既存oracleと照合。
  変更が共有kernelに及ぶ場合は20のradius/幅設定とchiral経路も回帰する。
- Rust/Python/Node/WASMのbit・status一致、detail説明の保持、A3のk=1/10/100・
  threshold・ordered IDs・unrounded scoresの回帰。native既定bitは変更しない。
- 分子編集後のprepared状態の失効、設定混用拒否、異なる分子へのscratch再利用、
  空/巨大/高次数/不正入力を検査。panic、hang、無制限cache、unsafe導入を速度の代償にしない。
  性能用高速pathでも既存の計算量上限・typed errorを保持する。

### 速度・資源・公開範囲

- 段階目標は「現行比改善」→「同条件RDKitと同等」→**paired speedupの95% CI下限>1.0**。
  主corpusの全セッションで点推定も>1.0、process-mean p95もRDKit以下を求める。
  同じ設定の独立確認集合でも通過するまで「RDKit超え」としない。
- 各変更の採用は同条件source baselineへの有意な主経路改善と、既存supported scopeの
  維持が必要。化学クラス別p95、native、detail APIの5%超の退行は調査・採用停止。
  5%以内のノイズもraw値に残す。標本不足なら合格でなくinconclusive。
- 初期化、WASM raw/gzip、測定可能なpeak memoryはbaseline比5%超の増加で採用停止・
  設計見直し。既存T3の安全性/メモリ絶対予算は緩めない。測れないmemoryを0や合格にしない。
- PF5で二つ目のhostとFirefox/WebKitを別集計。全browser/hostが通るまでは
  成功した環境に限定した表現を使う。Node/Python/nativeの勝利は各laneの同じ出口で判断。
  Python laneは両方のPython API、native laneは両方のnative kernelを比較する。
- CIは通常PRでdeterministic correctness、固定専用runnerで性能回帰を判定する。
  共有CIの時間は参考値。公開packageの最終再測定前は「candidate結果」と明記する。

PF5のCI gateは [`scripts/check_parse_morgan_rdkit_speed_gate.py`](../scripts/check_parse_morgan_rdkit_speed_gate.py)
と [`parse-morgan-rdkitjs-gate.yml`](../.github/workflows/parse-morgan-rdkitjs-gate.yml) に固定する。
後者はPRでは`performance`ラベル時、またはmanual dispatch時だけ走らせる。各browserのraw
recordをartifactとして保持し、speedupの点推定ではなくpaired log-speedupの95%下限が1.0を
厳密に上回ることを要求する。

旧「さらに1.10x SMILES」目標は引き続き中止。この計画は2026-09-20の明示的な
Parse＋FP競争目標であり、その旧目標や全面的なRDKit精度優位を復活・達成したものではない。
速度未達でも正しさの改善は別途採用可能だが、性能目標はopenのままにする。

## 5. 既存資産・証跡・並行作業

runnerを増殖させず、次を拡張する。新しいCLI flagや合格artifactは未実装であり、
ここに挙げた計画を既存runnerが既に検証するとは扱わない。

- 測定: `scripts/bench_browser_wasm_vs_rdkit_isolated.py`。
  PF0でmanifest/同一出力/paired統計を追加し、実際の実行コマンドをraw recordに保存。
- bit gate: `scripts/check_browser_rdkit_ecfp4_parity.py`、
  `scripts/rdkit_ecfp4_cross_binding_parity.py` と既存sparse/bitInfo gate。
- kernel oracle: `scripts/ecfp_rdkit_morgan_ecfp4_parity.py`、
  `scripts/ecfp_rdkit_environment_parity.py`、`scripts/ecfp_rdkit_suppression_parity.py`。
- 出力: `validation/results/` に契約版・commit・環境・日付つきraw/summary、
  `benchmarks/` に再現コマンドと判定、benchmark indexにリンク。
  旧recordは上書きせず、最速runのみの抜粋をしない。

まずPF0→PF1を最優先の性能作業にする。T1.6凍結精度packetは候補とデータを分離して
並行維持する。silent corruption・panic・資源上限違反は常に優先修正。
他の大型機能・3D拡張よりPF2–PF5を優先するが、既存Trust RCの精度/安全性条件は免除しない。
