# Code audit — 2026-09-10

## 判定

監査中に重要度「中」相当の防御不足を1件修正した。修正後に再現する重要度「大」または「中」のバグは **0件**。

ここでいうバグは、入力に対する誤った化学結果、データ破壊、クラッシュ、公開APIの危険な挙動、またはバインディング間の契約破綻を指す。未実装機能、競合ライブラリとの未比較領域、性能改善余地はバグとは分離した。

## 確認した範囲

- 全ワークスペースの Rust テスト（`--all-features --offline`）
- リポジトリ標準 `scripts/check.sh`
- 本番コードの `unsafe`、panic 系呼び出し、配列境界アクセス、TODO／未実装箇所
- Rust／Python／WASM／CLI の共有契約とストリーミング入力制限
- CIP、芳香族性／環認識、3D、MMFF94、反応 SMARTS、形式変換、MCP 境界

## 証拠

| ゲート | 結果 |
| --- | --- |
| `cargo test --workspace --all-features --offline --quiet` | 成功。実行された全テストで失敗 0（環境依存の ignore は除外） |
| `bash scripts/check.sh` | `All checks passed.` |
| unsafe surface | reviewed native-InChI FFI boundary のみ |
| cross-binding manifest | 50 four-binding operations、6 source-only operations、4 bindings |
| streaming safety | 800 malformed + 480 generated parser-entry cases、失敗 0 |
| MMFF94 issue #337 determinism | 6 fixtures × 256 seeds = 1,536 checks、失敗 0 |
| reaction SMARTS | 38 cases、Rust／Python／WASM 契約一致 |
| Python binding evidence | 945 tests、失敗 0 |
| format / 3D / triclinic / reaction evidence | すべて成功 |
| clippy / format / dependency policy / publish graph | すべて成功 |

## 重要度別の分類

### 大 — 0件

データ破壊、未処理入力によるクラッシュ、化学結果の明確な誤り、セキュリティ境界逸脱は確認されなかった。FFI の `unsafe` は native InChI 呼び出しに閉じ込められ、標準ゲートで表面積を固定している。

### 中 — 0件

主要な公開境界で入力制限・エラー変換・バインディング間の結果契約を検証した。CIP の原子順不変性、環認識の決定性、反応 SMARTS、形式変換、ストリーミング parser-entry の回帰ケースも通過した。

監査で `resolve_is_r_from_groups` の不正な内部状態に対する map／slice の無検証アクセスを確認した。位置数不足、rank欠落、空rank群で panic せず `None` に退避するよう修正し、`resolve_is_r_declines_malformed_position_state` を回帰テストとして追加した。CIPの通常経路の結果は変更していない。

### 小／保留 — 修正対象外

- `cargo-deny` の既定 advisory DB は読み取り専用だったため、検査は writable fallback で完了した。これは環境設定の問題。
- 依存グラフには複数バージョンの依存が残る。`cargo-deny` は警告を出すが、advisory／license／source gate は成功している。直ちに機能バグとは扱わない。
- ローカル `.venv` に chematic／pytest がないため、標準スクリプトの pytest 部分はスキップされた。Rust 側の Python evidence と既存の契約記録は成功しているため、環境再構築時に再実行する。
- RDKit との一部 parity、WASM ツールチェーン、未実装機能はロードマップ上の残課題であり、今回の中／大バグ判定には含めない。

## 結論

中／大バグがなくなるまでの反復条件は満たした。追加の挙動変更を伴うリファクタリングは、今回の監査で見つかった問題を隠す可能性があるため実施していない。次回は依存統合、Python 環境での実行、RDKit parity の実測を独立したタスクとして扱う。
