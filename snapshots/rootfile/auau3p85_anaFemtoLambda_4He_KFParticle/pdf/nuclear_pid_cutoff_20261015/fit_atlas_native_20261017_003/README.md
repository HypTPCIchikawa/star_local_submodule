# 全低R区間の1D fit図（ヘッダー調整済み正式版）

2026-10-17。3.85 GeV・2021 P24iyの保存ALL fitを再表示する。
R=|p/q|=0.40–2.00 GeV/cの160 native区間（幅0.01）×4核種×4名目モデルを全て含む。
各冊はR幅0.20・20区間・80ページ、全32冊・2560ページ。
新fit、正規化、イベント再解析、PID変更は行わない。

**全32冊・2560ページの生成と検査は完了しています。** PDF合計は619,296,231 bytes（590.6 MiB）。各核種8冊・640ページです。[completed.json](completed.json)に全冊の閉PDF読み戻し・SHA・保存fit入力不変性を記録しました。

各区間4ページで、Gaussian／log-normalとlinear／exponential背景の4モデルを2つずつ表示する。
全dE/dx範囲（0–120 keV/cm）とfit範囲の両方をlinear／logで示す。
dE/dxはnative幅0.1 keV/cmのまま、表示rebinなし。
4モデル間の生count・誤差・flowの一致を各区間で検査する。

| R [GeV/c] | d | t | ³He | ⁴He |
|---|---|---|---|---|
| 0.40–0.60 | [PDF](d_R040_060.pdf) | [PDF](t_R040_060.pdf) | [PDF](he3_R040_060.pdf) | [PDF](he4_R040_060.pdf) |
| 0.60–0.80 | [PDF](d_R060_080.pdf) | [PDF](t_R060_080.pdf) | [PDF](he3_R060_080.pdf) | [PDF](he4_R060_080.pdf) |
| 0.80–1.00 | [PDF](d_R080_100.pdf) | [PDF](t_R080_100.pdf) | [PDF](he3_R080_100.pdf) | [PDF](he4_R080_100.pdf) |
| 1.00–1.20 | [PDF](d_R100_120.pdf) | [PDF](t_R100_120.pdf) | [PDF](he3_R100_120.pdf) | [PDF](he4_R100_120.pdf) |
| 1.20–1.40 | [PDF](d_R120_140.pdf) | [PDF](t_R120_140.pdf) | [PDF](he3_R120_140.pdf) | [PDF](he4_R120_140.pdf) |
| 1.40–1.60 | [PDF](d_R140_160.pdf) | [PDF](t_R140_160.pdf) | [PDF](he3_R140_160.pdf) | [PDF](he4_R140_160.pdf) |
| 1.60–1.80 | [PDF](d_R160_180.pdf) | [PDF](t_R160_180.pdf) | [PDF](he3_R160_180.pdf) | [PDF](he4_R160_180.pdf) |
| 1.80–2.00 | [PDF](d_R180_200.pdf) | [PDF](t_R180_200.pdf) | [PDF](he3_R180_200.pdf) | [PDF](he4_R180_200.pdf) |

## 品質注記・限界

各パネルに対象核種・R・モデル名、数値fitの`valid/FAILED`、応答品質の`pass/FAIL`を表示する。
ヘッダー用の余白を確保し、軸の指数や上枠との重なりを避けた。
未収束・未支持点も診断として省略せず載せ、purity=0へ置換しない。
成分の線が描かれていることと、核種帰属やpurityが支持されることは別である。
数値fit・応答品質がpassでも、それだけで最終採用や完全な系統誤差評価を意味しない。

[coverage.tsv](coverage.tsv)に全2560 model/native区間の状態・境界・支持判定を記録した。
採用範囲は[採用表](../adopted_ranges_001/adopted_ranges.tsv)、
purity判定は[元の比較PDF](../views_final_001/purity_cutoff_overview.pdf)を参照する。
本分冊はR<2のTPC-only研究で、高RのTOF要求後purityとは異なる観測試料である。

## 来歴と検査

- [plan.json](plan.json): 明示YAML、全分冊・期待頁数、保存ROOT/JSONのSHA。
- `source/`: driver・新renderer・再利用元renderer・wrapper・mainconf・YAML。
- `quality_*.tsv`、`models_*.tsv`、各冊TSV: 実際の表示入力。
- `*.command.json` / `*.log`: STAR SL24y／ROOT5・batch-matched環境での実行記録。
- `*.readback.log` / `*.completed.json`: 閉じたPDFのGhostscript読戻し・80頁・SHA・サイズ・所要時間。
- 全冊後の`completed.json`: 全32冊・2560頁と保存fit入力の不変性・source/table hash。

初冊He4 R=0.80–1.00のR=0.855 Gaussian+exponential未収束例を目視確認してから、
同時描画2worker以下で残りを作成した。各冊を独立したROOTプロセスで閉じる。
未完成ファイルや過去出力は上書きせず、[001](../fit_atlas_native_20261017_001/README.md)・
[002](../fit_atlas_native_20261017_002/README.md)の予備試験PDFも保持する。

再現driverは`tools/nuclear_pid_cutoff/render_native_atlas.py`、
表示値は`tools/nuclear_pid_cutoff/native_atlas_settings.yaml`に明示。
新規ディレクトリへ`prepare --settings YAML --output-dir NEW_DIR`、
次に`run --output-dir NEW_DIR --book he4_R080_100`で初冊、
`run --output-dir NEW_DIR --book all`で全冊（1worker）、最後に`check --output-dir NEW_DIR`を実行する。
完了済み冊はPDFのSHA一致を確認してskipする。

本実行は追加の`tools/nuclear_pid_cutoff/run_native_atlas_parallel.py`から最大2workerで起動し、全冊成功後の`check`まで終了しました。同launcherと起動引数も記録しています。
