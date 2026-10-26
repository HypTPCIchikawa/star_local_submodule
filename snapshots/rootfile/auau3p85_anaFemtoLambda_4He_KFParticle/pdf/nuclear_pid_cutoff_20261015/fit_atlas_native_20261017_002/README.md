# 全低R区間の1D fit図：品質注記付き正式版

追記：この002はHe4 R=0.80–1.00（80頁）だけの表示試験として保持する。
status注記が軸指数・上枠へ重なったため、[003正式版](../fit_atlas_native_20261017_003/README.md)でヘッダー余白を調整した。
下表の全冊は準備時の計画で、002で全冊が完成したことは意味しない。

2026-10-17。3.85 GeV・2021 P24iyの保存ALL fitを再表示したもの。
新fit・正規化・イベント再解析・PID変更はしない。
R=|p/q|=0.40–2.00 GeV/cの160 native区間（幅0.01）を、4核種×4名目モデルで全て載せる。
各冊はR幅0.20の20区間・80ページ。全32冊・2560ページ。

1区間当たり4ページでGaussian／log-normal×linear／exponential背景を2モデルずつ示し、
全dE/dx範囲（0–120 keV/cm）とfit範囲の両方をlinear／logで表示する。
dE/dxはnative幅0.1 keV/cmで、表示rebinはしない。
全モデル間の生count・誤差・flow一致を各native区間で検査する。

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

## 品質注記の意味

各パネルの`fit: valid/FAILED`は数値fitの判定、`response: pass/FAIL`は保存応答品質判定。
`FAILED`や`FAIL`も省略せず診断用に描画する。線があるだけで核種帰属やpurityを支持するものではない。
未支持をpurity=0へ置換していない。名目4モデルだけで、別のfitへ差し替えない。

[coverage.tsv](coverage.tsv)は全2560 model/native区間のfit状態・境界・応答支持を記録する。
最終的な採用範囲は[採用表](../adopted_ranges_001/adopted_ranges.tsv)、
purityの判定は[元の比較PDF](../views_final_001/purity_cutoff_overview.pdf)を参照する。
本分冊はR<2のTPC-only解析で、R≥2のTOF要求後のpurity推定とは別の試料である。

## 再現性と確認

- [plan.json](plan.json): 明示YAML、全分冊・頁数、保存ROOT/JSONのSHA。
- `source/`: ROOT5 renderer・再利用元renderer・driver・wrapper・mainconf・YAMLのsnapshot。
- `quality_*.tsv`、`models_*.tsv`、各冊TSV: 実際に使った品質・保存モデル・native index。
- `*.command.json` / `*.log`: batch-matched SL24y／ROOT5実行の来歴。
- `*.readback.log` / `*.completed.json`: 閉じたPDFのGhostscript読み戻し、80頁確認、SHA、サイズ、時間。
- 全冊確認後の`completed.json`: 32冊・2560頁、保存fit入力の不変性、source/table hash。

最初のHe4 R=0.80–1.00を試験し、R=0.855のGaussian+exponential数値未収束表示を目視確認する。
その後は同時描画2worker以下で進める。未完成PDF・既存成果物は上書きしない。
[001予備版](../fit_atlas_native_20261017_001/README.md)は品質注記導入前のd R=0.40–0.60だけの試験として保持する。

再現用driverは`tools/nuclear_pid_cutoff/render_native_atlas.py`。
`prepare --settings tools/nuclear_pid_cutoff/native_atlas_settings.yaml --output-dir NEW_DIR`で新規準備し、
`run --output-dir NEW_DIR --book he4_R080_100`で初冊、`--book all`で全冊を作る。
最後に`check --output-dir NEW_DIR`で総合確認する。
途中で完了した冊はPDFのSHA一致を検査してskipし、未完了冊を無断上書きしない。
