# 低R・全native区間の1D fit図（保存fitの再表示）

注記：この001は品質注記導入前の最初のd R=0.40–0.60（80頁）だけの予備試験として保持する。
全冊は[品質注記付き002正式版](../fit_atlas_native_20261017_002/README.md)へ移行した。
以下の全冊計画は001準備時の記録で、001で32冊全てが完成したことを意味しない。

2026-10-17作成。3.85 GeV・2021 P24iyの2026-10-15研究で保存した
ALLスペクトルと予測countだけを再表示する。新fit、イベント再解析、PID係数変更は行わない。
対象はR=|p/q|=0.40–2.00 GeV/c、幅0.01の160区間×4核種×4名目モデル。
下表の各冊はR幅0.20、20 native区間、80ページで、全32冊・2560ページ。

各native区間の4ページは、Gaussian／log-normalとlinear／exponential背景の4モデルを
2モデルずつ表示し、全dE/dx（0–120 keV/cm）とfit範囲をそれぞれlinear／logで示す。
dE/dxはnative幅0.1 keV/cmのままで、正規化・表示rebin・再fitはしていない。
ALLの生count/error/flowが4モデル間で完全一致することをrendererで確認する。

## 分冊

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

描画された成分が物理的に支持されたことは意味しない。
数値未収束・応答支持なし・採用外の区間も省略せず、診断用に含める。
各native区間・モデルの`numericalValid`、`fitFailure`、`atBoundary`、
`fitResponseQuality`、`targetFreeSupported`は[coverage.tsv](coverage.tsv)に記録する。
Purityの採用判定は[元の比較PDF](../views_final_001/purity_cutoff_overview.pdf)と
[採用範囲表](../adopted_ranges_001/adopted_ranges.tsv)を参照する。
この分冊は低RのTPC-only研究で、高RのTOF要求後purityとは別の観測試料である。

## 再現・完了判定

リポジトリ直下で、明示した新しい出力先を用いる。

```bash
python3 -B tools/nuclear_pid_cutoff/render_native_atlas.py prepare \
  --settings tools/nuclear_pid_cutoff/native_atlas_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/fit_atlas_NATIVE_NEW
python3 -B tools/nuclear_pid_cutoff/render_native_atlas.py run \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/fit_atlas_NATIVE_NEW \
  --book d_R040_060
```

最初の1冊を確認後、`--book all`で残りを1workerで作成する。
既に完了した冊はSHA一致を確認して再描画しない。途中の不完全なPDFは上書きしない。
全冊後に同driverの`check`を実行すると、既存fit入力の不変性と全冊のSHA・頁数を確認し
`completed.json`を新規作成する。未完了時にはこのファイルは存在しない。

- [plan.json](plan.json): 明示YAML・32冊の範囲と期待頁数・入力ROOT/JSONのSHA。
- `source/`: driver・renderer・wrapper・mainconf・表示YAMLのsnapshot。
- `*.tsv`: 実際に描画するnative index・モデルROOTパス・表示設定。
- `*.command.json` / `*.log`: ROOT5/STARでの実行記録。
- `*.readback.log` / `*.completed.json`: 閉じたPDFのGhostscript読み戻し・80頁確認・SHA。

未収束点を除外して図を良く見せたり、名目以外のfitへ差し替えたりしていない。
native幅による射影の広がりは元studyのままで、自由fit幅をintrinsic detector幅と断定しない。
