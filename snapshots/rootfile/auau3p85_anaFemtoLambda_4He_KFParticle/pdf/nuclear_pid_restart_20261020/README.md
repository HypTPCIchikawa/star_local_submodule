# 2026-10-20：Git元PIDから再開した核PID study

対象はAu+Au 3.85 GeV、2021 P24iyの保存済みraw track QA。
過去のfit・校正点・cutoffを初期値として使わず、Git `8f7d539` の中心3次・幅2次多項式から再開した。
元ROOTは変更せず、再解析・job投入・haddは行っていない。

**今回の結果は[現在の応答・fit範囲比較のまとめ](../../../../Todo/2026/10/20/summary_current_responses_20261020.md)を参照。**
応答をさらに増やさず、ここまでの比較を整理した。新係数・純度・cutoffは本番採用していない。

## 現在の出力

| Directory | 内容 | 判定 |
|---|---|---|
| `tof_tag_smoke_001` | 初回準備のパス処理エラー時の途中snapshot | 未実行・不採用。履歴保持 |
| `tof_tag_smoke_002` | 4核種×R=0.875／1.775のGaussian試験とGit曲線COLZ | 計算・保存・PDF検査完了。物理採用ではない |
| `tof_tag_survey_001` | 各核種95点、計380のGaussian診断 | 190ページPDF・readback完了。純度・係数・cutoffは未認定 |
| `tof_tag_local15_smoke_001` | 同じ8点でfit範囲だけ±2.5σ→±1.5σ | 4ページPDF。範囲を狭めても不一致あり |
| `response_gaussian_smoke_001` | 新しい応答比較エンジンのGaussian基準 | 8行保存・readback完了 |
| `response_{gaussian,left_a1,left_a25,right_a1,right_a25}_smoke_002` | Gaussian／左右tail、初期αを2通り比較 | 40 fit保存。モデル採用ではない |
| `response_comparison_001` | 最初の3モデル比較 | 履歴保持。表示は002で改善 |
| `response_comparison_002` | 同じ保存fitを横長PDFと短い注記へ変更 | 履歴保持。α識別性の追加診断前 |
| `response_comparison_003` | α未拘束の解を別判定した最新比較 | 8ページ・48パネル、データ同一性／PDF検査成功 |
| `native_storage_export_001` | 代表8スライスの全9616 bin count／errorの書出し | 元ROOTとの完全一致・再読込成功 |
| `pico_streamer_001` | 同productionの実PicoDst 2本の保存形式 | Float16既定12bit。全入力のschema監査ではない |
| `storage_math_001` | ROOT5保存処理との量子化モデル照合 | 57075項目成功。実データfitの採用ではない |
| `storage_default12_smoke_001` | 同じ8スライスで予測積分に量子化を反映 | count不変、ROOT／JSON／bin表のreadbackと4ページPDF検証成功 |
| `storage_recenter25_smoke_001` | 実ピーク中心±2.5元σ、Gaussian＋量子化 | 8 fitと4ページPDFを保存・検証 |
| `storage_recenter15_smoke_001` | 実ピーク中心±1.5元σ、Gaussian＋量子化 | 8 fitと4ページPDFを保存・検証 |
| `split_normal_smoke_001` | 共通mode・左右幅を分ける応答、Git ROI | 8 fitと4ページPDF、人工試験6773項目成功 |
| `split_normal_survey_001` | 同じ応答、明示抽出95点×4核種、Git ROI | 380 fit、bin対応16654行、190ページPDFを保存・検証 |
| `split_recenter15_smoke_001` | 同じ応答、実ピーク中心±1.5元σ | 8 fitと4ページPDFを保存・検証。全Rへ一般化した結果ではない |
| `native_storage_export_survey_001` | 380スライスのnative count／error全456760行 | 元ROOTとの完全一致・再読込成功。Local ROI全surveyは未実施 |

すべてのfit設定、native R bin、dE/dx範囲、入力・実行ソースSHAは各 `plan.json` と `source/` に固定。
`fit.log.completed.json` は技術的な計算・保存の完了であり、fit品質・純度・PIDの採用記録ではない。
各描画先の `summary.json`／`completed.json`（ツールにより異なる）と `provenance.json` にPDF検査と入力不変性を記録する。
`quality_002`はtarget modeがfit範囲内かも検査した後処理の表。元fitのフラグを変更しない。

## 最初に見る図

- [Git元曲線とraw ALL／TOFのCOLZ（8ページ）](tof_tag_smoke_002/git_pid_baseline.pdf)
- [Gaussianの代表8 fit、linear／log（4ページ）](tof_tag_smoke_002/views_001/tag_peak_fit_diagnostics.pdf)
- [Gaussian全380 fit、linear／log（190ページ）](tof_tag_survey_001/views_001/tag_peak_fit_diagnostics.pdf)
- [狭い局所fitとの比較（4ページ）](tof_tag_local15_smoke_001/views_001/tag_peak_fit_diagnostics.pdf)
- [収束・χ²・統計・初期値→fit値の集計](tof_tag_survey_001/diagnostic_quality_summary.json)
- [最新：Gaussian／左tail／右tailの同一ページ比較](response_comparison_003/response_comparison.pdf)
- [最新：モデル・seed差・α識別性・χ²・誤差の表](response_comparison_003/comparison.tsv)
- [Float16保存精度を扱ったGaussian診断（4ページ）](storage_default12_smoke_001/views_001/tag_peak_fit_diagnostics.pdf)
- [Local ROI＋非対称応答の代表8点（4ページ）](split_recenter15_smoke_001/views_001/tag_peak_fit_diagnostics.pdf)
- [非対称応答の380 fit集（190ページ）](split_normal_survey_001/views_001/tag_peak_fit_diagnostics.pdf)
- [非対称応答surveyの最終品質表](split_normal_survey_001/quality_002/quality.tsv)

## 入力と解釈の注意

- raw cacheの `ALL`、`TAG_D`、`TAG_T`、`TAG_HE3` のみを使用。旧fitの係数等は使用しない。
- ⁴Heとdは全く同じm²窓 `3.268773 < m² < 3.693147` の `TAG_D` を共有する。
  この窓は純粋な⁴Heタグではない。tは `7.2564275 < m² < 8.2816925`、
  ³Heは `1.78285 < m² < 2.06485`。
- raw QAは電荷符号・centralityを積分したtrack分布。候補数をevent数や信号効率と呼ばない。
- 最初のsurvey下端R=0.405は全4核種の元式が定義される範囲を調べたもの。
  新しいPID下限ではなく、d・tのさらに低Rは別途調べる。
- 380行中376行は数値flagが良好だが、χ²による科学的採用ではない。
  125行は他核種の元中心がfit範囲内にあり、対象1ピークのまま純度fitへ転用しない。
- 高χ²や範囲依存があるため、Gaussian＋指数tailの比較をユーザーと合意して実施した。
  tailの幅はGaussian core幅であり、全応答のRMSとは区別する。
- tailが全native-R nodeでfit範囲の外にある場合、有限のHesse誤差が出てもαを決められない。
  003ではこの識別性を保存fitから後処理で判定し、元fitを変えず比較判定へ追加した。
- 3応答でも低Rの高χ²が残るため、新係数・純度・cutoffはまだ認定していない。
- Float16の保存丸めは予測側の積分にだけ反映し、raw値・bin幅は変更していない。
  細かい周期構造と整合するが、残るピーク不一致を全て解消したわけではない。
- Local ROI＋非対称応答で小さなχ²となる³He・⁴Heの代表点があるが、全域の校正・純度を保証しない。
  背景と信号の肩・中心・幅の相関は、ALL fitと範囲依存も含めて確認する。
- χ²≈1でもtarget yield=0となるGaussian fitを発見した。校正点や低純度・vetoの根拠にしない。
- 全surveyのLocal ROI追試と多項式校正は次手順。新しい準備拡張を保存・検証済みとはしていない。

最新の検討・未完了部分は [Todoの作業レポート](../../../../Todo/2026/10/20/report_pid_restart_20261020.md)、
再現手順は [studyツール](../../../../tools/nuclear_pid_restart_20261020/README.md) を参照。
