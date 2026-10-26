# 2026-10-21：10/17の方法を基準に、滑らかな従来PID関数で純度を再評価

**最新試行：** [全核種Gaussianのみ・全区間fit PDFと再評価](he4_r061_20261021.md)。

**前回採用（比較履歴）：** [全核種Gaussian＋linearのみのPurityと各核種3代表点を含む15頁PDF](gaussian_linear_report_20261021.md)。

[当日のTodo](STAR_Todo20261021.txt)は2026-10-20に計画として保存し、10/21のユーザー指示で実行した。
**再現・微調整・低R＋高R Purity・全R fit集・最終検証が完了。本番設定は変更していない。**

[今回の実施報告・全成果物](report_smooth_pid_20261021.md)を参照。

**最新の追加更新：** [4He R>0.61・低Rピーク位置の修正と更新Purity](he4_r061_20261021.md)。低R参照を補正し、0.70以上は前回の拡張範囲fitと数値一致。


## 最優先

1. [10/17の方法・結果](../17/README.md)を設定と実装から再現する。d・tを最初の回帰対照にする。
2. ガタガタした補間は使わず、元の`log10(R)`中心3次・幅2次式で係数を確認する。
3. 10/17のpurity推定手順を保って、ALL低R＋TOF高RのPurity vs Rを再評価する。

PID中心・幅のR依存曲線と、各スライスの信号応答・背景fitは別に点検する。
10/20のsplit-normal固定・GLS・追加採否maskを、そのまま再開の既定にしない。
今回の失敗だけで「多項式では推定できない」と結論しない。

## 参照資料

- [10/17の実施記録](../17/report_highrigidity_purity_20261017.md)：当時のモデル・TOF窓・集約方法・設定。
- [10/17の低R＋高R統合PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_tpc_tof_final_20261017_001/purity_tpc_tof_combined.pdf)：再現確認の比較基準。
- [10/20の結果](../20/summary_polynomial_all_purity_20261020.md)：変更点と失敗条件の履歴。明日の標準手法にはしない。
- [10/20の日報](../../../../analysisnote/20261020/summary20261020.md)：作業全体と終了時の方針。

本番はGit元PIDを保持する。新cutoffは未確定で、過去のROOT・図・設定は上書きしない。


完了：[10/17再現と滑らかなPID再評価](report_smooth_pid_20261021.md)。候補16,640スライス、元係数比較3,400スライス、代表24点の再現、PDF全1,610頁の読戻しがPASS。

高RのR=2からの連続支持域はd<3.05、t<3.40、3He<3.55、4He<3.20 GeV/c。域内の最小モデル下側Purityは99.13%、98.60%、97.90%、98.39%。未評価を低純度とせず、穴をつなぐ上限や本番vetoは作成していない。

追加表示修正完了：[linear/logの横軸幅を約2倍にした新版fit集・全8冊](fit_atlas_wide_20261021.md)。全1,600頁の読戻しがPASS。数値結果・Purity・001版PDFは保持。

Purity図の追加表示修正完了：[1ページ目のt/3He/4Heの横軸を点の範囲へ調整した最新版PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/purity_axis_002/purity_vs_p_over_q.pdf)。d・2ページ目・数値は維持。

- [背景なしfit・3He低側半値までのfit範囲拡張（完了）](no_smooth_he3window_20261021.md)

- [4He R>0.61・指定ピーク位置での再評価（完了）](he4_r061_20261021.md)

- [Gaussian＋linearのみ・図入りPDFレポート（完了）](gaussian_linear_report_20261021.md)
