# 2026-10-17 高R TOF選択後purityと全p/q fit集

依頼を [STAR_Todo20261017.txt](STAR_Todo20261017.txt) に保存しました。日付は日本時間です。

- [実施計画](../15/plan_highrigidity_tof_purity_20261017.md)：実際のm²窓・選択順序を監査し、低R結果と区別して統合。
- [実施状況・結果](report_highrigidity_purity_20261017.md)：入力、最終purity、採用区間、新旧PID比較、本番適用・検証。
- [前回の低R結果・設定検証](../15/README.md)：モデル低側50%以上の範囲、独立設定、最終PDF。
- [本日の日報](../../../../analysisnote/20261017/summary20261017.md)：作業全体、作成資料、検証結果、次回の残作業。

低Rの全native fit集は[32分冊・2560ページとして完成](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/fit_atlas_native_20261017_003/README.md)しました。低Rの80実行と既存PDFは保持しています。高Rの最終表も確定しました。

- [低R＋TOF後高Rの統合PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_tpc_tof_final_20261017_001/purity_tpc_tof_combined.pdf)：新3ページ＋元の低R3ページ、計6ページ。
- [本番採用区間と除外理由](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/production_ranges_20261017_001/README.md)：低Rの凍結範囲＋モデル検証済み高R。
- [共通本番PID設定](../../../../config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_production_20261017.yaml)：5入口から参照。各入口にprePurity20261017の戻し用設定を用意。
- [高R全区間fit集](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/highrigidity_tof_20261017_001/fit_atlas_001/README.md)：24冊840ページ。
- [新旧PID・実fit比較PDF（23ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/pid_parameter_comparison_final_20261017_001/raster_report_001/new_old_pid_final.pdf)：中心・幅・差・比・TOF窓を比較。
- [全5入口の20event試験](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/production_purity50_20261017_001/smokes_001/README.md)：閉ROOTとmetadataの読戻しPASS。

低R／高R評価、本番cutoff、追加のk*–核PID TH2／PID前対照／8D sparseまで実装・検証を完了しました。

- [新ヒストの意味・軸・射影](../../../../mdfiles/femto/implementation/femto_lambda_histograms_20260930.md)。
- [QA追加の検証索引](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/nuclear_pid_pair_qa_20261017_001/README.md)：実20event off/on・合成閉包・非ゼロmixing pool・容量stress・toy hadd。既存大統合ROOTは再解析していません。

追加指令の新旧dE/dx中心・幅・選択帯とTOF窓の比較、および本番でのpurity 0.5 cutoff適用も完了しました。ユーザー確認済みの方針に従い、d・t・³He・⁴Heの本番入口を共通校正・cutoffへ揃え、旧専用入口との差の記録と入口ごとの検証を行いました。本番jobは未投入で、次回の大規模統合前にはSparse監査の追加と実ファイル規模の容量確認が必要です。
