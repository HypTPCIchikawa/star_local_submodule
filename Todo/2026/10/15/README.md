# 2026-10-15 原子核PIDの使用可能rigidity下限

依頼原文は [STAR_Todo_20261015.txt](STAR_Todo_20261015.txt) です。2026-10-19の整理依頼により、旧`Todo/20161015/`から`Todo/2026/10/15/`へ移動しました。

- [実施計画・判断基準](plan_nuclear_pid_cutoff_20261015.md)：TPC-onlyでの核種別purity下限、全4核種の高rigidity TOF条件、検証方法。
- [実施状況・結果](report_nuclear_pid_cutoff_20261015.md)：80条件の集計、核種別使用候補範囲、最終設定・図と検証完了記録。
- [2026-10-15の終了・引継ぎログ](../../../../analysisnote/20261015/summary20261015.md)：完了した内容と明日の再開手順。

前日の低R応答関数の全域最適化を続けるのではなく、CFに使える範囲を先に決めます。過去のfit・設定は履歴として保持し、新しい条件は独立した設定として検証します。

**終了時点：** fit・感度集計と下限／除外域のJSON保存は完了しました。最終範囲のYAML作成、同条件の実イベント検証、最終図は明日へ持ち越しです。本番設定は未変更です。

**2026-10-17（日本時間）の再開分は完了：** [最終条件の独立mainconf](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_purity50_20261015.yaml)と[核PID YAML](../../../../config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_20261015.yaml)を作成し、batch一致build・先頭20イベント・閉ROOT検査・domain 1359項目・Python 17件を確認しました。[検証記録](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/purity50_20261015/README.md)、[最終PDF一覧](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_final_001/README.md)を参照してください。既存設定・校正係数は保持し、本番mainconfの切替やjob投入はしていません。

**追加依頼を実施中：** [高Rのm²要求後purityと低R結果の統合計画](plan_highrigidity_tof_purity_20261017.md)。保存済みヒストグラムが実際の高R TOF窓に対応するかをまず監査し、低R TPC-onlyと区別した統合版を作成します。

追加の「各p/qの1次元fit集」も含む今日の依頼・進捗は [Todo/2026/10/17](../17/README.md) にまとめています。
