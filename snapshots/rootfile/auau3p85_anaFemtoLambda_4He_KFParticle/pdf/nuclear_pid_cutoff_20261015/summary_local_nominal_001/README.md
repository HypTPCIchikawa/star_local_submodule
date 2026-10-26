# 集約判定の旧版（使用しない）

このsummaryの支持判定は、4核種×4モデルの全16実行数を、各核種の4モデル数と
誤って比較していたため、全区間が`unassessed`となっていた。
fit値・モデルpurity・元ROOTは変更していない。

修正済みの正式集約は[summary_local_nominal_002](../summary_local_nominal_002/summary.json)。
この旧版は問題の記録として保持するが、cutoff判断や正式図には使用しない。
