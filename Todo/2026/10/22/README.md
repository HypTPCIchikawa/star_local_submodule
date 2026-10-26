# 2026-10-22：4HeのみのPID係数調整・Purity推定

**最新改訂：** [低R Gaussianのみ・中心固定（002版、11頁report）](report_he4_fixedmean_20261022.md)。R≥0.8とPID係数は前版保持。

[当日のTodo](STAR_Todo20261022.txt)。**作成・検証完了。** [実施結果とPDF索引](report_he4_20261022.md)。⁴He全199区間を評価し、図入り10頁report・Purity図・候補PIDを保存。

d・t・³HeのPID関数形と係数はGit元版（参照8f7d539）を維持する。⁴Heのみ従来関数形の係数を調整し、Purityと図入りPDF reportを作成する。

| R [GeV/c] | 試料 | ⁠⁴He Purity評価モデル |
|---|---|---|
| 0.610～0.80 | ALL、TOF要求なし | Gaussianのみ。中心は前版Gaussian＋linear値に固定 |
| 0.80～2.0 | ALL、TOF要求なし | Gaussianのみ、smooth背景なし |
| >2.0 | 既存TOF m²選択、dと同じ窓 | ⁠⁴He Gaussian＋linearのみ。他核種成分なし |

R=0.80は中心自由側に含める。詳細なbin境界・選択・初期値・固定条件は実施時に明記する。低Rの混入モデルと、高Rで他核種を除外する指定を区別する。

[前日の全区間Gaussianのみ試行・目視確認資料](../21/he4_r061_20261021.md)と[旧linear版の図入りreport](../21/gaussian_linear_report_20261021.md)は比較履歴として保持する。

**本番採用作業：** [確定⁴He係数・全核種TOF条件の反映](pid_adoption_20261022.md)。

本番採用完了：確定⁴He係数を対象mainconfへ反映し、全核種TOF境界R=2と³He k* controlを検証。STAR再build・14,147＋1,359チェック通過。[適用設定・検証](pid_adoption_20261022.md)。

[Git版からのヒストグラム追加・変更一覧](histogram_changes_vs_git_20261022.md)：³He背景control、⁴He PID/TOF QA、ペアPID QA・Sparseの用途と定義数を整理。

[Λ・⁴Heのβ–η／β–Yアクセプタンス追加](beta_acceptance_20261022.md)。

[⁴HeアクセプタンスのR・TOF別TH2](he4_acceptance_samples_20261022.md)。

⁴HeのR・TOF別TH2は実装・検証完了。低RはTOF不要／有効TOF／TOF+m²、高RはTOF+m²の4区分。[区分名と比較方法](he4_acceptance_samples_20261022.md)。

**全ヒストグラムの総合説明：** [Λ–原子核Femto解析のヒストグラム一覧（現行⁴He設定へ更新）](../../../../mdfiles/femto/implementation/femto_lambda_histograms_20260930.md)。既存分布を含む用途・候補条件、全1,787キー・型・軸、確定PID/TOF、βおよびR/TOF別アクセプタンスを収録。

[更新済み⁴He解析の少数実データ検証と本番再投入](production_relaunch_20261022.md)（作業中）。

⁴He更新版は実データ4,027入力で正常終了、全追加TH2充填を確認。本番再投入はSUMS `A76FAF721CDAB36920D53B9E8EEA9230`で全件job生成中（26,994予定、完了前）。[検証・投入状況](production_relaunch_20261022.md)。

⁴He全job終了後にyudai.ichikawa.d3@tohoku.ac.jpへ1通通知する監視を設定済み。Git pushはjob終了後まで保留。[通知・投入記録](production_relaunch_20261022.md)。

**本日の作業終了ログ：** [summary20261022.md](../../../../analysisnote/20261022/summary20261022.md)。投入・全job終了メール監視は継続。Git pushは解析終了後まで保留。

[⁴He改訂report](report_he4_fixedmean_20261022.md)はTOF m²選択後R–dE/dx＋確定PID曲線を6頁に追加した12頁版へ更新。Purity・fit結果は維持。

⁴He再投入処理終了（自動記録）：SUMS `A76FAF721CDAB36920D53B9E8EEA9230`、終了コード0。[検証・投入記録](production_relaunch_20261022.md)。farm解析の完了は別途確認。

**本番再投入A76は失敗：全26,994件が入力形式不一致でイベント処理前停止、ROOT0件。** [原因・修正・farm再検証](batch_input_fix_20261022.md)。Condor exit0は解析成功ではない。

A76失敗原因の修正後、farm 1jobで2000イベント・全1787ヒスト・ROOT回収を確認済み。修正版の全件再投入と終了通知を002として開始。[詳細](batch_input_fix_20261022.md)。

修正版の全件再投入SUMS ID：`04BF2DAE47164B6F040E3A460E539289`。投入進行中。終了通知もこの新IDを対象に設定。旧A76は失敗記録として保持。


**通知条件を更新：** 修正版SUMS `04BF2DAE47164B6F040E3A460E539289` は、全件解析・ROOT確認後にhaddし、結合結果検証・個別ROOT保存後にメールする。先の「全job終了時点で通知」はこの設定に置換。[実施記録](batch_input_fix_20261022.md#全件終了後の検査付きhaddとメールへ変更)。現在は投入継続・runmeta待ち。

⁴He再投入処理終了（自動記録）：SUMS `04BF2DAE47164B6F040E3A460E539289`、終了コード0。[検証・投入記録](production_relaunch_20261022.md)。farm解析の完了は別途確認。
