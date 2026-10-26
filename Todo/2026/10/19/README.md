# 2026-10-19 核PID関数形と高運動量上限の確認

**終了時の方針変更：以後は[10/20の新Todo](../20/README.md)を優先します。**
Gitの元PIDへ本番関数を戻し、m²選択での校正から仕切り直す計画です。
過去の補間校正・fit採否・cutoffは新studyの前提にせず、このページ以下は履歴として保持します。
本日は計画と[日報](../../../../analysisnote/20261019/summary20261019.md)の保存までで、Git式への本番復帰はまだ実行していません。

[STAR_Todo20261019.txt](STAR_Todo20261019.txt)に本日の依頼を記録しました。

1. 元の対数多項式を基準に、PID可能な支持範囲内で中心・幅を更新する。低純度の低Rまで無理に延ばさない。
2. 全核種でR<2 GeV/cはTPC-only、R≥2はTOF m²併用。⁴Heはdと同じm²窓を使う。
3. 純度P=S/(S+B)≥50%を基準に、核種ごとのp/q vs purityを再評価する。
4. 局所fit図・残差からχ²基準を提案し、不良fitを系統差から除外する。
5. 飛び地の除外は当該・両隣のfit図を示して再検討する。未評価を低純度と取り違えない。
6. dE/dxの縦軸と局所表示を見直し、linear／logの両方でピークを見せる。

状態：独立study設定・切り戻しと診断ツールを作成し、22試験PASS。
保存fit100実行・14840行の品質監査が完了。48頁の診断図は全頁readback済みで、
veto区間前後を含む追加版68頁も全頁readback済みです。本番PID切替、新gateの純度再fit、χ²閾値採用は未実施です。
固定Gaussian core＋tailの48実行64fit比較は完了。追加のfree背景pilotは新指示により未実行で保留。
背景³He・tの中心・幅を既知のPID多項式で固定し、実peakを含む局所fitの整合を優先します。

元PID関数固定の4モデル・28fitと比較PDF28頁も完了しました。
[最新：He・tの実peakと固定PID曲線の比較](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/peak_anchored_views_001/peak_anchored_diagnostics.pdf)
は各RでHe領域／t領域をlinear・log表示します。背景核種の中心・幅は動かしません。
25試験・ROOT/PDF再読込み・入力SHA検査済みですが、元曲線とpeakのずれでχ²は不良です。
純度・cutoff採用には使わず、次は各核種の実peakを別校正して多項式係数を検証します。
[固定PID比較の再現手順](../../../../tools/nuclear_pid_cutoff/peak_anchored_20261019/README.md)を追加しました。

注意：このPDFは元多項式のみの参照studyであり、現在のanaFemto本番入口とは同一ではありません。
現在の³He・⁴Heは10/13の校正点の線形補間が有効です。元多項式への本番切替は未実施です。
R=.875の数値照合、元ROOT生成時との区別、20–28付近がt確認ページであることを
[レポート冒頭](report_pid_quality_20261019.md)へ追記しました。コード・設定・PDFは変更していません。
Git HEADとの照合も追記しました。元多項式ファイルは完全一致ですが、10/13・10/17の
未コミットの校正導入・参照変更で実効PIDは変わっています。今回rollbackはしていません。

追加指示の左linear表示を修正しました。全linearが核種近傍、右logは広域／局所を維持。
[最新⁴He図（004・28頁）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/representative_views_004/fit_diagnostics_he4.pdf)
を参照してください。003は保持し、全68頁readbackと表示回帰11試験が通っています。

[実施報告・診断PDFへのリンク](report_pid_quality_20261019.md)に、旧vetoの未評価理由と次の比較をまとめました。

- [旧fit品質監査・描画の再現手順](../../../../tools/nuclear_pid_cutoff/quality_20261019/README.md)
- [元多項式study設定検査](../../../../tools/nuclear_pid_cutoff/logpoly_20261019_config_audit.py)

承認障害時は[AI運用手順](../../../../docs/ai/skills/approval-recovery.md)に従い、未実行の操作を明示し、手動承認への切替を速やかに案内します。承認設定の変更はユーザーが行います。
今回のfit・PDFの実行承認と検査が済んだ時点で、自動承認へ戻してよいとユーザーへ連絡しました。
[実施計画](implementation_plan.md)と[10月のTodo索引](../README.md)を参照してください。

10/14・10/15の記録は、2026-10-19の整理で同じ`Todo/2026/10/`配下へ移動しました。
旧設定でのpurity50判定は、PID窓更新後の認証として引き継ぎません。

- [対象PDF（p.5）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_tpc_tof_final_20261017_001/purity_tpc_tof_combined.pdf)
- [元の多項式](../../../../StMaker/common/NuclearIdDeDxVsMom.h)
- [前回の実施報告](../17/report_highrigidity_purity_20261017.md)
- [前回の日報](../../../../analysisnote/20261017/summary20261017.md)
