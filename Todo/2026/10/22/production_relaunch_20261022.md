# 更新済み⁴He解析の少数実データ確認と本番再投入

2026-10-22。ユーザー指示：更新済みanaFemtoLambda_4He_KFParticleを少ない統計で実行し、正常動作とヒストグラム充填を確認後、jobを再投入する。今回の指示に検証後の本番投入まで含まれる。

作業中：2021 P24iyの明示入力、STAR再build、少数実データ、閉ROOTの正常終了・既存/追加分布・TOF区分/軸/全体一致を検証する。その後mainconfからjoblistを再生成し、2021 P24iy catalog・全イベント・5ファイル/jobと設定bundleを確認して投入する。

入口：config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml。監査先：rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/production_relaunch_20261022_001/。旧ROOT・投入済みsnapshotは保持する。

## 少数実データ検証：通過

実入力1本は4027 Entries。最初の10000指定はavailable不足を検出して処理前に停止（黙って短縮しないガードが動作）。その後同じ入力を全件指定(-1)で実行し、4027入力・3858再構成・選択Λ62、signalアクセプタンスΛ53、最終4He250を確認。2021 run/event/candidate ledgerと終了completedが一致。

既存・追加を含む全1787ヒストグラムの型・軸とROOT読戻しが通過。追加12TH2すべて非ゼロ。区分別4Heは低R TOF不要247、有効TOF73、TOF+m²64、高R TOF+m²3。全bin/flowで全体=低R TOF不要+高R TOF+m²、低R包含関係、Lambda/核のβとpTアクセプタンス件数、SE/MEのsignal上下分割の一致を検証。4He signal SE9/ME251、3He signal SE29も非ゼロ。jointSparse346962 entries。

現行mainconfから本番joblistを再生成し、catalogの2021/P24iy・nFiles=all・maxEvents=-1・5ファイル/jobを確認。現在本番投入へ進む。

## 本番投入：進行中

SUMS ID：`A76FAF721CDAB36920D53B9E8EEA9230`。XML・preflight（SL24y/sl73、ELF64、依存ライブラリ）通過。カタログ149,817行から重複LFNを除き80,982ファイル、SUMS分割は26,994 process。設定は最大5ファイル/jobで、今回SUMSは3ファイル単位へ分割している。catalogの件数は実測イベント総数ではなく、全統計の解析完了はまだ確認していない。

現在は全件のjobスクリプトを生成中で、**全件投入完了ではない**。`submit.sh`はtmux `he4_relaunch_20261022_001` 内で継続。終了時に標準スクリプトがconfiglog・joblistlog・git差分・runtime bundle・SUMS artifacts・runmetaを保存する。

補助tmux `he4_relaunch_monitor_20261022_001` は30秒間隔で`submission_status.json`を更新し、submit終了コード確定後に本記録・日別README・当日ログへ自動追記する。farmの全job終了・haddやメール送信はこの補助処理の対象外。別の再投入や旧ROOTとの統合は行わない。

- [投入ログ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/production_relaunch_20261022_001/submit.log)
- [進捗JSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/production_relaunch_20261022_001/submission_status.json)
- [少数実データ検証結果](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/production_relaunch_20261022_001/local_validation_passed.json)
- [検証ROOT](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/production_relaunch_20261022_001/local4027.root)

ローカル検証はwall time24.60秒、最大RSS578,932 KiB、終了コード0。解析本体はrequested/completed=4027、reconstructed=3858、status=0。入力は2021 run 22123040で、元カタログURLもproduction_3p85GeV_fixedTarget_2021/P24iyと照合。実入力SHAと元URLを監査先に保存した。

## Git push保留・全job終了メール

ユーザー指示によりgit commit/pushは解析job終了後まで保留。mainのHEADとindexは変更せず、fetchと別の一時worktreeの準備までで停止した。解析本体・設定・ライブラリの変更なし。自動git pushは設定していない。

今回のSUMS ID `A76FAF721CDAB36920D53B9E8EEA9230` の全job終了後、`yudai.ichikawa.d3@tohoku.ac.jp`へ1回だけメールを送る設定を追加した。投入中のため、まず標準runmetaの完成を待ち、その後既存`script/watch_job_notify.py --watch`を起動する。全26994 jobのsubmit/terminationログ一致、対象queueの消失、5分の安定待ちを確認して送信する。5分間隔で確認し、失敗・削除・未回収ROOT・完了ログ不明はメールへ明記する。hadd・再投入・即時テストメールは行わない。

tmux：`he4_completion_mail_20261022_001`。設定・待機状態：監査先の`completion_email_settings.json`、`completion_email_status.json`。起動ログ：`completion_email_launcher.log`。本監視起動後は`completion_email_watch.log`および`job/run/notify/auau3p85_anaFemtoLambda_4He_KFParticle_A76FAF721CDAB36920D53B9E8EEA9230/`へ記録。送信済み記録と排他lockで重複通知を防ぐ。

設定時のread-only確認ではsubmit済み2100、正常終了1840、実行中260。残りは投入継続中なので通知条件は未成立。sendmail実行可能であることと、tmux待機プロセスの稼働を確認した。

## 本番投入処理の終了（自動追記）

2026-10-22T15:24:19.611367+09:00。SUMS ID：`A76FAF721CDAB36920D53B9E8EEA9230`。submit.sh終了コード：0。生成process数：26958。
投入処理が終了した。各farm jobの実行完了・統合完了を意味しない。
runmeta：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_KFParticle_A76FAF721CDAB36920D53B9E8EEA9230.json`。監査：`submission_status.json`、`submit.log`、`submit.exit`。

## 重要訂正：A76の解析は全件失敗

Condorの26994 exit0とjob終了メールは解析成功ではない。全件がSUMS二列入力とpath-only parserの不一致でイベント処理前停止し、ROOTは作成されなかった。全統計解析完了とは扱わず、統合・Git pushを行わない。[原因・テンプレート修正・farm再試験](batch_input_fix_20261022.md)。旧ログ・runmetaはそのまま保持する。

修正版farm試験17DCはROOT回収・解析正常完了まで通過。全件の修正版再投入はtests/production_relaunch_20261022_002/で継続。新ID・終了メール監視も別管理。[再試験結果](batch_input_fix_20261022.md)。


## 全件終了後の検査付きhaddとメールへ変更

ユーザー指示により、SUMS `04BF2DAE47164B6F040E3A460E539289` の通知条件を「全subjob終了」から「全件解析・ROOT回収を確認し、haddと結合結果の検証が終了」に変更した。`script/watch_job_merge_notify.py` とbatch同一ROOT5の `script/singularity_root_batch.sh` を使用する。全件終了・空queue・5分のsettle後、ROOTヘッダーと解析completedログを照合し、全入力の読み取り検査、最大100ファイルずつの段階hadd、結合ROOT再読取り・ヒストグラム加算照合を行う。部分結合は許可しない。失敗時は成功とは通知せず、原因を記載した要確認メールとする。

成功時は出力先直下へ `auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root` を上書きせず公開し、結合済み個別ROOTを `subjobs/04BF2DAE47164B6F040E3A460E539289/` に保持する。その後 `yudai.ichikawa.d3@tohoku.ac.jp` へ1通通知する。旧・終了時点だけの通知workerを停止し、同名tmux `he4_completion_mail_20261022_002` を新設定で再起動した。解析投入workerは継続。現時点はrunmeta生成待ちであり、haddやメールは未実行。

検証：merge21＋archive17の38テスト通過、batch ROOT5.34/38起動、実起動予定コマンド・宛先・対象ID・部分結合無効を照合。変更前設定を `.before_hadd` として保持。監査：`tests/production_relaunch_20261022_002/hadd_notification_update.json`、監視状態：同ディレクトリ `completion_email_status.json`。hadd開始後の状態・結果・送信記録は `job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289/`。Git pushは引き続き保留。

## 本番投入処理の終了（自動追記）

2026-10-22T19:33:40.521743+09:00。SUMS ID：`04BF2DAE47164B6F040E3A460E539289`。submit.sh終了コード：0。生成process数：26969。
投入処理が終了した。各farm jobの実行完了・統合完了を意味しない。
runmeta：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289.json`。監査：`submission_status.json`、`submit.log`、`submit.exit`。
