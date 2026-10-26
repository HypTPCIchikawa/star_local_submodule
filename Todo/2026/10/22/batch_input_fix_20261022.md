# 2026-10-22 本番全job未解析：入力リスト形式不一致の修正

ユーザーから26994 exit0/ROOT0件の連絡。全26994 stderrでERROR: malformed input-list entryを確認。SUMSのpath+catalog events二列を、厳密なpath-only parserへ渡し、イベントループ前に停止。保存前joblistには整形と出力guardがあったが、generate_joblistが使用する共有テンプレートにはなく、再生成後のXMLで失われた。ローカル1列4027イベント試験は通るがfarm経路を検証していなかった。

SUMS wrapperが子解析の失敗を外側のexit0へ隠したため、Condor終了は解析成功を意味しない。今回は解析0イベント、ROOTは回収失敗以前に作られていない。旧runmeta/logは保持してfailureと明記。Git pushは引き続き保留。

修正中：テンプレートで二列を正規化し、入力/SCRATCH/コピー/実行/出力を検査する。通知文もbatch全終了と解析正常完了を区別する。同形式のlocal試験と少数farmでROOT回収まで確認し、全件再投入はその結果に基づく。監査：rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/batch_input_fix_20261022_001/。

## 原因の確定と修正

全26994 stderrに同一の入力形式エラー。入力数のcatalog列が付いたままのリストを渡し、Prepare段階でrejectしており、入力イベント解析は0。ローカル4027イベント試験はpath-onlyリストであったため差異を検出できなかった。本番XMLを再生成した際、事前の手修正版の整形・guardが共有テンプレートに反映されていなかったことを見落とした。

共有job_template_from_conf.xmlを修正し、SUMS二列リストをpath-onlyへ正規化、数字でない件数・余分な列・未staging/空リスト・過剰ファイル数を拒否する。SCRATCHとコピー・実行status・期待ROOT存在も検査。固定root4star/toolchainを明示し、scratch内にのみコピーする。generated XMLにも処理が残ることを検査した。

通知のsubjectは元から「終了・要確認」だったが、本文先頭が「解析jobの全件終了」だったため誤解を招いた。失敗時の先頭を「バッチjobは全件終了しましたが、解析の正常完了は確認できません」に変更し、真の成功時だけ「解析の正常完了」と表示する。既送信メールの記録は変更せず、追加メールは現時点では送っていない。

workerの入力処理を実際のcshで試す7テストと通知10テストを通過。新farm試験は2021 P24iy catalog先頭3ファイル、1 job、2000イベント、physics/mainの係数・cutは同じ。

farm試験SUMS ID：`17DC1E7C434E4D58D6E3E234832FFD4D`。mainconf：`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_batchcheck20261022.yaml`。解析名をbatchcheck付きにして旧結果と分離。Condor cluster：543888。ROOT回収・completed・イベント数・追加TH2を確認するまで全件再投入を行わない。

## Farm試験：ROOT回収・解析・ヒスト確認まで通過

17DC pilotは2000入力・1908再構成・選択Λ30、最終⁴He121（高R TOF m²後5）、signal Λアクセプタンス23。ROOT3548321 byteを指定出力先へ回収し、completedとledger、全1787キー/型/軸、追加β・R/TOF区分の整合性、SE/ME信号上下分割を閉ROOTから確認。2000イベントのsignal ⁠⁴He SEは0、ME85であり、少数試料で全histの全binが非ゼロという意味ではない。

この実farm経路の確認後、修正版の全件XMLを通常mainconfから再生成・guardを照合し、tmux he4_relaunch_20261022_002で全件再投入を開始した。新しい監査先：tests/production_relaunch_20261022_002/。新IDと進捗はsubmission_status.jsonへ保存する。終了メールも新IDだけを対象に設定し、解析の正常完了とbatch終了を区別する。全件投入完了や本番解析完了とはまだ扱わない。旧A76のログ/出力記録は保持、Git pushは保留。

修正版の全件再投入SUMS ID：`04BF2DAE47164B6F040E3A460E539289`。投入進行中。終了通知もこの新IDを対象に設定。旧A76は失敗記録として保持。


## 全件終了後の検査付きhaddとメールへ変更

ユーザー指示により、SUMS `04BF2DAE47164B6F040E3A460E539289` の通知条件を「全subjob終了」から「全件解析・ROOT回収を確認し、haddと結合結果の検証が終了」に変更した。`script/watch_job_merge_notify.py` とbatch同一ROOT5の `script/singularity_root_batch.sh` を使用する。全件終了・空queue・5分のsettle後、ROOTヘッダーと解析completedログを照合し、全入力の読み取り検査、最大100ファイルずつの段階hadd、結合ROOT再読取り・ヒストグラム加算照合を行う。部分結合は許可しない。失敗時は成功とは通知せず、原因を記載した要確認メールとする。

成功時は出力先直下へ `auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root` を上書きせず公開し、結合済み個別ROOTを `subjobs/04BF2DAE47164B6F040E3A460E539289/` に保持する。その後 `yudai.ichikawa.d3@tohoku.ac.jp` へ1通通知する。旧・終了時点だけの通知workerを停止し、同名tmux `he4_completion_mail_20261022_002` を新設定で再起動した。解析投入workerは継続。現時点はrunmeta生成待ちであり、haddやメールは未実行。

検証：merge21＋archive17の38テスト通過、batch ROOT5.34/38起動、実起動予定コマンド・宛先・対象ID・部分結合無効を照合。変更前設定を `.before_hadd` として保持。監査：`tests/production_relaunch_20261022_002/hadd_notification_update.json`、監視状態：同ディレクトリ `completion_email_status.json`。hadd開始後の状態・結果・送信記録は `job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289/`。Git pushは引き続き保留。


2026-10-23再確認：修正版04BFは全batch終了・queue0。26955件が全予定entry正常処理（約21.01億event）、39件が入力ROOT関連で未成功。strict全件条件でhadd未開始。[確認結果](../23/he4_completion_audit_20261023.md)。
