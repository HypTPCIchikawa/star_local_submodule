# Au+Au 13p5 GeV KFParticle Λ–4He解析の出力

## KFParticle版summary（2026-10-03）

- [まとめレポート：summary_4He_KFParticle.pdf](pdf/summary_4He_KFParticle.pdf)
- [各PDFの内容・計算条件・再現手順](pdf/README.md)

下記の **2020 P24iy・542/545件の部分統合ROOT**を読み取り専用で使用した。
`tools/` の旧summaryと同じ構成で、Λ質量分布、centrality別SE/ME、sideband差引き、
purity補正、質量fit方式のCF、探索的LL fitを `pdf/` に出力している。
旧固定S/Nは流用せず、このROOTのinclusive質量分布から純度proxyを計算する。
これはペア依存の純度を確定するものではなく、fit不成立・背景モデル・欠損入力の注意点をレポートに記載した。
元のROOT、旧解析のマクロ／summary、解析cutは変更していない。

## 保存場所（2026-10-03整理）

整理完了: **2026-10-03 10:03:08 JST**。統合ROOTのサイズは241,812,808 bytes。

- 利用する統合ROOT:
  [auau13p5_anaFemtoLambda_4He_KFParticle_2020P24iy_partial_merge.root](auau13p5_anaFemtoLambda_4He_KFParticle_2020P24iy_partial_merge.root)。
  **2020年P24iy、542 / 545 subjobsの検査済み部分統合**。未統合は `149, 250, 251`。
- 元の名前 `auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053_merge.root`
  から内容を変えず改名した。SUMS IDは `1213058CE62BCF113B2E3C75ACBD7053`。
- 統合済みの個別ROOT: [subjobs/1213058CE62BCF113B2E3C75ACBD7053/](subjobs/1213058CE62BCF113B2E3C75ACBD7053/)。
  2026-10-01 23:29（日本時間）に542本を移動済み。削除・内容変更なし。
- 移動先のREADMEと `archive_manifest.json` が、旧パスと現在の保存先の対応表。
- hadd途中の7 ROOT: [archive/hadd_work/1213058CE62BCF113B2E3C75ACBD7053/](archive/hadd_work/1213058CE62BCF113B2E3C75ACBD7053/)。
  初回失敗時の空作業ディレクトリも同じ階層の `_failed_path_20261001_1926/` に保持。
- 試験結果: [`tests/pilot_20260930/`](tests/pilot_20260930/) に分離。本投入と混ぜてmergeしない。
- 試験の詳細: [`tests/pilot_20260930/README.md`](tests/pilot_20260930/README.md)。
- [今回の整理完了・SHA256記録](../../job/run/merge_notify/auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053/layout_20261003/result.json)
- [今回の改名・作業ディレクトリ移動の対応表](../../job/run/merge_notify/auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053/layout_20261003/relocation_manifest.json)

過去のresult・監査TSV・投入記録は書き換えていないため、旧名・旧 `.merge_staging_...` パスが残っている。
現在の場所は上記の対応表を参照する。解析の再実行・再hadd・ROOTデータの削除は行っていない。
統合ROOTと `subjobs/`・`archive/hadd_work/` を同時にhaddしないこと（二重加算になる）。
13.5 GeVと3.85 GeVは、それぞれ別の結果として保持している。

既存の試験ROOTはファイル名を変えず、旧 `rootfile/auau13p5_anaFemtoLambda_4He_KFParticle_highpurity_pilot_20260930/`
から上記サブディレクトリへ移動済み（2026-10-01）。ROOTの移動前後のSHA-256は一致。
投入時のXML/config/runmetaの保存記録は改変していない。
今後も試験・cut比較は `tests/<条件または日付>/` のように分離し、条件をREADMEに残す。

## 本投入の解析条件

| 項目 | 条件 |
|---|---|
| 対象 | Au+Au 13p5 GeV fixed-target、2020、P24iy |
| ファイル選択 | `st_physics`、non-ADC、`daq_reco_picoDst`、`storage=local` |
| 規模 | FileCatalogに該当する全ファイル、各ファイルの全イベント |
| 分割 | 最大5ファイル/job |
| 実行環境 | SL24y、ROOT 5.34/38、STAR Singularity |
| 解析 | `run_anaFemtoLambda_4He.C` → `anaFemtoLambda_4He.C` |
| Λ選択 | KFParticle Finder/Topo、高purity設定。成功した試験の物理設定を維持 |
| Λ質量ゲート | 平均1.11596、σ=0.00296 GeV/c²、±3σ |
| 原子核 | 4He、既存NuclearId互換選択 |

設定入口は `config/mainconf/main_auau13p5_anaFemtoLambda_4He_KFParticle.yaml`。
Λのcutは `config/cuts/kf/kf_auau13p5_anaLambda_KFParticle_highpurity.yaml`。
`highpurity`は出力ディレクトリ名から省いているが、cutを緩めたわけではない。
旧Helix解析や他の核種・衝突エネルギーは今回の本投入対象ではない。

入力はXML内のFileCatalog検索で決める。少数イベントのvalidationリストは使用しない。
実際に選ばれたファイルはSUMSのdataset/listおよびrunmetaの保存物で確認する。
FileCatalogのイベント数は実ROOTのエントリ数と一致しない場合があるため、
処理件数の確定には解析ログの `requested/completed` を用いる。

2026-10-01の事前catalog確認では、non-ADC **1,635ファイル・22 runs・約1.85 TB**。
カタログ記載イベント数の合計は43,242,488だが、実際の処理数ではない。
最大5ファイル/jobの指定に対し、SUMSが今回は3ファイルずつに分割し、実際は**545 jobs**となった。
ADCの186ファイルは今回の対象外。成功したfarm試験と従来の本投入と同じstream条件。

出力にはΛ質量、アクセプタンス、SE/MEのk*、質量窓の低側・高側、sideband、
centrality別k*–mTを含む。mTの定義とヒスト名は
`mdfiles/femto/implementation/femto_lambda_histograms_20260930.md` を参照。

## 投入・再現性の記録

2026-10-01（日本時間）、`starsub03`から本投入済み。
`submit.sh`は終了コード0、SUMSは `Scheduling successful`。

| 項目 | 記録 |
|---|---|
| SUMS jobid | `1213058CE62BCF113B2E3C75ACBD7053` |
| subjob番号 | `0`～`544` |
| 登録数 | 545 jobs（各3ファイル） |
| Condor cluster | `320168`、`320169`、`320170`、`320171`、`320172`、`320173` |
| cluster別job数 | 順に100、100、100、100、100、45 |
| 実入力 | 全1,635ファイル。subjobリスト合計1,635件、ユニーク1,635件、重複なし |
| 受付確認 | 2026-10-01 13:49頃（日本時間）に実行中303、待機242、合計545 |

受付時の検索結果は1,686行だったが、SUMSが同一LFNの51行を除いて1,635ファイルとした。
生成された545個のsubjobリストについて、各3件、番号の欠落なし、全入力の重複なしを確認した。
上記キュー状態は受付直後の一時点であり、今後変化する。

本投入ROOT名は
`auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053_<subjob>.root`。
投入時はこのディレクトリ直下へ回収した。現在は統合検査済みの542本を
`subjobs/1213058CE62BCF113B2E3C75ACBD7053/` に整理済み。
`tests/`内のROOTは本投入のmerge対象外。

- [本投入mainconf](../../config/mainconf/main_auau13p5_anaFemtoLambda_4He_KFParticle.yaml)
- [本投入XML](../../job/joblist/joblist_auau13p5_anaFemtoLambda_4He_KFParticle.xml)
- [投入時設定](../../job/run/configlog/config_auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053.txt)
- [投入時XML](../../job/run/joblistlog/joblist_auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053.xml)
- [runmeta](../../job/run/runmeta/runmeta_auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053.json): コード・ライブラリ・生成物の保存先。

投入コマンド（実行済み。確認なしに再実行しないこと）:

```bash
cd job/run
./submit.sh ../joblist/joblist_auau13p5_anaFemtoLambda_4He_KFParticle.xml
```

今回は本投入の受付までを確認し、全job完了待ち・再投入・自動mergeは行っていない。
試験と同じ8項目の物理設定、ROOT5/Singularityビルド、XML/csh構文、入口テスト、
ヒストグラム契約テストを投入前に確認した。物理カットや解析ソースは変更していない。

stdout/stderrは `log/auau13p5_anaFemtoLambda_4He_KFParticle/` と
`err/auau13p5_anaFemtoLambda_4He_KFParticle/` に保存する。
コード・設定・ライブラリの記録は `job/run/runmeta/`、設定一覧は `job/run/configlog/`、
投入XMLは `job/run/joblistlog/`。同じSUMS jobidで照合する。
投入完了は解析完了ではない。全jobの正常終了・ROOT回収は別途確認する。

生成XMLには、試験で成功したROOTSYS明示、scratch作成・書込み検査、
SUMS入力の第1列正規化、1～5ファイルの件数検査、入力転送失敗時の停止、
ROOT出力の存在・非空検査を追加した。
XMLを再生成するとこれらの追加処理が失われるので、その場合は再検証する。
各jobは開始時にソース・設定・ライブラリをscratchへコピーする。
解析条件を揃えるため、本投入のjobが終わるまでは該当ファイルを変更しないこと。

## 終了確認とメール通知（2026-10-01）

2026-10-01 15:25（日本時間）の確認では、全545 jobのCondor終了記録が揃い、
キューは0件、終了コードは全件0だった。ただし、ROOTの回収・ヘッダーと解析完了ログを
確認できたのは542件で、subjob `149`・`250`・`251` の3件は両方とも未確認。
したがって「解析全件成功」ではなく、未回収分の確認が必要な状態である。
ROOT内部の全ヒストグラム検査・物理QA・merge・再投入は実施していない。

今回のSUMS jobidについて、`yudai.ichikawa.d3@tohoku.ac.jp` 宛てに
「終了・要確認」の集約通知を1通送信し、**2026-10-01 15:36:55（日本時間）**に
ローカルsendmailの受付成功を確認した。相手の受信箱への到着は未確認。
既に全job終了済みだったため、`script/watch_job_notify.py --watch --settle-sec 0` で
終了状態・ROOT・ログを再確認して通知した。処理は終了し、常駐監視は残していない。

- [通知・出力確認の記録](../../job/run/notify/auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053/sent.json)
- [通知設定](../../job/run/notify/auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053/settings.json)
- [通知スクリプトの使い方](../../job/run/README.md#完了メールのみmerge再投入なし)

送信済み記録により、同じjobidへの自動重複通知を防止する。
今後の新しい投入には、その投入のrunmetaとイベントログを指定して別途起動する。

## 回収済みROOTのhadd（2026-10-01）

ユーザーの依頼により、回収済み542件を1ファイルにまとめ、
2026-10-01 19:39:56（日本時間）にROOT再読取と282ヒストの加算検査まで成功した。
統合ROOTは241,812,808 bytes。
未回収のsubjob `149`・`250`・`251` は現在も欠落しているため、
**全545件の統合ではなく542件の部分統合**である。これら3件の再投入はしていない。
入力ROOTの合計サイズは982,160,691 bytes。元のファイルは保持する。
23:29:51に542本を `subjobs/1213058CE62BCF113B2E3C75ACBD7053/` へ移動し、
全ファイルのsize/inode/mtimeが移動前と一致することを検証した。
親ディレクトリ直下のROOTは統合ファイル1本のみ。

未回収3件はこの移動でも補完されておらず、統合は引き続き542/545件の部分統合である。
今後も検査付きhadd成功後に同じ形式で整理し、再検証可能な元ROOTは削除しない。

統合出力名は `auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053_merge.root`。
ROOTの読み取り検査、段階hadd、結合後の再読取り・ヒストグラム加算検査を通過してから公開する。
既存mergeは上書きしない。入力の正確な一覧・欠落番号・処理ログを保存する。

- ホスト: `starsub03`、tmux socket: `star-lambda-production`。
- session: `merge13p5_1213058`。2026-10-01 19:31頃（日本時間）再開。
- [処理ログ](../../job/run/merge_notify/auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053/monitor.log)
- [現在の状態](../../job/run/merge_notify/auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053/status.json)
- [最終結果（処理終了後に生成）](../../job/run/merge_notify/auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053/result.json)
- [入力一覧](../../job/run/merge_notify/auau13p5_anaFemtoLambda_4He_KFParticle_1213058CE62BCF113B2E3C75ACBD7053/inputs.list)

19:26の初回処理は、Singularity内で `/star/u/...` の別名パスが見えず、
最初の入力を開く段階で停止した。ROOT破損ではなく、hadd自体も未実行だった。
読み書き先を `/gpfs/...` 実パスへ解決するよう修正し、失敗時の状態記録・空の作業領域を
`_failed_path_20261001_1926` 付きの別名へ移して保持した上で再開した。
このhaddでは追加メールを送らない（13p5のjob終了通知は15:36に送信済み）。
3p85のhadd完了通知は別途その解析の監視処理が送信する。
