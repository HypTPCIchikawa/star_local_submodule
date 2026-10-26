# Au+Au 3.85 GeV KFParticle Λ–4He出力

## 2026-10-17：核PID purity統合と本番cutoff

- [低R＋TOF後高Rの統合PDF](pdf/nuclear_pid_cutoff_20261015/views_tpc_tof_final_20261017_001/purity_tpc_tof_combined.pdf)：新概要3頁＋元低R3頁。採用域と除外・未評価域を表示。
- [本番採用するR区間](pdf/nuclear_pid_cutoff_20261015/production_ranges_20261017_001/README.md)：4核種、低Rの穴と高R TOF条件。
- [低R全区間のfit集](pdf/nuclear_pid_cutoff_20261015/fit_atlas_native_20261017_003/README.md)：32冊2560頁、未支持fitも掲載。
- [当日のTodo・結果・新旧比較・検証](../../Todo/2026/10/17/README.md)。

本番入口を共通校正と検証済みcutoffへ接続したが、この既存統合ROOT自体を再解析・上書きしたものではない。
元ROOT、2019年の隔離データ、旧study/PDFは保持。今回job投入・hadd・メールは行っていない。

## KFParticle版summary（2026-10-03）

- [まとめレポート：summary_4He_KFParticle_3p85.pdf](pdf/summary_4He_KFParticle_3p85.pdf)
- [各PDFの内容・計算条件・再現手順](pdf/README.md)

下記の **2021 P24iy・26,956/26,994件の部分統合ROOTだけ**を読み取り専用で使用した。
2019年、pilot、subjob個別ROOTは混ぜていない。
旧summaryと同じ構成で、Λ質量分布、centrality別SE/ME、sideband差引き、purity補正、
質量fit方式のCF、探索的LL fitを `pdf/` に保存している。
旧固定S/Nを流用せず、このROOTのinclusive質量分布から純度proxyを計算した。
fitが不成立・不安定な場合は図に警告し、purityとモデルの仮定・欠損38件の扱いもレポートに明記した。
元のROOT、旧解析のマクロ／summary、解析cutは変更していない。

## 利用する統合ROOT（2026-10-03整理）

**2021年P24iy、26,956 / 26,994 subjobsの部分統合**。
欠損38件は含まない。「完成」はこの回収済み部分の統合・整理を意味し、全入力の解析完了ではない。
整理完了: **2026-10-03 10:05:24 JST**。統合ROOTのサイズは10,522,646,853 bytes。

- 利用ファイル:
  [auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root](auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root)
- SUMS ID: `0BEFC64E62F19BBA534A129CABC9B106`。
- 元の個別ROOT26,956本: [subjobs/0BEFC64E62F19BBA534A129CABC9B106/](subjobs/0BEFC64E62F19BBA534A129CABC9B106/)。
- hadd途中の274 ROOT: [archive/hadd_work/0BEFC64E62F19BBA534A129CABC9B106/](archive/hadd_work/0BEFC64E62F19BBA534A129CABC9B106/)。
- 試験結果は `tests/`、2019年の旧結果は `archive/2019P24iy/`、中断した入力復旧作業は
  `recoveries/20261003/` に保持。相互に混合しない。

全282ヒストグラムの名前・型・bin数・Entriesは一致し、ROOT読取り検査も完了。
`hRefMultWeight` と `hVz` の集約量の相対差のみ、ユーザーが **0.000015%以下を今回に限り許容**した。
最大差は約0.000014654%。一般の検査基準・解析cutは変更していない。
丸め誤差だけが原因だと証明したものでも、物理QAの完了認定でもない。

最新の完了判定と対応表:

- [整理の完了記録](../../job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_0BEFC64E62F19BBA534A129CABC9B106/manual_acceptance_20261003/result.json)
- [明示的な許容・SHA256の記録](../../job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_0BEFC64E62F19BBA534A129CABC9B106/manual_acceptance_20261003/acceptance.json)
- [旧→新パスの対応表](../../job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_0BEFC64E62F19BBA534A129CABC9B106/manual_acceptance_20261003/relocation_manifest.json)

以前の自動処理の `result.json` / `status.json` は **failedの履歴をそのまま保存**している。
今回の承認後の状態は、上記 `manual_acceptance_20261003/result.json` を参照する。
古い監査ログにある `.merge_staging_...` は移動前のパスで、対応表の `stagingArchive` が現在の場所。
旧予定名 `..._0BEFC64E62F19BBA534A129CABC9B106_merge.root` ではなく、上記の年・production・partial付きの名前を使用する。

**ROOTデータは削除していない。統合ROOTと `subjobs/`・`archive/hadd_work/` を一緒にhaddしないこと**
（同じイベントの二重加算になる）。欠損jobの再投入は引き続き停止中。

以下は整理前の作業履歴。現在の保存先・状態はこの節と完了記録を優先する。

## 2026-10-03：回収済み分の部分haddと欠損分の復旧

### 終了時の方針変更：再投入は後日

ユーザー指示により、欠損分の復旧用tmux `recovery3p85_0BEFC64` とその処理を停止した。
**再投入はまだ行っていない。** `submission.started`、投入ログ、再投入jobidは未作成。
`job/run/recovery_3p85_0BEFC64_20261003/` にカタログ照会・入力検査の記録と
復旧用スクリプト／設定を保持した。取得した入力も
`recoveries/20261003/input_cache/` に保持し、共有の元PicoDstは変更していない。
`.attempt_*` は途中取得・未検証・不良の可能性があり、そのまま解析へ使用しない。
次回は記録を確認して明示的に再開する。開始済みガードを削除してlauncherを再実行しない。

**部分hadd用のtmux `partial3p85_0BEFC64` は停止していない。**
終了確認時点は26,956本を対象とする統合前のROOT内容検査中であり、その後の
段階hadd・加算検査・個別ROOT整理・メール通知へ進む。SSH／VS Code切断後も継続する。

対象SUMS IDは `0BEFC64E62F19BBA534A129CABC9B106`（2021 P24iy）。
全26,994jobは10月2日08:20:36 JSTに終了したが、回収ROOTは26,956本で38本不足している。
SUMSは全件dispatch後、Java heap不足で失敗終了し、自動hadd監視へ引き継がれなかった。

ユーザー承認により、10月3日00:11:26 JSTから **`--allow-partial` を明示した検査付き部分hadd** を
`starsub03`のtmux socket `star-lambda-production`、session `partial3p85_0BEFC64` で開始した。
開始時点はROOTヘッダー・解析完了ログの全件照合中。照合とROOT内部読取り検査が済んでから
100本ずつの段階hadd、統合後の再検査、個別ROOT整理、メール通知へ進む。
起動したことと、haddが完了したことは区別する。SSH／VS Codeを閉じても継続する。

- [監視設定](../../job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_0BEFC64E62F19BBA534A129CABC9B106/settings.json)
- [現在の段階](../../job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_0BEFC64E62F19BBA534A129CABC9B106/status.json)
- [処理ログ](../../job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_0BEFC64E62F19BBA534A129CABC9B106/monitor.log)
- 成功時の公開名: `auau3p85_anaFemtoLambda_4He_KFParticle_0BEFC64E62F19BBA534A129CABC9B106_merge.root`。
  **この名前にpartialが含まれなくても、全件統合ではない。** 正確な統合・除外番号は処理後の
  `result.json` と `input_manifest.json` を根拠にする。
- 検査後の元ROOT保存先: `subjobs/0BEFC64E62F19BBA534A129CABC9B106/`。削除しない。
- 通知先: `yudai.ichikawa.d3@tohoku.ac.jp`。部分統合であること・除外番号を明記する。

不足38件:

```text
362, 1967, 2001, 5273, 8944, 10263, 10340, 10470, 10784, 10786,
11005, 11006, 11494, 12357, 21109, 21385, 21400, 21888, 22032, 22258,
22762, 23201, 23362, 23486, 23536, 23539, 23587, 23590, 23598, 23628,
23672, 23673, 23678, 23924, 24043, 24063, 25445, 25555
```

ログ上は、入力ROOT読取不能34件、PicoDst tree欠落2件、空PicoDst1件、途中のbasket読取エラー1件。
元ファイル・replicaとstaging後の問題を切り分け、正常な同一2021 P24iy入力を確保できた分を
別jobid・別出力先で再投入する調査を開始した。共有の元PicoDstは書き換えない。
再投入ROOTは元jobidのディレクトリへ差し戻さず、部分統合との二重加算を避ける。
以下の本投入開始時の記録は履歴として保持する。

## 現在の方針：2021年P24iyのみ（2026-10-01更新）

今後の解析には `production_3p85GeV_fixedTarget_2021` / P24iyのみを使用する。
本投入analysis_infoとXML、4核種のKF高purity設定を2021年へ変更した。
物理cut・centrality・ヒストグラム設定は変更していない。
2021年3ファイル限定試験は54,914イベントを処理し、回収・検査付きhadd・
イベント台帳／ヒスト検査が成功した。QA PDFも生成済み。purityの定量評価や
全データの網羅確認は未実施。本投入処理は2026-10-01 23:51:57（日本時間）に開始した。
試験の最新状態は
[tests/2021P24iy_pilot_20261001/README.md](tests/2021P24iy_pilot_20261001/README.md)を参照。

この下に記録されたSUMS ID `A8A883FA3E75C40C723A9456C8FCE429` のROOTは
**2019年の旧結果**で、2021年の解析結果ではない。2026-10-01にユーザー指示で
[archive/2019P24iy/A8A883FA3E75C40C723A9456C8FCE429/](archive/2019P24iy/A8A883FA3E75C40C723A9456C8FCE429/)
へ72本全てを移動した。サイズ・inode・mtimeの一致と親直下にROOTが残っていないことを検証済み。
削除はしていない。今後の2021年結果への混入やhaddを禁止する。ROOT選択はディレクトリ全体の
ワイルドカードではなく、runmetaとjobidごとの入力一覧で行う。

2021年のカタログ照会（P24iy、non-ADC `st_physics`、`storage=local`）では、
filename重複除去後80,982ファイル、約72.96 TB、カタログ記載2,184,032,750 events。
これは実PicoDstのEntriesを数えた値ではなく、希望する23 billion Entriesとの
対応・全データの収録範囲・各ファイルの読取可否はまだ確認できていない。

## 2021 P24iy本投入（2026-10-01）

ユーザーの本投入指示を受け、成功した2021年pilotと同じ物理設定で開始した。
対象は2021 P24iy / non-ADC `st_physics` / catalog `storage=local` に一致する
全ファイル・全イベント（`nFiles: all`、`maxEvents: -1`、最大5ファイル/job）。
2019年データを追加・再投入しない。23B実Entriesを網羅したとの確認ではない。

- ホスト: `starsub03.sdcc.bnl.gov`
- tmux: socket `star-lambda-production`、session `full3p85_2021_1001`
- [一回限定launcher](../../job/run/launch_production_3p85_2021P24iy_20261001.sh)
- [投入ログ](../../job/run/production_3p85_2021P24iy_20261001_submit.log)
- [投入から監視への引継ぎログ](../../job/run/production_3p85_2021P24iy_20261001_workflow.log)
- 再投入防止記録: `job/run/production_3p85_2021P24iy_20261001.started/`。削除・再実行禁止。

SUMS ID: `0BEFC64E62F19BBA534A129CABC9B106`。
SUMSはカタログ149,817行から重複LFN 68,835行を除外し、**80,982ファイル**を選択した。
分割予定は**26,994 subjobs**（各3ファイル、subjob `0`～`26993`）。
2026-10-01 23:59:54 JST時点は600件受付、209件実行中、391件待機、holdなし、
終了記録0件。生成済み637入力リストの1,911パスは全て2021 P24iyで範囲不一致なし。
ROOT回収も1本開始したが、全解析の成功やROOT検査完了の意味ではない。
大量のjobを末尾番号から順次生成・投入中であり、全件受付完了ではない。

全件の正常受付後に新しいrunmetaと実入力リスト（年・production・重複・分割数）を検査し、
そのSUMS IDだけの終了監視へ自動で引き継ぐ。監視記録は
`job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_<jobid>/` に作成する。
現在は投入段階なので、hadd監視の `settings.json` / `status.json` はまだ作成されていない。
runmetaも `submit.sh` が全件の投入を終えた後に保存する。

終了後は全出力・解析完了ログの確認 → ROOT読取り検査 → 段階hadd → 加算検査 →
統合ROOT公開 → 個別ROOTを `subjobs/<jobid>/` に移動 →
`yudai.ichikawa.d3@tohoku.ac.jp` への通知、の順に処理する。
統合ROOTはこのディレクトリ直下、試験ROOTは `tests/`、2019年の旧結果は `archive/2019P24iy/`。
削除・自動再投入・無断の部分統合はしない。失敗時は警告通知する。
SSH切断後も継続するが、サーバー再起動には耐えない。
投入または投入後の入力監査で停止した場合はログを確認し、jobを二重投入しないこと。

投入前にfarmと同じSingularity/SL24y/ROOT5で増分ビルド成功、
XML構文・出力ディレクトリ書込み・hadd監視21件・個別ROOT整理17件のテスト成功。
実行中はjobが開始時にコピーする解析ソース・設定・ライブラリを編集／再ビルドしない。

## 以下は2019年投入の履歴（今後の解析対象外）

## 解析条件

2026-10-01（日本時間）に、既存の3p85高purity条件で本投入。
3.85 GeVはfixed-targetのビームエネルギー表示で、colliderの√sNNではない。

| 項目 | 条件 |
|---|---|
| データ | 2019 FXT、`production_3p85GeV_fixedTarget_2019`、P24iy |
| stream | `st_physics`、non-ADC、`daq_reco_picoDst`、catalog `storage=local` |
| 入力範囲 | 上記catalog検索の全ファイル・全イベント。小統計validation listは使わない |
| 分割指定 | 最大5ファイル/job |
| 解析 | `run_anaFemtoLambda_4He.C` → `anaFemtoLambda_4He.C` |
| 環境 | SL24y、ROOT 5.34/38、farmと同じSingularity |
| Λ選択 | `kf_auau3p85_anaLambda_KFParticle_highpurity.yaml` |
| 原子核 | 4He、既存NuclearId互換選択 |

mainconfのanalysis項目以外の参照は、検証済み
`main_auau3p85_anaFemtoLambda_4He_KFParticle_highpurity.yaml` と全て同一。
event、centrality、KF cut、nuclear PID、maker、mixing、hist設定を変更していない。
過去の別年・別productionの3p85全データをまとめて投入するものではない。

事前catalog確認は225 uniqueファイル、約137 GB。
catalog記載の5,744,126 eventsは実PicoDst entriesの確定値ではない。
処理件数は解析完了ログの `requested/completed` で確認する。

## 保存先と自動処理

- 2019年の個別ROOT: `archive/2019P24iy/A8A883FA3E75C40C723A9456C8FCE429/` 以下。
  ファイル名は `auau3p85_anaFemtoLambda_4He_KFParticle_<jobid>_<subjob>.root` のまま。
- 統合ROOT: `auau3p85_anaFemtoLambda_4He_KFParticle_<jobid>_merge.root`。
- 今後の試験結果は `tests/<条件または日付>/` に分離する。
- 全job終了後、全出力と解析完了ログを確認してからhaddし、ROOT再読取りと
  ヒストグラム加算検査を通過した時点で統合ファイルを公開する。
- その後 `yudai.ichikawa.d3@tohoku.ac.jp` へ1通通知する。
  出力欠落・job失敗・hadd失敗があれば成功扱いせず警告し、自動再投入や部分統合はしない。
- 元のROOTは保持し、既存mergeは上書きしない。物理QAは別途実施する。

[本投入mainconf](../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml)、
[本投入XML](../../job/joblist/joblist_auau3p85_anaFemtoLambda_4He_KFParticle.xml)、
[監視・結合手順](../../job/run/README.md)を参照。

## 再開時の事前検証

Singularityで `make --no-clean -j4 all` 相当の増分ビルド成功。
入口テスト6件、ヒストグラム契約テスト、監視通知テスト、
ROOT5のTH1D・TH2D・TTree・TNamedを含む人工3ファイルの2段階haddと加算検査が成功。
解析ソース・物理cutは今回変更していない。

## 投入と終了後処理の記録

2026-10-01 19:26頃（日本時間）に `starsub03` から投入し、SUMSの
`Scheduling successful` と投入時XML・config・runmeta保存を確認した。

| 項目 | 記録 |
|---|---|
| SUMS jobid | `A8A883FA3E75C40C723A9456C8FCE429` |
| Condor cluster | `330143` |
| job数 | 75（subjob `0`～`74`） |
| 実入力 | 225 uniqueファイル、各jobに3ファイル。重複・番号欠落なし |
| 19:30:54の状態 | 67件終了コード0、8件実行中。途中時点の状態であり最終結果ではない |
| tmux socket | `star-lambda-production`（ホスト `starsub03`） |
| 監視session | `merge3p85_A8A883F`（19:31頃起動） |
| 通知先 | `yudai.ichikawa.d3@tohoku.ac.jp` |

- [runmeta](../../job/run/runmeta/runmeta_auau3p85_anaFemtoLambda_4He_KFParticle_A8A883FA3E75C40C723A9456C8FCE429.json)
- [投入ログ](../../job/run/production_3p85_20261001_submit.log)
- [終了監視・haddログ](../../job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_A8A883FA3E75C40C723A9456C8FCE429/monitor.log)
- [現在の状態](../../job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_A8A883FA3E75C40C723A9456C8FCE429/status.json)
- [最終結果（処理終了後に生成）](../../job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_A8A883FA3E75C40C723A9456C8FCE429/result.json)

統合出力名は `auau3p85_anaFemtoLambda_4He_KFParticle_A8A883FA3E75C40C723A9456C8FCE429_merge.root`。
5分ごとに終了を確認し、全件終了確認からさらに5分待機した後、検査・hadd・再検査を行う。
成功・失敗いずれも1通通知し、通知後プロセスは終了する。
SSH/VS Codeを閉じても継続するが、サーバー再起動の場合は記録を確認して復旧が必要。

```bash
tmux -L star-lambda-production list-sessions
```

投入初回のtmux起動は既定tcshにbash用redirectを渡したため、投入前に終了した。
bash明示後に上記1回だけ投入されており、二重投入はしていない。

## 2026-10-13: ³He/⁴He PID校正と³He背景CF用出力

現行の3p85⁴He mainconf（通常・highpurity・tofExtended）は新しい校正と³He controlを参照する。
上記の過去の投入記録は当時の条件であり、既存の統合ROOTには新校正・³He controlは含まれない。
昔のconfig/runmeta・入力・統合ROOTは保持し、新設定の小試験は
[tests/helium_calibrated_he3Control_20261013/](tests/helium_calibrated_he3Control_20261013/README.md)
へ分離した。以前のROOTと同じglobで再haddしない。

fit結果を使う[校正YAML](../../config/calibration/dedx_helium_auau3p85_2021P24iy_from040.yaml)、
[新しいTPC選択窓のPurity比較PDF](pdf/he4_purity_20261013/calibrated_selection/purity_calibrated_windows.pdf)、
[Task3実装・検証報告](../../Todo/20261013/task3_helium_pid_cf_implementation_report.md)を参照。
各現行mainconfの`_pidLegacy`版で従来条件へ戻せる。今回本投入はしていない。
