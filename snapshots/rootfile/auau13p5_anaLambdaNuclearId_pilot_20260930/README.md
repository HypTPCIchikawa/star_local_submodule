# 13p5 GeV 旧LambdaNuclearId farm試験（2026-09-30）

## 目的・条件

旧Helix解析 `anaLambdaNuclearId` の入力転送・実行・出力回収を試験する。
本投入ではなく、完了と物理QAの確認は後日行う。

| 項目 | 内容 |
|---|---|
| 入力 | Au+Au 13p5 GeV fixed-target、2020、P24iy、run 21033026 |
| ファイル | `st_physics_21033026_raw_7500026.picoDst.root`（non-ADC） |
| 制限 | 1ファイル、1 job、その全3,202入力イベント（maxEvents: -1） |
| mainconf | `config/mainconf/main_auau13p5_anaLambdaNuclearId_pilot_20260930.yaml` |
| XML | `job/joblist/joblist_auau13p5_anaLambdaNuclearId_pilot_20260930.xml` |
| 実行環境 | SL24y、ROOT 5.34/38、STAR Singularity |

物理カット・ヒストは元の `main_auau13p5_anaLambdaNuclearId.yaml` のまま。
変更は投入識別子、入力固定、イベント上限、出力先などの試験メタデータのみ。
入力記録は `config/picoDstList/auau13p5_farm_pilot_20260930.list`。
XMLではファイル名をFileCatalogの検索条件で固定し、最大1ファイルとする。
同じ入力で新KFParticle版4He解析も別jobで実行するが、旧解析は複数核種なので
ROOT内の全ヒストが一致する試験ではない。

## 出力の見方

`auau13p5_anaLambdaNuclearId_pilot_20260930_<SUMS jobid>_<subjob>.root`
がこの試験の解析ヒストグラム。既存本番出力は上書きしない。

- stdout: `log/auau13p5_anaLambdaNuclearId_pilot_20260930/stdout.<JOBID>.out`
- stderr: `err/auau13p5_anaLambdaNuclearId_pilot_20260930/stderr.<JOBID>.err`
- 設定/XML/コードの記録: `job/run/configlog/`、`job/run/joblistlog/`、`job/run/runmeta/` の同じanaName・jobidのファイル。

投入前に本ディレクトリとstdout/stderrの出力先を作成する。
実行ノードの `$SCRATCH` はSUMSが割り当てるため、必要な解析サブディレクトリを
job開始時に作成し、存在・書込み可能性を検査してから解析する。

## 投入記録

2026-09-30 11:51:24 EDT、`starsub03` から投入済み。SUMSの受付結果は `Scheduling successful`、
1 process / 1 file / `Success YES`。Condorへの登録も確認した。

| 項目 | 値 |
|---|---|
| Condor ID | `294294.0` |
| SUMS jobid | `8DDB3616BC656FBC9AD7E286A0026D51` |
| subjob | `8DDB3616BC656FBC9AD7E286A0026D51_0` |
| 出力予定ファイル | `auau13p5_anaLambdaNuclearId_pilot_20260930_8DDB3616BC656FBC9AD7E286A0026D51_0.root` |

投入したXML: `job/run/joblistlog/joblist_auau13p5_anaLambdaNuclearId_pilot_20260930_8DDB3616BC656FBC9AD7E286A0026D51.xml`。
設定記録: `job/run/configlog/config_auau13p5_anaLambdaNuclearId_pilot_20260930_8DDB3616BC656FBC9AD7E286A0026D51.txt`。
コード・ライブラリ・SUMS生成物の記録:
`job/run/runmeta/runmeta_auau13p5_anaLambdaNuclearId_pilot_20260930_8DDB3616BC656FBC9AD7E286A0026D51.json` と同じ識別子のsidecar/tar。

投入コマンド（実行済み。重複投入しないこと）:

```bash
cd job/run
./submit.sh ../joblist/joblist_auau13p5_anaLambdaNuclearId_pilot_20260930.xml
```

出力・stdout・stderrディレクトリは投入前に作成・検査済み。
SL24y/Singularityの増分 `make all`、XML構文、実commandの `tcsh -n` を確認した。
試験XMLには、生成後にROOTSYS明示、scratch検査、入力の第1列正規化・1件検査、
入力転送失敗時の停止、ROOT出力の存在・非空検査を加えている。
XMLを再生成する場合は、これらの試験用追加処理を失わないよう確認すること。

FileCatalogのイベント数欄は **3,852** だったが、投入前に同じXRootDファイルを
読み出して確認した `PicoDst` treeの実エントリ数は **3,202**。
SUMSの入力リストから第2列のカタログ件数を除き、実ファイル全体を読む設定。
Event / Track / TrackCovMatrix / BTofPidTraits branchの存在を確認済み。

**ここでは投入受付までを確認した。正常終了・ROOT回収・ヒストグラムQAは未確認。**
自動再投入・本投入・完了待機・自動mergeは行っていない。
後日は上記Condor/SUMS IDのstdout・stderrと、このディレクトリのROOTファイルを確認する。
SUMSのruntimeコピーはjob開始時なので、試験が終わるまでは関連ソース・設定・ライブラリの
変更を避ける。投入時の状態はrunmeta側にも保存されている。
