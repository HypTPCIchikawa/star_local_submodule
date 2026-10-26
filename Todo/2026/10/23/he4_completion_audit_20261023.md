# 2026-10-23：修正版⁴He解析の終了状態確認

対象SUMS：`04BF2DAE47164B6F040E3A460E539289`。11:37 JSTに監視スクリプトのread-only --checkで再確認：submitted/terminal各26994、現在queue0。バッチは全件終了しているが、全件正常解析・hadd完了ではない。

26994 stdoutの末尾を照合し、26955件がstatus0、requested=completedを確認。成功分の実処理イベント数合計2,100,726,543、再構成2,067,647,474。成功件の解析RealTime中央値417.499秒・最大782.56秒。全解析summaryのCPU合計2917.256時間。Condor開始・終了イベントから最大同時実行2156件、最初submitから最後terminateまで4時間46分06秒。以前2日との条件・待ち時間の差は未比較で、短縮理由の断定はしない。今回の設定は2021 P24iy、maxEvents=-1。ログによる全予定entry処理の確認であり、全ROOTの再帰読み取り・物理QAは未実施。

39件は未成功：34件unreadable PicoDst、2件missing PicoDst tree、1件empty PicoDst、2件はzlib/TBasket展開エラーで途中停止。途中停止はsubjob362（103029予定/102449処理）と24813（54828予定/53765処理）。元ファイル破損とstaging時の問題の切り分けは未完了。Condor exit0でもSUMS外側wrapperが解析失敗を隠す既知の問題は継続。

メールの「統合入力0/26994」は解析イベント0という意味ではない。39件の不足・未完了によりstrict全件結合条件でhadd開始前に停止したため、haddに実際に渡した入力が0。残り26955 ROOTはヘッダー確認済みで回収済み。欠損を自動除外せず、全件hadd・元ROOT整理・成功通知を実行していない。今回の調査では再投入・部分hadd・追加メール・Git pushは行っていない。待つだけでは39件は解消しない。

記録：`job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289/result.json`、`/tmp/yichikawa/sched04BF2DAE47164B6F040E3A460E539289.condor.log`、対象stdout/stderr。既存の失敗結果・送信記録は変更しない。


続報：既存2021統合ROOTとの入力イベント数差は−0.00261%。ユーザーの条件付き許可に従い39件を除く部分haddを開始した。[比較・実行記録](he4_event_comparison_partial_hadd_20261023.md)。
