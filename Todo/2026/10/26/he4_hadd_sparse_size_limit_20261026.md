# 2026-10-26：最終haddの巨大Sparse書き出し上限

対象SUMS：04BF2DAE47164B6F040E3A460E539289、2021 P24iy。元26955 subjobの全270区間と、大きな3中間ROOTは読み取り/Entries/加算検証を通過。最終3本のhaddが03:48 JSTに終了した際、ROOTのTBufferFile::WriteByteCountが bytecount too large (more than 1073741822) を3回報告し、workerはfailed・03:49 JSTに要確認メールを送った。09:08 JSTの問い合わせ時点で最終ROOTは未公開、元ROOT移動も未実施。これは前回のTH1F丸め差の問題とは別である。

全キーの大きさをTKeyから読み、唯一大きなオブジェクトが8次元hNuclearPidPairJoint（THnSparseT<TArrayD>）であることを確認。未完成final.rootの当該キー非圧縮長は1,209,524,077 bytesで、単一オブジェクトのbyte count上限1,073,741,822を超える。中間0/1/2はそれぞれ941,829,809 / 941,472,648 / 841,684,938 bytesで上限内。各中間ROOTは正常読み取り検証済み。エラーの出たfinal.rootを完成品として公開しない。

根拠：job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289/root_size_diagnosis_20261026/keys.log、およびpartial_merge_native_verified_20261025/final/hadd.log。前回のTH1F問題は独立ROOT Addとの全bin一致検査を導入し解決済みであり、巨大オブジェクトを無視して公開する権限は与えられていない。

元データ・中間ROOTを全て保持したまま保存形式の選択をユーザーへ依頼。推奨は通常1786ヒストを既定名の1 ROOTへ統合し、巨大Sparseだけ内容・誤差を保った3個の別ROOTへ保存する方式。既存のhNuclearPidPairJointの読み取り先が変わるため、無断で削除・型変更・誤差配列削除をしない。代替は全ヒストを含む検証済み3中間ROOTを分割最終出力として保存する方式。現時点では返答待ちで、形式依存の再結合/公開は未実施。

Sparseを合算したキーは取得時にエラーを出す可能性があるため、未完成final.rootの内容読取りや再利用を試さず、TKeyメタデータだけを診断に使用した。解析コード・PID・閾値・元ROOTは未変更。追加メール・再投入・Git pushなし。


## 2026-10-26 09:24 JST：通常ヒスト＋3 Sparseの方式を実行

ユーザーが「通常1786ヒストを1本、巨大PIDヒストは3分割の別ROOT」を明示選択。保存方式の承認を受け、ROOT5 TFileMergerのkSkipListedでhNuclearPidPairJointのみを除いて3検証済み中間ROOTを統合する。Sparseは中間0/1/2からそのまま別ROOTへ書き出し、全binの座標・内容・誤差²、全軸/ラベル、Entries、Sumw2とmomentsの完全一致を再読取りで検査する。

6段階の小規模実ROOT試験（nested通常ヒストを残してSparseだけを除く、通常ヒスト加算、3 Sparseのweighted fill・flow・誤差保存）を通過。09:23:49 JSTにsplit_sparse_final_20261026/finish_split.pyをnohup・独立sessionで起動し、09:24に通常統合とSparse抽出開始を確認。最大同時ROOTは2、nice10。現時点は実行中で完了ではない。

通常ROOTの読取り・1786キー・Entries・入力/再構成イベント数、7差承認＋追加単精度差のROOT Add全bin再現検査を維持。元26955 ROOTのfingerprint不変を確認し、全4ファイルを検証してから公開する。個別ROOTの整理は通常ROOTと3 Sparseの保存を確認した後だけ実施。完了メールは「部分hadd完了・Sparse分割保存」とし全保存先を記載する。旧failed/送信記録、元データ、旧中間ROOTは保持。Git push・再解析なし。

出力形式・8次元ヒストの説明・射影方法：rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/README_outputs_04BF2DAE47164B6F040E3A460E539289.md。最終manifestと状態recordで完成を判定する。


## 自動完了記録 — 2026-10-26T11:13:22.408564+09:00

状態：完了（通常ROOT＋3 Sparse、個別ROOT保存済み）。メール：ローカルsendmail受付済み。
通常ROOT：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root`（14064010756 bytes）。
Sparse：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_pid_pair_joint_part000.root`。
Sparse：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_pid_pair_joint_part001.root`。
Sparse：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_pid_pair_joint_part002.root`。
検証・公開・保存・送信記録：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289/split_sparse_final_20261026`。
