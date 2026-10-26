# Task 3再校正: ³He/⁴Heの共通log-polynomial係数fit

## 現在の状態

再検証・再校正中。[計画](task3_helium_logpoly_recalibration_plan.md)に沿い、表示修正と生分布・係数fitを分ける。現在の補間校正を支持された最終校正と仮定しない。

## 確認した問題

前回の補間providerはR<0.41で従来値、R≥0.41で最初のknotsを使用するため、境界で⁴He中心が104.765→82.021 keV/cm、約−22.745 keV/cm跳ぶ。これはデータ由来の物理的な急変でなく、domainの切り替えである。R=.405で従来106.058へ戻ることが、前回PDF p.4の最初の跳びを作っている。

205点の元独立fitは6品質通過・199参考で、raw局所最大があるのは31点。核種の位置・幅が正しく決まったとするには、生分布との一致・背景と重複の識別性を再確認する必要がある。欠落fit gapのlegacy切り替えも段差を作る。

## 進捗

| 項目 | 状態 |
|---|---|
| 上下で同じ物理範囲の低剛性図 | 完了、10ページ・40パネル、元fitは変更なし |
| 生p/q×dE/dx COLZと旧/現行曲線の重ね図 | 検証済み、別名の可読性改善版も完成 |
| ³He/⁴Heの共通係数によるfit | ROOT5 compile・数値80 tests PASS、20区間prototypeは実行完了／未収束 |
| 新係数の適用・ROOT5検証 | 未実施、候補の支持を確認してから判断 |

元ROOT・旧fit/校正/PDFは保持した。今回の再校正ではまだ現行mainconfや係数を変更していない。本投入・git操作は行わない。

## 下段dE/dxの表示修正

[全範囲を表示した低剛性fit PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/views_full_dedx/he4_low040_components_full_dedx.pdf)
を別名で保存した。上下とも従来校正のnσ=-8〜+4を物理dE/dxへ変換した範囲を使い、
保存軸0〜120 keV/cmの範囲に制限する。tritonと³Heの中心は10区間すべて下段で表示される。
元の`views/he4_low040_components.pdf`は保持した。

保存bin・count・error・座標変換の監査に加え、10ページを目視した。
ROOT5の表示上端の丸めでoverflowが含まれるケースはnormal binへの明示制限で解消した。
初回表示は`views_full_dedx_pre_axis_bin_guard`へ保持し、最終表示と区別する。

再描画では保存studyの旧mainconfと現在のmainconfのSHAが異なるため、
`--allow-runtime-mainconf-drift`を明示し、変更前後のSHAとruntime設定snapshotを記録した。
使用したfit数値・設定・入力fingerprintの一致検証は維持し、通常はSHA差異を拒否する。

## 独立検証に使う分布

既存の`true/hDedxP_3He_m2`は、TPC核種pullやbest-speciesによらずTOFの
raw `(m/q)^2`を1.78285〜2.06485で選んだ³HeタグのdE/dx分布である。
R=0.4〜0.5にも約144万countがあり、³Heの位置・幅のcross-checkに利用する。
ALL分布の部分集合なので、ALLとタグ分布を独立データとして同時に足すことはしない。
TOF受容性・qualityに由来するshape biasを確認する。

一方、この統合ROOTには独立した⁴He TOFタグ分布はない。低剛性の⁴He成分を
滑らかな式で描けても、それだけで核種帰属が確定したとは判断しない。

## ユーザーの終了指示による中断・明日の再開箇所

2026-10-13、ユーザーの「今日はここで終了、明日再開」に従い、新しいfit・解析を開始せず中断した。
全範囲PDFは完成済み。生分布の描画処理もexit 0で終了した。共通係数fitは準備のみで、
候補係数を新たに決定したとは扱わない。本番校正・mainconf・元ROOTは今回の再検証で変更していない。

### 保存した生分布の確認用出力

- [低剛性のALL COLZと旧/現行曲線](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_baseline/raw_low.png)
- [TOFで³Heをタグした低剛性COLZ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_baseline/tof3_low.png)
- [R=.425のALL/タグshape重ね図](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_baseline/all_tag_shape_R_425.png)
- [描画studyフォルダ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_baseline)
- [COLZ比較PDF・4ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_baseline/raw_dedx_calibration_overlay.pdf)
- [生分布とTOF射影PDF・25ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_baseline/low_rigidity_raw_projections.pdf)
- [ALL/TOF shape比較PDF・25ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_baseline/tof3_vs_all_shape_overlay.pdf)

COLZ4ページ・低剛性射影25ページ・ALL/タグshape重ね25ページを生成した。
`projection_audit.tsv`・`tof3_moment_summary.json`に全210区間のタグmode・窓内momentを保存する。
mode付近の窓内mean/RMSはGaussianの中心・幅そのものではなく、タグ・窓への条件付きの診断量である。
ALLとタグのcore近傍積分比でshapeを表示するscaleも記録する。これはTOF効率やCF背景係数ではない。
翌日は実画像、subset監査、保存精度、元データ不変・completionの記録を確認してから解釈する。
中断時の終端確認では`completed.json`の`phase=complete`・`inputPreserved=true`を確認した。
生成時に元ROOTのfingerprint不変と製品/sourceのSHAを記録している。追加の一括再hashは未実施。
今回のROOT確認・描画のsessionはすべてexit 0で終了し、共通fitの実行sessionはない。
描画用新ツール3件の実行説明は、翌日のREADME追記に残した。

### 共通係数fitの保存状態

- engine: `tools/he4_purity/global_logpoly/GlobalHeliumFit.C`。
  He3/He4共有14係数、dE/dx bin積分・R quadrature、AMPのKKT確認・step=0対策・非収束profile拒否、
  近傍核種の連続R域での中心順序確認まで保存した。まだcompile・実fitをしていない。
- launcher: `tools/he4_purity/global_logpoly/run_global_fit.py`、
  設定: `global_logpoly_settings.yaml`、単体テスト: `test_run_global_fit.py`。
  6単体テスト、実入力のprepare-only、3初期値の共通係数準備、既存outputを拒否する試験はPASS。
- [準備だけのstudy](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/launcher_prepare_001)
  はfit結果ではない。legacy・低R側shift・旧knotsの多項式回帰を初期値として保存した。
  old-knots回帰は真値・事前分布として拘束しない。feasibleなlegacyとのblend率を記録する。
- 初回prototype設定はR幅0.02、dE/dx幅0.5 keV/cm、低Rを密に含む20区間・legacy初期値・outer4回。
  先に数値テスト・少数区間の挙動を確認し、全210区間への実行はその後に判断する。

### 明日、実fitの前に修正・確認すること

1. engineへunit-count・variance・ALL/TAG軸一致と単独呼び出し時の既存出力保護を追加する。
   追加patchは審査timeoutで未反映。設定側の`countIntegerTolerance`、
   `countVarianceRelativeTolerance`、`axisTolerance`は保存済みだが、現engineではまだ使っていない。
2. 採用した最終係数と共分散の評価点を一致させ、AMP/KKT・local/global・outerの収束状態を出力する。
   conditional covarianceを、近傍核種・背景shapeまでprofileした完全な誤差と誤記しない。
3. launcherでrun後に元summaryのSHAと、実際にcompileしたsource snapshotのSHAも再確認する。
   現版は元source・入力の前後照合まで。既存outputは空でも拒否する。
4. 数値selftest・ROOT5 compile → 少数R prototype → 複数seed・背景・bin幅・R quadrature・held-out感度へ進む。
5. TOF³HeのtagとNOT_TAGを用いる場合はALL+TAGを二重計上しない。タグと全trackのshape一致も仮定しない。
   低Rの⁴He Gaussianが³Heの非Gaussian tailを拾っていないか、resolved Rからの外挿と比較する。
6. 十分な支持が得られた場合にだけ、別の係数YAMLを保存して直接式評価の解析profileを検討する。
   今の補間校正・legacy・従来profileへの切り戻しを保持する。

README/REFERENCEの共通係数fit入口の追記は途中。既存の全範囲表示説明・再検証の注意書きは保存した。
中断前の読み取り確認も一部審査timeoutで終了したため、未確認部分を「確認済み」にしない。

## 2026-10-14: 再開時の確認と実行環境の問題

再開し、raw studyの69製品・5 source snapshotのSHA一致と元ROOTのfingerprint不変を確認した。
低RのALL/TAG画像25区間も目視した。R=.4〜.5の前回⁴He成分は独立local最大でなく肩であり、
TOF³Heの非Gaussian tailとの重なりが残る。R≈.525以降は高側bumpが見える。
TOFの約144万countは全dE/dxの総数で、³He近傍診断窓は33,655、coreは17,318。
R=.425は窓917・core430である。core規格化後の高側TAG/ALL積分比は.425で.626、
.4〜.5合計で.508。高側excessをすべて⁴Heと断定しない。

launcherへ元summaryと消費するsource/settings snapshotの実行前後hash確認を追加した。
engineは入力監査等のhelper追加までで、最終mainループ・共分散/収束flagの反映は未完了。
現在のengine SHAは`d041e5349c34fa2565098d3fd8dadae1f40e0e3677fdf39eec69895828a8c89e`。
未compile・未fitで、現行校正・mainconf・元ROOTは今回変更していない。

通常コマンドは`bwrap: Creating new namespace failed: No space left on device`、
承認付きコマンド・編集も自動承認のtimeoutが続き、無害な`true`も完了できなかった。
安全拒否ではないが、実行環境が安定しないためprototypeを開始しない。
次回は編集の適用状態確認→入力guard/mainへの接続・採用gと共分散/収束flagの整合→
ROOT5 selftest/compile→少数区間prototypeから再開する。

### 同日の追加診断: 承認付き実行で再開可能

ユーザーから原因調査と実行の明示承認を受け、ホスト設定を変更せず切り分けた。
通常sandboxではシェル開始前に失敗する一方、承認付きでは同じ作業ディレクトリのコマンドが成功した。

| 確認 | 実測結果 |
|---|---|
| `user.max_user_namespaces` | 1500 |
| `user.max_net_namespaces` | 0 |
| `unshare --user /usr/bin/true` | exit 0 |
| `unshare --user --net /usr/bin/true` | exit 1、ENOSPC |
| `bwrap --unshare-user --ro-bind / / /usr/bin/true` | exit 0 |
| 上記bwrapに`--unshare-net`を追加 | exit 1、通常sandboxと同じエラー文 |
| STAR／ROOT5コンテナの`root-config --version` | exit 0、5.34/38 |
| 共通係数launcherの単体テスト | 6件PASS |

network namespaceの作成制限で同じエラーが再現するため、ROOTや入力PicoDstの変化ではなく、
通常実行ツールのsandbox作成側の問題と判断する。Codex実行ファイルにもuser/netのunshareフラグが含まれる。
ただし、この確認で内部の実起動argvを取得したわけではない。
承認timeoutは別の層の問題で、ログにdeadline超過が残るが内部遅延原因は不明。
明示承認後の呼び出しは正常終了しており、現在は承認付きの経路で作業を再開できる。

ディスク容量は/tmp空き約89 GB、GPFS空き約197 TB。GPFSの`df -i`／statfsのinode表示は
値が整合しないため「空き18」と断定しない。GPFSネイティブの最大inode値は896,925,696、
解析出力先への一意な空ディレクトリ作成は成功し、同じ空ディレクトリのみrmdirで除去した。
sysctl・既存プロセス・既存データには変更を加えていない。

engineの未完了部分を完成してから、同じmainconf・STAR/ROOT5で数値selftestとprototypeへ進む。

## 10月14日: 承認付きで行った数値検証と実データprototype

承認付き実行で、engineの未完了main接続と収束・共分散の整合を完成した。
hostのnamespace上限やCodex設定・既存プロセスを変更せず、同じユーザーでSTAR/ROOT5を実行した。

- launcherの8単体テストPASS。元ROOT・元summary・実compile source/settings snapshotの改変検知とnative R-edge整列を含む。
- [ROOT5数値selftest・80項目](validation/global_helium_numeric_20261014_001/selftest.json)と
  [来歴](validation/global_helium_numeric_20261014_001/provenance.json)を保存した。
  AMP/KKTのゼロ開始・ゼロ境界・近接成分・高count・反復上限、Gaussian尾積分、GL1/3/5、
  count/variance・軸・TAG subset・正値/核種順序を確認した。これは核種同定の検証ではない。
- frozen engine SHA: `3536a4466bc58070224a98b57f9c2d08747de57f21aad0618897e4bd399b2b66`。
  元の中間sourceの記録は上記に保持し、各studyはその時のsource snapshotを持つ。

### prototypeの結果と採用しない理由

[20区間prototype](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/prototype_20261014_001)
はexit 0／実行完了。R幅0.02、dE/dx幅0.5 keV/cm、legacy初期値、outer4回である。
元ROOTのfingerprint、元summaryとsource・snapshotの前後SHA、入力整数性・varianceを確認した。
`completed.json`は`physicsValidationComplete=false`であり、実行完了をfit収束と混同しない。

| 検証・量 | 結果 |
|---|---|
| AMP/KKT | 20区間すべて収束 |
| joint deviance | 14.527M → 11.539M → 8.753M → 6.847M |
| 最終global/local/outer stationarity | 全て未成立 |
| 最終scale正規化係数変化 | 0.681、許容0.001を超える |
| 最終global gradient/count | 0.00715、許容0.000001を超える |
| R=.41の⁴He中心／幅 | 約93.856／4.992 keV/cm、未収束の参考値 |
| R=2.49の⁴He幅 | 約0.127 keV/cm、fit bin幅0.5未満で要再検証 |
| 新係数の本番適用 | 行わない |

高Rのtritonはfit窓外で、R≥1.61の6区間では窓内期待収量がほぼゼロになり、
freeのT中心・幅が識別不能になる。低Rでもlocal/globalの改善余地が大きい。
従って、Hesse警告を単に無視して「fit成功」とすることも、今回の係数を正しいPIDとすることもできない。

次の検証は、Tshape固定の別variantによる数値原因の切り分け、反復・global予算と係数基底のconditioning、
低RのTshape感度、resolved-Rからの予測・保留R・背景tail/model感度である。
元の205独立点を平滑化して正解とみなす方法へは戻らない。

### 生分布の可読性改善版

[新しいCOLZ重ね図・4ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_readable_20261014/raw_dedx_calibration_overlay.pdf)、
[低RのALL/TOF射影・25ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_readable_20261014/low_rigidity_raw_projections.pdf)、
[ALL/TOF形状比較・25ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_overlay_readable_20261014/tof3_vs_all_shape_overlay.pdf)
を別名のstudyへ保存した。caption・軸title・指数の重なりを解消した。
元の曲線・窓・moment・subset数値は旧baselineとbyte-identicalで、描画だけを変更した。
新旧各69製品・5 source hashと元ROOTのfingerprint不変を確認した。
未収束のprototype係数をこの図や本番校正へ追加していない。

TOF³Heの条件付き窓momentを共通式へ回帰した
[参考cross-check](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/raw_validation_20261014/conditional_moment_crosscheck/conditional_moment_crosscheck.json)
も保存した。Rmin=.6／.75の回帰からR=.425へ外挿すると中心は72.718／79.515 keV/cmと変わる。
range内残差が小さくても低Rへ約6.8 keV/cmの外挿差があり、これもGaussian中心・幅の確定ではない。

次revisionではglobalの片側stationarity probeについて、真の関数制約境界とAMP profile失敗を
区別して記録する。現在のprototypeはstationarity=falseなので採用判断に影響しないが、
将来profile非収束を境界扱いしてtrueと判定しないための改善事項である。
