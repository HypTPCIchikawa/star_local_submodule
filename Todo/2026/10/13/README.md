# 2026-10-13 STAR Todo の実施記録

依頼元：[STAR_Todo_20261013.txt](STAR_Todo_20261013.txt)。
各Taskの計画、結果、再現方法をこのディレクトリにまとめる。

本日の全体・設定・検証・MD/PDF索引・次回の再開箇所は[日次ログ](../../../../analysisnote/20261013/summary20261013.md)を参照。

| Task | 計画 | 結果・報告 |
|---|---|---|
| 1. Git管理・中間生成物の整理 | [計画](task1_git_cleanup_plan.md) | [報告](task1_git_cleanup_report.md) |
| 2. 高剛性⁴HeのTOF選別 | [計画](task2_he4_tof_pid_plan.md) | [報告](task2_he4_tof_pid_report.md) |
| 3. ⁴HeのTPC純度推定 | [計画](task3_he4_purity_plan.md) | [報告](task3_he4_purity_report.md) |
| 3追補. 低剛性ピークの再探索・linear表示 | [計画](task3_lowrigidity_refinement_plan.md) | [改訂結果](task3_lowrigidity_refinement_report.md) |
| 3再追補. 0.40 GeV/c開始のfit | [計画](task3_from040_refit_plan.md) | [3.85 GeVの結果](task3_from040_refit_report.md) |
| 3実装追補. ³He/⁴He校正とCF背景用スペクトル | [実装計画](task3_helium_pid_cf_implementation_plan.md) | [進捗・報告](task3_helium_pid_cf_implementation_report.md) |
| 3再校正. 広いdE/dx表示・共通式の係数fit | [計画](task3_helium_logpoly_recalibration_plan.md) | [再検証の進捗](task3_helium_logpoly_recalibration_report.md) |
| 3方針変更. 0.40開始free fitから反復して共通式へ | [計画](task3_helium_iterative_from040_plan.md) | [最新報告](task3_helium_iterative_from040_report.md) |

## 出力の置き場所

- Task 1：元リポジトリの外に復元可能な隔離領域を設け、報告にmanifestを記載する。
- Task 2：新しいmainconf・核PID YAMLは`config/`、検証ログは本ディレクトリの`validation_*.log`、実イベント検証ROOTは解析出力内の専用サブディレクトリに保存する。
- Task 3：処理マクロ・設定は`tools/he4_purity/`、PDF・数値表・入力来歴は各解析ROOTディレクトリの`pdf/he4_purity_20261013/`に保存する。

⁴Heの追加選別はYAMLで明示的に有効化する。3.85 GeVの入力は2021 P24iyに限定し、保存済みの統合ROOTを更新しない。

## 実施結果

1. 中間生成物・旧単発コード63件を復元可能な外部領域へ隔離した。既追跡15件の整理はGitHubへ通常push済み（cleanup `a9d63d6`、報告 `8f7d539`）。
2. `_tofExtended` mainconfで、R<2.0はTPCのみ、R≥2.0はTPC＋TOFを使う⁴He選別を追加した。STAR環境のビルド、両エネルギー各845 checks、2021 P24iy先頭1000イベントの出力検証が成功した。小統計では高Rの有効TOF候補がなく、実データでの高R収量増加はまだ確認できていない。
3. 旧設定では各エネルギーで210剛性区間の純度fit・背景比較・数値表を作成した。13.5 GeVは168区間でfit品質基準を満たすが、候補数coverageは34.94%。3.85 GeVは純度を参考値として示した。その後の低剛性・bin幅の見直し結果は下記追補を参照。TOF併用後・最終best-species後の確定純度ではない。

### 主要なPDF

| 対象 | 剛性対純度 | 背景モデル比較 | 全射影fit |
|---|---|---|---|
| 3.85 GeV・2021 P24iy | [1ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_purity.pdf) | [2ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_background_models.pdf) | [53ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_fits.pdf) |
| 13.5 GeV・2020 P24iy | [1ページ](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_purity.pdf) | [2ページ](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_background_models.pdf) | [53ページ](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_fits.pdf) |

今回pushしたのはTask 1の整理とその記録。Task 2/3のコード・設定・報告は作業ツリーに保存している。

## Task 3追補: 旧fitのlinear表示

旧fitを変更せず、保存済みの射影・fit関数から表示だけ追加した。各エネルギーの`views/`に独立保存している。

| 対象 | 全210区間のlinear表示 | 低剛性の物理dE/dx・log/linear比較 |
|---|---|---|
| 3.85 GeV・2021 P24iy | [53ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/views/he4_fits_linear.pdf) | [8ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/views/he4_low_rigidity_dedx.pdf) |
| 13.5 GeV・2020 P24iy | [53ページ](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/views/he4_fits_linear.pdf) | [8ページ](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/views/he4_low_rigidity_dedx.pdf) |

元のPDF・ROOTは保持する。新studyは`refined_lowrigidity/`に保存し、元校正の窓・fitによる新しい窓・数値的な収束・fit品質・ピークと肩を区別して記録する。

## Task 3追補: 改訂fitの結果

低剛性の信号平均制約、複数初期値、³He・tritonのtailを見直し、Float16保存による振動を抑えるためdE/dxを0.5 keV/cmへrebinした。R分割と元PID窓の観測countは維持した。両エネルギーで0.50–0.51 GeV/cの候補成分をfitできたが、3.85 GeVでは独立した生ピークではなく背景の肩として分離している。

| 対象 | 改訂fit・log | 改訂fit・linear | 低剛性dE/dx・log/linear |
|---|---|---|---|
| 3.85 GeV・2021 P24iy | [53ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_fits.pdf) | [53ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_fits_linear.pdf) | [8ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/views/he4_low_rigidity_dedx.pdf) |
| 13.5 GeV・2020 P24iy | [53ページ](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_fits.pdf) | [53ページ](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_fits_linear.pdf) | [8ページ](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/views/he4_low_rigidity_dedx.pdf) |

3.85 GeVは数値的な参考fitが200/210区間に増えたが、品質基準を通る区間は0。13.5 GeVは169品質通過・30参考・11不採用で、品質通過区間の候補数coverageは65.51%。背景関数間で純度が大きく違う区間もあるため、covariance誤差だけを精度と解釈しない。

各`refined_lowrigidity/`内に元窓／fit中心窓の比較CSVとエネルギー別のPID曲線候補YAMLを保存した。候補はstudy用で、productionのPIDヘッダー・設定・元ROOTは変更していない。fit後の校正候補保存で発生したPython依存エラーは標準ライブラリ化で解消し、再fitせず完成させた。その際も全fit ROOT・PDF・数値のSHA256不変を確認し、旧snapshotを`revisions/postprocess_repair_001/`に保持した。

## 0.40 GeV/c開始の再fit 3.85 GeV

最新の再fitは3.85 GeVの2021 P24iyのみを対象とし、`refined_lowrigidity_from040/`へ別出力した。triton・³Heの2つの背景Gaussianを含め、0.50未満の信号位置を探索する。0.42–0.48の6区間が品質基準を通過したが、最初の0.40–0.41は信号平均が下端に残り不採用、0.41–0.42と0.48–0.50は参考値として区別した。raw分布では背景上の肩としての推定であり、核種の独立確認や全sampleの確定純度ではない。

- [背景2成分・全分布と候補付近のlog/linear図 10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/views/he4_low040_components.pdf)
- [全210区間のlinear図 53ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/views/he4_fits_linear.pdf)
- [条件・全10低剛性区間の数値・採用可否・検証](task3_from040_refit_report.md)

元PID窓の候補数・全210区間のrebin countとvarianceは保存され、0.50以上の400 fitは前回と完全一致した。この再fit段階では旧結果・本番PIDを保持し、13.5 GeVは再処理していない。

## Task 3実装追補: ³He/⁴He校正とCF背景用スペクトル

今回の依頼により、fitで得た³He・⁴Heの中心値とGaussian幅を、3.85 GeV・2021 P24iyの現行⁴He解析設定に導入する。グローバルな従来PID係数は変更せず、元設定へ戻すmainconfも保存する。校正済み±2σ窓とbest-species条件でTPC純度を再評価し、旧窓の純度とは区別する。

⁴Heの既存SE/ME出力に加え、同じ解析ROOTへ実質量のΛ–³Heと、同じ³He候補を⁴He質量で計算する比較用スペクトルを追加する。後者は誤識別時のk*座標を調べる入力であり、自動的にCFを差し引くものではない。実装・検証の進捗と現在有効な設定は[実装報告](task3_helium_pid_cf_implementation_report.md)へ逐次記録する。

実装・検証は完了。新PID8,840項目、旧TOF両energy各845項目、人工pair・Maker policyguardと新旧同一1,000イベントの試験がPASS。最終新ROOTの各1,999,865項目の読み戻しでΛ選別不変・³He実質量/⁴He質量仮説・全cent/質量/m_T・mixing QAを検証した。

- [更新したPurity対p/q・背景関数・中心/幅の4ページPDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/calibrated_selection/purity_calibrated_windows.pdf)
- [3つの質量表現のCF入力QA（20 MeV/c、2ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/helium_calibrated_he3Control_20261013/helium_control_inputs.pdf)
- [新旧試験ROOTの条件・来歴](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/helium_calibrated_he3Control_20261013/README.md)

現行3p85⁴He設定には新校正と³He controlを適用済み。`_pidLegacy` mainconfを別に保存して従来条件へ戻せる。旧統合ROOTと元fit成果物は保持し、今回本投入・git pushはしていない。

## 最新の再検証: 独立点の補間から共通式の係数へ

ユーザーの指摘を受け、下段dE/dx表示を全範囲へ広げ、p/q×dE/dxの生分布で中心・幅を確認する。R=0.41の切り替えで⁴He中心に約22.7 keV/cmの人工的段差があり、前回の補間曲線を最終校正と扱わない。解析コードと同じlog-polynomialの係数を複数R区間で共通fitする再校正を開始した。[最新報告](task3_helium_logpoly_recalibration_report.md)を参照。旧結果と日次ログは当時の事実として保持する。

### 本日の中断・翌日の再開

ユーザーの終了指示に従い中断。全範囲の低剛性PDFと生分布の重ね描画は生成済み。
共通係数fitは実装・初期値準備までで、未compile・未実fit。本番PIDは今回変更していない。
入力精度/軸guard・収束/共分散の記録・来歴照合を完成させてからprototypeを始める。
詳細・新しい図のリンク・未解決項目は[中断時の再校正報告](task3_helium_logpoly_recalibration_report.md)へ保存した。

### 10月14日の再開状況

生図・元ROOTの来歴を確認し、低Rの³He/⁴He重なりとTOFタグ統計の解釈を再検証した。
launcherの来歴保護を追加、engineの入力監査helperを保存したが、mainへの接続・収束/共分散確認は未完了。
通常実行のnamespaceエラーと承認処理timeoutが続くため、ROOT5 compile・共通係数fitは未実行。
元ROOT・本番校正は変更していない。[最新状態・環境問題](task3_helium_logpoly_recalibration_report.md)を参照。

ユーザーの明示承認後に追加診断し、network namespace上限0で通常sandboxと同じ失敗を再現した。
ホスト設定を変えず、承認付きでROOT5.34/38の起動とlauncher6単体テストが成功。
現在はこの承認済み経路でengine修正・数値検証を再開する。承認timeout自体の内部原因は不明として区別する。

その後、engineのROOT5 compile・数値80項目／launcher8単体テストPASS、20区間prototypeの実行完了を確認した。
prototypeはAMPが収束してもglobal/local/outerは未収束のため、新係数は本番未採用。
生COLZ等は旧数値を保持した可読性改善版も別保存済み。[検証結果・次のstudy](task3_helium_logpoly_recalibration_report.md)を参照。

### 0.40開始free fitからの反復study（最新）

ユーザー指定の旧`he4_low040_components.pdf`は存在し、その保存fitを初期値にした。
旧legacy相対の制約を外した区間別free、soft平滑化、共有14係数のraw fitを実装・実行した。
低Rのnσ／元dE/dxのlog・linear 4面図と、全区間の各軸linear図を別studyへ保存した。
Python11、ROOT数値1,390、人工描画5 testsがPASS。実データの収束・核種同定とは区別する。
全域の共有式5追加条件は低R再現が悪く未採用。TOF対照だけでも⁴He中心は確定できない。
同じ式の低domain診断は成分を保持するが、120反復まで増やしても残差・数値収束が未達で未採用。
追加依頼の元dE/dx全域0〜120と補助45〜120の表示、生COLZへの中心・±σ／±2σの重ねも別出力へ保存・確認済み。
全曲線はdomain0.40〜0.70のみ。滑らかに見えても核種同定・purity・校正式の確定とは扱わない。
今回も旧PDF・元ROOT・本番PIDは変更しない。
図・数値・現在の採用判断は[反復study報告](task3_helium_iterative_from040_report.md)にまとめる。
