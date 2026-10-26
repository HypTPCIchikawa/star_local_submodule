# 2026-10-19 原PID関数固定・peak-anchored局所fit

[当日Todo・進捗](../../../Todo/2026/10/19/README.md)に対応する独立studyです。
本書は設計と実行手順を説明します。実行成功・数値収束・物理的妥当性・purityの認定は別の判定です。
出力フォルダが存在するだけでは完了と判断せず、各段階の`completed.json`と検査結果を確認してください。

## 方針変更と目的

旧局所fitでは、対象ROIから離れたtの裾をROI内で正規化し、その振幅だけを自由にすると、
実際のtピークと整合しない巨大な全線yieldでも対象領域の背景を表現できる問題がありました。
単に成分を増やす自由core＋tail実データpilotは、この問題の解決策としては採用せず、
[従来study](../coretail_20261019/README.md)から本方式へ方針を変更しています。
旧fit・設定・準備済み出力は履歴として残し、自動的に再開しません。

新方式では、同じALLヒストグラム内の**tピークそのもの、³Heピークと左側の肩、⁴Heピークと右裾**を
含む連続した領域を同時fitします。tの振幅は同じ観測ピークのcountによって制約され、
TOFタグの絶対収量をALLへ移すことはしません。

## 核種の中心・幅は全固定

`shapePolicy: fixed_original_pid_functions`を必須とします。
`StMaker/common/NuclearIdDeDxVsMom.h`にある元の`log10(R)`多項式
（中心：3次式、幅：2次式）を、保存済みの独立baselineと照合して使用します。

- d・t・³He・⁴Heの中心・幅は固定し、背景核種の中心や幅を都合よく動かしません。
- 各native R bin内の3点Gauss–Legendre nodeで、同じ元関数を評価します。
- TAG fitで得た各R独立の中心・幅や、折れ線補間した校正表は使いません。
- fitする核種はt・³He・⁴Heで、非負の振幅と明示した滑らかな背景を調整します。
- Gaussian／log-normalとlinear／exponential背景の4モデルを比較します。

log-normal対照では、元の数値をmode・局所曲率幅として入力します。
この幅はGaussianの標準偏差と同義ではなく、元PID関数がlog-normal応答を独立に検証したという意味でもありません。
実ピークとの不一致が残れば別途ピーク校正を検討し、滑らかな多項式係数として評価します。
このpilotから各Rの中心・幅を直接本番へ転用しません。

## 入力と領域

対象はAu+Au FXT 3.85 GeV・**2021 P24iyのみ**です。
既存部分統合ROOT由来の、count・error・flow保存を検査済みのnative ALL cacheを使用します。
元ROOTを再解析・変更せず、23B全統計を新たに処理したという主張もしません。

基準は`logpoly_baseline_001/baseline.json`に保存されたTPC-only・元多項式の独立baselineです。
準備時にbaseline、元header、runtime reader、設定、cache、抽出記録のSHAを確認します。
元多項式からbest-species gateを再計算しますが、今回のgate積分はraw候補数・残差の診断用です。

設定は[study_settings.yaml](study_settings.yaml)に明示しています。
代表Rは`0.825, 0.875, 0.925, 0.975, 1.725, 1.775, 1.825 GeV/c`の7点です。
各Rの3 node全体を覆うように領域を求め、0.5 keV/cm単位のnative bin境界へ外向きに揃えます。

| 領域 | 下端 | 上端 | 目的 |
|---|---|---|---|
| 共通fit領域 | tと³Heそれぞれの中心−2幅の小さい方 | ⁴He中心＋4幅 | 実際のtピークからHeピークまで同じ振幅でfit |
| He表示・診断領域 | ³He中心−2幅 | ⁴He中心＋4幅 | ³Heの左肩と両Heピークを確認 |
| t表示・診断領域 | t中心−2幅 | t中心＋2幅 | tピークと元PID関数の整合を確認 |

nodeの最小・最大を使うため、表の「中心・幅」は単一のbin中心値だけを指すものではありません。
4モデル×7点＝28 fitを保存し、失敗や不良fitも省略しません。
joint fitの実際の自由度と、重なり得る各診断領域のPearson寄与を区別します。
各領域に独立した自由度を割り振ってreduced χ²を作ることはしません。

## d除外の診断とpurityの限界

dは明示的にfit対象から除外し、その元関数から有限軸・各ROI・⁴He gateへの寄与上限を診断します。
ALL軸内観測countに`5√N`を加えたmarginを用いますが、これは**近似的・記述的なmargin**であり、
厳密な信頼区間でもfit振幅への強制制約でもありません。
寄与の大きさを確認する基準`materialOmissionFraction: 0.01`は診断用です。
元の応答関数自体の不確かさを網羅したsystematic errorや、purityの認定条件ではありません。

除外dの寄与が無視できない場合やfitが合わない場合は、その理由を保存します。
単に背景へ吸収されたことをもって良い物理fitとは判断しません。
今回の`selection.available`は`false`、purity等は`null`であり、未評価を0%に置き換えません。
新しいcutoff・PID校正・CF補正の本番採用には、応答・背景・帰属・不確かさの検証が別途必要です。

## ツールと出力

| ファイル | 内容 |
|---|---|
| [prepare_peak_anchored.py](prepare_peak_anchored.py) | `prepare`／`run`／`check`。入力とsourceをsnapshotし、生成TSVもSHAで保護。 |
| [FitPeakAnchored.C](FitPeakAnchored.C) | 元PID関数の中心・幅を固定した共通ROI fit、領域別残差・省略d診断・ROOT読み戻し。 |
| [render_peak_anchored.py](render_peak_anchored.py) / [PlotPeakAnchored.C](PlotPeakAnchored.C) | 保存済み結果からHe領域とt領域をlinear／log表示。再fitしない。 |
| `test_*.py` | 設定契約、元多項式、native領域、固定shape、上書き防止、描画等の検査。 |

study直下の`plan.json`に設定・7点の領域・入力SHAを保存します。
4つのモデル別サブフォルダには`settings.tsv`、`shapes.tsv`、`gates.tsv`、`regions.tsv`、
ROOT／JSON、実行ログと個別完了記録を保存します。
元count/error/flow・成分和・領域mask・残差の閉じたROOT読み戻しを検査します。
描画は別フォルダへ保存し、線形図はHeまたはtのピーク付近に限定します。
0–120 keV/cm全軸のlinear図でピークが潰れる表示へ戻しません。

## 実行手順

リポジトリ直下から実行します。`NEW`は未使用の名前へ置き換えてください。
既存出力への上書きや失敗jobの同じ場所での自動再実行はしません。
`run`は明示mainconfのbatch相当ROOT5 wrapperを使用し、farm jobは投入しません。

```bash
python3 -B tools/nuclear_pid_cutoff/peak_anchored_20261019/prepare_peak_anchored.py prepare \
  --settings tools/nuclear_pid_cutoff/peak_anchored_20261019/study_settings.yaml \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_logpolyLowRReference20261019.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/peak_anchored_pilot_NEW

python3 -B tools/nuclear_pid_cutoff/peak_anchored_20261019/prepare_peak_anchored.py run \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/peak_anchored_pilot_NEW

python3 -B tools/nuclear_pid_cutoff/peak_anchored_20261019/prepare_peak_anchored.py check \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/peak_anchored_pilot_NEW

python3 -B tools/nuclear_pid_cutoff/peak_anchored_20261019/render_peak_anchored.py \
  --study-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/peak_anchored_pilot_NEW \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_logpolyLowRReference20261019.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/peak_anchored_views_NEW
```

本番mainconf・本番PID係数・R veto・元ROOT・過去PDFは変更しません。
完了manifestは数値実行と保存検査の記録であり、物理QAの代わりではありません。

## 2026-10-19の実行結果

`peak_anchored_pilot_001`で4実行28fitを完了しました。数値収束28/28、全shape固定、
ROOTのcount/error/flow・成分和・mask・残差の再読込みと入力SHA不変を検査済みです。
`peak_anchored_views_001/peak_anchored_diagnostics.pdf`は28頁・112panelの表示検査と全頁readback済み。
prep8件・描画17件の計25 Python試験もPASSしました。

元曲線を固定するとt実peakからの規格化で旧遠方tailの巨大な外挿は避けられましたが、
joint χ²/ndfは約1887–6928で、元の中心・幅とデータの不一致が明瞭に残りました。
28件の数値収束を良好な物理fitと取り違えず、全件purity未認定のままです。
この出力をpurityの系統差やcutoff採用へ使いません。
次は核種別実peakを独立に校正し、同じ多項式の係数として検証します。
He4 fit内でt／³Heの中心・幅を自由にする方式には戻しません。

[結果・PDFへのリンク](../../../Todo/2026/10/19/report_pid_quality_20261019.md)を参照してください。
