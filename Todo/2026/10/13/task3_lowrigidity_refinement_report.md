# Task 3追補: 低剛性⁴Heの探索・校正候補・linear表示

現行コードへの³He・⁴He校正導入と³He CF入力の追加は[実装追補](task3_helium_pid_cf_implementation_report.md)で管理する。本報告内の多項式候補は比較用のstudy成果物のまま保持し、解析では0.40開始の採用205点の補間校正を使う。

この報告は探索開始0.50 GeV/cの結果を記録している。3.85 GeVの0.40開始の最新再fitは[別の追補報告](task3_from040_refit_report.md)を参照。以下の旧結果・図・設定は比較用に保持した。

前回の低剛性での「有効なfitなし」は「⁴Heがいない」という意味ではなかった。信号平均の狭い制約、近傍核種の境界、fit範囲外に出たtritonのtailの除外に加え、PicoDstのdE/dx保存量子化と0.1 keV/cm binの不整合が原因に含まれていた。本追補では現在のPID曲線の近傍を探索し、元の±2校正gateを維持した推定と、fit平均・幅に移動したgateの参考推定を分けて保存する。

## 保存場所

旧Task 3のファイルは保持し、各エネルギーで以下のサブフォルダを使用する。

- [3.85 GeV追補フォルダ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity)
- [13.5 GeV追補フォルダ](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity)

各フォルダに純度PDF、linear/exponential背景比較PDF、全射影のlog PDFとlinear PDF、CSV、監査JSON、ROOT、設定・ソースsnapshotをまとめる。raw dE/dxのnative射影も保持する。3.85 GeVは2021 P24iyのみ、13.5 GeVは2020 P24iyで、元の部分統合ROOTを読み取り専用で使用した。

| 対象 | Purity（1ページ） | Background比較（2ページ） | log fit（53ページ） | linear fit（53ページ） |
|---|---|---|---|---|
| 3.85 GeV | [PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_purity.pdf) | [PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_background_models.pdf) | [PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_fits.pdf) | [PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_fits_linear.pdf) |
| 13.5 GeV | [PDF](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_purity.pdf) | [PDF](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_background_models.pdf) | [PDF](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_fits.pdf) | [PDF](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_fits_linear.pdf) |

- 3.85 GeV: [全modelの平均・幅・両gate CSV](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_locus_fit_candidates.csv) / [校正候補YAML](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_pid_curve_candidate_3p85.yaml)
- 13.5 GeV: [全modelの平均・幅・両gate CSV](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_locus_fit_candidates.csv) / [校正候補YAML](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/he4_pid_curve_candidate_13p5.yaml)

低剛性を物理dE/dx軸で見やすくしたlog/linear対比（各8ページ）も追加した：[3.85 GeV](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/views/he4_low_rigidity_dedx.pdf) / [13.5 GeV](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity/views/he4_low_rigidity_dedx.pdf)。保存済みfitの再描画で、再fitしていない。元の旧fitにも`views/he4_fits_linear.pdf`を追加し、旧PDF・数値は保持した。

## 旧fitの問題と今回の対処

### 平均制約と近傍核種

旧fitは⁴He平均を元校正の±0.8 nSigma内に制限していた。しかしR≈0.50–0.55 GeV/cでは実データの成分平均がこの範囲を越える。新studyではR=0.50–0.85のみ平均の許容範囲を−3.2～+1.5元校正nSigmaに広げ、−3.0～+1.0を0.4刻みの初期値で走査する。幅・初期値・許容範囲は新しいstudy YAMLに記録した。

³Heの大きいピークを⁴Heと取り違えないよう、信号平均を近傍核種の平均より最低0.8元校正nSigma上に保つ。各componentの核種・平均・幅・初期値を保存し、信号なし背景modelとのAIC差も計算する。

近傍Gaussianは平均がfit範囲外でも、そのtailが入る場合に残す。例えばR≈0.70ではtriton平均は範囲−7～+4より外にあるが、無視できないtailが範囲内へ入る。tail-only背景成分の平均が境界へ寄る場合は、信号の平均・幅が内部で数値的に安定していれば参考component値を記録する。ただし品質基準通過にはしない。信号境界や範囲内の近傍ピーク境界は引き続き採用しない。

### dE/dxの保存格子とalias

ROOT Float16 round-tripと2021 PicoDstの782 trackで、70–90 keV/cm付近のdE/dxは1/64=0.015625 keV/cmの格子に載ることを確認した。元の0.1 keV/cm binには6個または7個の格子点が入り、平坦な分布でも約16.7%の周期差を生む。

- 人工一様分布の6/7 binの比は1.1664。
- 実際のR=0.505のnative射影の比は1.168。
- 格子点容量を補正した診断では周期ripple RMSが7.77%→0.81%。
- 0.5 keV/cmまたは1.0 keV/cm幅では人工一様分布のこのaliasが消える。

証拠は [Float16監査JSON](validation/float16_dedx_alias_002/float16_alias.json) に保存した。新studyではfactor 5、すなわち0.5 keV/cm幅にまとめる。期待bin数は解析的Gaussian・背景積分で計算し、元のcandidate数は保持する。観測gate候補数はnative binで求め、連続な元校正±2の積分gateは変更していない。

R=0.505の3.85 GeV代表試験ではχ²/ndfが約111→14.6に改善した。ただし他の剛性でも背景・Gaussian形状が完全に説明できるようになったわけではない。native binの大きいχ²を全て物理的な形状不一致と見なす旧解釈は不十分だったが、量子化だけで全問題が解決するとも言えない。

## ピーク・肩・fit品質を分ける

raw射影に局所最大があるか、Gaussian成分を含むモデルが信号なしモデルより適合するか、fitが数値的に収束するか、残差の品質基準を満たすかを別々に記録する。

3.85 GeVのR=0.50–0.51は、rebin後のraw射影に独立した局所最大がなく、背景上の肩として分離する。Gaussian成分がfitできたことを「独立したraw peakを観測した」とは書かない。条件付きの振幅／covariance誤差やnull-model AIC差も、厳密な粒子発見の有意度ではない。

原子核ラベルは校正近傍・平均の順序を使ったcomponentの解釈であり、低剛性でも背景modelへの依存や³Heとの混合を考慮する。新しいTOF選択後や最終best-species選択後の確定純度には使用しない。

## 全剛性fitの到達点

| 対象 | 品質基準通過 | 数値的に成立する参考fit | 不採用 | 品質通過区間の元gate候補数coverage |
|---|---:|---:|---:|---:|
| 3.85 GeV | 0/210 | 200/210 | 10/210 | 0% |
| 13.5 GeV | 169/210 | 30/210 | 11/210 | 65.51% |

3.85 GeVでは0.50 GeV/cから低剛性のcomponentを記録できるようになったが、高統計の残差・tail-only背景の制約が残るため参考値である。13.5 GeVは区間数だけでなく、品質通過区間が覆う候補数も前回の34.94%から65.51%へ増えた。主表示に選ばれた数値的に有効なmodelは、3.85 GeVでlinear 76／exponential 124、13.5 GeVでlinear 71／exponential 128（参考fitも含む）。

保存母集団と元gateは同じで、全210剛性区間のnative gate候補数が旧結果と一致した。今回の改善によってproduction解析の選択track数を変更したわけではない。元gate候補数の合計は3.85 GeVで198,592,594、13.5 GeVで734,686。

代表低剛性区間の純度は以下。いずれも元校正±2窓内のPである。

| 対象 | R [GeV/c] | linearのP | exponentialのP | 主表示 | 解釈 |
|---|---|---:|---:|---|---|
| 3.85 GeV | 0.50–0.51 | 65.85% | 62.61% | exponential | 独立raw最大なし、肩のcomponent参考値 |
| 3.85 GeV | 0.55–0.56 | 58.71% | 43.43% | exponential | model差15.28 pp、参考値 |
| 3.85 GeV | 0.60–0.61 | 49.07% | 53.85% | linear | 参考値 |
| 3.85 GeV | 0.65–0.66 | 73.10% | 76.01% | linear | 参考値 |
| 3.85 GeV | 0.70–0.71 | 91.86% | 90.38% | linear | tail-only tritonの境界あり、参考値 |
| 3.85 GeV | 0.75–0.76 | 95.98% | 95.87% | linear | tail-only tritonの境界あり、参考値 |
| 13.5 GeV | 0.50–0.51 | 68.69% | 67.14% | linear | 今回の品質基準通過 |
| 13.5 GeV | 0.55–0.56 | 64.99% | 62.90% | exponential | 今回の品質基準通過 |
| 13.5 GeV | 0.60–0.61 | 63.87% | 38.70% | exponential | 両modelとも基準内だがmodel差25.17 pp |
| 13.5 GeV | 0.65–0.66 | 51.78% | 54.43% | linear | 今回の品質基準通過 |
| 13.5 GeV | 0.70–0.71 | 65.54% | 68.20% | linear | 今回の品質基準通過 |
| 13.5 GeV | 0.75–0.76 | 80.44% | 80.84% | linear | 今回の品質基準通過 |

「品質基準通過」はχ²/ndf≤5等の今回のmodel検査基準を指す。13.5 GeVの0.60–0.61のように2関数とも基準内でも純度が大きく異なるため、covariance誤差だけを精度と解釈できない。P≈1で背景振幅が非負制約の下端へ寄る区間では、背景の片側上限もcovariance誤差から十分評価できない。

## 元gateと校正候補gate

主表示の純度は引き続き元のPID曲線のnSigma±2内のS/(S+B)。fit平均・幅が移動しても、この窓を黙って移動していない。

別列でfit平均±2×fit幅のgateも計算する。これは校正を変更した場合の参考値で、元gateの純度と同じ物理選択ではない。新gateがfit範囲を越える場合は外挿による純度を採用しない。

代表試験の13.5 GeV、R=0.50–0.51では、旧校正の平均約84.65 keV/cmに対し、Gaussian成分の平均は約71.56～71.63 keV/cm、幅は約5.75～5.81 keV/cm。信号平均と³He平均の分離は約2.57元校正nSigmaで、平均制約の内部に収まった。元gateでの純度は約67～69%である一方、移動したgateでは約45%。単に校正窓を動かすと同じ純度を保ったまま効率だけ改善する、とは限らない。

`he4_locus_fit_candidates.csv` に各modelの平均・幅・初期値、raw peak判定、null-fit比較、元gateと提案gateの端・純度・model signal数を保存する。`he4_pid_curve_candidate_<sample>.yaml` はエネルギー別の経験的平均・幅とlog10 polynomialの**study候補**であり、`productionEnabled: false`。元のMaker、共通PIDヘッダー、production YAMLは変更しない。候補曲線の観測範囲外への外挿も承認していない。

R=0.50–0.51における両energyの校正候補とgateの違いは次のとおり。

| 対象 | model | fit平均 [keV/cm] | fit幅 [keV/cm] | 元gateのP | fit平均±2fit幅のP |
|---|---|---:|---:|---:|---:|
| 3.85 GeV | linear | 74.479 | 6.036 | 65.85% | 39.73% |
| 3.85 GeV | exponential | 74.657 | 5.899 | 62.61% | 38.78% |
| 13.5 GeV | linear | 71.560 | 5.808 | 68.69% | 45.45% |
| 13.5 GeV | exponential | 71.629 | 5.748 | 67.14% | 45.18% |

元曲線の平均は両energyとも84.65 keV/cm付近であり、low-Rでは実componentの平均と差が大きい。平均・幅を移動するとsignalと背景の両方の収量が変わる。この候補を「旧cutと同じ純度・効率のまま正しい校正へ置き換えられる」と扱わず、曲線・生データ・model差を確認する必要がある。

## 検証

- [全210区間×2エネルギーの独立rebin監査](validation/rebin_audit_final_manifest.json)で、各fit binがnative binの厳密count和・variance和に等しく、flow・entries・内部名／ROOT keyも正しいことを確認した。旧／新JSON／独立TSVのgateCounts・fitCountsは全区間で完全一致し、検査エラーは0。
- 旧結果、改訂結果、追加図の完成manifestに記載した全productのSHA256一致と原本input不変を確認した。両エネルギーの改訂fitはlog53ページ＋linear53ページ、追加図はlinear53ページ＋物理dE/dx対比8ページ。低剛性0.505／0.555／0.605の画像を確認した。
- 低剛性R=0.55/0.70、linear/exponential背景、元平均制約を越えるGaussian信号の4つの人工例で、平均・幅、元gate・移動gateの純度、近傍核種の順序、信号なし比較を検証。
- signalなし例と弱い単調肩の例では、raw局所最大・component evidenceの誤検出がないことを確認。
- 旧機能を無効にした既存の2つのclosureも再検証し、旧数値を回収。

結果は [新6caseのclosure](validation/lowrigidity_toy_final)、[旧closure回帰](validation/legacy_toy_after_lowrigidity) に保存した。人工例の成立は実データの背景関数の正しさや同位体ラベルの完全な同定を保証しない。

最終C++ fitソースのSHA256は `6b664ce5ea97e62314ea7786ac85d4a444e643ef7f74cb8e95291ece885a8c70`。両energyのROOT fitをこのソースで完了した後、host NumPyのFortran依存による校正候補出力のimport失敗を修正した。小さい候補多項式fitを標準Pythonの再直交化QRへ置き換え、既知の3次・2次多項式を1e−9以内で回収した。ROOT fit・全fit PDFは再実行せず、入力・設定・ソースhashを確認した `--finish-existing` で後処理だけを再開し、全fit生成物のhash不変を確認した。修復前のcontext・wrapper・候補出力コードは `revisions/postprocess_repair_001/` に保存した。

原本ROOTは変更していない。解析jobやproduction PIDの再投入・書換えは行っていない。

設定は [purity_lowrigidity_settings.yaml](../../../../tools/he4_purity/purity_lowrigidity_settings.yaml)、実装前計画は [task3_lowrigidity_refinement_plan.md](task3_lowrigidity_refinement_plan.md)。
