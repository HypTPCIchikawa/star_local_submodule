# ⁴Heの低剛性fit: 0.40 GeV/cからの再解析（3.85 GeV）

2021 P24iyの3.85 GeV保存ROOTに対し、⁴He候補Gaussianと低いnSigma側のtriton・³Heの2つの背景Gaussianを用い、位置探索を0.40 GeV/cから行った。0.41–0.50区間のcomponentを回収でき、0.42–0.48では今回の残差品質基準を通るfitも得られた。一方、最初の0.40–0.41は信号平均が許容範囲の下端に寄るため、fit曲線を保存するが純度を採用しない。

独立したrawピーク、背景上の肩をGaussianで分離したcomponent、品質基準を通る純度を区別する。今回の低剛性raw分布では⁴He候補側に独立した局所最大はなく、肩のcomponentとしての推定である。

## 設定・原本と保存場所

- [実装前計画](task3_from040_refit_plan.md)
- [新study設定](../../../../tools/he4_purity/purity_lowrigidity_from040_settings.yaml)
- [新結果フォルダ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040)
- [前回の0.50開始結果](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity)

入力は2021 P24iyのみの元の部分統合ROOT。解析済みのpre-PID `true/hDedxP`を使用する。対象はcharge-inclusive・accepted centrality全体のTPC窓であり、最終best-species、TOF後、Λとのペアに条件づけた⁴He純度ではない。13.5 GeVは今回再解析していない。

| 内容 | PDF |
|---|---|
| 元校正±2窓の純度 | [he4_purity.pdf](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/he4_purity.pdf) |
| linear/exponential背景比較 | [he4_background_models.pdf](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/he4_background_models.pdf) |
| 全210区間のlog fit | [he4_fits.pdf](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/he4_fits.pdf) |
| 全210区間のlinear fit | [he4_fits_linear.pdf](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/he4_fits_linear.pdf) |

- [低剛性2背景ピークの個別表示（10ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/views/he4_low040_components.pdf): 各区間で全nSigma分布のlog／linearと、⁴He候補付近のraw dE/dx拡大のlog／linearを4面表示。triton・³He・⁴He候補・背景和を別の線で表示する。
- [保存fitを読み直したlinear全区間表示（53ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/views/he4_fits_linear.pdf): 数値を再fitせず保存済みfitから描画。最初の不採用modelと保存上端によるfit範囲clipも明示する。

[平均・幅・元gate／移動gateのCSV](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/he4_locus_fit_candidates.csv)、[校正候補YAML](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040/he4_pid_curve_candidate_3p85.yaml)も保存する。この多項式候補はstudy専用のまま保持する。今回の依頼に伴う解析コードへの導入は、別途作る205点の線形補間校正を使い、[校正・CF入力の実装報告](task3_helium_pid_cf_implementation_report.md)で管理する。

## 前回との変更

前回は0.50未満で位置探索を行わず、信号平均を元校正の±0.8 nSigmaに制限していた。triton・³Heは既に含まれていたが、信号が低い側へずれると、³He側の幅や平均がそのcomponentを吸収して境界へ寄っていた。

今回、探索開始を0.40へ変更した。同じ−3.2下端での最初の試験では0.405–0.425の平均が下端へ寄った一方、0.435以降の平均が−3.1から−2.2へ連続的に動くことを確認した。そのため0.40≤R<0.50だけ、平均の許容範囲を−3.8～+1.5、初期値走査を−3.6～+0.8（0.4刻み）とした。raw局所最大の診断窓もこの範囲全体を調べる。

0.50以上の平均制約・初期値・幅・背景の扱いは前回と同一。0.40–0.50にも既存の低剛性探索で用いた幅0.35～2.2元校正nSigmaを適用する（前回の探索対象外では上端1.8だった）。近傍成分の平均移動係数・幅係数はそのまま使い、平均順序を保つため³He等の平均上端は今回の信号平均の許容下端から制限する。fit要求範囲−7～+4、Rの0.01 GeV/c幅、dE/dxの0.5 keV/cm幅、元校正±2gateは変更していない。平均をさらに広げて最初の区間を強制採用することはしていない。

⁴He候補平均と全近傍成分平均の間は最低0.8元校正nSigmaを要求する。0.40–0.50のモデルにはtriton・³Heの2つの主背景ピークが含まれ、deuteron Gaussianはfit範囲へのtail条件に該当せず含まれない。smoothな背景関数は別にlinear/exponentialを比較する。

## 低い側の2背景ピークと⁴He候補

以下のuは元の⁴He校正で定義したnSigmaであり、fitで中心を移動した後のnSigmaではない。backgroundのピーク位置と⁴He候補の平均を別に保存する。主表示モデルはこの低剛性範囲ではexponential。

| R [GeV/c] | triton平均u | ³He平均u | ⁴He候補平均u | linearの元gate純度 | exponentialの元gate純度 | exponential χ²/ndf | 判定 |
|---|---:|---:|---:|---:|---:|---:|---|
| 0.40–0.41 | −6.549 | −5.193 | −3.800 | 不採用 | 不採用 | 16.67 | 信号平均が下端、採用不能 |
| 0.41–0.42 | −6.608 | −5.084 | −3.507 | 64.32% | 53.50% | 7.46 | 参考値 |
| 0.42–0.43 | −6.632 | −5.013 | −3.279 | 66.62% | 59.81% | 4.33 | exponentialのみ品質基準通過 |
| 0.43–0.44 | −6.593 | −4.944 | −3.083 | 67.68% | 62.64% | 3.29 | 両modelが品質基準通過 |
| 0.44–0.45 | −6.543 | −4.865 | −2.924 | 67.88% | 63.44% | 3.08 | 両modelが品質基準通過 |
| 0.45–0.46 | −6.511 | −4.821 | −2.783 | 67.64% | 63.81% | 2.72 | 両modelが品質基準通過 |
| 0.46–0.47 | −6.491 | −4.784 | −2.668 | 67.68% | 64.28% | 2.94 | 両modelが品質基準通過 |
| 0.47–0.48 | −6.471 | −4.748 | −2.572 | 67.92% | 64.86% | 2.84 | 両modelが品質基準通過 |
| 0.48–0.49 | −6.482 | −4.741 | −2.403 | 67.68% | 64.66% | 6.08 | 参考値 |
| 0.49–0.50 | −6.486 | −4.727 | −2.224 | 66.90% | 63.73% | 8.79 | 参考値 |

0.41–0.42では背景model差が10.82 percentage pointsあり、covariance由来の小さい誤差だけを精度と見なせない。0.42–0.43もlinearではχ²/ndf=5.21で今回の基準を満たさず、両関数の差6.81 ppはmodelの評価を分けて読む必要がある。品質基準通過はmodelが唯一正しいことや核種ラベルの独立確認を保証しない。

0.42–0.43のexponentialでは⁴He候補平均約81.42 keV/cm、幅約7.84 keV/cm。一方、元校正平均は101.05 keV/cmで、元gate下端89.08 keV/cmは候補Gaussianの平均より高い。この結果は元gateが信号の高いdE/dx側を選んでいることを示す。fit平均±2fit幅へgateを移すと参考純度は約42.04%となり、元gateの59.81%と同じ選択ではない。

## 保存上限と最初の区間

dE/dx保存範囲は0–120 keV/cm。最初の区間中心R=0.405では元校正平均106.06、幅6.265 keV/cmで、元±2gateは93.53–118.59 keV/cm。元gate全体は保存範囲内だが、fit上端は+4ではなく実軸上端u=2.225までになる。

R=0.425の軸上端はu=3.164。厳密なR=0.400では元gate上端が120をわずかに超えるが、今回の最初の点は0.40–0.41 binの中心0.405であり、両者を混同しない。120以上のoverflowにはdE/dxごとの情報がないため、その部分の形状fitや外挿は行わない。

0.405では平均が−3.8下端に寄るので、signalなしmodelとの差や数値的収束があっても採用不能である。背景上の肩を「独立した⁴Heピークを確認した」と断定せず、さらに物理的な背景形状・校正を調べる余地を残す。

## 全区間と検証

最終再解析は完了した。

| 項目 | 結果 |
|---|---:|
| 品質基準通過 | 6/210 |
| 数値的に成立する参考fit | 203/210 |
| 不採用 | 1/210（0.40–0.41） |
| 主表示linear／exponential | 76／133 |
| 元gate候補数の合計 | 198,592,594 |
| 品質基準通過区間の元gate候補数 | 10,850,976（全体の5.46%） |

前回の3.85 GeVは品質基準通過0、参考200、不採用10だった。今回は新しい低剛性で6区間が品質基準を通るが、それでも候補数の5.46%だけであり、全sampleの確定純度や平均純度は出せない。今回の結果をcompleted CFへ一律の⁴He純度として差し替えない。

独立監査で全210射影のnative→factor 5 rebinについてcount・variance・under/overflow・entries・保存名が厳密に保存されることを確認した。元のnative gate count・fit countは旧refined結果・今回のJSON・独立再計算TSVで全210区間が完全一致。0.50以上の200区間×2modelでは、共有fit field 15,200項目と主表示の選択が前回とexactに一致した（追加の`lowerRigiditySearchApplied`監査fieldだけを除く）。結果は [独立監査manifest](validation/from040_audit_final_manifest.json) に保存した。実行した位置探索によって本番のtrack選択数を変えたわけではない。

R=0.405・0.425でtriton・³Heの2ピークと元校正から移動した⁴He信号を含む人工例を使い、平均・幅、元gate純度、独立な数値積分、実保存軸でのfit高端clip、低剛性raw診断範囲を検証した。移動gateが保存範囲を超える例ではその純度を不採用とし、厳密R=0.400で元gate自体が不足する場合はfitを拒否した。さらに旧平均下端−3.2より低い平均−3.5の分離した人工ピークを追加し、新しい平均制約とraw探索窓で回収できることを確認した。全7例はPASS。結果は [closure provenance](validation/from040_toy_final_v4/provenance.json) と [ROOT実行log](validation_from040_toy_final_v4.log) に保存した。

最終C++ソースのSHA256は `ab25f839121320a72b7864f1eb2a64e2d25fb3c3185875d49745c8119185abe3`。raw診断窓の最終補正前にも一度全fitを実行したため、その成果物23ファイルは削除せず新フォルダの `revisions/pre_peakwindow/` へ保存し、移動前後のhash一致を確認した。最終結果は補正後source snapshotと設定を用いて再実行している。

追加図は保存済みfitを使った描画だけであり、数値fitの再実行ではない。成分和、nSigmaとraw dE/dxの座標変換、count・誤差の同値、PDFページ数、旧結果と新結果の全生成物hashを検証した。

この再fit時点では元ROOT、production PIDヘッダー／YAML、以前のfit成果物を保持した。その後の明示的な依頼により3p85の現行⁴He mainconfのPID参照を更新し、旧条件への切り戻し設定を保存する。共有PIDヘッダー、追加keyがない他のpreset、元ROOTとfit成果物は変更しない。現行の変更・検証状態は[校正・CF入力の実装報告](task3_helium_pid_cf_implementation_report.md)を参照。
