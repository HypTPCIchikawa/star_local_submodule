# Task 3: ⁴HeのTPC純度推定 — 完了報告

2026-10-13の最新状況：0.40開始のfitから³He・⁴Heの中心／幅を3p85の現行⁴He解析へ導入し、³Heの背景CF入力を追加した。[実装・検証報告](task3_helium_pid_cf_implementation_report.md)と[更新したPurity対p/q（4ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/calibrated_selection/purity_calibrated_windows.pdf)を参照。以下は比較のため残す旧設定の結果であり、新しい選別の純度として流用しない。

この報告の数値は**旧設定の結果**です。低剛性のfit不採用はピーク不存在の判定ではありません。[低剛性再探索の計画](task3_lowrigidity_refinement_plan.md)に沿って追補を実施しています。dE/dxのFloat16保存と0.1 keV/cm bin幅による振動も[独立検証](validation/float16_dedx_alias_002/float16_alias.json)で確認され、旧fitの大きいχ²を応答形状の不一致だけに帰すことはできません。旧PDF・ROOT・数値表は比較用に保持します。

既存ROOTから剛性を0.01 GeV/cごとに分割したdE/dx射影をfitし、元の⁴He校正nSigma±2窓内の純度を推定する処理とPDFを作成した。13.5 GeVでは210区間中168区間が今回のfit品質基準を満たす。一方、3.85 GeVではGaussian＋単純な背景関数では高統計の形状を統計誤差の精度で説明できず、純度は**背景モデルに依存する参考推定**である。全⁴Heサンプルの確定純度や、TOF追加後の純度として使用しない。

## 出力ファイル

| 対象 | 剛性対純度（1ページ） | 背景モデル比較・差（2ページ） | 全射影fit（53ページ） | 数値表 |
|---|---|---|---|---|
| 3.85 GeV FXT、2021 P24iy | [purity PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_purity.pdf) | [background PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_background_models.pdf) | [fit PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_fits.pdf) | [CSV](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/purity.csv) |
| 13.5 GeV FXT、2020 P24iy | [purity PDF](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_purity.pdf) | [background PDF](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_background_models.pdf) | [fit PDF](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/he4_fits.pdf) | [CSV](../../../../rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/purity.csv) |

それぞれの `pdf/he4_purity_20261013/` にREADME、`purity.json`、`he4_fits.root`、設定・ソースsnapshot、実行ログ、SHA256検査記録をまとめた。途中のコンパイル・描画試行は内部 `revisions/` に保持し、移動前後のhashを確認した。元ROOTの変更、再hadd、新規job投入は行っていない。

## 対象と実際の選択条件

3.85 GeVは2021 P24iyの26,956/26,994 subjob部分統合、13.5 GeVは2020 P24iyの542/545 subjob部分統合を使用した。3.85 GeVの2019年データやpilotは混ぜていない。両エネルギーともfixed-targetのビームエネルギー表記である。

fit入力は `true/hDedxP`。`FemtoLambdaLegacy::CollectNuclei` で核種PID・TOF・最終剛性上限より前に保存されるtrack分布で、`nHitsDedx >= 15`、`gPt >= 0.1 GeV/c`、dE/dx>0を満たす。正負電荷と全採用中央度を合算し、Λ・原子核のペア数で重み付けしていない。

既存の `true/hDedxP_4He` は既にTPC±2のgateで切断されているため、背景fitに使わなかった。保存時の原子核設定は `maxNSigmaNuclear: 2.0`、`maxPOverQ: 2.5`、`m2_selection: false`、best-species選択あり。Todoに記載された上限2.0とは異なり、既存ROOTの実際の上限は2.5だった。今回のfit範囲は0.40–2.50 GeV/cである。

## fitと純度の定義

保存TH2は剛性幅0.01 GeV/c、dE/dx幅0.1 keV/cm。各剛性binを射影し、その中央剛性Rで `NuclearIdDeDxVsMom.h` の⁴He校正平均μcal(R)・幅σcal(R)を直接使用した。fit図の横軸は元のdE/dx射影を線形変換した校正nSigma `(dE/dx−μcal)/σcal`。元のkeV/cm射影も `he4_fits.root` の `raw_000` 等に保存している。

fitモデルは次の2種類で比較した。

- ⁴He Gaussian ＋ 非負のlinear背景 ＋ fit範囲内の近傍核種Gaussian。
- ⁴He Gaussian ＋ exponential背景 ＋ fit範囲内の近傍核種Gaussian。

³Heのピークは⁴Heから約−5～−3校正nSigmaにあり、そのtailも⁴He窓内の背景に含める。低剛性で範囲に入るtritonも背景Gaussianとして扱う。平滑背景だけに全ての核種ピークを吸収させていない。

各保存dE/dx binの期待数はGaussianのerf積分、linear/exponentialの解析積分で計算し、Poisson likelihoodでfitする。収束・covariance・物理的パラメータ・境界を確認したモデルのうち、Poisson AICが小さい方を主表示に選ぶ。fit範囲・幅・平均の制限・振幅上限等は全てYAMLに記録する。

純度は **校正nSigmaの−2～+2内**の `S/(S+B)`。fit Gaussianの幅を使って±2σの窓を移動した値ではない。⁴He Gaussianの積分がS、平滑背景＋近傍核種Gaussianの積分がBである。

誤差棒はfit covarianceを伝播した誤差である。背景モデル間差は別途表示する。この2関数の差は完全な系統誤差ではなく、片方のモデルが明らかにpoor fitの区間もあるため、無条件に「対等な良いモデルの誤差」と解釈しない。

## 結果とfit品質

| 対象 | 今回の品質基準を通過 | 収束するが品質基準外の参考fit | 純度を採用できない区間 | 品質通過区間が覆うTPC窓内候補数の割合 |
|---|---:|---:|---:|---:|
| 3.85 GeV | 0/210 | 175/210 | 35/210 | 0% |
| 13.5 GeV | 168/210 | 14/210 | 28/210 | 34.94% |

今回の品質基準はfit成功、covariance status≥2、信号等が境界にない、非負で有限な期待数に加え、Pearson χ²/ndf≤5、純度の相対covariance誤差≤20%。これは今回設定したモデル検査基準であり、純度の物理的正確さを保証するものではない。

13.5 GeVで168/210区間が通っても、候補数でのcoverageは34.94%である。低剛性の高収量区間にはfit不採用・形状不一致があるため、全体の純度を168区間の平均で代表させることはできない。数値的に有効なfitが主表示に選ばれた数は、品質基準外の参考14区間も含めてexponential 158、linear 24（計182区間）。3.85 GeVではexponential 135、linear 40で、いずれも参考値である。

Pがほぼ1となる区間には、背景振幅が非負制約の下端付近へ寄る場合がある。この境界ではcovarianceだけによる対称誤差が背景の片側上限を過小評価し得るため、非常に小さい誤差棒を「背景がゼロと証明された」と解釈しない。今回はprofile-likelihoodによる片側上限までは計算していない。

代表区間の両モデルの結果は以下。差の単位ppはpercentage pointである。ここに記す数値は剛性binごとのTPC窓の推定で、最終⁴He候補全体の純度ではない。

| 対象 | R [GeV/c] | linearのP | exponentialのP | 主表示 | 主表示のcovariance誤差 [pp] | モデル差 [pp] | 主表示χ²/ndf |
|---|---|---:|---:|---|---:|---:|---:|
| 3.85 GeV | 1.20–1.21 | 80.82% | 90.89% | exponential・参考 | 0.077 | 10.07 | 85.23 |
| 3.85 GeV | 1.60–1.61 | 65.41% | 84.43% | exponential・参考 | 0.171 | 19.01 | 76.72 |
| 3.85 GeV | 2.00–2.01 | 68.49% | 87.65% | exponential・参考 | 0.168 | 19.16 | 48.11 |
| 3.85 GeV | 2.40–2.41 | 62.82% | 64.73% | linear・参考 | 0.184 | 1.91 | 52.57 |
| 13.5 GeV | 1.20–1.21 | 86.00% | 92.86% | exponential | 1.267 | 6.86 | 1.01 |
| 13.5 GeV | 1.60–1.61 | 65.88% | 83.30% | exponential | 5.783 | 17.41 | 0.94 |
| 13.5 GeV | 2.00–2.01 | 71.77% | 83.18% | exponential | 6.186 | 11.41 | 0.98 |
| 13.5 GeV | 2.40–2.41 | 51.27% | 54.40% | linear | 6.522 | 3.13 | 1.00 |

3.85 GeVの代表4区間で、TPC窓内のbinごとの絶対残差合計／観測数は順に4.76%、3.91%、3.28%、2.63%。図上では主ピークに近くても、非常に高い統計ではこの程度の形状差が大きなχ²となる。covariance誤差が0.1 pp程度という結果だけを精密な純度測定と解釈できない。背景モデル差は10～19 ppになる区間もある。

背景モデルをAICで切り替えるため、主表示の曲線にはモデル切替由来の段差も生じる。物理的な純度の不連続と断言せず、2つのモデルを並べたPDFも確認する。

## この結果から言えること・言えないこと

剛性の大きい領域では³Heとの分離が小さくなり、TPC±2だけで十分に高い純度を維持できるとは限らない。例えば2.40–2.41 GeV/cでは、2関数の推定は3.85 GeVで約63～65%、13.5 GeVで約51～54%となる。Task 2の高剛性TOF選択を検討する動機と整合するが、TOF追加の改善量を今回のTH2から測ったわけではない。

以下は既存ROOTから確定できない。

- 正電荷⁴Heだけの純度（保存ヒストに電荷軸がない）。
- 4核種間のbest-species判定後の厳密な純度。
- R≥2でTOFを併用した新しい選択後の純度（TPCとTOFのjoint情報が保存されていない）。
- 中央度・ペア・k*別の⁴He純度。

Task 2で追加した `true/hDedxP_4He_TOFControl` / `true/hNSigmaHe4VsP_TOFControl` は、高剛性ではTOF条件後・TPC gate前の母集団を保存するため、将来の再解析でこのfitを同じ考え方で適用できる。今回の旧ROOTのTPC純度をCFの⁴He purity係数へ自動代入していない。

剛性bin中央の校正gateも近似である。bin内のtrackの連続な剛性位置は復元できず、0.40–0.41 GeV/cの校正平均は約0.417 σcalだけ変化する。gate端の変化をCSV/JSONに残した。最も細かい保存binで計算したが、特に低剛性でtrack単位の厳密なgate再適用と一致するとは限らない。

## 検証と再現

- ROOT 5.34/38・SL24yのSingularity環境で両エネルギーのfitを実行し、原本ROOTのfingerprint不変を確認した。
- 完成PDFは各サンプルで主表示1ページ、背景比較2ページ、全射影fit53ページ。ソース・設定・最終出力のhashを検査した。
- fitの状態・covariance・Gaussian幅・全成分の積分・gate内残差を `purity.json` に保存した。
- 人工Gaussian＋linear/exponential＋³He混合分布のclosure testで、既知純度約0.990413、0.995375をそれぞれ回収し、fit status 0・covariance 3を確認した。解析的積分と独立したTF1数値積分も1e−8以内で一致した。このtestは積分・fit処理の検証であり、実データの背景モデルの正しさを保証しない。

最終ソースを `.C++` で強制再compileした [closure結果](validation/task3_toy_final/toy_closure.json)、[ソースSHA256記録](validation/task3_toy_final/provenance.json)、[実行ログ](validation_he4_purity_toy_final.log) も保存した。最終closureの純度回収の絶対差はlinear約7.03×10⁻⁹、exponential約2.71×10⁻⁹である。

再現コマンド:

```bash
python3 tools/he4_purity/run_purity.py --sample 3p85 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --output-dir /tmp/he4_purity_3p85_reproduction

python3 tools/he4_purity/run_purity.py --sample 13p5 \
  --mainconf config/mainconf/main_auau13p5_anaFemtoLambda_4He_KFParticle.yaml \
  --output-dir /tmp/he4_purity_13p5_reproduction
```

出力先が非空の場合は上書きしない。設定は [purity_settings.yaml](../../../../tools/he4_purity/purity_settings.yaml)、手順は [tools README](../../../../tools/he4_purity/README.md)、実装前計画は [task3_he4_purity_plan.md](task3_he4_purity_plan.md) に記録している。
