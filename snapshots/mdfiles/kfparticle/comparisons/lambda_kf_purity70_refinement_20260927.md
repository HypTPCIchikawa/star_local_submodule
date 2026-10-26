# KFParticle Lambda high-purity refinement — 2026-09-27

Study started 2026-09-27 EDT; report finalized 2026-09-28 EDT. Filenames retain the study start date.

## Conclusion

「標準KF Imp5に対し信号数が約70%になってもよい」という方針で、同じ先頭10,000入力イベントの保存済みKF候補に3,756通りの追加cutを適用した。今回のgrid内で**default sidebandによる名目purityが最大**だった条件は次の3条件。

```yaml
# Additional/tighter cuts on the original KF Imp5-selected candidates
nSigmaProton: 2.5                 # abs(Pico nSigmaProton) <= 2.5
maxDistanceToPv: 0.45             # cm; final KF parent-to-PV distance
minDecayLengthSignificance: 10.0  # final PV-constrained decay length / error
```

この `nominal70` はNwin=892、背景推定B=15、信号推定S=877、purity **98.32%**、信号保持率 **71.49%**。ただし、背景曲率を許す別モデルではpurity推定が約95.86–98.82%となり、98%を保証できたわけではない。

background領域の選び方と少数countの不確かさを考慮した事前の順位付けでは、別の `robust70` が最良だった。こちらを**安定性を重視する実用候補**として併記する。

```yaml
# Alternative full-data exploratory candidate; not a production default switch
imp5MinDCAPion: 1.75              # cm; original pion track gDCA to PV
maxChi2Ndf: 5.0                  # final Lambda fit chi2/NDF, not Finder threshold
minDecayLengthSignificance: 10.5
```

`robust70` はNwin=905、B=16.5、S=888.5、purity **98.18%**、信号保持率 **72.43%**。モデルを変えたpurity範囲は約96.76–98.43%。これも98%確定ではないが、今回の名目最良条件より背景曲率への変動が小さかった。

前回の高purity条件は98.08%・信号保持率74.91%なので、今回の改善は小さい。約70%まで減らすことによって劇的にpurityが伸びた、とは結論しない。**最終的なproduction cutへの切替は行っていない。**

## What was actually executed

- 元PicoDstへのbounded再接続はexit51 / Connection errorで失敗し、以前のcacheも不在。[入力確認](../../../rootfile/auau13p5_anaLambda_KFParticle/purity70_refinement_20260927/provenance/input_status.md)を参照。
- 使用したのは以前の10,000入力 / 9,744再構成eventの[保存済みcandidate export](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_study_20260927/baseline_export/metadata.txt)。**新規PicoDst再構成ではない。**
- Au+Au、13.5 GeV beam-energy FXT。sqrt(sNN)は約5.2 GeVで、collider13.5 GeVではない。
- 元の `selected==1` に追加cutをANDした。元KF/Imp5で失った候補の回復、daughter primary chi-square、hits、Finder内部cut、TOFの再評価は行っていない。
- 前回のtraining / validation両方の結果を既に見ているため、今回は全10,000eventを用いた**探索的refinement**。旧splitは安定性の診断のみで、新しい独立holdoutではない。
- 質量窓、event選別、mass constraintなし、binning、pT/eta受容域を維持した。背景の少ない質量窓へ狭めてpurityを上げたものではない。

## Metrics and predeclared ranking

| Quantity | Definition |
|---|---|
| Signal window | [1.110, 1.122) GeV/c2 |
| Default sidebands | [1.098, 1.106), [1.126, 1.134) GeV/c2 |
| Near sidebands | [1.100, 1.108), [1.124, 1.132) GeV/c2 |
| Far sidebands | [1.094, 1.102), [1.130, 1.138) GeV/c2 |
| Nwin | All candidates in the fixed signal window |
| Background B | 0.75 x (Nleft + Nright) |
| Signal S | Nwin - B |
| Purity | S / Nwin |
| Signal retention | S / Sbaseline with the same sideband definition |
| Candidate retention | Nwin / 1370; different from signal retention |

9軸のsingle / pairと、崩壊長有意度を含むtripleを走査した。3,756点の内訳はbaseline1、single37、pair572、triple3,146。閾値は[QA YAML](../../../config/qa/lambda_kf_purity70_refinement_20260927.yaml)、事前規則は[plan](../plans/plan_lambda_kf_purity70_refinement_20260927.md)に保存した。

70/75/80%の保持率条件は、default/near/farそれぞれのS/Sbaselineの**最小値**に課した。70%で1,536点、75%で866点、80%で408点が適格だった。

robust順位は各sidebandの `1 - 0.75 * PoissonUpper90(Nsideband) / Nwin` の最小値を最大化する。背景0件でも片側90%上限はln(10)>0なので、誤差0のpurity100%と扱わない。前回のsideband>=8という足切りは外した。同scoreでは追加条件数が少ない、Sが大きい、固定point IDの順で選ぶ。

`robust70` のrobustという名称はこのsideband / Poisson上限の順位付けに限定した意味で、あらゆる背景曲率に対して安定という保証ではない。

このscoreは**順位付け用proxy**で、Nwinの不確かさ、event内相関、多重探索、背景形状を含む厳密なpurity信頼下限ではない。`nominal70` は同じ保持率適格条件の中でdefaultの名目purityを最大化した別の選択。背景fitは順位付け後のcross-checkであり、fitを見て選び直していない。

## Complete cut comparison

各列は標準KF Imp5候補からの選別。`Disabled`は追加のfinal cutなしを意味し、元のFinder/Topo品質条件まで無いという意味ではない。

| Parameter | KF Imp5 baseline | old_aggressive | robust70 | nominal70 | robust75 | robust80 |
|---|---:|---:|---:|---:|---:|---:|
| nSigmaProton, absolute maximum | 3.0 | 3.0 | 3.0 | 2.5 | 3.0 | 3.0 |
| nSigmaPion, absolute maximum | 3.0 | 3.0 | 3.0 | 3.0 | 3.0 | 3.0 |
| imp5MinDCAProton [cm] | 0.7 | 0.7 | 0.7 | 0.7 | 0.7 | 0.7 |
| imp5MinDCAPion [cm] | 1.0 | 2.0 | 1.75 | 1.0 | 2.5 | 1.5 |
| maxDaughterDistance [cm] | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |
| maxDistanceToPv [cm] | 0.5 | 0.5 | 0.5 | 0.45 | 0.5 | 0.5 |
| minCosPointing | 0.998 | 0.998 | 0.998 | 0.998 | 0.998 | 0.998 |
| minNHitsFit | 15 | 15 | 15 | 15 | 15 | 15 |
| minNHitsRatio | 0.52 | 0.52 | 0.52 | 0.52 | 0.52 | 0.52 |
| imp5MaxPathLength [cm] | 100 | 100 | 100 | 100 | 100 | 100 |
| Final maxChi2Ndf | Disabled | Disabled | 5.0 | Disabled | 5.0 | 5.0 |
| Final maxTopoChi2Ndf | Disabled | Disabled | Disabled | Disabled | Disabled | Disabled |
| Final minDecayLength [cm] | Disabled | Disabled | Disabled | Disabled | Disabled | Disabled |
| Final minDecayLengthSignificance | Disabled | 10.0 | 10.5 | 10.0 | 9.5 | 8.5 |
| Final minVertexLineSignificance | Disabled | Disabled | Disabled | Disabled | Disabled | Disabled |
| Final maxMassError | Disabled | Disabled | Disabled | Disabled | Disabled | Disabled |

nHitsRatioは元のlambda_imp5規則どおりnHitsMax>0の場合に適用。helix pathは両daughterの絶対値に対する上限。TOF PIDは使用せず、HFT trackであることも要求しない（HFT trackを排除する意味ではない）。

全列でinterfaceChiPrimaryCut=3、finderChiPrimary2D=3、finderChi2Ndf2D=10、finderLCut=1 cm、finderLdL2D=3、finderMaxDaughterDistance=1.5 cm、topoChi2NdfCut=3およびcovariance品質条件を継承する。`maxChi2Ndf` はfinal Lambdaの値、`minDecayLengthSignificance` はPV constraint後のparentの崩壊長/誤差であり、upstreamのdaughter primary chi-squareやFinder L/σLとは別の量。

## Results on all 10,000 inputs

| Selection | Point ID | Nwin | B | S | Purity [%] | Signal retention [%] | Minimum retention across sidebands [%] | Candidate retention [%] |
|---|---:|---:|---:|---:|---:|---:|---:|---:|
| KF Imp5 baseline | 0 | 1370 | 143.25 | 1226.75 | 89.54 | 100.00 | 100.00 | 100.00 |
| old_aggressive | 75 | 937 | 18.00 | 919.00 | 98.08 | 74.91 | 74.20 | 68.39 |
| robust70 | 2005 | 905 | 16.50 | 888.50 | 98.18 | 72.43 | 71.86 | 66.06 |
| nominal70 | 2379 | 892 | 15.00 | 877.00 | 98.32 | 71.49 | 70.88 | 65.11 |
| robust75 | 1983 | 949 | 18.00 | 931.00 | 98.10 | 75.89 | 75.11 | 69.27 |
| robust80 | 1947 | 1013 | 24.00 | 989.00 | 97.63 | 80.62 | 80.09 | 73.94 |

「70%」は信号保持率であり、nominal70の背景込み候補数の保持率は65.11%。解析した入力event数を7,000に減らしたのではない。

nominal70は標準baseline比で背景推定を143.25から15へ約90%減らす一方、信号を約28.5%失う。old_aggressiveとの比較では、さらに信号42（919→877）を失い、背景推定は3（18→15）減る。前回に比べた改善の規模はこの程度である。

## Statistical uncertainty and comparison with the previous cut

全10,000event（候補0件を含む）への共通Poisson(1) weightを使った400回のpaired event bootstrap。選んだcutを固定した場合の統計変動のみで、探索による偏りや背景モデル系統誤差は含まない。

| Selection | Purity +/- bootstrap SD [%] | Signal-retention 16–84% interval [%] | Purity gain vs old_aggressive [percentage points] | Gain 16–84% interval [percentage points] |
|---|---:|---:|---:|---:|
| old_aggressive | 98.079 +/- 0.386 | 73.44 to 76.45 | Reference | Reference |
| robust70 | 98.177 +/- 0.385 | 70.90 to 74.07 | +0.098 | -0.013 to +0.215 |
| nominal70 | 98.318 +/- 0.373 | 69.98 to 73.09 | +0.239 | +0.029 to +0.391 |
| robust75 | 98.103 +/- 0.382 | 74.49 to 77.35 | +0.024 | -0.078 to +0.122 |

新旧のpurity差は相関した候補から計算しており、別々の誤差を単純に二乗和したものではない。nominal70の差のbootstrap SDは0.177 percentage points。多数のcutを同じデータで選んだ後の結果なので、前回から有意に改善したという強い主張はしない。70%の信号保持も点推定の条件であり、不確かさを含めた保証ではない。

旧splitの診断ではnominal70はpurity98.57% / 98.09%、信号保持率67.53% / 75.46%。robust70はpurity98.62% / 97.76%、保持率69.98% / 74.89%。全体での保持率条件を部分標本ごとに満たす保証はなく、これらを独立validationとは扱わない。

## Background-model dependence: the main remaining limitation

前回の単一Gaussian+pol2 fitは負の背景などで不適格だった。今回はsignalピークをfitせず、sidebandだけへ非負Bernstein background degree0/1/2をintegrated-bin Poisson likelihoodでfitした。degree0は一定、degree1は線形、degree2は曲率を許す。54通りすべてで数値解とprofile上限計算が成功した。

| Selection | Default constant/linear purity [%] | Near constant/linear purity [%] | Far constant/linear purity [%] | Default curved purity [%] | Near curved purity [%] | Far curved purity [%] |
|---|---:|---:|---:|---:|---:|---:|
| old_aggressive | 98.08 | 97.68 | 98.16 | 97.63 | 96.54 | 97.95 |
| robust70 | 98.18 | 97.93 | 98.43 | 97.74 | 96.76 | 97.98 |
| nominal70 | 98.32 | 97.90 | 98.49 | 96.94 | 95.86 | 98.82 |
| robust75 | 98.10 | 97.71 | 98.10 | 97.66 | 96.58 | 98.06 |

対称窓では一定・線形背景の積分推定は一致する。左右のsideband countが異なるだけで単純sideband法が破綻するわけではない。今回の大きな不確かさは、ピーク下へ外挿する背景の曲率である。

nominal70のdefault sidebandは合計20候補のみ。単純推定ではB=15だが、同じdefault領域のcurved fitはB=27.27となる。near領域のcurved fitではB=36.96、purity95.86%。nominal70が他の背景仮定でも必ず最大purityになるとは言えない。

profileではwindow内Bを固定してshapeを再最適化した。purityへの換算では観測Nwinを固定しているため、Nwinも含めた同時のpurity信頼区間ではない。DeltaNLL=0.821187207574908を名目片側90%の漸近閾値として使ったが、少数count・境界・cut選択を含めた厳密なcoverageではない。nominal70のdefault curved modelのprofile上限はB=42.75（purityに換算95.21%）、near curved modelはB=51.03（94.28%）。robust70のnear curved modelはB=44.76（95.05%）。これらも「真のpurityの保証下限」とは呼ばない。

**従って現時点では、約98%というsideband点推定を持つ条件は見つかったが、真のpurityが98%以上と確定したとは言えない。** positive background modelは負の多項式を避けるだけで、peaking reflection、signal tailのsideband漏れ、cutによる背景形状変化を除外する証拠ではない。fitの妥当性とモデルの正しさも別である。

## Practical choice and next validation

名目値を最大にするならnominal70、背景領域依存と少数countへの配慮を優先するならrobust70が今回の候補。追加の数%の収量を失うメリットが小さいと判断する場合は、98.10%・保持率75.89%のrobust75も合理的である。ユーザーが70%程度を許容しても、無理に70%ぎりぎりまで落とす必要はない。

productionへ反映する前に、同じ元PicoDstでの再構成と、別event/runでの確認が必要。特にnSigma/gDCAは今回post-KF候補に適用しており、upstream YAMLを変えた再処理と完全同一だとは未検証。元PicoDstが回復すれば、以前の[analysis note/qhu1比較](kfparticle_lambda_qhu1_cut_comparison_20260908.md)にあるdaughter primary chi-squareやhits/Finder条件も調べられるが、今回は未実施。

signalには実Lambdaのfeeddownも含むため、これはprimary-Lambda purityではない。DCA・崩壊長のcutはpT/rapidity/lifetime受容域にも影響し得る。S保持率は検出器・再構成効率の代わりにはならない。

## Files and verification

- [Main mass overlay](../../../share/figure/auau13p5_Lambda_KFParticle_purity70_20260927/lambda_purity70_overlay_00.pdf)：baseline / old_aggressive / robust70 / nominal70。1.10–1.14 GeV/c2、左counts / 右全[1.05,1.25)範囲Normalize。
- [Higher-retention overlay](../../../share/figure/auau13p5_Lambda_KFParticle_purity70_20260927/lambda_purity70_overlay_01.pdf)：baseline / old_aggressive / robust75 / robust80。
- [nominal70 background diagnostics](../../../share/figure/auau13p5_Lambda_KFParticle_purity70_20260927/lambda_purity70_background_05.pdf)、[robust70 background diagnostics](../../../share/figure/auau13p5_Lambda_KFParticle_purity70_20260927/lambda_purity70_background_02.pdf)。すべてPNG版も保存。
- [Background CSV](../../../share/figure/auau13p5_Lambda_KFParticle_purity70_20260927/lambda_purity70_background.csv)、[background audit](../../../share/figure/auau13p5_Lambda_KFParticle_purity70_20260927/lambda_purity70_audit.txt)。
- [Output index / reproduction](../../../rootfile/auau13p5_anaLambda_KFParticle/purity70_refinement_20260927/README.md)、[selection validation](../../../rootfile/auau13p5_anaLambda_KFParticle/purity70_refinement_20260927/provenance/selection_validation.md)、[independently selected ROOT](../../../rootfile/auau13p5_anaLambda_KFParticle/purity70_refinement_20260927/verified_selections.root)。

Python7tests、ROOT背景計算の人工データtests、全6条件の元ROOTからの200bins+2flows照合、33,804行の別計算によるmetric/ranking照合、10,000event IDの一致を確認した。ROOT検証はsingularity-local-build-run手順に従いSL24y / ROOT5.34/38で実施した。標準YAML、両Maker、`.current_mainconf`、元ROOT、前回studyのコード・結果が不変であることをSHA256で確認した。farm job、default切替、git commit/pushは行っていない。
