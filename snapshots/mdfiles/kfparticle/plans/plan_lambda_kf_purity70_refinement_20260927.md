# Plan: KFParticle Lambda high-purity refinement — 2026-09-27

## User objective and scope

前回の比較を踏まえ、ユーザーは信号数が約70%になってもpurity約98%を優先する方針を示した。ここでは「70%」を標準KF Imp5に対する同一質量窓の背景差引き信号 S の保持率として扱い、背景込み候補数 Nwin の保持率も併記する。標準cutは自動変更しない。

2026-09-27 EDTに元PicoDstの最初のファイルを再確認したが、xrdstar.rcf.bnl.gov:1095のbounded statはexit51、Connection errorとなった。以前のcacheも不在。従って、今回も[前回の同じ10,000入力イベントの保存済みKF候補](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_study_20260927/README.md)を追加選別する。新規PicoDst再構成・upstream track/Finder cutの変更ではない。

前回のtraining / validation両方を既に見ているため、今回は**全10,000eventを用いた探索的refinement**とする。既存splitは安定性の診断としてのみ再利用し、新しい独立holdoutとは呼ばない。母数は同じ10,000入力 / 9,744再構成event、Au+Au13.5 GeV beam FXT。

## Fixed definitions

| Item | Definition |
|---|---|
| Signal window | [1.110, 1.122) GeV/c2 |
| Default sidebands | [1.098, 1.106), [1.126, 1.134) GeV/c2 |
| Near sidebands | [1.100, 1.108), [1.124, 1.132) GeV/c2 |
| Far sidebands | [1.094, 1.102), [1.130, 1.138) GeV/c2 |
| Background | B = 0.75 x (Nleft + Nright) |
| Signal | S = Nwin - B |
| Purity | S / Nwin |
| Signal-retention denominator | KF Imp5 baseline S with the same sideband definition |
| Display range | 1.10–1.14 GeV/c2; unchanged1MeV/c2 native bins |
| Baseline | Original selected==1, with all inherited KF and Imp5 cuts |

mass constraint、質量窓最適化、pT/eta受容域cut、TOF、massError cutは追加しない。共分散等の元のKF品質条件は維持する。

## Predeclared bounded grid

| Parameter | Values |
|---|---|
| minDecayLengthSignificance | 7, 8, 8.5, 9, 9.5, 10, 10.5, 11, 11.5, 12, 13 |
| imp5MinDCAPion [cm] | 1.5, 1.75, 2, 2.25, 2.5, 2.75, 3 |
| nSigmaProton, absolute maximum | 2.75, 2.5, 2.25, 2 |
| minCosPointing | 0.999, 0.9995, 0.9998 |
| imp5MinDCAProton [cm] | 0.8, 0.85, 0.9 |
| maxDistanceToPv [cm] | 0.45, 0.4 |
| maxTopoChi2Ndf | 2.5, 2 |
| maxDaughterDistance [cm] | 0.45, 0.4, 0.35 |
| maxChi2Ndf | 5, 3 |

baseline、全single、全pair、minDecayLengthSignificanceを含む全tripleを列挙する。合計3,756点（上限6,000）。追加cutは最大3条件で、baseの選別を緩めない。前回の `minDecayLengthSignificance>=10 && imp5MinDCAPion>=2 cm` を固定referenceとして含める。

## Ranking fixed before refined results

1. full sampleの3種類のsidebandそれぞれでS/Sbaselineを計算し、その最小値が70/75/80%以上であることを要求する。mass窓候補数の下限もQA YAMLで固定する。
2. 少数sideband count kにはPoisson平均の片側90%上限 muUpper を使う（CDF(k;muUpper)=0.10）。k=0でもmuUpper=ln(10)>0となるため、背景0件を誤差0のpurity100%として優遇しない。
3. 各sidebandの `1 - 0.75 * muUpper / Nwin` を計算し、その最小値を最大化する。これは少数背景とsideband位置への感度を考慮した**順位付けproxy**。Nwinの不確かさ、event内相関、多重探索を含む厳密なpurity信頼下限ではない。
4. old scanの「sideband>=8」の打切りは新rankingに用いない。すべての少数カウントを上限付きで評価する。
5. robust70/75/80の代表と、同じ70%保持率条件下でdefaultの名目purityが最大の点を保存する。両者が同じなら重複させない。僅差で複雑なcutを過大評価しないよう、同scoreでは少ない追加条件・高いS・固定ID順を用いる。
6. 代表点を決めた後に、既存split別結果、paired event bootstrap、background-model cross-checkを評価する。モデルを見てcutを再最適化したことにはしない。

## Background and implementation checks

前回は単一Gaussian+pol2 fitの負の背景とsignal tailのモデル依存が問題となった。今回はsignal形状に依存しないsideband-onlyのintegrated-bin Poisson fitを追加し、非負Bernstein background degree0/1/2で比較する。degree1は傾き、degree2は曲率感度の診断。sidebandの左右非対称自体は問題ではなく、対称窓の単純sideband推定は線形背景でも正しい。

window Bを固定してshapeをprofileした上限を診断として出し、DeltaNLL=0.821187207574908（名目片側90%の漸近閾値）を使う。sparse count / boundaryでの厳密なcoverageを保証しない。positive backgroundであっても未観測のpeaking background、誤同定反射、signal leakageがない証明にはならない。

独立ROOT計算で元treeに選定cutを再適用し、各条件の全200bins+2flowsを照合する。全10,000eventを含むpaired event bootstrapは固定cutの統計誤差のみ。既存コード・標準YAML・元ROOT・前回結果は不変とし、新規script/config/出力領域へ保存する。

Outputs: `rootfile/auau13p5_anaLambda_KFParticle/purity70_refinement_20260927/` and `share/figure/auau13p5_Lambda_KFParticle_purity70_20260927/`. ROOT検証はsingularity-local-build-run手順のROOT5/SL24y環境。farm job、production切替、Git書込みは行わない。
