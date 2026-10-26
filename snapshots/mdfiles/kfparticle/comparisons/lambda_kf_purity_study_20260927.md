# KFParticle Lambda purity study — 2026-09-27

## Result and scope

先頭10,000入力イベントの保存済みKF Imp5候補に追加cutを適用し、1,020条件を比較した。統計を優先する場合、`minDecayLengthSignificance >= 5` が有望である。単独でsideband推定purityは89.54%から94.46%、背景差引き信号の保持率は94.95%となる。事前に定めたtraining側の選定規則で選ばれた組合せでは、これに `abs(nSigmaProton) <= 2.5` を加えるとpurity94.77%、信号保持率93.09%となった。

ユーザーの「結果の比較を見て決める」という方針に従い、production設定は変更していない。約98%のpurityを狙う条件では信号を約25%失うため、単一の最適設定とは結論しない。また、以下のpurityは背景モデルに依存する推定値であり、質量fitでは数ポイント低くなる。BG freeや真の再構成効率の測定とはしない。

重要な実行範囲：前回のPicoDst cacheが消失し、元XRootDへの接続も失敗したため、今回は**同じ10,000イベントから既に再構成した候補treeの再選別**を行った。PicoDstからの新規再構成ではない。元trackのhitsやdaughter primary chi-square等、treeに保存されていないupstream変数の変更は未実施。[入力・接続確認](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_study_20260927/provenance/input_status.md)に区別して記録した。

## Dataset and baseline

- Au+Au、13.5 GeV **beam-energy FXT**（sqrt(sNN) approximately 5.2 GeV）。collider sqrt(sNN)=13.5 GeVではない。
- 入力リストの先頭10,000イベント、再構成9,744イベント。保存済みordered run/event IDと候補の所属を確認した。
- 元ROOT：[default_regression.root](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_kf_specific_scan_20260914/default_regression.root)。raw treeは42,871候補、標準選別後は質量全範囲で3,135候補。
- `backend: finder_topo`、`selectionProfile: lambda_imp5`、Lambdaのみ。mass constraintなし。
- 今回のbaselineは**geometry cutsを含む標準KF Imp5**。以前のgeometry-off studyではない。すべての走査点は元の `selected == 1` に追加cutをANDしており、upstreamで落とした候補を復活させない。
- baselineの元選別判定と200 bins + underflow/overflowを完全再現してから走査した。

## Fixed estimators and selection protocol

質量窓は `[1.110, 1.122)` GeV/c²で固定し、表示範囲だけを1.10–1.14 GeV/c²とした。元の1 MeV/c² binを維持する。

| Quantity | Definition |
|---|---|
| Nwin | Candidate count in [1.110, 1.122) GeV/c2 |
| Default sidebands | [1.098, 1.106) and [1.126, 1.134) GeV/c2 |
| Background B | 0.75 x (Nleft + Nright) |
| Signal S | Nwin - B |
| Purity | S / Nwin = S / (S + B) |
| Signal retention | S / Sbaseline in the same split and sideband definition |
| Candidate retention | Nwin / Nwin_baseline; not signal retention |
| Near sidebands | [1.100, 1.108) and [1.124, 1.132) GeV/c2 |
| Far sidebands | [1.094, 1.102) and [1.130, 1.138) GeV/c2 |

runId/eventIdのhash（seed20260927）でevent単位にtraining4,958 / validation5,042へ分割した。同じeventの候補を両側に分けない。全12軸のsingle、すべてのpair、事前指定した4種類のtriple familyを調べ、合計1,020点とした。[走査YAML](../../../config/qa/lambda_kf_purity_study_20260927.yaml)が閾値の完全な一覧である。

代表4条件は**training / default sidebandのみ**で、信号保持率95/90/80/70%以上のそれぞれについて `purity - statistical_error` が最大の点を選んだ。Nwin>=50、sideband合計>=8も要求し、背景0件だけで最良とする選択を避けた。これは探索の順位付けであり、多重探索補正済みの信頼下限ではない。validation・全体・別sidebandを見て代表点を選び直していない。以下の保持率ラベルはtrainingで課した下限であって、全体・validationでの保証ではない。

## Selected cut comparison

各列は元のKF Imp5選別を維持した上での追加・強化条件。全列でTOF PIDは使用せず、HFT trackであることも要求しない（HFT trackの排除ではない）。nHitsFit>=15、nHitsMax>0の場合のnHitsRatio>=0.52、両daughterの元helix pathの絶対値<=100 cmを維持する。

| Parameter | KF Imp5 baseline | choice_01 | choice_02 | choice_03 | choice_04 |
|---|---:|---:|---:|---:|---:|
| Training signal-retention floor | Reference | 95% | 90% | 80% | 70% |
| nSigmaProton, absolute maximum | 3.0 | 2.5 | 2.5 | 3.0 | 3.0 |
| nSigmaPion, absolute maximum | 3.0 | 3.0 | 3.0 | 3.0 | 3.0 |
| imp5MinDCAProton [cm] | 0.7 | 0.7 | 0.7 | 0.7 | 0.7 |
| imp5MinDCAPion [cm] | 1.0 | 1.5 | 1.0 | 3.0 | 2.0 |
| maxDaughterDistance [cm] | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |
| maxDistanceToPv [cm] | 0.5 | 0.5 | 0.5 | 0.5 | 0.5 |
| minCosPointing | 0.998 | 0.998 | 0.998 | 0.998 | 0.998 |
| minNHitsFit | 15 | 15 | 15 | 15 | 15 |
| minNHitsRatio | 0.52 | 0.52 | 0.52 | 0.52 | 0.52 |
| imp5MaxPathLength [cm] | 100 | 100 | 100 | 100 | 100 |
| Final maxChi2Ndf | Disabled | Disabled | Disabled | Disabled | Disabled |
| Final maxTopoChi2Ndf | Disabled | Disabled | Disabled | Disabled | Disabled |
| Final minDecayLength [cm] | Disabled | Disabled | Disabled | Disabled | Disabled |
| Final minDecayLengthSignificance | Disabled | Disabled | 5.0 | 5.0 | 10.0 |
| Final minVertexLineSignificance | Disabled | Disabled | Disabled | Disabled | Disabled |
| Final maxMassError | Disabled | Disabled | Disabled | Disabled | Disabled |

`Disabled`は**追加のfinal cutが無い**ことを表す。KF upstream cutまで無いわけではない。全列でinterfaceChiPrimaryCut=3、finderChiPrimary2D=3、finderChi2Ndf2D=10、finderLCut=1 cm、finderLdL2D=3、finderMaxDaughterDistance=1.5 cm、topoChi2NdfCut=3とcovariance品質条件を継承する。完全なeffective設定は[baseline metadata](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_study_20260927/baseline_export/metadata.txt)に保存した。

`minDecayLengthSignificance`はPV constraintを付けたparentのdecay length / errorで、共分散を使用するKF-specific量。元track gDCAやnSigmaとは異なる。`minVertexLineSignificance`はraw particleからPVへのvertex-line量で、同じ名前のように扱わない。これらはFinderのSIMD L/σL cutとも同一の選別ではない。

## Results on all 10,000 inputs

以下は全体を用いた記述的な比較。purity誤差は400回のpaired Poisson event bootstrapの標準偏差（統計のみ）。候補0件のeventを含む全10,000eventを対象に、各条件とbaselineへ共通event weightを使った。

| Selection | Nwin | B | S | Purity [%] | Signal retention [%] | Candidate retention [%] | S/B |
|---|---:|---:|---:|---:|---:|---:|---:|
| KF Imp5 baseline | 1370 | 143.25 | 1226.75 | 89.54 +/- 0.83 | 100.00 | 100.00 | 8.56 |
| choice_01 | 1307 | 117.00 | 1190.00 | 91.05 +/- 0.77 | 97.00 | 95.40 | 10.17 |
| choice_02 | 1205 | 63.00 | 1142.00 | 94.77 +/- 0.58 | 93.09 | 87.96 | 18.13 |
| choice_03 | 1088 | 49.50 | 1038.50 | 95.45 +/- 0.57 | 84.66 | 79.42 | 20.98 |
| choice_04 | 937 | 18.00 | 919.00 | 98.08 +/- 0.39 | 74.91 | 68.39 | 51.06 |

choice_02は背景推定を143.25から63.00へ約56%減らし、信号減少は約7%。choice_03はchoice_02に比べpurityの上昇が約0.68 percentage pointsにとどまり、信号保持率はさらに約8.44 points下がる。choice_04はpurity優先の別選択肢であり、統計を保つ目的とは明確なtrade-offがある。

## Held-out validation

| Selection | Training purity [%] | Training signal retention [%] | Validation Nwin | Validation S | Validation purity [%] | Validation signal retention [%] |
|---|---:|---:|---:|---:|---:|---:|
| KF Imp5 baseline | 91.31 | 100.00 | 697 | 612.25 | 87.84 | 100.00 |
| choice_01 | 92.66 | 97.11 | 663 | 593.25 | 89.48 | 96.90 |
| choice_02 | 95.19 | 91.70 | 613 | 578.50 | 94.37 | 94.49 |
| choice_03 | 95.64 | 80.31 | 572 | 545.00 | 95.28 | 89.02 |
| choice_04 | 98.49 | 71.81 | 489 | 477.75 | 97.70 | 78.03 |

| Selection | Validation purity gain [percentage points] | Paired-bootstrap 16–84% interval [percentage points] | Validation signal retention, 16–84% interval [%] |
|---|---:|---:|---:|
| choice_01 | +1.64 | +1.19 to +2.13 | 95.86 to 97.99 |
| choice_02 | +6.53 | +5.65 to +7.44 | 92.58 to 96.24 |
| choice_03 | +7.44 | +6.47 to +8.49 | 87.16 to 90.94 |
| choice_04 | +9.86 | +8.74 to +10.95 | 75.90 to 80.25 |

同じ10,000event内のholdoutなので、別run/別データによる確認ではない。bootstrapも固定したcutの統計変動のみで、cut選択手順そのものの不確かさや背景モデル系統誤差を含まない。

## Single-cut responses

主に何が効いているかを見るため、同じ事前grid内の単独cutを示す。この表は探索結果の説明であり、代表4条件の再選定や追加の未開封validationとはしない。

| Added cut | Scan point | Full Nwin | Full S | Full purity [%] | Full signal retention [%] | Validation purity [%] | Validation signal retention [%] |
|---|---:|---:|---:|---:|---:|---:|---:|
| abs(nSigmaProton) <= 2.5 | 1 | 1338 | 1207.50 | 90.25 | 98.43 | 88.59 | 97.96 |
| imp5MinDCAPion >= 1.5 cm | 11 | 1338 | 1209.00 | 90.36 | 98.55 | 88.78 | 98.90 |
| imp5MinDCAPion >= 3.0 cm | 13 | 1118 | 1051.25 | 94.03 | 85.69 | 93.38 | 89.83 |
| maxChi2Ndf <= 3.0 | 24 | 1351 | 1212.25 | 89.73 | 98.82 | 88.19 | 98.82 |
| minDecayLengthSignificance >= 5.0 | 35 | 1233 | 1164.75 | 94.46 | 94.95 | 94.05 | 96.77 |

この範囲では崩壊長の有意度が最も有望な単独cut。point35とchoice_02の差はproton nSigmaを3.0から2.5へ強化するかどうかで、全体ではpurity約+0.31 points、信号保持率約-1.85 pointsである。統計を優先するなら、point35の単独cutも十分検討価値がある。point35のtraining信号保持率は93.12%なので「95% training floor」の代表には入らない。全体の94.95%と混同しない。

## Background-model checks

| Selection | Near-sideband purity [%] | Default-sideband purity [%] | Far-sideband purity [%] | Gaussian + pol2 purity [%], fit range 1.090–1.145 |
|---|---:|---:|---:|---:|
| KF Imp5 baseline | 89.49 | 89.54 | 90.47 | 86.30 |
| choice_01 | 90.99 | 91.05 | 91.74 | 87.78 |
| choice_02 | 94.65 | 94.77 | 94.90 | 91.73 |
| choice_03 | 95.24 | 95.45 | 95.93 | 92.63 |
| choice_04 | 97.68 | 98.08 | 98.16 | Invalid |

sideband位置を変えても改善の方向は変わらない。一方、Gaussian + quadratic backgroundの積分bin Poisson fitではpurityの絶対値が約3 points下がる。choice_02のfit由来信号保持率は全体93.46%、validation94.59%で、収量減少の程度はsideband法と整合するが、purity絶対値は確定しない。

fitは独立cross-checkでありcut選定には使っていない。狭いfit範囲でも全体のdeviance/NDFはbaseline81.88/49、choice_02は89.69/49（漸近p約0.0022、0.00035）となり、単一Gaussian + pol2のモデル適合には懸念がある。sparse binではこのp値自体も近似。`accepted=1`は数値的・物理的チェック通過であって、モデルの正しさの保証ではない。

広いfit範囲1.080–1.150では全体の全5条件が負の背景等で不適格。choice_04は狭い範囲も不適格。これらのS/B/purity結果はCSVで空欄にし、採用していない。無効fitのdeviance/p値も解釈しない。現段階で「98%を確定」「BG free」とは言えない。[fit CSV](../../../share/figure/auau13p5_Lambda_KFParticle_purity_20260927/lambda_purity_fit_metrics.csv)と[audit](../../../share/figure/auau13p5_Lambda_KFParticle_purity_20260927/lambda_purity_audit.txt)を参照。

## Reference-note connection and untested cuts

参考：[analysis note PDF](../../../ref_analysisnote/Analysnote_local_lambda_hyperons_polarization_BESII_v3_20260601.pdf)、[qhu1 comparison](kfparticle_lambda_qhu1_cut_comparison_20260908.md)。PDF p16 / Section 2.3はdaughter primary chi-square>10、hits>15等、p33 / Section 3はhits・primary chi-squareの変動を示す。Fig.9はcollider7.7 GeV、20–50%で、今回のFXT全対象とは違う。

原PicoDstが再取得できれば、daughter primary chi-square3→5/10/20、hits15→20、Finder L1→3/5 cm、Finder L/σL3→5/7/10が次の候補。ただし今回は**未実施**。parentに付け直したdaughter chi-square、単独trackのprimary chi-square、final decay-length significanceを同一視しない。TOFも元Imp5 treeでは未評価placeholderであるため、今回のtreeから効果を判定していない。

purityは質量ピークの純度であり、feeddown由来の実Lambdaを除くprimary-Lambda purityではない。DCA・崩壊長・pointingを変えるとpT/rapidity/lifetime受容域も変わり得る。今回の積分信号保持率はMC効率や微分効率の代わりではなく、最終physics解析では別途補正が必要。

## Figures and saved products

- [Purity versus signal retention](../../../share/figure/auau13p5_Lambda_KFParticle_purity_20260927/lambda_purity_retention.pdf)：全5条件、training / validationを分離。
- [Main mass overlay](../../../share/figure/auau13p5_Lambda_KFParticle_purity_20260927/lambda_purity_overlay_00.pdf)：baseline + choice_01/02/03、左counts / 右全[1.05,1.25)範囲でunit-area Normalize。
- [High-purity mass overlay](../../../share/figure/auau13p5_Lambda_KFParticle_purity_20260927/lambda_purity_overlay_01.pdf)：baseline + choice_04。5線を重ねず別図にした。
- [Full-sample fit cross-check](../../../share/figure/auau13p5_Lambda_KFParticle_purity_20260927/lambda_purity_fit_range00_full.pdf)。PNG版も同じ名前で保存。
- [Output index and reproduction commands](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_study_20260927/README.md)。全走査CSV、候補ID、bootstrap、source/config snapshot、SHA256を保存。
- [Independent selected ROOT histograms and candidate trees](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_study_20260927/verified_selections.root)。選定4条件 + baselineを元ROOTから再計算し全202 bins一致。
- [Plan](../plans/plan_lambda_kf_purity_study_20260927.md)。初期計画の入力取得確認中という状態は、現在は上記の取得失敗・保存tree再選別に確定。

## Validation and preservation

Python単体・人工データintegration testは9件成功。ROOT5で独立選別検証・描画・fitを実行し、元treeからの5条件全bin一致、histogram全体=train+validation、質量窓集計、Normalizeの分母を検査した。全走査9,180行とbootstrap135行は別計算でも検査済み。

`singularity-local-build-run`の手順に従い、STAR SL24y / ROOT5.34/38の既存環境を用いた。source/config/図の新規追加のみで、StLambdaMaker、StLambdaKFParticleMaker、標準Imp5 YAML、`.current_mainconf`、既存ROOT2本が不変であることをSHA256で確認した。farm jobの投入、default切替、git commit/pushは行っていない。
