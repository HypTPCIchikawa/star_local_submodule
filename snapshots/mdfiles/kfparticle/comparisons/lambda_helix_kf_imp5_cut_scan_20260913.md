# AuAu13p5：Imp5とKFParticleの候補数比較（2026-09-13／09-14）

**最新の2026-09-14 study:** Finder内部の幾何cutを外し、KF-specific cutを広く緩めても、今回の探索での最接近点は1,613候補（Helix Imp5: 1,775）だった。同数条件は未達。結果・全cut表・デフォルトへの復帰は末尾の [Geometry-off KF-specific scan](#geometry-off-kf-specific-scan-2026-09-14) を参照。以下の従来の4つのKF列は2026-09-13の距離cut scanであり、新studyとは区別する。

## 2026-09-13 distance-cut scan

**2026-09-14補足：KFParticle固有のcutを外した比較ではない。** 共通のTPC PID・track条件をImp5に合わせ、KF固有の共分散・Finder/Topo選別を維持したまま、追加の距離・pointing cutを調整した。調整後は対応する距離cutの閾値自体も旧Imp5とは異なる。

KF側の `maxDistanceToPv: 0.5 → 0.8 cm` と `maxDaughterDistance: 0.5 → 0.6 cm` を緩め、`minCosPointing: 0.998` を維持すると、Λ近傍の候補数は旧Imp5の約99.1%まで戻る。ただし、その条件ではsideband推定S/Bも旧Imp5とほぼ同じになり、明確な改善は確認できなかった。

同じ先頭10,000入力イベントの保存済み出力を使用。質量窓は **1.110 ≤ M(pπ⁻) < 1.122 GeV/c²**、1 MeV/c² bin。数は背景を含むΛ候補数であり、衝突イベント数・真のΛ収量・絶対効率ではない。

| Selection | maxDCAV0 [cm] | maxDaughterDCA [cm] | minCosPointing | Nwin | Nwin / Helix | S/B (sideband) |
|---|---:|---:|---:|---:|---:|---:|
| Helix Imp5 | 0.5 | 0.5 | 0.998 | 1775 | 100% | 5.20 ± 0.35 |
| KF baseline | 0.5 | 0.5 | 0.998 | 1370 | 77.2% | 8.56 ± 0.74 |
| KF: maxDCAV0 | 3.0 | 0.5 | 0.998 | 1653 | 93.1% | 5.02 ± 0.35 |
| KF: maxDCAV0 + maxDaughterDCA | 0.8 | 0.6 | 0.998 | 1759 | 99.1% | 5.22 ± 0.35 |
| KF: maxDCAV0 + maxDaughterDCA + minCosPointing | 0.75 | 0.6 | 0.997 | 1783 | 100.5% | 4.94 ± 0.33 |

誤差は非重複の質量窓・sidebandの独立Poisson近似による統計目安。背景形状、共有イベント・候補間相関、同じデータでcutを選んだ影響を含まない。従って差の有意性をこの誤差の単純な合成から主張しない。

## NotionとImp5設定の照合

公開後の [7月1日ページ](https://yellow-wolverine-c1e.notion.site/2026-07-01-3902284fda2880c4b97cd2fa9798edb3) の「各フェーズの具体的なCut条件パラメータ」を確認し、現行Imp5との一致を確認した。[9月8日ページ](https://yellow-wolverine-c1e.notion.site/2026-09-08-3d52284fda2880b4a510c34e49009007) も主要本文・比較表・設定整理・出力移動・ログを確認した。一部子ブロックや画像は未確認であり、全要素の完全取得とはしていない。

### Cut parameters

2026-09-14追記。添付表の10項目を同じ順序で、同値も省略せず全列に記入した。列名の `KF:` に続く英語パラメータ名は、`KF baseline` からの調整対象を示す。数値は設定された閾値で、`min` は下限、`max` は上限、`nSigma` は絶対値の上限。

| Parameter | Helix Imp5 | KF baseline | KF: maxDCAV0 | KF: maxDCAV0 + maxDaughterDCA | KF: maxDCAV0 + maxDaughterDCA + minCosPointing |
|---|---:|---:|---:|---:|---:|
| `nSigmaProton` | 3.0 | 3.0 | 3.0 | 3.0 | 3.0 |
| `nSigmaPion` | 3.0 | 3.0 | 3.0 | 3.0 | 3.0 |
| `minDCAProton` [cm] | 0.7 | 0.7 | 0.7 | 0.7 | 0.7 |
| `minDCAPion` [cm] | 1.0 | 1.0 | 1.0 | 1.0 | 1.0 |
| `maxDaughterDCA` [cm] | 0.5 | 0.5 | 0.5 | 0.6 | 0.6 |
| `maxDCAV0` [cm] | 0.5 | 0.5 | 3.0 | 0.8 | 0.75 |
| `minNHitsFit` | 15 | 15 | 15 | 15 | 15 |
| `minNHitsRatio` | 0.52 | 0.52 | 0.52 | 0.52 | 0.52 |
| `minCosPointing` | 0.998 | 0.998 | 0.998 | 0.998 | 0.997 |
| `maxPathLength` [cm] | 100.0 | 100.0 | 100.0 | 100.0 | 100.0 |

`paired_cut_overlay.png` は `Helix Imp5`、`KF baseline`、`KF: maxDCAV0 + maxDaughterDCA` の3曲線。`pv_only_overlay.png` の追加曲線は `KF: maxDCAV0`、`count_matched_overlay.png` の追加曲線は `KF: maxDCAV0 + maxDaughterDCA + minCosPointing`。

表を添付と同じ名称で読めるよう、KF側の一部の名称をHelix側に対応させている。**表の `maxDCAV0` をそのままKF YAMLのキーとして使わないこと。** 実装上の対応は以下。

| Parameter in table | Helix YAML key | KF YAML key | Definition |
|---|---|---|---|
| `minDCAProton` | `minDCAProton` | `imp5MinDCAProton` | Original PicoTrack `gDCA(pv.X(),pv.Y(),pv.Z())`; identical definition in both methods |
| `minDCAPion` | `minDCAPion` | `imp5MinDCAPion` | The same original PicoTrack gDCA definition |
| `maxDaughterDCA` | `maxDaughterDCA` | `maxDaughterDistance` | Helix distance between the two closest-approach points / KF daughter `GetDistanceFromParticle`; not identical reconstructed quantities |
| `maxDCAV0` | `maxDCAV0` | `maxDistanceToPv` | Helix distance to the PV calculated from the V0 midpoint and momentum / KF parent `GetDistanceFromVertex` |
| `minCosPointing` | `minCosPointing` | `minCosPointing` | The same pointing-cosine criterion, evaluated using each method's reconstructed decay vertex and momentum |
| `maxPathLength` | `maxPathLength` | `imp5MaxPathLength` | Absolute values of the original daughter-helix `pathLengths` in both methods; applied to both daughters |

`nSigmaProton`、`nSigmaPion`、`minNHitsFit`、`minNHitsRatio` は両YAMLで同名。PIDは保存済みTPC nSigmaで、両方式ともTOF・kaon PIDを使わない。hit ratioは旧Makerと同様に `nHitsMax > 0` のときだけ検査する。追加のpT・eta・nHitsDedx・dEdxError acceptance、HFT限定、anti-Λ再構成は今回のImp5比較では使わない。

### Retained selection in the original four KF columns

2026-09-14分類更新。上の4つのKF列は、以下の条件をすべて引き継いでいるが、**すべてをKF-specific cutとは呼ばない。** 今回は「既存Helixにはない共分散伝播・距離誤差・vertex fit等の追加計算を必要とする選別」を **KF-specific** と定義し、幾何・運動学的な条件、入力・fitの妥当性確認と分ける。これは今回の比較での実務上の呼び方であり、KFParticle以外のライブラリでは原理的に実装できないという意味ではない。

#### KF-specific cuts

以下は、既存Helixに閾値だけを追加しても適用できず、共分散伝播やfit計算の実装が必要な条件。`N/A` は旧Helixに対応する計算・判定がないことを表す。

| Parameter / check | Helix Imp5 | KF (all 4 columns) | Stage / meaning |
|---|---|---|---|
| `interfaceChiPrimaryCut` | N/A | 3.0 | Interface: scalar `chiPrimary < 3.0` classifies a track as primary. Lambda daughters must be secondary, requiring `chiPrimary ≥ 3.0` at this stage |
| `finderChiPrimary2D` | N/A | 3.0 | Finder: recomputed SIMD `chiPrimary > 3.0` for both secondary Lambda daughters |
| `finderChi2Ndf2D` | N/A | 10.0 | Finder: parent geometric-fit `chi2/NDF < 10.0` |
| `finderLdL2D` | N/A | 3.0 | Finder: decay-length significance `ldlMin = L/σL > 3.0`; distance uncertainty calculated from the PV/SV covariances |
| `chi2/NDF (PV-constrained copy)` | N/A | < 3.0 | TopoReconstructor: evaluated on a separate copy with a PV constraint; internal upstream condition |

`chiPrimary` は共分散を使ったdaughter–PVのχ²型適合度で、距離cmである元trackのgDCAとは別。Interfaceは `GetDeviationFromVertex` の返値を再度二乗せず分類に使う。InterfaceとFinderでは別の時点で計算するため、上表の `≥ 3` と `> 3` を区別する。閾値3を「3σ」や「χ² > 9」と読み替えない。

旧HelixのSVはdaughterの最近接点の中点であり、共分散付きの共通頂点fitではない。従って `L` は計算できるが、上表のfit `chi2/NDF` やSV誤差に基づく `L/σL` は、そのままでは得られない。

#### Geometry / kinematics (not KF-specific)

以下は位置・運動量から計算できるため、今回の定義では **KF-specificではない**。Helixに同種の選別を導入できるが、現行Imp5での適用状況は項目ごとに異なる。

| Parameter / check | Helix Imp5 | KF (all 4 columns) | Stage / meaning |
|---|---|---|---|
| `finderMaxDaughterDistance` [cm] | No Finder preselection; the corresponding final cut already requires daughter DCA ≤ 0.5 (`maxDaughterDCA`) | 1.5 | Finder: daughter distance after transport `dr < 1.5`; also distinct from the additional final KF `maxDaughterDCA` cut |
| `p1p2 > -p12` and `p1p2 > -p22` | This dot-product condition is not applied | required | Finder: daughter momentum dot products must satisfy `p1·p2 > −p1²` and `p1·p2 > −p2²`; internal upstream condition |
| `finderLCut` [cm] | `L` is calculated, but this lower-bound cut is not applied | 1.0 | Finder: PV–SV distance `lMin = L > 1.0` |
| `lMin < 200` [cm] | This upper-bound cut is not applied | required | Finder: PV–SV distance `L < 200`; internal upstream condition on a quantity distinct from `maxPathLength` |
| `isParticleFromVertex` | No identical API check; outward pointing is already required by cos(pointing angle) ≥ 0.998 (`minCosPointing`) | true | With the current `L/σL > 3` requirement, effectively `(SV−PV)·p > 0`; classified as geometric only under this combined selection |

`lMin` / `ldlMin` はFinderの `GetDistanceToVertexLine` に基づく量であり、元Helixのpath lengthや、PV制約後の `GetDecayLength` と混同しない。今回の外部PVは1個で、Finderの `L` は `L = |SV−PV|`、`σL` はPV/SV共分散から計算する距離誤差である。**`L` の上下限は幾何cut、`L/σL` の下限はKF-specific cut** と区別する。

`isParticleFromVertex` の一般形は「有効な分散かつ `L < 3σL`、または `(PV−SV)·p < 0`」で、一般には距離誤差を含む。今回はさらに `L/σL > 3` を課すため前半の許容条件は使われず、実効的にPVから外向きの `(SV−PV)·p > 0` が必要になる。この方向条件は、今回の全KF列の追加最終 `minCosPointing`（0.998または0.997）でも満たされる。旧Helixの `minCosPointing: 0.998` も、そのHelix再構成値で外向きを要求している。

daughter間距離やpointingにはこのように両方式で対応するcutがあるが、評価するSV・運動量・transport後の状態は同一ではない。同じ閾値を指定しても、同じ候補を選ぶことまでは保証しない。

#### KF input / fit validity checks

以下は上の物理的なKF-specific cutsとは分けて記載する。共分散品質は入力前処理、`chi2 > 0` はfitの数値的妥当性確認であり、距離やPIDの追加物理cutとは区別する。

| Parameter / check | Helix Imp5 | KF (all 4 columns) | Stage / meaning |
|---|---|---|---|
| `rejectBadCovariance` | Not applied | true | BuildTrack: reject Pico covariance matrices with the bad flag set |
| `maxPositionVariance` [cm²] | Not applied | 100.0 | BuildTrack: require `0 ≤ variance < 100.0` for each Cartesian position component |
| `maxMomentumVariance` [(GeV/c)²] | Not applied | 1.0 | BuildTrack: require `0 ≤ variance < 1.0` for each Cartesian momentum component |
| `chi2 > 0` | No corresponding parent fit | required | Finder: parent chi2 must be positive and finite; KF fit validity check |

入力共分散の3条件は、Helix側でも `TrackCovMatrix` を読み、必要なCartesian変換を追加すれば適用でき、KFParticleによるΛのfitは必要ない。従って独立した「KF input quality」として扱う。`chi2 > 0` はfit計算に依存するが、ここでは調整対象の物理cutとは分けている。

加えてKF入力の有限な状態・共分散、非ゼロ運動量、daughter ID/PDGの整合、正のNDFとmass/length error、候補量の有限性などの妥当性確認も維持している。これはTOF/PID等の追加物理cutではないが、旧Helixと同じ候補集合を保証するものでもない。

**追加Maker cutの `maxChi2Ndf` / `maxTopoChi2Ndf` 等が無効でも、KF-specific cuts表のFinder `chi2/NDF < 10` やTopo `chi2/NDF < 3` は残る。** `primaryProbCut` はこの外部Pico PV経路ではFinder閾値を後から `finderChiPrimary2D` で上書きするため、独立した有効cutとして数えていない。今回そろえた共通cutとは別にこれらのKF-specific cutsが有効なため、Imp5とすべての選別条件が同じ比較ではない。

### Mass range and shared event selection

| Parameter / check | Helix Imp5 | KF (all 4 columns) |
|---|---|---|
| `minMass` / `maxMass` (Maker) [GeV/c²] | no Maker mass cut | 1.05 / 1.25 |
| `massMin` / `massMax` (comparison QA) [GeV/c²] | 1.05 / 1.25 | 1.05 / 1.25 |
| `window` (Nwin) [GeV/c²] | [1.110, 1.122) | [1.110, 1.122) |
| `parent mass constraint` (plotted mass) | none | none |

旧Helixには同じMaker mass cutはないが、比較QAは両方式とも `1.05 ≤ M < 1.25` の同じbin範囲を使用。KF Makerの境界判定自体は `1.05 ≤ M ≤ 1.25`。主窓のNwin・sideband積分は上限を含まない共通bin定義である。

Finder内部のPDG massからの幅による `saveMother` 選別やmass constraintは、二次解析用の別候補コンテナに対するもの。今回のraw Λは `saveParticle` を満たして先に `Particles` に保存されたものを取得しており、そこにそのPDG mass窓や親mass constraintを適用しているとはしない。ただし、上表のFinder/Topo選別は適用済み。

eventは前回と同じImp5 policy（`maxNTr: 0`、追加Vz/Vr/RefMult/VPD cutなし、共通centrality選別）で、両方式とも10,000入力・9,744再構成イベント。今回そろえたのは共通PID/track条件と主窓の候補数であり、**全cutが同一、KF固有cutがゼロ、真の再構成効率が同一という意味ではない。** Finder/Topo以前に失われた候補は保存済みtreeから回復できないため、KF固有cutを外す検証には別設定での再処理が必要になる。

### Sources and scope of this update

- [Helix cut YAML](../../../config/maker/maker_auau13p5_anaLambda.yaml)、[KF Imp5 cut YAML](../../../config/cuts/kf/kf_auau13p5_anaLambda_KFParticle_Imp5.yaml)、[scan QA YAML](../../../config/qa/lambda_kf_imp5_cut_scan_20260913.yaml)。
- 調整後の値はYAML本番設定を書き換えた値ではなく、[保存済みscanの代表条件](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_cut_scan_20260913/main/scan_summary.txt) と [CSV](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_cut_scan_20260913/main/scan_metrics.csv) の値。
- 定義の照合：[StLambdaMaker](../../../StMaker/StLambdaMaker/StLambdaMaker.cxx)、[StPicoKFParticleInterface](../../../StMaker/kfparticle/StPicoKFParticleInterface.cxx)、[KfParticleHelper](../../../StMaker/kfparticle/KfParticleHelper.cxx)、[StLambdaKFParticleMaker](../../../StMaker/StLambdaKFParticleMaker/StLambdaKFParticleMaker.cxx)、[KFParticleFinder](../../../StRoot/KFParticle/KFParticleFinder.cxx)、[KFParticleTopoReconstructor](../../../StRoot/KFParticle/KFParticleTopoReconstructor.cxx)、[KFParticleSIMD](../../../StRoot/KFParticle/KFParticleSIMD.cxx)。
- 2026-09-14の先行する分類・表更新では文書のみを変更した。その後の実装・再処理は末尾の新study節に分離して記載する。

7月1日の歴史的Imp5 fit S/N=5.217等は716,241イベントの結果であり、今回の10,000イベントのS/Bとは混同しない。9月8日時点の「Notion本文未照合」は今回の指定表の確認で解消した。過去ノート・既存設定の履歴コメント自体は変更していない。

## 図・詳細・出力

2026-09-13の表示更新：`paired_cut_overlay`、`count_matched_overlay`、`pv_only_overlay` は左を無規格化の候補数、右を各曲線の通常bin合計が1になるNormalize表示とした。左右とも1.05 ≤ M < 1.25 GeV/c²、1 MeV/c² bin。右は `Scale(1/Integral(1,Nbins))`（bin幅係数なし、underflow/overflow除外）で、元Count histogram・cut・質量窓・S/Bの数値は変更しない。補助fit図とPV距離scan曲線は今回の表示更新の対象外。

- [主要比較図：pointingを維持した2距離cut緩和（PNG）](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/paired_cut_overlay.png) / [PDF](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/paired_cut_overlay.pdf)
- [PV距離だけを緩めたときの候補数とS/B](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/pv_only_curve.png)
- [3条件scanの最接近点の比較図](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/count_matched_overlay.png)
- [詳細解析ノート：手法、fit cross-check、制約、再現手順](../../../analysisnote/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/note.md)
- [360条件のCSV](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_cut_scan_20260913/main/scan_metrics.csv) / [ヒストグラムとcanvasのROOT](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_cut_scan_20260913/main/scan.root)
- [出力の条件・用途索引](../../../rootfile/auau13p5_anaLambda_KFParticle/mdfiles/cut_scan_20260913.md)

更新前の図・成果物の復元先は [before_normalized_overlays_20260913.tar.gz](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/before_normalized_overlays_20260913.tar.gz)。元provenanceは履歴として保持し、表示更新分は [provenance_normalized_20260913.tar.gz](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/provenance_normalized_20260913.tar.gz) を参照。

2026-09-13の距離cut scanで追加したのは保存済み候補を再選別するQAコード・QA設定・成果物のみ。この時点では既存 `StLambdaMaker` / `StLambda` / KF再構成コード・本番YAML・元ROOT・`.current_mainconf` は変更していない。PicoDstの再処理、farm job投入、git commit/pushはしていない。

全raw候補が保存されているので、今回の追加最終cutの変更はPicoDstを読み直さずに評価できる。一方、保存前の入力PID・共分散・Finder/Topoで落ちた候補は回復できない。**今回の範囲でS/B改善が見えなかったことは、KFParticle全体の性能限界を意味しない。** qhu1との比較は依頼どおり中断した。

## Geometry-off KF-specific scan (2026-09-14)

### Result and interpretation

依頼に従い、上記の **Geometry / kinematics (not KF-specific)** に挙げたFinder内部の5条件をΛ/anti-Λについて無効化できる設定を追加した。この節は、保存済みtreeの最終距離cutを変えた前節とは別に、同じPicoDstの先頭10,000イベントを上流から再処理したstudyである。全runで9,744イベントを再構成し、最終Imp5の距離・pointing・PID・track条件は固定した。

**幾何cutの無効化だけではNwinは1,370のまま。KF-specificの5閾値もかなり広く緩めた最接近点で1,613（90.9%）、Helixの1,775より162少なく、同数条件は今回の探索では見つからなかった。** 最接近点のsideband S/Bは5.98 ± 0.43、Helixは5.20 ± 0.35。ただし候補数がそろっていないので、「同数でS/Bが改善した」とは結論しない。

| Selection | applyLambdaGeometryCuts | Nwin | Nwin / Helix | S/B (sideband) |
|---|---|---:|---:|---:|
| Helix Imp5 | N/A | 1775 | 100.0% | 5.20 ± 0.35 |
| KF baseline | true | 1370 | 77.2% | 8.56 ± 0.74 |
| KF: applyLambdaGeometryCuts=false | false | 1370 | 77.2% | 8.51 ± 0.73 |
| KF: topoChi2NdfCut | false | 1491 | 84.0% | 8.75 ± 0.73 |
| KF: finderChi2Ndf2D + topoChi2NdfCut | false | 1506 | 84.8% | 8.34 ± 0.68 |
| KF: chiPrimary + chi2/NDF + L/sigmaL (moderate) | false | 1606 | 90.5% | 6.09 ± 0.44 |
| KF: chiPrimary + chi2/NDF + L/sigmaL (loose) | false | 1613 | 90.9% | 5.98 ± 0.43 |

ここでの最接近点は **実行した5つのgeometry-off条件の中** で選んだもので、全パラメータ空間の最適点・理論的な収量上限ではない。選択基準は主窓の候補数差、同点なら変更KF閾値数、さらに同点なら一覧順。S/Bを最大にする点を選んだものではない。

`1e6` は非常に緩い有限の上限であり、cutの数学的な完全削除ではない。`interfaceChiPrimaryCut = 1e-6` も正の閾値を残し、`finderChiPrimary2D = 0` と `finderLdL2D = 0` のFinder判定は厳密な `> 0`。従って「全KF cutをゼロにした」とは表現しない。入力共分散・fit妥当性、再構成量の定義の違いも残るため、この結果を「KFParticleパッケージを使うと必ず効率が落ちる」と一般化しない。

Nwinは背景を含む候補数であって真のΛ効率ではない。背景を含む主窓の一致が未達のため、追加のcutを勝手に緩めて一致させることはしていない。今回の通常Imp5距離cutは、以前のscanで候補数を回復させた `0.8 / 0.6 cm` へ変更せず、`0.5 / 0.5 cm` のまま保持した。

### Shared Imp5 cuts — all values repeated

`KF: KF-specific cuts (loose)` は上の最接近点 `all_specific_loose` を指す。表の別名とYAMLキーの対応は前節と同じ。

| Parameter | Helix Imp5 | KF baseline | KF: applyLambdaGeometryCuts=false | KF: KF-specific cuts (loose) |
|---|---:|---:|---:|---:|
| `nSigmaProton` | 3.0 | 3.0 | 3.0 | 3.0 |
| `nSigmaPion` | 3.0 | 3.0 | 3.0 | 3.0 |
| `minDCAProton` [cm] | 0.7 | 0.7 | 0.7 | 0.7 |
| `minDCAPion` [cm] | 1.0 | 1.0 | 1.0 | 1.0 |
| `maxDaughterDCA` [cm] | 0.5 | 0.5 | 0.5 | 0.5 |
| `maxDCAV0` [cm] | 0.5 | 0.5 | 0.5 | 0.5 |
| `minNHitsFit` | 15 | 15 | 15 | 15 |
| `minNHitsRatio` | 0.52 | 0.52 | 0.52 | 0.52 |
| `minCosPointing` | 0.998 | 0.998 | 0.998 | 0.998 |
| `maxPathLength` [cm] | 100.0 | 100.0 | 100.0 | 100.0 |

### KF-specific thresholds tested

各列は実際の再処理設定。幾何cutは最初の `KF baseline` だけ有効で、右の5列では無効。Helixにはこの5閾値に対応する共分散・fit選別はない。新設 `topoChi2NdfCut` は、従来の内部固定条件 `chi2/NDF (PV-constrained copy) < 3` をΛだけ設定可能にしたもの。

| Parameter | KF baseline | geometry_off | topo_only_loose | fit_loose | moderate_specific_loose | all_specific_loose |
|---|---:|---:|---:|---:|---:|---:|
| `interfaceChiPrimaryCut` | 3 | 3 | 3 | 3 | 1 | 1e-6 |
| `finderChiPrimary2D` | 3 | 3 | 3 | 3 | 1 | 0 |
| `finderChi2Ndf2D` | 10 | 10 | 10 | 100 | 100 | 1e6 |
| `finderLdL2D` | 3 | 3 | 3 | 3 | 1 | 0 |
| `topoChi2NdfCut` | 3 | 3 | 1e6 | 1e6 | 1e6 | 1e6 |
| `applyLambdaGeometryCuts` | true | false | false | false | false | false |

### Geometry gates disabled in this study

これは **Finder内部** の切替であり、通常Imp5に対応するMaker最終cutの無効化ではない。位置・運動量・共分散の計算やtransport自体は残し、以下の選別マスクだけをΛ/anti-Λについてバイパスする。このstudyはΛのみを再構成し、他PDGのマスクを直接バイパスしない。

| Parameter / check | KF baseline | All geometry-off trials | Classification |
|---|---|---|---|
| `finderMaxDaughterDistance` [cm] | dr < 1.5 | Disabled | Finder geometry preselection |
| `p1p2 > -p12` and `p1p2 > -p22` | Required | Disabled | Daughter momentum dot products |
| `finderLCut` [cm] | L > 1.0 | Disabled | Lower PV-SV distance bound |
| `lMin < 200` [cm] | Required | Disabled | Upper PV-SV distance bound |
| `isParticleFromVertex` | Required | Disabled | Vertex-line direction / compatibility gate |

`finderMaxDaughterDistance` と `finderLCut` はsnapshotに元の値が残るが、`applyLambdaGeometryCuts: false` のΛ選別には適用されない。前節の `isParticleFromVertex` を純粋な外向き条件とみなす説明は、元の `L/σL > 3` との組み合わせに対するもの。今回の緩和studyではその有意度も変えるので、同じ簡略化を適用せず、当該gateを明示的にバイパスしたと記載する。

入力・fit妥当性は全runで固定：`rejectBadCovariance: true`、Cartesian位置対角分散 `0 ≤ variance < 100 cm²`、運動量対角分散 `0 ≤ variance < 1 (GeV/c)²`、正で有限なfit χ²、正のNDF・error、各候補量の有限性などは維持した。親の質量制約は図のmassに適用しない。共通event/centrality、TPC-only PID、anti-Λ無効も変更しない。

### Windows and cross-checks

主窓は前節と同じ `[1.110, 1.122)`、既定sidebandは `[1.098, 1.106)` と `[1.126, 1.134)`。局所一定密度を仮定し、`B = 0.75 × (Nleft + Nright)`、`S = Nwin − B`、`S/B` とした。最接近点は `Nleft = 145`、`Nright = 163`、`B = 231`、`S = 1382`。Helixは `B = 286.5`、`S = 1488.5`。これも背景モデルに依存する推定であり真の効率測定ではない。

| Cross-check | Helix Imp5 | KF: KF-specific cuts (loose) |
|---|---:|---:|
| Nwin, [1.110, 1.122) | 1775 | 1613 |
| S/B, default sidebands | 5.20 ± 0.35 | 5.98 ± 0.43 |
| S/B, near sidebands | 5.43 ± 0.37 | 5.87 ± 0.42 |
| S/B, far sidebands | 5.70 ± 0.39 | 6.49 ± 0.48 |
| Nwin, [1.112, 1.120) | 1627 | 1469 |
| S/B, narrow window / default sidebands | 7.52 ± 0.48 | 8.54 ± 0.60 |

誤差は前節と同じ独立Poissonの目安。共有イベント・候補間相関、背景形状、同じデータで条件を選んだ影響を含まず、この数値だけで差の有意性は主張しない。

### Reversibility and validation

新しい [study mainconf](../../../config/mainconf/main_auau13p5_anaLambda_KFParticle_KFSpecificStudy.yaml) は、**geometry-offだけを行い、5つのKF-specific閾値は元の値を維持する入口** として残した。非常に緩い探索点を本番デフォルトにはしない。最接近点を再現する場合は [all_specific_loose_main.yaml](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_kf_specific_scan_20260914/configs/all_specific_loose_main.yaml) と対応するKF snapshotを使用する。

元のcutに戻すには、新しいROOTプロセスで [元のImp5 mainconf](../../../config/mainconf/main_auau13p5_anaLambda_KFParticle_Imp5.yaml) を指定するだけでよい。再ビルド不要。新規キー省略時の既定値は `applyLambdaGeometryCuts: true`、`topoChi2NdfCut: 3`。緩和YAMLでgeometryだけtrueに戻しても、ほかの閾値は元に戻らないため、完全復帰には元mainconfを使う。

- ROOT 5.34/38・STAR SL24y・GCC 4.8.5でKF関連ビルドと人工データのfull-chain / Pico adapterテスト成功。ROOT 6は導入していない。
- importerの再現性検査は39 files / 11 translation unitsで一致。ローカル拡張を生成手順とPROVENANCEに記録。
- 元Imp5で10,000イベントを再実行し、旧KF出力と68ヒストグラムの2,536,396 cells、および候補treeの42,871 entries × 42 leavesが完全一致。新規設定文字列があるためROOTファイル自体のbyte一致ではない。
- 各study出力で10,000入力・9,744再構成、完成した非recovered ROOT、構造・候補量・counter QAを確認。元出力との有効な共通cut・event設定一致も検査した。入力順は別途保存した7ファイルのリスト、SHA256、10,000件のrun/event ID監査で追跡する。
- 旧 `StLambda` / `StLambdaMaker` / Helixライブラリ、元Imp5 YAMLを変更せず、`.current_mainconf` も元の内容を維持。farm job投入・git commit/pushなし。

幾何プリフィルタを外すと処理コストは増える。今回のROOT内計測ではdefault約32秒、geometry-off約60秒、最も広い緩和約1,506秒。並列実行中の実測であり厳密なbenchmarkではないが、広い緩和を通常運転の推奨設定とはしない。

### New figure and outputs

- [New comparison PNG](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/kf_specific_cut_overlay.png) / [PDF](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/kf_specific_cut_overlay.pdf)：Helix Imp5、KF baseline、実行点中の最接近条件。左は実候補数、右は `[1.05, 1.25)` の全通常bin積分でNormalize。以前の図を上書きしていない。
- [Metrics CSV](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_kf_specific_scan_20260914/comparison/kf_specific_cut_scan_metrics.csv) / [Summary](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_kf_specific_scan_20260914/comparison/kf_specific_cut_scan_summary.txt) / [Comparison ROOT](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_kf_specific_scan_20260914/comparison/kf_specific_cut_scan.root)
- [Output index and run conditions](../../../rootfile/auau13p5_anaLambda_KFParticle/mdfiles/kf_specific_scan_20260914.md) / [Reproduction and provenance](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_kf_specific_scan_20260914/provenance/README.md)
- [New QA macro](../../../common/macro/compareLambdaKfSpecificCuts.C) / [Default-output equality test](../../../tests/compare_kfparticle_default_outputs.C)
