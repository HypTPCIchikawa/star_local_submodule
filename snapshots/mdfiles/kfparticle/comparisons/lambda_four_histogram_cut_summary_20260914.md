# AuAu13p5: Four-histogram cut and result summary (2026-09-14)

対象：[4系列の統合図（PDF）](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/combined_cut_overlay_readable_1p08_1p16.pdf) / [PNG](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/combined_cut_overlay_readable_1p08_1p16.png)。

同じ先頭10,000入力イベント（9,744再構成イベント）、Λ → pπ⁻のみ。共通のcutも省略せず記載した。最後の2列は**別々の調整条件**であり、両方の変更を同時適用した条件ではない。

## Cuts and results

| Parameter / result | Helix Imp5 | KF baseline | KF: maxDCAV0 + maxDaughterDCA + minCosPointing | KF: chiPrimary + chi2/NDF + L/sigmaL (geometry off) |
|---|---:|---:|---:|---:|
| Plot style | Charcoal circles | Blue dashed line | Vermilion squares | Green triangles |
| `nSigmaProton` | 3.0 | 3.0 | 3.0 | 3.0 |
| `nSigmaPion` | 3.0 | 3.0 | 3.0 | 3.0 |
| `minDCAProton` [cm] | 0.7 | 0.7 | 0.7 | 0.7 |
| `minDCAPion` [cm] | 1.0 | 1.0 | 1.0 | 1.0 |
| `maxDaughterDCA` [cm] | 0.5 | 0.5 | 0.6 | 0.5 |
| `maxDCAV0` [cm] | 0.5 | 0.5 | 0.75 | 0.5 |
| `minNHitsFit` | 15 | 15 | 15 | 15 |
| `minNHitsRatio` | 0.52 | 0.52 | 0.52 | 0.52 |
| `minCosPointing` | 0.998 | 0.998 | 0.997 | 0.998 |
| `maxPathLength` [cm] | 100.0 | 100.0 | 100.0 | 100.0 |
| TOF PID | Not used | Not used | Not used | Not used |
| `interfaceChiPrimaryCut` (KF-specific) | N/A | 3.0 | 3.0 | 1e-6 |
| `finderChiPrimary2D` (KF-specific) | N/A | 3.0 | 3.0 | 0 |
| `finderChi2Ndf2D` (KF-specific) | N/A | 10.0 | 10.0 | 1e6 |
| `finderLdL2D` (KF-specific) | N/A | 3.0 | 3.0 | 0 |
| `topoChi2NdfCut` (KF-specific) | N/A | 3.0 | 3.0 | 1e6 |
| `applyLambdaGeometryCuts` (not KF-specific) | N/A | true | true | false |
| `rejectBadCovariance` (input quality) | N/A | true | true | true |
| `maxPositionVariance` [cm²] (input quality) | N/A | 100.0 | 100.0 | 100.0 |
| `maxMomentumVariance` [(GeV/c)²] (input quality) | N/A | 1.0 | 1.0 | 1.0 |
| **Nwin [1.110, 1.122) GeV/c²** | **1775** | **1370** | **1783** | **1613** |
| **Nwin / Helix** | **100.0%** | **77.2%** | **100.5%** | **90.9%** |
| **S/B (sideband)** | **5.20 ± 0.35** | **8.56 ± 0.74** | **4.94 ± 0.33** | **5.98 ± 0.43** |

## Reading the table

- `maxDCAV0` / `maxDaughterDCA` は比較用のHelix側名称。KF YAMLではそれぞれ `maxDistanceToPv` / `maxDaughterDistance`。両方式で評価する再構成量は同一ではない。KFのdaughter gDCAとhelix path lengthは `imp5MinDCAProton/Pion`、`imp5MaxPathLength` に対応する。
- KF-specific判定は、Interface `chiPrimary ≥ cut`、Finder `chiPrimary > cut`、親fit `chi2/NDF < cut`、`L/σL > cut`、PV拘束した別copyの `chi2/NDF < topoChi2NdfCut`。緑の `1e6` は緩い有限上限、`0` は厳密な `> 0` を残すため、KF固有cutの完全削除ではない。
- `applyLambdaGeometryCuts: false` はFinder内部の5項目（daughter距離 `dr < 1.5 cm`、daughter運動量内積条件、`L > 1 cm`、`L < 200 cm`、`isParticleFromVertex`）を無効化する。通常Imp5に対応する最終距離・pointing cutは、緑でも表の値を維持する。
- 全KF系列で入力共分散とfit妥当性の確認を維持。位置・運動量の対角分散は0以上かつ表の上限未満。図の親massに質量制約はかけない。
- 主窓のNwinは**背景を含む候補数で、真のΛ収量や効率ではない**。S/Bはsideband `[1.098,1.106)` と `[1.126,1.134)` GeV/c²から `B = 0.75 × (Nleft + Nright)`、`S = Nwin − B` として推定。誤差は独立Poisson近似の統計目安で、系統誤差・共有候補相関・同じデータでcutを選んだ影響は含まない。
- 図の表示範囲は1.08–1.16 GeV/c²、1 MeV bin。右panelの規格化範囲は従来どおり `[1.05,1.25)`。表の候補数・S/BはNormalize前の値。

橙はHelixとほぼ同数だが、S/Bの改善は確認できていない。緑はS/Bの推定値がHelixより高いものの候補数は90.9%であり、同数条件での改善を確認した結果ではない。

## Sources

- Cuts: [Helix Imp5 YAML](../../../config/maker/maker_auau13p5_anaLambda.yaml), [KF baseline YAML](../../../config/cuts/kf/kf_auau13p5_anaLambda_KFParticle_Imp5.yaml), [KF-specific loose snapshot](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_kf_specific_scan_20260914/configs/all_specific_loose_kf.yaml).
- Results: [Distance/pointing scan CSV](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_cut_scan_20260913/main/scan_metrics.csv) — `point=86, window=primary, sidebandVariation=default`; [KF-specific scan CSV](../../../rootfile/auau13p5_anaLambda_KFParticle/Imp5_kf_specific_scan_20260914/comparison/kf_specific_cut_scan_metrics.csv) — `helix / baseline / all_specific_loose, primary, default`.
- [Detailed comparison and definitions](lambda_helix_kf_imp5_cut_scan_20260913.md). 橙の値は保存済みscanの選択点であり、本番YAMLを書き換えた値ではない。この要約作成では解析・cut・図を変更していない。
