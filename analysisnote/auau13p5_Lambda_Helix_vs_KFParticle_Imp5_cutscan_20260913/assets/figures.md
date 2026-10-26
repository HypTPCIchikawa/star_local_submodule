# 図の索引

図の実体は `share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/`。

- [paired_cut_overlay.png](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/paired_cut_overlay.png)：Helix、KF baseline、PV0.8/daughter0.6/cos0.998。
- [count_matched_overlay.png](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/count_matched_overlay.png)：3条件gridで候補数が最も近い点との比較。
- [pv_only_overlay.png](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/pv_only_overlay.png)：PV距離だけを3 cmへ緩めた比較。
- [pv_only_curve.png](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/pv_only_curve.png)：PV距離だけを変えたNwinとsideband推定S/B。x配置はカテゴリで、数値間隔は等間隔ではない。
- [fit_mass_range00.png](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/fit_mass_range00.png)、[range01](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/fit_mass_range01.png)：gaus+pol2の2つのfit範囲。数値収束とモデル適合性は別であり、詳細ノートの注意を参照。

3種類のoverlayは、左が無規格化の候補数、右が各曲線の通常bin合計を1にするNormalize表示。左右とも1.05 ≤ M < 1.25 GeV/c²、1 MeV/c² bin。右は `Scale(1/Integral(1,Nbins))`（bin幅係数なし、underflow/overflow除外）を表示用histogramに適用する。元Count histogram・cut・質量窓・S/Bは変更しない。PV距離scan曲線と無規格化の補助fit図は今回の表示更新の対象外。

各PNGと同名のPDFも保存。scan ROOTに元histogramとcanvas、fit ROOTにfit入力histogram・TFitResult・canvasを保存した。更新前の図・成果物は [before_normalized_overlays_20260913.tar.gz](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/before_normalized_overlays_20260913.tar.gz) から復元できる。元provenanceは履歴として保持し、表示更新分は [provenance_normalized_20260913.tar.gz](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/provenance_normalized_20260913.tar.gz) を参照。
