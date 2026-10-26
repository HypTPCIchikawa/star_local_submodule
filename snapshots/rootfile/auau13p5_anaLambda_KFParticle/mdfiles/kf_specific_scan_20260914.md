# 2026-09-14：Imp5幾何条件を固定したKF-specific study出力

出力ディレクトリは [Imp5_kf_specific_scan_20260914/](../Imp5_kf_specific_scan_20260914/)。同じ先頭10,000入力イベントを再処理し、通常Imp5の追加最終cutを固定したまま、上流KF条件を調べる。既存ROOTを置換せず、各runを別名で保存する。

主窓は **1.110 ≤ M(pπ⁻) < 1.122 GeV/c²**。比較対象の [旧Helix Imp5](../../auau13p5_anaLambda/local_compare_10000_20260907.root) は `Nwin = 1775`。今回探索したgeometry-offの5条件で最接近した `all_specific_loose` は **1613候補（旧Helixの90.9%、162候補不足）で、同数には届かなかった**。これは背景を含む候補数で、真のΛ収量や効率ではない。探索した条件内での最接近点であり、数学的な候補数上限は主張しない。

## Run索引

下表の全runが正常完了し、すべて10,000入力／9,744再構成イベント。同名の設定snapshotは [configs/](../Imp5_kf_specific_scan_20260914/configs/) の `<run>_main.yaml` と `<run>_kf.yaml`。ただし `default_regression` は元のImp5 mainconfを使う。

表の `Interface` は `interfaceChiPrimaryCut`、`Finder primary` は `finderChiPrimary2D`、`Finder fit` は `finderChi2Ndf2D`、`L/σL` は `finderLdL2D`、`Topo` は `topoChi2NdfCut` の設定値。`Geometry` は `applyLambdaGeometryCuts`。

| Run / ROOT output | Geometry | Interface | Finder primary | Finder fit | L/σL | Topo | Nwin | Status |
|---|---|---:|---:|---:|---:|---:|---:|---|
| [default_regression](../Imp5_kf_specific_scan_20260914/default_regression.root) | true | 3 | 3 | 10 | 3 | 3 | 1370 | Complete; exact baseline regression |
| [probe_loose](../Imp5_kf_specific_scan_20260914/probe_loose.root) | true | 1e-6 | 0 | 1e6 | 0 | 3 | 1468 | Complete; original upstream Topo limit retained |
| [geometry_off](../Imp5_kf_specific_scan_20260914/geometry_off.root) | false | 3 | 3 | 10 | 3 | 3 | 1370 | Complete |
| [topo_only_loose](../Imp5_kf_specific_scan_20260914/topo_only_loose.root) | false | 3 | 3 | 10 | 3 | 1e6 | 1491 | Complete |
| [fit_loose](../Imp5_kf_specific_scan_20260914/fit_loose.root) | false | 3 | 3 | 100 | 3 | 1e6 | 1506 | Complete |
| [moderate_specific_loose](../Imp5_kf_specific_scan_20260914/moderate_specific_loose.root) | false | 1 | 1 | 100 | 1 | 1e6 | 1606 | Complete |
| [all_specific_loose](../Imp5_kf_specific_scan_20260914/all_specific_loose.root) | false | 1e-6 | 0 | 1e6 | 0 | 1e6 | 1613 | Complete; closest among the five geometry-off trials |

`default_regression` は [既存Imp5 KF出力](../Imp5/local_compare_10000_20260908.root) と、68ヒストグラムの2,536,396 cells、および42枝×42,871 entriesの候補treeが完全一致した。ROOTファイル全体のbyte一致を意味しない。

## 候補数とS/B

正式比較の主窓・default sidebandによる結果。`Raw candidates` はTopo後の保存候補、`Selected candidates` は追加Maker選別後の全候補であり、主窓の `Nwin` とは異なる。

| Sample | Raw candidates | Selected candidates | Nwin | Sideband S/B |
|---|---:|---:|---:|---:|
| Helix Imp5 | N/A | 4914 | 1775 | 5.195 ± 0.349 |
| geometry_off | 74131 | 3137 | 1370 | 8.514 ± 0.733 |
| topo_only_loose | 1356206 | 3326 | 1491 | 8.745 ± 0.727 |
| fit_loose | 2811947 | 3386 | 1506 | 8.340 ± 0.681 |
| moderate_specific_loose | 3699047 | 4314 | 1606 | 6.091 ± 0.445 |
| all_specific_loose | 9423004 | 4405 | 1613 | 5.983 ± 0.434 |

代表点は主窓候補数とHelixとの差が最小になる条件で選び、S/Bは最適化していない。S/Bは局所的に一定の背景密度を仮定するsideband推定で、過去のgaus+pol2 fit値とは異なる。誤差は独立Poisson近似の統計目安で、共有イベント・候補間相関、背景形状、同じデータで条件を選んだ影響を含まない。上表は表示用に丸めており、全精度と窓・sideband変更時の結果は正式CSVに保存している。

## 固定条件と今回の切替

全runで追加Maker選別は `maxDistanceToPv = 0.5 cm`、`maxDaughterDistance = 0.5 cm`、`minCosPointing = 0.998` を維持する。TPC PID、hits、元trackのproton/pion gDCA下限0.7/1.0 cm、元helixの絶対path length上限100 cm、event/centrality条件、mass範囲1.05–1.25 GeV/c²、anti-Λ無効も元Imp5設定を維持する。

`applyLambdaGeometryCuts: false` はΛ/anti-ΛのFinder内部に限って、transport後daughter距離、娘運動量の内積、PV–SV距離の下限・200 cm上限、`isParticleFromVertex` の5項目をバイパスする。上記の通常Imp5追加Maker cutを外す設定ではない。

共分散品質は緩めない。`rejectBadCovariance: true`、位置対角分散 `0 ≤ variance < 100 cm²`、運動量対角分散 `0 ≤ variance < 1 (GeV/c)²`、有限値・正のfit/errorなどの妥当性確認を維持する。`L/σL` は独立したKF-specific条件として設定値で選別する。geometry切替が有意度や共分散計算を無効化することはない。

## Study設定と元Imp5への復帰

Study入口は [main_auau13p5_anaLambda_KFParticle_KFSpecificStudy.yaml](../../../config/mainconf/main_auau13p5_anaLambda_KFParticle_KFSpecificStudy.yaml)。この標準study設定は **geometry-offだけを適用し、KF-specificの5閾値は元の3、3、10、3、3を保持**する。最接近点の極端に緩い閾値を標準設定には採用しない。

最接近点を再現するときは [configs/all_specific_loose_main.yaml](../Imp5_kf_specific_scan_20260914/configs/all_specific_loose_main.yaml) と、その参照先 [all_specific_loose_kf.yaml](../Imp5_kf_specific_scan_20260914/configs/all_specific_loose_kf.yaml) を使う。各出力の条件確認では設定snapshotとROOT内 `KFParticleEffectiveConfiguration` を優先する。

元の挙動へ戻すには、runnerに [main_auau13p5_anaLambda_KFParticle_Imp5.yaml](../../../config/mainconf/main_auau13p5_anaLambda_KFParticle_Imp5.yaml) を指定する。この元mainconfと参照YAML、旧Helix解析は変更していない。新規2キーを省略した既定値は `topoChi2NdfCut = 3`、`applyLambdaGeometryCuts = true`。元Imp5 mainconfを使えば他のKF設定も元に戻る。

## 正式な比較出力と再現記録

- 比較図：[PNG](../Imp5_kf_specific_scan_20260914/comparison/kf_specific_cut_overlay.png) / [PDF](../Imp5_kf_specific_scan_20260914/comparison/kf_specific_cut_overlay.pdf)
- 比較ヒストグラムとmetadata：[kf_specific_cut_scan.root](../Imp5_kf_specific_scan_20260914/comparison/kf_specific_cut_scan.root)
- 数値と判定：[kf_specific_cut_scan_metrics.csv](../Imp5_kf_specific_scan_20260914/comparison/kf_specific_cut_scan_metrics.csv) / [kf_specific_cut_scan_summary.txt](../Imp5_kf_specific_scan_20260914/comparison/kf_specific_cut_scan_summary.txt)
- 比較対象5run：[runs.tsv](../Imp5_kf_specific_scan_20260914/runs.tsv)
- 保存した入力・ログ・補助処理の索引：[provenance/README.md](../Imp5_kf_specific_scan_20260914/provenance/README.md)。同配下の [input/](../Imp5_kf_specific_scan_20260914/provenance/input/)、[logs/](../Imp5_kf_specific_scan_20260914/provenance/logs/)、[helpers/](../Imp5_kf_specific_scan_20260914/provenance/helpers/) に入力順序・イベント同一性の記録、実行／回帰／テストログ、補助処理を保存。

比較には実際に再処理した各runの選別後ヒストグラムを使用した。保存済みtreeから上流で失われた候補を回復した結果ではない。正式比較の固定cut・event metadata照合は通過している。

前回の追加最終geometry cut scanは [2026-09-13出力索引](cut_scan_20260913.md)、cutの定義と既存比較は [比較本文](../../../mdfiles/kfparticle/comparisons/lambda_helix_kf_imp5_cut_scan_20260913.md) を参照。
