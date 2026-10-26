# KF Lambda purity optimization: AuAu13p5 / AuAu3p85 — 2026-09-28

## 1. Conclusion

両energyで、同じ先頭10,000入力イベントを使った上流再構成cutを含む探索を行い、最適条件を固定した。その後、探索と重複しない100,000入力イベントでdefault / optimizedを実際に再構成し、比較を完了した。

| FXT beam energy | Default purity | Optimized purity | Signal retention | Raw window retention |
|---|---:|---:|---:|---:|
| 13.5 GeV | 90.42% | 98.67% ± 0.12 pp | 57.42% ± 0.50 pp | 52.62% |
| 3.85 GeV | 81.69% | 97.73% ± 0.39 pp | 53.99% ± 1.39 pp | 45.13% |

誤差は400回のpaired-event bootstrapによる統計誤差。`Signal retention` は背景差引き `S/Sdefault` であり、検出効率ではない。3p85の候補数比45.13%には除去した背景も含まれるため、信号保持率53.99%と区別する。

今回の「信号約半分を許容してpurityを優先」という方針に沿った候補が得られた。ただし、探索で見えた3p85の99.39%を最終結果に使わず、**独立100kの97.73%**を採用する。背景モデルを変えた感度範囲は13p5で98.14–98.79%、3p85で96.24–98.06%であり、「常に98%以上」「BG free」の保証ではない。

既存の標準設定、Helix/KF Maker、共有ライブラリは変更していない。新規study用設定・QA補助コード・結果だけを追加した。ROOT5.34/38 / STAR SL24yでローカル実行し、farm job、commit、pushは行っていない。

## 2. Final figures and ROOT files

表示は前回の指定を維持し、横軸1.10–1.14 GeV/c²、左がraw counts、右が[1.05,1.25)全積分でNormalize。元の1 MeV/c² binを保ち、rebin・smoothing・parent mass constraintは加えていない。両曲線にPoisson error barsを表示した。

| Sample | PDF | PNG | Default ROOT | Optimized ROOT |
|---|---|---|---|---|
| AuAu13p5, 100k | [PDF](../../../share/figure/auau13p5_Lambda_KFParticle_purity_upstream_20260928/validation100k_default_vs_best.pdf) | [PNG](../../../share/figure/auau13p5_Lambda_KFParticle_purity_upstream_20260928/validation100k_default_vs_best.png) | [ROOT](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/validation100k/default.root) | [ROOT](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/validation100k/best50.root) |
| AuAu3p85, 100k | [PDF](../../../share/figure/auau3p85_Lambda_KFParticle_purity_upstream_20260928/validation100k_default_vs_best.pdf) | [PNG](../../../share/figure/auau3p85_Lambda_KFParticle_purity_upstream_20260928/validation100k_default_vs_best.png) | [ROOT](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/validation100k/default.root) | [ROOT](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/validation100k/best50.root) |

## 3. Input samples and event selection

XRootDのport1095を使わず、読み出せるNFSのPicoDstを使用した。13p5の探索は以前の先頭7ファイルと同じ入力順。3p85は今回アクセス可能だったrun20160024を使用した。

| Dataset | Production | NFS mount | Run | Discovery input / reconstructed | Validation input / reconstructed |
|---|---|---|---|---:|---:|
| AuAu13p5 | production_13p5GeV_fixedTarget_2020, P24iy | /star/data22 | 21033026 | 10,000 / 9,744 | 100,000 / 97,278 |
| AuAu3p85 | production_3p85GeV_fixedTarget_2019, P24iy | /star/data24 | 20160024 | 10,000 / 9,505 | 100,000 / 95,601 |

default / optimizedは各energy内で同じ入力順・event policy・再構成イベント数。100kは探索に使用したファイルを丸ごと除いた別リストであり、全run/event ID照合で探索との重複0、各リスト内の重複0を確認した。**同一run内の独立eventであり、独立runでの検証ではない。**

入力台帳・catalog query・実際のROOT entriesは [13p5 input provenance](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/provenance/input/README.md)、[3p85 input provenance](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/provenance/input/README.md) に保存した。catalogのevents欄ではなく、実際のTTree entries / GetEntryと解析完了数で入力件数を保証した。

13p5は**13.5 GeV beam-energy FXT**であり、collider sqrt(sNN)=13.5 GeVではない。両energyとも既存 `lambda_imp5` / `fxtmult` 方針を維持し、明示的なVz/Vr/RefMult/VPD cutを新たに加えていない。mode別中心(-0.4,-2)cmはこのprofileではQA用であり、vertex acceptance cutとは異なる。

重要な既存制約: `StRefMultCorr/Param.h` の拡張済みindex0のrun範囲19151029–22179999が両runに一致し、`setParameterIndex` は最初の一致を採用する。他energy用indexを覆う既存の挙動を引き継いでいる。今回は較正・event selectionを変更せずenergy内の比較を固定したため、**正しく較正されたcentrality百分率や他runへの普遍性を主張しない**。別途検証が必要。

## 4. Search procedure and additional upstream cuts

[事前計画](../plans/plan_lambda_kf_upstream_purity_20260928.md)と[探索設定](../../../config/qa/lambda_kf_upstream_purity_20260928.yaml)に沿い、各energyで16種類の実再構成を行い、各出力に3,256通りの追加cutを適用した。合計**52,096点 / energy**。full Cartesian productではなく、single / pairとdecay-length significanceを含むtripleの事前定義gridであり、普遍的な最適解ではなくこの探索範囲でのbestである。

今回PicoDstから初めて比較できた上流条件は、`minNHitsFit`、`minNHitsRatio`、original daughterのInterface / Finder primary threshold、`finderLCut`、`finderLdL2D`、`finderChi2Ndf2D`、`topoChi2NdfCut`。保存済みcandidate treeだけのpost-selectionでは失われた候補を含むこれらの比較はできなかった。

[qhu1 / analysis-note比較](kfparticle_lambda_qhu1_cut_comparison_20260908.md)で優先項目としたdaughter primary偏差・飛行距離・有意度を含めた。13p5のbestにはqhu1側に近いprimary10、Finder L5、L/sigma5の組み合わせが選ばれた。ただし、再フィット後daughterのprimary χ²、TOF、acceptance等までqhu1と同一にしたものではなく、Fig.9の再現ではない。

順位付けは3種類のsidebandに対する最小信号保持率が50%以上の点の中で、背景のPoisson上限を用いた最も保守的なpurityスコアを最大化した。背景0を誤差0のpurity100%として優遇しない。同点は変更parameter数、信号数、固定point IDで決めた。70%・80%保持候補とnominal purityだけのbestも参考として保存した。

両energyのbest50を**100kのmass分布を見る前に固定**し、同じ10kを全YAML条件で再構成した。娘粒子ペアの多重集合、202質量bin（under/overflowを含む）、各対応ペアのmassが完全一致した。13p5は873、3p85は144のselected候補が一致し、mass差の最大値はともに0だった。その後100kを処理し、そこでcutを選び直していない。

## 5. Frozen cut settings

`Default` は今回の基準である **KF Imp5**。Helixとの再比較ではない。表はparameter設定値であり、Finderのstrict不等号とMakerのinclusive不等号の違いを保持している。

| Parameter | Unit / meaning | Default | Optimized AuAu13p5 | Optimized AuAu3p85 |
|---|---|---:|---:|---:|
| nSigmaProton | absolute TPC nSigma maximum | 3.0 | 3.0 | 3.0 |
| nSigmaPion | absolute TPC nSigma maximum | 3.0 | 3.0 | 3.0 |
| imp5MinDCAProton | cm, original daughter gDCA minimum | 0.7 | 0.7 | 1.0 |
| imp5MinDCAPion | cm, original daughter gDCA minimum | 1.0 | 1.0 | 2.5 |
| maxDaughterDistance | cm, final KF daughter distance maximum | 0.5 | 0.5 | 0.5 |
| maxDistanceToPv | cm, final KF parent-PV distance maximum | 0.5 | 0.4 | 0.5 |
| minNHitsFit | daughter TPC fitted hits minimum | 15 | 20 | 25 |
| minNHitsRatio | applied when nHitsMax > 0 | 0.52 | 0.52 | 0.52 |
| minCosPointing | final pointing cosine minimum | 0.998 | 0.998 | 0.998 |
| imp5MaxPathLength | cm, original helix absolute path maximum | 100 | 100 | 100 |
| interfaceChiPrimaryCut | original daughter primary threshold | 3 | 10 | 3 |
| finderChiPrimary2D | transported daughter primary threshold | 3 | 10 | 3 |
| finderLCut | cm, Finder vertex-line distance minimum | 1 | 5 | 1 |
| finderLdL2D | Finder vertex-line significance minimum | 3 | 5 | 3 |
| finderMaxDaughterDistance | cm, Finder daughter distance maximum | 1.5 | 1.5 | 1.5 |
| finderChi2Ndf2D | Finder geometric chi2/NDF maximum | 10 | 10 | 10 |
| topoChi2NdfCut | upstream PV-topology chi2/NDF maximum | 3 | 3 | 3 |
| minDecayLength | cm, final PV-constrained decay length minimum | Disabled | 10 | Disabled |
| minDecayLengthSignificance | final PV-constrained significance minimum | Disabled | 12 | 10 |
| maxChi2Ndf | additional final geometric chi2/NDF maximum | Disabled | Disabled | Disabled |
| maxTopoChi2Ndf | additional final topology chi2/NDF maximum | Disabled | Disabled | Disabled |
| minVertexLineSignificance | additional final vertex-line significance minimum | Disabled | Disabled | Disabled |
| maxMassError | additional final mass-error maximum | Disabled | Disabled | Disabled |
| minMass / maxMass | GeV/c2, stored selected mass range | 1.05 / 1.25 | 1.05 / 1.25 | 1.05 / 1.25 |
| applyLambdaGeometryCuts | upstream Finder geometry gates | true | true | true |
| useTof / strictTofPid / cleanKaonsWithTof | TOF gates | false | false | false |
| useHftTracksOnly | HFT requirement | false | false | false |
| reconstructAntiLambda | anti-Lambda reconstruction | false | false | false |
| rejectBadCovariance | covariance safety gate | true | true | true |
| maxPositionVariance | cm2 | 100 | 100 | 100 |
| maxMomentumVariance | (GeV/c)2 | 1 | 1 | 1 |

`backend: finder_topo`、`selectionProfile: lambda_imp5`、`pidProfile: pico_nsigma`は共通。`maxDistanceToPv`が以前の比較の`maxDCAV0`に対応するKF側の量。Interface / Finderのprimary偏差、Finderのvertex-line有意度、final PV-constrained decay-length significanceは異なる量なので相互置換しない。`Disabled`はoptional final boundが負値という意味で、上流のFinder / Topo条件まで無効という意味ではない。

| Energy | Default mainconf | Optimized mainconf | Optimized KF YAML |
|---|---|---|---|
| AuAu13p5 | [YAML](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/default_validation100k_main.yaml) | [YAML](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/best50_validation100k_main.yaml) | [YAML](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/best50_kf.yaml) |
| AuAu3p85 | [YAML](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/default_validation100k_main.yaml) | [YAML](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/best50_validation100k_main.yaml) | [YAML](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/best50_kf.yaml) |

標準設定へ戻す操作は不要で、既存defaultはそのまま残っている。新しい最適条件は上記mainconfを明示した場合だけ使われる。`.current_mainconf`も変更していない。

## 6. Fixed definitions and numerical results

質量単位はGeV/c²。固定signal windowは `[1.110,1.122)`。default sidebandsは `[1.098,1.106)` / `[1.126,1.134)`。`B=0.75*(Nleft+Nright)`、`S=Nwin-B`、`purity=S/Nwin`、`signalRetention=S/Sdefault`。

### Independent 100k validation

| Energy | Selection | Nwin | B | S | Purity | S / Sdefault | Nwin / NwinDefault | S / B | S / sqrt(Nwin) |
|---|---|---:|---:|---:|---:|---:|---:|---:|---:|
| AuAu13p5 | Default | 13,357 | 1,279.50 | 12,077.50 | 90.4208% | 100% | 100% | 9.439 | 104.50 |
| AuAu13p5 | Optimized | 7,029 | 93.75 | 6,935.25 | 98.6662% | 57.4229% | 52.6241% | 73.976 | 82.72 |
| AuAu3p85 | Default | 2,269 | 415.50 | 1,853.50 | 81.6880% | 100% | 100% | 4.461 | 38.91 |
| AuAu3p85 | Optimized | 1,024 | 23.25 | 1,000.75 | 97.7295% | 53.9924% | 45.1300% | 43.043 | 31.27 |

purity上昇は13p5で8.25±0.25 percentage points、3p85で16.04±0.85 percentage points（paired bootstrap）。一方、信号を減らすため単純な`S/sqrt(Nwin)`は低下している。今回の目的はpurity優先であり、統計的有意度最大化とは異なる。

### Discovery 10k and alternative retention choices

| Energy | Choice | Nwin | B | S | Purity | Signal retention | Validated on 100k |
|---|---|---:|---:|---:|---:|---:|---|
| AuAu13p5 | Default | 1,370 | 143.25 | 1,226.75 | 89.54% | 100% | Yes |
| AuAu13p5 | best50 | 730 | 6.75 | 723.25 | 99.08% | 58.96% | Yes |
| AuAu13p5 | best70 | 898 | 9.75 | 888.25 | 98.91% | 72.41% | No |
| AuAu13p5 | best80 | 1,010 | 18.00 | 992.00 | 98.22% | 80.86% | No |
| AuAu13p5 | nominal50 | 678 | 4.50 | 673.50 | 99.34% | 54.90% | No |
| AuAu3p85 | Default | 269 | 40.50 | 228.50 | 84.94% | 100% | Yes |
| AuAu3p85 | best50 | 123 | 0.75 | 122.25 | 99.39% | 53.50% | Yes |
| AuAu3p85 | best70 | 173 | 2.25 | 170.75 | 98.70% | 74.73% | No |
| AuAu3p85 | best80 | 192 | 4.50 | 187.50 | 97.66% | 82.06% | No |
| AuAu3p85 | nominal50 | 150 | 0.75 | 149.25 | 99.50% | 65.32% | No |

`best50`はnominal purityだけではなく、sidebandを変えたときの弱い条件も考慮した順位。`nominal50`の見かけのpurityが高くても採用しなかった。70% / 80%候補は統計重視の場合の参考であり、独立100kで検証済みとは扱わない。詳細cutは各 `scan/chosen_points.json` に保存した。

## 7. Background sensitivity and limitations

near sidebandsは `[1.100,1.108)` / `[1.124,1.132)`、farは `[1.094,1.102)` / `[1.130,1.138)`。追加で各sidebandに非負Bernstein degree0/1/2のPoisson fitを行い、signal window下の背景を外挿した。これはcut再選択のためには使わず、背景モデルの感度確認に限定した。

| Validation sample | Selection | Three-sideband purity range | Nine-model purity range |
|---|---|---:|---:|
| AuAu13p5 | Default | 90.03–91.02% | 88.35–91.02% |
| AuAu13p5 | Optimized | 98.48–98.79% | 98.14–98.79% |
| AuAu3p85 | Default | 81.09–82.55% | 80.97–82.70% |
| AuAu3p85 | Optimized | 97.44–98.02% | 96.24–98.06% |

9モデルはすべて数値収束検査を通過した。3p85の一部にはshape parameter boundaryの警告があり、範囲は信頼区間ではない。3 sidebandすべてにおける最小信号保持率は13p5が57.12%、3p85が53.59%だった。

背景診断CSV: [13p5](../../../share/figure/auau13p5_Lambda_KFParticle_purity_upstream_20260928/validation_background_checks_background.csv)、[3p85](../../../share/figure/auau3p85_Lambda_KFParticle_purity_upstream_20260928/validation_background_checks_background.csv)。同じディレクトリにモデル別の図・ROOTも保存した。

注意点:

- 10k探索にはcut選択による楽観バイアスがある。最終判断は独立100kを優先し、同じ100k上で再調整していない。
- purityはmass peakのsideband推定であり、MC truth、primary-Lambda fraction、feeddown除去を測ったものではない。peaking背景・signal tail・cutによるmass shape変化を完全には排除できない。
- L / DCA / hits等の強いcutはpT・rapidity・崩壊位置の受容を変える。検出効率、物理収量、偏極等への補正は別途MC / embeddingが必要。
- 同一run内の結果であり、他run・他centrality・異なるproductionへの適用は追加検証が必要。既存centrality較正の制約は§3のとおり。
- TOF PID、parent-refitted daughter primary χ²、HFT、pT/eta/centrality最適化、covariance安全条件の緩和、mass constraintは今回試していない。TOFは既存Imp5 profileでは禁止され、再フィット後daughter χ²は現在未実装なので、実施済みと混同しない。

## 8. Verification, preservation, and new files

- 新しい13p5 default10kは前回 `default_regression.root` と68個のTH1/TH2/TH3、2,536,396セル（flow / error / axisを含む）、42,871行×42 leafが完全一致した。
- 32本の探索run、2本の最適10k再実行、4本の100k runが完了。要求入力数への暗黙clampは許さず、全runのinput / completed / event policy / 出力構造を確認した。
- 最適10kは候補ペア・全202 mass bin・各pairのmassが完全一致。source別の独立ROOT選択照合も合格した。
- 最終100kは別実装のstdlib numeric auditorで実candidate CSVから再計算し、Nwin/B/S/purity/保持率・全bin・正規化・event非重複・凍結設定とhashを照合した。
- Pythonテストはpurity関連32件、frozen replay8件、freeze jobs3件、独立numeric auditor2件が合格。ROOT5コンパイル、合成図、背景ゼロ時の正の上限、上書き拒否も確認した。
- `provenance/before_SHA256SUMS` にある既存mainconf / cut / analysis info / Maker / Interface / 共有ライブラリ / `.current_mainconf` は変更なし。
- 新規補助コードは `common/macro/*LambdaKf*{Study,Purity,Comparison,Replay,Background}*` 等、試験設定は `config/qa/lambda_kf_upstream_purity_20260928.yaml`。各studyの `provenance/code_snapshot` に実行コードを保存し、`run_inventory.csv` に19本ずつの再構成出力とROOT checksumを記録した。コードsnapshotは完全な独立STAR実行環境ではなく、実行版を追跡するための保存。

| Artifact | AuAu13p5 | AuAu3p85 |
|---|---|---|
| Study README | [README](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/README.md) | [README](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/README.md) |
| Frozen decision | [JSON](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/frozen_decision.json) | [JSON](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/frozen_decision.json) |
| Final metrics | [CSV](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/validation100k/comparison/chosen_metrics.csv) | [CSV](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/validation100k/comparison/chosen_metrics.csv) |
| Paired bootstrap | [CSV](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/validation100k/comparison/paired_event_bootstrap.csv) | [CSV](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/validation100k/comparison/paired_event_bootstrap.csv) |
| Final independent audit | [JSON](../../../rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/validation100k/final_numeric_audit.json) | [JSON](../../../rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/validation100k/final_numeric_audit.json) |

初回の補助設定ではConfigManagerが絶対参照へ `config/` を付ける問題があり、`*_main_resolved.yaml` を新規作成して解決した。初回CINT exportのcast解釈エラーは既存exporterをACLiCでコンパイルして解決。独立ROOT checkerに必要なsource別native baselineも新規manifestに追加した。失敗ログ・途中成果は監査用に残し、最終集計には使用していない。これらは再構成アルゴリズムや既存解析の変更ではない。

## 9. References

- [前日のpurity検討](lambda_kf_purity_study_20260927.md)
- [前日の70%保持率を含む追加検討](lambda_kf_purity70_refinement_20260927.md)
- [qhu1 / analysis-note cut比較](kfparticle_lambda_qhu1_cut_comparison_20260908.md)
- [今回の事前計画](../plans/plan_lambda_kf_upstream_purity_20260928.md)
