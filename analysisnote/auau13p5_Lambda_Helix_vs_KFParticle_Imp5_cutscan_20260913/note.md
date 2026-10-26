# AuAu13p5 Λ：Imp5 / KFParticle、候補数整合cut scan

更新日：2026-09-13。状態：保存済み10,000入力イベントの探索的QA完了。本番cutの決定、絶対効率測定、公開用の有意差評価ではない。

## 1. 目的と結論

ユーザー指定によりqhu1比較を中断し、旧Helix Imp5とKF Imp5のΛ近傍の**背景込み候補数**をそろえた場合にもS/Bの向上が残るかを調べた。これは絶対効率の代用と主張せず、指定された比較指標として扱う。

KF追加最終cutのPV距離上限を0.8 cm、daughter距離上限を0.6 cm、cos pointing下限を0.998とすると、主窓候補数は1759（旧1775の99.1%）。sideband推定S/Bは5.22 ± 0.35で、旧5.20 ± 0.35とほぼ同じ。今回調べた範囲では、候補数を戻した状態で明確なS/B改善を確認できなかった。

[結果要約・Notion Imp5 cut照合表](../../mdfiles/kfparticle/comparisons/lambda_helix_kf_imp5_cut_scan_20260913.md) を参照。

## 2. データ・選別・出典

入力はAuAu13p5の前回の先頭10,000イベント。同一性は前回の入力ファイル内容・順序・ループ範囲の監査に基づく。今回のread counter一致だけでイベントIDの同一性を証明したとはしない。

- [Helix元ROOT](../../rootfile/auau13p5_anaLambda/local_compare_10000_20260907.root)：read=10000、reconstructed=9744、図示全質量範囲4914候補。
- [KF Imp5元ROOT](../../rootfile/auau13p5_anaLambda_KFParticle/Imp5/local_compare_10000_20260908.root)：同じread/reconstructed、全質量範囲3135 selected、valid raw tree 42871候補。
- [前回の入力監査・実行記録](../../mdfiles/kfparticle/comparisons/lambda_helix_kf_imp5_comparison_10000_20260908.md)。そこにある旧KFディレクトリ名は後日整理されており、現在の実体は上のリンク。
- [前回の入力・実行provenance](../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5/comparison_provenance_20260908.tar.gz)。入力ID一覧SHA256は `aeb0e8c9da2d62e1c3fed9dd7e1897d473af766887a38a1d05539ed5618da4d3`。
- [Helix cut YAML](../../config/maker/maker_auau13p5_anaLambda.yaml)、[KF Imp5 cut YAML](../../config/cuts/kf/kf_auau13p5_anaLambda_KFParticle_Imp5.yaml)。

共通centrality選別（0–80%）を維持。qhu1の20–50%等へ変更していない。Imp5のevent policyは現行旧Makerに合わせたもの：maxNTr=0、追加Vz/Vr/RefMult/VPD cutなし、共通centrality不正256イベントを除外。collider/fixed-targetのmode-aware QA座標設定を今回変更していない。PIDはproton/pionの保存済み |nSigma| ≤ 3だけで、TOF・kaon仮説を使わない。両daughterのnHitsFit ≥ 15、ratio ≥ 0.52、元gDCA ≥ 0.7/1.0 cm、元helix経路長絶対値 ≤ 100 cmを固定した。

KF固有の既存設定を固定：共分散の有限性・妥当性とvariance上限100/1、interfaceChiPrimaryCut=3、Finder daughter距離1.5 cm、L cut=1 cm、chiPrimary2D=3、chi2/NDF=10、L/dL=3等。採用済みupstreamのΛ PV χ²/NDF < 3も変更しない。最終親massはPDG値へのconstraintなし。anti-Λは含まない。

## 3. 保存済み候補による再選別

`StLambdaKFParticleMaker.cxx` はvalid raw候補をtreeへ記録した後でselected判定により選別する。このtree全42871候補を使い、元のpath・mass・その他最終cutを維持し、次の3種類だけを置換した。3つ以外の固定最終cutを通過する候補は21362。

| パラメータ | scan値 |
|---|---|
| maxDistanceToPv [cm] | 0.5, 0.6, 0.75, 0.8, 1, 1.5, 2, 3, 5, 無効化 |
| maxDaughterDistance [cm] | 0.5, 0.6, 0.75, 1, 1.25, 1.5 |
| minCosPointing | 0.998, 0.997, 0.995, 0.99, 0.98, 0.95 |

360 unique条件。負のPV距離設定は**追加Maker cutだけ**の無効化を意味し、upstreamのPV選別を外さない。0.8 cmはNotionの以前のphaseにも現れる比較anchor、それ以外は緩和範囲と飽和を確認する今回のQA選択である。

代表点は、(a) PVのみ、(b) PVとdaughterでpointing固定、(c) 3条件全体のそれぞれについて、旧主窓候補数との絶対差が最小の点。同差なら変更パラメータ数が少ない点、さらに同じならmanifestの順序。**S/Bは代表点選択に使わない。** 同じデータでのcut調整なので、独立検証済みの最適cutとはしない。

## 4. 質量窓と背景の定義

2026-09-13の表示更新により、scanのmass overlay（`paired_cut_overlay`、`count_matched_overlay`、`pv_only_overlay`）は左を無規格化の候補数、右を各曲線の通常bin合計が1になるNormalize表示とした。左右とも1.05 ≤ M < 1.25 GeV/c²、1 MeV/c² bin。右は表示用histogramに `Scale(1/Integral(1,Nbins))` を適用し、bin幅係数を使わずunderflow/overflowを除外する。元Count histogram・cut・質量窓・以下の候補数とS/Bは変更しない。補助fit図は無規格化の候補数表示を維持し、その表示範囲は別途指定するfit範囲に対応する。補助fit図とPV距離scan曲線は今回の表示更新の対象外。

| 用途 | 質量範囲 [GeV/c²] |
|---|---|
| 主窓 W | [1.110, 1.122) |
| 狭い窓 | [1.112, 1.120) |
| 基準sideband左右 | [1.098, 1.106), [1.126, 1.134) |
| 近いsideband左右 | [1.100, 1.108), [1.124, 1.132) |
| 遠いsideband左右 | [1.094, 1.102), [1.130, 1.138) |

これらは今回のQA設定であり、Notionの歴史的fit手順を再現する値ではない。

窓内実測候補数をNwin、左右sideband数をNL/NRとし、局所背景密度を一定と近似する。α = 窓幅 / 左右sideband合計幅、B = α(NL+NR)、S = Nwin−B、S/B = Nwin/B−1。主窓ではα=0.75。負のSを0へ切り上げない。

統計目安は Var(B)=α²(NL+NR)、Var(S)=Nwin+Var(B)、Var(S/B)=Nwin/B² + Nwin²Var(B)/B⁴。窓が非重複の独立Poisson近似であり、1イベントからの複数候補の相関や共有daughterは含まない。両手法は同じイベントを共有しscan間も候補を共有するため、手法差・scan差の誤差や有意性は未評価。

## 5. 結果

| 条件（PV距離, daughter距離, cos） | 全質量範囲N | 主窓N | N/Helix | 基準sideband S/B |
|---|---:|---:|---:|---:|
| Helix Imp5 | 4914 | 1775 | 1 | 5.195 ± 0.349 |
| KF baseline (0.5, 0.5, 0.998) | 3135 | 1370 | 0.7718 | 8.564 ± 0.739 |
| KF PVのみ (3.0, 0.5, 0.998) | 4585 | 1653 | 0.9313 | 5.022 ± 0.348 |
| KF 2距離 (0.8, 0.6, 0.998) | 4838 | 1759 | 0.9910 | 5.221 ± 0.353 |
| KF 3条件 (0.75, 0.6, 0.997) | 5276 | 1783 | 1.0045 | 4.943 ± 0.329 |

基準sidebandでHelixはNL/NR=182/200、B=286.5、S=1488.5。KF baselineは84/107、B=143.25、S=1226.75。2距離緩和では168/209、B=282.75、S=1476.25。従って背景込み候補数だけでなく、この単純背景推定下の信号数も旧に近づくが、S/B改善は残らない。

PV距離だけを0.5, 0.6, 0.75, 0.8, 1, 1.5, 2, 3, 5 cm, 無効へ緩めた主窓数は、順に1370, 1488, 1564, 1590, 1628, 1650, 1651, 1653, 1653, 1653。**PV距離cutだけでは、この他条件を固定した保存候補集合で旧1775には届かない。** 3条件gridで旧候補数の±1%に入るのは上表の1759と1783の2点だけだった（事後の範囲確認であり新たな選択規則ではない）。

### 質量窓・sideband位置のcross-check

狭い窓のNはHelix1627、KF baseline1285、PVのみ1512、2距離1613、3条件1631。2距離緩和の旧比はここでも約99.1%。

| 条件 | 基準sideband S/B | 近いsideband | 遠いsideband |
|---|---:|---:|---:|
| Helix | 5.195 | 5.431 | 5.704 |
| KF baseline | 8.564 | 8.514 | 9.498 |
| KF PVのみ | 5.022 | 4.862 | 5.352 |
| KF 2距離 | 5.221 | 5.172 | 5.663 |
| KF 3条件 | 4.943 | 4.899 | 5.240 |

2距離緩和の明確な改善がないという解釈は、このsideband変動で覆らない。値の変動を厳密な系統誤差として換算・統合はしていない。

### gaus+pol2の補助確認

全5条件に同じgaus+pol2関数族（数値安定化のため中心をずらした2次多項式）をfit。ROOT5 ACLiC、`RSQLI0N`：Poisson likelihood、bin積分。同じ主窓でGaussian積分Sと背景多項式積分Bを計算。fit範囲は[1.090,1.150]と[1.095,1.140] GeV/c²。σ探索範囲は[0.0003,0.008] GeV/c²、最低fit範囲entries=50。

| 条件 | S/B：広いfit範囲 | S/B：狭いfit範囲 | χ²/NDF：広い / 狭い |
|---|---:|---:|---|
| Helix | 4.660 ± 0.286 | 4.342 ± 0.315 | 126.90/54 / 98.66/39 |
| KF baseline | 6.622 ± 0.466 | 6.461 ± 0.561 | 98.47/54 / 85.56/39 |
| KF PVのみ | 3.890 ± 0.229 | 3.754 ± 0.262 | 85.93/54 / 72.25/39 |
| KF 2距離 | 4.127 ± 0.237 | 3.920 ± 0.267 | 77.51/54 / 66.51/39 |
| KF 3条件 | 3.925 ± 0.223 | 3.784 ± 0.256 | 79.99/54 / 62.36/39 |

全10 fitでstatus=0、covariance quality=3、正背景・境界非接触等の数値検査を通過した。しかし**数値的成功はモデルの適合性を保証しない**。特にHelixのχ²/NDFは良好とはいえず、単純gaus+pol2のモデル依存がある。誤差は固定モデル・fit範囲での全共分散（S–B相関を含む）だけ。fitのS+Bを窓内実測数へ拘束していないため、Nwinと完全一致しない。

この補助fitでも候補数整合後にKFのS/Bが旧を上回らないが、「KFが本質的に劣る」「差が有意」とは結論しない。sidebandではほぼ同じ、fitでは低いという背景推定依存を残す。[fit CSV](../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/fit_mass.csv)、[fit図](../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/fit_mass_range00.png)。

## 6. 検証・未確定点

- 元の保存出力QA成功。再評価baselineの全候補selected判定、元bin幅/共通bin幅のmass histogram内容を完全再現した。
- 元Helix/KF histogramから独立に直接積分し、主窓1775/1370を確認。
- ROOT5 ACLiCで360条件scan成功、CSVは2窓×3sideband×(360+Helix)=2166行。初回と最終版のCSV完全一致。
- 今回のNormalize表示への更新前には、表示だけを修正した最終scanと初回scan、fit入力の5 histogramを読み直し、flow binを含む1010 cellsの内容・誤差・bin下端を比較、差分0。
- 更新前のoverlayの軸・凡例・無規格化の重ね描きを目視確認済み。fitは上記モデル品質の注意付き。
- YAML generatorは未知/重複キー、不正数値、不整合の窓等を拒否し、manifestと成果物の上書きを拒否。QA設定・入力ROOT・QAソースのSHA256を保存。
- ROOT5 CINTでscanを直接解釈する方法はSTL dictionaryの問題で失敗したため、検証済み経路は `.C+` によるACLiC。ROOT6導入や既存解析の再buildは不要だった。

保存前のPID・共分散・Finder/Topoの不採用候補は今回のtreeから復元不能。運動学/centrality依存、候補共有による相関、より十分な背景/ピークモデル、独立データでのcut検証は未評価。MC/embedding効率は今回ユーザー指定で対象外。[Open Questions](open-questions.md) に意図的な省略と次の検証を記録した。今回の比較結果を公開用の物理測定や全KF機能の最適化完了とはしない。

## 7. 再現性・変更範囲

元ROOT SHA256（今回も不変）：

```text
1d8bf45a0446ce7ac7353d0f55d7e6386a1f000ab6adb8f65dba30a16ba580b6  Helix
925a856d05692889a9d2920ec2e4202ae19af52ce9245fa5f3aea9d657f06b67  KF Imp5
```

HEADは `1c79dfc78cdb587c1ff5b32cfe1a6ebd8a921ed5`、既存dirty worktreeを保持。今回の新規QA：

- [scanLambdaKfImp5Cuts.C](../../common/macro/scanLambdaKfImp5Cuts.C)
- [prepareLambdaKfImp5Scan.py](../../common/macro/prepareLambdaKfImp5Scan.py)
- [fitLambdaKfScanMass.C](../../common/macro/fitLambdaKfScanMass.C)
- [QA YAML](../../config/qa/lambda_kf_imp5_cut_scan_20260913.yaml)

既存のMaker/StLambda、Helix/KF本番YAML、元ROOT、`.current_mainconf` を変更していない。farm job、PicoDst再処理、git commit/pushなし。`config/qa` は既存 `.gitignore` の `config/*` で無視されるため、将来明示的にcommitする場合はこのQA YAMLの含め方に注意する（今回ignore/indexも変更していない）。

検証環境：既存Singularity STAR image、SL24y、ROOT 5.34/38、GCC 4.8.5。host manifest生成はPython3/PyYAML。環境wrapperは既存のLD_LIBRARY_PATHを保持する。setupによるmainconf切替を行わず、ACLiC生成物は/tmpへ分離。

プロジェクトrootから、未使用のmanifest/出力名を指定する：

```bash
python3 common/macro/prepareLambdaKfImp5Scan.py \
  config/qa/lambda_kf_imp5_cut_scan_20260913.yaml /tmp/new_scan_manifest.tsv

# 既存STAR/ROOT5環境内。ACLiC生成先を作業用へ設定するdriver使用を推奨。
root4star -b -q 'common/macro/scanLambdaKfImp5Cuts.C+("/tmp/new_scan_manifest.tsv","new-scan-output")'

root4star -b -q 'common/macro/fitLambdaKfScanMass.C+("new-scan-output/scan.root","new-fit-prefix","HelixLambdaCounts,KfBaselineCounts,KfNearestPvCounts,KfNearestPvDaughterCounts,KfNearestCounts","1.09:1.15,1.095:1.14",1.110,1.122,0.0003,0.008,50)'
```

実際のwrapper、ACLiC driver、最終manifest、baseline/scan/fit/一致検証ログ、ソース・設定・SHA256を [元provenance archive](../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/comparison_provenance_20260913.tar.gz) に保存。これはNormalize表示への更新前の履歴として保持する。表示更新分の再現用記録は [provenance_normalized_20260913.tar.gz](../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/provenance_normalized_20260913.tar.gz)、更新前の図・成果物の復元先は [before_normalized_overlays_20260913.tar.gz](../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5_cutscan_20260913/before_normalized_overlays_20260913.tar.gz)。wrapper内の環境pathは実行時のもの。再実行時は新しい作業先を選ぶ。

成果物は [KF出力索引](../../rootfile/auau13p5_anaLambda_KFParticle/mdfiles/cut_scan_20260913.md)、図は [assets索引](assets/figures.md) を参照。Notionの取得JSONは認証不要の公開ページの一時作業領域に保持し、解析データとともに外部へ送信したりgitへ追加したりしていない。
