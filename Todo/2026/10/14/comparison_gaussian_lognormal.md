# Gaussian／log-normal比較のまとめ

## 今回行ったこと

3.85 GeV・2021 P24iyの高統計統合ROOTに保存されたALLとTOF3タグを使用。
R=0.30〜2.50 GeV/cをnative幅0.01で分け、各応答modelについて
TAGd／TAGt／TAG³Heの計660スペクトルを、linear／exponential背景それぞれでfitした。
各モデルの入力bin・保存count・fit窓は同一で、比較対象として照合した。
未収束・境界・zero-yield・軸外・未確認の核種帰属は表と図へ残した。

ユーザー指定どおりGaussian／log-normalの関数モデルに絞り、実測TOF形状テンプレートは追加していない。
原ROOT・旧fit・既存PDF・本番PID／解析コードは今回の作業で変更していない。

## 正式な図と数値

出力基点は`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/`。

- [中心・幅・model適合の比較（4ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/comparison_gaussian_lognormal_004/response_parameters.pdf)。
- [同一データへの重ね描画（24ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/comparison_gaussian_lognormal_004/saved_spectra_overlays.pdf)。
  R=0.405／0.605／0.825を、全dE/dx域0〜120とピーク付近、linear／logで対照表示。
- [全1760行の数値・警告](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/comparison_gaussian_lognormal_004/comparison.tsv)。
- [ALL／TOFタグの生COLZと中心・幅ガイド（5ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/nominal_gaussian_guides_004/nominal_gaussian_colz_guides.pdf)。
  失敗・未支持点は赤色で、連続した検証済み曲線として接続していない。
- [GaussianのR=0.80〜0.90の残差（30ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/tag_native_linear_core_views_002/tag_fit_residuals.pdf)。
- [log-normalの同じR域の残差（30ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/tag_lognormal_native_smoke_linear_views_002/tag_fit_residuals.pdf)。

紙面外へ軸や注記が出るROOT5の表示問題を修正した別名の正式版を案内する。
旧表示版・コンパイル失敗時のログ・source snapshotも保持している。数値は描画修正で変更していない。

## 全域での比較結果

Poisson deviance（観測分布とmodelのずれ）について、両モデルが数値収束した区間だけで集約した比。
小さい方が当該保存点でデータへ近い。ただし大域的な最小値・検証合格・有意度を意味しない。

| TOFタグ | linear背景でlog-normal／Gaussian | exponential背景でlog-normal／Gaussian |
|---|---:|---:|
| d／⁴He互換TAGd | 1.763 | 1.859 |
| TAGt | 2.124 | 2.456 |
| TAG³He | 1.705 | 1.716 |

全域でlog-normalへ一括切り替える結果ではなかった。
TAGdの指標は同じd＋名目⁴Heのfitであり、dと⁴Heの独立測定として二重加算しない。

一方、R=0.80〜0.90のTAGtはlog-normalでdevianceがGaussianの約0.28倍に改善した。
R=0.825の例:

| スペクトル | Gaussianのχ²/ndf | log-normalのχ²/ndf |
|---|---:|---:|
| TAGt | 93.62 | 15.90 |
| TAG³He | 438.21 | 771.39 |
| TAGd | 493.76 | 799.79 |

tの局所改善を全Rへ一般化しない。改善例も、現在の応答モデルが分布を十分説明できたとは判断していない。

## なぜ関数・purity・本番PIDをまだ確定しないか

### 高統計域の分布形

R=0.80〜0.90のGaussian fitで、ピーク±2幅の観測とmodelの絶対差を計数で割った量は
d約7%、t約6%、³He約12%、名目⁴He約6%。
1 keV/cmへ再集約しても約5〜11%残るため、単にnative binが細かいことだけでは説明できない。
粗いR幅0.05の平均化だけでもなく、R幅0.01でも応答形の不一致が残る。

中心点を滑らかにするだけで、この不一致や信号／背景の分解誤差はなくならない。
現行と同じlog(R)中心3次／幅2次、連続接続した区分3次を試験するコード・候補は保存したが、
全4核種で検証済み関数を得た状態ではない。

### 低Rの⁴He帰属

同じTAGd、同じ入力・窓でも背景modelを変えると、名目⁴Heの結果が大きく変わる。

| R（GeV/c） | linearの中心／幅／窓収量 | exponentialの中心／幅／窓収量 |
|---|---|---|
| 0.405 | 86.64／5.33 keV/cm／51.5 | 71.22／13.12 keV/cm／256.5 |
| 0.455 | 79.59／7.35 keV/cm／95.6 | 71.53／11.99 keV/cm／198.2 |

χ²が約1、数値statusが正常でも、これだけでは⁴Heピークの同定を保証しない。
zero-yieldや、他の帯を幅広い名目⁴Heへ割り当てた例もある。
**R=0.4〜0.5を信頼できる⁴He校正下限として認定する根拠は、今回の比較では得られていない。**
これは低Rの⁴Heが存在しないという意味ではない。TOF受理・保存軸・統計・成分識別の制約を含む。
各核種の信頼できる下限・上限の確定は、残差と帰属の検証が済むまで未完了とする。

### パラメータの意味

- Gaussian: 中心は平均、幅は標準偏差。
- log-normal: 中心はmode、幅は局所曲率から定義した量。平均／標準偏差ではない。
- 両方ともnative R binの射影応答の比較で、intrinsic幅の検証結果ではない。
- log-normalのmode±2曲率幅の積分確率に、Gaussianの95.45%をそのまま代入しない。
- TAGはALLの部分集合。タグを独立の追加データとして尤度へ足さない。
  TAG→ALLの転用と核種帰属も未検証のため、確定purityを報告しない。

## 完了・未完了

完了: 入力監査、両関数モデル・両背景の全域fit、モデル比較、全軸表示・COLZガイド、数値表、再現コード。

未完了: 信頼できる全4核種の関数・支持域、そのshapeを用いたALL purityの検証、最終校正YAML・本番PIDの更新。
ALLのrunnerはROOT5数値検証と入出力準備試験までで、今回の未検証候補を確定値として適用していない。

次の検討候補は、核種ごとの応答をGaussian 2成分などの関数で記述する比較。
それを試すか、現状の単一応答modelによる参考値の範囲でまとめるかは、今回の結果を示して相談する。
実測TOF形状テンプレートを許可されたものとして追加することはしない。
