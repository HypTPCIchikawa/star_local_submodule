# Task 3: ⁴HeのTPC純度評価 — 実装計画

## 目的と入力

既存のKFParticle Λ–⁴He解析の `true/hDedxP` を細かい剛性 `R=|p/q|` 区間に分け、dE/dx射影をfitして⁴Heの校正nSigma窓内の `S/(S+B)` を推定する。対象は3.85 GeVの2021 P24iy統合ROOTと13.5 GeVの2020 P24iy統合ROOT。原本ROOTは読み取り専用とし、新しいイベント解析・job投入は行わない。

- 3.85 GeV: `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root`
- 13.5 GeV: `rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/auau13p5_anaFemtoLambda_4He_KFParticle_2020P24iy_partial_merge.root`

## 確認した既存実装

`FemtoLambdaLegacy::CollectNuclei` の `hDedxP` は核種PID・TOF・最終剛性上限より前にfillされ、`nHitsDedx >= 15`、`gPt >= 0.1 GeV/c`、正のdE/dxを満たすtrackを含む。正負電荷を区別せず、選択済みΛの有無とも独立である。`hDedxP_4He` は既に±2のTPC窓を適用した分布なので、背景推定用fitには使用しない。

保存TH2は剛性600 bin / 0–6 GeV/c（幅0.01）、dE/dx 1200 bin / 0–120 keV/cm（幅0.1）。既存最終選択の実際の上限は `maxPOverQ: 2.5`（Todoの説明の2.0とは異なる）。核種間の最小絶対nSigmaによるbest-species選択もある。

## 実装手順

1. `tools/he4_purity/` にROOT5用fit macro、YAML設定、mainconfを引数にする実行wrapperを追加する。
2. `StMaker/common/NuclearIdDeDxVsMom.h` の⁴Heおよび近傍核種の校正関数を直接使用し、係数を別実装で複製しない。YAMLで剛性区間、fit範囲、許容fit状態、背景候補、plot設定を管理する。
3. 剛性0.40–2.50 GeV/cを保存TH2の最小幅0.01 GeV/c刻みで射影する。0.05幅では低剛性側のローカスの変化が校正幅と同程度になるため、fitは最も細かい保存binで行う。各区間のfit座標は中央剛性での校正nSigmaへ変換し、fitするのは元の未切断dE/dx射影とする。剛性bin内のローカス移動・校正gate端の変化は監査に残す。完全に連続なtrack単位のgateは保存済みTH2から復元できないので、これはbin中央の校正gateで評価した推定である。
4. Gaussianの⁴He信号に、pol1相当の正値線形背景とexponential背景をそれぞれ加えてfitする。³Heなど近傍核種のローカスがfit範囲に入る場合は、そのGaussianを背景側に加え、近傍ピークを⁴He信号と取り違えない。実分布・残差・fit状態・モデル間差を比較する。
5. Purityの積分窓はfit Gaussianの±2σではなく、解析が使う校正関数の `mean(R) ± 2 * sigma(R)` とする。fit幅・fit平均は独立に記録し、背景modelやGaussianの破綻した区間には有効な純度を付けない。
6. fit covarianceから純度の誤差を伝播し、背景モデルの違いは別のモデル差として残す。低統計・境界・fit失敗・大きな残差等はCSVとplotで識別する。
7. 各入力ROOTのフォルダ内 `pdf/he4_purity_20261013/` に、剛性対純度のPDF、全射影fitをまとめたPDF、CSV、JSON、設定・ソースsnapshot・実行ログを保存する。
8. 完了報告を `Todo/2026/10/13/task3_he4_purity_report.md` に日本語で作り、入力・結果・fitの限界・出力リンクをまとめる。

## 解釈の限界

この評価は中央度を合算したtrackのTPC ±2校正窓の純度である。保存TH2に電荷・TOF・核種best判定の追加軸がないため、正電荷⁴Heだけの純度、best-species後の厳密な純度、新規R≥2のTOF併用選択後の純度は既存ファイルから復元できない。Task 2で追加するQAヒストにより将来の解析で確認する。

dE/dxの幅には非Gaussian tailがあり得る。単一Gaussianと簡単な背景で説明できない区間は、その事実を報告し、目標純度の数値を無理に作らない。保存TH2Fの大きなbin値には浮動小数点加算の丸めもあるが、ここでは原本を修正しない。

## 検証

- ROOT5 / SL24yのSingularity wrapperでコンパイル・実行する。
- 元TH2の軸・射影sum・integral窓を監査JSONに記録する。
- fit結果のvalidity、covariance、境界、非負の期待bin数、残差を確認する。
- 入力ROOTのsize・mtime・inodeが処理前後で変化しないことを確認する。
- PDFを目視し、⁴Heのピーク同定とpurity plotの表示を確認する。

## 実装後の到達点

両エネルギーで210区間×2モデルを処理し、主表示PDF、両背景モデルの比較・差、全射影fit、CSV・JSONを保存した。13.5 GeVの168区間は今回のモデル品質基準を満たすが、候補数coverageは34.94%。3.85 GeVの高統計では全区間で単純Gaussianモデルの残差が基準を満たさず、収束した175区間を参考推定として明確に区別した。品質基準を緩めて数値を採用することはしていない。結果・限界・出力リンクは [完了報告](task3_he4_purity_report.md) を参照。
