# 現在のピーク応答・fit範囲比較のまとめ（2026-10-20）

## 結論と今回の区切り

Git元PIDから仕切り直し、Gaussian・指数tail・左右非対称Gaussian、保存精度とfit範囲の影響を比較した。
**新しい応答モデルをさらに増やさず、ここまでの結果を先にまとめる**というユーザー方針に従う。

- 本番6入口の中心・幅はGit元式へ復帰済み。中心の3次式・幅の2次式は変えていない。
- 保存丸めとfit範囲のずれが、大きなχ²の一部に寄与していた。左右非対称応答との組合せで改善する点がある。
- 代表点の³He・⁴Heでは局所的に小さなχ²を得たが、d・tに不一致が残る。
  **新しい係数・純度・cutoffの本番採用はしていない。**

## 同じ代表8スライスの比較

Au+Au 3.85 GeV、2021 P24iyの保存済みtrack QA。native R幅0.01 GeV/c、dE/dx幅0.1 keV/cm。
d／⁴Heは同じm²窓、t／³Heはそれぞれのm²窓を使用。電荷符号・centralityを積分した分布である。
値はχ²/ndfであり、純度ではない。すべて同じraw cacheを使い、観測count・誤差・binは変更していない。

| R [GeV/c] | Species | Gaussian, Git ROI | Gaussian + Float16, Git ROI | Split-normal + Float16, Git ROI | Gaussian + Float16, local ROI | Split-normal + Float16, local ROI |
|---|---|---:|---:|---:|---:|---:|
| 0.875 | d | 3190.343 | 3152.430 | 531.537 | 305.765 | 19.574 |
| 0.875 | t | 106.762 | 95.056 | 40.289 | 24.175 | 8.856 |
| 0.875 | ³He | 391.914 | 362.861 | 36.143 | 7.804 | 1.942 |
| 0.875 | ⁴He | 23.987 | 21.088 | 3.746 | 2.271 | 1.050 |
| 1.775 | d | 284.932 | 234.801 | 211.803 | 6.667 | 3.818 |
| 1.775 | t | 196.309 | 197.715 | 34.881 | 3.955 | 0.594 |
| 1.775 | ³He | 18.696 | 11.376 | 6.379 | 1.071 | 0.851 |
| 1.775 | ⁴He | 1.926 | 1.305 | 1.254 | 1.074 | 1.043 |

- **Git ROI**：元式の中心±2.5元σをnative bin端に揃えたfit範囲。
- **Local ROI**：元式付近でraw 3bin移動和の最大を探し、そのbin中心±1.5元σとした範囲。
  移動和は範囲を探すためだけに使い、fit dataは平滑化しない。初期値・nodesはGitのまま。
- **Float16**：信号と背景の予測積分へ保存丸めを含める。観測値の補正・再binではない。
- **Split-normal**：1つのピークで左右Gaussian幅を分ける。別核種のfreeピークを追加したものではない。
  中心はmode、σは左右幅の幾何平均。全RMSやGaussianの標準偏差と同じとは扱わない。

local ROIではbin数・自由度も変わるので、χ²/ndfの低下だけでモデル選択や物理的な採用を決めない。
対象ピーク＋局所線形背景の診断であり、全核種を分離した純度fitではない。
特に、背景が信号の肩を吸収する可能性、中心・幅と背景の相関、範囲依存が残る。
小さなχ²だけで、得られた幅を本番nσへ採用したり、信号帰属が確定したと判断したりしない。

最初に見る図：

- [Local ROI＋非対称応答の代表8点（4ページ、linear／log）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/split_recenter15_smoke_001/views_001/tag_peak_fit_diagnostics.pdf)
- [同じLocal ROIのGaussian比較（4ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/storage_recenter15_smoke_001/views_001/tag_peak_fit_diagnostics.pdf)
- [Git ROI＋非対称応答の代表8点（4ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/split_normal_smoke_001/views_001/tag_peak_fit_diagnostics.pdf)
- [以前のGaussian／左右指数tailの比較（8ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_003/response_comparison.pdf)

## 95点×4核種のsurveyと品質検査

R=0.405〜4.995 GeV/cの明示抽出95点／核種を、Git ROIでGaussianと非対称応答により比較した。
全600 native R binやd・tの0.4未満まで調べ終えたという意味ではない。
非対称応答は380スライスを保存し、raw・成分予測・16,654行のbin対応表を再読込検証した。
linear／logの190ページPDFも全ページを検査した。

表は収束・誤差等の必要条件を満たし、他核種の元ピーク中心がfit範囲内になく、
target yield／有効統計が正で、target modeがfit範囲内にある点数。
**2／3／5は比較用のχ²水準で、採用基準を決めたわけではない。**

| Species | Gaussian ≤2 / ≤3 / ≤5 | Split-normal + Float16 ≤2 / ≤3 / ≤5 |
|---|---:|---:|
| d | 2 / 2 / 2 | 2 / 3 / 3 |
| t | 2 / 2 / 5 | 13 / 20 / 28 |
| ³He | 4 / 9 / 15 | 18 / 21 / 25 |
| ⁴He | 15 / 27 / 31 | 27 / 34 / 42 |

改善しても、点数だけでは1本の多項式を全域に採用できない。
Rの支持範囲・空白、誤差の大きさ、ピーク両側の範囲、曲線の残差・微分・正の幅を別途確認する。
高Rだけの良好点から、未支持の低Rへ外挿して埋めない。

小さなχ²の落とし穴として、Gaussianの`t_b0061`と`he4_b0051`はtarget yield=0なのに、
status=0・covariance良好・χ²/ndf≈1だった。巨大だが有限な形状誤差も返していた。
これらはピークを同定できておらず、校正・純度・低R vetoの根拠に使わない。
品質検査は保存fitのフラグを書き換えず、別のJSON／TSVへ記録した。

- [Gaussian surveyの品質表](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/tof_tag_survey_001/quality_002/quality.tsv)
- [非対称応答surveyの品質表](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/split_normal_survey_001/quality_002/quality.tsv)
- [非対称応答surveyの190ページfit集](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/split_normal_survey_001/views_001/tag_peak_fit_diagnostics.pdf)

## 検証と未完了の区別

- Python回帰試験128件、Float16のROOT5照合57,075項目、非対称応答の人工データ試験6,773項目が成功。
- 同じproductionの実PicoDst代表2本でFloat16保存形式を確認。ただし全入力のschema監査ではない。
- 本番6入口は[診断後の静的再監査](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_git_restore_20261020_001/restore_audit_after_diagnostics.json)も成功。
- 元ROOT・過去PDFを保持。新係数の適用、job投入、hadd、git commit/pushは実施していない。

**次に残すもの**：現在の応答を用いた全surveyのLocal ROI追試、支持点からの多項式再校正とCOLZ重ね描き、
ALLでの共通曲線を使う核種分離、純度S/(S+B)と核種別cutoff／高R上限の評価。
全surveyのLocal ROI追試は未実施で、準備ツールの拡張も保存・検証済みとはしていない。
現時点で新しいPurity vs R図や新cutoffを確定した、という扱いにはしない。
