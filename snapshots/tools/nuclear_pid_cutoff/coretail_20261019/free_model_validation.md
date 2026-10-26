# 自由2Gaussian応答の検証方針（2026-10-19）

## この文書の位置付け

Gaussian core＋tailの固定profile比較に続き、同一核種を2本のGaussianで表す自由モデルを
調べる際の独立した検証案をまとめる。**科学的な採用基準の提案であり、本番採用・新しい数値閾値の決定ではない。**
旧ALL native count、旧局所fit窓、旧PID gateを保持した応答比較と、元多項式gateのpurity再評価は区別する。
固定profileで得た図・良いχ²・数学self-testだけでは、自由2Gaussianの核種帰属やpurityを認証しない。

対象はAu+Au FXT 3.85 GeV・2021 P24iyの既存部分統合ROOT由来のinclusive track QA。
本書は本番main、PID係数、R veto、元ROOT、既存fit/PDFを変更しない。
現在の実行状況・出力先は[当日Todo](../../../Todo/2026/10/19/README.md)と
[固定profile比較の説明](README.md)を参照する。

## 1. 全線yieldと局所窓yieldを混同しない

各Gaussianを数学的な全実数軸で正規化したdensity `g_j(x)` とし、保存fit窓を `W` とする。

```text
A_j = integral_W g_j(x) dx
Y_j = fit窓内のcolumn yield
N_j = Y_j / A_j                  （モデルによる全線yield）
f_tail = N_tail / (N_core + N_tail)
```

したがって、`Y_tail / (Y_core + Y_tail)` は局所窓で観測される割合であり、全線tail fractionではない。
特にtailが広い／窓の外へずれる場合、2つの割合は大きく異なる。
`A_j` は実際にfitした中心・幅・R node・重みから計算し、seed値や描画後のrebinから推定しない。
R nodeを含む場合は、各Gaussian列の実際の正規化kernelの窓積分を用いる。

全線yieldは観測窓外へのモデル外挿を含む。検出器の有限dE/dx軸のcountや、物理的に正のdE/dx範囲の
yieldと同一ではない。有限軸acceptance、局所窓acceptance、窓外外挿への感度を別々に記録する。
非常に小さな `A_j` によって巨大な全線yieldが得られても、データがそのtailを支持したとはいえない。

必要な保存値は、少なくとも各列の `Y_j`、`A_j`、`N_j`、有限軸acceptance、
全線fraction、局所窓fraction、中心・幅・その不確かさである。

## 2. core／tailのラベルと核種のラベルは別の問題

2Gaussianには成分ラベルを交換しても同じ予測になる自由度がある。
このstudyで「core＝狭く、全線yieldでも多数の成分」と呼ぶなら、次を同時に確認する。

- coreの幅がtailより狭い。
- coreの**全線**yieldがtailより大きい。peakの高さや局所窓yieldの大小では判定しない。
- 極端に近い2つの中心・幅、ほぼゼロの片方のyield、強い相関による非識別性を記録する。
- seed・初期値・最適化経路を変えても、同じ物理的な核種に帰属する解へ戻るか比較する。

狭い成分と多数成分が一致しない解を、無言のswap・clip・再fitによって合格にしてはいけない。
2Gaussianとして記述できても、上記のcore＋tail解釈を支持しない場合がある。
その場合は「ラベル契約不成立」等の理由を残し、応答校正への採用を保留する。
閾値や許容差が必要な比較は、明示設定と感度比較を用意してから判定する。

さらに、同一核種内のcore／tail交換と、d・t・³He・⁴He間の核種ラベル混同を分離する。
ある核種のtailが隣の核種のpeakを取り込むだけでもχ²は改善し得る。
各Rについて、対象と競合核種の中心・幅・全線yield・選択gate内寄与を並べ、
旧Gaussian／log-normal、TOFで観測したpeak位置、隣接Rとの整合を比較する。
対象だけ自由にする場合と、局所競合核種も自由にする場合を混ぜずに保存する。

## 3. 数値収束と物理的な支持を分離する

| 判定層 | 確認する内容 | 不成立時の扱い |
|---|---|---|
| 数値 | 最適化終了状態、有限予測、尤度・勾配、境界、同一初期値以外からの再現性 | 数値不成立として保存し、省略しない |
| 不確かさ | 対象shape/yieldの共分散、正定性、誤差、強相関、外挿fractionの安定性 | 対象の評価が不可能なら `unassessed` |
| 核種帰属 | core/tailラベル、他核種peakとの交差、対象／競合成分の識別 | 帰属未支持として `unassessed` |
| 応答の妥当性 | wide/local残差、gate内残差、窓・背景・初期値への依存、R方向の整合 | 理由付きで保留し、良いχ²だけで採用しない |
| purity | 同じ対象gateに対する核種総寄与と総候補、複数妥当モデル間の比較 | 未評価と支持された低purityを分ける |

共分散が存在することは必要な情報だが、それだけで完全な不確かさが検証済みとはならない。
shape固定時の条件付きyield共分散と、自由shapeを含む不確かさは同じではない。
対象2成分や、そのgateへ無視できない寄与を持つ競合成分の誤差・共分散が評価できない場合、
小さな形式的誤差やpurityの中心値だけを使って採用しない。
ゼロyieldの非寄与nuisanceに未定義shapeが残る場合も、寄与の上限を明示して区別する必要がある。

片方のGaussianがゼロyield境界へ達したとき、その中心・幅は識別できないことがある。
その2Gaussian解を共分散検証済みと呼ぶのではなく、別に保存した1Gaussian controlとの関係を記録する。
不確かさの検証が必要なら、profile・再標本化・toy等の追加試験を明示的に設計する。

## 4. χ²、自由度、ΔDの読み方

図のPearson χ²とconditional ndfは実保存値をそのまま用いる。
ndfでは、そのfitで実際に自由なshape、背景、内部yieldを考慮し、固定profileの値を自由parameterとして引かない。
境界・強い縮退・profile選択がある混合モデルでは、表示したχ²/ndfから
単純な標準分布に基づくp値や独立な適合度検定を自動的に主張しない。

2Gaussianで同一核種に2列を割り当てた場合、片方の列を落としたときの `ΔD` は
**その核種の有意度ではない**。もう片方が同じ核種のyieldを保持できるためである。
必要なら同一核種の両列を同時に除き、他核種・背景を再最適化した群としての比較を別途行う。
それでも境界や非識別性のため、`ΔD` をそのままGaussianのσ有意度へ変換しない。
単列の `ΔD` はその名前のまま記録し、旧1列核種の有意度screenを無条件で流用しない。

## 5. purityで最良モデルを選別しない

1つの核種の信号寄与は、その核種のcoreとtailを合わせて計算する。
tailを別の「背景」として捨てると、核種応答のモデル変更とpurityの定義変更が混ざる。
一方、隣接核種のcore／tailは対象核種に対する背景として残す。

```text
S_gate = 対象核種の両Gaussianを、実際の固定PID gate上で積分したcount
T_gate = 同一gateを満たす全候補count
P_gate = S_gate / T_gate
```

`T_gate` が実測gate countか、全fit成分の予測和かを明記し、両者とgate内残差を比較する。
モデルの積分と観測countが不整合なまま、比を都合よく切り詰めない。
信号和が負／非有限、分母が未支持、purityが物理範囲外の場合は理由を保存し保留する。

採用候補モデルは、purityの高さを見る前に応答・背景・核種帰属の妥当性で決める。
妥当なGaussian／log-normal／2Gaussianと背景の比較を保存し、
高いpurityを返したprofile・R点だけを残さない。
全Rのモデル比較、失敗点、旧vetoの位置と両隣を表示し、穴の新規採用を自動化しない。
ユーザーが指定したpurity基準との比較は、妥当なモデルの差を考慮した下側で行う。
未支持モデルを無理に数値envelopeへ入れることと、悪いモデルを無言で除くことの両方を避け、
除外理由と未評価領域を明示する。

旧gateに対する参考purityは、元多項式で変更されたgateのpurityへ転用できない。
本番へ反映する前には、対象とする新gateで候補count・各成分積分・モデル差を評価し直す。

## 6. 滑らかな校正関数への採用は別段階

自由2Gaussian fitの中心・幅は各Rの応答記述であり、そのままPID中心・nσ幅の新しい表ではない。
core幅、混合RMS、有限R投影の幅、検出器のintrinsic幅は区別する。
固定Rでの混合分散には、成分内分散に加えて中心間の分離が寄与する。

```text
Var_mix = (1-f)*sigma_core^2 + f*sigma_tail^2
          + f*(1-f)*(mu_tail-mu_core)^2
```

R nodeを積分した投影では、さらにRに伴う中心の変化が幅へ寄与する。
どの中心・幅を校正対象にするかを決め、元の対数多項式等の滑らかな関数でR方向の整合を検証してから判断する。
COLZへの関数重ね書きと各Rの1次元fitの両方を確認し、独立fit点の折れ線へ無言で戻さない。

## 7. 描画互換性

[共通painter](../PlotCutoffSpectra.C)と
[品質注記macro](../quality_20261019/PlotQualityDiagnostics.C)は、
ROOT内に保存された次のhistを読むため、2Gaussianを4核種へfoldして保存すれば利用できる。

- `raw_native`：元count・error・flowを保持。
- `prediction_total`、`prediction_smooth`。
- `prediction_d`、`prediction_t`、`prediction_he3`、`prediction_he4`：各核種の2成分和。
- `prediction_fit_mask`：窓外を予測ゼロとして描かないためのmask。

component予測とtotalの加算一致、閉じたROOTのreadbackを先に検査する。
窓外の全線外挿yieldを、観測済みのcount曲線としてwide panel全域へ描き足さない。
実数値の `numericalValid`、Pearson χ²、conditional ndfを保持し、科学的支持が未評価なら
`supported=0, supportKnown=0` として「not evaluated」を表示する。

ただし[pilot renderer](render_pilot.py)のPython準備処理は、既定の `--mode fixed` では
48 jobs・4 fixed profiles・12冊64頁のplanを明示的に検査する。
自由2Gaussianでは `--mode free_coretail` を明示し、専用adapter・engine・responseの宣言と、
24 jobs・Gaussian control＋自由2Gaussian seedの2組・12冊32頁を検査する。
job数だけから自由モデルと推定せず、既存の固定profile結果も上書きしない。
どちらの場合も失敗fitを省略せず、全入力SHA、全PDFの閉じたreadback、ページ数、各panelの軸記録を保存する。

## 8. 次の結果に必要なチェックリスト

- [ ] 全Gaussian列の窓yield・acceptance・全線yield・fractionを記録する。
- [ ] core/tailラベルと他核種への帰属を独立に検証する。
- [ ] 初期値・窓・背景・自由にする競合核種を比較し、失敗も保存する。
- [ ] 共分散・境界・非識別性を確認し、未評価と低purityを分ける。
- [ ] 単列ΔDを核種の有意度として流用しない。
- [ ] 旧gate参考比較と新多項式gateの再評価を分離する。
- [ ] purityが高いこと自体をモデル採用基準にしない。
- [ ] COLZ、wide/local・linear/log、R方向の変化を同時に確認する。
- [ ] 滑らかな校正の確定前に、本番係数・R veto・本番mainを切り替えない。

このtrack QAから、Λとのpair条件付きpurity、k*別purity、最終CF補正、Λや原子核の真の再構成効率を
直接測定したとはいえない。図・レポートでも、そのような表現へ置き換えない。
