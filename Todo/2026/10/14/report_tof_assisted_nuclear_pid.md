# TOF-assisted原子核PID study: 進捗・結果

## 現在の状態

最新版の[Todo](STAR_Todo_20261014.txt)と[計画](plan_tof_assisted_nuclear_pid.md)を確認した。
対象は3.85 GeV・2021 P24iyの既存高統計ROOT。
入力cacheの作成と独立監査、粗いR分割とnative分割でのTAG初期fitが完了した。
単一Gaussianとlog-normalの全域比較、局所残差・COLZ・比較PDFの作成まで完了した。
全4核種の校正／ALL purityを確定できる適合性には達していない。
校正係数の採用・本番解析の変更はまだ行っていない。
ユーザーは「必要なら低R／高Rを連続接続する関数も比較する」ことを許可した。

指定TXTは初回読取り時には10/13と同内容だったが、保存し直された後の最新版を現在の依頼として使う。
古いGit整理・TOF拡張タスクの点検は停止し、移動・削除・Git操作は実施していない。

## 入力監査

入力は`auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root`。
4本ともTH2F、Rは600 bins・0〜6 GeV/c（0.01幅）、dE/dxは1200 bins・0〜120 keV/cm（0.1幅）。
Sumw2は未保存。誤差は保存countからのPoisson型で、stored countをそのまま使う。

| Histogram | 全保存binのcount合計（核種収量ではない） |
|---|---:|
| `true/hDedxP` | 316,515,592,807 |
| `true/hDedxP_d_m2` | 3,555,393,696 |
| `true/hDedxP_t_m2` | 528,861,766 |
| `true/hDedxP_3He_m2` | 252,915,567 |

全723,604セルで軸一致、有限非負count、`ALL−TAGd−TAGt−TAG3>=0`がPASS。
ALLのGetEntriesは316,515,591,143で、stored countとの差1664はTH2Fの丸めと整合する。
これらはtrack QAの数であり、入力イベント数や各原子核の真の数ではない。

## 低Rの初期確認

R=0.40〜0.50のタグ分布は低いdE/dxの一般粒子／誤TOF等の成分が多い。
例えばTOFtは全域522,534 countだが、55〜120 keV/cmは2,101だけ。
同じ核領域のTOFdは3,924、TOF³Heは36,083 count。
これだけで⁴He中心やpurityを確定しない。

粗い非重複R分割で、R=0.8〜0.9のd／t／³Heタグの全域modeはそれぞれ約12.75／23.25／35.25 keV/cm。
これはピークfitの中心測定ではなく、帯の分離を確認する手掛かりである。
dタグ同域の45〜120 keV/cmには3,633,646 countあり、高い側の候補を調べられる。
dタグを純粋な⁴Heとせず、他タグのshapeも使って核種の分担を検証する。

0.005刻みの隣接R指定は同じnative binを重複利用するため使わない。
非重複0.01以上の区間と、明示した隣接bin合併でstudyする。

## 未実施・今後の作業

1. TAGのGaussian＋背景fit、Rごとの支持域と関数の比較。
2. TAGのshapeを使うALL fit、実選択窓のpurityと系統差。
3. 検証済み校正候補YAML、支持域外の扱いの相談、必要に応じた解析profile更新。

正式PDF・数値表・結論は[Gaussian／log-normal比較のまとめ](comparison_gaussian_lognormal.md)へ記録した。
旧studyと現行PIDは保持する。

## 入力cache・表示の完成

新しい出力は`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001/`。

- [TOFタグとALLのCOLZ（12ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001/source_colz_atlas.pdf)。
- [Rごとの全dE/dx、log／linear（24ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001/coarse_full_dedx_atlas.pdf)。
- [native bin監査](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001/native_audit.json)。

ALL／3タグ／非重複remainingの5個のTH2D cacheを作成。
remainingの誤差は部分集合の相関を考慮して`Var=Nremaining`とし、ALLとTAGの分散を独立に加算していない。
26製品・9ソース・17保護対象のSHA、元ROOT fingerprint、全723,604セルを独立に照合してPASS。
全軸0〜120 keV/cmの図を保存し、低dE/dx側を隠していない。

## 実装の数値検証（物理的な核種識別の検証ではない）

- ROOT5のTAG core／macro: 2,088項目PASS、コンパイルPASS。
- ROOT5のALL／purity macro: 15,024項目PASS、コンパイルPASS。
- 現在のC++係数・校正knot・strict `R<2.5`を読み取り、実際のPIDを再現する6テストPASS。
- 曲線・支持域の判定を含むPython回帰試験は最終17テストPASS。
- 別実装のlog-normal core／macroはROOT5数値3,093項目PASS、コンパイルPASS。
- Gaussianはbin積分、yieldは非負、背景は正値。boundary・zero-yield・未支持成分を別途保存する。
- 新しい候補関数は中心3次／幅2次のlog(R)式と、滑らかなlog(R)の区分3次式を比較する。
  未支持域・未測定のgapを検証済み校正として外挿しない。

dタグにはd帯以外の高いdE/dx帯も見えるため、未知成分を⁴Heへ強制的に帰属しない。
最初に0.80〜0.90 GeV/cで背景modelを比較し、その後0.30〜2.50 GeV/cのnative分割へ進んだ。

## Gaussian応答による初期fit

初期fitは採用結果ではない。TAGd中のt／³Heは初期値の旧関数を固定した未検証のnuisance成分。

| 分割 | 背景 | スペクトル数 | 数値収束した数 |
|---|---|---:|---:|
| R=0.80〜0.90、幅0.05 | linear | 6 | 6 |
| R=0.80〜0.90、幅0.05 | exponential | 6 | 3 |
| R=0.30〜2.50、幅0.01 | linear | 660 | 646 |
| R=0.30〜2.50、幅0.01 | exponential | 660 | 604 |

native分割のR=0.805で、linearの名目ピークは次のとおり。
幅は射影Gaussianの幅で、最終的なintrinsic幅や核種帰属を確定した値ではない。

| 名目成分 | 中心（keV/cm） | 幅（keV/cm） | 当該スペクトルのχ²/ndf |
|---|---:|---:|---:|
| d（TAGd） | 13.791 | 1.082 | 468.35 |
| t（TAGt） | 25.637 | 1.642 | 76.81 |
| ³He（TAG³He） | 37.510 | 2.103 | 404.63 |
| ⁴He候補（TAGd） | 50.936 | 3.007 | 468.35 |

中心はRに対して滑らかで背景model間でも概ね近いが、数値収束だけでは残差を説明できていない。
TAGdはd収量が大きいため全体χ²をdが支配する。⁴He近傍の局所残差を別途診断する。
目視ではdピークの高dE/dx側の裾と、単純背景で再現できない谷がある。
粗いR幅の平均化だけではなく、非Gaussian応答・背景構造・保存粒度を切り分ける。

出力は`tag_smoke_{linear,exponential}_001/`と`tag_native_{linear,exponential}_001/`。
粗い分割の図は`tag_smoke_linear_views_001/`へ保存した。
紙面内へ収めた正式[全軸・残差の図](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/tag_native_linear_core_views_002/tag_fit_residuals.pdf)は
native分割のR=0.80〜0.90。表示追加だけで、fit値は変更していない。

## 関数候補の初期判定（未採用）

`curves_native_001/support_audit.json`に背景依存・数値収束・境界・χ²・統計・誤差の判定を保存。
数値shapeとして条件を通った点はd=12、t=16、³He=18、⁴He候補=2。
この数は真の核種として支持された点数ではない。
⁴He候補の2点はR=0.465／0.475で、中心約66.8／70.1、幅約13.6／10.4と広く、
同域の³He成分と重なるため核種帰属は未確定。最低6点の関数化条件も満たさず、⁴He曲線は作成していない。
低統計でχ²が小さいことを理由に、これを信頼できる⁴He下限とはしない。

`logpoly.json`／`pchip_logR.json`は比較用候補で、`productionEnabled=false`、核種帰属・転用・幅の校正は未検証。
ガタガタの点補間を新本番設定へ置き換える段階にはまだ達していない。

## log-normal比較の結果と次の検討

Gaussian結果を保持したまま、log-normal応答も両背景で全660スペクトルずつ比較した。
パラメータはピークmodeと局所曲率の幅であり、分布の平均／標準偏差と混同しない。
その定義ではGaussianのR-variance差引きを無条件に流用しない。
fit窓・native count・背景の比較条件を固定した。tのR=0.80〜0.90は改善したが、全域では両モデル収束区間でも
log-normalのdevianceが大きく、一括採用は支持されない。低R⁴Heの同定も未確定。
ユーザーはまずGaussian／log-normalを比較する方針を選択したため、実測TOF形状テンプレートは追加していない。

1320組の両モデル入力照合、PDF紙面内の軸・注記、生成物SHA、原ROOT fingerprintがPASS。
次の検討は応答をGaussian 2成分等の関数で記述するかどうか。
現段階の不良適合を隠してALLの確定purityや本番PIDへ進めない。詳しい表・図と未完了項目は比較まとめを参照。
