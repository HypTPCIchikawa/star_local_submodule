# Task 3: 0.40開始free fitを起点とする反復校正

## 現在の状態

[計画](task3_helium_iterative_from040_plan.md)を保存した。
旧0.40開始PDF・fit ROOT/JSONは存在し、出発点に使う。旧PDF・本番校正・元ROOTは保持する。
新しいfree→平滑化→共通関数fitと、nσ／元dE/dx両軸の表示を実装した。
低剛性10区間と全210区間・14段階の実fit／描画は完了。
最初の全域共通式fitは未収束で低R再現も悪化したため、係数は未採用。
背景初期値・数値条件・係数初期値の5通りの追加比較も完了したが、良好な全域共通式は得られていない。
既存TOF制御分布との比較は完了したが、低Rの⁴He中心を独立に決める統計・純度は不足した。
同じ式で0.40〜0.70に限定し、最大120反復まで追加した診断も完了した。
成分の消失は避けられたが、残差と数値収束が未達のため本番の係数採用は行っていない。
追加依頼のraw表示0〜120の図と、生COLZへの低R関数重ねも別出力へ作成・確認した。
旧PDF・fit数値は保持し、表示専用図による追加fit・本番PID変更はない。

旧fitの探索制約、今回の数学的guard、後段の平滑化penaltyは分けて記録する。
低Rで良好な形状fitが得られても、³Heのtailと⁴Heの核種帰属は別の検証事項である。

## 実装・数値検証の進捗

- 元の0.40開始fitの全210区間から、旧TF1の高さ・中心・幅をraw単位と窓内収量へ変換した。
  Gaussian全収量への換算には`σ_raw sqrt(2π) / 0.5 keV/cm`を用いる。
  旧±2nσ窓の信号countを独立に復元するテストを含むPython 11 testsがPASSした。
- `smoothing.py`は正規化したlogR基底のweighted QRと、二次式の解析最小値による正値幅を実装した。
  正値だがBernstein係数の一部が負になる式も保持する数値テストがPASSした。
- 新しい描画器は人工5 testsがPASS。counts/error/flow、raw↔nσの同bin性、成分和とfit-maskを検査した。
  全区間のnσ／raw dE/dx linearと、低Rの4面log/linearを独立PDFに出す。
- ROOT5／SL24y／GCC4.8.5のfitコア数値テスト1,390項目がPASSした。
  [検証記録](validation/iterative_helium_20261014_002/validation_manifest.json)。
  初回の自己テストnamespace compile失敗は`001`へ保持し、修正後を別`002`へ保存した。
- 共通係数候補は未採用。人工テストの成功を実データfitや核種同定の成功とは扱わない。

旧`he4_low040_components.pdf`・旧fit ROOT・現在の本番校正YAMLは保持している。

## 低剛性10区間の区間別free fit

`iterative_from040/low_rigidity_free_20261014_001/00_local_free/`へ保存した。
R=0.40〜0.50の10区間すべてでstatus 0、AMP KKT収束、有限予測、幅・中心の境界なしを確認した。
Gaussian入れ替えはなく、最大AMP KKT残差は約1.93×10⁻¹¹。
同じ窓での旧採用modelのdata deviance合計6,231.025から6,149.750へ改善した。
改善は最初の区間が主で、0.415〜0.495では旧modelの局所解をほぼ再現する。

R=0.405の⁴Heラベル成分は中心81.387、幅9.232 keV/cmとなった。
次のR=0.415は中心82.021であり、従来providerのlegacy境界による約22.7 keV/cmの跳びは
このstudyのfit点には含まれない。ただしfree点には小さな非単調が残るため、以降の共通関数fitを行った。

- [低R全10区間・nσ／元dE/dxのlog/linear 4面](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_rigidity_free_20261014_001/00_local_free/views/low_rigidity_components.pdf)
- [nσのlinear図](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_rigidity_free_20261014_001/00_local_free/views/fits_nsigma_linear.pdf)
- [元dE/dxのlinear図](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_rigidity_free_20261014_001/00_local_free/views/fits_raw_dedx_linear.pdf)

数値収束と残差品質は別で、Pearson χ²/ndfは0.405で約15.9、0.415で約7.46、0.425で約4.33。
低Rの緑Gaussianは幅広の肩の成分であり、独立したraw第3ピークや⁴Heの真の収量を確定するものではない。
図にはpoor/reference・品質未確認・核種帰属のmodel依存を明示する。

全域の反復studyは`iterative_from040/full_iterative_20261014_001/`へ別保存した。
free段階の全210区間ではAMPは収束したが、2区間はminimizerの提案stepを非改善として保持した。
全域の数値収束・残差品質を10区間の成功だけで断定せず、後段の結果を確認した。

## 全域の反復結果: 滑らかさだけでは採用できない

| 段階 | Raw data deviance合計 | Penalty合計 | 判定 |
|---|---:|---:|---|
| 区間別free | 1,920,542 | 0 | 全AMP収束。2非改善stepを保持、3境界sliceを区別 |
| λ=0.05、3反復 | 2,247,515 | 649,907 | 数値収束するが低Rの形状・成分分担が悪化 |
| λ=0.2、3反復 | 2,762,810 | 893,994 | 低Rで³He窓内収量が0となる例あり |
| λ=1、3反復 | 4,000,678 | 897,621 | 平滑化を優先した形状変化 |
| λ=5、3反復 | 5,349,449 | 372,867 | 局所数値収束でもデータ再現は保証しない |
| penaltyなし、共有14係数・GL3 | 5,541,915 | 0 | 未収束、status 3、stationarity評価に失敗点あり。未採用 |

全段階でfit窓・観測countは同一。raw devianceとpenaltyを混同しない。
R=0.425ではfreeのχ²/ndf≈4.33からλ=0.05第3反復で≈172.41へ悪化した。
T Gaussianの中心が61.34→69.82、幅4.37→6.75へ移り、³He窓内収量が4.285M→0.265Mに低下した。
背景成分が³Heの帯を吸う別の分担に移っており、曲線の滑らかさを正しさとは評価できない。

最終候補のR=0.425の³He中心は78.60、⁴He中心は89.32 keV/cm。
freeの中心71.03／81.42とは異なる。TOF条件付き³He窓のmean≈69.10ともずれる。
TOF mean自体も純粋な核種のGaussian中心ではないが、生COLZでも低Rのずれが明瞭である。
最終共有fitのEDM・covarianceは未収束のため、係数精度として使用しない。

- [free→反復→最終の中心・幅・残差まとめ 3ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/full_iterative_20261014_001/iterative_calibration_summary.pdf)
- [最終候補のnσ／元dE/dx・低R 4面 10ページ（未収束参考）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/full_iterative_20261014_001/13_common_raw_refit/views/low_rigidity_components.pdf)
- [生COLZへの候補重ね 4ページ（未収束参考）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/common_candidate_raw_overlay_20261014_001/raw_dedx_calibration_overlay.pdf)
- [元dE/dx射影・TOF帯との比較 25ページ（未収束参考）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/common_candidate_raw_overlay_20261014_001/low_rigidity_raw_projections.pdf)

各14段階に全区間のnσ linear・元dE/dx linear各53ページと、低R 4面10ページを保存した。
全描画でcounts/error/flow・rebin・同じ表示bin・予測成分和がPASS。最大相対成分和差は約3.8×10⁻¹⁶。
描画の整合とfitの良さは別である。初回candidate YAMLはstudy-onlyで、PIDへ未適用。

### 関数domainの参考切り分け

同じ3次／2次式を0.40〜0.70のfree点だけへQR fitする診断も別保存した。
R=0.425の³He中心は、全域QRで74.23、低R限定QRで71.11、freeで71.03となった。
⁴Heは全域QR83.48、低R限定QR80.81、free81.42。
全域を一本の低次数式で合わせることと、局所分解の自由度の両方を切り分ける必要がある。
この低R限定式は投影幅の初期化診断だけで、R-bin平均raw fit・核種同定・本番のpiecewise校正ではない。
[数値表](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/domain_diagnostic_20261014_001/domain_comparison.json)。

## 追加の原因切り分け

同じ最終係数を初期値に、T／背景の開始値をfreeと最後のregularizedで変えて、penaltyなしcommon refitを比較した。
さらにTの中心・幅をfreeデータの値へ固定する条件付き比較を行った。
legacy式への固定ではなく、低Rの背景が³Heを吸う自由度の診断である。
AMP profile失敗と、窓acceptance等の数学条件による不連続を分けて確認する。
追加処理は別studyに保存し、初回14段階の結果は保持する。

### 共通fitの数値不安定さ: 実データで確認した原因

[独立probe](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/diagnostic_case_20261014_001)で
最終係数・T／背景・窓内振幅を復元し、基点のdevianceが保存値と完全一致した。
正規化した14座標の±probeでは、厳しいAMP条件で6/28 probeが失敗した。
原因はすべてKKT残差の小さな超過（約1.05×10⁻⁸〜1.93×10⁻⁸）で、形状・acceptance・templateの失敗はなかった。
反復上限ではなく、4〜8反復でline-searchが止まる極近傍の数値精度問題である。
無効点を1×10²⁰とする目的関数に人工的な急変が入り、Minuitを不安定にする。

AMP許容残差を1×10⁻⁶へ変えると全29 case・全210点が有効となり、
厳しい条件と有効caseを比べた最大deviance差は約0.000377。
一方、元の最終係数の勾配/countは約0.000775で、設定したstationarity閾値の約775倍。
数値判定を緩和しただけで収束済みと扱わず、この条件で最適化を再実行した。

共通式の再fit専用入口は`tools/he4_purity/iterative_logpoly/run_common_refit.py`。
元studyの全製品・source・入力fingerprintを照合して別studyへコピーcompileする。
数値設定は`common_refit_numerical_settings.yaml`（AMP 1×10⁻⁶、最大1000反復、共通block最大12反復）。
物理cut・fit窓・背景model・R平均・本番PIDは変えない。係数初期値を含む5通りの比較は完了した。

### 数値条件修正後の5通り: 低Rの問題は残る

| 共通係数の初期値 | T／背景の開始値・条件 | Total raw D | R=0.425のraw D | 状態 |
|---|---|---:|---:|---|
| 区間別free（比較基準、共通式ではない） | 各R自由 | 1,920,542 | 481 | 低Rの形状参考 |
| 初回最終係数 | free初期値、T可変 | 5,825,030 | 44,327 | ³He成分消失、広いTが吸収 |
| 初回最終係数 | 最後のregularized初期値、T可変 | 5,002,930 | 45,222 | 同様、block反復未収束 |
| 初回最終係数 | free初期値、T中心・幅固定 | 15,486,032 | 517,456 | 下側構造を背景が吸収 |
| free点の初期QR係数 | free初期値、T可変 | 4,993,955 | 44,908 | 初期値変更でも³Heが消失 |
| free点の初期QR係数 | free初期値、T中心・幅固定 | 15,486,032 | 517,450 | 最終係数開始とほぼ同じ不良形状 |

全比較でobserved count・fit窓・model・GL3は同一。
AMPは全210区間で収束、globalの軸方向probeも通り、最初の極近傍KKT判定の障害は解消した。
ただしoverallの数値収束は全5候補でfalse。条件付きglobal最小化だけを全blockの収束と扱わない。
曲線の正値・核種名の順序・滑らかさも、観測ピークの正しい分担を保証しない。

たとえばfree初期QRからのT可変caseでも、R=0.425の³He窓内収量は0、
⁴He中心は90.70（freeは81.42）へ移る。
T固定caseは³He中心77.12となり、freeの71.03とは合わない。
固定・初期値変更だけで解決したとはいえず、全域3次／2次式と背景分解のモデルをさらに切り分ける必要がある。
「原式では絶対に不可能」とはまだ断定しないが、今回の候補の本番採用はしない。

| Study | 低Rの元dE/dx・nσ 4面PDF |
|---|---|
| free nuisance・最終係数開始 | [10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/free_nuisance_common_20261014_001/views/low_rigidity_components.pdf) |
| regularized nuisance・最終係数開始 | [10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/regularized_nuisance_common_20261014_001/views/low_rigidity_components.pdf) |
| free T形状固定・最終係数開始 | [10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/free_fixedT_common_20261014_001/views/low_rigidity_components.pdf) |
| free nuisance・初期QR開始 | [10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/initialqr_free_common_20261014_001/views/low_rigidity_components.pdf) |
| free T形状固定・初期QR開始 | [10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/initialqr_fixedT_common_20261014_001/views/low_rigidity_components.pdf) |

### 独立TOF制御による追加確認

保存済み`true/hDedxP_d_m2`はTPC PIDとは独立したTOF条件である。
TAGはALLの部分集合であり、ALLと統計的に独立した標本という意味ではない。
設定のm²/q²窓はdeuteron中心3.48096±1.5×0.141458。
⁴Heの`(3.72742/2)²≈3.473`も通るため、D／⁴He互換のTOF制御として高dE/dx部分を調べる。
純粋な⁴He選別とは呼ばず、deuteron・誤TOF対応・他核種等の混入を区別する。
この分布と³He制御・ALLを比べた。R=0.40〜0.50でDタグの全域countは1,356,992だが、
55〜120 keV/cmの核領域は3,924 countだけ。³Heタグの同領域は36,083 countである。
R=0.425のDタグ核領域は311 count、0.5 keV/cmの各binは概ね1〜15 count。
80／100 keV/cm付近に再現性ある独立ピークを決められず、どちらが⁴Heの正しい中心かは確定できない。
これは「⁴Heがない」という判定ではなく、TOFの受理率・純度・統計の制限がある対照標本である。
modeや窓内meanをGaussian中心・purityとして採用しない。

- [核領域55〜120 keV/cm、ALL／TOF³He／D・⁴He互換タグのlog/linear対照10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/tof_compatible_control_20261014_002/tof_controls_nuclear_dedx.pdf)
- [全dE/dx範囲0〜120 keV/cmの対照10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/tof_compatible_control_20261014_002/tof_controls_full_dedx.pdf)
- [生COLZの対照3ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/tof_compatible_control_20261014_002/tof_controls_nuclear_colz.pdf)

001を保持し、002は線形y軸の桁表示だけを指数表記へ修正した。
001／002の60 hist・43,320 binのcount/error/flow/axisが完全一致、数値TSVもbyte一致。
元ROOT・旧fit/PDFの不変を確認した。新fitや本番PID変更はない。

## 低domainのraw共通fit: 全域より改善するが12反復では未収束

上のQR初期化から、R=0.40〜0.70の同じ30区間をraw/R-bin平均で共通fitした。
各Rのfit窓・observed count・背景modelはfreeと同一。関数は中心3次・幅2次のままである。
このdomain選択は診断用で、解析入力cut・本番のpiecewise式ではない。
比較する全observed countは646,768,724。

| 同じ30区間での比較 | Raw D合計 | R=0.425のD | R=0.425のχ²/ndf | 判定 |
|---|---:|---:|---:|---|
| 区間別free | 169,477.494 | 481.070 | 4.33 | 比較基準、各R自由 |
| 共通式、T形状可変、最大12反復 | 316,704.793 | 3,443.891 | 29.86 | 成分は保持、block未収束 |
| 共通式、free T形状固定、最大12反復 | 502,223.710 | 19,553.365 | 167.05 | 下側の分担が変化、未採用 |

比較基準の30区間freeにもR=0.535の非改善proposalを保持した1区間が含まれる。
その値を最良解の保証とは扱わず、各caseで同じ保存値を基準にする。
R=0.425のfreeは数値収束している。

T可変では全30区間で³He／⁴Heの窓内収量が正となり、全域fitで見られた³He消失は避けられた。
R=0.425の³He中心69.771、幅5.467、⁴He中心82.189、幅7.534 keV/cm。
freeの³He71.031／4.293、⁴He81.419／7.839とは中心・幅・収量が変化し、残差は増えている。
最大12反復の終点でも係数driftは0.05054、local改善/countは2.666×10⁻⁶で、overall収束はfalse。
global status 0、covariance 3だけを全fitの収束とは扱わない。

T固定では最初の3区間のT収量が0となる。R=0.425の³He幅は9.247となって下側の2構造を吸い、
⁴He中心は87.176へ移る。T固定だけで解決したとはいえない。
全原本・snapshot SHA各456件と製品31件ずつを独立監査し、全30区間の入力窓/count一致を確認した。

- [T可変・低Rのnσ／元dE/dx、log/linear 4面10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_domain_free_common_20261014_001/views/low_rigidity_components.pdf)
- [T形状固定・同じ4面10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_domain_fixedT_common_20261014_001/views/low_rigidity_components.pdf)

各studyにはnσ linear・元dE/dx linearを各8ページでも保存した。
raw/nσ同count、rebin/error/flow、成分和監査はPASS。最大相対成分和差は約4.2×10⁻¹⁶。
T可変はDが単調に改善し続けていたため、同じ初期値・modelで上限だけ120反復へ増やす
`low_domain_convergence_settings.yaml`による別studyを追加し、数値予算の打切りとmodel品質を区別した。

### 最大120反復まで増やした結果

| 同じ30区間 | Raw D合計 | R=0.425のD | R=0.425のχ²/ndf |
|---|---:|---:|---:|
| 区間別free | 169,477.494 | 481.070 | 4.33 |
| 共通式、T可変、最大12反復 | 316,704.793 | 3,443.891 | 29.86 |
| 共通式、T可変、最大120反復 | 314,565.202 | 4,370.826 | 37.91 |

120反復のtraceは単調下降し、最初の12点は元の12反復caseと完全一致した。
D合計の追加改善は約0.676%だが、R=0.425の局所残差は逆に増えた。
全域共通fitより低Rの成分分担を保持できても、データ再現の改善とは同一ではない。
最終R=0.425の³He中心／幅は69.309／5.835、⁴Heは82.367／7.536 keV/cm。
全30区間で³He・⁴Heの収量は正だが、核種の独立した同定を意味しない。

global status 0、global軸方向probeはPASSだが、overallはfalse。
係数drift=5.734×10⁻⁴は閾値1×10⁻⁴を超え、R=0.665／0.695でlocal status 3が残る。
最後のΔDは約0.266でほぼ平坦になった。反復上限だけ増やせば解決するとは判断しない。
条件付きcovarianceを本番の係数精度へ採用しない。

- [120反復候補・低Rのnσ／元dE/dx、log/linear 4面10ページ（未収束参考）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_domain_convergence_common_20261014_001/views/low_rigidity_components.pdf)
- [同候補の参考係数YAML・domain（本番未適用）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_domain_convergence_common_20261014_001/analytic_candidate.yaml)

独立監査でsource snapshot15件・製品31件のSHA不一致0、全30区間のfreeとのmask/count一致を確認した。
元ROOT・旧成果物・本番PIDは保持している。現行の点補間校正も低Rの確定校正とは扱わない。
結論は「今回の同一関数family・背景model・初期値では十分な共通fitが得られていない」。
原式が数学的に不可能だという証明ではない。今後は関数domain／連続接続と背景shapeの仮定を見直すが、
本番の関数familyや選択条件の変更は今回自動実施しない。

## 追加表示1: 元dE/dxの全体形を見える範囲へ

`iterative_from040/wide_views_20261014_001/`にfree10と最新120反復候補の図を別保存した。
元dE/dxは**保存軸0〜120 keV/cmの全normal bin**、nσも同じbin・count・errorを表示する。
全域のlinear図では通常粒子の巨大ピークに核領域が埋もれるため、全域図を残した上で
T／³He付近を見やすくした45〜120 keV/cmの補助図も併記した。
fit窓は広げず、窓外はデータだけ。灰色線はfit境界であり、PID cutではない。

| 元のfit | 全域0〜120のnσ／dE/dx log・linear 4面 | 低Rの45〜120補助4面 |
|---|---|---|
| 区間別free、R=0.40〜0.50 | [10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/wide_views_20261014_001/free10/wide_low_rigidity_components.pdf) | [10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/wide_views_20261014_001/free10/supplemental_nuclear45_120_components.pdf) |
| 共通式120反復、低R部分 | [10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/wide_views_20261014_001/common120_30/wide_low_rigidity_components.pdf) | [10ページ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/wide_views_20261014_001/common120_30/supplemental_nuclear45_120_components.pdf) |

各フォルダには`wide_raw_dedx_log.pdf`／`wide_raw_dedx_linear.pdf`と対応するnσの各軸別図も保存した。
freeは各3ページ、共通式は30区間を各8ページ。R=0.425の全域・補助PDF本体を実画像で確認した。
全域のnormal count合計はfree10が10,655,283,430、共通30が29,518,164,728。
これはpre-PIDの**track数で、イベント数や³He／⁴He収量ではない**。
全native bin12,000／36,000個が表示され、低R補助は7,500個。
raw/nσ同bin・同error、rebin/count/variance/flow保存、成分和監査はPASS、最大相対成分和差4.26×10⁻¹⁶。
元source／製品SHAと入力ROOT fingerprintを前後で照合した。

補助PNGにはROOT5のexportによる文字のclipが見られるが、PDF本体は正常でありPDFを参照する。
PNGも旧製品として保持し、数値や元PDFを上書きしていない。

## 追加表示2: p/q vs dE/dx COLZに今回の関数を重ねる

[低R関数の生COLZ重ね・4ページPDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_domain_raw_overlay_20261014_001/raw_colz_low_domain_overlay.pdf)
を別出力した。元ROOTの`true/hDedxP`（PID前、正負電荷を合算）にTF1を`Draw("same")`する表示専用処理である。

1. 保存dE/dx全域0〜120 keV/cmのCOLZと、最新の³He／⁴He中心・±1σ・±2σ。
2. 同じCOLZの核領域20〜120 keV/cm。色軸はlogで、低側の複数の帯も見える。
3. 最新120反復・12反復・初期free点QRの幅比較。
4. 同じ生COLZへの3条件の中心比較。

式は解析コードと同じ`L=log10(|p/q|)`、`μ=c0+c1 L+c2 L²+c3 L³`、`σ=d0+d1 L+d2 L²`。
評価域は0.40〜0.70 GeV/cのみで、外へ延長しない。表示のR範囲は0.35〜0.80である。
青は³He、赤は⁴He。実線が中心、破線が±1σ、点線が±2σ。
ここでσは共通fitのGaussian intrinsic幅で、係数の不確かさの帯や核種同定・purityの保証ではない。
初期QRの幅は投影幅のseedなので、共通fitのintrinsic幅と同じ量として評価しない。

1,806点でPythonの元係数評価とROOT TF1を照合し、最大相対差6.45×10⁻¹⁶。
原式の幅正値をdomain端・頂点で検査、curve域外描画falseを記録した。
表示用のdE/dx再binは0.1→0.5 keV/cm。R表示範囲のY flowも含むcount合計
43,006,175,713が前後一致し、元ROOT／旧study／source／製品の保全検査もPASS。

COLZでは中心関数は連続で、³Heラベルは低Rの主な帯の近く、⁴Heラベルは高側の肩・帯を辿る。
特に⁴Heについて、形の重なりだけで独立した⁴Heピークやそのpurityを確定したとは判断しない。
最新候補は未収束・未採用と図中にも記載した。滑らかな関数へつなぐ処理と表示は実装できたが、
**採用可能な校正式の確定は未完了**であり、本番PIDは今回変更していない。
