# 0.40開始の区間別fitから共通PID式へつなぐstudy

旧`refined_lowrigidity_from040`の保存済みGaussianを初期値として、まず区間ごとの
中心・幅を自由にfitし、その結果を段階的なsoft penaltyで滑らかな関数へつなぐ。
旧legacy中心・幅に相対的な制約や振幅capは引き継がない。
元のPDF・ROOT・現在の解析PIDを上書きしない、postprocessing専用の入口である。

[計画](../../../Todo/20261013/task3_helium_iterative_from040_plan.md)と
[結果・最新状態](../../../Todo/20261013/task3_helium_iterative_from040_report.md)を参照。

## 実行

入力は保存済み3.85 GeV・2021 P24iyのみ。`mainconf`はROOT5/SL24yの選択に使い、
現在のPID cutをヒストへ再適用するものではない。
全条件は`iterative_settings.yaml`にある。新しい、まだ存在しない出力先を指定する。

```bash
python3 -B tools/he4_purity/iterative_logpoly/run_iterative_fit.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/NEW_STUDY
```

先頭10区間の自由fitだけを先行確認する場合：

```bash
python3 -B tools/he4_purity/iterative_logpoly/run_iterative_fit.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --indices 0,1,2,3,4,5,6,7,8,9 --stop-after-free \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/NEW_PROTOTYPE
```

`--prepare-only`は来歴・設定・初期値の保存だけ、`--no-render`は数値fitだけ。
`--indices`は元cacheのindexを保つ部分試験であり、全域を観測した校正とは区別する。
出力済みstudyは再開・上書きせず、新しいstudyに実行する。

## 計算と表示の意味

- 生の0.1 keV/cm射影を保存し、0.5 keV/cm fit binとのcount/error/flowを照合する。
- 固定した旧fit窓、旧採用linear/exponential背景modelを全phaseで用いる。
  fit窓外には別の核種もあるため、広い表示範囲とfit対象を区別する。
- Gaussian templateはfit窓で正規化し、fitする振幅は**非負の窓内count**。
  旧TF1のGaussian高さは`H σ_raw sqrt(2π)/Δdedx`で全countへ換算してから、窓acceptanceを掛ける。
- 区間別free/regularizedの幅は投影幅。最後に共有14係数でR-bin平均を含むraw likelihoodを再fitし、
  intrinsic幅へ再最適化する。投影幅をそのままintrinsic幅として二重に広げない。
- 中心3次・幅2次の`log10(R)`式は解析コードと同じ形式。数値計算では中心化・正規化した基底を使う。
  幅の正値は端点・二次式の頂点で確認し、Bernstein係数の非負に関数familyを限定しない。
- penalty情報量は`sum(component_i^2 / total_i)`、scaleはfree幅から固定する。
  これは混合model内の重みで、核種の確定・独立測定の誤差ではない。
  同じ重み・中心/幅の差の二乗でQRと局所penaltyを対応させる。
- `fit_result.json`はraw data devianceとpenaltyを分離。実行完了、局所収束、共有係数の収束、
  物理的な核種の検証は別の状態として保存する。
- 低Rの肩、zero-yield、入れ替わり、非識別・低品質は隠さない。failed binでlegacyへ突然戻す関数は作らない。

各phaseの`views/`には`fits_nsigma_linear.pdf`、`fits_raw_dedx_linear.pdf`、
`low_rigidity_components.pdf`を保存する。rawとnσは同じbin countを表示し、y軸のJacobian変更はしない。
最後の`analytic_candidate.yaml`はstudy用で、本番設定へ自動適用しない。

## 数値テスト

```bash
python3 -B -m unittest discover -s tools/he4_purity/iterative_logpoly -p 'test_*.py'
```

ROOT側の`SelfTestIterativeHelium.C`と`TestIterativeRenderer.C`は、明示したSL24y mainconfで
`script/singularity_root_batch.sh`からACLiC compileして実行する。
既存のfit・PDF・元ROOTの不変と、使用したソースsnapshotをstudyごとに記録する。

## 共通fitだけの数値・背景初期値の感度

最初の全域fitはAMP判定による目的関数の急変と、低Rの背景分担の変化を確認したため未採用。
`run_common_refit.py`は完成済み全210区間studyを検証し、元studyを変更せず新しいsiblingへ再fitする。
同じ係数から、`--nuisance-start free`と`last_regularized`を比べられる。
`--fix-triton-shape`はfreeデータのT中心・幅に条件付ける診断で、legacy式への固定ではない。

```bash
python3 -B tools/he4_purity/iterative_logpoly/run_common_refit.py \
  --study-folder rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/full_iterative_20261014_001 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --nuisance-start free --coefficient-start final \
  --settings tools/he4_purity/iterative_logpoly/common_refit_numerical_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/NEW_COMMON_CHECK
```

数値override YAMLは厳密schemaで、物理cut/modelの変更を拒否する。
`--prepare-only`は照合・snapshot・初期値準備まで。`--no-render`は数値結果だけ保存する。
fixed-Tの結果・誤差はfreeのT modelに条件付きであり、その固定自体で核種同定を証明しない。

`diagnose_common.py --study STUDY --phase 13_common_raw_refit --mainconf MAINCONF
--output-dir NEW_DIAGNOSTIC --variant strict:1e-8:250 --variant loose:1e-6:1000
--timeout-seconds 3600`は、保存された最終点と±14座標で、数学的support・template・AMP profileの失敗理由を分離する。
独立診断の精度・計算数と実行例はそのsource／出力`settings.tsv`、最新Todo報告に記録する。

## 低剛性domainだけの同一関数形の診断

`run_low_domain_common.py`は、同じ全域free結果の先頭30区間（R=0.40〜0.70）だけを取り出す。
中心3次・幅2次の形式は変えず、低Rのfree点のQR係数から共通raw fitを開始する。
`low_domain_common_settings.yaml`は診断domainと数値条件だけを定義する。
解析の入力範囲・PID cutを変更する入口ではなく、本番へpiecewise校正を適用しない。

```bash
python3 -B tools/he4_purity/iterative_logpoly/run_low_domain_common.py \
  --study-folder rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/full_iterative_20261014_001 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --domain-source rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/domain_diagnostic_20261014_001/domain_comparison.json \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/NEW_LOW_DOMAIN_CHECK
```

`--fix-triton-shape`あり／なしを別studyへ保存する。比較基準は全域Dではなく同じ30区間のfree D。
全原本・source・入力fingerprintと各区間のfit窓・countを前後に照合する。
各軸のlinear図8ページと、0.40〜0.50の4面図10ページを保存する。
診断が実行完了しても、収束・核種同定・本番の適用可否は別に判定する。
`low_domain_convergence_settings.yaml`は同じ初期値・modelで反復上限だけ12から120へ増やす診断。
上のコマンドへ`--settings tools/he4_purity/iterative_logpoly/low_domain_convergence_settings.yaml`を追加し、
別の新しい出力先を指定する。物理条件を緩めたり、未収束の候補を本番へ適用する設定ではない。

## 表示専用: 元dE/dx全範囲とCOLZへの関数重ね

`draw_wide_iterative_helium.py`／`WideIterativeHelium.C`は、完成したstudyの保存済みphase ROOTを読む。
`wide_iterative_settings.yaml`で全保存軸0〜120 keV/cmと低Rの45〜120補助範囲を指定し、
rawとnσの同じbinをlog/linearで表示する。fit窓は変えず窓外はデータのみ。
元のimmutable rendererをsourceへコピーして使い、新fit・PID変更は行わない。

```bash
python3 -B tools/he4_purity/iterative_logpoly/draw_wide_iterative_helium.py \
  --study rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_domain_convergence_common_20261014_001 \
  --phase-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/low_domain_convergence_common_20261014_001 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --settings tools/he4_purity/iterative_logpoly/wide_iterative_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/iterative_from040/NEW_WIDE_VIEWS \
  --phase-label 'COMMON / UNCONVERGED / UNADOPTED'
```

`draw_low_domain_raw_overlay.py`／`DrawLowDomainRawOverlay.C`は、元`true/hDedxP`のCOLZへ
保存済み係数のTF1を`Draw("same")`で重ねる。表示条件は`low_domain_overlay_settings.yaml`。
`--latest-study`、`--reference-study`、`--qr-source`、`--mainconf`、`--output-dir`を明示する。
今回の最新studyは上の120反復、referenceは`low_domain_free_common_20261014_001`、QR sourceは
`domain_diagnostic_20261014_001/domain_comparison.json`。同じ入力の完成study以外は拒否する。
4ページPDFに全域COLZ・核領域COLZ・幅・中心比較を保存し、評価domain0.40〜0.70の外へは描かない。
±σはGaussian intrinsic幅であり、係数誤差・核種同定・purityを保証する帯ではない。

`draw_tof_compatible_helium.py`／`DrawTofCompatibleHelium.C`は`tof_compatible_settings.yaml`に従い、
元ALLとTOF³He／Dタグの未切断分布を対照表示する。Dタグは⁴Heも通るm²/q²域だが、純粋な⁴He標本ではない。
これも表示専用で、fit・PIDは変更しない。全保存軸図と核領域図を別PDFに保存する。

全入口は新しい出力先だけを使い、元のstudy製品／source SHA・元ROOT fingerprintを前後に照合する。
広範囲図の補助PNGはROOT5 exportで文字がclipする場合があるため、確認済みPDFを参照する。
図・数値・採用判断の最新リンクは[Todo報告](../../../Todo/20261013/task3_helium_iterative_from040_report.md)にある。
