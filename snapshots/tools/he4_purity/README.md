# ⁴HeのTPC純度推定

## 最新の依頼: 0.40開始free fitから反復して共通式へ

旧`refined_lowrigidity_from040/views/he4_low040_components.pdf`は存在する。
これを初期値に、区間別自由fit→soft平滑化→共有14係数のraw/R-bin平均fitを
新しい`iterative_logpoly/`と別studyに実装する。legacyからのcold global fitは今回の出発点にしない。
nσだけでなく元dE/dxで、全区間のlinear図と低Rの4面log/linearを保存する。
旧成果物と本番PIDは保持する。
全域の共有式は低Rの成分再現が悪化したため未採用。`iterative_logpoly/run_low_domain_common.py`で
同じ関数形・同じfree点のR=0.40〜0.70診断を別保存し、domainと背景分担を切り分ける。

[実行方法](iterative_logpoly/README.md)と
[最新報告](../../Todo/20261013/task3_helium_iterative_from040_report.md)を参照。

## 再検証中: 低剛性の中心・幅を共通係数で求める

現行の点ごとの補間校正は最終校正ではない。2026-10-13の再確認で、
R=0.410 GeV/cのlegacy値への切り替え境界に、⁴He中心の約22.7 keV/cmの段差を確認した。
元の校正・fit・PDFは保持し、生のp/q×dE/dx分布との比較と、
解析コードと同じ`log10(R)`の中心3次式・幅2次式の共通係数fitを別studyとして進めている。
以下の「導入」は前回の作業用校正の記録であり、低剛性の正しさを保証しない。

[再校正計画](../../Todo/20261013/task3_helium_logpoly_recalibration_plan.md)と
[再校正・検証結果](../../Todo/20261013/task3_helium_logpoly_recalibration_report.md)を参照。

## 最新: ³He/⁴He校正の導入と選択窓の純度比較

2026-10-13の明示的な依頼により、0.40開始の採用205点の中心・Gaussian幅を
3p85 / 2021 P24iyの現行⁴He解析へ導入した。共有の多項式係数は変更せず、
専用の`config/calibration/dedx_helium_auau3p85_2021P24iy_from040.yaml`を使う。
旧gate、新gate、新しいTPC best-speciesの比較は
`pdf/he4_purity_20261013/calibrated_selection/purity_calibrated_windows.pdf`へ別保存した。
このYAMLは多くの参考fit点を含む作業用校正で、TOF後・ペア別の確定純度ではない。

校正exportと保存済み背景modelの再積分を再現する入口は、以下。
出力済みYAML・QAは上書きしないので、再実行は必ず新しい名前を指定する。

```bash
python3 tools/he4_purity/export_runtime_calibration.py \
  --fit-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040 \
  --calibration-output UNUSED_CALIBRATION.yaml \
  --output-dir UNUSED_QA_DIRECTORY \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml
```

fitの再実行はしない。校正YAMLの補間・支持half-bin・legacy復帰を使う同じ窓で、
linear/exponentialの2背景modelを積分する。保存されていない新windowの共分散誤差は作らない。
`check_runtime_calibration.py`の独立検証は、選択区間・数値積分と実C++値を照合する。
最終結果・実装箇所・切り戻し・³He CF入力・再解析の必要性は
[Task3実装報告](../../Todo/20261013/task3_helium_pid_cf_implementation_report.md)を参照。

## 保存済みヒストの元fit workflow

既存の未PID `true/hDedxP` を剛性で細かく分け、校正nSigma±2内のGaussian信号と背景の積分比から純度を推定します。イベント解析やjob投入は行いません。

```bash
python3 tools/he4_purity/run_purity.py --sample 3p85 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml
python3 tools/he4_purity/run_purity.py --sample 13p5 \
  --mainconf config/mainconf/main_auau13p5_anaFemtoLambda_4He_KFParticle.yaml
```

設定は `purity_settings.yaml` にあり、入力は既存summaryとfingerprintを照合します。ROOT処理は `script/singularity_root_batch.sh` のROOT5環境を使います。出力済みのフォルダは上書きせず、再実行する場合は `--output-dir` で別の場所を指定します。

`rootfile/<anaName>/pdf/he4_purity_20261013/` に剛性対純度（1ページ）、linear/exponential背景モデル比較・差（2ページ）、全射影fit（53ページ）、CSV、監査JSON、設定・ソースsnapshotを保存します。fit図は元dE/dx射影を線形変換した校正nSigma軸で表示し、元のkeV/cm射影もfit ROOTに保存します。正負電荷を合算したTPC窓の推定であり、最終best-speciesやTOF条件後の純度とは区別してください。品質基準を通らないfitは参考値として区別します。

旧設定の結果では、13.5 GeVの品質通過区間は168/210ですが、候補数coverageは34.94%です。3.85 GeVの旧fitは品質基準を通る区間がありませんでした。追補ではdE/dxの保存精度によるbin間の振動と低剛性の成分制約を別々に再検証します。参考fitのcovariance誤差だけを純度の精度として使用しないでください。

[実装計画](../../Todo/20261013/task3_he4_purity_plan.md)

[結果・出力リンク・解釈の限界](../../Todo/20261013/task3_he4_purity_report.md)

## 低剛性ピークの再探索とlinear表示

低剛性の信号平均の移動範囲・複数初期値・近傍核種のtailを見直すstudy設定は
`purity_lowrigidity_settings.yaml`です。元の設定・元の±2校正窓と比較して使います。
`dedxRebinFactor: 5`でdE/dxの0.1 keV/cm binを0.5 keV/cmにまとめ、Float16保存の離散格子とbin幅による振動を抑えます。Rの0.01 GeV/c bin、全count、元窓の観測countは維持します。元の射影は`raw_native_*`、fit用は`raw_*`と`projection_*`として保存し、元の校正窓の純度とfitされた中心・幅の窓の純度を区別します。キーを省略する旧設定ではrebinしません。

```bash
python3 tools/he4_purity/run_purity.py --sample 3p85 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --settings tools/he4_purity/purity_lowrigidity_settings.yaml
```

改訂結果は元の出力内`refined_lowrigidity/`に保存します。
全射影の`he4_fits.pdf`（log）に加え、`he4_fits_linear.pdf`（linear）を作ります。
13p5はsample・mainconfを13p5へ変更します。非空の出力先は上書きしません。

保存済みfitを再fitせず表示方式だけ追加する場合は、以下を使用します。

```bash
python3 tools/he4_purity/render_fit_views.py \
  --study-folder rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml
```

`views/he4_fits_linear.pdf`に全区間のlinear図、
`views/he4_low_rigidity_dedx.pdf`に低剛性の物理dE/dx軸（keV/cm）でlog/linearを並べた図を保存します。
表示設定は`view_settings.yaml`。元のfit ROOT・PDF・数値は変更せず、曲線の座標変換と
保存binのcount/errorの一致を検査します。
`No accepted fit`は不採用の意味で、ピーク不存在の判定ではありません。

[再探索の計画](../../Todo/20261013/task3_lowrigidity_refinement_plan.md)

[再探索の結果・低剛性の表・解釈](../../Todo/20261013/task3_lowrigidity_refinement_report.md)

新studyは`he4_locus_fit_candidates.csv`に元の窓とfit中心・幅の窓を別列で保存し、`he4_pid_curve_candidate_<sample>.yaml`にstudy用の曲線候補を出します。productionへ自動適用しません。候補曲線の小さい最小二乗計算は標準ライブラリで行い、hostのNumPyやlibgfortranには依存しません。

`--finish-existing`は、ROOT fitが正常終了した後に候補保存が失敗したstudyだけを検査して後処理を完成させる修復入口です。入力・設定・ソースsnapshot・fit完了ログが一致し、未完了の場合に限って使えます。旧publisher・contextを`revisions/postprocess_repair_001/`へ保持し、fitのROOT/PDF/数値は再実行・上書きしません。すでに`completed.json`がある出力には使えません。

## 0.40 GeV/cからの再fit

3.85 GeVだけを対象に探索開始を0.40へ下げる設定は`purity_lowrigidity_from040_settings.yaml`です。0.50未満の信号平均・初期値を別キーで指定し、0.50以上と元のPID窓は維持します。背景にはtritonと³Heの2つのGaussianピークを含めます。

```bash
python3 tools/he4_purity/run_purity.py --sample 3p85 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --settings tools/he4_purity/purity_lowrigidity_from040_settings.yaml
```

出力は`pdf/he4_purity_20261013/refined_lowrigidity_from040/`。最初の0.40–0.41区間は中心R=0.405で評価し、fit高端は保存されたdE/dx上限120 keV/cmで制限します。範囲端に張り付く信号は採用せず、数値的なfitと純度の品質を区別します。既存の出力先は上書きしません。

試験用の`run_lowrigidity_prototype.py`は`--indices 0,1,2,3,5,7,9,10,20`でcached射影の対象binを指定でき、`--cache PATH`で入力cacheを明示できます。既に要求幅へrebinされたcacheを再度rebinしません。指定したindicesとcacheのhashを`context.json`へ記録します。

表示は`render_fit_views.py`に`--settings tools/he4_purity/view_settings_from040.yaml`を指定します。低剛性の全分布と⁴He候補付近を別panelにしてlog・linearで描き、triton・³Heの背景成分も個別表示します。

[0.40開始の計画](../../Todo/20261013/task3_from040_refit_plan.md)

[0.40開始の結果と採用可否](../../Todo/20261013/task3_from040_refit_report.md)。新`views/he4_low040_components.pdf`は10ページ・各4面（全分布と候補付近のlog/linear）で、実fit範囲と保存上限を明記します。`views/he4_fits_linear.pdf`は全210区間・53ページです。最初の0.40–0.41の曲線は診断表示で、採用した純度とは区別します。

### 下段の物理dE/dxも全体を表示する

`view_settings_from040_full_dedx.yaml`では、上下とも従来nσ=-8〜+4に対応する
同じ範囲を表示する。triton/³Heの中心が下段に入ることと、上下のnormal bin・count・errorを検査する。
表示だけの変更で再fitはしない。旧viewsは上書きしない。

```bash
python3 tools/he4_purity/render_fit_views.py \
  --study-folder rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --settings tools/he4_purity/view_settings_from040_full_dedx.yaml \
  --output-dir UNUSED_VIEWS_DIRECTORY \
  --allow-runtime-mainconf-drift
```

最後のオプションは、保存studyと同じpathのmainconfがその後変更された場合だけ明示する。
元studyのfit/source/settings/inputの検証は省略せず、旧・現行SHAとruntime設定snapshotを保存する。
現行設定はROOT5の実行環境の指定だけに使い、保存fitの解析条件には置き換えない。
オプションなしの通常動作ではSHA差異を拒否する。
今回の完成出力は`refined_lowrigidity_from040/views_full_dedx/he4_low040_components_full_dedx.pdf`。

## 生分布への曲線重ね描画

`draw_raw_calibration.py`は元ROOTのpre-PID `true/hDedxP`とTPC選択に依存しない
TOF³Heタグ分布を読み、生COLZ・射影・ALL/タグ形状比較を新studyへ保存する。
表示設定は`raw_overlay_settings.yaml`で、margin・字体・captionの位置も記録する。
表示を変更しても既存studyを上書きしない。

```bash
python3 tools/he4_purity/draw_raw_calibration.py \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --calibration config/calibration/dedx_helium_auau3p85_2021P24iy_from040.yaml \
  --tof-nuclear-config config/cuts/nuclearid/nuclearid_femtoLambda_legacy.yaml \
  --settings tools/he4_purity/raw_overlay_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/UNUSED_RAW_VIEWS
```

最新の表示だけを改善した出力は`raw_overlay_readable_20261014/`。
元の数値・曲線・診断窓・subset監査は`raw_overlay_baseline/`と同じで、元ROOTは変更しない。
TOFタグ総数は³Heピークの統計ではない。タグはALLの部分集合で、規格化後の高側差を
そのまま⁴Heの収量やCF背景係数とみなさない。

任意のstudy-only係数候補は`--analytic-coefficients CANDIDATE.yaml`（旧alias
`--candidate-calibration`）で追加表示できる。schemaVersion=1、functionalForm=log10_polynomial、
productionEnabled=false、candidateOnly=true、単位GeV/c・keV/cm、outsideRangePolicy=legacyを要求する。
species.he3/he4はdomain、meanCoefficients（4係数）、sigmaCoefficients（3係数）を持つ。
表示したことは採用・PID検証を意味しない。今回の未収束prototype係数は表示・採用していない。

## 複数R区間で共有する係数fit

`global_logpoly/GlobalHeliumFit.C`は、³He/⁴Heの中心3次式と幅2次式の合計14係数を
全R区間で共有する。dE/dx binにはGaussian積分、R binにはGauss–Legendre平均を使い、
Rごとの規格化・背景・近傍核種のshapeだけをnuisanceとして扱う。
各区間に独立したHe3/He4の中心・幅を持たせない。

```bash
python3 tools/he4_purity/global_logpoly/run_global_fit.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/global_logpoly_recalibration/UNUSED_STUDY
```

`global_logpoly_settings.yaml`は少数20区間のprototype設定で、最終校正ではない。
`--prepare-only`は実fitせず入力・初期値と来歴を保存する。`--seeds legacy legacy_low_shift old_knots_poly`
は初期値感度の比較で、旧独立fit点を真値として拘束しない。
feasibleなlegacyとのblend率も保存する。既存outputは空でも拒否する。
元ROOT・summary・sourceと実際にcompileするsnapshotの前後一致を検査する。

`completed.json`のexecutionCompleteはfit収束・核種識別・本番PIDの検証完了ではない。
JSONのAMP/KKT・local/global/outerとstationarity、conditional covarianceの評価点を別々に確認する。
Gaussian amplitudeは全線収量／fit窓の観測countで、窓内fractionや総和1の制約ではない。
conditional covarianceは近傍核種・非線形背景shapeを固定しており、完全な系統・校正誤差ではない。
[最新の数値・実データ検証](../../Todo/20261013/task3_helium_logpoly_recalibration_report.md)を参照。
