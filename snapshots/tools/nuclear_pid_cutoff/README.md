# TPC原子核PIDのpurity下限とR cutoff study

## 2026-10-19追加：元多項式studyと旧fitの品質監査

[当日Todo](../../Todo/2026/10/19/README.md)と[品質監査の実行手順](quality_20261019/README.md)を参照。
`logpoly_20261019_config_audit.py`は元のlog10(R)多項式に戻す独立study設定と切り戻しを検査。
`quality_20261019/audit_fit_quality.py`は旧集約の独立gate countを再現した上でχ²を比較する。
共通`PlotCutoffSpectra.C`の表示縦軸は、見えている範囲のraw count+errorと全モデルから設定する。
`PlotQualityDiagnostics.C`は旧holeと隣接Rの広域／局所・linear／log図を別名保存する。
本番設定の切替・新gate純度の再認証は、これらの診断完了と別の手順で行う。
[core＋tail比較の実行手順](coretail_20261019/README.md)に、追加許可された局所応答比較と数学試験をまとめた。

## 2026-10-17追加：高R TOFと本番採用

[当日Todo・統合PDF・新旧比較・本番検証](../../Todo/2026/10/17/README.md)と[高R再現手順](highrigidity/README.md)を参照。
旧TPC-onlyの80実行とPDFは保持。新本番設定は nuclearid_femtoLambda_auau3p85_purity50_production_20261017.yaml。
4核種の5本番入口が共通参照し、対応する prePurity20261017 のmainconfで切り戻せる。
prepare_production_ranges.py は凍結低R域＋高R全モデル合格表から採用域と理由を出力する。
combine_tpc_tof.py は明示settings/output-dirから新3頁＋旧3頁を生成し、入力SHAと全PDF頁を検査する。
render_native_atlas.py と highrigidity/render_highrigidity_atlas.py は全R・4nominalモデルの保存fitを失敗点も省略せず表示。
compare_pid_selections.py と compare_fitted_pid_parameters.py は新旧設定と実fitの中心・幅・差・比率を比較。
旧study用 purity50_20261015 は不変更。新本番の⁴He TOFはd同一窓であり、旧3<m²<4窓のpurityは未評価。
本番は検証済みR域だけを使用。fit幅は有限R投影幅で、校正係数へ直接採用していない。

以下は旧TPC-only studyと、その独立profile検証の記録。

対象はAu+Au FXT 3.85 GeV・2021 P24iyの既存部分統合ROOT。
低RではTOFを選択条件に使わず、`ALL`の実TPC ±2σ・best-species窓を評価する。
TOFタグは応答形の手掛かりに限り、TOF purityをTPC purityへ置き換えない。
今回の[Todo](../../Todo/2026/10/15/STAR_Todo_20261015.txt)に対応する後処理で、
新しいPicoDst解析・job投入・CF補正は行わない。

## 入力と保護

- 旧`nuclear_pid_tof_20261014/inputs_001/cache.root`を読み取り専用で使用する。
- R幅0.01、dE/dx幅0.1と元TH2F保存count・flowを保持し、数値fitと表示rebinを分ける。
- `freeze_baseline.py`で実PID係数・作業用ヘリウム校正・YAML・adapterを保存する。
- 今日別途許可された本番adapter／profile変更は旧baselineとの差として記録する。
  旧入力保護を無言で回避して新旧の設定を混ぜることはしない。
- 旧ROOT・PDF・fit・本番設定を上書きしない。出力はすべて新しいstudy名に保存する。

## ツール

| ファイル | 内容 |
|---|---|
| `study_settings.yaml` | 入力4組、R分割、数値条件、支持条件、0.5基準、感度・表示設定。 |
| `local_settings.yaml` | 核種ごとの局所窓、広窓での誤った背景制約を切り分ける数値guard、対象パラメータ自身の支持判定。 |
| `free_local_settings.yaml` | TOF-linearを初期値・帰属の参照とし、ALL内で対象の中心・幅をfitする比較。PID gateは固定。 |
| `freeze_baseline.py` | 現TPC選択の不変snapshotを作成。 |
| `run_cutoff.py` | Gaussian／log-normal、linear／exponentialの1組をROOT5で実行。 |
| `run_suite.py` | 宣言した複数組を限定したworker数で実行。既存出力へ再実行しない。 |
| `FitNuclearCutoff.C` | ALLだけの尤度で固定応答をfitし、実best-species窓の成分積分・観測count・残差を保存。 |
| `FitLocalNuclearCutoff.C` | d→⁴He全域の単一背景fitと分離した、対象核種局所窓のfit。各対象の結果を独立尤度として加算しない。 |
| `summarize_cutoff.py` | nativeの信号と全候補を合算して0.05幅のpurity、model envelope、支持範囲を出力。 |
| `prepare_views.py` | 保存済みfitの比較図に必要な明示TSVを作成。 |
| `ExtractNuclearGateCounts.C` | fit成否と独立に、元ALLと固定gateだけからnative候補count・bin境界の上下限を取得。 |
| `export_adopted_ranges.py` | 完了した感度比較から、採用候補の下限・明示veto・描画用範囲TSV・raw保持率を派生出力。 |
| `PlotNuclearCutoff.C` / `PlotCutoffSpectra.C` | purity／生COLZと固定gate、代表Rの全dE/dxスペクトルを描画。 |
| `test_cutoff_*.py` | variation、支持判定、集約、未評価・穴の保護の回帰試験。 |

## 読み方

TAGのnative投影幅はR内の変化を含むため、そのままintrinsic幅と断定しない。
Gaussianは平均／標準偏差、log-normalはmode／局所曲率幅であり、積分も各モデルを使う。
同じnative応答を3つのGL nodeに置き、固定TPC gateは実校正関数をnodeごとに評価する。
これは転用感度の比較であり、TOF受理の条件付けによるbiasが消えたという意味ではない。

失敗したGaussianの負の中心も元の値のまま保持するが、支持済み核種とは扱わない。
暗黙のlegacy置換、peakの強制生成、未評価をpurity=0に置き換える操作はしない。

purityは区間ごとに`sum(signalSelected)/sum(totalSelected)`で計算する。
Rごとのpurityの単純平均は使わない。下側envelopeは応答・背景・明示した転用感度の最小値で、
統計的な信頼下限や全系統誤差を包含する保証ではない。
固定応答のyield誤差は別記し、核種帰属・転用bias・共通のmodel誤記述を完全には含まない。

fitに失敗した旧数値出力の`rawGateCount=0`は、実測候補がゼロという意味ではない。
正式なfree-stage集約では独立gate extractorの640行を必須入力とし、
gateのSHA・native coverage・測定可能な各modelとのcount一致を検証する。
raw countは保存histのdE/dx bin部分重なりとRのGL3による近似で、
event-level再選別や真の再構成効率ではない。full-bin／touching-binの範囲も保持する。

支持には数値収束、対象peakの非境界・支持、gate内の局所／積分残差、
寄与のある未支持nuisanceの割合を別々に確認する。高統計χ²だけを理由に全点を無効にしない。
`unassessed`（未評価）と`below_purity_threshold`（支持された評価が0.5未満）を分ける。
支持された連続区間は`safeIntervals`へ保存し、未支持の穴をつながない。
単一の下限cutは、安全区間がTPC-only上端まで連続するときだけ`lowerCutoff`として出力する。

結果はinclusive accepted-eventのcharge-inclusive track purityであり、
Centrality別・Λとのpair条件付きpurityや最終CFの補正係数を直接与えるものではない。

### 全核種共通窓と対象局所窓

最初の全核種共通窓では、大量のd側countと全域の背景近似がrareな⁴Heへ影響し、
高Rで⁴Heのfit収量がゼロとなった。これは⁴He不存在の証明ではない。
exponential背景のnormalized-window slope上限8も、狭いTAG fit窓の数値guardを
広いALL窓へそのまま使うと強すぎる制限になり得る。
旧4実行は問題切り分けの基準として保持し、本番cutoffへ採用しない。

新しい局所設定では`--target d|t|he3|he4`を明示し、対象応答と実gateを含む共通窓を作る。
各対象の中では4モデル・全感度条件が同じ窓となる。他核種のGaussian／log-normal tailも保持し、
窓受理が明示した数値floor未満のものだけ非活性とし理由を保存する。
また、TAG_D内の⁴Heや背景が境界でもdの中心・幅まで境界と判定しない。
元TAGのnamed mean／log-widthと記録されたmath guardを対象ごとに照合し、
全rowの境界・適合度・数値収束と分けて保存する。

### ALLによる対象中心・幅のfit

固定TAG応答の転用でgate近傍の残差が残る場合、`targetShapeMode: free_target_mean_sigma`で
対象核種だけの中心／幅をALLで再最適化する。他核種と現PIDのgateは区別する。
TAG-linear応答は初期値・核種帰属の参照であり、ALLの背景linear／exponentialは両方比較する。
旧TAG-exponentialの数値境界を両方のALL背景の初期値へ持ち込まない。
TAG-exponentialとの差は削除せず診断へ記録する。

対象の許容中心shift・幅ratioはfitを押し込む制約ではなく、fit後の帰属screenである。
数値収束・shape covariance・対象非境界・screen通過をすべて区別して保存する。
このmodeのmean／width感度variationはnuisanceだけを変え、対象は再最適化する。
最適化が同じ解へ戻る単なる初期値変更を、対象の形状系統誤差と称しない。
TAG側の誤差だけでなく、ALLでfitした対象の中心誤差／幅と幅誤差／幅、対象窓収量も
設定した支持基準で確認する。Hessianが形式的に正でも誤差が大きい点は支持済みにしない。
native parameter表の`fitResponseQuality`はこの実ALL応答の確認であり、
最終の核種帰属・全モデルenvelope・残差判定を一つで代用するflagではない。

## 再実行

各CLIの`--help`で明示引数を確認する。数値設定はYAMLに置き、出力名は必ず未使用名にする。
runtime入口は`script/singularity_root_batch.sh MAINCONF`、STAR SL24y／ROOT5。
hostのROOTや暗黙mainconfは使わない。

```bash
python3 -B -m unittest discover -s tools/nuclear_pid_cutoff -p 'test_*.py'
```

全16組の新しいfree-target比較を実行する例（`NEW`は未使用suffix）：

```bash
python3 -B tools/nuclear_pid_cutoff/run_suite.py \
  --settings tools/nuclear_pid_cutoff/free_local_settings.yaml \
  --baseline rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/baseline_001/baseline.json \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --variation nominal --target d --target t --target he3 --target he4 \
  --suffix NEW --workers 2
```

全16組の成功後に、同じ設定・baselineを渡して集約する：

```bash
python3 -B tools/nuclear_pid_cutoff/summarize_cutoff.py \
  --settings tools/nuclear_pid_cutoff/free_local_settings.yaml \
  --baseline rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/baseline_001/baseline.json \
  --raw-gates rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/independent_gate_counts_001/gate_counts.tsv \
  --runs rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/all_{d,t,he3,he4}_*nominal_NEW \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/summary_free_NEW
```

wide／fixed-local／free-localは設定ファイルで明示的に分ける。既存runの途中へ
別設定を足したり、窓やgateの違う実行を一つのenvelopeへ混ぜたりしない。
各runにcommand・設定・source snapshot・完了manifestを、各summaryにも派生成果の
source snapshot・入力SHA・生成物SHAを保存する。

本studyの実行完了と、cutoff／新PID関数の本番採用は別の判断である。

## 2026-10-15の集約版

`nuclear_pid_cutoff_20261015/summary_local_nominal_002/`が局所窓16実行の修正済み集約。
`summary_local_nominal_001`では、各対象の4モデル数を全16実行数と誤って比較して
全点を未評価にしていたため、支持フラグの判定に使用しない。元fit値は不変。
回帰試験に4核種×4モデルの独立性と、高Rの穴を無視して下限を認証しない試験を追加した。

最終の数値入力は`summary_free_sensitivity_001/`。
nominal 16組とnuisance感度64組の計80組を、全R・全対象について比較している。
geometry・実gate・fit窓・元ROOTを固定し、全runのROOT／表／数値sourceのSHAを確認した。
各対象内の20条件で利用可能な値がそろわない区間を、残りの都合のよいモデルだけで採用しない。

範囲を採用する場合は`adopted_ranges_001/adopted_ranges.json`の`adoptedMin`と
`adoptedVeto`を組にする。`summary.json`の単一suffix用`lowerCutoff`を、穴を許容した
最小Rと取り違えない。描画用は同フォルダの`adopted_ranges.tsv`。
R>=2のTOF条件は別profileの設定であり、このTPC-only studyから高Rのpurityを推定しない。

`numerical_recovery_settings.yaml`によるt／⁴HeのGaussian+exponential全R再実行は、
AMP精度を厳しくしても既存の数値穴を解消せず、追加の無効点が増えた。
これら`recovery001`は診断として保持し、最終80条件へ混ぜたり、点ごとによい方へ置換したりしない。

## 2026-10-17（日本時間）の再開・最終表示と設定検証

前回完了したfitの再実行は不要。最終表示は
[views_final_001](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_final_001/README.md)
で、purity・COLZ・候補数の3ページ、応答比較5ページ、核種ごとの代表fitを保存した。
`PlotNuclearCutoff.C`は省略可能な第8引数に使用範囲TSVを受け取り、穴をつながず除外域と使用域を示す。
`PlotNuclearParameters.C`はnominalの実ALL応答検査を通った点だけを表示する。
`PlotCutoffSpectra.C`は保存したfit countをそのまま表示し、再fitや表示規格化は行わない。
未支持fitを診断として描いたことと、purityの使用候補と認めたことを区別する。

[最終mainconf](../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_purity50_20261015.yaml)と
[核PID YAML](../../config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_20261015.yaml)は、
保存した最終下限とvetoを実解析へ明示的に適用する独立プロファイル。
後処理fitとは別に、同mainconfでbatch一致build・2021ローカル入力の先頭20イベント・閉ROOTを検証した。
[検証記録](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/purity50_20261015/README.md)を参照。
ここでのイベント試験は接続確認だけであり、新たなpurity・効率測定や本番job投入ではない。
旧本番mainconf・中心／幅の係数は変更せず、baselineを新しい選択窓へ黙って置き換えていない。
