# 4Heの背景なしfitと3He低側への範囲拡張

最新の依頼に従い、**実際のfit下限を³Heの低側半値まで広げ、背景なし・linear背景・exponential背景を同一範囲で再fitした**。Gaussian/lognormalの各160区間、R=0.4–2.0 GeV/c、ALL試料を対象とする。既存の最終Purityと本番設定は保持。

## 図と数値

- [代表3区間×2応答、6頁の比較PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/no_smooth_he4_low_he3window_001/representative_fits.pdf)
- [Gaussian全160区間](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/no_smooth_he4_low_he3window_001/fit_atlas_gaussian_he4.pdf)
- [lognormal全160区間](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/no_smooth_he4_low_he3window_001/fit_atlas_lognormal_he4.pdf)
- [各区間のχ²・NDF・deviance・Purity・品質フラグ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/no_smooth_he4_low_he3window_001/metrics.tsv)
- [0.05 GeV/cごとのΣS/Σtotal集約と採否](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/no_smooth_he4_low_he3window_001/aggregate.tsv)
- [範囲と半値位置の根拠](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/no_smooth_he4_low_he3window_001/fit_windows.json)
- [全比較の機械可読記録](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/no_smooth_he4_low_he3window_001/comparison.json)

図は上から背景なし、linear、exponentialの3段、左linear・右log。共通の広い表示軸で生データを示し、fit曲線は実fit範囲内のみ描画。各パネルにfit範囲、χ²/NDF、deviance、numerical・identity・historical qualityを記す。灰色の縦点線がfit下限。背景なしのsmoothは厳密にゼロ（描画用ゼロ成分だけ残る）。

## Fit下限の規則

各Rで、従来の4He fitに含まれる固定³He参照と、数値・free-response支持条件を通る既存ALL ³He fitを参照する。Gaussian/lognormal双方の低側半値位置と旧下限の最小値を取り、0.5 keV/cm刻みで切り下げた。追加の経験的余白は0。Gaussianはμ−√(2ln2)σ、lognormalはmode×exp[−√(2ln2)width/mode]（保存ソースの定義）。参照選択を含め全値をfit_windows.jsonに保存。

160区間のうち143区間を拡張、17区間は元から十分広く維持。上限不変更。たとえばR=0.805では38→34.5、1.205では25.5→22.0、1.605では19→17.0 keV/cm。全モデルが同じ下限・上限を使用し、すべての参照半値位置を含む。

## 背景除去の比較結果

smooth backgroundのテンプレート列・振幅・傾き自由度を除き、d/t/³Heの核種成分は残した。背景なしにするとPurityはほぼ1になる区間があるが、背景をモデルから除いた結果でもあり、清浄な試料である根拠にはならない。代表区間でもχ²/NDFが数十〜数百と大きく、絶対的なfit品質には課題が残る。既存の支持フラグは従来条件を適用したもので、χ²が良いという意味ではない。

| 応答 | 対照背景 | 数値比較可能区間 | 背景なしでχ²/NDF改善 | χ²/NDF比の中央値（なし/あり） | 双方identity通過で改善 |
|---|---|---:|---:|---:|---:|
| gaussian | linear | 151 | 0 | 1.656 | 0/136 |
| gaussian | exponential | 151 | 0 | 4.396 | 0/136 |
| lognormal | linear | 150 | 1 | 1.535 | 0/134 |
| lognormal | exponential | 149 | 1 | 2.610 | 1/134 |

ほとんどの区間で背景なしは悪化する。lognormalではχ²/NDFが改善する少数例もあるが、全体として背景除去を支持する結果ではない。Purityの最終採用値には反映していない。

背景ありは従来初期値と背景なしの最適化結果を用いた初期値の双方からfitし、数値的に成立する解のうちdata devianceが小さい解を比較に使用した。ピークidentity判定の基準は元のまま保持。ゼロ背景を含むモデルのdevianceが背景なしより悪くなる入れ子関係違反（許容0.01）は全比較で0。Pearson χ²・χ²/NDFはdevianceとは別に集計し、NDFは背景パラメータを除いた数を使用。

旧範囲の試行no_smooth_he4_low_001は診断用履歴として保持し、本報告の比較表には使っていない。異なるfit範囲間のχ²の大小を、そのまま背景効果とは扱わない。

## 実装・検証と再現

- 実装・設定：tools/nuclear_pid_smooth_20261021/no_smooth/{run.py,analyze.py,PlotNoSmooth.C,render.py,ValidateNoSmooth.C,finish.py,settings_he3window.yaml}。本番ヘッダは変更せず、10/17保存コードの私有コピーに背景なし分岐を加えた。
- 拡張範囲で1,600 fit：背景なし320、従来初期値の背景あり640、背景なし初期値の背景あり640。両応答で各160区間を完走。
- 事前検証：旧範囲の代表6fitは元のyield・形状・Pearson・devianceを相対1e-8で再現。Gaussian/lognormal既知収量toyで背景なしモデルの回収、4核種列のみ、smooth=0、NDFを確認。
- 今回のROOT検証：背景なし全320区間でsmooth全binが0、total=核種成分和。1,280対照でraw_nativeの全bin内容・誤差が完全一致。全モデルでPID gateのSHA、fit範囲を検証。
- PDF全320頁と代表6頁をGhostscriptで読み戻し、頁数とSHAを保存する。元の数値・PDF・表示修正版・最終PurityのSHA不変も検証する。最終証跡はcompleted.json。
- mainconf=config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml。script/singularity_root_batch.sh経由のSTAR SL24y/ROOT5。jobid/configlogなし（ローカルstudy）。job投入・hadd・commit/pushなし。
- 初回同時ACLiCビルドが衝突したため失敗ログと生成物を保存し、単独ビルド後に全fitを実行。runnerは初回ビルドを直列化するよう修正。計算済みPDFの検証は、この環境にあるGhostscriptに統一した。

再現手順：新しい空出力先をYAMLに指定し、run.pyのprepare→none→baseline→warm、analyze.py analyze、render.py render、ValidateNoSmooth.C、finish.py finishの順。Python各コマンドの第2引数は明示的なsettings_he3window.yaml。コマンド・実効設定・私有ソースと差分・最終runnerを出力ディレクトリ内へ保存。

## 最新追加：4He R>0.61・低Rピーク位置

[更新fit図・Purity・実施報告](he4_r061_20261021.md)。ユーザー指定の³He/⁴Heピークを低R参照へ反映し、4He適用下限をR>0.61にした独立study。0.70以上は直前の半値範囲fitと数値差ゼロ。旧成果物は保持。最新の4He統合Purityはhe4_r061_001内に保存している。
