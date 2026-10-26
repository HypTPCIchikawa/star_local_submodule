# 2026-10-21：Git元の滑らかな核PIDによるPurity再評価

完了：候補16,640スライス、元係数比較3,400スライス、代表24fitを実行。統合図・COLZ・全R fit集の計1,610ページを全頁読み戻し、元ROOT・本番設定・過去入力の保護も確認した。

## 成果物

**Purity図の最新版：** [1ページ目の核種別横軸を調整したPDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/purity_axis_002/purity_vs_p_over_q.pdf)。dは維持、2ページ目もそのまま。

**1次元fit図の最新版：** [表示幅2倍の002版・全8冊](fit_atlas_wide_20261021.md)。001版も保持。

- [最終Purity vs p/q PDF（全域＋高R、2頁）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/purity_vs_p_over_q.pdf)
- [最終Purity PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/purity_vs_p_over_q.png)
- [ALL／各TOFと元・候補のPID中心／±2σ、8頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/dedx_vs_p_over_q_curves.pdf)
- [全Purity・S/総数・条件付き誤差・モデル別採否JSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/candidate_summary.json)
- [Purity一覧TSV](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/candidate_summary.tsv)
- [元／候補係数と調整根拠YAML](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/coefficients.yaml)
- [独立study結果YAML（本番無効）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/study_results.yaml)
- [支持され50%以上の区間JSON（自動vetoではない）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/assessed_ranges.json)
- [Git元ヘッダと旧baselineのSHA監査](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/git_original_audit.json)
- [最終検証manifest](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/completed.json)

## 方法と10/17の再現

- 対象はAu+Au FXT 3.85 GeV、2021 P24iyの部分統合ROOT由来。charge-inclusive、centrality積分のinclusive track QA。全23B統計・Λ pair・k*条件付きPurityではない。
- 入力は10/14保存cacheのALL=true/hDedxP、TAG_D=true/hDedxP_d_m2、TAG_T=true/hDedxP_t_m2、TAG_HE3=true/hDedxP_3He_m2。元ROOTは10,522,646,853 bytes。保存TH2F count・error・flowとnative軸を保持し、Entriesへの規格化や整数丸めを行わない。
- 元投入時のtrack QA選択を継承し、現在のMakerから条件を再定義しない。来歴は[元cache context](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001/context.json)と元jobid 0BEFC64E62F19BBA534A129CABC9B106のconfiglogを参照。
- 10/17の保存C++・header・shapes・fit窓・settingsを新出力へ複製。全4核種×低R/高R×3代表点=24fitでS、総数、raw count、Purityの96比較が相対1e-8以内で一致した。
- 元のPIDヘッダは`git show 8f7d539:StMaker/common/NuclearIdDeDxVsMom.h`から直接取り出し、旧baselineのSHAと一致した。現在の解析ソースや補間校正は計算に使用していない。mainconfはSTAR SL24y/ROOT5の環境選択に使用した。
- 10/17の低R80実行、高Rのd/t/3He native計12実行と4He coarse4実行を同じ数学で再実行。既存96実行13,640スライスではgate以外のfit窓・最適化パラメータ・yield・適合指標が完全一致した。d/tの大幅不一致を式の限界へ帰属していない。

## 元の関数形と微調整

R=|p/q| [GeV/c]、x=log10(R)。中心はa0+a1*x+a2*x²+a3*x³、幅はb0+b1*x+b2*x²のまま。

- d/tの係数はGit元のまま。
- Heは旧Gaussianの支持されたALLおよび一致するTOFのtarget-freeピークを参照。R=0.8–2.5で元係数との差を同じ次数へ回帰した。Gaussian標準偏差を使い、10/20のsplit-normal幾何平均幅を移植していない。
- 高統計の一点が支配しないよう、元σで規格化した残差へ等しい重みを与えた。全0.4–5の1,001点で中心の変更を0.5元σ以内、幅の変更を20%以内へ縮小した。
- 4核種とも中心・幅は正。中心に内部極値はなく、元の式に不自然な振動を加えていない。
- **PIDのR依存曲線と、スライスの信号／背景応答を分離した。** Gaussian/lognormalの各R応答を候補多項式へ固定していない。背景核種は10/17のTOF応答参照に固定し、背景ピークを局所的な自由位置へ追加していない。

## Fit・積分・採否

native R幅は0.01 GeV/c、dE/dx幅は0.1 keV/cm。表示のみ0.5 keV/cmへ加算。
対象中心／幅はALLまたは対応TOF分布自身でfree、他核種応答は旧参照。Gaussianはmean/std、lognormalはmode/曲率幅の意味を保つ。
fit範囲は10/17と同じ対象局所窓、全応答参照とgateの±4幅の和集合を0.5 keV/cm境界へ合わせたもの。初期値、自由度、制約は各runのsettings.tsv・shapes.tsv・fixed_windows.jsonに保存。

| 条件 | 10/17と同じ設定 |
|---|---|
| 応答／背景 | Gaussianまたはlognormal × linearまたはexponential |
| 低R | 4モデル×nominal・nuisance中心±0.1σ・幅±10%=20条件 |
| 高R | 4 nominalモデル |
| TPC選択 | strict |nσ|<2、d/t/3He/4Heの最小|nσ|、同値はこの順 |
| 積分 | native RのGL3各点でwinner窓を計算し、各成分を積分 |
| 集約 | ΣS/Σ(S+B)。Purityの単純平均や曲線高さの比は不使用 |
| 対象帰属screen | TOF参照からの中心差≤1参照σ、幅比0.5–2、label crossingなし |
| 応答支持 | yield≥50、effective≥25、drop ΔD≥25、中心誤差/σ≤0.25、幅誤差/σ≤0.25等の旧条件 |
| 選択窓残差 | absolute≤15%、integral≤10%、未支持nuisance寄与とomission bound≤1% |

旧screen・数値収束・共分散・選択窓coverageを区別して保存。新しいχ²閾値、split-normal、GLS、複数モデルの追加条件や10/20のmaskは導入しなかった。
低Rは20条件すべて、高Rは全4条件が旧基準を通過した区間を評価可能とする。モデル帯へ不良fitは混入させない。
目視で見つかった低Rの別ピークへの移動は、全20条件を通過した支持点には混入していなかった。

誤差は旧方式の条件付きyield誤差の加重上界と、モデル／nuisance感度の範囲を別々に記録。帯は信頼区間でも全系統誤差でもない。
旧quality通過を、応答や物理的帰属の完全な検証と同一視しない。

## TOFと高Rの追加評価

R<2はALL、R≥2は対応TOFタグだけを使う。TAGはALLの部分集合なので、両者を独立標本として足し合わせない。
TOF不合格時に次点核種へ繰り上げない。

| 核種 | 保存されているstrict m²窓 |
|---|---|
| d / 4He | 3.268773 < m² < 3.693147 |
| t | 7.2564275 < m² < 8.2816925 |
| 3He | 1.78285 < m² < 2.06485 |

4Heはdと同じ窓内でdE/dxにより分離した。旧3<m²<4は使用していない。

10/17ではd/t/3Heの計算が2.5で止まっていたため、今回はその3核種も2.5–5の全native区間を12実行、3,000fit追加した。
保存TOF窓を再照合し、元Git応答は対象の初期値として使い、対象は分布中でfreeにした。
他核種のprovider-only応答は未検証として記録し、無視できない寄与があれば旧基準で未評価とした。
高Rの4Heは全2–5を0.05幅で独立fitする旧手順を維持。有限R投影幅を新しいintrinsic PID幅へ移植していない。

| 核種 | R=2からの連続支持域 [GeV/c] | 域内の最小モデル下側Purity | 全高R評価可能区間 |
|---|---|---:|---:|
| d | [2.00, 3.05) | 99.13% | 44 / 60 |
| t | [2.00, 3.40) | 98.60% | 42 / 60 |
| he3 | [2.00, 3.55) | 97.90% | 31 / 60 |
| he4 | [2.00, 3.20) | 98.39% | 24 / 60 |

d/tは連続域より上にも支持点がある。未評価の穴をつながず、最終図に各点を示す。
上表の端は統計限界の確定値ではなく、応答・帰属・残差・数値収束・統計の旧支持基準による連続域の端。
5 GeV/cも元He PIDの共通評価端であり、統計上限と呼ばない。

## 低Rと50%基準

| 核種 | 評価可能区間 | うちモデル下側≥50% | 評価可能点での最小モデル下側 |
|---|---:|---:|---:|
| d | 21 / 32 | 21 | 88.48% |
| t | 23 / 32 | 23 | 81.09% |
| he3 | 25 / 32 | 24 | 48.41% |
| he4 | 24 / 32 | 23 | 44.42% |

- 3He [0.75,0.80)のモデル下側は48.41%、4He [0.65,0.70)は44.42%。低Rを無理に広げる再調整はしなかった。
- 50%を超える評価可能点と未評価を区別した。3Heは[0.80,2.00)が連続して50%以上。
- d/t/4Heの支持域内には未評価の穴もある。未評価から自動的に局所vetoや単一cutoffを作らない。本番への変更はない。

## 10/17との比較

下表は両方で評価可能な低R区間だけの、nominal Gaussian+linear Purity差。微調整後と10/17の両方が20条件を通過する点で比較した。

| 核種 | 共通評価区間数 | 最大の|ΔP| |
|---|---:|---:|
| d | 21 | 0.0000 percentage points |
| t | 23 | 0.0000 percentage points |
| he3 | 25 | 3.5666 percentage points |
| he4 | 24 | 11.3081 percentage points |

d/t自身の式とfit応答は不変だが、競合するHeのpullが変わるとbest-speciesの窓が変わる。これをd/tの再校正と混同しない。
未調整Git元係数でも名目32実行3,400fitを別途実行し、[比較JSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/original_summary.json)へ保存。
4 nominalだけの比較を低R20条件と同じ認定結果にはしていない。
最終PDFは10/17の20条件／4条件と今回の同じ条件を重ねて比較する。

## 全R fit集

各ページで同じRの全4モデルを並べ、左linearは対象ピーク近傍、右logはfit全域を表示。
成分、fitの数値成立／旧支持判定を明記。不良fit・未評価点も省いていない。
COLZの±2σ線は単独帯のguideであり、実Purity積分にはbest-species境界を使った。

| 核種 | 低R | 高R |
|---|---|---|
| d | [160頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_low_001/fit_atlas_low_d.pdf) | [300頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_high_001/fit_atlas_high_d.pdf) |
| t | [160頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_low_001/fit_atlas_low_t.pdf) | [300頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_high_001/fit_atlas_high_t.pdf) |
| he3 | [160頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_low_001/fit_atlas_low_he3.pdf) | [300頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_high_001/fit_atlas_high_he3.pdf) |
| he4 | [160頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_low_001/fit_atlas_low_he4.pdf) | [60頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_high_001/fit_atlas_high_he4.pdf) |

## 検証と保護

- 低Rはfit独立の候補count 640行を元／候補の両方で測定し、モデルごとのcountと照合した。
- 低R32区間、高R60区間、核種別20/4モデルの網羅と非重複を検査した。
- 高R拡張の初回集約では区間原点の取り違えを検出し、最終図の前に修正した。修正前はsummary_extension_offset_diagnostic_001へ使用禁止と明記して保存。native fitとfit集のR座標は不変。最終集約はcoverageとduplicate-modelをfail-closedで検査する。
- 元cache・過去fit／source・本番設定など1,339ファイルのSHAを監査。元10GB ROOTはfingerprint不変を照合する。
- STAR SL24y / ROOT5.34/38のbatch一致環境。farm job投入、hadd、commit/push、メールはなし。本番PID・既存CF・Lambda/KFParticle・13p5は変更していない。
- Source／settings／command／各fitのROOT読戻しマーカーとcompletedを各runに保持。
- 再現手順は[tools README](../../../../tools/nuclear_pid_smooth_20261021/README.md)。

最終検証：低R4冊640頁、高R4冊960頁、Purity2頁、COLZ8頁、合計10冊1,610頁をGhostscriptで全頁読み戻し、頁数・SHAを検査してPASS。
候補16,640、Git元係数比較3,400、代表再現24、計20,064のfitスライスを処理した。元と共通のd/t評価可能区間のPurity差は最大3.33e-16。
元10GB ROOTのfingerprint、監査対象1,339ファイル、本番設定不変を最終照合した。最終ソース・設定は出力内final_sourceへ固定した。

[本日の日報](../../../../analysisnote/20261021/summary20261021.md)に実行環境・出力・未評価の意味を記録した。


## 追加依頼：1次元fit図の横軸を約2倍へ拡張（完了）

ユーザーが解析結果を了承した後、linear/logの横軸表示幅を約2倍へ広げる依頼を受けた。
fitの再実行、PID係数、Purity、採否、元ヒストグラムは変更しない。
旧atlas_low_001/atlas_high_001を保持し、新しいatlas_low_002/atlas_high_002へ全4核種を出力する。
各パネルは元の表示範囲の中点を保って幅2倍とし、保存ヒストの0–120 keV/cm軸でclipする。
旧図で表示範囲が逆転していたものは、数値を変えずヒスト全軸表示とする。
fit範囲の外側は生ヒストのみを表示し、未計算のfit曲線をゼロとして描いたり外挿したりしない。
4Heの代表3区間をプレビュー・目視・全頁読戻し済み。全1,600頁の新版を再生成し、全頁読戻し・頁数・SHA検証がPASS。fit結果・Purity・元001版PDFの不変性も確認した。

[表示幅2倍の新版・全8冊の索引](fit_atlas_wide_20261021.md)。
設定と手順：tools/nuclear_pid_smooth_20261021/wide_atlas_settings.yaml、render_wide_atlas.py。

## 追加依頼：Purity図1ページ目の横軸を核種別に調整（完了）

[更新版PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/purity_axis_002/purity_vs_p_over_q.pdf)／[1ページ目PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/purity_axis_002/page1.png)。
元PDFを保持し、purity_axis_002へ表示更新版を保存した。dは元の表示のまま。
t/3He/4Heは実際に描画する評価可能点（10/17比較点を含む）のbin端から0.05 GeV/cずつ余白を付けた。
未評価マーカーだけが続く領域で横軸を広げない。

| 核種 | 1ページ目の横軸 [GeV/c] |
|---|---|
| d | 0.40–5.00（元のまま） |
| t | 0.75–4.90 |
| 3He | 0.70–3.60 |
| 4He | 0.60–3.25 |

再fit・Purity再計算なし。入力表・数値成果物・元PDFのSHA不変を確認。
全2ページを読み戻し、2ページ目は同条件で画像化した元PDFと画素単位で一致した。
設定・ソース・command・検証manifestをpurity_axis_002へ保存。
再現手順はtools/nuclear_pid_smooth_20261021/render_purity_axes.py。

## 追加比較：背景なし・3He半値までの実fit範囲拡張

[新しい比較図と報告](no_smooth_he3window_20261021.md)。ユーザー追加依頼で低R 4Heの実fit下限を拡張し、背景あり・なしを同じ範囲で再fitした。表示だけの変更とは別のstudy。元の確定成果物は保持し、最終Purityへの採用変更は行っていない。

## 最新追加：4He R>0.61・低Rピーク位置

[更新fit図・Purity・実施報告](he4_r061_20261021.md)。ユーザー指定の³He/⁴Heピークを低R参照へ反映し、4He適用下限をR>0.61にした独立study。0.70以上は直前の半値範囲fitと数値差ゼロ。旧成果物は保持。最新の4He統合Purityはhe4_r061_001内に保存している。

## 最終採用モデルの更新：Gaussian＋linear

ユーザーの目視判断・確認により全核種Gaussian＋linearの名目fitだけを採用。[最新Purity図と図入りPDF報告](gaussian_linear_report_20261021.md)。本ページの旧モデル範囲は比較履歴として保持し、最新の採用図では他モデルの点・範囲を混ぜていない。
