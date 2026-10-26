# Gaussian＋linearのみのPurityと図入りPDFレポート

**最新版：** [各核種3代表点を含む改訂版PDF（15頁）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_002/pid_purity_report_20261021.pdf)。TOFなしfitは8–11頁。

**ユーザー確認済み：全核種Gaussian応答＋linear背景の名目fitだけを採用。作成・検証完了。**

- [図入り日本語PDFレポート（13頁）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_001/pid_purity_report_20261021.pdf)
- [更新Purity図（全体・TOF拡大の2頁）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_001/purity_vs_p_over_q.pdf)
- [選択した点・誤差・採否のTSV](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_001/purity.tsv)、[元値と統計・理由のJSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_001/summary.json)
- [PID式・係数](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_001/coefficients.json)、[study用PID選択設定](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_001/pid_selection.yaml)、[再現設定](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_001/settings.yaml)
- [最終検証manifest](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_001/completed.json)、[364区間の元値照合](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/gaussian_linear_report_001/point_validation.json)

## 採用条件

全d/t/³He/⁴HeでGaussian＋linear、nominalのみ。BGなし、lognormal、exponential、mean/width感度変化の点や包絡線はPurity図へ含めない。誤差はGaussian＋linearの条件付き収量誤差で、native区間の集約には分母重みの線形和による保守的上限を使う。応答・背景選択を含む完全な系統誤差ではない。

採否もGaussian＋linear自身のsupportedを使い、以前の全モデル共通maskを残さない。支持されたPurity<50%の点も除かず表示。曲線で未評価点の穴をつながない。横軸は支持点の範囲に合わせ、dの全体図は従来の0.4–5を維持。

4Heはhe4_r061_001の最新低R fit、他核種と高Rは既存のモデル別結果を使用。新しいfit、イベント解析、本番設定変更は行っていない。Gaussian＋linearの初期値比較で選ばれた保存解を使い、他モデルのPurityを混合したものではない。

## 結果

| 核種 | 領域 | 支持区間/全区間 | 支持点R最小–最大 | 最小Purity |
|---|---|---:|---:|---:|
| d | low | 24/32 | 0.40–1.75 | 88.50% |
| d | high | 55/60 | 2.00–4.75 | 99.13% |
| t | low | 24/32 | 0.80–2.00 | 85.08% |
| t | high | 57/60 | 2.00–4.85 | 98.60% |
| he3 | low | 27/32 | 0.65–2.00 | 29.33% |
| he3 | high | 31/60 | 2.00–3.55 | 97.92% |
| he4 | low | 28/28 | 0.61–2.00 | 45.04% |
| he4 | high | 27/60 | 2.00–3.35 | 98.08% |

⁴Heの0.61–0.65は63.70 ± 0.07%（条件付き集約誤差上限）。最初のnative [0.61,0.62)は59.89%。支持された集約点で50%未満なのは³He [0.65,0.70)=29.33%、[0.70,0.75)=34.20%、⁴He [1.95,2.00)=45.04%。新しいvetoは追加していない。

## PDFの内容

日本語本文に加え、最新Purity全体/TOF拡大、PID中心/幅のGit元との比較、ALL dE/dxマップ、⁴He低Rの2点較正、全4核種の低R/TOF代表fitを収録。全10頁分の図を含む13頁。PID中心3次・幅2次の式と係数、strict 4He R>0.61、TOF m²窓、ΣS/Σtotal、誤差・採否・制限と再現方法を説明する。

PID選択窓の多項式と各Gaussian fitの中心・幅を区別し、低R³He参照中心の較正がTPC PID多項式の局所補間ではないことを明記。本番未反映。³He中心を独立再測定したと扱わず、指定較正値を参照に固定したことを記載。

図の元vector PDFを保持し、レポートには200 dpiのPNGを使用。元図の白い余白だけをtrimして大きく掲載し、データや軸ラベルは変更していない。日本語は既存フォントで組版し、欠字・本文はみ出しゼロを確認。

## 再現と検証

- 364区間のPurity、S、総数、生計数、条件付き誤差、採否をモデル別元JSONへ完全一致照合。支持点273、未評価点91。掲載代表fitのresponse/backgroundメタデータもGaussian/linearを確認。
- 旧入力23ファイルのSHA不変。新しいfitなし。図10頁＋PDFレポート13頁を全頁読み戻し・頁数・エラーなし・SHA検証。代表頁と全頁配置の目視確認も実施。
- 実装：tools/nuclear_pid_smooth_20261021/linear_report/{settings.yaml,prepare.py,PlotLinearReport.C,render.py,build_report.py,finalize.py}。prepare→render→build_report→finalizeの順。source/へ最終版を保存。
- mainconf=config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml。描画はSTAR SL24y/ROOT5、script/singularity_root_batch.sh経由。PDF組版はXeLaTeX。jobid/configlogなし（ローカルの抽出・描画・組版）。job投入・hadd・commit/pushなし。
- 旧モデル比較図・旧報告は履歴として保持。今回の採用モデル図はgaussian_linear_report_001へ独立保存。


## 代表点3点への追加改訂（002版）

各核種d/t/³He/⁴HeのTOF m²選択なしALL fitをR=[0.80,0.81)、[1.20,1.21)、[1.60,1.61) GeV/cの3点で掲載する。全12点はGaussian＋linearの支持条件を通過。旧001版を保持し002版を作成済み。再fitとPurity変更は行わない。

改訂版の8–11頁はd、t、³He、⁴Heの順に各3点×linear/logの6パネル。低Rピーク較正2点と高R TOF代表図も保持。364区間のPurity・誤差・支持条件・表示軸が001版と完全一致、purity.tsvもSHA一致。全12代表点の元メタデータでALL/Gaussian/linearを確認。

再現：tools/nuclear_pid_smooth_20261021/linear_report_three_points/ のprepare.py → render.py → build_report.py → finalize.py。設定はsettings.yaml、検証証拠は002版completed.jsonとpoint_validation.json。

検証完了：レポート15頁＋図12頁を全頁読戻し、欠字・はみ出しゼロ。旧入力と001版の28ファイルSHA不変。代表頁と全15頁の配置を目視確認済み。


## 後続のGaussianのみ試行

ユーザーの追加指示でlinearを全て除いた全核種・全Rの再評価へ進んだ。[Gaussianのみの最新結果・全区間PDF](he4_r061_20261021.md)。本ページのlinear版15頁レポートは前回結果として保持。
