# 高R・TOF要求後の全1D fit図

2026-10-17：全24冊・840頁、100,043,218 bytesを作成し、閉PDF読戻し・頁数・SHA・保存fit入力の不変性を確認した。[completed.json](completed.json)参照。新fit・正規化・PID変更は行わない。

| 核種 | sample | R [GeV/c] | R幅 | 冊数 | 頁数 |
|---|---|---|---|---:|---:|
| d | TAG_D | 2.0–2.5 | 0.01 native | 3 | 200 |
| t | TAG_T | 2.0–2.5 | 0.01 native | 3 | 200 |
| ³He | TAG_HE3 | 2.0–2.5 | 0.01 native | 3 | 200 |
| ⁴He | TAG_D | 2.0–5.0 | 0.05 coarse直接fit | 15 | 240 |

| R [GeV/c] | d | t | ³He |
|---|---|---|---|
| 2.0–2.2 | [80頁](d_R200_220.pdf) | [80頁](t_R200_220.pdf) | [80頁](he3_R200_220.pdf) |
| 2.2–2.4 | [80頁](d_R220_240.pdf) | [80頁](t_R220_240.pdf) | [80頁](he3_R220_240.pdf) |
| 2.4–2.5 | [40頁](d_R240_250.pdf) | [40頁](t_R240_250.pdf) | [40頁](he3_R240_250.pdf) |

⁴Heは各16頁：[2.0–2.2](he4_R200_220.pdf)、[2.2–2.4](he4_R220_240.pdf)、[2.4–2.6](he4_R240_260.pdf)、[2.6–2.8](he4_R260_280.pdf)、[2.8–3.0](he4_R280_300.pdf)、[3.0–3.2](he4_R300_320.pdf)、[3.2–3.4](he4_R320_340.pdf)、[3.4–3.6](he4_R340_360.pdf)、[3.6–3.8](he4_R360_380.pdf)、[3.8–4.0](he4_R380_400.pdf)、[4.0–4.2](he4_R400_420.pdf)、[4.2–4.4](he4_R420_440.pdf)、[4.4–4.6](he4_R440_460.pdf)、[4.6–4.8](he4_R460_480.pdf)、[4.8–5.0](he4_R480_500.pdf)。表示範囲すべてが採用範囲ではない。

1区間4頁でGaussian／log-normal×linear／exponential背景の全4モデルを、全dE/dx（0–120 keV/cm）とfit範囲のlinear／logで表示する。dE/dx幅はnative 0.1 keV/cm。凡例はTAG observed countsで、低RのALLではない。

各パネルに核種・R区間・sample・モデル・R幅・TOF窓を表示する。fitは数値判定、responseは応答品質、saved runQはモデル別保存品質で、最終4モデル包絡や完全な系統誤差ではない。失敗・未支持を省略せず、purity=0へ置換しない。[coverage.tsv](coverage.tsv)は全840 model/区間の状態・理由を記録する。

raw TOF (m/q)²の厳密開区間[GeV²/c⁴]はD=(3.268773,3.693147)、T=(7.2564275,8.2816925)、³He=(1.78285,2.06485)。⁴Heの[study]はfit保存時点の旧本番(3,4)窓との差を表す来歴で、後日の本番採用状態とは別。D窓はd／⁴He互換で純⁴He truthではない。

⁴Heは元5 native binの厳密加算後に別stageで直接fitし、native purityを平均・補間していない。射影fit幅をintrinsic detector幅とは断定しない。[低R全fit集](../../fit_atlas_native_20261017_003/README.md)や元大ROOT・既存fitは変更していない。

[plan.json](plan.json)に計画と元入力SHA、source/に表示source・YAML・runtime、inputs/に元fit context16個と応答品質表を保存した。各PDFの検査は*.completed.json、全体はcompleted.json。ROOT5コンパイル・空入力拒否は[検証ログ](../atlas_renderer_validation_001/compile_guard.log)、目視確認は[He4例](he4_highR_fit_page3.png)。元大ROOTを再読していない。

再現入口はtools/nuclear_pid_cutoff/highrigidity/render_highrigidity_atlas.py。新規出力へprepare --settings tools/nuclear_pid_cutoff/highrigidity/atlas_settings.yaml --output-dir NEW_DIR、初冊run --output-dir NEW_DIR --book he4_R480_500 --workers 1、目視後run --output-dir NEW_DIR --book all --workers 2で全冊と最終checkを行う。既存PDFは上書きしない。
