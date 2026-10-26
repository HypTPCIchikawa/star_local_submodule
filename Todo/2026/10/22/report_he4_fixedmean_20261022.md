# 4He低R：Gaussianのみ・中心固定の再評価（002版）

**作成・検証完了。** ユーザー訂正に従い、0.61<R<0.80の19 native区間だけ、前版Gaussian＋linearの⁴He中心に固定してlinear背景を除去した。⁴He幅と核種収量は自由。R≥0.8の180 native区間・84集約点、全PID係数、TPC選択・TOF要求は前版と一致。

## 成果物

- [改訂版の図入りPDFレポート（11頁）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_he4_20261022_002/he4_pid_purity_report_20261022.pdf)。6頁が低R中心固定・Gaussianのみの3代表点、11頁が変更前後のPurity・幅と解釈。
- [低R中心固定・Gaussianのみ、全19区間fit PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_he4_20261022_002/fit_low_fixed.pdf)
- [⁴He全199区間fit集](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_he4_20261022_002/fit_all_regions.pdf)。先頭19頁が今回の再fit、以後は前版と同じデータ・モデル・fit値。
- [統合Purity図（2頁）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_he4_20261022_002/purity_vs_p_over_q.pdf)
- [全集約値と支持条件](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_he4_20261022_002/summary.json)、[固定中心の照合と幅](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_he4_20261022_002/fixed_center_validation.json)、[低R変更前後](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_he4_20261022_002/low_comparison.json)、[設定](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_he4_20261022_002/settings.yaml)、[最終検証manifest](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_he4_20261022_002/completed.json)

## 現在のモデル

| R [GeV/c] | 試料 | ⁴He応答中心 | fit成分 |
|---|---|---|---|
| 0.61–0.80 | ALL、TOFなし | 前版Gaussian＋linear値に固定 | 4核種Gaussianのみ、smoothなし |
| 0.80–2.00 | ALL、TOFなし | 自由、前版結果を保持 | 4核種Gaussianのみ、smoothなし |
| 2.00–5.00 | 既存TAG_D | 自由、前版結果を保持 | ⁠⁴He Gaussian＋linearのみ |

⁴Heはstrict R>0.61。0.80は中心自由側に含める。d/t/³HeはGit 8f7d539のPID関数形・係数・適用域のまま。⁴Heの候補PID係数も001版から再調整しない。変更したのは低Rのスペクトルfitモデルと⁴He応答中心の自由度であり、TPCのPID関数の中心を固定したという意味ではない。

今回の最終条件は**Gaussianのみ**。最初の発言を受けたGaussian＋linear維持の案は、ユーザーの直後の訂正で取り下げ、実行していない。⁴He中心だけは前版Gaussian＋linearの値を使い、各GL3点の応答変化は前版のnative内形状を維持する。低R³Heの較正中心・0.70未満の幅自由などの既存条件も保持する。

## 再評価の結果

| 集約R | 前版Gaussian＋linear | 中心固定Gaussianのみ | 新モデルの支持 |
|---|---:|---:|---|
| 0.61–0.65 | 63.6581% | 99.6185% | 不通過 |
| 0.65–0.70 | 87.0207% | 99.8125% | 通過 |
| 0.70–0.75 | 75.8713% | 99.9516% | 通過 |
| 0.75–0.80 | 87.4978% | 99.9413% | 通過 |

固定した中心は19区間すべて元の値と全桁一致。例えばR=0.61–0.62ではμ=67.06310514545 keV/cm、幅は前版3.63034から6.58988 keV/cmへ広がる。0.64–0.65の中心は64.38844196478 keV/cmで、幅は3.78698→6.37176。

最初の2 native区間では幅比が従来のidentity参照に対して2.08694、2.12677となり、既存上限2を超える。従って最初の0.61–0.65集約点は支持しない。図では下端の×で位置を示し、計算値とfit図は表・全区間PDFに残した。採否閾値は緩めていない。

支持点は低R中心固定3/4、中R24/24、高R48/60、合計75/88。R≥0.8の採否・Purity・誤差も前版と完全一致。背景を除いたためPurityが1に近づくが、データの背景不存在やfit品質の改善を証明したものではない。中心を固定しても幅や残差の問題が残る。

固定中心のmeanError=0は条件付きであり、前版で同じデータから推定した中心の較正不確かさを含まない。自由パラメータから中心を除き、幅は自由のまま、NDFが残る形状・収量自由度と一致することを検証した。縦誤差もモデル・形状に条件を置いた収量誤差で、完全な系統誤差ではない。

## 検証と再現

- 新しいfitは19区間のみ。他の180区間のJSON/ROOT等は前版の保存値とSHA一致、84集約点は辞書の全値が一致。全PID係数・候補ヘッダ・PID曲線TSVもSHA一致。
- 固定中心19件の全桁一致、meanError=0、幅自由、shapeパラメータ数1減とNDFを検証。
- 独立生計数556件を照合。全低R139区間の³He低側半値を417 GL3点で包含。
- 全199区間のraw内容・誤差不変、成分和、R<2でsmooth=0、高R d/t/³He=0をROOTで検証。
- 前版と過去入力236ファイルのSHA不変。各PDF405頁＋report11頁、計416頁を全頁読み戻し、頁数・エラーなしを確認。レポートは欠字・はみ出しゼロ。全頁配置、代表fit、変更前後の表を目視確認した。
- 実装：tools/nuclear_pid_he4_fixedmean_20261022/。run.py prepare → run.py fit → analyze.py → validate.py → render.py → build_report.py → finish.py。設定、ソース差分、実行コマンドと検証を出力先へ保存。
- mainconf=config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml。STAR SL24y/ROOT5、script/singularity_root_batch.sh。既存2021 P24iy部分統合cacheを使用し、元ROOT fingerprintとcache SHAを照合。新規jobid/configlogなし（ローカルstudy）。本番反映・job投入・hadd・commit/pushなし。

[前版001の記録](report_he4_20261022.md)とPDFは履歴として保持。今回の改訂結果は002版を参照する。

## ユーザー確定後の適用（2026-10-22追記）

Purity確定の指示を受け、⁴HeのPID係数を対象解析へ反映済み。d/t/³HeのGit PIDは保持し、全核種でR<2はTOF不要、R≥2はTOF必須。[適用設定・検証記録](pid_adoption_20261022.md)。本資料の「本番未反映」はstudy完成時点の履歴として保持する。

## TOF選択後R–dE/dx図の追加（作業中）

既存TAG_D分布と今回の確定PID係数を使い、TOF m²選択後のマップをreportへ追加する。旧PIDを重ねた過去図の関数だけを今回の係数へ揃え、fit/Purity/入力データは変更しない。旧11頁reportはtests/report_tof_map_20261022_001/に保持。

### TOFマップ追加：完了

指定のreportを同じパスで12頁へ更新し、6頁へTOF m²選択後のR–dE/dx＋PID曲線を追加。既存TAG_D試料（strict 3.268773<(m/q)²<3.693147、TPC PID前）を使い、R=0.61–5の全体と2–5の拡大を掲載。過去図のPID係数は古いため、今回のcoefficients.json/curves.tsvに揃えて重ね直した。d/t/³HeはGit、⁴Heは元と確定係数、中心と±2σを表示。低R PurityはTOF不要ALLのまま。

report全12頁と追加図を読戻し、追加頁と全頁配置を目視確認。欠字・はみ出しなし。保護対象のPurity/係数/全fit PDF等のSHA不変を確認。旧11頁reportと更新前tex/manifestは解析rootfile配下tests/report_tof_map_20261022_001/へ保持し、現manifestを更新。解析本体・稼働job・通知設定・Git pushには変更なし。
