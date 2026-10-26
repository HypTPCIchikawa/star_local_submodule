# 横軸表示範囲を約2倍にした1次元fit集

完了：全8冊・1,600ページを生成し、全頁読戻し・頁数・SHA検証がPASS。fit数値・Purity・元PDFの不変性も確認した。

各R・各モデルの左linearと右logの横軸を、元の表示中点を保って幅2倍へ拡張した。
保存ヒストの0–120 keV/cmを超える部分は軸端でclipする。
旧図で表示範囲が成立していなかった低Rの5パネルはヒスト全軸表示とした。

fit範囲・関数・パラメータ・Purity・採否・生ヒストグラムは変更していない。
広げた領域では生ヒストを表示し、fit曲線は保存済みfit範囲内にだけ描く。
001版のPDFと元の数値成果物は保持している。

| 核種 | R<2 GeV/c | R≥2 GeV/c |
|---|---|---|
| d | [低R：160頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_low_002/fit_atlas_low_d.pdf) | [高R：300頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_high_002/fit_atlas_high_d.pdf) |
| t | [低R：160頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_low_002/fit_atlas_low_t.pdf) | [高R：300頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_high_002/fit_atlas_high_t.pdf) |
| he3 | [低R：160頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_low_002/fit_atlas_low_he3.pdf) | [高R：300頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_high_002/fit_atlas_high_he3.pdf) |
| he4 | [低R：160頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_low_002/fit_atlas_low_he4.pdf) | [高R：60頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_smooth_20261021/atlas_high_002/fit_atlas_high_he4.pdf) |

低R4冊640頁、高R4冊960頁、合計8冊1,600頁。
新版出力はatlas_low_002とatlas_high_002。
各ディレクトリのcompleted.jsonへ、全頁読戻し・頁数・SHA・数値結果と旧PDFの不変性を記録する。

設定・再現手順：tools/nuclear_pid_smooth_20261021/wide_atlas_settings.yaml、render_wide_atlas.py、PlotWideAtlas.C。
元の[解析結果・条件](report_smooth_pid_20261021.md)はそのまま使用する。

## 追加比較：背景なし・3He半値までの実fit範囲拡張

[新しい比較図と報告](no_smooth_he3window_20261021.md)。ユーザー追加依頼で低R 4Heの実fit下限を拡張し、背景あり・なしを同じ範囲で再fitした。表示だけの変更とは別のstudy。元の確定成果物は保持し、最終Purityへの採用変更は行っていない。
