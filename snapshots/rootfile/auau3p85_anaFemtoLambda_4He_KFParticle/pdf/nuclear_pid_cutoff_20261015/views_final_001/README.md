# 最終集計の図（2026-10-17 JST再開）

10月15日に完了した`summary_free_sensitivity_001`と`adopted_ranges_001`から作成した表示用PDFです。fit・イベント処理・purity再計算は実施していません。旧図は保持しています。

| PDF | ページ数 | 内容 |
|---|---:|---|
| [purity_cutoff_overview.pdf](purity_cutoff_overview.pdf) | 3 | purityのモデル包絡と使用範囲／除外範囲、ALLの広いdE/dx COLZと現行PID窓、実候補数と残す範囲。 |
| [response_parameters.pdf](response_parameters.pdf) | 5 | nominal 4モデルのALL free fit中心・幅と現行関数。最後はALL COLZへの重ね描き。 |
| [fit_spectra_d.pdf](fit_spectra_d.pdf) | 12 | dのR=0.605、0.655、1.795 GeV/c。 |
| [fit_spectra_t.pdf](fit_spectra_t.pdf) | 12 | tのR=0.755、0.805、1.395 GeV/c。 |
| [fit_spectra_he3.pdf](fit_spectra_he3.pdf) | 12 | ³HeのR=0.705、0.755、1.995 GeV/c。 |
| [fit_spectra_he4.pdf](fit_spectra_he4.pdf) | 16 | ⁴HeのR=0.655、0.705、0.855、1.995 GeV/c。 |

代表スペクトルは各Rにつき、全範囲0〜120 keV/cmと保存されたfit窓の2表示、Gaussian／log-normal × linear／exponential背景の4モデル、linear／logの両軸を表示します。native dE/dx bin幅0.1 keV/cmを保持し、表示の規格化・再fitは行いません。未検証点も診断として表示するため、曲線が描かれていることをfit支持の証拠にしないでください。

overviewの斜線は使用しないR範囲、緑線は今回書き出した使用候補範囲です。除外にはpurity不足だけでなく、数値共分散・モデル整合性が未確認の区間も含みます。淡青帯は比較したモデル・感度の包絡であり、信頼区間や全系統誤差ではありません。応答比較はnominalかつown-fit応答検査を通った点であり、すべてのpurity採用判定を通過した点だけを意味しません。

`source/`に描画ソース、`inputs/`に集計表のコピー、`*.command.txt`と`*.log`に実行記録を保存しました。`run_views.sh`はbatch一致ROOT5の描画入口で、既存PDFを上書きしません。`check_views.py`がGhostscript全ページ読み戻し、ページ数、入力・保存ROOT・PDF・ログのSHAを検査し、[completed.json](completed.json)へ保存します。目視QA PNGの`parameters_page_1/2.png`は応答PDFの4/5ページに対応します（Ghostscriptの出力番号は1から始まります）。

PIDの中心・幅の係数、元ROOT、既存mainconfはこの描画で変更していません。使用候補範囲を解析へ適用するYAMLの状態・実イベント検証は[Todo報告](../../../../../Todo/2026/10/15/report_nuclear_pid_cutoff_20261015.md)を参照してください。
