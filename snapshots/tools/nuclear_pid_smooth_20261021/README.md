# 2026-10-21 Git元の滑らかな核PIDとPurity

指定Todo `Todo/2026/10/21/STAR_Todo20261021.txt` に対応する独立study。
本番Maker・PID・cutoff・CFは変更しない。

入力は2021 P24iy部分統合由来の保存native cache。数学は10/17の保存C++/headerを複製して実行する。
現在の解析コードをfitへ読み込まない。mainconfはSL24y/ROOT5環境の選択に使う。
Git `8f7d539` のPIDヘッダを直接取り出し、当時baselineのSHAと一致させた。

`settings.yaml`はJSON互換YAMLで、参照先・出力先・微調整制限・TOF窓を明示する。
別の再実行はoutputBaseを新しい未使用先へ変更する。元出力や既存runへ異なる設定を混ぜない。
同じ設定の完了runだけstudy.pyが再利用でき、未完了runは上書きしない。

## 再現順序

リポジトリ直下で実行する。最初の回帰だけはホストの混在PYTHONPATHを外してシステムNumPyを使う。
それ以外は標準Pythonと既存ROOT5 wrapper。解析条件を環境変数から供給しない。

```bash
env -u PYTHONPATH /usr/bin/python3 tools/nuclear_pid_smooth_20261021/calibrate.py
python3 -B tools/nuclear_pid_smooth_20261021/study.py reproduce --workers 2
python3 -B tools/nuclear_pid_smooth_20261021/audit.py
python3 -B tools/nuclear_pid_smooth_20261021/study.py candidate --workers 2
python3 -B tools/nuclear_pid_smooth_20261021/count_gates.py candidate
python3 -B tools/nuclear_pid_smooth_20261021/study.py original --workers 2
python3 -B tools/nuclear_pid_smooth_20261021/count_gates.py original
python3 -B tools/nuclear_pid_smooth_20261021/extend_high.py
python3 -B tools/nuclear_pid_smooth_20261021/summarize.py original
python3 -B tools/nuclear_pid_smooth_20261021/summarize.py candidate
python3 -B tools/nuclear_pid_smooth_20261021/prepare_views.py
python3 -B tools/nuclear_pid_smooth_20261021/prepare_atlas.py low
python3 -B tools/nuclear_pid_smooth_20261021/prepare_atlas.py high
```

`PlotStudy.C+(OUTPUT,CACHE,0)`で統合図とCOLZを出力。
`PlotStudy.C+(OUTPUT/atlas_low_001,"",1)`とhigh版でfit集を出力する。
必ず`script/singularity_root_batch.sh MAINCONF root4star -l -b -q ...`経由で呼ぶ。
実際の絶対パス・引数は各runのcommand.jsonと日報に記録。

```bash
python3 -B tools/nuclear_pid_smooth_20261021/validate_pdfs.py low
python3 -B tools/nuclear_pid_smooth_20261021/validate_pdfs.py high
python3 -B tools/nuclear_pid_smooth_20261021/validate_pdfs.py summary
python3 -B tools/nuclear_pid_smooth_20261021/finalize.py
```

## 解釈

- PID中心はlog10(R)の3次、幅は2次。d/tは固定、Heだけ微調整。
- 各スライスのGaussian/lognormal応答は旧target-free fitのまま。PID gateの多項式へ固定しない。
- 低RはALL、R≥2は対応TOFタグ。4Heはd窓。ALLとTAGを加算しない。
- 低R20条件・高R4条件すべてを旧基準で検査。新しいchi²閾値・GLS・split-normal・採否maskは導入しない。
- 積分したSと総数を加算後に割る。モデル帯と条件付きyield誤差の上界は別。
- 2.5超のd/t/3Heは旧high-R engineで全域を追加fit。対象はデータ内free、他核種のprovider-only応答は未検証とし、既存のmaterial-nuisance条件を適用。
- 未評価をP=0や50%未満としない。高Rの孤立した支持点から穴をつないだ上限は作らない。
- assessed_ranges.jsonはstudyの支持区間であり、本番cutoff/vetoではない。

[当日報告](../../Todo/2026/10/21/report_smooth_pid_20261021.md)を正式索引とする。

## 表示範囲を広げた1次元fit集（追加依頼）

`render_wide_atlas.py`は001版の保存atlas manifestとfit ROOTだけを読み、
linear/logの表示幅を`wide_atlas_settings.yaml`のwidthScale=2へ広げる。
再fitやPurity再計算はしない。0–120 keV/cmでclipし、旧表示範囲が不成立だったものは全軸表示。
fit曲線は保存済みfit範囲内だけに描き、外側のdataへ外挿しない。

```bash
python3 -B tools/nuclear_pid_smooth_20261021/render_wide_atlas.py low
python3 -B tools/nuclear_pid_smooth_20261021/render_wide_atlas.py high
```

出力はatlas_low_002とatlas_high_002。各runのcompleted.jsonへ全PDFの頁数・SHA、
数値成果物・元PDFの不変性確認を記録する。既存出力には上書きしない。
別の再描画はwide_atlas_settings.yamlのoutputSuffixを未使用値へ変える。

## Purity図1ページ目の核種別横軸

`python3 -B tools/nuclear_pid_smooth_20261021/render_purity_axes.py`は保存表から
実際の支持点のbin端を集計し、purity_axis_settings.yamlの余白を付ける。
dは元の範囲を維持する。出力purity_axis_002は既存なら上書きせず停止。
2ページ目は元PDFとの画素一致を検査し、全2頁の読戻しと元の数値・PDFの不変性を記録する。

## 4He R>0.61・低Rピーク修正

`he4_r061/settings.yaml`に適用下限、ユーザー指定3He/4Heピーク、低R補正域、半値範囲規則を保存。
3Heの中心は指定ピーク参照、幅は低RだけALL fit、4He中心/幅は自由fit。TPC PIDの元の中心3次/幅2次式は保持し、4He適用下限を厳密なR>0.61とする。

新しい空の出力先で`he4_r061/run.py prepare`→`full`→`warm`、`analyze.py`、`ValidateR061.C`、`render.py`、`finish.py`の順。Pythonは各ディレクトリのスクリプトをproject rootから実行。初回ACLiCビルドは直列化する。
設定した2点は初期値・参照形状に反映し、4He fitを2点へ強制固定しない。低R補正は0.70以上に影響させず、前回の同じ半値範囲fitとの数値一致を確認する。
詳細な成果物・制限は`Todo/2026/10/21/he4_r061_20261021.md`。最終数値は旧結果を置き換えず`he4_r061_001`に保存する。

## Gaussian＋linear採用と図入りreport

`linear_report/settings.yaml`で全核種gaussian/linear/nominalを選択。prepare.py→render.py→build_report.py→finalize.py。既存モデル別JSONの点・誤差・支持条件を直接抽出し、他モデルの共通mask・包絡線を使わない。ROOT vector図、200dpi画像、白余白を除いたレポート用画像、日本語13頁PDF、元値照合をgaussian_linear_report_001へ保存。新fit/本番変更なし。
