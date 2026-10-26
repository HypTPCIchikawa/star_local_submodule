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
