# TOFタグを手掛かりにした原子核PID study

対象は3.85 GeV・2021 P24iyの保存済み高統計ROOT。
ユーザー指定の[今回Todo](../../Todo/2026/10/14/STAR_Todo_20261014.txt)・
[計画](../../Todo/2026/10/14/plan_tof_assisted_nuclear_pid.md)・
[最新結果](../../Todo/2026/10/14/report_tof_assisted_nuclear_pid.md)を参照。
2026-10-19の整理で、Todo記録を`Todo/2026/10/14/`へ移動した。過去の実行snapshotは保持する。

## 元データを保護する方針

- `true/hDedxP`とTOFのd／t／³He窓のヒストを、native binと保存countのまま使う。
- TOFdはd／⁴He互換のタグであり、純粋な⁴He試料ではない。
- ALLとTAGは部分集合として相関する。独立尤度として重複加算しない。
- 古いROOT・PDF・旧fit・共通PID・現行YAMLを上書きしない。新しい出力先だけを使う。
- 実行成功、fit収束、核種の帰属、TAGからALLへの転用、校正の採用は別の判定。
  ここでの候補JSONは`productionEnabled: false`で、本番設定ではない。

## 処理の入口

| ファイル | 役割 |
|---|---|
| `prepare_inputs.py` / `ExtractTofInputs.C` | native入力cache、ALL−TAGのremaining、監査と全軸COLZ・射影図。 |
| `input_settings.yaml` | 入力図・保護対象・軸の明示設定。 |
| `run_tag_study.py` / `FitTofTagShapes.C` | 非重複R分割のTAG fit、native Gaussian bin積分、linear／exp背景。 |
| `tag_study_settings.yaml` | R分割、初期値、数学的設定、支持域判定、関数比較条件。 |
| `runtime_selection.py` | 実際のC++係数と現行YAMLを読み、現在のTPC窓を再現（変更しない）。 |
| `build_curves.py` / `nuclear_curves.py` | 背景依存・不良fitを分離し、滑らかな中心／幅の候補を作る。 |
| `run_all_study.py` / `FitAllNuclearPurity.C` | 固定TAG候補のALL fit、実際のbest-species窓での条件付きpurity。 |
| `render_study.py` / `PlotNuclearStudy.C` | 保存fitの全軸linear／log、fit窓の補助図。 |
| `run_lognormal_study.py` / `FitLognormalTofTagShapes.C` | 原Gaussianと分離したlog-normal応答の比較。mode／局所曲率幅の定義で、本番への採用は禁止。 |
| `lognormal_study_settings.yaml` / `lognormal_full_native_settings.yaml` | 同入力・同窓のR=0.80〜0.90試験／R=0.30〜2.50全域の明示設定。 |
| `audit_tag_shapes.py` / `AuditTagShapes.C` | ピーク±2／±3幅と全fit窓の局所残差、0.1／0.5／1 keV/cm集約の診断。再fitなし。 |
| `compare_responses.py` / `CompareNuclearResponses.C` | 同観測データのGaussian／log-normal重ね、中心・幅・適合指標のPDFと数値表。 |
| `render_guides.py` / `PlotNuclearGuides.C` | 生COLZへ名目fit点、失敗点、実際の現PID・legacy参照、支持域内だけの関数を重ねる。 |
| `display_settings_readable.yaml` / `overlay_settings_readable.yaml` | ROOT5の紙面内に軸・注記を収める表示専用設定。旧表示設定・出力は保持。 |
| `SelfTest*.C` / `test_*.py` | ROOT5の数値処理・実選択再現・曲線連続性・未支持域の回帰試験。 |

## 入力抽出の再現例

リポジトリ直下で実行する。新しい未使用の出力名を指定する。

```bash
python3 -B tools/nuclear_pid_tof/prepare_inputs.py \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --settings tools/nuclear_pid_tof/input_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_NEW
```

各runnerは`--help`で入力を確認し、明示したYAMLと新しい出力先を渡す。
SL24y／ROOT5は`script/singularity_root_batch.sh MAINCONF`を通して使う。
Python補助の回帰試験:

```bash
python3 -B -m unittest discover -s tools/nuclear_pid_tof -p 'test_*.py'
```

## 数値の読み方

粗いR幅0.05のfitは初期確認用。最終候補の中心／幅はnative幅0.01の非重複fitを使う。
射影の幅にはR内の中心変化が含まれるため、単純なfree-fit幅をintrinsic幅と呼ばない。
GL3によるR-smearing variance差引きは候補初期値の近似で、R平均モデルによるスペクトル再検証が必要。
連続曲線の間に未測定binがあっても、そのbinを自動的に支持済みへ変更しない。

purityは名目Gaussian分解とTAG shape転用を仮定した`S/(S+B)`。
現在のTPC±2σ・best-species窓と、新しい候補窓・単独Gaussian窓を区別して保存する。
固定shapeのyield誤差だけでは、核種帰属・背景model・TAG転用を含む最終誤差にはならない。
未支持の低Rやpoor fitの値を確定purityとしてつなげない。

Gaussianとlog-normalの最新比較・未完了の範囲は
[今回の比較まとめ](../../Todo/2026/10/14/comparison_gaussian_lognormal.md)を参照。
log-normalのmode／曲率幅をGaussian平均／標準偏差へ無条件に置き換えない。
ユーザー方針はまず両関数モデルの比較であり、TOF実測形状テンプレートは追加していない。
全域での不良適合・低R⁴Heの背景依存が残るため、本番PIDの変更・確定purity算出は未実施。
