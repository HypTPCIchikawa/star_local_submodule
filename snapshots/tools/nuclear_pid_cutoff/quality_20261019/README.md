# 2026-10-19 保存fitの品質監査と診断図

[当日Todo](../../../Todo/2026/10/19/README.md)に対応する、旧fitの読み取り専用診断です。
元ROOT・fit・旧PDF・本番PIDを上書きしません。新しい多項式gateのpurity認証ではありません。

## ファイル

| ファイル | 役割 |
|---|---|
| `audit_fit_quality.py` / `audit_settings.yaml` | 保存済み低R80実行・高R20実行の品質再集計。独立gate countのSHAと640行coverage、旧集約の再現を検査。 |
| `PlotQualityDiagnostics.C` | 共通painterを使い、広域／対象核種近傍、linear／logを表示。失敗fitも省かず、χ²と旧支持理由を注記。 |
| `render_quality_diagnostics.py` / `render_quality_settings.yaml` | 旧除外域の当該・隣接RについてROOT5で診断PDFを生成。全頁readback・入力SHA不変を確認。 |
| `test_audit_fit_quality.py` | 品質をpurityより先に判定し、モデル不足を未評価とする回帰試験。 |
| `test_render_quality_diagnostics.py` | 局所表示範囲、失敗fit、欠損χ²の表示を確認。 |

χ²は保存値の`Pearson / conditionalNdf`です。自由shapeのパラメータを完全に差し引いた
厳密な適合度検定のndfと同じとは限りません。現在の閾値scanは基準案比較用で、採用値ではありません。
native支持は各行の診断、旧0.05 GeV/c集約の支持は`grouped_quality.tsv`を参照します。

## 実行

リポジトリ直下で実行します。`NEW`は未使用の出力名に置き換えます。

```bash
python3 -B -m unittest discover -s tools/nuclear_pid_cutoff/quality_20261019 -p 'test_*.py' -v
python3 -B tools/nuclear_pid_cutoff/quality_20261019/audit_fit_quality.py \
  --settings tools/nuclear_pid_cutoff/quality_20261019/audit_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/audit_NEW
python3 -B tools/nuclear_pid_cutoff/quality_20261019/render_quality_diagnostics.py \
  --settings tools/nuclear_pid_cutoff/quality_20261019/render_quality_settings.yaml \
  --audit-table rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/audit_NEW/native_quality.tsv \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/representative_views_NEW
```

実行はSTAR SL24y／ROOT5、`script/singularity_root_batch.sh`経由です。
描画は保存count・モデルのまま行い、再fit・表示normalizationをしません。
モデル線は保存済みfit mask内だけに描き、窓外の未計算binを予測ゼロとして結びません。
局所t図はt近傍だけを拡大し、大きなdピークを含む全体図と分けます。
2026-10-19追加指示により`plots.linearFocusOnly: 1`を明示し、左linearは全ページで対象核種近傍を表示します。
広域ページの右logは0–120 keV/cmを維持し、局所ページの右logは対象近傍を表示します。
再fit・count変更はなく、パネルごとのfocus方針をログと回帰試験で検査します。
最新の`representative_views_004`はこの表示方針で68頁readback済み。003以前はそのまま保持しています。
現設定はd28頁、t12頁、⁴He28頁の計68頁を予定し、実完了は`completed.json`で判定します。
旧veto内のnative隣接点だけでなく、veto区間そのものの前後0.05-R区間を含めます。
dの[1.70,2.00)はTPC-only上端に達するため、同じ試料条件の右隣は存在しません。
途中出力だけでは完了としません。

## 関連する独立study設定

`../logpoly_20261019_config_audit.py`は元のlog10(R)多項式の係数、幅の正値、
5入口のstudy／切り戻し設定を監査します。calibration補間と旧vetoを外すのは独立studyだけです。
係数の算術検査は、C++runtimeや新gateのpurity検証を代用しません。
