# 再開時の独立読み取り監査（2026-10-17 JST）

前回の最終集計と、新しく作成した最終cutプロファイルを独立に照合する記録です。新しいfit・PicoDst解析・ROOTの変更は行いません。この監査自身が作成するファイルは、この専用ディレクトリ内の新規レポートだけです。

## 再現コマンド

リポジトリ直下で実行します。既存レポートの上書きは拒否するため、再検査時は新しいファイル名を指定してください。

```bash
python3 rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/validation_resume_audit_20261017_001/replay_audit.py \
  --project /gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer \
  --report rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/validation_resume_audit_20261017_001/report.json
```

## 検査内容

- `adopted_ranges_001`と`summary_free_sensitivity_001`の保存済み生成物SHA。
- 全80実行の完了manifest SHAとcritical生成物SHA、4核種×4モデル×5条件のcoverage。
- 160 native区間／実行、集約128行・モデル詳細2560行。
- 保存された区間ごとの候補数から保持率を再計算し、モデル低側purityの最小値と使用候補JSON／TSVを照合。
- 新しい`purity50_20261015`核PID YAMLの全4種の下限・vetoをJSONと照合。
- 全4核種R≥2のTOF窓、旧R<2.5上限と明示的な⁴He拡張、TPC-onlyとの競合がないことを確認。
- 新mainconfがcanonicalな2021 P24iyの解析情報と既存出力名を使い、参照先が存在することを確認。
- 旧本番mainconf・核PID YAML、PID係数header・校正YAMLが研究開始時の凍結SHAと一致することを確認。

元の約10GB ROOTは再ハッシュしません。各実行に保存された入力保存・保護snapshot検査の完了情報を検証します。100MBを超える予期しない監査対象も拒否します。

これは保存物と設定の整合性監査です。実行境界の検証は別のdomain試験に依存し、モデル低側50%以上を全系統誤差やCFペアごとのpurity保証へ拡張しません。候補保持率を真の信号効率と解釈しません。

実測結果は[report.json](report.json)を参照してください。
