# 保存fitの品質監査（2026-10-19）

旧補間gateのfitを読取専用で再集計した診断。元多項式gateのpurity再認証ではない。
元ROOT・fit・PID設定・旧PDFは変更していない。新規fitなし。

対象: 100実行、14840 fit行。低R80実行、高Rnative16実行、He4 coarse4実行。
χ²/ndfは有限R投影のtarget-local窓のPearson値。品質screenをPの判定より先に適用。
0.05区間は全native fitの最大χ²/ndfを使う。支持済み条件のPはΣS/Σ(S+B)。
高RのHe4 nativeとcoarseは別stageで比較し、良い点の継ぎ接ぎはしない。

## 低R nominalのχ²/ndf分布

| Species | N | Median | 90% quantile | Max |
|---|---:|---:|---:|---:|
| d | 634 | 7559.405 | 221350.266 | 6871590.993 |
| he3 | 634 | 535.354 | 1708.269 | 17304.298 |
| he4 | 634 | 82.531 | 1218.512 | 2333366013.565 |
| t | 634 | 811.769 | 1715.413 | 49335.288 |

## 主な注意

不良modelを除外しても、応答2family・背景2family・nominal2条件以上が残らなければ未評価。
2条件による両family被覆は2×2全組の保持とは異なる。各行のretainedLabelsを必ず確認する。
χ²cutの採用値をこの監査だけで決定しない。隣接binと除外fitの局所図を確認する。
飛び地はfit未収束・帰属・応答精度・残差等で発生し、真の低purityと区別する。

## 出力

- native_quality.tsv: 失敗を省かない全fitのχ²・残差・旧支持理由。
  native支持は行単位の診断。旧0.05区間の集約支持はgrouped_quality.tsvを参照。
- focus_native_quality.tsv: 下限、飛び地とその両隣。
- grouped_quality.tsv: 0.05区間ごとの元支持とχ²。低Rの旧集約数値・理由を照合済み。
- threshold_scan.tsv: 1.5/2/3/5/10/20の品質限定envelope（旧gate参考）。
- distributions.json: 低R、高Rnative、高Rcoarseを分けた分位値。

注目域の理由集計（同区間・感度条件の重複を含み、独立事象数ではない）:

```json
{
  "fitted_meanError_over_sigma": 36,
  "fitted_sigmaError_over_sigma": 37,
  "free_target_shape_or_identity_screen": 428,
  "numericalValid": 106,
  "gate_integral_residual": 218,
  "local_gate_absolute_residual": 251,
  "fitted_window_yield_too_small": 17,
  "target_yield_boundary": 7,
  "available": 10,
  "axisCovered": 10,
  "gateActive": 10,
  "identitySupport": 10,
  "purity_unavailable": 10,
  "selectionCovered": 10,
  "unsupported_material_nuisance": 10
}
```
