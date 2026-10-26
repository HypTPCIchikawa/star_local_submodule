# 再評価fitパラメータと既存PID関数の比較：native診断版

2026-10-17に作成した18ページの途中診断版です。低R ALLの4核種、高R TOFタグの4核種を比較し、旧設定比較9ページへ追加しています。新fitや係数変更は行っていません。

⁴He高Rの0.05 GeV/c幅による全域独立fit、および採用済み低R範囲だけの再表示は、後続の`../pid_parameter_comparison_final_20261017_001/`に追加しています。通常はそちらのREADMEと全面配置PDFを参照してください。この診断版は来歴保持のため上書きせず残しています。

`source/`、`context.json`、`fit_points.tsv`、`coverage.tsv`、`completed.json`に固定入力とSHA256を保存しています。旧設定は以前の比較reportのsnapshotから取得しており、現在の本番設定を旧設定として再利用していません。fitの初期値は既存タグfitやprovider seedであり、旧PIDの選択中心・幅そのものではありません。

Gaussianの平均・標準偏差とlog-normalのmode・局所曲率幅は異なる定義です。有限R幅の投影応答なので、検出器固有の校正幅として採用していません。色付き／灰色の意味、TOF条件付き試料の制限、過去の⁴He `[3,4]` 窓と今回の狭いD窓の違いは最終版READMEに説明しています。
