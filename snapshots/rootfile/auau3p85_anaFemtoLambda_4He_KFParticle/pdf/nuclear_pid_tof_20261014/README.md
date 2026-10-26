# 2026-10-14: TOFタグを用いた原子核PID study

対象は3.85 GeV・2021 P24iyの保存済み高統計ROOT。
原ROOT・旧PDF・現在の解析PIDは変更していない。
[計画と最新報告](../../../../Todo/2026/10/14/README.md)を参照。

| 保存先 | 内容・採用状況 |
|---|---|
| `inputs_001/` | ALL／3タグ／remainingのnative cache、軸・計数・誤差・来歴監査、全dE/dxのCOLZ・射影図。完成。 |
| `tag_smoke_prepare_001/` | 最初のGaussian試験の準備記録。fitは実行していない。 |
| `tag_smoke_linear_001/` / `tag_smoke_exponential_001/` | R=0.80〜0.90、幅0.05のGaussian応答＋背景比較。未採用。 |
| `tag_smoke_linear_views_001/` / `tag_smoke_exponential_views_001/` | 上記をlinear／log・全dE/dx・zoom・残差で表示。再fitなし。 |
| `tag_native_linear_001/` / `tag_native_exponential_001/` | R=0.30〜2.50、幅0.01のGaussian TAG fit。各660スペクトル。未採用。 |
| `tag_native_linear_core_views_001/` / `tag_native_exponential_core_views_001/` | 上記のR=0.80〜0.90、30スペクトルの表示専用PDF。 |
| `curves_native_001/` | 数値支持域のauditとlog(R)多項式／滑らかな区分3次の候補。⁴Heは関数化条件を満たさず曲線なし。本番PIDではない。 |
| `native_roi_residual_audit_001/` | ピーク近傍±2／±3幅と全fit窓の残差。0.1／0.5／1 keV/cmの再集約による診断で、再fitなし。 |
| `all_gaussian_contract_prepare_001/` | ALL fit用の入出力ガード・窓・GL3の準備試験。ROOT fitとpurity算出はまだ実行していない。 |
| `plot_compile_001/` | ROOT5での新表示マクロのコンパイル試験。物理fitではない。 |
| `tag_lognormal_native_smoke_{linear,exponential}_001/` | Gaussianと同じR=0.80〜0.90の30 nativeスペクトルでlog-normal応答を比較。未採用。 |
| `tag_lognormal_native_full_{linear,exponential}_001/` | R=0.30〜2.50の全660スペクトルをlog-normalで比較。未採用。 |
| `tag_native_{linear,exponential}_core_views_002/` | Gaussianの紙面内修正版。数値は001版と同じ。 |
| `tag_lognormal_native_smoke_{linear,exponential}_views_002/` | log-normalの紙面内修正版。数値は001版と同じ。 |
| `nominal_gaussian_guides_004/` | 生ALL／TOFタグに名目fit点・失敗点・旧／現PID参照を重ねた正式COLZ 5ページ。⁴He支持曲線は未確定。 |
| `comparison_gaussian_lognormal_004/` | 同じnativeデータの正式比較PDF 4＋24ページ、1760行の表、入力照合と来歴。以前の001〜003も保持。 |

## 現段階の重要点

TOFdはd／⁴He互換のタグで、純粋な⁴Heではない。TAGはALLの部分集合であり、独立尤度として足さない。
高統計のピーク近傍にGaussian／log-normalからの形のずれが残る。
低統計でχ²が小さいことや、fitが数値収束したことだけでは核種帰属・purity・PID校正を確定しない。

ユーザーの方針は、まずGaussian／log-normalの関数モデルを比較すること。
TOF形状の実測テンプレートへの拡張は行わない。
log-normal試験は別source・別出力に保存し、mode／局所曲率の幅とGaussianの平均／標準偏差を区別する。
全域では一括log-normal採用を支持せず、低R⁴Heも背景依存が大きい。
モデル検証が済むまでは、本番PIDへの反映は行わない。最終関数・支持域・ALL purity検証は未完了。
詳しい比較と次の検討事項は[今回比較まとめ](../../../../Todo/2026/10/14/comparison_gaussian_lognormal.md)を参照。
