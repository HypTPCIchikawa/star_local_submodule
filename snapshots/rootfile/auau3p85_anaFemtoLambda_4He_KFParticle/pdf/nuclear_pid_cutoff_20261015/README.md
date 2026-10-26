# 原子核PID cutoff研究（2026-10-15）

3.85 GeV、2021 P24iyの保存済みALL／TOF-tag QAから、低RのTPC-only PIDが使える範囲を調べています。[計画・実施報告](../../../../Todo/2026/10/15/README.md)を参照してください。

**2026-10-17（日本時間）の再開分は完了しました。** 80実行の保存済みfit・感度集計から最終PDFを作成し、最終範囲の独立YAML・mainconfを追加して20イベント接続検証を実施しました。旧本番mainconfは切り替えていません。[実施報告](../../../../Todo/2026/10/15/report_nuclear_pid_cutoff_20261015.md)と、履歴として[10/15終了・引継ぎログ](../../../../analysisnote/20261015/summary20261015.md)を参照してください。

| 出力 | 内容・状態 |
|---|---|
| `baseline_001/` | 元のmainconf、実際のTPC中心・幅、参照ソースの凍結。後続の新設定と分離。 |
| `all_{gaussian,lognormal}_{linear,exponential}_nominal_001/` | d〜⁴Heを跨ぐ広い窓での条件付きALL fit。背景モデル不適合の診断結果でありcutoff採用用ではない。 |
| `summary_nominal_001/` | 上記のpurity／残差／支持判定。`purity_overview.pdf`は3ページ。 |
| `spectrum_views_nominal_readable_002/spectra.pdf` | 全4広窓モデルの生スペクトル・成分・背景。全dE/dxとfit窓、linear／log、24ページ。旧表示版を保持。 |
| `all_{d,t,he3,he4}_{gaussian,lognormal}_{linear,exponential}_nominal_001/` | 核種ごとの局所窓。応答はTOF-tagから固定転用。省略tail上限も記録。 |
| `summary_local_nominal_001/` | **支持判定に集計バグがあった旧版**。purity値を変更せず002へ修正。採用判定に使わない。 |
| `summary_local_nominal_002/` | 修正後の固定応答・局所fitの集計とpurity／ALL COLZ／収量の3ページPDF。最終cutoffではない。 |
| `local_contract_prepare_001/` | 局所窓の入出力準備。これ自体は実行完了の証拠ではない。 |
| `validation_domain_001/` | batch一致build、PID domain・既存TOF・校正の回帰ログと人工fixture。初回の誤引数による失敗ログも保持。 |
| `validation_fit_001/` | 広窓・局所固定・局所freeのROOT5数値テスト、ソース、実行ログ。 |
| `independent_gate_counts_001/` | fitとは独立した全4核種×160区間のTPC候補数と保存bin境界による上下限。fit不成立でも統計を0にしないための共通入力。 |
| `all_<species>_<response>_<background>_nominal_free001/` | ALLで目的核種の中心・幅を自由にした比較。選択窓自体は現行のまま。採用判断は全モデルの集計を参照。 |
| `summary_free_sensitivity_001/` | 最終nominal 16＋感度64実行の集約。核種ごと20条件でモデル低側・数値支持・残差を評価。 |
| `adopted_ranges_001/` | 下限・除外域・保持率のJSONと描画用TSV。研究上の使用候補範囲であり`productionEnabled=false`。 |
| `validation_domain_veto_001/` | 相対出力パスの二重解決による旧テストharness失敗ログ。保持。 |
| `validation_domain_veto_002/` | harness修正後のbatch一致build、絶対／相対出力で各1359 checks PASS。 |
| [views_final_001/](views_final_001/README.md) | 最終purity・使用範囲・COLZ・候補数、中心・幅、代表fit。全6 PDF、全頁読み戻し・目視QA・SHAを保存。 |
| [validation_resume_audit_20261017_001/](validation_resume_audit_20261017_001/README.md) | 80実行のcoverage、1062件のSHA、最終YAMLとの一致、旧本番設定・校正不変の独立監査。 |

`views_final_001/`は10/15中断時には未作成でしたが、10/17に完成しました。既存のnominal／途中経過の図は履歴として保持し、最終版はこの新しい図を参照してください。

低Rのpurityは実際のTPC±nσ＋best-species受理域での目的核種／全成分を計算しています。TOF-tag内のpurityではありません。R幅0.01が基本で、0.05幅の表は信号・総数を加算してから比を求めています。単純な平均ではありません。

モデルの低側と残差・数値収束・核種帰属を分けて保存します。未評価はpurity=0ではありません。TAGはALLの部分集合なので独立尤度として二重加算していません。ここでのtrack countは入力イベント数、CFペア数、真の効率とは異なります。

高Rの全核種TOF実装の20イベント検証は [別のtestsフォルダ](../../tests/allNucleiTof20261015/README.md) に保存しています。本研究の元ROOT、過去のPDF、既存本番設定は保持しています。新規job投入やhaddは実施していません。

最終下限・vetoまで組み合わせた独立プロファイルのbuild・20イベント・閉ROOT metadata・回帰検証は [tests/purity50_20261015](../../tests/purity50_20261015/README.md) です。前の高R TOFのみの検証と区別してください。
