# チーム共有Gitの整理・本番更新の公開

実施日：2026-10-26（日本時間）。ユーザーの実行指示と追加指定に従い、HypTPCIchikawa / Yudai Ichikawa が追加した個人資料等を履歴から照合して整理した。ローカル実体は保持し、過去のGit履歴は書き換えていない。

## 公開結果

- リポジトリ：`star-tohoku/star-analyzer`、branch：`main`。
- コミット：`31b7d126830b978c9ad95f23b5deb730ed68e472`。通常のfast-forward push：`62c5f4c → 31b7d12`。リモートmainのSHAも照合済み。
- 最新チームmainを基に独立worktreeで検証した。CandidateFemto/BDT等のチーム更新を保持。
- 元作業場所のmain/indexも同期した。hard resetは使わず、無関係な変更・個人資料・ROOT/PDFを保存した。

## 共有した内容

- 採用4HeのTPC中心3次・幅2次のlog10(R)多項式、確定係数YAML、厳密な0.61<R<5 GeV/c、核種別R/TOF条件。
- d/t/3Heの `NuclearIdDeDxVsMom.h` は係数・関数とも変更なし。全核種でR<2はTOF不要、R≥2は有効TOF+m²。4Heはdのm²窓。
- 物理3Heと3He-as-4He control、Λ/4Heのβ–η・β–Y、4HeのR/TOF区分別β–η・β–Y・pT–Y。
- 任意のpair PID TH2/8次元Sparse、その設定依存・回帰テストと全ヒスト説明。
- SUMS入力path/countの正規化、入力・終了・出力の検査、汎用の監視・統合・保存helper。
- 校正YAMLの出典pathは同定用metadata。実行時には出典ROOT・fitコードを開かず、必要係数は配布YAMLに含む。

## 管理解除した個人資料等（40件）

追加コミットのauthorがYudai Ichikawaであることを全件確認した。

- `Todo/20261013/task1_git_cleanup_plan.md`
- `Todo/20261013/task1_git_cleanup_report.md`
- `Todo/20261013/task1_quarantine_manifest.json`
- `Todo/20261013/task1_stage_candidates.txt`
- `analysisnote/20260904/summary20260904.md`
- `analysisnote/20260907/summary20260907.md`
- `mdfiles/configuration/unused_lambda_cut_cleanup_20260908.md`
- `mdfiles/femto/plans/plan_femto_lambda_nuclei_kfparticle_20260929.md`
- `mdfiles/kfparticle/comparisons/kfparticle_lambda_qhu1_cut_comparison_20260908.md`
- `mdfiles/kfparticle/comparisons/lambda_four_histogram_cut_summary_20260914.md`
- `mdfiles/kfparticle/comparisons/lambda_helix_kf_comparison_10000_20260907.md`
- `mdfiles/kfparticle/comparisons/lambda_helix_kf_imp5_comparison_10000_20260908.md`
- `mdfiles/kfparticle/comparisons/lambda_helix_kf_imp5_cut_scan_20260913.md`
- `mdfiles/kfparticle/comparisons/lambda_kf_purity70_refinement_20260927.md`
- `mdfiles/kfparticle/comparisons/lambda_kf_purity_study_20260927.md`
- `mdfiles/kfparticle/comparisons/lambda_kf_upstream_purity_20260928.md`
- `mdfiles/kfparticle/plans/plan_kfparticle_lambda.md`
- `mdfiles/kfparticle/plans/plan_kfparticle_lambda_full_reconstruction.md`
- `mdfiles/kfparticle/plans/plan_kfparticle_pro_integration.md`
- `mdfiles/kfparticle/plans/plan_lambda_kf_purity70_refinement_20260927.md`
- `mdfiles/kfparticle/plans/plan_lambda_kf_purity_study_20260927.md`
- `mdfiles/kfparticle/plans/plan_lambda_kf_upstream_purity_20260928.md`
- `rootfile/auau13p5_anaLambda_KFParticle/mdfiles/README.md`
- `rootfile/auau13p5_anaLambda_KFParticle/mdfiles/rootfile_move_map_20260908.tsv`
- `rootfile/auau13p5_anaLambda_KFParticle/mdfiles/rootfiles_SHA256SUMS_20260908`
- `rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/best50_replay10k_analysis.yaml`
- `rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/configs/hist.yaml`
- `rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/configs/frozen/best50_replay10k_analysis.yaml`
- `rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/configs/hist.yaml`
- `rootfile/femto_lambda_kf_validation_20260930/README.md`
- `rootfile/femto_lambda_kf_validation_20260930/acceptance_split/README.md`
- `rootfile/femto_lambda_kf_validation_20260930/pair_mt/README.md`
- `tools/summary_4He.tex`
- `share/figure/femto_lambda_kf_validation_20260930/README.md`
- `tools/add_missing_hists.py`
- `tools/apply_patch.py`
- `tools/create_fit_macro.py`
- `tools/fix_missing_hists.py`
- `tools/merge_hist_yaml.py`
- `tools/patch.py`

必要な旧ヒスト軸JSON fixture 1件は `tests/fixtures/femto_lambda_legacy_histogram_manifest.json`へ移動し、テスト・技術資料の参照を更新した。

## 除外・保存

`.gitignore`にTodo、日報、生成物、個人報告/計画、今回のPurity fit/作図専用ディレクトリ、監視状態、日付付きの単発jobスクリプト等を明記した。trackedだったものは別途indexから管理解除した。共有コードやMarkdown全体の一括除外は行っていない。

今回のPurity推定エンジン・fit図・report PDFは非公開。本番係数・runtime設定と必要なQA/回帰テストは配布対象とした。

元の個人資料は元位置に保持。旧Todo4件は以前に移動済みの現在位置を照合した。更新前の資料・diff・公開差分・由来確認・同期計画・検証ログも下記のローカルauditに保存した。

`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/git_publish_20261026_001/`

主な記録：`before_tracked.patch`、`before_status.txt`、`backup_manifest.json`、`local_backup/`、`removed_paths_authors.json`、`published.patch`、`published_name_status.txt`、`validation.json`、`primary_sync_result.json`。

保持した無関係な変更：AGENTSのローカル規則、AI規則の手元追加（最新チーム規則へ3-way統合）、4核種highpurity別名の手元study設定、5つの既存Lambda作図ツール変更。その他の未追跡studyソース・データも変更していない。

## 検証

- 最新mainベースの新規checkoutでSL24y / sl73_x8664_gcc485 / ROOT5.34.38 / GCC4.8.5の全ビルドが成功。
- 設定/入口：9 tests中1件は生成job XMLがないcheckoutのためskip。旧ヒスト契約632比較、k*–mT契約が成功。
- batch入力7、終了通知10、統合監視21、ROOT保存17 testsが成功。検証では実際のメールを送らない。
- ROOT：PID 14,147、β 53、acceptance区分734,505、Sparse 53,234,314、pair QA 129,387,488、pool 65,063,426 checksが成功。
- QA OFF fixtureの参照表記をONと統一した。また過去のHe4領域を前提にした期待値を、現行0.61<R<5へ更新して成功を確認した。物理ソースのカット変更ではない。
- 2021P24iyの実20イベントは20/20処理・再構成し正常完了。閉ROOT再読取りで通常1,786ヒストとSparse、ledger20、Λ質量/β各1 entry、全acceptance区分のbin整合を確認。
- この20イベントの採用4Heは0件。4Heが非ゼロでfillされることとTOF境界は上記の候補fixture/acceptance検証で確認した。

元作業場所の本番ライブラリも再ビルド成功（exit 0）。未公開の実験Makerを含めず、追跡済みソースを明示した。最初の再ビルドはwrapper内部の引数展開で停止したため、空白入りmake変数のquoteを修正して再実行した。公開checkoutのbuild成功とは別に記録する。

記録時刻：2026-10-26T11:53:19.102792+09:00。

最終確認：2026-10-26T11:55:15.719033+09:00。リモートmainと元作業場所の公開ファイルが全件一致。Todo/analysisnote/rootfile/shareの追跡0件。ローカルバックアップ132件のSHA/サイズ一覧を `complete_backup_manifest.json` に保存。無関係なWIPは保持。全作業完了。
