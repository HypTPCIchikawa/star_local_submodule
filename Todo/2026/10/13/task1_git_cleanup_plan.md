# Task 1: Git 管理整理の実装計画（2026-10-13 JST）

[依頼](STAR_Todo_20261013.txt)に従い、内容・参照を確認した不要な中間生成物と単発デバッグ用ファイルを整理する。解析ソース・設定・ROOT/PDF成果物・投入履歴・再現用studyツールは保持する。

## 隔離と対象の判断

隔離先は `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/_quarantine/star-analyzer_20261013/`。元の相対パスと内容を維持し、SHA-256、Git追跡状態、symlinkのリンク文字列をmanifestに保存する。現在変更済みのキャッシュも現在の内容を保存する。復元時は既存ファイルを上書きしない。

| 分類 | 対象・根拠 |
| --- | --- |
| Python生成物 | 既追跡 `.pyc` 3件と `__pycache__`。対応する `.py` は保持。 |
| ROOT生成物 | root直下 `AutoDict_vector_LksPoint_...`。元studyソースは保持。 |
| TeX生成物 | `tools/summary_4He*.aux` のみ。`.tex`・PDFおよび `rootfile/` のレポート一式は保持。 |
| 偶発空ファイル | `0.52`、`1`、`Get`、`Write`、`⚠️`。全て0 byte、参照なし。 |
| 旧単発コード | root直下 `.gen_kstar_yaml.py`、`append_*.py`、`check_yaml.py`、`debug_run.C`、`list_keys.C`、`test_config.C`、`test_dawson.C`。固定パスの旧YAML追記/ローカルテストで通常フローから参照なし。 |
| 未追跡の旧単発コード | `make_femto_hist.py`、`tools/debug_crash.sh`、`tools/fix_yaml_th2f*.py`、`tools/update_hist_yaml.py`、`tools/val.py`、`tools/generate_3p85_joblists.py`。固定設定への一度きりの書換/検査で通常フローから参照なし。最後の生成scriptは旧2018/2019設定も含み現在方針と不一致。複数scriptには末尾 `EOF`/shell行混入あり。 |
| 生成debug wrapper | `script/singularity_run_debug.sh`。旧 `tools/debug_crash.sh` の `sed` 生成物、10イベント固定。通常Singularity wrapperは保持。 |
| 旧記録・backup | `diff_summary.txt`、`analysis/anaLambdaNuclearId.C.bak`。履歴/旧ソース内容を隔離保存。 |
| 個人symlink | `tmp_Nuclear_id -> ../../tmp_ana/Nuclear_id`。参照なし。リンクだけ移動し、リンク先のデータは触らない。 |

曖昧なものは保持する。特に `common/macro/` のKF比較/cut scan/purity study、`tools/kfparticle_summary/`、新規Makerソース、`job/run/` の監視/統合処理/投入記録は保持する。

## 実装と検証

1. 対象の追跡状態・ハッシュ・リンク文字列を保存する。
2. 元相対パスを保ち隔離し、移動後ハッシュ/リンク文字列一致を検証する。
3. `.gitignore` にキャッシュ、dictionary、TeX補助ファイル、backupを追加する。`.py`/`.C`/`.tex` は一括ignoreしない。
4. `git diff --check` と追跡削除対象を確認する。物理コード変更はなく、本TaskのSTAR再ビルドは不要。
5. cleanup削除・`.gitignore`・本Task記録だけを明示stageし、既存dirty tree全体を `git add -A` しない。
6. 親担当と調整して通常commit/pushでGitHub現行treeから不要物を除く。履歴改変・force pushはしない。

## 出力

- [計画](task1_git_cleanup_plan.md)
- 隔離/復元manifest: `task1_quarantine_manifest.json`
- 結果報告: `task1_git_cleanup_report.md`
- GitHub整理用stage候補: `task1_stage_candidates.txt`
