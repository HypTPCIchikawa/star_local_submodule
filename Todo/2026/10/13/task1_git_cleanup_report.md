# Task 1: Git 管理整理の結果（2026-10-13 JST）

[計画](task1_git_cleanup_plan.md)に従い、63件（通常ファイル62件・symlink1件、通常ファイル合計737,455 byte）をリポジトリ外へ隔離した。既追跡15件はGitHub現行treeから除く対象、残り48件は未追跡またはignore済みのローカル生成物である。内容は削除していない。

隔離先: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/_quarantine/star-analyzer_20261013/`

[全対象・元パス・SHA-256・復元先manifest](task1_quarantine_manifest.json)。隔離先にも同じ `manifest.json` を保存した。

| 分類 | 隔離数 | 既追跡数 |
| --- | ---: | ---: |
| Pythonキャッシュ | 32 | 3 |
| 旧単発生成/検査/debugコード | 18 | 10 |
| 旧diff snapshot | 1 | 1 |
| 個人作業場所へのsymlink | 1 | 1 |
| 偶発空ファイル | 5 | 0 |
| ROOT自動dictionary/生成物 | 3 | 0 |
| 旧解析backup | 1 | 0 |
| TeX補助ファイル | 2 | 0 |
| 合計 | 63 | 15 |

## GitHubから除く既追跡ファイル

- `.gen_kstar_yaml.py`
- `append_hMult.py`
- `append_hPvsY.py`
- `append_kstar_mass.py`
- `append_mixed_hist.py`
- `check_yaml.py`
- `debug_run.C`
- `diff_summary.txt`
- `list_keys.C`
- `test_config.C`
- `test_dawson.C`
- `tmp_Nuclear_id`
- `script/time_NYT_to_JST.pyc`
- `script/__pycache__/analysis_info_helper.cpython-38.pyc`
- `script/__pycache__/time_NYT_to_JST.cpython-38.pyc`

単発コードは旧固定設定への追記/検査/debugで通常build/run/docsからの参照がない。`test_dawson.C` 等には `EOF` /shell行混入もあった。symlinkはリンク自体だけ移動し、リンク先には触れていない。

## 再混入防止と保持対象

[.gitignore](../../../../.gitignore)へPythonキャッシュ、ROOT dictionary、TeX補助ファイル、backup、既知の偶発空ファイル/個人symlinkを追加した。`.py` / `.C` / `.tex` は一括ignoreしていない。

現行Maker・analysis・config、ROOT/PDF成果物、job履歴、KF比較/cut scan/purity study、summaryの元ソースを保持し、既存dirty treeは破棄していない。物理解析の計算内容を変更していないため、本TaskのSTAR再ビルドは不要。

## 検証と復元

62通常ファイルのSHA-256・サイズ・mode・mtimeが移動前後で一致。symlink文字列 `../../tmp_ana/Nuclear_id` が一致。63元パスが移動直後に空いていることを確認。`git diff --check` と `git check-ignore --no-index -v` は成功した。GitHub main は開始時HEAD `4324ceeba8b5b160390c49cb55f94caa8ac5e8de` と一致（ls-remote）。

Python実行後のcache再生成は許容し、以後ignoreする。現在変更済みだった `analysis_info_helper.cpython-38.pyc` の内容も保存した。

復元は隔離manifestの `path` を元repository相対パスとして使う。復元先に既存ファイルがあれば上書きせず比較する。隔離先 `README.md` に手順を記載した。Gitの過去commitにも既追跡旧ファイルが残る。

## GitHub反映

GitHubへの通常commit/pushは完了した。cleanup commit は [a9d63d6](https://github.com/star-tohoku/star-analyzer/commit/a9d63d62cdf0601db720de5e6b01965a9fd5b171)（`a9d63d62cdf0601db720de5e6b01965a9fd5b171`）。`origin/main` のremote HEADがこのcommitと一致することを `git ls-remote` で確認した。

開始時stagingが空であることを確認し、[stage候補](task1_stage_candidates.txt)の20パス（既追跡15件のuntrack/削除、`.gitignore`、Task 1記録4件）だけをcommitした。再生成したPythonキャッシュは `git rm --cached` で追跡を外し、現在のローカルファイルは保持する。既存変更やTask 2/3は含めていない。force push/履歴改変は行っていない。

実行logは `Todo/2026/10/13/task1_git_push.log` に保存する（ローカル記録、`.gitignore` の `*.log` 対象）。本節はcleanupのpush成功後に追補し、報告だけの通常commitで反映する。
