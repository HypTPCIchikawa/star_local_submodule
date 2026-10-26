# STAR personal workspace

`star-tohoku/star-analyzer`とは別に、Todo・日報・個人資料・評価コード/設定・小さなPDF/図を管理するリポジトリです。解析ROOTと大きな生成物は解析用ファイルシステムに保持します。

## 配置

- `Todo/`、`analysisnote/`：管理する原本。親解析ディレクトリの同名パスはここへの相対symlinkで、普段どおり編集できます。内容の変更はこのGitで直接認識されます。
- `snapshots/<元の相対パス>`：個人の評価コード・設定・報告・小さな成果物のコピー。既存プログラムの実行パスを保つため、元ファイルは親ディレクトリに残しています。編集したらcaptureを実行します。
- `metadata/snapshots.json`：元パス、親Gitの基準commit、種類、サイズ、SHA256。
- `metadata/excluded_large_files.json`：サイズ超過等で保存しないファイルの一覧。元ファイルは削除しません。
- `policy.json`：保存対象とサイズ制限。通常1ファイル1 MiB以下。ROOT、バイナリ、ビルド生成物、ログを除外します。出力ディレクトリの大きな数値配列は保存せず、設定・要約・由来情報を選択します。より大きい完成PDFを保存する場合は`extraLargeReports`に元パスと上限byteを明示します。

`rootfile/`、`config/`、`tools/`等の直下symlinkはローカル参照用です。Gitには入りません。これによりTodoの相対リンクから親の実データ・共有コードを引き続き参照できます。

## 取得と接続

親リポジトリのこのサブモジュールは`update = none`です。チームの通常の再帰取得では自動取得せず、解析ビルドも依存しません。利用者だけ明示的に取得します。

```bash
git -c submodule.star_local.url=git@github.com:HypTPCIchikawa/star_local_submodule.git submodule update --init --checkout star_local
git -C star_local switch main
git -C star_local config core.hooksPath hooks
python3 star_local/scripts/sync_workspace.py mount
```

既存のTodo/日報が実ディレクトリとしてある場合、mountは上書きせず停止します。バックアップと内容照合をしてから移行してください。rootfileはこのGitから復元できません。

## 更新

親解析ディレクトリから実行します。Todo・日報は直接管理、他の個人コード等はcaptureで同期します。captureは親の解析ファイルを変更しません。

```bash
python3 star_local/scripts/sync_workspace.py plan
python3 star_local/scripts/sync_workspace.py capture
git -C star_local status --short
git -C star_local add -- README.md policy.json .gitignore scripts hooks verification Todo analysisnote snapshots metadata
python3 star_local/scripts/sync_workspace.py check-index
git -C star_local diff --cached --stat
git -C star_local commit -m 'Update personal analysis records'
git -C star_local push origin main
```

親で共有するのはサブモジュールのcommit参照だけです。必要な時に`git add star_local`で参照を更新して通常commit/pushします。共有解析コードの未commit変更を一括stageしないでください。

## 復元と保護

`python3 star_local/scripts/sync_workspace.py restore`はsnapshotの元ファイルがない場合だけ復元します。既存の異なる内容やSHA不一致は、書込み前の全件確認で停止します。共有コードへの手元overrideも含むため、既存checkoutとの相違は手動で比較・統合してください。復元後、コードは親の元パスから実行します。

`hooks/pre-commit`はstage済み内容を調べ、ROOT・大きなファイル・バイナリ・symlinkの混入を拒否します。追加cloneでも`core.hooksPath`を設定してください。完全なROOT/ログ等の保管はこのGitの役割には含めません。

```bash
python3 -m unittest discover -s star_local/verification -v
```
