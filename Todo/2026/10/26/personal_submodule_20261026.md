# 個人資料のサブモジュール管理

対象：`https://github.com/HypTPCIchikawa/star_local_submodule.git`。

## 構成

- 親解析の`star_local/`を個人リポジトリのサブモジュールにする。
- Todo/日報の原本を`star_local/Todo`・`star_local/analysisnote`へ移し、親の従来パスは相対symlinkで保持する。内容258件のSHA一致を確認済み。移動前backupも保持。
- 個人評価コード・設定・報告・小さな図は`star_local/snapshots/<元パス>`へ保存する。コードの原本は実行パスの依存を壊さないよう元位置に保持し、captureで同期する。
- 元の`rootfile`や大きな生成物はそのまま。サブモジュール直下のrootfile等は参照用symlinkで、Gitに入れない。
- 通常1ファイル1 MiB以下のPDF/PNG等を含める。出力JSONの要約は256 KiB以下。ROOT・ビルド生成物・ログ・大量数値配列を除外する。上限超過は消さず、一覧へ記録する。
- 親はURLとcommit参照のみ管理。`update=none`により通常の再帰cloneに個人リポジトリを要求せず、解析の実行依存にも加えない。

## 日常の操作

```bash
python3 star_local/scripts/sync_workspace.py capture
git -C star_local status --short
git -C star_local add -- Todo analysisnote snapshots metadata
git -C star_local commit -m 'Update personal records'
git -C star_local push origin main
```

Todo/日報は直接編集をGitが認識する。その他は元ファイルを編集後、captureして保存する。必要な時だけ親の`star_local`参照を更新してcommitする。詳細はサブモジュールのREADME。

## 検証・進捗

- 空のリポジトリへの接続を確認し、初期管理ツールを`main`へpush済み。
- サイズ・ROOT/バイナリ除外・既存パス保護・復元時全件事前検査・任意submoduleの再帰cloneの9 testsが成功。
- Todo/日報の元ファイル保持と元パスでの参照を確認済み。
- snapshot収集と照合が完了。9,179ファイル・284,013,346 bytes（270.86 MiB）、PDF 541件・PNG 1,024件。全SHA一致を確認し、今回の4He reportとβ–η PDFも登録した。サイズ等による除外1,191件の一覧を保存。続いて資料のpushと親のcommit参照を確定する。

開始・進捗記録：2026-10-26T12:57:44.931351+09:00。

## 資料の公開・最終設定

個人repoの`main`へ資料コミット `096eb5feb8cc07bf9178eaf01c3dbdcdca400ded` をpush済み。原本・snapshot・管理情報は計9363件、PDF 549件、PNG 1024件。ROOT・ログ・ビルド生成物は追跡0件。共有repoには`.gitmodules`、任意submoduleのcommit参照、既存パスaliasのignoreと利用説明だけを追加する。

初回のGit登録はGPFSへの多数オブジェクト書込みが遅かったため、専用の一時checkoutでcommit/pushし、元のsubmoduleへpackで取り込む。原本・snapshotは保持。Git上のsnapshotは改行変換を無効にしてSHAの再現性を維持する。

9 tests成功、snapshot全SHA一致、元の11件のtracked変更は完全一致。以降、Todo・日報は直接編集、その他は `python3 star_local/scripts/sync_workspace.py capture` で同期してから個人Gitへcommit/pushする。1 MiB超の図は`policy.json`で個別登録できる。

最終公開のSHA/時刻・親mainとの照合は `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/personal_submodule_20261026_001/completed.json` に保存する。

資料push・設定確認時刻：2026-10-26T13:17:09.766558+09:00。
