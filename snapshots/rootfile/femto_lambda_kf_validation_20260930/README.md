# KFを用いたΛ–原子核解析の実装検証（2026-09-30）

各解析の比較では、エネルギーごとに同じ先頭10,000入力イベントを指定しています。
元解析の基準は**analysis/anaLambdaNuclearId.C**であり、
未完成だった旧anaFemtoLambda_4He.Cではありません。

| ディレクトリ／ファイル | 内容 |
|---|---|
| `13p5/{d,t,3He,4He}.root` | 新しいStFemtoMaker。ビームエネルギー13.5 GeVのFXT。保存済みの高purity KFカットを使用 |
| `3p85/{d,t,3He,4He}.root` | 新しいStFemtoMaker。ビームエネルギー3.85 GeVのFXT。保存済みの高purity KFカットを使用 |
| `<energy>/legacy.root` | 変更していない元のHelix版anaLambdaNuclearIdを、元のmainconfで実行 |
| `<energy>/standalone.root` | 同じKFプリセットとイベント選別方針を用いた、単独KF Λ解析の対照試験 |
| `13p5/smoke100.root` | 事前の100イベント動作確認。10k比較には使用しない |
| `audit*.csv` | 解析全体の検証結果と収量の集計。完了状態は検証ログを参照 |
| `provenance/` | 実行・ビルド・テストのログ、設定の保持確認、ソースのスナップショット |

## Λ信号候補の質量選択

信号チャネルで原子核とのペアに使用するΛは、既存の4核種設定で
**質量中心から±3σ**に選択されています。
`lambdaSignalMean: 1.11596`、`lambdaSignalSigma: 0.00296`、
`lambdaSignalNSigma: 3.0`から、信号質量窓は
**[1.10708,1.12484] GeV/c²（両端を含む）**です。

包括的なΛ質量ヒストグラムと質量対k*は、信号窓で狭める前の
全質量範囲1.05–1.25 GeV/c²を保持しています。
比較図で使う**[1.10,1.14) GeV/c²**は表示・面積正規化の窓であり、
±3σの信号選択窓とは異なります。全質量のentries、比較図の表示窓内カウント、
信号窓内のペア数を区別してください。
この記述は既存設定の説明であり、解析・図・YAMLの変更は行っていません。

入力：`config/picoDstList/auau{13p5,3p85}_femtoLambda_validation_20260930.list`。
既存の出力は上書きしていません。初期のビルド・テスト・QAで失敗した試行もログに残し、
成功した実行は実装報告で明示しています。

- [実装報告](../../mdfiles/femto/implementation/femto_lambda_nuclei_kfparticle_20260930.md)
- [QA図と統計のまとめ](../../share/figure/femto_lambda_kf_validation_20260930/README.md)

これらは実装検証の出力であり、新たなpurity最適化の結果ではありません。
旧Helix解析と新しい高purity KF解析の収量が一致する必要はありません。
単独KF解析との候補の厳密な一致、および元のNuclearId選別の保持は、
それぞれ別の検証として確認しています。
