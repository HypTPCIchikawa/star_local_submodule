# 2026-10-19 Gaussian core＋tail局所比較

[当日Todo・進捗](../../../Todo/2026/10/19/README.md)に対応する独立study。
旧Gaussian／log-normal・本番PID・旧ROOT/PDFを変更しません。

## ツール

| ファイル | 内容 |
|---|---|
| `NuclearCoreTailFitCommon.h` | 元Gaussian engineの重み付きnodesをcore／tailへ展開。元engineは不変更。 |
| `SelfTestCoreTail.C` / `coretail_math_settings.tsv` | Gaussian極限・正規化・bin積分・R重み・幅変換の数学試験。 |
| `run_selftest.py` | 明示mainconfのROOT5で数学試験。階層を保つsnapshot・SHA・ログ・完了manifestを保存。 |
| `prepare_study.py` | 元多項式・TPC-only R<2の独立入力設定を作成。fitや本番切替はしない。 |
| `FitCoreTailPilot.C` | 同じ旧ALL・局所窓・旧gateで応答を比較。失敗fitでも独立raw countを保存。 |
| `run_coretail_pilot.py` / `pilot_settings.yaml` | 代表4R・48実行64fitのprepare／run／check。生成入力もSHA保護。 |
| `run_pilot_suite.sh` | 残りの未実行pilotを順次実行し、最後にcheck。tmux用の限定launcher。 |
| `render_pilot.py` | 完了済み48実行64fitから12冊64頁を描画。χ²と数値収束を表示し、帰属・purity認証は未評価と明示。 |
| `test_*.py` | 元多項式参照設定、相対パスsnapshot、改変検出、旧窓保存、二重実行防止等の回帰試験。 |

## 応答の定義と制限

核種のdensityは`(1-f) Gaussian(core) + f Gaussian(tail)`。
最初の比較では`f<0.5`、幅比`rho>1`、seed幅で定義する初期offsetを固定profileとする。
offsetはseed幅から設定したdE/dx差のまま、fit幅変更時にも固定。
核種yieldはcore／tail共通の1列で、coreの中心・幅を調整する。

core幅は本番PID校正幅や混合RMSではない。profile走査の最良値は独立な適合度検定ではない。
同じ旧gateによる参考purityは、元多項式の新gateのpurity認証に使わない。
対象のみfreeと、局所競合核種もfreeの場合を分け、応答変更とnuisance変更を混同しない。

## 数学試験

```bash
python3 -B tools/nuclear_pid_cutoff/coretail_20261019/run_selftest.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_logpolyLowRReference20261019.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_coretail_20261019/math_selftest_NEW
```

`math_selftest_002`はROOT5で5489項目PASS。001は名前空間の曖昧性によるcompile失敗記録として保持。
数学試験は実データfit・PID校正・purityを認証しない。

## pilot実行

`NEW`は未使用名に置き換える。まずprepareし、tailなしcontrolの旧fit再現を確認する。

```bash
python3 -B tools/nuclear_pid_cutoff/coretail_20261019/run_coretail_pilot.py prepare \
  --settings tools/nuclear_pid_cutoff/coretail_20261019/pilot_settings.yaml \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_purity50_20261015.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/coretail_pilot_NEW
python3 -B tools/nuclear_pid_cutoff/coretail_20261019/run_coretail_pilot.py run \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/coretail_pilot_NEW \
  --job d_linear_target_only_gaussian_control
bash tools/nuclear_pid_cutoff/coretail_20261019/run_pilot_suite.sh \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/coretail_pilot_NEW
```

最後のlauncherはtmux内で実行できる。`suite.log`と各実行`run.log`を保存する。
`pending`は完了済みのSHAを確認してskip。失敗・途中実行を同じ出力へ自動retryしない。
完了はpilot直下の`completed.json`で判定。raw count/error/flow・成分和・閉じたROOT読み戻し・
旧fitとのGaussian再現・全入力SHAを検査する。farm job・メール・hadd・Git pushは実行しない。

```bash
python3 -B tools/nuclear_pid_cutoff/coretail_20261019/render_pilot.py \
  --pilot-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/coretail_pilot_NEW \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_purity50_20261015.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/coretail_pilot_views_NEW
```

001の固定profile pilotは48実行64fit・8点のGaussian厳密再現検査まで完了。
対応するviews001の12冊64頁は全頁readback済み。
固定tailでχ²の改善はあるがd/t等では不十分であり、次にtail自体もfitする独立pilotを比較する。

## 方針更新：自由core＋tail実データpilotはsuperseded

上記の「次にtail自体もfitする」は当初の検討方針です。
2026-10-19の追加指示を受け、自由core＋tailの実データpilotはprepareまでで停止し、
以後の解析方針としては[原PID関数固定・peak-anchored study](../peak_anchored_20261019/README.md)へ移行しました。
準備済みの24実行32fitを完了済みの実データ結果とは扱わず、自動的に再開しません。
既存の固定profile結果・数学試験・設定・source・過去PDFは履歴として保持します。

理由は、ROIから離れたtの裾を局所窓内で正規化して振幅を自由にすると、
実際のtピークと整合しない巨大な全線yieldでも背景を説明し得るためです。
単に自由成分を増やすのではなく、t・³He等の背景核種の中心・幅を元PID多項式へ固定し、
実際のtピークと³Heの左肩を含む連続領域で振幅を同時にfitします。
核種の関数更新が必要な場合は、別途ピーク校正と滑らかな係数の検証を行います。
新方式も現段階では応答・残差の診断であり、purityや新しい本番cutの認定ではありません。
