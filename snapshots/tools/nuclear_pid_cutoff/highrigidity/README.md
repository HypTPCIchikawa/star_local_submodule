# 高rigidityのTOF条件付き核PID purity（2026-10-17）

低Rの既存80実行・source・数学・集約結果を変更せず、新しい実行部だけを分離する。
`FitHighRigidity.C`は既存local-fitの数学・TPC gate積分・帰属診断・ROOT読み戻しを再利用する。

## 入力と条件

- 2021 P24iy部分統合ROOT由来の不変nativeキャッシュを使用する。
- d／t／³Heは保存済みpre-TPCの`TAG_D`／`TAG_T`／`TAG_HE3`を直接fitし、保存m²窓と現行設定の一致を検証する。
- ⁴Heは今回のユーザー指定によりdと同じ`TAG_D`窓で評価する。現行本番設定の`3<m²<4`とは異なり、本番窓のpurityとは称さない。
- 4核種の元TPC winnerを先に決め、targetのm²窓・R範囲はその後に適用する。TOF不合格の競合核種を除いて再割当しない。
- d／t／³Heは`2<=R<2.5`、⁴He studyは`2<=R<=5 GeV/c`のprovider有効範囲を調べる。
- 2.5超には既存TAG応答fitがない。providerは初期値としてのみ使用し、対象は選択スペクトル自身で再fitする。固定nuisanceは未検証と記録し、寄与が無視できない場合は未評価とする。

同じデータに由来する旧TAG fitは初期値／帰属参照だけに用い、独立priorや追加likelihoodにしない。
対象中心と幅はfree、選択TPC gateは現行関数に固定する。Gaussianの平均／標準偏差とlog-normalのmode／局所曲率幅を区別する。

## 実行と出力

`settings.yaml`にruntime mainconf、native R、m²対応、数学・支持条件を明示する。
`run_highrigidity.py`は単独実行、`run_suite.py`は明示target／モデル／感度variationの最大2並列実行。
`--stage pilot`はnative3点、`upper-pilot`は⁴Heの高R3点、`native`は宣言範囲の全native R binをfitする。
`--self-test-only`で同じsnapshotのROOT5 C++試験を実行できる。

モデルROOTの`raw_native_%03d`、`prediction_total|smooth|d|t|he3|he4|fit_mask_%03d`は既存描画と互換。
`native_summary.tsv`に数値purityと支持flag、fitした中心・幅・誤差を保存する。
`summarize_highrigidity.py`は`sum(S)/sum(total)`を使って0.05 GeV/c区間へ集約し、
一つでもnative fit／モデルが未支持なら区間全体を未評価とする。失敗点を補間しない。

```bash
python3 -B -m unittest discover -s tools/nuclear_pid_cutoff/highrigidity -p 'test_*.py'
python3 -B tools/nuclear_pid_cutoff/highrigidity/run_suite.py \
  --settings tools/nuclear_pid_cutoff/highrigidity/settings.yaml \
  --stage pilot --targets d t he3 he4 --variations nominal --suffix NEW --workers 2
```

数値結果はinclusive track purityであり、イベント数・真の信号効率・CF pair条件付きpurityではない。
TOFタグをpurity=1と仮定しない。高統計χ²だけで全点を無効とせず、対象共分散・誤差・収量・帰属・gate残差・未支持nuisanceを確認する。
model envelopeは信頼区間・完全な系統誤差ではない。既存出力を上書きせず、source／設定SHAと元ROOT fingerprintを前後検証する。

## 最終結果とcoarse独立再fit

出力baseは`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/highrigidity_tof_20261017_001/`。

| 成果物 | 内容 |
|---|---|
| `summary_001/` | 全1800 native fitの診断。既存結果として保持。 |
| `coarse_inputs_002/` | 元.01-Rキャッシュの隣接5binsを厳密加算。全cell・保存Sumw2・flow・閉じたROOTの読戻しを検証。001はROOT5 API compile失敗ログのみ。 |
| `coarse_summary_001/` | He4全2–5を.05-R幅で独立に再fit。240fit、良いnative点の摘み取りやpurity平均ではない。 |
| `summary_final_001/` | d/t/He3は5native fitのΣS/Σtotal、He4は全域coarse独立fitを使う最終表。各手法を明記。 |

高RはGaussian／log-normalとlinear／exponentialの**4 nominal条件**のenvelope。
低Rは4モデル×5 nuisance感度の**20条件**であり、帯の条件数と意味が同一ではない。
He4のcoarse幅は.05-R投影分布自身の幅で、intrinsic幅や新しい本番PID係数として流用しない。

高Rの支持区間はd/t/He3が`[2,2.5)`、He4が`[2,3.2)`。
各区間中の最小model-low purityは順に99.43%、99.51%、98.68%、98.52%。
He4の`[3.2,5)`は未評価であって、purityが50%未満と確定した領域ではない。
特に3.95–4.00では.05化でabsolute残差が約7%まで改善しても未検証nuisance寄与が残り、支持済みにはしない。
4.95–5.00では低統計・残差・対象幅の帰属screenも残る。判定基準は緩和していない。

He4高Rのraw候補は約356万、支持域は約338万（94.79%）。これはinclusive候補数であり真の信号効率ではない。
詳細な未支持理由・nuisance割合は`summary_final_001/summary.json`、中心／幅／誤差は`native_parameters.tsv`。

native全域の再実行例（未使用suffixを指定）：

```bash
python3 -B tools/nuclear_pid_cutoff/highrigidity/run_suite.py \
  --settings tools/nuclear_pid_cutoff/highrigidity/settings.yaml \
  --stage native --targets d t he3 he4 --variations nominal --suffix NEW --workers 2
python3 -B tools/nuclear_pid_cutoff/highrigidity/run_coarse_suite.py \
  --settings tools/nuclear_pid_cutoff/highrigidity/coarse_settings.yaml \
  --stage native --variations nominal --suffix NEW --workers 2
```

`summarize_highrigidity.py`／`summarize_coarse.py`の`--help`で明示入力と新規出力先を確認する。
前者はnative診断表、後者は独立coarse表と最終混合表を別directoryへ作成する。
元cache・fit結果・本番設定を自動変更／上書きする手順ではない。
