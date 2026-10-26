# pT–Y_Lab・信号質量窓二分割の検証（2026-09-30）

既存のKF Λ–原子核解析にヒストグラムだけを追加した再解析です。
各エネルギーについて、前回と同じ入力リストの先頭10,000イベントを全4核種で処理しました。
高purityカット・±3σの信号窓・原子核選別・mixingは変更していません。
比較元の `../{13p5,3p85}/{d,t,3He,4He}.root` は保持しています。

## 保存内容

- `{13p5,3p85}/{d,t,3He,4He}.root`：新ヒスト12個を追加した出力。1ファイル192ヒスト。
- `provenance/<energy>_<species>.log`：10,000入力イベント処理と正常終了のログ。
- `provenance/audit_<energy>_<species>.log`：既存180ヒスト・eventLedger・selectedLambdaの新旧一致検証。
- `provenance/qa_new_<energy>_<species>.log`：新QAの実行ログ。
- `provenance/build.log`、`toy_retry_*.log`：SL24y/ROOT5のビルドと4核種の人工データ検証。
- `provenance/toy_d.log`：最初のテスト用ヘッダーがROOT5の辞書生成で失敗した記録。本番較正は変更せず、テスト実装だけをCINTから隠して再試験しました。

ヒストの定義・運動量源・質量境界は[日本語の一覧](../../../mdfiles/femto/implementation/femto_lambda_histograms_20260930.md)を参照してください。

## 実イベント検証結果

全8条件で、以下が成功しました。

- 既存180ヒストの型・軸・全bin内容・誤差・entriesが完全一致（underflow／overflowも含む）。
- 10,000イベントの処理記録と採用Λの元イベント・娘index・質量・運動量が一致。
- SE／ME、1次元／cent9別の各binで `signalLow + signalHigh = signal`。誤差の二乗も加算して一致。
- ΛのpT–Y_Labをcandidate ledgerから再構成して一致。娘2枚とΛのmapのentriesが±3σ内Λ数に一致。

| ビームエネルギー（FXT） | 核種 | 採用Λ（全質量） | ±3σ内Λ／各娘mapのentries | 最終採用原子核mapのentries | 新QA PDF |
|---|---|---:|---:|---:|---|
| 13.5 GeV | d | 873 | 741 | 20344 | [図](../../../share/figure/femto_lambda_kf_validation_20260930/acceptance_split/new_13p5_d.pdf) |
| 13.5 GeV | t | 873 | 741 | 3419 | [図](../../../share/figure/femto_lambda_kf_validation_20260930/acceptance_split/new_13p5_t.pdf) |
| 13.5 GeV | ³He | 873 | 741 | 1122 | [図](../../../share/figure/femto_lambda_kf_validation_20260930/acceptance_split/new_13p5_3He.pdf) |
| 13.5 GeV | ⁴He | 873 | 741 | 175 | [図](../../../share/figure/femto_lambda_kf_validation_20260930/acceptance_split/new_13p5_4He.pdf) |
| 3.85 GeV | d | 144 | 124 | 36398 | [図](../../../share/figure/femto_lambda_kf_validation_20260930/acceptance_split/new_3p85_d.pdf) |
| 3.85 GeV | t | 144 | 124 | 6504 | [図](../../../share/figure/femto_lambda_kf_validation_20260930/acceptance_split/new_3p85_t.pdf) |
| 3.85 GeV | ³He | 144 | 124 | 3919 | [図](../../../share/figure/femto_lambda_kf_validation_20260930/acceptance_split/new_3p85_3He.pdf) |
| 3.85 GeV | ⁴He | 144 | 124 | 878 | [図](../../../share/figure/femto_lambda_kf_validation_20260930/acceptance_split/new_3p85_4He.pdf) |

これらはpair数や検出効率ではありません。原子核mapは同じイベントのΛの有無によらず最終採用trackを数えます。
旧ROOTには娘運動量・原子核候補の独立ledgerがないため、その運動量源と1track1回の記録は人工Pico fixtureで別途検証しました。
低統計の分割CFの物理解釈を保証するものではありません。

## 再実行例

既存出力への上書きは拒否されるため、再実行では未使用の出力ファイル名を指定してください。

```bash
bash script/singularity_run_anaFemtoLambda_d.sh \
  config/mainconf/main_auau13p5_anaFemtoLambda_d_KFParticle_highpurity.yaml \
  config/picoDstList/auau13p5_femtoLambda_validation_20260930.list \
  rootfile/femto_lambda_kf_validation_20260930/acceptance_split/13p5/d_retry.root \
  acceptance_split_13p5_d_retry 10000
```

実行前に `script/singularity_make.sh <mainconf> --no-clean -j4 all` で再ビルドします。
図は `script/singularity_checkHistAnaFemtoLambda.sh <ROOT> <mainconf> <新しいPDF>` で作成します。
新形式は7ページで、5ページ目がpT–Y_Lab、6ページ目が低質量側CF、7ページ目が高質量側CFです。
旧ROOTに新ヒストは遡及追加されません。旧形式のQAは警告付きで元の4ページを出力します。

## QAの追加確認

全8 PDFは7ページで正常に開けることを確認しました。13.5 GeV・Λ–dのpT–Y_Labページと低質量側CFを画像化して目視確認しています。旧ROOTは警告付きで4ページを生成して終了コード0、新12ヒストの一部欠落または型の不一致を与えた試験は終了コード1でPDFを作らず停止しました。ログは `qa_old_13p5_d.log`、`qa_negative_partial.log`、`qa_negative_wrong_dimension.log` です。
