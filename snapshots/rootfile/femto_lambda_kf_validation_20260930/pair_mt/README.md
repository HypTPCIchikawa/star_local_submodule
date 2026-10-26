# centrality別k*–m_Tの追加検証（2026-09-30）

既存のKF Λ–原子核Femto解析へ、SE／ME・cent9=0〜8・5質量領域ごとの90 TH2Dを追加した出力です。
高purityカット、±3σ、sideband、原子核選別、mixing、既存192ヒストは変更しません。
比較元は `../acceptance_split/{13p5,3p85}/{d,t,3He,4He}.root` で、元ROOT・PDFは保持しています。

## m_Tの定義

ユーザー確認済みの式（c=1）：

```text
m_T = sqrt((|pT_vector_Lambda + pT_vector_nucleus|/2)^2
           + ((m_Lambda + m_nucleus)/2)^2)
```

ペア全体の横質量ではありません。質量は既存YAMLの固定質量、横運動量はペア用の物理運動量です。
HeのZ補正は候補作成時の1回だけです。式と使用質量はROOTの `FemtoLambdaPairMtDefinition` に記録します。

名前・軸・射影例は[ヒストグラム一覧](../../../mdfiles/femto/implementation/femto_lambda_histograms_20260930.md)を参照してください。
X=k*（200 bin、0〜1 GeV/c）、Y=m_T（200 bin、0〜10 GeV/c²、幅0.05）。
centralityをSE／ME別に加算した後、同じm_T bin範囲をXへ射影し、その範囲内でCFを正規化します。

## 実行・保存先

すべてローカル実行で、各条件とも前回と同じ先頭10,000入力イベントです。farm投入はありません。

- mainconf：`config/mainconf/main_auau<energy>_anaFemtoLambda_<species>_KFParticle_highpurity.yaml`
- 入力：`config/picoDstList/auau<energy>_femtoLambda_validation_20260930.list`
- `<energy>`：`13p5,3p85`、`<species>`：`d,t,3He,4He`
- job ID：`pair_mt_<energy>_<species>`
- 出力：`<energy>/<species>.root`（1ファイル282ヒスト＋Tree／metadata）
- 実行ログ：`provenance/<energy>_<species>.log`
- 監査ログ：`provenance/audit_<energy>_<species>.log`
- QAログ：`provenance/qa_<energy>_<species>.log`
- ビルドログ：`provenance/build.log`（SL24y／ROOT5、make all）
- SUMS job ID／configlog：なし（ローカル実行）。KF実効設定はROOTへ保存。

入力リスト・ROOT・PDF・実行ログはローカルファイルで、新規cloneには含まれません。
再実行時は入力を用意し、既存ROOTを上書きしない新しい出力名を指定してください。

## QA PDF

従来7ページを保ち、8〜12ページ目にsignal／SBPos／SBNeg／signalLow／signalHighのk*–m_Tを追加しています。
各ページの左がSE、右がMEです。図はmakerの `cfCent9Min..cfCent9Max` を合算しますが、ROOTはcent9ごとに保存しています。

| データ | d | t | ³He | ⁴He |
|---|---|---|---|---|
| 13p5 | [QA](../../../share/figure/femto_lambda_kf_validation_20260930/pair_mt/new_13p5_d.pdf) | [QA](../../../share/figure/femto_lambda_kf_validation_20260930/pair_mt/new_13p5_t.pdf) | [QA](../../../share/figure/femto_lambda_kf_validation_20260930/pair_mt/new_13p5_3He.pdf) | [QA](../../../share/figure/femto_lambda_kf_validation_20260930/pair_mt/new_13p5_4He.pdf) |
| 3p85 | [QA](../../../share/figure/femto_lambda_kf_validation_20260930/pair_mt/new_3p85_d.pdf) | [QA](../../../share/figure/femto_lambda_kf_validation_20260930/pair_mt/new_3p85_t.pdf) | [QA](../../../share/figure/femto_lambda_kf_validation_20260930/pair_mt/new_3p85_3He.pdf) | [QA](../../../share/figure/femto_lambda_kf_validation_20260930/pair_mt/new_3p85_4He.pdf) |

## 検証方法

- `tests/test_femto_lambda_hist_contract.py`：旧632契約項目を保ち、新90キー×4核種、軸、型、282ヒストを検査。
- `tests/femto_lambda_legacy.C`：平均質量・横運動量ベクトル和の数値、縦運動量に依存しないこと、He補正、全cent9／質量領域／SE・ME、flow、誤差を人工データで確認。
- `tests/check_femto_lambda_mt.C`：旧192ヒストが全bin・誤差・entriesまで不変、全イベント／採用Λ記録が一致、90 TH2のflow込みProjectionXが既存k*と一致、18組の低側＋高側＝signalを検査。
- QA：旧ROOTは7ページを維持、新ROOTは12ページ。部分欠落・異次元・軸不一致はPDFを作らず失敗する。

実イベントROOTには原子核ペアごとの独立ledgerがないため、m_Tの数値そのものは人工データ試験で検証します。
生成したm_T別CFの物理的な解釈・fitや、低統計領域のpurity評価は今回行っていません。

## 検証結果

2026-09-30にすべて完了しました。

| 検証 | 結果 |
|---|---|
| SL24y／ROOT5.34/38／GCC4.8.5のmake all | 成功 |
| 4核種の人工データ試験 | 全て成功。初回d試験の明示的Sumw2配列要求はテスト側の過剰条件として修正し、bin誤差・射影誤差の検証は維持 |
| 13p5／3p85 × 4核種の実イベント | 全8実行、各10,000イベント、status=0 |
| 旧192ヒスト・イベント／候補記録 | 全8条件で完全一致 |
| 新TH2の全m_T射影と誤差 | 90個×8条件で既存k*と一致（flow込み） |
| signalLow＋signalHigh＝signal | SE／ME×9cent＝18組×8条件、全二次元bin・誤差二乗和が一致 |
| QA PDF | 全8ファイルで12ページ。13p5 Λ–dの追加ページを目視確認 |
| 旧ROOT／異常ROOT | 旧4／7ページを維持。部分欠落・異次元・軸不一致はexit1、PDF未生成 |
| 保護対象コード・設定 | 289ファイル不変 |
| 静的テスト | 入口6、旧ヒスト632項目＋新90キー／核種が成功 |

13p5では入力10000、再構成9744、採用Λ873、3p85では入力10000、再構成9505、採用Λ144で、追加前と一致しました。
人工データのログは `provenance/toy_retry_d.log`、`toy_t.log`、`toy_3He.log`、`toy_4He.log`。
