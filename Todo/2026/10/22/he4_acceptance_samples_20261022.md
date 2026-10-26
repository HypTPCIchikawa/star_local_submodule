# ⁠⁴HeアクセプタンスのR・TOF別TH2

2026-10-22（JST）。ユーザー依頼により、⁴Heのβ–η、β–Y、pT–YにR・TOF別のTH2Dを追加した。低RのTOFありは「有効なTOFあり／TOF＋m²選択の両方を保存」との回答に従い分けた。

## 保存する4区分

Rは元トラックの|p/q| [GeV/c]。境界は既存nuclearid設定`highRigidityTof.thresholdPOverQ: 2.0`を使用する。

| 名前の末尾 | R範囲 | 条件 |
|---|---|---|
| `_lowR_noTof` | R < 2 | TOF不要。TOFがある候補も含む低Rの基準試料 |
| `_lowR_tofMatch` | R < 2 | 有効TOFあり。m²窓は要求しない |
| `_lowR_tofM2` | R < 2 | 有効TOF＋⁴Heのm²窓 |
| `_highR_tofM2` | R ≥ 2 | 有効TOF＋⁴Heのm²窓 |

有効TOFは既存解析と同じ、有効なTOF traits参照・有限かつ正のβ_TOFから有限m²を計算できる条件。m²窓は既存設定の`3.268773 < (m/q)² < 3.693147`（dと共通、両端除外）。区分判定は元トラックから得たdouble精度のRとm²で行い、候補へのFloat保存後の丸め値は使わない。

これらはすべて最終選別済み⁴Heの部分試料で、TPC best-species、R>0.61、R<5、品質条件は共通。低RにTOFを追加した比較であり、解析本体の候補選別やCF用候補storeは変更しない。

## ヒストグラム名と軸

以下の3つの基底名に、上表の4種類の末尾を付けた計12個をROOTトップ階層に保存する。

| 基底名 | 横軸 | 縦軸 |
|---|---|---|
| `hNucleus_BetaVsEta_he4` | η_lab | β=|p|/E |
| `hNucleus_BetaVsYLab_he4` | Y_lab | β=|p|/E |
| `hNucleus_PtVsYLab_he4` | Y_lab | pT [GeV/c] |

例：`hNucleus_PtVsYLab_he4_lowR_noTof`、`hNucleus_PtVsYLab_he4_lowR_tofMatch`、`hNucleus_PtVsYLab_he4_lowR_tofM2`、`hNucleus_PtVsYLab_he4_highR_tofM2`。

各軸のbin・範囲は対応する既存全体分布と同じ。βは候補4元運動量から求め、TOF測定βではない。電荷補正後の⁴Heの物理運動量P=2Rを使う。実験室系のη/YでCMシフトなし。1候補につき1回、ペア数による重みなし。

低Rでは`tofM2 ⊆ tofMatch ⊆ noTof`で重複するので、3試料を足し合わせない。TOF要求による受入れの比較には対応する分布を同じbinで比べられる。各変数で`lowR_noTof + highR_tofM2 = 既存全体分布`が成立する。Λのβ分布は既存のsignal候補ごとの全体分布を保持し、今回のR/TOF別分割は⁴Heを対象とする。

## 適用先・保存情報

- [対象mainconf](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml)だけの`femtoHist`参照を[専用YAML](../../../../config/hist/hist_auau3p85_femtoLambda_he4_he3_control_acceptanceSamples20261022_kf.yaml)へ変更。
- 既存133個の定義を保持して12個を追加し、KFヒストグラム群は145個。旧共有YAML・前回β用YAMLは保持。
- [FemtoLambdaLegacy.cxx](../../../../StMaker/common/FemtoLambdaLegacy.cxx)の候補収集処理で、最終⁴He候補だけを区分して充填。既存の全体アクセプタンス充填はそのまま保持。
- 12個の一部欠落・TH2でない設定・既存全体分布との軸不一致、または明示的なhighRigidityTof設定がない場合はInit失敗。区分分布を持たない旧プロファイルは従来動作。
- ROOTの`FemtoHe4AcceptanceSamples`（TNamed）へ試料の定義、R境界、m²窓、重複・全体一致の規則を保存。

## 検証・再現

STAR SL24y / sl73_x8664_gcc485 / ROOT5で再build成功。[専用テスト](../../../../tests/test_he4_acceptance_samples.C)734,505チェック通過。PicoDst fixtureで、低R TOFなし・窓内・窓外・無効β、高R TOFなし・窓内・窓外、R=2および直前・直後、R下限・上限を評価した。

全12 TH2について、各トラックから独立に計算した期待bin内容、entries、全体分布とのbin単位の一致、低Rの包含関係をunderflow/overflowまで照合。ROOT読み戻しで12個のTH2と区分メタデータを確認。既存βアクセプタンスの53チェックも再実行して通過。

```bash
./script/singularity_make.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml --no-clean -j2
bash script/singularity_root_batch.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml root4star -b -q 'tests/run_he4_acceptance_samples_test.C("config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml","tmp/he4_acceptance_samples_recheck")'
```

テスト出力先は未使用ディレクトリを指定する。監査先：`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/he4_acceptance_samples_20261022_001/`。`before/`、`build.log`、`test_002.log`、`runtime/he4_acceptance_samples.root`、`beta_regression.log`、`completed.json`を保存。初回`test.log`はROOT5のGetNcells非対応によるテスト側コンパイル診断で、修正後の最終結果は`test_002.log`。

今後新たに解析するROOTへ追加される。既存実データROOTは変更せず、今回実データ新規event処理・farm投入・hadd・commit/pushは行っていない。
