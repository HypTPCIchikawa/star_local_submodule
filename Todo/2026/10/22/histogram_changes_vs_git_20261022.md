# anaFemtoLambda_4He_KFParticle：Git版からのヒストグラム追加・変更

作成日：2026-10-22（JST）。比較基準：Git HEAD `8f7d539`。

現在の入口は `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml`。この名前のmainconfはGit HEADには存在しないため、Git版の `main_auau3p85_anaFemtoLambda_4He_KFParticle_highpurity.yaml` が参照する⁴He KF解析のヒストグラム定義と比較した。以下は現在の設定・実装の比較であり、過去に作成したROOTファイルへの収録を保証するものではない。

## 追加・変更の概要

| 追加・変更 | 主なヒストグラム名 | 簡単な説明 |
|---|---|---|
| Λ–³Heのk*分布 | `hKstarSE_lambda_he3`、`hKstarME_lambda_he3` | 選別した³Heを物理的な³He質量で扱う相関分布。中心度別・Λ信号窓・左右sideband別も用意。 |
| ³Heを⁴He質量で扱う背景分布 | `hKstarSE_lambda_he3_as_he4`、`hKstarME_lambda_he3_as_he4` | 同じ³He候補の運動量を使い、質量を⁴Heとしてk*を再計算。³He混入が⁴He相関に与える影響を調べる。 |
| ³He関連の詳細ペア分布 | `hKstarMass_3He*_CentBinN`、`hKstarMt_3He*_CentBinN`、`hQlab_3He*` | k*とΛ不変質量、ペアmT、実験室系相対運動量の分布。`3HeAs4He`版とMixed-event版も追加。 |
| ³Heのアクセプタンス・近接ペアQA | `hNucleus_PtVsYLab_he3`、`hNucleus_PtVsYLab_he3_as_he4`、`hDphiDeta_proton_3He*` | 選別した核のpT–rapidity分布と、Λ娘陽子との近接ペアを確認。 |
| ⁴HeのPID選択前後のQA | `hDedxP_4He_TOFControl`、`hNSigmaHe4VsP_TOFControl`、`hDedxP_4He_selected` | TPC PID前の対照試料と最終選択後の分布を比較。対照試料もR≥2ではTOF条件を要求。 |
| ⁴HeのTOF・選別段階QA | `hM2P_4He_TPC`、`hM2P_4He_selected`、`hHe4PidStage` | TOF選択前後のm²–R分布と、TPC判定から最終選択までの段階別カウント。 |
| k*に依存したPID QA | `hKstarNSigma_*`、`hKstarPOverQ_*`と各`PrePid`版 | k*ごとのnσ・R分布。⁴He、³He、³HeAs4Heについて、選択後／TPC PID前、SE／ME、中心度・Λ質量領域別に保存。 |
| PIDの多次元分布 | `hNuclearPidPairJoint`（THnSparseD） | k*・R・nσ・中心度・核種表現・SE/ME・Λ質量領域・PID段階を同時に保存。Rとnσを組み合わせて後から投影できる。 |
| Mixing QAの拡張 | `hMixSamplerQA` | ³Heチャネル追加に合わせ、チャネル軸を4から8へ拡張。 |

## 名前・選択条件の読み方

- `SE`：同一イベント、`ME`または`Mixed`：混合イベント。
- `CentBinN`：cent9の0～8。大きい番号ほど中心衝突側。
- `signal`：Λ信号窓、`leftSB`／`SBNeg`：低質量側sideband、`rightSB`／`SBPos`：高質量側sideband。
- `signalLow`／`signalHigh`：既存のΛ信号窓を中心質量で上下に分割。Lowは中心未満、Highは中心以上。
- `3HeAs4He`／`he3_as_he4`：³He候補の運動量を保ち、⁴He質量を仮定した表現。別の独立した核候補やイベント試料ではない。nσは元の³He仮説を使用する。
- `PrePid`：トラック品質、核種別R範囲・veto、TOF条件を適用し、TPCのnσ選択・best-species判定を課す前のペアQA。
- 現在の最終候補は全核種でR<2 GeV/cはTOF不要、R≥2 GeV/cは核種別TOF窓を要求。³Heの物理質量版と⁴He質量仮定版は同じ選別済み候補を使う。既存のd/t/³Heの最小R・vetoも適用される。
- ペアPID QAのRは保存された核4元運動量からP/2で求めた値。⁴Heと³Heの電荷数2を考慮している。

⁴Heのk*–mT分布、単粒子アクセプタンス、Λ信号窓の`signalLow/High`分割は既にGit版に含まれる。今回の新規追加として数えているのは、主にその³He／³HeAs4Heへの展開とPID QAである。

## 定義数の比較

| 定義群 | Git版 | 現在 | 追加 | 既存定義の変更 | 削除 |
|---|---:|---:|---:|---:|---:|
| KF・相関・アクセプタンス | 79 | 129 | 50 | 1 | 0 |
| NuclearId・詳細ペア分布 | 203 | 561 | 358 | 0 | 0 |
| ペアPID QA（TH2） | 0 | 1,080 | 1,080 | 0 | 0 |
| 合計（通常ヒストグラム） | 282 | 1,770 | 1,488 | 1 | 0 |

これに8次元THnSparseDを1個追加。ペアPID QAの1,080個は、2変数（nσ／R）×3表現×2（SE／ME）×5質量区分（信号全体・上下半分・左右sideband）×9中心度×2 PID段階に対応する。信号全体と上下半分は重複を含む表示であり、独立試料として足し合わせない。

変更1個は`hMixSamplerQA`のチャネル軸（4→8）。名称・軸定義が同じ既存ヒストグラムも、今回採用したPIDとTOF条件により入力候補や内容は変わる。

## 定義・実装ファイル

この文書からの相対パスはリポジトリ内の対応ファイルを指す。

- [現在のmainconf](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml)
- [相関・アクセプタンス定義](../../../../config/hist/hist_auau3p85_femtoLambda_he4_he3_control_kf.yaml)
- [核PID・詳細ペア分布定義](../../../../config/hist/hist_auau3p85_femtoLambda_he4_he3_control_nuclear.yaml)
- [ペアPID QA・Sparse定義](../../../../config/hist/hist_auau3p85_femtoLambda_he4_he3_pid_pair_qa.yaml)
- [選別・ヒストグラム充填処理](../../../../StMaker/common/FemtoLambdaLegacy.cxx)
- [Sparse実装](../../../../StMaker/common/FemtoNuclearPairPidSparse.cxx)
- [確定PID採用と検証記録](pid_adoption_20261022.md)

比較元はGit `8f7d539`の`config/hist/hist_femtoLambda_he4_kf.yaml`と`hist_femtoLambda_he4_nuclear.yaml`。現在の同名作業ファイルではなく、`git show`で取得した内容を比較に使用した。

ROOT出力では、KF・相関ヒストグラムとSparseはトップ階層、NuclearIdおよびペアPID QAの同一イベント側は`true/`、混合イベント側は`mix/`へ保存される。

## 2026-10-22追記：Λ・⁴Heのβアクセプタンス

後続依頼により`hLambda_BetaVsEta_signal`、`hLambda_BetaVsYLab_signal`、`hNucleus_BetaVsEta_he4`、`hNucleus_BetaVsYLab_he4`の4個のTH2Dを追加した。横軸は実験室系η/Y、縦軸はβ=|p|/E。Λ信号候補または最終選別⁴He候補につき1回充填する。TOF測定βではない。

現在の対象mainconfは`hist_auau3p85_femtoLambda_he4_he3_control_betaAcceptance20261022_kf.yaml`を参照。本文の旧共有YAMLは変更前の定義として保持した。KF群は133個、通常ヒストグラム合計は1,774個、Git比較の追加は1,492個となる。既存変更1個・削除0・Sparse追加1個は同じ。[実装・検証詳細](beta_acceptance_20261022.md)。

## 2026-10-22追記：⁴HeのR・TOF別アクセプタンス

β–η、β–Y、pT–Yの各分布に、低R TOF不要・低R有効TOF・低R TOF＋m²・高R TOF＋m²の4区分（計12 TH2D）を追加。低RはR<2、高RはR≥2。対応する全体分布を保持し、低R3区分の重複を明示した。現在のmainconfは`hist_auau3p85_femtoLambda_he4_he3_control_acceptanceSamples20261022_kf.yaml`を参照。KF群145個、通常ヒストグラム合計1,786個、Git比較の追加1,504個。既存変更1個・削除0・Sparse追加1個は同じ。[区分名・条件・検証](he4_acceptance_samples_20261022.md)。
