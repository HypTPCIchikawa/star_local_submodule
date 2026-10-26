# Task 3実装追補: ³He/⁴He校正とCF背景用スペクトル

## 追補: 低剛性校正は再検証中

以下は前回の点ごとの補間校正を導入・動作検証した記録である。
その後の生分布・表示確認で、R=0.410 GeV/cの境界に⁴He中心の約22.7 keV/cmの
人工的な段差を確認した。動作テストの成功は、核種ピークの位置・幅の正しさを保証しない。
現在は、解析コードと同じ関数形の共通係数fitを別studyで検討している。
[再校正・検証結果](task3_helium_logpoly_recalibration_report.md)を優先して参照する。
元の係数・出力・切り戻し設定は比較用に保持している。

## 現在の状態

PID校正・³He出力の実装、STARビルド、同一1,000入力の実イベント試験、ROOT読み戻し監査・人工pair閉包を完了し、最終PDF・来歴・切り戻し設定を保存した。計画は[実装計画](task3_helium_pid_cf_implementation_plan.md)、元fitの採用可否・数値・検証は[0.40 GeV/c開始の報告](task3_from040_refit_report.md)を参照。

対象は3.85 GeV・2021年P24iyのΛ–⁴He解析。13.5 GeV、既存Helix/Phi解析、元の統合ROOT、過去の投入設定は変更しない。現行設定の更新前に従来PID用mainconfを保存する。

## 採用する校正と純度の定義

- 元fitで採用した205個の剛性knotsを保存し、中心値とGaussian幅を区分線形補間する。別背景モデルへの切替えで有効bin数を増やすことは行わない。
- 校正の対象域は0.41–2.50 GeV/c。0.40–0.41では⁴Heの信号位置がfit境界に残るため採用せず、対象域外は明示的に従来校正を使用する。
- 有効な隣接fit bin間は補間し、欠落binをまたぐ広いgapでは、その有効binの半幅0.005 GeV/c以内のみ端点値を使う。未支持の中央区間は従来校正に戻す。
- 中心・幅を変更するが、選択のnσ上限は従来の2.0を維持する。d・tの校正は変えない。
- Purity対p/qは、旧校正±2σ、新校正±2σ、新校正＋best-speciesのTPCモデル推定を分けて再計算する。背景モデル差を併記し、TOF後やペア・centrality別の確定純度とは扱わない。
- 低剛性⁴Heは背景の肩をmodelで分離した結果で、³Heも同じfitの背景成分からの校正である。大部分のbinは参考fitであり、この校正を独立に検証したPIDと主張しない。

### 解析に適用する関数

[校正YAML](../../../../config/calibration/dedx_helium_auau3p85_2021P24iy_from040.yaml)に、各核種205個の`rigidity / mean / sigma`をkeV/cmで保存した。source fitは6品質通過・199参考fitで、各点のmodel・参考区分・境界・誤差などを来歴に保存する。大域多項式候補は今回の低剛性構造を平滑化し過ぎるため導入しない。

隣接する支持点の間では中心と幅をそれぞれ線形補間し、`nσ_s = (dE/dx − μ_s(R))/σ_s(R)`を使う。半bin以内の端点clamp、未支持gapとdomain外のlegacy復帰はYAMLとproviderの共通契約である。独立したPython実装と実際のC++値・校正適用flagを5,370点で照合した。

`NuclearIdDeDxVsMom.h`の全体係数を書き換えるのではなく、`FemtoLambdaLegacy`の所有する新providerへ導入した。従って³He/⁴HeのTPC仮説QA・best-species・最終候補が同じ新pullを使用し、d/tのpullと他エネルギーの従来profileは保持される。明示された校正が不明・不正ならInitは失敗する。実行設定にもauau3p85・year2021・P24iyを要求する。

### 更新したPurity対p/q

- [新旧窓・best-species・背景関数・中心/幅の4ページPDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/calibrated_selection/purity_calibrated_windows.pdf)
- [1ページ目PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/calibrated_selection/purity_calibrated_windows.png)
- [全model・選択窓の数値CSV](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/calibrated_selection/purity_windows.csv)
- [校正/Purityの再現と独立監査](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/calibrated_selection/validation/runtime_final_001/validation.json)

下表は主表示の背景modelで同じ新窓を積分した比較。全てcharge-inclusive・cent積分のTPCモデル推定であり、数値をTOF後やCFのペア純度へ置き換えない。

| p/q中心 [GeV/c] | 元校正±2σ | 新校正±2σ | 新校正＋TPC best-species |
|---:|---:|---:|---:|
| 0.425 | 59.81% | 42.04% | 67.33% |
| 0.505 | 62.61% | 38.78% | 57.01% |
| 2.495 | 53.94% | 43.58% | 58.94% |

中心と幅を更新すると窓の母集団が変わり、単独の±2σ純度は下がる区間もある。実解析はbest-speciesを要求するため、新窓の³Heとの競合も反映した青曲線を別に示す。元曲線は元±2σのみであり、旧解析のbest-species込みと新解析のefficiencyを比較した図ではない。signal数の元窓比もPDF内に示すが、fit積分の比であって再構成効率ではない。

背景2関数の差は同じ固定窓で評価し、covarianceが保存されていない新しい窓の統計誤差を作らない。図の斜線範囲はmodel spreadで、校正されたsystematic errorでもない。0.405と成立しない主fitの点は純度を採用せず、runtimeのlegacy復帰と図の不採用を区別する。PDF・元fitを上書きせず`calibrated_selection/`へ別保存した。

## 追加するCF入力

実質量のΛ–³Heと、同一の³He候補の運動量に⁴He質量を与えた比較用出力を、⁴He出力と同じROOTへ追加する。SE/ME、native cent9、Λ信号・左右sideband・信号の低/高質量側、質量×k*、k*×m_Tを保存する。

³Heを選ぶ窓の候補は、⁴He窓へ漏れ込む³Heのtailと同一ではない。従ってinclusive TPC純度を差し引き係数として自動適用せず、ペア条件や誤識別の移行を別途調べられる入力だけを作る。

新ROOTは690ヒスト（最上位129、`true/`300、`mix/`261）を持つ。新しい仮想核種をcandidate storeへ追加するのではなく、物理的³Heの運動量をZ=2補正した後に⁴He質量でエネルギーを計算し直す。mapped表現のために再びZ=2を掛けたりmixingを実行したりしない。k*とm_Tは再計算されるが、q_labと同一cent/質量区分のpair数は実質量表現と同じ。

ヒストの名前・軸・対象は[ヒストグラム一覧](../../../../mdfiles/femto/implementation/femto_lambda_histograms_20260930.md)に追記した。`FemtoNuclearDedxCalibration`には実際に消費したYAML全文と絶対パスを、`FemtoLambdaHe3ControlDefinition`等には質量仮説・m_Tの定義を保存する。

## 現行設定と切り戻し

| 入口 | 現行mainconf | 従来PID・⁴Heのみへ戻すmainconf |
|---|---|---|
| 通常 | [KFParticle](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml) | [pidLegacy](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_pidLegacy.yaml) |
| highpurity | [highpurity](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_highpurity.yaml) | [highpurity_pidLegacy](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_highpurity_pidLegacy.yaml) |
| 高R TOF追加 | [tofExtended](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_tofExtended.yaml) | [tofExtended_pidLegacy](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_tofExtended_pidLegacy.yaml) |

各々のanalysis名・output stemは同じなので、新旧比較は別output filenameで行う。共有の旧maker/hist/PID YAMLを変更せず、3p85専用YAMLへ参照を変更した。過去のpilot/recovery mainconf、13p5、他核種の旧入口は変更していない。通常設定ではR<2.5のTPC選択、高R TOF版ではTask2と同じR≥2.0でTOFも要求し、2.5より上の⁴Heは従来校正を使う。

## 変更・検証の進捗

| 項目 | 状態 |
|---|---|
| 実装計画・Todo索引・継続更新ルール | 保存済み |
| 校正YAMLと純度図の再計算 | 完了、5,370点のC++/Python照合・選択窓の独立積分もPASS |
| エネルギー限定のPID適用と切り戻し設定 | 完了、3現行mainconfの参照変更と3切り戻し設定を保存 |
| ³He／³He-as-⁴HeのSE/MEと同等ヒストグラム | 完了、690件の定義・名前・元ヒスト維持の静的テストPASS |
| STAR/SL7環境ビルド・人工PID・実イベント検証 | 最終make all成功、新PID8,840項目、旧TOF各energy845項目PASS、3条件各1,000入力の解析成功 |
| ROOT読み戻し監査 | 最終新通常・新TOFの各1,999,865項目PASS、同一Λ・SE/ME・cent/質量/m_T/half閉包とmixing QAのflowを検証 |
| 人工pair閉包・Makerガード | 各表現・SE/ME45pair、全centの独立boost/m_T oracle・count/error・不正pair排除PASS、旧単一核種profileのInit/failure/ownership試験PASS |
| 最終来歴整理・入力QA | source/config snapshot・manifestと2ページPDFを保存 |

### 小さい実イベント試験

入力は2021 P24iyの復旧cache `st_physics_22123040_raw_5000023.picoDst.root`（4027 Entries）の先頭1,000件。解析出力は[専用testsフォルダ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/helium_calibrated_he3Control_20261013/README.md)へ保存した。

| 条件 | 完了入力 | 再構成イベント | 採用Λ | 最終⁴He track | 最終³He track | Λ–³He signal SE / ME |
|---|---:|---:|---:|---:|---:|---:|
| 新校正・通常 | 1,000 | 965 | 16 | 103 | 370 | 14 / 354 |
| 新校正・TOF追加 | 1,000 | 965 | 16 | 101 | 370 | 14 / 354 |
| 従来PID・⁴Heのみ | 1,000 | 965 | 16 | 本比較では集計対象外 | pairを保存しない | 保存しない |

全3条件はexit 0 / status completed。event ledger・Λ候補ledger・Λ QA・未PID dE/dxは新旧で完全一致した。³Heの実質量／⁴He質量仮説で同一pair数、q_lab、pTが一致し、全centの質量×k*・m_Tの射影、signalLow+High、canonical/legacyのaliasも一致した。小統計では物理CFの形状や³He誤識別係数を確定できない。

最終ログは本Todoの`validation_helium_controls_build_v2.log`、`validation_helium_pid_runtime_final_v2.log`、`validation_helium_events_{new_base,new_tof}_v2.log`（旧条件は`validation_helium_events_pid_legacy.log`）、`validation_helium_outputs_{new_base,new_tof}_final_v2.log`、`validation_femto_helium_{pair,policy_guards}_v2.log`。

独立レビューで、channel名と`partB`の核種対応のチェック不足、およびHe3 base channel index=4が旧mixing QAの4 bin軸のoverflowへ入る問題を見つけた。対応を厳密に検証するInitガードを追加し、新profileだけ`hMixSamplerQA`を8 channel軸へ拡張した。旧profileは4 binを維持し、元の物理ヒストの定義は保持する。取り違えたbase/信号channelはInitで拒否されるnegative testを追加した。最終実イベントROOTではHe3 QAのin-range加重和6,791、channel under/overflow=0を確認した。

初回の`new_{base,tof}_1000.root`とログは履歴として保持し、最終利用対象は`new_{base,tof}_1000_v2.root`。初回・最終で物理結果は変わらず、QA軸を修正して再出力した。途中の監査マクロにROOT6向け`GetNcells`が含まれたため最初の監査compileも失敗したが、ROOT5で使える軸bin数からの計算へ変更して再実行した。元統合ROOTやfit成果物には変更していない。

[3つの質量表現のSE/ME・左右SB・m_T入力QA（20 MeV/c、2ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/helium_calibrated_he3Control_20261013/helium_control_inputs.pdf)は小統計の接続確認。最終ROOT・ソース・参照YAML・ライブラリ・ログ・入力cacheのSHAとsource/config snapshotは[validation manifest](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/helium_calibrated_he3Control_20261013/validation_manifest.json)へ保存する。独立Python照合の元8,826項目fixtureと最終8,840項目fixtureの校正TSVはSHA256まで一致し、レビュー修正による校正数値の変更はない。

## 残る作業と運用

- 新設定で必要な統計を再解析してから³He control CFを作る。今回、本投入はしていない。
- ³He core→⁴He窓のtailの移行、TOF・centrality・k*・m_Tに条件づけた混入率・背景形状を検証し、背景CFの係数を決める。今回のinclusive TPC純度をそのまま補正係数へ使わない。
- 現在はdirty worktreeへの保存で、git pushはしていない。通常`config/`等はignore対象なので、将来のcommit時には新しい校正・専用maker/PID/hist・切り戻しmainconfも明示的に含め、参照先が欠けないようにする。
- `docs/ai/AGENT_RULES.md`へ、Todoの進捗・現行設定・出力リンク・検証・未完了点を分析中から更新する継続ルールを追加した。

本投入・git pushは行っていない。今回追加するヒストグラムは既存の統合ROOTへ後付けされず、再解析で保存される。
