# Task 2 実装結果：高剛性⁴HeのTOF選別

[計画](task2_he4_tof_pid_plan.md)に従い、Λ–⁴Heの核選別へTPC＋TOFを追加した。
新profileは`_tofExtended`、元profileもそのまま使用できる。

| 条件 | 元profile | 新⁴He profile |
|---|---|---|
| TPC | 較正nσの絶対値 < 2、d/t/³He/⁴Heのbest species | 同じ |
| R < 2.0 GeV/c | TPCのみ | 同じ |
| R ≥ 2.0 GeV/c | R < 2.5はTPCのみ、2.5以上は除外 | 有効TOF＋3.0 < (m/q)² < 4.0を要求 |
| ⁴He上限 | R < 2.5 | capを外す。TPC較正0.4～5.0の範囲内 |
| d/t/³He上限 | R < 2.5 | 同じ |

R=2.0は高剛性側、TOF窓の両端は除外。R=|p/q|、TOF量は`R²(1/β²−1)=(m/q)²`。
⁴Heの中心は`(3.72742/2)²=3.4734149641`。Z=2補正は候補四元運動量で一度だけ掛ける。
初期窓3.0～4.0は変更可能で、TOF分解能fitによる最適化値ではない。
TPCで他核種が近い場合、TOFを理由に⁴Heへ繰り上げない。
未知/欠損YAMLキー、不正窓、非有限値は検査し、不正設定でInitを失敗させる。

## 設定・変更箇所

- [3.85 GeV新mainconf](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_tofExtended.yaml)
- [13.5 GeV新mainconf](../../../../config/mainconf/main_auau13p5_anaFemtoLambda_4He_KFParticle_tofExtended.yaml)
- [核PID YAML](../../../../config/cuts/nuclearid/nuclearid_femtoLambda_he4_tofExtended.yaml)
- [新核QA YAML](../../../../config/hist/hist_femtoLambda_he4_nuclear_tofExtended.yaml)
- [選別コード](../../../../StMaker/common/FemtoLambdaLegacy.cxx)
- [ヒスト一覧・意味](../../../../mdfiles/femto/implementation/femto_lambda_histograms_20260930.md)

従来へ戻す場合は、元の`main_auau{3p85,13p5}_anaFemtoLambda_4He_KFParticle.yaml`を指定する。
旧核PID・旧ヒストYAML・ΛのKF条件・±3σ窓は保持。3.85 GeVは2021 P24iyのみ。
新6ヒストは`true/`の`hDedxP_4He_TOFControl`、`hNSigmaHe4VsP_TOFControl`、
`hM2P_4He_TPC`、`hDedxP_4He_selected`、`hM2P_4He_selected`、`hHe4PidStage`。
TOFControlはTOF後・TPC PID前なのでピークがTPC窓で切断されず、次回の純度fitに使える。
旧ROOTから求めるTask 3のTPC純度を、新TOF後の純度として扱わない。

## 検証

STAR SL24y／ROOT 5.34/38／GCC 4.8.5のSingularityでmake all相当がexit 0。
新mainconfで`bash script/singularity_make.sh MAINCONF --no-clean -j4`を実行した。
[ビルドログ](validation_build.log)。
[ROOT5テスト](../../../../tests/test_he4_highrigidity_tof.C)は両エネルギー各845 checks PASS。
[3.85 GeVログ](validation_he4_tof.log)、[13.5 GeVログ](validation_he4_tof_13p5.log)。
旧PID組合せ回帰、境界・非有限値・³He競合・他種上限、11不正設定、閾値変更、
14人工Pico trackのβ変換・Z=2・6QA収量を確認した。
人工入力で旧⁴He3候補→新4候補、stage counts=[8,1,7,4,3,4]と期待値に一致。

2021 P24iyの検証済みローカルPicoDst
`st_physics_22123040_raw_5000023.picoDst.root`（4027 Entries）から先頭1000入力を解析。
completed=1000、reconstructed=965、selectedLambda=16、status=0、エラー/致命エラー0、exit 0。
[実イベントログ](validation_he4_tof_events.log)。
出力は[local1000.root](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/he4_tofExtended_20261013/local1000.root)。
ROOTの新6ヒストを読み取り検査し、最終⁴He84候補、selected-dE/dx84、acceptance84が一致した。
stage countsは[92,84,8,0,0,84]。この1000イベントでは高R TPC候補8件に有効TOFがなく、
高R採用0件なので、高Rでの実データ収量増加までは検証できていない。
[出力監査JSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/he4_tofExtended_20261013/output_audit.json)。
旧エントリの静的テスト8件、旧ヒスト契約632比較も成功した。
人工候補の差は実データの収量・純度改善率ではない。本番ジョブは投入していない。
