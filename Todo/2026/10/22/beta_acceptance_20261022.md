# Λ・⁴Heのβ–η／β–Yアクセプタンス追加

2026-10-22（JST）。ユーザー依頼により、anaFemtoLambda_4He_KFParticleへ4個のTH2Dを追加した。

| ヒストグラム | 横軸 | 縦軸 | 充填対象 |
|---|---|---|---|
| `hLambda_BetaVsEta_signal` | η_lab | β | Λ信号窓内の候補 |
| `hLambda_BetaVsYLab_signal` | Y_lab | β | 同上 |
| `hNucleus_BetaVsEta_he4` | η_lab | β | 最終PID選択後の⁴He候補 |
| `hNucleus_BetaVsYLab_he4` | Y_lab | β | 同上 |

β=|p|/Eは候補4元運動量から求める運動学的な速度（c単位）で、TOF測定βではない。Λは既存アクセプタンスと同じ候補4元運動量、⁴Heは電荷数2の補正済み物理運動量と⁴He質量を用いる。Yとηは実験室系で、CMシフトなし。η/Yの横軸は200 bin、−5～5、βの縦軸は200 bin、0～1。軸設定はYAMLに置く。

Λは既存signalアクセプタンスと同じ候補、⁴Heは既存の最終選別済み核候補を1候補につき1回充填する。Λ sidebandは含めず、SE/MEペア数による重みもかからない。βまたは各横軸の値が非有限の場合は該当する追加分布への充填を行わず、既存の候補選択は変更しない。

## 実装

- [対象mainconf](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml)の`femtoHist`を[専用ヒストグラムYAML](../../../../config/hist/hist_auau3p85_femtoLambda_he4_he3_control_betaAcceptance20261022_kf.yaml)へ変更。
- 旧共有YAMLを保持し、既存129個の定義を同じ内容で引き継ぎ、4個を追加した。対象以外のmainconfは変更していない。
- [FemtoLambdaLegacy.cxx](../../../../StMaker/common/FemtoLambdaLegacy.cxx)の単粒子アクセプタンス充填処理に追加。4個の一部だけを設定した場合やTH2でない場合はInit失敗。追加分布がない旧プロファイルは従来通り動作する。
- ROOTのトップ階層へ保存。更新後に新たに解析して作るROOTへ収録され、既存ROOTへの追記・書換えは行っていない。

## 検証

STAR SL24y / sl73_x8664_gcc485 / ROOT5のコンテナ環境で再build成功。

[専用テスト](../../../../tests/test_femto_beta_acceptance.C)の53チェックが通過した。4個のTH2・軸、Λの正負方向の運動量、β=|p|/Eの該当bin、⁴HeのP=2R、R下限・高R TOFの選別済み候補だけの充填、Λ sideband除外、ペア充填を繰り返しても単粒子entriesが増えないこと、ROOT読み戻しを検証した。

既存PID採用テストも新mainconfで再実行し14,147チェック通過。d/t/³HeのGit係数、⁴He確定係数、TOF境界、³He背景controlは維持される。

監査出力：`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/beta_acceptance_20261022_001/`。`before/`、`build.log`、`test_002.log`、`runtime/beta_acceptance.root`、`pid_regression.log`、`completed.json`を保存。最初の`test.log`はテスト側TH2ヘッダ不足によるコンパイル診断で、修正後の結果は`test_002.log`。

```bash
./script/singularity_make.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml --no-clean -j2
bash script/singularity_root_batch.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml root4star -b -q 'tests/run_femto_beta_acceptance_test.C("config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml","unused","unused","tmp/beta_acceptance_recheck")'
```

テスト出力先には未使用ディレクトリを指定する。中間の2引数は共通形式のテストrunner用の未使用引数。実データの新規event処理・farm投入・hadd・commit/pushは行っていない。
