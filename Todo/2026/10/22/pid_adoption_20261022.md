# 2026-10-22 確定⁴He PIDの解析コードへの採用

ユーザー承認：002版Purityを確定し、anaFemtoLambda_4He_KFParticleへ反映。d/t/³Heの関数・係数はGit 8f7d539のまま。追加指示：全核種でR<2はTOF不要、R≥2はTOF必須。³He背景評価用k*候補にも同一条件を適用する。

## 実装・検証完了

- 解析専用nuclearid YAMLと⁴He多項式YAMLを作成。共有Git PIDヘッダと他mainconfは変更しない。
- ⁠⁴HeのR>0.61、R<5、旧⁴He vetoなしを採用。d/t/³Heの既存最小R・vetoは維持。TOF窓は核種別の既存値を維持し、⁴Heはd窓。
- 低R純度fitの固定Gaussian中心はfitの条件であり、PIDには確定したlog10(R)多項式を使う。
- STAR SL24y / sl73_x8664_gcc485 / ROOT5コンテナで再build完了。専用回帰14,147チェック、旧設定互換性1,359チェックが通過。

確定資料：[002版report索引](report_he4_fixedmean_20261022.md)。監査出力：rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_adoption_20261022_001/。既存結果は保持。

## 適用先と選別

- 入口：`config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml`。
- 核種設定：`config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_he4_accepted_20261022.yaml`。
- 係数：`config/calibration/dedx_he4_auau3p85_2021P24iy_accepted_20261022.yaml`。`x=log10(R)`、mean3次・sigma2次の関数形を維持。
- mean = `[40.07075429018867, -104.6984111474443, 140.39957818755835, -64.63448390958877]`。
- sigma = `[2.3688048759048863, -6.957225467583047, 7.124880203400031]`。

| 核種 | R<2 GeV/c | R≥2 GeV/c のTOF (m/q)²窓 |
|---|---|---|
| d | TOF不要 | 3.268773 < m² < 3.693147 |
| t | TOF不要 | 7.2564275 < m² < 8.2816925 |
| ³He | TOF不要 | 1.78285 < m² < 2.06485 |
| ⁴He | TOF不要 | dと同じ窓 |

これはTOF要件の表であり、他のトラック品質・最小R・既存veto条件も引き続き適用される。d/t/³Heの旧veto（特にR≥2.5）は今回の変更対象に含めず保持。⁴Heの旧vetoだけは確定studyの選択に合わせて解除した。全核種のbest-species判定を先に行い、下限/TOF等で不合格のwinnerを他核種へ再割当てしない。⁴He係数変更に伴う競合境界の変化は意図した変更。

`FemtoLambdaLegacy`が候補の選別・全4核種pull・selected QAとpre-PID pair QAの最終domain条件を共通管理する。³Heの物理質量によるk*と、同じ³He候補を⁴He質量で扱う背景controlの双方は選別済み候補を使う。SE/MEとも同じ候補storeに基づく。未選別raw QAは従来通り別に保持。

## 検証と再現

確定studyから独立に出力した4核種×1,001点のmean/sigma表と実行時pullを照合。d/t/³Heは共有Git表とbit単位一致。共有PIDヘッダはGit 8f7d539とファイル全体が完全一致し、他mainconf・共有旧cutは作業前SHAと一致。確定studyの係数・settings・summaryとfitソースのSHAも照合。

R=0.61の厳密除外と直上、R=2と直上/直下のTOF要件、strict m²窓、R=5除外、旧⁴He veto域の受入れ、winner再割当て禁止を検証。PicoDst fixtureで⁴He候補4件・³He候補2件を確認し、³HeのTOFなし高R候補と誤ったd窓候補が落ち、両k* controlへ各2件のみ入ることを検証。ROOTへの正確な消費YAML・strict下限・veto/TOFメタデータ保存を読み戻した。欠落ファイル・不正係数・異なる年・負の幅は明示的に失敗し、Init失敗も理由付きで確認。

```bash
./script/singularity_make.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml --no-clean -j2
bash script/singularity_root_batch.sh config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml root4star -b -q 'tests/run_femto_he4_pid_adoption_test.C("config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml","config/calibration/dedx_he4_auau3p85_2021P24iy_accepted_20261022.yaml","rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_adoption_20261022_001/expected_coefficients.tsv","tmp/he4_adoption_recheck/config_output")'
```

テスト出力先は未使用のディレクトリを指定する。実行済み出力は監査先の`runtime_002/`、`test_002.log`、`legacy_domain.log`、`build.log`、`completed.json`。`runtime_001`は初期診断で保持し、最終結果は002。旧設定回帰は変更していない`main_auau3p85_anaFemtoLambda_4He_KFParticle_highpurity.yaml`を入力として`tests/run_femto_nuclear_pid_domain_test.C`を実行した。

実データの新規event処理・farm投入・hadd・commit/pushは実施していない。採用元Purityの統計的/モデル上の制限は002版reportの通り。今回の完了は設定・解析選別処理への反映とローカル回帰検証を指す。
