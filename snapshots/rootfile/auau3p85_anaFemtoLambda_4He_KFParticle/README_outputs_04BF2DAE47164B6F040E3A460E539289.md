# 2021 P24iy：通常ヒスト統合ROOTと8次元PID分割ROOT

この保存形式は2026-10-26にユーザー承認。データはSUMS 04BF2DAE47164B6F040E3A460E539289の正常終了26,955/26,994 subjobs。解析コード・PID・選択条件は変更していない。ファイルの完成判定は同じ出力先のmanifestと監視記録result.json（phase=merged）、archive_result.json、sent.jsonを確認する。

- 通常の1,786ヒスト：`auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root`。巨大な`hNuclearPidPairJoint`はこのROOTから除き、代わりに`NuclearPidPairJointSplitManifest`を保存する。
- 巨大Sparse：`auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_pid_pair_joint_part000.root`、`part001.root`、`part002.root`。各ファイルには同じ名前・型・8軸の`hNuclearPidPairJoint`を保存し、bin内容・座標・誤差²・Sumw2・Entries・統計量を元の中間ROOTと完全一致させる。
- `part000/001/002`は入力イベント群の分割（正常subjobs 10,000 / 10,000 / 6,955件）。Rや核種の領域分割ではない。同じ座標binが複数partに現れるので、必要な射影を各partで作り、射影を加算する。
- 出力manifest：`auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_outputs.json`。Sparse一覧：`auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_pid_pair_joint_parts.list`。射影helper：`auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_project_joint_parts.C`。
- 通常ROOTと元subjob・中間ROOTを一緒にhaddしない。3 Sparseを再び8Dの1オブジェクトとして書き出すと約1GBのROOT上限に再度達する。元データ・中間ROOTは保持する。

## hNuclearPidPairJointの意味

Λ–原子核ペアのk*、原子核R=p/|q|、dE/dx nσの同時分布を保持するペアPID QA。Purityのfit結果は保存していない。後からR・nσを選び、k*分布やPID背景の影響を調べるための生カウント。

| 軸 | 内容 |
|---:|---|
| 0 | Kstar：200bin、0–1 GeV/c |
| 1 | POverQ：120bin、0–6 GeV/c |
| 2 | NSigma：120bin、−12–12 |
| 3 | Cent9：0–8 |
| 4 | 0=⁴He、1=³He、2=³Heを⁴He質量として扱う背景表現 |
| 5 | 0=SE、1=ME |
| 6 | Λ質量窓：0=左SB、1=signalLow、2=signalHigh、3=右SB |
| 7 | 0=selected、1=prePID |

prePIDは核種自身のTPC nσ条件とbest-species競合を外すが、対象R領域やTOF条件等は保持。he3as4のRとnσは元の³Heを使い、k*を⁴He質量で計算し直す。各記録はweight=1でSumw2付き。signalはMassClassの1+2で復元し、独立カテゴリとして二重に保存していない。selected/prePIDやhe3/he3as4は別区分なので、全カテゴリを足したEntriesは独立な物理ペア数ではない。

## 射影例

ROOT5でhelperを明示ロードする。ファイル名と関数名が異なるため`.L`してから関数を呼ぶ。次はcent9=8、⁴He、SE、selected、Λ信号窓のk*を全3partから加算する例。出力先は新しいファイル名にする。

```cpp
.L auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_project_joint_parts.C+
project_joint_parts("auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_pid_pair_joint_parts.list", "he4_se_selected_cent8.root", 0, 8, 0, 0, 0, 1, 2);
```

引数順：partsList, output, 射影軸（0=k*,1=R,2=nσ）, cent9（−1=全cent）, representation, SE/ME, PID stage, Λ massClass低端, 高端。R・nσなど物理軸のflowを含めて射影し、内容と誤差²を加算する。

監視記録：`job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289/split_sparse_final_20261026/`。


## 自動完了記録 — 2026-10-26T11:13:22.408564+09:00

状態：完了（通常ROOT＋3 Sparse、個別ROOT保存済み）。メール：ローカルsendmail受付済み。
通常ROOT：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_merge.root`（14064010756 bytes）。
Sparse：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_pid_pair_joint_part000.root`。
Sparse：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_pid_pair_joint_part001.root`。
Sparse：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289_pid_pair_joint_part002.root`。
検証・公開・保存・送信記録：`/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/job/run/merge_notify/auau3p85_anaFemtoLambda_4He_KFParticle_04BF2DAE47164B6F040E3A460E539289/split_sparse_final_20261026`。
