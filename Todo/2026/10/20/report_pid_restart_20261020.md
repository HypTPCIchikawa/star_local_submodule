# 核PIDをGitの元式から再開する作業記録（2026-10-20）

**今回の区切りは[現在の応答・fit範囲比較のまとめ](summary_current_responses_20261020.md)。**
以下は段階ごとの履歴。最新の比較表・品質表・PDFと未完了事項は上記に整理した。
ユーザー方針により応答モデルをさらに増やさず、実施済みの比較までで一旦まとめる。

## 依頼と方針

[本日のTodo](STAR_Todo20261020.txt) に従い、10月13〜19日に得た校正・fitの結論を
新しい初期値として引き継がず、Gitの元PIDから仕切り直す。
過去のファイルは削除せず、今回の採用条件や検証結果と区別して保持する。

中心・幅の関数形は元の多項式を維持する。TOF選択した分布、TOFを要求しない分布の順に検証し、
`P = S/(S+B)` と不良fitを区別して核種別の低p/q下限を決める。
高p/qではTOFを組み合わせ、統計とfit成立範囲を確認して上限を決める。

## 完了：解析入口のPIDを元式へ復帰

対象は3.85 GeV / 2021 P24iy用の6入口。追加校正がGitのPID値を上書きしないよう、
新規の `gitpid_restart` / `gitpid_tofExtended` 設定を用意してmainconfから参照させた。
元の日付付き設定は編集で潰さず、変更前の設定一式も別途保存した。

元式は [`NuclearIdDeDxVsMom.h`](../../../../StMaker/common/NuclearIdDeDxVsMom.h)。
中心値の `log10(p/q)` 3次式、幅の2次式を使う。
追加校正の `enabled` を無効化して元式を参照する復帰であり、
以前の純度研究に由来するcutoff・除外区間・TOF条件を新たに認定したものではない。
それらの再評価はこれから行う。

| 確認対象 | 結果 | この結果から言えること |
|---|---|---|
| 設定・providerの静的回帰試験 | 16件成功 | 復帰設定と参照経路の整合 |
| SL24y / ROOT 5実行時provider | 6入口、864行、2592項目の照合成功 | 実際に返すPID中心・幅等が元式の基準と一致 |
| TOFタグのGaussian診断fit | 380スライス保存、readback成功 | fit品質の採用や新係数の認定ではない |
| 新しい多項式係数・純度・cutoff | 未認定 | 不良fitから本番条件を変更しない |

検証・保存物：

- [変更前設定の保存](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_git_restore_20261020_001/before/)
- [復帰監査JSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_git_restore_20261020_001/restore_audit.json)
- [実行時provider検証JSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_git_restore_20261020_001/runtime_provider.json)

特別な `tofExtended` 入口は、復帰時点では旧 `3 < m² < 4` を保持している。
これは今回の新しいstudyの選択条件とは別であり、新studyでは⁴Heにdと同じm²窓を使う。
新study用の入力・条件と、既存入口に残る条件を混同しない。

## 完了：Git曲線の重ね描きと新しいGaussian診断

生のALL、d／t／³Heのm²タグを持つnative cacheのSHAを照合して使用した。
以前の係数・fit点・purity・除外区間は入力にしていない。
対象は2021 P24iyの既存部分統合ROOTから保存された、電荷符号・centralityを分けていないtrack QA。
TOFタグはALLの部分集合なので、別の独立標本として重複加算しない。

| 保存物 | 内容 |
|---|---|
| [Git元曲線と生COLZ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/tof_tag_smoke_002/git_pid_baseline.pdf) | ALL／3種m²タグ × 低R／高R、8ページ。元中心と±2σ |
| [全RのGaussian fit集](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/tof_tag_survey_001/views_001/tag_peak_fit_diagnostics.pdf) | 各核種95点、計380スライス。局所linear/log、190ページ |
| [代表8点のfit集](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/tof_tag_smoke_002/views_001/tag_peak_fit_diagnostics.pdf) | R=0.875、1.775 GeV/c、4ページ |
| [狭い範囲の比較](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/tof_tag_local15_smoke_001/views_001/tag_peak_fit_diagnostics.pdf) | 同じ8点でfit範囲だけ±2.5σから±1.5σへ変更。4ページ |
| [品質・統計の集計](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/tof_tag_survey_001/diagnostic_quality_summary.json) | 収束状態、χ²感度表、初期値→fit値、入力・出力SHA |

surveyのRは0.405〜4.995 GeV/cのnative R bin（幅0.01 GeV/c）から明示抽出した95点。
全4核種の元曲線が定義される範囲を最初のsurveyとしただけで、**0.4を新しい下限にしたのではない**。
d・tのさらに低Rの評価は未実施。native dE/dx binは0.1 keV/cmのままで、fit時のrebinはしない。
有限R bin内の元曲線の変化は3点の積分nodeで扱う。

380行すべての計算・保存・native count/error/flow・成分加算の再読取が完了した。
うち376行が数値収束等のflagを満たすが、これはχ²による科学的採用ではない。
125行では別核種の元ピーク中心が対象fit範囲内にあり、対象ピークだけのfitを純度評価には使えない。
この段階では他核種の中心を自由に動かしたり、純度を計算したりしていない。

### Gaussianのfit範囲比較

表の範囲はGitの初期σを基準にnative bin端へ揃えたもの。数値はχ²/ndf。

| R [GeV/c] | Species | ±2.5σ range | ±1.5σ range |
|---|---|---:|---:|
| 0.875 | d | 3190.34 | 209.00 |
| 0.875 | t | 106.76 | 21.84 |
| 0.875 | ³He | 391.91 | 93.22 |
| 0.875 | ⁴He | 23.99 | 9.73 |
| 1.775 | d | 284.93 | 29.79 |
| 1.775 | t | 196.31 | 7.49 |
| 1.775 | ³He | 18.70 | 12.89 |
| 1.775 | ⁴He | 1.93 | 2.22 |

例えばR=0.875の³He中心は広いfitで34.273、狭いfitで34.626 keV/cm、
⁴Heは47.130→47.517 keV/cmとなり、単に範囲を狭めるだけでは形と中心の安定性を保証できない。
これらの不良fitから係数を更新したり、低純度・孤立vetoを認定したりしない。

ユーザーにこの結果を示して相談し、**ピーク形のモデル比較を追加する**方針で合意した。
次はGaussianと片側の指数tail付きGaussianを同じnative bin・同じ範囲で比較する。
PID曲線は元の中心3次・幅2次多項式のまま。tail応答のcore幅と全分布のRMSを区別し、
今回の実データから得たfitの品質・安定性を確認するまで本番へ反映しない。

検証：Python単体テスト104件、ROOT5 native bin等の人工データ試験12件が成功。
元ROOTの再解析、job投入、hadd、git commit/pushは行っていない。

## 完了：Gaussianと左右の指数tail応答の比較

- [同じ生データで3モデルを比較したPDF（8ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_003/response_comparison.pdf)
- [比較表：中心・core幅・α・誤差・χ²・収束・初期値差](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_003/comparison.tsv)
- [詳細JSONと入力不変性の記録](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_003/summary.json)

Gaussianと左右の指数tail付きGaussianを、4核種×2つのRで比較した。
tailのα初期値は1.0／2.5の2通りを別出力で実行し、合計5条件40 fitを保存した。
同じモデルの初期値2通りを、独立した物理モデルとは数えない。
描画する解は数値条件を確認した上で目的関数から選ぶが、**科学的な採用ではない**。
各ページは3モデルを横3列、linear／logを上下2行とし、8ページ・48パネル。
同じRのnative bin、count、error、flowが3列で完全に一致することをROOTで照合した。
全ページをGhostscriptで再読取し、元データ・保存fitは変更していない。

表はχ²/ndf。「α未拘束」は、全native-R nodeのtail／core境界がfit範囲内で見えず、
その範囲のデータではαを識別できない解。χ²が低純度を示すわけではない。

| R [GeV/c] | Species | Gaussian | Left exponential tail | Right exponential tail |
|---|---|---:|---:|---:|
| 0.875 | d | 3190.34 | α未拘束 | α未拘束 |
| 0.875 | t | 106.76 | α未拘束 | α未拘束 |
| 0.875 | ³He | 391.91 | 218.07 | 396.22 |
| 0.875 | ⁴He | 23.99 | 25.50 | α未拘束 |
| 1.775 | d | 284.93 | 303.77 | α未拘束 |
| 1.775 | t | 196.31 | 181.11 | α未拘束 |
| 1.775 | ³He | 18.70 | 17.66 | 15.69 |
| 1.775 | ⁴He | 1.93 | 1.95 | 1.92 |

tailの最適化（MIGRAD）と誤差計算（Hesse）は別々に保存し、Hesse成功でMIGRAD不良を上書きしない。
さらに、Hesseが有限の誤差を返しても、tailがfit範囲の外ならαは実際には決められない。
この追加診断は保存node・fit値から後処理で求め、元fitのフラグや予測は書き換えず、
比較の数値OK判定へ別条件として加えた。数学的にtailを見られない解を校正へ採用しない。

例えば³HeのR=0.875は左tailで改善するが、χ²/ndf≈218が残る。
⁴Heの同じRでは十分な改善がない。**今回の3形状から低Rの校正や純度を確定するには不十分**。
一方、R=1.775の⁴Heは3モデルともχ²/ndf≈1.9で、形状により成否が変わる例と区別できる。
これらを低純度・純度ゼロ・途中のvetoという判断へ置き換えない。

中心は応答のmode、σはGaussian core幅として比較した。全RMSとは異なるため、
新しいnσの本番選択帯へそのまま転用せず、曲線・選択帯・純度を合わせて検証する。
新しい多項式係数、ALLの純度、核種別の下限／上限はまだ認定していない。

tail数学試験はROOT5で4585項目成功した。これは実データfitの採用ではなく、
積分・規格化・人工データ回収・誤差・未収束状態の取り扱いの検証。
本番6入口の復帰は[コメント整理後の再監査](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_git_restore_20261020_001/restore_audit_after_comments.json)も成功。
入口の古い「校正を使用」とするコメントは現在のGit元PIDに合わせて直したが、追加の物理値変更はない。

## 完了：保存精度による周期構造の切り分け

入力とモデルの取り違えを調べ直し、raw QA／m²タグは核種nσ選択より前であること、
R bin・dE/dx単位・projectionの対応を再確認した。
代表8スライスの全9,616 bin（flowを含む）を元ROOTから書き出し、count・errorの完全一致を確認した。
最大binは1,961,877、すべて整数countでerror²=countと整合する。
この代表範囲ではTH2Fの整数表現上限による飽和は見られない。

`StPicoTrack::mDedx`は`Float16_t`であり、同じ2021 P24iyのローカル実PicoDst 2本を
イベントを読まずに検査した。writer ROOT 5.34/38、StPicoTrack v9、checksum 500008820、
range／factorなしで、既定12bit mantissa保存に対応する。
この2本は欠損job復旧用の代表キャッシュであり、部分統合に入った全PicoDstのschema監査ではない。
ROOTの丸め処理は[公式実装](https://root.cern/doc/master/TBufferFile_8cxx_source.html)でも確認できる。

この保存形式をモデル化し、**観測値・native binを変えずに、信号と背景の予測積分範囲だけを
保存丸めの逆像へ変換**した。2の冪を跨ぐ丸め・clamp、ROOT自身とのroundtrip、native bin所属等を
ROOT5で57,075項目検証した。再binやデータの平滑化ではない。
nominal／effective両方のbin端をTSVへ記録し、保存ROOT・JSON・TSVを再読込した。

| R [GeV/c] | Species | Continuous prediction χ²/ndf | Float16-aware prediction χ²/ndf |
|---|---|---:|---:|
| 0.875 | d | 3190.34 | 3152.43 |
| 0.875 | t | 106.76 | 95.06 |
| 0.875 | ³He | 391.91 | 362.86 |
| 0.875 | ⁴He | 23.99 | 21.09 |
| 1.775 | d | 284.93 | 234.80 |
| 1.775 | t | 196.31 | 197.72 |
| 1.775 | ³He | 18.70 | 11.38 |
| 1.775 | ⁴He | 1.93 | 1.30 |

同じGaussian＋局所線形背景、同じGit初期値・同じfit範囲で比較した。
細かい周期構造には説明が付くが、**巨大χ²のすべてを保存精度で説明できたわけではない**。
中心・幅の変化は小さく、滑らかな非対称ピーク形状とfit範囲の依存を引き続き確認する。
保存精度の効果を根拠にχ²の閾値だけを緩めたり、本番PID係数を変えたりしない。

- [量子化を扱った8点のlinear／log PDF（4ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/storage_default12_smoke_001/views_001/tag_peak_fit_diagnostics.pdf)
- [native／effective bin対応表](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/storage_default12_smoke_001/native_effective_bins.tsv)
- [ROOT5との数学照合ログ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/storage_math_001/test.log)
- [実PicoDstの保存形式確認](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/pico_streamer_001/onfile_mDedx.json)
- [元count／errorの書出し・再読込確認](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/native_storage_export_001/completed.json)

## 完了：非対称応答と実ピーク中心のfit範囲比較

共通modeの左右Gaussian幅を分けた応答を診断用に追加し、ROOT5人工データ試験6,773項目に成功。
同じGit初期値で380スライスの計算・保存・raw／成分予測／bin対応表readbackと、190ページPDFを検査した。
代表8点では、実ピーク中心±1.5元σへROIだけを動かす比較も実施した。
rawを平滑化・再binせず、初期値・元曲線・nodes・診断gateを保持した。

R=0.875の³He／⁴Heは組合せによりχ²/ndf≈1.94／1.05となる一方、d／tは≈19.57／8.86。
小さなχ²でも信号・背景の帰属や新しい幅が確定したわけではなく、全域へ採用しない。
幅は左右幅の幾何平均で、元Gaussianのσ・全RMSと同じとはしない。

元Gaussianの`t_b0061`と`he4_b0051`はtarget yield=0でもstatus=0、χ²≈1であった。
こうした見かけ上の良好fitを校正へ使わないため、収束・正のtarget統計・誤差・他ピーク混入・
modeがfit範囲内かを後処理で検査し、別表へ保存した。χ²の2／3／5は記述的な比較値である。
最終Python回帰試験128件と、本番6入口の診断後静的再監査が成功。
全survey Local ROI追試・多項式再校正・ALL純度・cutoffは未完了。

## 続けて行うstudy

1. raw cacheのd・t・³Heのm²タグ、dと同じm²タグ中の⁴Heについて、p/qスライスのdE/dxをfitする。
   元の多項式を基準に、必要な場合は次数を変えずに係数だけを更新する。
2. TOFを要求しないALL分布を、更新した各核種のPID中心・幅と整合する成分でfitする。
   特に低p/qで、³He・⁴Heと他核種の実ピークを確認する。
3. 良好なfitで得た純度が0.5未満になる低p/q範囲を打ち切る。
   不良fitを純度ゼロへ置き換えず、未評価とする。孤立した除外区間は元のfitと前後の区間を提示する。
4. 低p/qのALLと `p/q >= 2 GeV/c` のTOF選択後の純度を一枚にまとめる。
   高p/q上限は、実統計とfitの成立状況から判断する。

新しい入口は [studyツールの説明](../../../../tools/nuclear_pid_restart_20261020/README.md)、
[準備・実行スクリプト](../../../../tools/nuclear_pid_restart_20261020/prepare_tag_study.py)、
[設定YAML](../../../../tools/nuclear_pid_restart_20261020/tag_study.yaml) を参照。
TOF Gaussian診断と3応答の代表点比較までは実施済み。
続きは、残るピーク形・局所残差の不一致を解決し、採用可能な応答・fit品質基準を確認するところから。
不良fitを使って先に係数を滑らかにしたり、純度図の空白を埋めたりしない。

## 完了：10/13 Todoのディレクトリ整理

`Todo/20261013/` と `Todo/STAR_Todo_20261013.txt` を
[`Todo/2026/10/13/`](../13/README.md) へ移動した。

- 移動したファイル：155件。
- 現行Markdownのナビゲーション修正：17ファイル。
- 移動先・相対リンク等の確認：223件、成功。
- 移行スクリプトの回帰試験：11件成功。
- 非Markdown、元のTodo TXT、validation内のsnapshot・ログ等は内容を保持した。
  ROOT、PDF、ジョブ提出時のsnapshotは移動・変更していない。

[移行manifest](../todo13_migration_20261020.json) に移動前後のSHA-256と検証結果を保存した。
これは**移行した時点の記録**である。今後の正当なMarkdown進捗更新はその時点のSHA-256と異なり得るが、
元TXT・ログ・snapshotを変更しなかった根拠は保持されている。

## まだ完了していないこと

- 残るピーク応答・残差の不一致の解決と、良好なTOF fitからの多項式係数の再校正。
- ALL分布の再fit、低p/q cutoffと高p/q上限の再決定。
- 新しいp/q対純度図と、各p/qのdE/dx fit集。
- 新係数・cutoffを本番へ適用する前の回帰試験と適用判断。

元PIDへ戻せたことだけをもって、純度改善や新しい使用範囲が検証されたとはしない。
