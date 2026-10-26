# GitのPID多項式から再開する核PID study（2026-10-20）

## 最新の実行結果

[共通多項式とALL純度の結果まとめ](../../Todo/2026/10/20/summary_polynomial_all_purity_20261020.md)
を参照。全380 TOF局所fit、同じ次数の候補係数、全352点×3条件のALL/TOF分離、
旧新COLZ・全点linear/log・Purity vs Rを出力した。制約付きGLSでも同じ比較を実施した。
ただし低Rのfit品質が足りず、新しいcutoff/高R最終上限は本番へ採用していない。

### 追加CLI

すべてプロジェクト直下から `python3 -B tools/nuclear_pid_restart_20261020/` に続けて指定する。

| Entry | Usage |
|---|---|
| `recenter_survey.py` | `--settings recenter_survey.yaml --output-dir NEWDIR`。全raw exportから局所ROIだけを準備 |
| `polynomial_calibration.py` | `--fit-dir DIR --settings polynomial_calibration.yaml --output-dir NEWDIR`。候補係数・支持点・感度比較 |
| `shared_curve_study.py prepare` | `--settings shared_curve_study.yaml --coefficient-dir DIR --variant linear/constant/wider --output-dir NEWDIR`。`--smoke`で代表8点 |
| `shared_curve_study.py fit` | `--output-dir DIR`。凍結した設定・係数でROOT5 fit、保存/readback |
| `shared_curve_report.py tests` | `--fit-dir DIR --output-dir NEWDIR`。ROOT5人工試験 |
| `shared_curve_report.py overlay` | 上記に`--coefficient-dir DIR`を追加。旧新COLZ |
| `shared_curve_report.py fits` | 全核種の保存予測を局所linear/log表示。refit/rebinしない |
| `shared_curve_report.py compare` | `--oracle-dir DIR --all-dir DIR`を追加。エンジン同値性／曲線圧縮の切り分け |
| `shared_curve_report.py summary` | `--other-fit-dir DIR`を複数指定。良好モデルのみの純度・未評価区間・採否YAML |

`--diagnostic-mode tof_fixed/tof_oracle` はTOF同値性・圧縮誤差の診断専用で、
全RにTOFを要求する。純度を採用せず、通常のALL低R/TOF高R studyと混同しない。
正の曲線や数値収束だけで本番採用はしない。GLSの制約付き候補ではLOOを省略し、
forward fitを別途検証したことをJSONに明示している。

新規出力のみを許可し、元ROOTと古いPDFは上書きしない。
`shared_curve_study.py` はTOF支持範囲外へ曲線を黙って外挿しない。
誤差は固定曲線での条件付き尤度曲率による近似で、PID曲線そのものの系統誤差ではない。

## この作業の位置づけ

10月13〜19日のfit・補間値を新しい校正の初期値として継承せず、Gitに保存された
[`NuclearIdDeDxVsMom.h`](../../StMaker/common/NuclearIdDeDxVsMom.h) のPID関数へ戻して進める。
中心値は `log10(p/q)` の3次式、幅は同じ変数の2次式を維持し、必要な場合だけ係数を更新する。
独立なスライスごとの中心値・幅を、そのまま折れ線の本番PIDへ採用しない。

初期段階で**6つの3.85 GeV本番入口の元PIDへの復帰と実行時照合、380スライスのGaussian診断fit、
190ページPDF、局所範囲の比較**を完了した。高χ²と範囲依存が残るため、
ユーザーと相談してGaussian＋左右指数tailも代表8点で比較・保存した。
その後の候補係数・純度評価は冒頭の最新結果へ記載した。cutoffの本番認定は未完了である。
旧研究の結果は削除せず、今回の採用済み結果とは区別して保存する。

## 復帰と検証

- 本番入口は新規 `gitpid_restart` / `gitpid_tofExtended` 設定へ付け替えた。
  日付付きの旧設定と、変更前mainconf一式は保持した。
- 新設定では追加dE/dx校正を無効化し、Gitの中心値・幅を使用する。
  **既存のcutoff、除外区間、TOF選択まで検証済みの新条件へ置き換えたわけではない。**
- SL24y / ROOT 5 の実行時providerで、6入口を通した864行・2592項目の照合に成功した。
  静的回帰試験は16件成功した。これはPID値の復帰確認であり、純度や信号効率の認定ではない。
- 特別入口 `tofExtended` の旧 `3 < m² < 4` は今回の復帰で勝手に変更しない。
  **新しい共通studyでは、⁴HeのTOF選択をdと同じm²窓にした。**

保存先：

- [変更前の設定](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_git_restore_20261020_001/before/)
- [復帰監査](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_git_restore_20261020_001/restore_audit.json)
- [実行時providerの検証](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_git_restore_20261020_001/runtime_provider.json)

## 新しい4段階のstudy

1. **TOFタグから係数を確認する。** d・t・³Heをm²で選び、⁴Heにはdと同じm²窓を使う。
   各p/qスライスのdE/dxをGitのPID関数からfitし、式の次数を変えずに係数を調整する。
   dと⁴Heはm²だけでは分離できないため、両者のdE/dxピークを別成分として扱う。
2. **TOFを要求しないALL分布をfitする。** 段階1の校正を基準とし、低p/qを含む各核種の
   ピークを確認する。他核種の成分も共通のPID関数と整合させ、任意に動く背景ピークで吸収しない。
3. **低p/qの使用下限を決める。** `P = S/(S+B)` が0.5未満になる範囲をPID対象から外す。
   不良fit・未収束・統計不足は低純度と同一視せず「未評価」に分ける。
   前後と矛盾する孤立した除外区間は自動採用せず、元のdE/dx fitを提示して判断する。
4. **p/q対純度を統合表示する。** `p/q >= 2 GeV/c` ではTOF選択後の純度を算出する。
   高p/q上限は保存データの実際の統計とfitの成立状況から決め、以前の上限をそのまま流用しない。

全段階で、各p/qのlinear/log dE/dx fitと、二次元分布へのPID曲線重ね書きを残す。
対象ピークと隣接ピークが分かる局所表示とし、空白の大きい縦軸や不要な0–120 keV/cm表示で
ピークが潰れないようにする。

## 実行入口

新しいstudyの設定は [`tag_study.yaml`](tag_study.yaml)、準備・実行の入口は
[`prepare_tag_study.py`](prepare_tag_study.py)。処理の順序は次のとおり。

```text
prepare --settings <tag_study.yaml> --output-dir <新規出力先> --smoke
baseline --output-dir <同じ出力先>
fit --output-dir <同じ出力先>
```

各サブコマンドは `python3 -B tools/nuclear_pid_restart_20261020/prepare_tag_study.py`
の後に指定する。上記はCLIの形式であり、未実行段階の成功を表すものではない。
入力は保存済みrawヒストグラムのcacheのみとし、過去のfitパラメータ・補間値は使用しない。
準備段階の設定と入力識別情報を確認してからfitへ進む。

`--smoke` を省略すると初回surveyの95点×4核種を準備する。
R=0.405〜4.995 GeV/cは全4元曲線の定義域を調べるための範囲であり、PID下限／上限の新設定ではない。
`tag_study_local15.yaml` はfit範囲だけを狭める比較用で、品質を改善したことにはしない。

保存済みfitを再fitせず描画する：

```text
python3 -B tools/nuclear_pid_restart_20261020/render_tag_study.py \
  --fit-dir <fit済み出力先> --output-dir <その直下の新規views先> \
  --settings tools/nuclear_pid_restart_20261020/view_settings.yaml
python3 -B tools/nuclear_pid_restart_20261020/summarize_tag_study.py \
  --fit-dir <fit済み出力先> --output-file <未使用の集計JSON>
```

描画ではnative count、誤差、保存された予測だけを使用し、正規化や曲線の外挿をしない。
χ²感度表の2／3／5／10は記述的な比較値であり、採用基準を自動設定するものではない。

## 新しいピーク応答の比較

設定は [`response_study.yaml`](response_study.yaml)。準備時に
`--response-variant gaussian`／`left_a1`／`left_a25`／`right_a1`／`right_a25` を明示する。
`--settings` は上のresponse設定、fit入口は同じ `prepare_tag_study.py` を使い、
各条件は別の新規出力先へ保存する。実行時は固定済みplanを使用するので、設定・variantの変更はprepare時だけ。

- 同じ4核種×2点、同じnative binとfit範囲で比較する。
- Gaussian、左指数tail付きGaussian、右指数tail付きGaussianを比較する。
  tail形状の初期値2通りは局所解の確認であり、別の物理モデルとして数えない。
- 新応答でも対象核種の中心とcore幅だけを校正対象とし、他核種を任意のfreeピークにしない。
  これは純度fitではなく、他核種の中心がfit範囲内にある場合は別途取り扱いが必要。
- 中心は応答のmode、σはGaussian core幅。非Gaussian応答の全RMSと同一視しない。
  本番のnσへ採用する前に、幅の意味と選択帯を改めて検証する。
- MinimizeとHesseの戻り値・statusを分離し、Hesse成功が最適化未収束を上書きしない。
  形状共分散と固定形状での条件付き収量共分散を、計算せず有効と表示しない。
- 新しいtail数学試験はROOT5で4585項目成功。意図的に未収束を作る試験も含む。
  この成功は実データのfit品質や純度の採用を意味しない。

実データ・PDFの一覧は[出力README](../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/README.md)。

モデルを横3列、linear／logを上下2行に比較する入口は
[`compare_response_studies.py`](compare_response_studies.py)。
`--fit-dir` を3〜5条件分指定し、`--output-dir` はそれらと同じ親の新規directory、
`--settings` は `view_settings.yaml` を明示する。再fit・ROOT再生成・正規化はしない。
3列でraw count／error／flow／軸の一致も検査する。
family内の最良初期値解は描画のために選ぶもので、科学的な採用とは区別する。

native-R nodeの実際のtail境界も保存fitから求め、tailとcoreの両方がfit範囲で見えない解を
「αを識別できない」として比較判定へ追加する。元fitフラグや予測を後から書き換えない。
過去のソースは実行snapshotのSHAで検査する。現在のruntime YAMLが変わっている場合、
コメントだけの訂正は値の完全一致を確認して許容するが、物理値・参照先の変更は拒否する。
wrapperは厳密なSHA一致を要求する。

104 Pythonテスト、12 native数学試験、4585応答数学試験が成功。
ただし実データには高χ²が残り、応答・係数・純度・cutoffの採用は未完了。

本番Makerや本番PIDへ新係数を適用する判断は、fit・純度・回帰試験の確認後に別途記録する。
現在の進捗は [10/20 Todo](../../Todo/2026/10/20/README.md) と
[再開作業レポート](../../Todo/2026/10/20/report_pid_restart_20261020.md) を参照。

## 保存精度・局所範囲・非対称応答の追加診断

すべて診断用であり、元ROOT・GitのPID多項式・本番cutoffを変更しない。
2021 P24iyの代表PicoDst 2本で、dE/dxのFloat16既定12bit保存を確認した。
native countを補正するのではなく、予測積分に保存丸めを含める。

| 入口／設定 | 役割 |
|---|---|
| `export_restart_spectra.py --fit-dir ... --output-dir ...` | 同じraw cacheの全bin count／errorを書出し。flow・ゼロbinも保持し再読込検証 |
| `verify_storage_mapping.py --output-dir ...` | ROOT5自身のFloat16処理と数学モデルを照合。新規出力先・buildのみ |
| `storage_study.yaml`／`FitStorageSlices.C` | 同じGaussianで保存丸めの有無を切り分ける。nominal／effective bin端を保存 |
| `split_normal_study.yaml`／`FitSplitNormalSlices.C` | 共通mode、左右のGaussian幅を分けた単一核種ピーク。背景核種のfreeピーク追加ではない |
| `verify_split_normal.py --fit-dir ... --output-dir ...` | 固定済みsnapshotと設定で人工データ試験。実データをfitしない |
| `recenter_tag_study.py --settings ... --variant local15 --output-dir ...` | 代表8点のfit範囲だけを実ピークへ移す準備。fitは既存prepare入口で実施 |
| `recenter_study.yaml`／`recenter_split_study.yaml` | Gaussian／非対称応答で、同じROI移動を比較する設定 |
| `assess_calibration_support.py --fit-dir ... --output-dir ... --chi2-levels 2 3 5` | 元fitを変えず、ゼロ信号・誤差・収束・χ²・範囲の支持を後処理で検査 |

`prepare_tag_study.py prepare --settings <storageまたはsplit設定> --output-dir <新規sibling> --smoke`
で8点、`--smoke`なしで380点を固定し、同じ入口の`fit`で実行する。
ROI移動は3bin移動和で範囲中心を探すだけで、fit dataは平滑化・rebinしない。
初期値・GL3 nodes・元曲線・診断gateはGitのまま保持する。

非対称応答の中心はmode、幅は`σ = sqrt(σLeft*σRight)`、追加の比は`σRight/σLeft`。
全分布のRMSや、元のGaussian標準偏差と同じだと扱わない。
本番nσへそのまま採用するのではなく、幅の定義・選択帯・純度を合わせて確認する。

χ²の2／3／5は支持点数を比較するための値であり、採用基準の決定ではない。
小さなχ²でもtarget収量がゼロならピークは測定できていない。
巨大だが有限なHesse誤差も「中心・幅を同定できた」ことを保証しない。
これらを低純度や孤立vetoへ読み替えず、未評価として分離する。

ユーザー方針：ここでは応答モデルをさらに増やさず、現在の応答と範囲比較を先にまとめる。
この段階の最終Python回帰試験128件、Float16照合57075項目、非対称応答人工試験6773項目が成功。
比較結果は[現在の応答のまとめ](../../Todo/2026/10/20/summary_current_responses_20261020.md)。
全surveyのLocal ROI追試・係数校正・ALL純度は未実施／未採用部分として残している。
