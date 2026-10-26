# 新旧核PID関数とTOF窓の比較（2026年10月17日）

- [比較PDF](pid_parameter_comparison.pdf)：9ページ。設定一覧、4核種の中心・幅・中心差・幅比、
  ³He/⁴He低R拡大、ALL dE/dx COLZ、TOF m² COLZ。
- [中心・幅の細かい数値表](response_curves.tsv)、[代表Rの表](reference_parameters.tsv)、
  [TOF窓と旧参照mean/sigma](tof_windows.tsv)。
- [入力・設定来歴](context.json)、[native TOF map検査](tof_m2_audit.json)、
  [完了・SHA256](completed.json)、[再現command](command.json)、[ソースsnapshot](source/)。

## 「旧設定」の区別

| 表示 | 実際に比較した核種設定 | 中心・幅 |
|---|---|---|
| `Old calibrated (canonical He4)` | 変更前canonical ⁴He / ⁴He highpurity入口の`calibrated.yaml` | ³He/⁴Heはfrom040校正、d/tは旧log-polynomial |
| `Old legacy (dedicated d/t/He3)` | 専用d/t/³He highpurity入口の`legacy.yaml` | 4核種とも旧log-polynomial |
| `New common calibrated` | `purity50_20261015.yaml` | 旧canonicalと同じ中心・幅。下限/veto/TOF要求が別途変わる |

**今回の比較は再fitではなく、設定されている関数の比較**。
黒（旧canonical）と赤（新common）の中心・幅は一致する。旧専用legacy（青）とは³He/⁴Heで異なる。
d/t自身の関数は変わらないが、競合する³He/⁴Heのpullが変わればbest-species境界は変わり得る。
同じ関数であっても、許可R範囲、veto、TOF条件が違えば選択される候補数は同じにならない。

各核種ページの±2σ線は**PIDゲートの端**であり、fit誤差や関数の信頼区間ではない。
R>2.5ではfrom040校正domain外なので、³He/⁴Heは既存のlegacy fallbackを表示している。
描画はR=0.4–5.0だが、d/t/³Heの最終候補にR≥2.5を許可する意味ではない。

## TOF m²について

| 核種 | 旧TOF参照mean | 旧sigma | 今回studyのstrict m²窓 | 以前のHe4 tofExtended／現在のHe4設定 |
|---|---:|---:|---|---|
| d | 3.48096 | 0.141458 | 3.268773–3.693147（±1.5σ） | d窓は同じ |
| t | 7.76906 | 0.341755 | 7.2564275–8.2816925（±1.5σ） | t窓は同じ |
| ³He | 1.92385 | 0.094 | 1.78285–2.06485（±1.5σ） | ³He窓は同じ |
| ⁴He | 再fitなし | 再fitなし | dと同じ3.268773–3.693147 | 3.0–4.0（studyと異なる） |

**TOFのmean・sigmaを新たに測定・再fitしたものではない**。
d/t/³Heの数値は既存YAMLの参照値をそのまま表示し、⁴Heの[3,4]をGaussian mean±sigmaとは扱わない。
旧canonical／旧legacyのCF核選択は`m2_selection:false`のTPC-onlyで、旧TOF窓はQAタグの条件だった。
`tofExtended`は⁴Heに限る別入口。今回の高R studyでは全核種に明示的なTOF条件を要求している。

TOF COLZは、元2021 P24iy部分統合ROOTから **`true/hM2P`の1キーのみ**を読み、
[tof_m2_native.root](tof_m2_native.root)にそのまま複製したもの。
bin内容・誤差・Entries・flowの全件読み戻し一致を確認した。元10.5 GB ROOT全体を読み直し・hashしない。

このTOF mapにはdE/dxとのjoint情報がないので、窓を重ねただけでTPC+TOF選択後のpurityを測れるわけではない。
⁴Heのstudy狭窓と、現在の広い3–4窓のpurityを同一視しない。
核種のpurity fit結果は別の高R研究成果物に示す。

## 再現と保護

元ROOT・既存低R fit・本番PID係数は変更していない。対象は2021 P24iy、26,956/26,994 subjobsの部分統合のみ。
このレポートは本番投入、欠損job再投入、全23B網羅の確認を含まない。

再現プログラムは`tools/nuclear_pid_cutoff/compare_pid_selections.py`。
保存した`command.json`は、snapshotしたROOTマクロをSL24y／GCC4.8.5／ROOT5.34/38で描画する実command。
既存出力への上書きを拒否するため、再生成には新しい出力ディレクトリを指定する。
