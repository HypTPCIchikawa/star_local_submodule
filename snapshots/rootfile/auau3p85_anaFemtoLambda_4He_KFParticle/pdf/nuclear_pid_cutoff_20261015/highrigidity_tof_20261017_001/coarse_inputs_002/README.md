# 高R研究用の厳密R集約入力

2026-10-17作成。3.85 GeV・2021 P24iyの検証済みnative cacheだけを入力とし、
隣接するR=|p/q|の5 native binを非重複に加算した。
これは高R・低統計の⁴Heを別stageで直接再fitするための入力で、fitやpurityの平均化ではない。
元の約10 GB ROOTは記録済みfingerprintをstat照合しただけで、開いたりSHA計算したりしていない。

## 入出力

| 項目 | 元のnative cache | 今回の集約cache |
|---|---|---|
| R範囲 | 0–6 GeV/c | 0–6 GeV/c |
| R分割 | 600 bin・幅0.01 GeV/c | 120 bin・幅0.05 GeV/c |
| dE/dx | 1200 bin・幅0.1 keV/cm・0–120 keV/cm | 変更なし |
| キー | ALL / TAG_D / TAG_T / TAG_HE3 / REMAINING | 同じ5キー |
| R境界 | 元のnative境界 | 元の5 binごとの境界を厳密継承 |
| flow | 各軸のunderflow / overflowと隅のflow | 全て保持・対応するRグループ内だけ加算 |
| Entries | 元のGetEntries | 各ヒストでそのまま保持 |
| variance | 元の保存Sumw2 | 保存arrayを直接加算し、誤差の平方から再構成しない |

全5ヒストで、入力3,618,020 cellから出力733,220 cell（146,644 cell/ヒスト、flow込み）を検証した。
別アルゴリズムのtarget側直接和とcount・保存varianceが全cellで厳密一致した。
閉じた出力ROOTを開き直し、全cellのcount・Sumw2・GetBinError・軸・Entriesの一致も確認した。

独立toyは非整数count、countと異なるvariance、ゼロcount、全flow、2通りの集約率、
端数が残る不正な集約率の拒否を含み、611チェックを通過した。
全実行は`COARSE_INPUTS_COMPLETE ... checks=20532676`で正常終了した。

## 計数の保存

| キー | 全保存count（flow込み） | Entries | flow count |
|---|---:|---:|---:|
| ALL | 316,515,592,807 | 316,515,591,143 | 381,505,788 |
| TAG_D | 3,555,393,696 | 3,555,393,696 | 322,296 |
| TAG_T | 528,861,766 | 528,861,766 | 891,161 |
| TAG_HE3 | 252,915,567 | 252,915,567 | 5,276 |
| REMAINING | 312,178,421,778 | 312,178,421,778 | 380,287,055 |

各ヒストのvariance総和も元のnative cacheと一致し、今回の総和比較で相対roundoffは0だった。
ALLのcountとEntriesの差1664は元TH2Fの保存丸めに由来する既存差で、Entriesへの再正規化はしていない。
TAGはALLの部分集合であり、独立な追加尤度として足し合わせてはいけない。
REMAININGは元cacheの値とvarianceをそのまま集約し、ALLとTAGの誤差を再加算していない。

## ファイルと来歴

- [cache.root](cache.root): 新しいR幅0.05の5ヒスト。1,571,842 bytes。
- [native_audit.json](native_audit.json): inputNativeAxisとcoarseAxis、全cell直接和・読戻し・計数監査。
- [context.json](context.json): 原ROOT fingerprint、元native cacheのSHA、TOF m²窓、集約率、明示runtime。
- [completed.json](completed.json): 実行完了、前後SHA/fingerprint不変、各出力・sourceのSHA。
- [root.log](root.log): ROOT5/STARでのtoy・全cell監査ログ。
- `provenance/`: 元native入力のcontext・audit・completedをコピーした不変の来歴。
- `source/`: 今回のC++/Python・mainconf・analysis_info・runtime wrapper snapshot。

`native_audit.json`という名前は後段との入出力契約に従う。
内容の`processingStage=exact_r_coarsening`および`coarseAxis`を必ず確認し、
幅0.05を元のnative幅0.01と取り違えないこと。
sourceのdタグm²窓は3.268773–3.693147 (GeV/c²)²で、⁴Heにも互換性があるが純⁴He試料ではない。
3–4の広い⁴He窓を再現した入力ではない。

最初の[coarse_inputs_001](../coarse_inputs_001/root.log)はROOT5 accessor非対応によるcompile失敗記録として保持した。
同値のbin総数式へ修正し、この002で成功した。001には成功cacheやcompletedはない。

## 再現コマンド

リポジトリ直下から、未使用の出力名へ実行する。

```bash
python3 -B tools/nuclear_pid_cutoff/highrigidity/prepare_coarse_inputs.py \
  --input-study rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_purity50_20261015.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/highrigidity_tof_20261017_001/coarse_inputs_NEW \
  --factor 5 --execution-timeout 1200
```

元cache・既存fit・本番PIDは変更しない。後段のcoarse fitは別source・別出力で管理する。
