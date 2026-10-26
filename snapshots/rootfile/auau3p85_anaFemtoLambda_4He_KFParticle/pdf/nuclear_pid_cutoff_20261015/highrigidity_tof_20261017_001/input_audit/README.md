# 高R・TOF選択後purityに使う保存ヒストグラムの監査

2026年10月17日（日本時間）。対象は **Au+Au FXT 3.85 GeV・2021 P24iyのみ**。
この記録は入力の選択条件・来歴の監査であり、新しいpurity fitや本番設定の更新ではない。

## 結論

| 核種 | 元ROOTのpre-TPCヒストグラム | cacheキー | 保存されたstrict `(m/q)^2` 窓 | 現行高R設定との対応 |
|---|---|---|---|---|
| d | `true/hDedxP_d_m2` | `TAG_D` | `3.268773 < m2 < 3.693147` | 一致。`2 <= R < 2.5`で現TPC選択後の推定に使用可能 |
| t | `true/hDedxP_t_m2` | `TAG_T` | `7.2564275 < m2 < 8.2816925` | 一致。`2 <= R < 2.5`で使用可能 |
| ³He | `true/hDedxP_3He_m2` | `TAG_HE3` | `1.78285 < m2 < 2.06485` | 一致。`2 <= R < 2.5`で使用可能 |
| ⁴He | 専用ヒストグラムなし | 該当なし | 現行の`3 < m2 < 4`は未保存 | このROOTだけでは現行窓のpurityを確定できない |

現行設定は `R >= 2` でTOFを要求する。⁴Heに限り`maxPOverQ=2.5`より高いRへの延長を認め、
d・t・³Heには`R < 2.5`上限を維持している。保存軸の`R < 6`とoverflowも区別する。

`TAG_D`はdと⁴Heがともに入り得るTOF窓であり、純粋なd/⁴Heの真値標本ではない。
⁴Heについてこの狭い窓の参考purityを求めても、**現行の広い`3 < m2 < 4`のpurityに置き換えてはならない**。

## 実際に照合した根拠

1. 元ROOTの全282ヒストグラムの既存監査表`merged_audit.tsv`。
   `hDedxP_4He_TOFControl`、`hNSigmaHe4VsP_TOFControl`、`hM2P_4He_TPC`、
   `hDedxP_4He_selected`、`hM2P_4He_selected`はなく、TH3もない。
   現行ソースにこれらの定義があることは、2026年10月1日の元ROOTに保存された証拠ではない。
   元の`true/hM2P`はR–m²だけなので、dE/dxとの同時選択を復元できない。
2. 提出時runtime bundleのSHA256を提出runmetaと照合した上で、
   当時の`FemtoLambdaLegacy.cxx/.h`、核種cut YAML、核種hist YAMLを抽出した。
   現行作業ツリーのC++を過去生成時のソースとして扱っていない。
3. 当時の`CollectNuclei`では、品質選択`nHitsDedx >= 15`、`gPt >= 0.1`、`dEdx > 0`後、
   `R = pMom().Mag()`とTOFの`m2 = R^2(1/beta^2 - 1)`を使用している。
   TOFには`isTofTrack`、正常なtraits index、`beta > 0`を要求する。
   QA充填はTPC best-speciesや最終R範囲に依存せず、その選択で`continue`しない。
4. 元YAMLのTOF中心±`1.5 sigma`から窓を再計算し、既存cacheの`context.json`と、
   `nuclearid_femtoLambda_auau3p85_purity50_20261015.yaml`の窓に一致することを数値検査した。
   TOF窓はtrack単位で適用済みなので、m²方向のヒストグラムbin境界近似ではない。
5. 小さい既存native cache、`context.json`、`native_audit.json`のSHA256を、
   既存`completed.json`の記録と照合した。元10.5 GB ROOTはstat fingerprintのみ再確認し、
   ファイル全読取り・全体SHA256再計算はしていない。

## 条件付きpurityとしての注意

- native軸は`R = 0–6`、幅`0.01 GeV/c`、`dEdx = 0–120`、幅`0.1 keV/cm`。
  現TPC nσ/best-species窓を再積分する際には有限binの近似を明示する必要がある。
- 元のTPC核種係数と現在の校正が異なっていても、pre-TPC QAから現ゲートを評価できる。
  元`hDedxP_4He`などのTPC選択済み標本を、そのまま現選択後の標本と扱ってはならない。
- TOFタグはALLのsubsetで、独立な追加観測としてALLとのlikelihoodへ二重加算しない。
- charge-inclusive・centrality積分のQAからの推定であり、CFペア、centrality、mT毎のpurityではない。
- 入力はSUMS `0BEFC64E62F19BBA534A129CABC9B106`の26,956/26,994 subjobsの部分統合。
  欠損38件は未収録。2019年は使用せず、全23 billion eventsを網羅したとの認定でもない。

## 保存内容と再現

- [audit.json](audit.json)：条件、窓照合、元ROOT fingerprint、キー一覧、限界の機械可読記録。
- [snapshot/](snapshot/)：過去runtimeから抽出した該当ソース/YAMLと、現在・既存cacheの証拠資料。
- [SHA256SUMS](SHA256SUMS)：snapshotのハッシュ。
- [completed.json](completed.json)：全監査成果物のハッシュと小さい元入力のハッシュ。
- [replay_audit.py](replay_audit.py)：非上書き再生成・検証用プログラム。

リポジトリ直下で検証する例（既存成果物を書き換えない）：

```bash
python3 rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/highrigidity_tof_20261017_001/input_audit/replay_audit.py \
  --project "$PWD" \
  --output rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/highrigidity_tof_20261017_001/input_audit \
  --verify
```

再生成は`--verify`を外し、`--output`に未作成の別ディレクトリを指定する。
既存audit・snapshotへの上書きは拒否する。再現処理もfit・元ROOT全体の読込みは行わない。
