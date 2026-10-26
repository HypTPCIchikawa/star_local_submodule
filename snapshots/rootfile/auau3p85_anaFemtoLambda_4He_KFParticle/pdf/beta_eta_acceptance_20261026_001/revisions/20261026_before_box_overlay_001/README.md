# Λ・⁴He：β vs η（2026-10-26）

[作図PDF：5頁](beta_vs_eta_lambda_he4.pdf)／[主図PNG](overview.png)／[低R TOF比較PNG](he4_low_tof_comparison.png)／[合算TH2 ROOT](beta_eta_summed.root)。

保存済みの検証済み3中間ROOTから対象6 TH2Dをそれぞれ1回だけ読み、Clone/Addした。SUMS 04BF2DAE47164B6F040E3A460E539289、2021 P24iy、正常終了26,955/26,994 subjobs、実入力2,100,726,543イベント・再構成2,067,647,474イベント。巨大Sparseを読み込まずに作図した。

## PDFの内容

1. Λ signal／⁴He 0.61<R<2（TOF要求なし）／⁴He 2≤R<5（有効TOF＋m²）の主図。
2. Λ signalの単独図。
3. ⁴Heの低R・高R比較。
4. 低RのTOF不要／有効TOF／有効TOF＋m²を並列比較。
5. 主図3試料を元のη全軸[-5,5]で表示。

主図のη表示は全6試料の非ゼロregular binの共通範囲に余白をつけた[-5,0.15]。β表示[0,1]、元200×200 binを保持。全図で同じ色スケール（log-z、raw candidate counts/bin）を使用。生成数等で割った効率や規格化した密度ではなく、選択候補の分布である。

## 選択とβの定義

βは運動学的|p|/E、横軸はη_lab。He運動量にはZ=2補正が一度だけ反映されている。Λは信号質量窓1.10708≤M≤1.12484 GeV/c²の候補で1候補1fill。核ペアの有無やペア数で重み付けしない。⁴Heは最終PID後、strict 0.61<R<5。境界R=2はTOF要求側に含める。TOF＋m²はstrict 3.268773<(m/q)²<3.693147、dと同じ窓。

低R noTofはTOFを要求しない試料で、TOFの付いたtrackも含む。低R3試料は重複する部分集合なので合算しない。既存TH2タイトル内セミコロンによる軸名誤解釈は描画用cloneで直し、元ROOT・本番設定は変更しなかった。

## 候補数と検証

| 試料 | 候補数（flow込み） | flow数 |
|---|---:|---:|
| Λ信号候補 | 26,618,096 | 5 |
| ⁴He 全R | 143,635,640 | 0 |
| ⁴He R<2、TOF要求なし | 140,067,476 | 0 |
| ⁴He R<2、有効TOF | 38,372,219 | 0 |
| ⁴He R<2、TOF＋m² | 30,836,103 | 0 |
| ⁴He R≥2、TOF＋m² | 3,568,164 | 0 |

Λのflow5件はNに含まれ、regular-binの色マップには描かれない。図中にもflow数を明記。⁴Heはflow0件。

全40,804 bin（flow込み）で、lowR_noTof+highR_tofM2=inclusive、lowR_tofM2⊆lowR_tofMatch⊆lowR_noTofが成立し、不一致0。全6試料でentries=全bin合計、有限・非負の候補数を確認。入力/再構成イベント数は承認済み試料と一致。元ROOTのfingerprintは作図前後で不変。5頁PDFはGhostscriptで全頁読取りを確認し、主図・低R比較を目視確認した。

## 再現用ファイル

- `PlotBetaEta.C`：ROOT5で対象TH2の合算・検査・描画。
- `plot.env`：入力リスト、図のラベル・余白等。
- `inputs.list`、`metadata.json`：3入力の実パス、fingerprint、選択条件。
- `counts.tsv`、`validation.tsv`、`completed.json`：候補数、flow、closure、完了・hash記録。
- `run.log`、`pdf_validation.log`：batchと同じSL24y/ROOT5実行、PDF検査。

プロジェクト直下から、`script/singularity_root_batch.sh`へ通常mainconfを指定し、root4starで`PlotBetaEta.C+("plot.envの絶対パス")`を実行する。これは保存済みヒストの後処理で、PicoDstイベント再解析・job再投入・hadd再実行は行わない。
