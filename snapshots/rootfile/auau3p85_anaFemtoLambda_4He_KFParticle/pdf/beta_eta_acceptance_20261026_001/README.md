# Λ・⁴He：β vs η（2026-10-26）

[作図PDF：6頁](beta_vs_eta_lambda_he4.pdf)／[主図PNG](overview.png)／[低R TOF比較PNG](he4_low_tof_comparison.png)／[合算TH2 ROOT](beta_eta_summed.root)。

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

全40,804 bin（flow込み）で、lowR_noTof+highR_tofM2=inclusive、lowR_tofM2⊆lowR_tofMatch⊆lowR_noTofが成立し、不一致0。全6試料でentries=全bin合計、有限・非負の候補数を確認。入力/再構成イベント数は承認済み試料と一致。元ROOTのfingerprintは作図前後で不変。当初5頁PDFはGhostscriptで全頁読取りを確認し、主図・低R比較を目視確認した。

## 再現用ファイル

- `PlotBetaEta.C`：ROOT5で対象TH2の合算・検査・描画。
- `plot.env`：入力リスト、図のラベル・余白等。
- `inputs.list`、`metadata.json`：3入力の実パス、fingerprint、選択条件。
- `counts.tsv`、`validation.tsv`、`completed.json`：候補数、flow、closure、完了・hash記録。
- `run.log`、`pdf_validation.log`：batchと同じSL24y/ROOT5実行、PDF検査。

プロジェクト直下から、`script/singularity_root_batch.sh`へ通常mainconfを指定し、root4starで`PlotBetaEta.C+("plot.envの絶対パス")`を実行する。これは保存済みヒストの後処理で、PicoDstイベント再解析・job再投入・hadd再実行は行わない。


## 改訂：1頁の文字サイズとBOX重ね描き

ユーザー指定により1頁の軸タイトル/ラベルをfont43のpixel指定（22px/17px、panel title24px）へ統一し、軸タイトルを中央配置した。ηの整数目盛と余白も調整した。元のη[-5,0.15]・β[0,1]、候補数・binning・selectionは維持。

[BOX図PDF](beta_eta_box_overlay.pdf)／[BOX図PNG](beta_eta_box_overlay.png)／[BOX図ROOT（3 TH2D＋canvas）](beta_eta_box_overlay.root)。Λ=青、低R⁴He=橙、高R⁴He=緑。ROOT描画optionは`BOX SAME`、fill style0の輪郭boxで重なりを見せる。共通の最大値202525.05でraw countの線形尺度を使い、試料ごとの規格化はしない。1頁のCOLZは従来通りlog-zである。

現行PDFは6頁。1頁を修正版へ置換し、旧PDFの2〜5頁を再利用（GSでPDF再出力）、6頁へBOX図を追加した。元5頁版は`revisions/20261026_before_box_overlay_001/`に保持。

`RenderBetaEta.C`＋`render_layout.env`は保存済み合算ROOTから1頁とBOX図を再描画する。`ValidateBoxOverlay.C`で新ROOTを再読取りし、3 TH2Dの全122,412 bin（flow込み）の内容・誤差、Entries、表示範囲、共通maxが合算ROOT/指定値と一致することを検証（不一致0）。色とcanvasも再読取りできる。元のbeta_eta_summed.rootはSHA256不変。現行hash/検証はcompleted.json、overlay_validation.tsv、render_layout.logへ記録。
