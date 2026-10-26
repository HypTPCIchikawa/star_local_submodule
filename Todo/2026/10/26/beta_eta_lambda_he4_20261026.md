# 2026-10-26：Λ・⁴Heのβ vs η分布

## 作業目的と対象

保存済みの検証済み3中間ROOT（SUMS 04BF2DAE47164B6F040E3A460E539289、2021 P24iy、正常26,955/26,994 subjobs）から、対象TH2Dだけを読み込みClone/Addしてβ–η分布を作る。巨大Sparseは読み込まない。全haddや元ROOT整理の処理は並行して継続し、今回の描画で解析コード・PID設定・イベント選択を変更しない。

主図はΛ signal、⁴He R<2（TOF要求なし）、⁴He R≥2（有効TOF＋dと同じm²窓）を別パネルで表示。補足として⁴He全Rと低Rの有効TOF／TOF+m²も表示する。境界R=2はコードに合わせ高R側へ含める。

## 作業項目・進捗

チェックリストは[当日のTodo](STAR_Todo20261026.txt)へ分離した。以下に定義・結果・検証・改訂履歴を記録する。

## 定義上の注意

X=η_lab、Y=運動学的β=|p|/E。TOFで測ったβではない。He運動量はZ=2を1回補正済み。Λは質量信号窓内候補（1.10708≤M≤1.12484 GeV/c²）、核ペア数に依存せず1候補1fill。⁴Heは最終PID後、strict 0.61<R<5 GeV/c。高RはR≥2かつ有効TOF、strict 3.268773<(m/q)²<3.693147。

低R noTofは「TOFを要求しない」であり、TOFの付いたtrackを除外しない。低R3種類は重複する部分集合なので足さない。候補カウントの分布であり、生成数等で割った検出効率そのものではない。

既存subset histogram titleに括弧内セミコロンがあり、ROOTが軸名を誤解釈している。元ROOTや本番YAMLは触らず、描画用cloneでη/βの軸名を明示する。binningは元の200×200、η[-5,5]、β[0,1]を保持し、全試料共通の表示η範囲と色スケールを使う。

## 保存先

`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/`。入力3本と設定は同ディレクトリmetadata.json、inputs.list、plot.envに保存。状態：作図・整合性検査・PDF/PNG目視確認まで完了。


## 完了結果

完成時刻：2026-10-26T09:47:32.276560+09:00。

- [図PDF（6頁・修正版）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/beta_vs_eta_lambda_he4.pdf)
- [主図PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/overview.png)
- [低R TOF条件比較PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/he4_low_tof_comparison.png)
- [合算TH2 ROOT](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/beta_eta_summed.root)
- [保存内容・再現手順](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/README.md)

| 試料 | 候補数（flow込み） | flow数 |
|---|---:|---:|
| Λ信号候補 | 26,618,096 | 5 |
| ⁴He 全R | 143,635,640 | 0 |
| ⁴He R<2、TOF要求なし | 140,067,476 | 0 |
| ⁴He R<2、有効TOF | 38,372,219 | 0 |
| ⁴He R<2、TOF＋m² | 30,836,103 | 0 |
| ⁴He R≥2、TOF＋m² | 3,568,164 | 0 |

主図1頁はΛ・⁴He低R／高R、4頁は低RのTOF3条件比較。Λはflow5件がNに含まれるため図中に注記し、regular図には描かれない。⁴Heは全試料flow0。

全40,804 bin/試料（flow含む）のpartitionとsubsetを確認し不一致0。元ROOTはfingerprint不変。PDF全5頁のGhostscript読取り・主図/低R比較の目視確認と、独立した定義・集計レビューを通過。表示調整前後でcounts.tsvのhash一致。PDF/図の作成完了は、本番haddの最終検証・メール送信完了とは別の状態として扱う。


## 1頁の文字サイズ修正・3色BOX重ね描き（完了）

ユーザー指定により1頁の軸/タイトルの大きさ・余白を修正し、同じη[-5,0.15]・β[0,1]でΛ、低R⁴He、TOF+m²要求高R⁴Heを3色で重ねる。ROOTの描画optionは`BOX SAME`を使用する。元の候補数・binning・選択・合算TH2は変更しない。元5頁PDF等は出力先revisions/20261026_before_box_overlay_001/へ保存。

軸/ラベル/タイトルはfont43のpixel指定でサイズを揃える。BOXは共通のraw count最大値を使い、試料ごとの規格化は行わない。重なりで隠れないよう色付き輪郭のboxと凡例を用いる。修正版1頁、別添BOX PDF/PNG、ROOT canvas/3 TH2を保存し、main PDFにもBOX頁を追加する。


改訂完了：2026-10-26T10:16:05.412521+09:00。1頁の文字サイズ（axis title22px / labels17px / panel title24px）と中央配置・整数η目盛を修正。η[-5,0.15]・β[0,1]を維持し、旧2〜5頁の内容をPDFへ再利用、6頁へ3色のBOX SAME図を追加した。

- [修正版6頁PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/beta_vs_eta_lambda_he4.pdf)
- [BOX図単独PDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/beta_eta_box_overlay.pdf)
- [BOX図PNG](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/beta_eta_box_overlay.png)
- [BOX図ROOT（3 TH2＋canvas）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/beta_eta_acceptance_20261026_001/beta_eta_box_overlay.root)

BOXは青Λ・橙低R⁴He・緑高R⁴He、共通raw-count線形尺度・輪郭box。旧1頁のCOLZ log色尺度とは別の描画であり、試料ごとの規格化は行わない。合算ROOTのhashは不変。保存後ROOTを再読取りし全122,412 binの内容・誤差²相当、Entriesと描画範囲・共通maxが一致（不一致0）。PDF全6頁＋単独BOXをGSで読取り、PNGとcropboxに沿ったPDF第1頁、ROOT canvas再読取り画像を確認し、独立レイアウトレビューを通過。初版はrevisionsへ保持。図の改訂で本番hadd・解析code・元中間ROOTを変更していない。
