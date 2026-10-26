# Plan: two-energy KF Lambda upstream-cut study — 2026-09-28

## Request and preservation

PicoDst読み出しが復旧を待たずNFS経由で可能になったため、AuAu13p5とAuAu3p85について、保存treeでは変更できなかった上流cutも再構成で調べる。各10,000入力イベントで条件を選び、各100,000入力イベントでdefaultと固定したbestを比較する。収量約半分を許容する方針を背景差引き信号 `S/Sdefault >= 0.50` として評価し、13p5の前回70%条件も比較に残す。背景込み候補数比は別に報告する。

既存Helix/KF Maker、標準YAML、元ROOT、`.current_mainconf` は変更しない。新規study設定・QA補助コード・出力だけを追加する。farm submit、commit/pushは行わない。実行は既存StChain/Maker経路のSTAR SL24y / ROOT5.34/38。探索結果はproduction defaultへ自動反映しない。

## Inputs and event policies

- 13p5: `production_13p5GeV_fixedTarget_2020`, P24iy, NFS `/star/data22/...`。探索は以前の先頭7ファイルから同じ10,000入力。保存SHA256との一致確認済み。
- 3p85: `production_3p85GeV_fixedTarget_2019`, P24iy, NFS `/star/data24/...`。探索はrun20160024の `st_physics_20160024_raw_4000004.picoDst.root` の先頭10,000入力。ROOT内17,941entries、先頭/中央/末尾の共分散読み出し確認済み。
- 100,000イベント比較には探索で使うファイルを丸ごと除いた別リストを確保する。default/bestには同じ入力順を使用。event headerから全run/event IDを記録し、重複・探索標本との非重複・正確な入力件数を検証する。
- 3p85 baselineは既存HelixのImp5共通cutと同じ数値に標準KF上流条件を加えたもの。13p5 baselineと同じ `lambda_imp5` / `pico_nsigma` / Finder+Topoを用いるが、専用analysis_infoに2019のproduction・run範囲・FXT modeを記録する。
- イベント・centrality条件はenergy内で固定する。既存Imp5は明示的なvertex radius/Vz cutを使わず、maxNTrとcentrality選別を使う。QA centerをvertex acceptance cutと混同しない。
- 既存StRefMultCorrの広いFXT run-range extensionとcalibration上の制約を引き継ぐ。`setParameterIndex` は最初の一致を採用するため、広げたindex0の範囲に他energyが入る点も記録する。今回centrality再較正や既存解析の修正を行わず、積分purity比較のbaseline/event選別を固定する。centrality百分率の精密測定を主張しない。

## Predeclared reconstruction grid

各energyで以下の16条件をPicoDstから10,000イベントずつ再構成する。未記載の値は標準KF Imp5を維持し、geometryを無効化しない。

| Variant | Upstream overrides |
|---|---|
| default | None |
| hits20 | minNHitsFit = 20 |
| hits25 | minNHitsFit = 25 |
| ratio060 | minNHitsRatio = 0.60 |
| primary5 | interfaceChiPrimaryCut = 5; finderChiPrimary2D = 5 |
| primary10 | interfaceChiPrimaryCut = 10; finderChiPrimary2D = 10 |
| primary20 | interfaceChiPrimaryCut = 20; finderChiPrimary2D = 20 |
| finderL3 | finderLCut = 3 cm |
| finderL5 | finderLCut = 5 cm |
| finderLsig5 | finderLdL2D = 5 |
| finderLsig7 | finderLdL2D = 7 |
| finderFit5 | finderChi2Ndf2D = 5 |
| topo1p5 | topoChi2NdfCut = 1.5 |
| combined_soft | minNHitsFit = 20; both primary thresholds = 5; finderLCut = 3; finderLdL2D = 5 |
| combined_tight | minNHitsFit = 20; both primary thresholds = 10; finderLCut = 5; finderLdL2D = 5 |
| detached_combination | both primary thresholds = 10; finderLCut = 5; finderLdL2D = 5 |

Interfaceのscalar primary分類、Finderのtransport後SIMD primary偏差、raw vertex-line有意度、PV-constrained final崩壊長有意度は異なる量。Finderのstrict不等号とMakerのinclusive不等号も区別する。hit ratioは元コードどおりnHitsMax>0時のみ。

今回の未検証項目: TOF PID（Imp5で禁止されており、別profileへ切替えると他のcutも変わる）、qhu1のparent-refitted daughter χ²（現在未実装）、covariance安全性の緩和、pT/eta/centrality受容域最適化、mass constraint。これらを実施したとは報告しない。上流再構成が必要だったhits・original daughter primary条件・Finder条件を今回の対象とする。

## Final-cut scan and fixed ranking

各upstream出力の `selected == 1` を出発点に、保存済み量への追加cutのsingle/pairとfinal decay-length significanceを含むtripleを走査する。gridは新規QA YAMLへ事前保存し、1出力あたり最大6,000点。nSigma p/pi、original daughter gDCA、final parent/daughter距離、pointing、final fit/topology χ²、崩壊長と有意度を含む。信号窓やmass errorは最適化しない。

| Quantity | Definition |
|---|---|
| Signal window | [1.110, 1.122) GeV/c2 |
| Default sidebands | [1.098, 1.106), [1.126, 1.134) GeV/c2 |
| Near sidebands | [1.100, 1.108), [1.124, 1.132) GeV/c2 |
| Far sidebands | [1.094, 1.102), [1.130, 1.138) GeV/c2 |
| B; S; Purity | 0.75*(Nleft+Nright); Nwin-B; S/Nwin |
| Retention denominator | Same-energy default reconstruction, same sideband convention |
| Eligibility | Nwin >= 20; positive finite S; minimum S/Sdefault over three sidebands >= 0.50 |
| Ranking | Maximum worst-sideband [1 - 0.75*PoissonUpper90(Nsideband)/Nwin] |
| Tie breaking | Fewer changed upstream plus downstream settings, larger S, fixed point identity |
| Additional comparisons | 70% and 80% signal-retention choices; nominal-purity best at 50% |

0背景でも上限は正なのでpurity100%を誤差0として優遇しない。順位付け量は厳密なpurity信頼下限ではない。探索は10,000全体を使う探索的なものとして扱い、以前の13p5 splitを独立holdoutとは呼ばない。chosen条件をevent単位の共通weight bootstrap（400回）と非負sideband backgroundモデルで確認し、モデルを見て100k側のcutを再最適化しない。

## Frozen-cut verification and 100k comparison

1. 13p5 default10kを以前の出力と比較し、3p85 default10kは新規baselineとして構造・候補・metadata QAを行う。
2. 全upstream runについてrequested/completed/input histogramが10,000、終了status正常、event policy不変であることを確認する。期待入力数への暗黙clampを成功扱いにしない。
3. 追加cutのexport/replayを独立ROOTで照合する。各energyのbest50を全YAML条件として固定し、PicoDstから10kを再処理する。post-KFのPID/gDCA scanと上流設定変更が一致するか検査し、差があれば隠さず記録・対処する。
4. defaultと固定bestを同じ独立100,000入力で再構成する。両energyのbestは100k結果を見る前に固定する。別eventで保持率が50%を下回っても隠さず報告し、100k上で選び直さない。
5. 質量比較は1.10–1.14 GeV/c2、左counts、右[1.05,1.25)の全積分でNormalize。default/best各100kのNwin/B/S/purity/S保持率、背景モデル依存、event数、cut表を記録する。

## Output layout

- `rootfile/auau13p5_anaLambda_KFParticle/purity_upstream_20260928/`
- `rootfile/auau3p85_anaLambda_KFParticle/purity_upstream_20260928/`
- 各配下に `configs/`, `scan10k/`, `validation100k/`, `exports/`, `scan/`, `provenance/input/`, `provenance/logs/`。
- 図: `share/figure/auau{13p5,3p85}_Lambda_KFParticle_purity_upstream_20260928/`。
- 最終報告: `mdfiles/kfparticle/comparisons/lambda_kf_upstream_purity_20260928.md`。

処理量・少数候補・入力障害で未完了となった項目があれば、完了済みと分けて記録する。purityは背景モデル依存のピーク純度であり、primary-Λ purity、真の検出効率、普遍的最適cut、BG-freeの保証とはしない。
