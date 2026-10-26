# Plan: KFParticle Lambda purity study — 2026-09-27

## Scope and fixed reference

今回はKF導入の有無を比較せず、KF Imp5を基準にpurityを高める追加・強化cutを探索する。既存解析・default YAML・旧出力は変更しない。先頭10,000入力 / 9,744再構成イベントを固定し、Λのみ、既存FXT/event/centrality条件を維持する。

ユーザーは許容する統計減少率を先に固定せず、結果を見て決める方針。従って単一の推奨設定へ自動変更せず、purity対signal retentionの比較と複数の候補を提示する。

Reference: [previous summary](../../../analysisnote/20260914/summary20260914.md), [four-curve table](../comparisons/lambda_four_histogram_cut_summary_20260914.md), [reference-note / qhu1 comparison](../comparisons/kfparticle_lambda_qhu1_cut_comparison_20260908.md).

## Estimators

- Primary mass window fixed: [1.110, 1.122) GeV/c².
- Default sidebands: [1.098, 1.106), [1.126, 1.134); B = 0.75 × sideband count, S = Nwin − B.
- Purity = S/(S+B) = S/Nwin. S/Bは別指標として併記する。
- Relative signal retention = S/Sbaseline。背景込みNwinの比、真の再構成効率とは区別する。
- Near/far sidebandsとmass fitをcross-checkし、差を背景モデル依存として記録する。mass窓・binを狭くしてpurityを上げたことにはしない。
- Baseline: completed saved KF Imp5 result, Nwin=1370, B=143.25, S=1226.75, sideband purity≈89.54%. Source identity and selected histogram replay must pass before scanning.

## Candidate cuts and reference limits

保存済み候補から、fit/topology χ²、崩壊長・有意度、daughter/PV距離、pointing、TPC nSigma、元track gDCAの追加cutを走査する。既存選別より厳しいcutに限定し、upstreamで捨てた候補を回復できるとはしない。

Reference PDF §2.3/p16・§3/p33は、daughter primary χ²10（変動20）、hits15（変動20）、qhu1現存実装のFinder L>5 cm、L/σL>5を追加候補の根拠とする。これらのupstream変数は旧treeに未保存なので、PicoDstが取得できれば独立したmainconfで再処理する。Scalar final cutsとSIMD upstream cuts、original daughterとparent-refitted daughterのχ²は別の変数として扱う。

TOFはImp5 tree内では未評価placeholderなので、このtreeから有効性を判定しない。pT/η/centralityのacceptance変更、massErrorによる選別、parent mass constraintは今回の初期gridに含めない。参考図のcollider7.7 GeV/20–50%と今回のFXT13.5 GeV beam/√sNN≈5.2の差を維持し、η・PV条件はコピーしない。

## Scan and validation design

1. 保存済み完成ROOTをREADで開き、既存structural/Imp5 QAを実行。全候補・metadata・元ヒストグラムを新規study領域へexport。
2. Gridは[QA YAML](../../../config/qa/lambda_kf_purity_study_20260927.yaml)に記録。single/pairと事前に指定した少数のphysics tripleを調べ、全次元の無制限gridは行わない。
3. runId/eventIdによる再現可能なhashでevent単位のtrain/validationへ分割。同じeventの候補を両方に分散させない。
4. Trainだけでsignal retention 95/90/80/70%以上の代表条件を選ぶ。Validationと全10,000イベントをその後に評価し、validationの結果でcutを選び直さない。選択点ごとにpaired event bootstrapでpurity差・signal retentionの統計的な幅を調べる。
5. Nwin、背景、S、purity、S/B、signal retention、cut差分を併記する。複雑な条件数と探索の統計的偏りに注意し、単純で安定した領域を優先して説明する。
6. 代表点を独立ROOT計算で照合。取得可能ならPicoDst上でも再実行し、最終YAMLとmetadataを保存する。取得不能なupstream条件を実施済みとは記さない。
7. 図は1.10–1.14 GeV/c²、左counts / 右全範囲[1.05,1.25) Normalize。purity対signal retentionも図示し、異なる保持率の選択をユーザーが比較できるようにする。

## Reproducibility and limitations

- Runtime follows singularity-local-build-run: STAR SL24y / ROOT 5.34/38 / GCC 4.8.5.
- Outputs: rootfile/auau13p5_anaLambda_KFParticle/purity_study_20260927/; figures: share/figure/auau13p5_Lambda_KFParticle_purity_20260927/.
- Existing .current_mainconf / default configuration / outputs remain intact; no farm jobs, commit or push.
- 2026-09-27開始時点で前回の/tmp PicoDst cacheは消失。保存済み10000event候補は存在しQA成功。元URIからの再取得可否を確認中。
- 10,000イベント内の探索結果を、普遍的な最適cut・真の効率測定・BG freeの証明とは呼ばない。独立データでの確定は次段階。
