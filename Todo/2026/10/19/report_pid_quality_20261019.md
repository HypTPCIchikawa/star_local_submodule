# 2026-10-19 核PID再検証：旧fit監査と独立studyの進捗

## 現時点の結論

旧設定の狭い除外区間は、純度が50%未満と測定された谷とは限らず、
少数のfit未収束や応答形状の支持失敗により区間全体が未評価となった例が確認できた。
本番PID・係数・本番ROOTは変更していない。新しい元多項式gateの純度認証はこれから行う。

最新追加指示では、背景³He・tの中心・幅は既知のPID多項式から固定参照する。
He4 fit内で自由に動かさず、必要な係数修正は各核種の実peakを用いた別の校正として行う。
He3頂上と左肩を含む局所範囲、および別位置にあるt実peakとの整合を優先する。

### 重要な区別：最新PDFは現行anaFemtoのPID設定ではない

ユーザーのピークずれの指摘を受け、入口から実効PIDまで読取り照合した。
`peak_anchored_views_001`は、専用`*_logpolyLowRReference20261019.yaml`で
`dedxCalibration.enabled=false`として元の`NuclearIdDeDxVsMom.h`を固定した参照studyである。
一方、現在の標準`main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml`は
`purity50_production_20261017.yaml`を通じて
`config/calibration/dedx_helium_auau3p85_2021P24iy_from040.yaml`を有効にしている。
現行³He・⁴Heは支持範囲内で校正点を線形補間し、d・tは元多項式を使用する。
現行校正は、今回新しく得た滑らかな多項式ではない。
元多項式へ戻す変更は独立studyのみで、本番入口への適用は未実施だった。

R=0.875の中心・幅は以下。単位はkeV/cm。これは関数値の照合であり、正しい校正の認定ではない。

| Species | PDF mean | PDF sigma | Current mean | Current sigma |
|---|---:|---:|---:|---:|
| t | 23.691733 | 1.715000 | 23.691733 | 1.715000 |
| ³He | 33.230633 | 1.903423 | 34.282979 | 1.856006 |
| ⁴He | 46.449511 | 2.679187 | 47.260798 | 2.909743 |

現行設定では、このRの⁴Heは別の`speciesPOverQVeto`の[0.85,0.95)で選択対象外である。
既存部分統合ROOTを生成したjob `0BEFC64E62F19BBA534A129CABC9B106`は、
投入時snapshotの`nuclearid_femtoLambda_legacy.yaml`を使用した。現行校正での再解析ではない。
従って「元ROOT生成時」「現在の本番設定」「今回の比較用設定」を混同してはならない。

R=0.875・Gaussian/exponentialの20–28付近の図はt controlを表示したp.7で、
表示範囲は20–27.5。He領域はp.5の29–57.5であり、共同fit窓は20–57.5である。
tとHeのページの用途をタイトルだけで区別しにくく、説明が不十分だった。
今回は確認結果のみ記録し、解析コード・設定・PDF・元ROOTは変更していない。

### Gitに保存された元PIDから何が変わったか

ローカルHEAD `8f7d539`と作業中ファイルを照合した。保存済み`origin/main`参照も同じcommitだが、
今回はfetch/pullをしていないので、リモートの最新状態を照会したという意味ではない。

**元の多項式係数は不変だが、有効なPIDはここ数日の作業で切り替わっている。**
`NuclearIdDeDxVsMom.h`の現在ファイルとHEADのGit blobは
`7f61dffb53d02742c81c59dad951af46dfdd1730`で完全一致する。
このファイルの最終コミットは2026-06-15の`f5ca317`である。
HEADの`FemtoLambdaLegacy.cxx`は、この元関数の`GetNSigma`を直接使用していた。
HEADの⁴He highpurity mainconfも`nuclearid_femtoLambda_legacy.yaml`を参照していた。
現在は未コミットのprovider導入と参照変更により、³He・⁴Heの中心・幅へ校正YAMLを優先適用する。
d・t自身の曲線が同じでも、競合するHeのpullやcutoff等の変更で最終選択は変わり得る。

| 作業記録上の日付 | 本番入口への変更 |
|---|---|
| 10/13 | ³He・⁴Heの中心/幅を線形補間する校正を導入し、⁴He系3入口へ接続。 |
| 10/17 | canonical⁴Heと⁴He/d/t/³He highpurityの5入口を共通校正・cutoffへ接続。 |
| 10/19 | 元多項式固定の独立study。本番を元PIDへ戻す変更は未実施。 |

根拠は[10/13の実装記録](../13/task3_helium_pid_cf_implementation_report.md)と
[10/17の接続記録](../17/report_highrigidity_purity_20261017.md)。
ファイルのmtimeだけで厳密な適用日時を推定してはいない。
新providerはuntracked、校正/PID YAMLは現在の`config/*` ignore規則によりGit管理外であり、
今回確認した校正の導入・有効化はHEADには未コミットである。

元関数と今回PDFの参照は一致しているため、そのPDFのピークずれ自体を
「最近の新校正で悪化した証拠」とは言えない。新校正を良いとも認定しない。
別設定を混ぜたfitの比較ではなく、同じ生分布上でGit元曲線と現行実効曲線の中心・幅を
直接比較することが必要である。今回は読取り監査とこの記録のみで、コード・設定・PDFは変更していない。
`_prePurity20261017`は10/13校正を残す切り戻しであり、元PIDへ戻す`_pidLegacy`とは区別する。
checkout/reset/rollbackは実行していない。

## 完了した検証

- 設定監査9試験、品質監査6試験、描画7試験、独立入力16試験の計38試験PASS。
  pilot runner9試験も通り、合計47試験PASS。
  core＋tail数学試験はROOT5で5489項目PASS（math_selftest_002）。
- 元の`log10(R)`多項式について、係数固定、正の幅、独立評価との算術一致160016行を検査。
  5入口のstudy／切り戻しと本番YAML不変更を確認。これはC++runtimeやpurityの保証ではない。
- 旧fit100実行・14840行を再集計し、独立gate count640行・gate SHA・旧集約の数値と理由を再現。
- 旧veto区間と前後の診断PDF68頁をROOT5で生成し、全頁読み戻し・入力SHA不変を確認。
  不良fitも表示。広域／局所、linear／logを分け、縦軸は表示内のcount+errorとモデルから設定。
  tの局所図では巨大なdピークを除き、モデル線は保存fit窓内だけに表示する。

## 旧vetoの監査例

以下の「不支持」は従来の20条件比較によるもの。新しいχ²条件を適用した認証結果ではない。

| 核種 | 旧除外R区間 [GeV/c] | 不支持の旧条件数 | 主な旧理由 |
|---|---|---:|---|
| d | [1.45, 1.50) | 1 / 20 | 数値未収束・対象shape／帰属screen |
| t | [1.35, 1.40) | 2 / 20 | 数値未収束・中心／幅の誤差・shape／帰属screen |
| ⁴He | [0.85, 0.90) | 3 / 20 | 数値未収束・中心／幅の誤差・shape／帰属screen |
| ⁴He | [0.90, 0.95) | 1 / 20 | 数値未収束・shape／帰属screen |
| ⁴He | [1.75, 1.80) | 1 / 20 | 数値未収束・中心／幅の誤差・shape／帰属screen |

これらは「除外不要」と確定したわけでもない。旧fitの前提が不十分なため、
良いモデルによる再評価が必要である。dの[1.70,2.00)には残差不良もあり、別に調べる。
この区間はTPC-only上端に達し、同じ試料条件の右隣はない。R≥2のTOF試料を同条件の隣とは扱わない。

## χ²から分かった問題

旧低R nominalのPearson χ²/conditionalNdf中央値は、d約7559、t約812、³He約535、⁴He約83。
これは対象局所fit窓内の値であり、表示している全dE/dx範囲を誤って加算した値ではない。
非常に高統計なので小さな形状不一致も強く検出される。旧「支持」は数値・帰属・gate残差等の
screenであり、良好なχ²を保証していなかった。

χ²上限1.5／2／3／5／10／20を比較すると、低Rでは必要な応答・背景familyの被覆が残らない。
したがって、悪いfitを落としても認証された低R purityは得られていない。
これは低purityの測定や、低R候補を全部本番から除く決定ではない。
χ²閾値だけを極端に緩めて認証することはしない。

## 追加回答と次の比較

ユーザーからGaussian core＋tail応答を局所fitへ追加して比較する許可を得た。
まず同じ保存ALL・旧fit窓・旧gateを使い、d R=1.475、t 1.375、⁴He 0.895／1.775の4点で比較する。
Gaussian control、対称tail、正／負のoffsetを持つtailと、linear／exponential背景を比較する。
対象だけのshapeを動かす場合と、局所の競合核種も動かす場合を分ける。

tail割合・幅比・初期offsetは最初のpilotでは固定profileとし、核種yieldは1列のまま扱う。
tail offsetはseed幅から設定したdE/dx差であり、fit後の幅に追従する量ではない。
core中心・幅は混合分布の平均・RMSや本番PID係数ではない。
profile走査後の最良χ²を、独立検定のp値や独立モデル数とは称さない。

pilotで応答と競合核種の扱いを確認した後に、元多項式の新gateへ進む。
旧gate countは新gateへ転用せず、新しい独立countとbest-species境界を評価する。
R<2はALL、R≥2は全核種でTOF m²を要求し、⁴Heはdと同じ窓でR<5まで調べる。
従来のTPC-only上限2.5をTOF併用時へ持ち込まない。

2026-10-19時点でd R=1.475のtailなしcontrolは旧fitのχ²・尤度・ndf・statusを完全再現。
続く48実行64fitの限定pilotはtmux `nuclear-pid-coretail-20261019-001`で開始した。
全48実行64fitが完了し、Gaussian control8点すべて旧χ²・尤度・ndf・statusを再現した。
対応する比較図12冊64頁も全頁readbackと入力SHA不変検査を通過した。
描画13試験を含め、Python試験は計60件PASS。

固定tailの例として⁴He R=1.775・exponential背景は、対象だけのGaussianでχ²/ndf=93.019、
競合核種もfreeのGaussianで12.067、同じfree条件のpositive tailで8.765となった。
93→8.8の全改善をtailだけの効果とは解釈せず、nuisanceと応答変更を分ける。
d/t等では依然大きなχ²が残り、固定profileだけでは十分なモデルと判断しない。
参考purityが高いことを採用理由にせず、次はtail自体の割合・中心・幅をfitして比較する。

### 追加指示による方針変更：自由な背景shapeではなくPID曲線との整合

上のfree比較案は追加指示で保留した。数学試験256項目はPASS、入力24実行32fitは作成したが、
実データfitは未実行。`superseded.json`を保存し、誤起動はrunnerで拒否する。

添付図の細い橙赤線はtモデルの右tailだった。R=0.875のlog-normal＋linearで、
t中心22.686 keV/cmに対してfit範囲35.5–59 keV/cm、窓内acceptanceは約3.0×10⁻¹⁰。
窓内だけで約10.9万の収量を自由にfitした結果、全線外挿は3.63×10¹⁴となり、
同じALL投影の6.94×10⁸計数の約52万倍を要求する。これはt実peakの測定ではなく、
極小tailの局所正規化が背景basisとして使われた病理であり、物理的混入やpurityに転用しない。

新しい初期比較は、既存のlog10(R)中心3次・幅2次多項式で全核種のshapeを固定する。
同じALLのt・³He・⁴He peakを含む連続した局所fit範囲を使い、t規格化を実peakで拘束する。
³Heの頂上と左肩も含める。ROI外のdを自由な遠方tailとして追加せず、寄与上限を診断する。
全域0–120 keV/cmの形をfitする方針ではない。各領域の残差・曲線と実peakのずれを表示する。
曲線がずれた場合、He4 fit内で背景中心を動かして隠さず、別途核種校正で多項式を更新する。
これは初期の整合性検証であり、純度認証・新cutoff採用・本番切替はまだ行わない。

### 元PID関数に固定した実peak同時fitの結果

`peak_anchored_pilot_001`の4実行28fitを、STAR SL24y／ROOT5で完了した。
R=0.825、0.875、0.925、0.975、1.725、1.775、1.825 GeV/cの7点について、
Gaussian／log-normalとlinear／exponential背景を比較した。
全核種の中心・幅は元の`log10(R)`多項式に固定し、R bin内の3 nodeでも同じ関数を参照する。
同じALLデータのt実peakから⁴He右裾までを連続した局所範囲でfitした。

全28件で数値収束した。閉じたROOTのraw count/error/flow、成分和、診断mask、残差、
入力SHAの不変を検査した。関連Python試験はprep8件・描画17件の計25件PASS。
要求R値とnative bin中心の最終桁の丸め差を修正したが、入力bin・窓・countは変えていない。

t規格化は全28件でALL軸内countの記述的margin以内だった。
R=0.875・log-normal/linearでは、t全線yieldは約1.288×10⁷、ALL軸内countは6.942×10⁸。
旧fitの3.63×10¹⁴という遠方tailの異常な外挿は発生しなかった。
ただし、これを物理的t混入量の測定と解釈してはならない。

重要なのは、**元PID曲線に固定するだけでは良いfitにならなかった**ことである。
joint Pearson χ²/conditionalNdfは全28件で約1887–6928と大きく、図にもpeakのずれが残る。
例えばR=0.875の元曲線は³He中心33.231、幅1.903、⁴He中心46.450、幅2.679 keV/cmで、
データの³He頂上とのずれを目視できる。これは単なる不収束点の飛び地とは異なる。
この28件はpurityの系統差・cutoff採用の根拠から除外し、`selection.available=false`、
`purity=null`のまま保存した。未評価を低purityや0%に置き換えない。

[最新比較PDF・28頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/peak_anchored_views_001/peak_anchored_diagnostics.pdf)
では、各RにHe領域2頁、t領域2頁を割り当て、4モデルを比較している。
³He頂上と左肩を含むHe領域、および別位置のtピークをそれぞれ表示する。
全112panelの表示条件と28頁のreadbackを検査し、R=0.875を目視確認した。
再fit・表示規格化はせず、旧PDFも保持した。

次はTAG_T／TAG_HE3の実peak近傍で独立に中心・幅を校正し、元と同じ中心3次・幅2次の
`log10(R)`多項式係数として比較する。He4 fit内で背景中心を動かしてずれを吸収しない。
TAG_D内の⁴He領域は純粋試料ではないので、校正した競合核種を固定して別途検証する。
既存`curves_native_001`はt／³Heの支持域が今回の0.8–1.8を覆わず、⁴Heの支持曲線もない。
過去のglobal／低R polynomial候補も未収束や旧free nuisance条件のため、そのまま再利用しない。
係数の更新・本番への採用はまだ行っていない。

## 成果物

追加指示の反映：004では左linearを全ページで核種近傍に統一。
⁴He図は³Heと⁴Heのピークを含む範囲をRごとに計算し、0–120 keV/cmのlinear図を廃止した。
右logの広域表示は残し、fit/countは不変更。旧003は保持。
全68頁readback・全linearのfocus方針検査・表示11試験PASS。

- [旧fit品質監査](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/audit_001/README.md)
- [dの最新診断図・28頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/representative_views_004/fit_diagnostics_d.pdf)
- [tの最新診断図・12頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/representative_views_004/fit_diagnostics_t.pdf)
- [⁴Heの最新診断図・28頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/representative_views_004/fit_diagnostics_he4.pdf)
- [独立study入力](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_coretail_20261019/study_inputs_001/preparation.json)
- [再現手順](../../../../tools/nuclear_pid_cutoff/quality_20261019/README.md)
- [core＋tail比較・再現手順](../../../../tools/nuclear_pid_cutoff/coretail_20261019/README.md)
- [最新：元PID関数固定・実peak同時fitの再現手順](../../../../tools/nuclear_pid_cutoff/peak_anchored_20261019/README.md)
- [最新：元PID関数固定・He／t比較28頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/peak_anchored_views_001/peak_anchored_diagnostics.pdf)
- [⁴He・exponential・競合核種もfitする比較](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/coretail_pilot_views_001/he4_exponential_local_neighbors.pdf)
- [t・linear・対象のみfitする比較](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019/coretail_pilot_views_001/t_linear_target_only.pdf)

旧PDF001／002、失敗試行のログも保持。本番job投入・実データhadd・メール・Git commit/pushは行っていない。
