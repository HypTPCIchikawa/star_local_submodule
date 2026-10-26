# qhu1さんのΛ選別と現在のAuAu13p5 KFParticle設定の比較

- 調査日・記録日: 2026-09-08
- 目的: 指定analysis noteのFig.9で背景が小さい理由を、現存コードの有効な選別条件から検討する。
- 対象: qhu1さんの現存コード、現在のKF標準設定、KF Imp5設定。既存Helix解析との比較は既存の検証メモを参照。
- 状態: コード・資料の調査結果。今回、新しいカットによる実イベント再解析やS/B測定は行っていない。
- 取扱い: 内部analysis noteを参照するローカル調査メモ。今回の作業ではコード・設定の変更、commit、pushは行わない。

## 1. 結論

qhu1さんのΛ選別は、現在のKF設定と同じではない。特に、daughterのprimary vertexからの統計的な隔たりを表すχ²、PV–崩壊点間距離、その有意度に、より厳しい条件がある。さらに、Treeを読み直す段階にもdaughterのχ²、pT、TOFの追加選別がある。**KFParticleを導入するだけで、Fig.9と同じ純度が得られるという比較にはなっていない。**

重要な相違は次のとおり。

1. qhu1側はoriginal daughterのχ²primary ≥ 10に加え、Λに関連付けて再フィットしたdaughterについても、下流でχ²primary ≥ 10を要求する。現在のKF設定はoriginal daughterに > 3で、後者の追加条件はない。
2. qhu1側のFinderはPV–DV距離 `l > 5 cm`、`l/σl > 5`。現在はそれぞれ `> 1 cm`、`> 3`。
3. qhu1側は下流でΛのpTを概ね0.5–5 GeV/cに限定し、TOF情報があるdaughterには固定m²窓も課す。現在のImp5はTOFを使用しない。
4. すべてがqhu1側で厳しいわけではない。parentのPV制約付きχ²/NDFはqhu1側が `< 5`、現在が `< 3`。Imp5にはさらにcos pointing ≥ 0.998、距離カット等がある。
5. Fig.9は**collider AuAu √sNN = 7.7 GeV、20–50%**。現在の入力は**Run20の13.5 GeV beam-energy FXT、√sNN ≈ 5.2 GeV**で、これまでの比較は0–80%。同じエネルギー・アクセプタンス・centralityの比較ではない。

以上から、daughter χ²と飛行距離・有意度は優先的に検証する価値がある。ただし、背景抑制の大きさと信号損失はまだ測定しておらず、「この変更でBG freeになる」とは結論できない。

## 2. 確認できたことと、確定していないこと

| 区分 | 内容 |
|---|---|
| PDFで確認済み | Fig.9のエネルギー、centrality、表示bin幅、本文のdaughter条件、sidebandによる背景評価の記述 |
| 現存コードで確認済み | Interface → Finder → Topo → Maker → Tree読出し → mass描画の各段階で適用される条件 |
| 未確定 | 現存コードのどの版・実行設定がFig.9の元ROOTファイルを作ったか |
| 未検証 | その条件を現在のFXTデータに適用した際のS/B、純度、信号収量、効率 |

Fig.9生成時のコード版・実行設定についてユーザーに確認したが、心当たりはないとの回答だった。現存の`analysis_set1.C`は9.2 GeV、`invariant_mass_cent2050.C`は19.6 GeV設定であり、これらをそのまま「Fig.9の7.7 GeV production設定」と扱わない。

また、[README_Alma9.txt](/star/u/qhu1/local_polarization/README_Alma9.txt:10)には、Tree作成後の解析をローカルUbuntu / ROOT 6.24/06で行い、マクロをRCFへコピーしてパスを変更した旨の記述がある。現存コピーと図作成時の実行環境が完全に一致する保証はない。

## 3. 調査対象と処理のつながり

### 3.1 資料・現在の設定

- [指定analysis note](../../../ref_analysisnote/Analysnote_local_lambda_hyperons_polarization_BESII_v3_20260601.pdf): PDFページ16、§2.3 / Fig.9。systematic variationsはページ33。データセット表はページ4。
- [現在のKF標準設定](../../../config/cuts/kf/kf_auau13p5_anaLambda_KFParticle.yaml)
- [現在のKF Imp5設定](../../../config/cuts/kf/kf_auau13p5_anaLambda_KFParticle_Imp5.yaml)
- [現在の13p5解析情報](../../../config/analysis/analysis_info_auau13p5_anaLambda_KFParticle.yaml): `production_13p5GeV_fixedTarget_2020`、P24iy、`mode: fxtmult`。
- [標準KFのevent設定](../../../config/cuts/event/event_auau13p5_anaLambda_KFParticle.yaml) / [Imp5のevent設定](../../../config/cuts/event/event_auau13p5_anaLambda_KFParticle_Imp5.yaml)

### 3.2 qhu1側の処理経路

```text
analysis_set1.C
  → StKFParticleInterface: track品質、TPC/TOF PID、primary/secondary分類
  → KFParticleFinder: daughter組合せ、崩壊点・飛行距離等の選別
  → KFParticleTopoReconstructor: parentのPV topology選別
  → StKFParticleAnalysisMaker: daughter再フィット、Tree保存
  → readTree_merge.C: daughter χ²・pT・TOF等の追加選別、mass histogram
  → invariant_mass_cent2050.C: 描画、sideband fit、背景評価
```

主要な参照箇所:

| ID | ファイル・確認箇所 |
|---|---|
| Q1 | [analysis_set1.C](/star/u/qhu1/local_polarization/lam_KF_BESII/analysis_set1.C:95): low-PV-track event除外、soft PID、`SetChiPrimaryCut(10)`、Λ/anti-Λ指定 |
| Q2 | [StKFParticleInterface.cxx](/star/u/qhu1/local_polarization/lam_KF_BESII/StRoot/StKFParticleAnalysisMaker/StKFParticleInterface.cxx:676): track受入れ。PIDはL313–412、primary分類はL491–495、再構成呼出しはL801 |
| Q3 | [KFParticleFinder.cxx](/star/u/qhu1/local_polarization/lam_KF_BESII/StRoot/KFParticle/KFParticleFinder.cxx:24): 初期値。実際のparent選別はL685–741、daughter距離はL1265 |
| Q4 | [KFParticleTopoReconstructor.cxx](/star/u/qhu1/local_polarization/lam_KF_BESII/StRoot/KFParticle/KFParticleTopoReconstructor.cxx:626): PV制約付きχ²/NDF `< 5`。Λ対象指定はL595、呼出しはL970 |
| Q5 | [StKFParticleAnalysisMaker.cxx](/star/u/qhu1/local_polarization/lam_KF_BESII/StRoot/StKFParticleAnalysisMaker/StKFParticleAnalysisMaker.cxx:709): parent mass、daughter再フィット、χ²計算、Tree保存（L789–809） |
| Q6 | [readTree_merge.C](/star/u/qhu1/local_polarization/histogram_macros/readTree_merge.C:470): daughter再構成・追加選別。default cut配列はL440–443、20–50% mass fillはL965–973 |
| Q7 | [invariant_mass_cent2050.C](/star/u/qhu1/local_polarization/analysis_macros/invariant_mass_cent2050.C:238): raw spectrum描画、rebin、背景fit、積分 |

現在側では、[StPicoKFParticleInterface.cxx](../../../StMaker/kfparticle/StPicoKFParticleInterface.cxx)のL103–231（track/PID）、L381–415（χ²分類・Topo実行）、[Finder](../../../StRoot/KFParticle/KFParticleFinder.cxx)のL825–869 / L1412–1416、[TopoReconstructor](../../../StRoot/KFParticle/KFParticleTopoReconstructor.cxx)のL719–730、[Maker](../../../StMaker/StLambdaKFParticleMaker/StLambdaKFParticleMaker.cxx)のL173–190を照合した。

## 4. 有効なカット条件の比較

以下のqhu1欄は、PDFの条件一覧だけではなく、**今回確認した現存コードのdefault経路**の条件である。Fig.9生成時と完全に同一と確定したものではない。

単位は、距離cm、運動量GeV/c、mass GeV/c²。χ²およびχ²/NDFは無次元。`l`の定義は§5参照。

### 4.1 KFの再構成・topology条件

| 条件 | qhu1現存コード | 現在KF標準 | 現在KF Imp5 |
|---|---|---|---|
| original daughterのχ²primary | p、πとも ≥ 10（Q1/Q2） | > 3 | > 3 |
| parentに関連付けたdaughterのχ²primary | 下流でp、πとも再度 ≥ 10（Q5/Q6） | この追加条件なし | この追加条件なし |
| FinderのPV–DV距離 `l` | 5 < l < 200（Q3） | 1 < l < 200 | 1 < l < 200 |
| Finderの `l/σl` | > 5（Q3） | > 3 | > 3 |
| Finderのdaughter間距離 | < 1（Q3） | < 1.5 | < 1.5。さらにMakerで ≤ 0.5 |
| parent geometric χ²/NDF | < 10（Q3） | < 10 | < 10 |
| parent PV-constrained topology χ²/NDF | **< 5**（Q4） | **< 3** | **< 3** |
| original daughterのgDCA下限 | このcm単位の条件は該当経路に見つからない | なし | p ≥ 0.7、π ≥ 1.0 |
| parentのPVまでの距離 | 上記topology条件。Imp5相当の固定上限は見つからない | 追加上限なし | KF状態で ≤ 0.5 |
| cos pointing | Imp5相当の ≥ 0.998 は見つからない | 追加下限なし | ≥ 0.998 |
| original Helixのpath長 | Imp5相当条件は見つからない | 追加条件なし | 両娘の絶対値 ≤ 100 |

補足:

- 現在のInterfaceでは `< 3` をprimaryと分類するが、その後のFinderにも `χ²primary > 3` の条件があるため、表の現在側の実効境界は `> 3`。
- qhu1のFinderには`fCuts2D[0]`を使うdaughter χ²条件があるが、L1241–1246ではコメントアウトされている。初期値3や`SetPrimaryProbCut(0.0001)`由来の約18.42を、そのまま有効なΛ daughter cutと読んではいけない。実効条件はInterface側の10と、下流の再フィット後の10。
- qhu1側もPV方向との整合性を要求する。Finderの`isParticleFromVertex`等を「pointing条件が一切ない」と表現してはいけない。ただし、固定の`cos ≥ 0.998`と同一条件ではない。
- parent topology条件はFinder後のTopoでも課される。Finderの候補保存だけで追跡を止めると、qhu1側の `< 5` を見落とす。

### 4.2 Track acceptance・PID・Λ acceptance

| 条件 | qhu1現存コード | 現在KF標準 | 現在KF Imp5 |
|---|---|---|---|
| daughter nHitsFit | ≥ 15（Q2/Q6） | ≥ 15 | ≥ 15 |
| nHitsFit/nHitsMax | Λ daughter経路では条件を確認できない | ≥ 0.52、nHitsMax > 0 | nHitsMax > 0の場合に ≥ 0.52 |
| nHitsDedx | 独立した下限は該当経路に見つからない | ≥ 5 | 追加条件なし |
| dEdxError | 0.04–0.12（Q2） | 0.04–0.12 | この品質窓は使用しない |
| daughter pT | 下流で0.15–5（Q6） | 0.15–10 | 追加acceptance cutなし |
| daughter η | PDFでは7.7 GeVで abs(η) < 1.5。ただし現存Q6に括弧の問題あり（§7） | lab η = −2.4–0 | 追加acceptance cutなし |
| proton / pionのTPC PID | `abs(dEdxPull(mode=1)) < 3`（Q2） | `dedx_pull`で3σ | PicoDst保存値の `abs(nSigma) ≤ 3` |
| TOF | 必須ではない。有効な場合は上流3σ判定に加え、下流m²窓（Q2/Q6） | 有効な場合にsoft 3σ判定。下流固定m²窓なし | 使用しない |
| 下流の固定TOF m²窓 | p: 0.5–1.5、π: −0.06–0.1、TOFがある場合（Q6） | なし | なし |
| parent pT | 保存parentで > 0.5、daughter四元運動量和で0.5–5（Q6） | 追加条件なし | 追加条件なし |
| anti-Λ | 再構成する | 再構成する | 再構成しない |

注意点:

- PDF本文の`hits > 15`、`pT > 0.15`という表記と、現存コードの`≥ 15`、`≥ 0.15`は境界が微妙に異なる。再現時はどちらに合わせるか明示する。
- qhu1のTPC PIDは保存済み`nSigmaProton/Pion`ではなく`dEdxPull`。両者に同じ「3」を設定しても、同一track集合になるとは限らない。
- qhu1の`mNSigmaDaughters = 3.0`はQ6で宣言されているだけで、その下流マクロの有効なTPC cutではない。実際のTPC判定はQ2まで追跡した。
- qhu1 Makerの`AcceptTrack()`にあるη、hits ratio等を、このΛ daughter経路で使われる条件と混同しない。
- 現在KF標準のTOF多項式係数はxwu2側の参照値であり、**Run20 AuAu13p5で校正・検証済みの係数ではない**。同じ3σという閾値だけでPIDの同等性を保証できない。
- Imp5の「追加acceptance cutなし」は、あらゆるtrackを使用するという意味ではない。PID、nHitsFit、gDCA、非零電荷、有限な運動量・共分散等の受入れ条件は残る。

### 4.3 Event・centrality条件

| 項目 | qhu1 / 指定Fig.9 | 現在KF標準 | 現在KF Imp5 |
|---|---|---|---|
| データセット | Fig.9はcollider 7.7 GeV、Run21 P22ib | Run20 13.5 GeV beam-energy FXT、P24iy | 同左 |
| mass比較のcentrality | Fig.9およびQ6の対象histは20–50% | これまでの比較は0–80% | これまでの比較は0–80% |
| vertex | 現存Makerは7.7 GeV相当分岐で abs(Vz) ≤ 145、原点からの半径 ≤ 2。Q6にもVz条件 | FXT Vz = 198–202、XY中心(−0.4, −2.0)、半径2 | 有限PVの安全検査。追加Vz/Vr cutなし。mode別中心はQA用のみ |
| その他event処理 | trigger、bad run、pileup、centrality、low-primary-track event除外等あり | 設定されたRefMult/VPD条件とcentrality等 | 既存Imp5のevent方針に合わせたmaxNTrとcentrality等 |

qhu1側のlow-primary-track除外は、使用trackに対するprimary track比が小さいeventを落とす処理で、現在のadapterには同じ処理を実装していない。event selection全体の同等性も、今回の調査では確立していない。

**13p5の意味を取り違えないこと。** 現在のファイル名の13p5はcolliderの√sNN = 13.5 GeVを意味しない。STARのbeam-use資料Table 8では、Run2020のFXT beam energy 13.5 GeVは√sNN = 5.2 GeVに対応する。[STAR Beam Use Request Runs 22–25](https://drupal.star.bnl.gov/STAR/files/STAR_Beam_Use_Request_Runs22_25.pdf)

colliderの対称なlab η条件や原点中心のvertex cutをFXTへそのままコピーすることは避ける。既存の`analysis.mode`に対応したvertexの扱いを維持し、acceptanceはFXTの運動学に合わせて検討する。

## 5. 同じ名前でも区別すべき量

### 5.1 original daughterと再フィット後daughterのχ²primary

上流のχ²primaryは、元trackから作成したKFParticleとPVとの整合性である。qhu1 Makerでは、その後daughterのコピーに対して:

```cpp
daughter.SetProductionVertex(parent);
daughter.TransportToProductionVertex();
```

を実行し、その状態から`GetDeviationFromVertex(PV)`を計算してTreeに保存する（Q5のL717–744）。Q6はその保存値に再度 ≥ 10を課している。

この2回のχ²条件は単なる重複ではない。現在の`interfaceChiPrimaryCut`を10へ変更するだけで、後段まで同等になるわけではない。再フィット後の値は現在の候補診断項目にも同じ形では保存されていないため、検証には明示的な計算・記録が必要になる。

また、ここで使うqhu1の`GetDeviationFromVertex`は共分散で規格化した二次形式 `dᵀ C⁻¹ d` に相当するχ²であり、平方根を取った「何σ」、TPC PIDのnSigma、cm単位のDCAとは異なる。参照: [KFParticleBase.cxx](/star/u/qhu1/local_polarization/lam_KF_BESII/StRoot/KFParticle/KFParticleBase.cxx:2873)。

### 5.2 Finderのl/σlと、PV制約後のdecay-length significance

今回の一様磁場経路でFinderが使用するのは、`GetDistanceToVertexLine`で求めたPV–DVの3次元距離とその誤差に基づく有意度である。参照: [KFParticleBaseSIMD.cxx](/star/u/qhu1/local_polarization/lam_KF_BESII/StRoot/KFParticle/KFParticleBaseSIMD.cxx:1406)、Q3のL701–741。

現在の候補には、このvertex-line系の量に加え、parentのコピーにPV制約を課した後の`decayLengthSignificance`も存在する（現在InterfaceのL287–300）。両者を同じ量として扱わない。

qhu1 Finderの `l/σl > 5` に合わせる第一候補は`finderLdL2D`であり、単に`minDecayLengthSignificance: 5`を追加することとは同一でない。QA・cut-flowにはどちらの定義かを残す。

### 5.3 parentのtopology条件とmassの定義

qhu1 Topoはparentの一時コピーに`SetProductionVertex(PV)`を課し、χ²/NDF `< 5`を判定する。これは元parentのmassをPDG Λ massへ固定する操作ではなく、この判定によって保存元parentのmassがPV制約付きコピーの値へ置き換わるわけでもない。

Finderには3σ mass選別を使う別の制約付き候補collectionがあるが、Makerへ返る元のparent候補一覧と区別する必要がある。「Fig.9がきれいなのはΛ massを固定しているから」とは、この経路からは言えない。

一方で、比較するmass observableには違いがある。

| 出力 | massの定義 |
|---|---|
| qhu1 Makerのraw invariant-mass hist | 元parentの`GetMass()`（Q5のL709 / L771） |
| qhu1 Q6の下流mass hist | Treeに保存した再フィット後daughterの運動量とp/π質量から四元運動量を作り、その和のmass（L470–499） |
| 現在のKF candidate mass | 元parentの`GetMass()`。parentのPDG mass constraintなし |

したがって、mass幅・収量の比較では、この定義の差も分離して確認する。どちらがどの程度背景を減らすかは未測定である。

## 6. Fig.9の「ほぼBG free」をどう読むか

- Fig.9上段のKFPは0.5 MeV/bin、下段のSDU “Topo”は0.1 MeV/bin。binごとの縦軸値をそのまま比較しない。また、この“Topo”という手法ラベルを、現在使用するクラス`KFParticleTopoReconstructor`との二者択一と解釈しない。
- Fig.9は背景が視覚的に小さいが、本文ではsidebandを二次多項式で評価し、背景を差し引いている。背景が数学的にゼロと示した図ではない。
- Q6では20–50%のmass histを、後続の狭いsignal window条件より先にfillする。したがって、このhistが狭いmass窓で切り取られたためにきれいに見えている、という処理ではない。
- 現存Q7もraw spectrumを描いてから、別コピーで背景を差し引く。積分範囲は1.1105–1.1205 GeV/c²、sideband fitは1.09–1.142から1.100–1.130を除外してpol2を使用する。ただし、これは現存マクロの設定でありFig.9当時との同一性は未確定。
- Q7の`signal to noise`という表示は計算上 `Nsignal/Ntotal = S/(S+B)`、すなわち**純度**。`S/B`と区別する。今後の比較では `S/B`、`S/(S+B)`、必要なら `S/√(S+B)`を名前と式を対応させて記録する。
- 現存Q6は0.5 MeV/binのhistを作るが、Q7には`Rebin(2)`もある。これも、現存マクロがFig.9そのものの再現設定と確定できない理由の一つ。

## 7. 現存コードをコピーする前の注意点

Q6のL480–495には、`fabs(Eta() >= 1.5)`のように、比較結果の真偽値へ`fabs`を適用する形の式がある。同様の判定が重複しており、意図された`fabs(Eta()) >= 1.5`と等価ではない。負側のη制限が意図通りに働かない可能性がある。

従って、PDFの「|η| < 1.5」を現存コードで正しく実行していると断定しない。また、この式を新実装へそのまま移植しない。現存コピーの問題であり、Fig.9生成時にも同じ不具合があったと断定する根拠はない。

同様に、宣言のみのcut変数、コメントアウトされたFinder条件、Λ経路で呼ばれないtrack選別関数を、有効な条件一覧へ混ぜないことが重要である。

## 8. 次に行う検証の優先順位（未実施）

既存`StLambdaMaker`、標準KF、Imp5の設定と出力を維持し、承認後に**独立したqhu1参考profile**を追加する方針が安全である。qhu1本番再現と、FXTに適したcut最適化は別の目的として扱う。

### P0: 比較の基準をそろえる

1. 同一PicoDstリスト・同一event範囲から出発し、event cutを通過した集合の差も記録する。まず20–50%の共通centralityで比較する。
2. parent/daughterのpT、FXTに適したηまたはrapidityの範囲、mass定義、bin幅、表示範囲、積分窓、背景fit方法をそろえる。
3. Fig.9の本当の再現を主張するには、7.7 GeV productionのコード版、実行マクロ、元Tree、最終hist/描画設定まで対応付ける。それまでは「現存qhu1コード参考」と明記する。

### P1: 主要なKF条件を一つずつ追加・走査する

1. original daughter χ²primary: 3 → 10。Interface分類と有効なFinder条件の両方を考慮し、境界も記録する。
2. Finderの `l/σl`: 3 → 5。
3. Finderの `l`: 1 → 5 cm。
4. qhu1と同じ手順で再フィット後daughter χ²primaryを計算・QAへ記録し、≥ 10の追加効果を分けて調べる。
5. parent/daughter pT条件、TOFありtrackへの固定m²窓を別々に検証する。TOFのRun20校正確認を伴わない単純な移植は避ける。
6. parent topology `< 3`と`< 5`、Imp5のpointing/DCA条件も比較する。すでに厳しい条件を残したままでは「qhu1と同じカット」にはならない。

各段階で、event/track/candidate cut-flow、mass spectrum、sideband背景、窓内SとB、統計誤差、S/B・純度、相対信号残存率を保存する。candidate数の減少だけを背景抑制と呼ばず、同じwindowでfitした信号数も確認する。embedding/MC等がない段階では、相対信号残存率を絶対再構成効率と呼ばない。

複数cutを同時に変更した最終profileだけでなく、一条件ずつの差分を残すことで、どの条件が有効だったか判断できる。

### 未解決事項

| 項目 | 現在の状態 | 解決に必要なもの |
|---|---|---|
| Fig.9の厳密なproduction版 | 未特定。ユーザーにも心当たりなし | 作成者のproduction macro / commit / ROOT出力 / 描画版の対応 |
| 現存Q6のη式とproduction版の関係 | 現存コードの問題を確認したのみ | 実際に図を作ったローカル版との照合 |
| FXTでのTOF PID | 現在の標準係数は参考値 | 対象runのTPC/TOF QAと校正確認 |
| BG抑制と信号損失の量 | 未測定 | 同一入力・共通条件でのcut scanと背景fit |
| publication相当の効率・系統誤差 | 今回の範囲外 | embedding/MC、acceptance/efficiency補正、cut/背景model変動 |

## 9. 再現性・既存結果との関係

今回の比較時点のstar-analyzer HEADは:

```text
1c79dfc78cdb587c1ff5b32cfe1a6ebd8a921ed5
Add full KFParticle Lambda reconstruction and Imp5 validation
```

これは現在側の基準であり、qhu1側のproduction commitではない。既存worktreeの無関係な変更は調査・メモ作成で変更していない。

関連する既存記録:

- [標準KFとHelixの10,000 event比較](lambda_helix_kf_comparison_10000_20260907.md)
- [Imp5を合わせた10,000 event比較](lambda_helix_kf_imp5_comparison_10000_20260908.md)
- [未使用cut設定の整理記録](../../configuration/unused_lambda_cut_cleanup_20260908.md)
- [KF ROOT出力の配置・条件一覧](../../../rootfile/auau13p5_anaLambda_KFParticle/mdfiles/README.md)

今回参照したファイルのSHA-256を以下に残す。**これらは調査時点のコピーを識別するためのもので、Fig.9 productionとの同一性を証明するものではない。** Q1–Q7は§3のリンクに対応する。

| 対象 | SHA-256 |
|---|---|
| 指定PDF | `660af49abac9f50711bbbc12a6090302946b66b9836e19a89b669e4a65b4a9ee` |
| 現在KF標準YAML | `93bc2d06afb18bd4efaaa5914c307240622400611d6764ec5a2bb3b2ff1f80fb` |
| 現在KF Imp5 YAML | `1fb8545991825c874709b5be70022024cfcdc4cced99058a7974ae6e6bcb85fa` |
| Q1 analysis_set1.C | `bf307f5ae2a83c566977e6aff18147cf830e552485fa0400b57383b55885c0d4` |
| Q2 Interface | `47c1c61ed8cec6c8addeb405b3dfdde7492ddee2281bbe9b9f426cc4328f9ce7` |
| Q3 Finder | `836d3231df9475aee368ca4647ea3a38a37291973eaeb17a6965a0bd8bddbcff` |
| Q4 TopoReconstructor | `0dbfc30ff1f6319c8e6794b0b7efc87694e61bdb66cc2cefbb5096bab21687e0` |
| Q5 Maker | `bdd2706a23d65f2cb1597ed39acfdfd23e5f198ce3da139fc24853ed15986533` |
| Q6 readTree_merge.C | `0c3860fe930584095a64238ba4f120967c0c1d8737ce1b07e92ca6c00dfc60a1` |
| Q7 invariant_mass_cent2050.C | `cc35d2aa7e3290dc549d067af8788e085fa74e528c4032a0e4a901a6b7792a7e` |

## 改訂履歴

| 日付 | 内容 |
|---|---|
| 2026-09-08 | PDFとqhu1現存コードの調査を保存。Finder後のTopo選別まで追跡した値、再フィット前後のdaughter χ²、mass定義、FXT/collider差、未確定事項と検証方針を記録。 |
