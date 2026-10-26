# Implementation plan: Lambda–nucleus femtoscopy with StFemtoMaker and KFParticle

作成日: 2026-09-29。状態: **実装前の計画（実装後の状態は11節）**。調査基準: `main` の `a11f942949da36ac4d7c5130860414e9c7c71f12` と、現在保持されているローカル変更・未追跡ファイル。

この文書を作成した段階では、解析コード・設定・ライブラリ・ROOT出力は変更していない。実装、ビルド、実イベント検証、job投入は未実施。

## 1. Goal and scope

`anaLambdaNuclearId.C` で行っていた Λ–原子核Femto解析を、原子核別の入口を維持しつつ **`StFemtoMaker` を利用する構成**へ整理する。Λは、今回採用した高purity条件で、既存のfull KFParticle Finder / TopoReconstructorにより再構成する。

| Entry macro | Bachelor species | Base channel |
|---|---|---|
| `anaFemtoLambda_d.C` | `deuteron` | `lambda_deuteron` |
| `anaFemtoLambda_t.C` | `triton` | `lambda_triton` |
| `anaFemtoLambda_3He.C` | `he3` | `lambda_he3` |
| `anaFemtoLambda_4He.C` | `he4` | `lambda_he4` |

- 各マクロは `#include "StMaker/StFemtoMaker/StFemtoMaker.h"` を用い、実際に `StFemtoMaker` を生成する。名前だけの変更や専用Makerの別名化では対応しない。
- 物理処理・pair loop・ヒストグラム出力はcompiled Maker / helper内へ置く。マクロは設定確認、StChain構築、`Init/Make/Clear/Finish` の制御だけにする。
- `anaLambdaNuclearId.C` の **対象原子核に関する全ヒストグラムと共通QA** を移植対象とする。特に全質量範囲の `k* × invariant mass`、SE/ME、signal/sideband、centrality別出力を必須とする。
- 3p85 / 13p5それぞれに設定を用意し、同じ実装で処理する。13p5も今回の対象は **beam energy 13.5 GeVのFXT** であり、colliderのsqrt(sNN)と混同しない。
- `StLambdaMaker`、旧 `anaLambda.C` / `anaLambdaNuclearId.C`、旧専用Makerのファイルは削除・置換しない。既存Phi/Kaon Femtoの選別・出力も維持する。
- 今回は新しいcut最適化やROOT6移行ではない。保存済みカットの移植とFemto出力の成立・回帰検証が目的。

## 2. Current implementation and gaps

### 2.1 Entry points and configuration

現行 `anaFemtoLambda_{d,t,3He,4He}.C` は、`StMaker/StLambdaNuclearMaker/StLambdaNuclearMaker.h` にある **`StFemtoMakerLambdaNuclear`** を使用している。runnerも `libStLambdaNuclearMaker.so` をロードしている。参照例の `anaFemtoPhiDeuteron.C` は `StFemtoMaker` を使用しており、両者はまだ別経路。

4Heの現行mainconfは `lambda:` / `femtoLambdaNuclearHist:` を参照し、内容もHelix Λ用cutである。`StFemtoMaker` に必要なspecies/channel定義がないため、includeだけの変更では動作しない。また `analysis_info_auau3p85_anaFemtoLambda_4He.yaml` の `baseAnaMacro: anaFemtoLambdaNuclear_4He` は実ファイル名と一致しない。全4種類で同様の整合性を検査する。

`ConfigManager` の `maker:` は現状 `PhiCutConfig` / `FemtoConfig` 等へ読み込まれる。新経路では `maker:` にFemto設定、`kf:` にKF Λ設定を置き、旧 `lambda:` の数値をKFへ暗黙適用しない。

### 2.2 Physics and histogram gaps

元の `anaLambdaNuclearId.C` は `StLambdaMaker`、`StNuclearIdMaker`、`StLambdaNuclearMixMaker` を利用し、一部のSE pair処理がマクロ内に残っている。ROOTのトップ階層にΛ QA、`true/` にNuclearId/SE、`mix/` にMEを保存する。

移植途中の `StLambdaNuclearMaker.cxx` には以下の相違がある。新実装ではこれらを互換性として再現しない。旧コード自体の修正は今回の対象外とする。

| Finding | Evidence | New implementation requirement |
|---|---|---|
| SE full-mass histogram is filled twice | `FillSameEventPairs` and `FillFemtoPair` both fill `hKstarMass_CentBin*` | One fill per accepted pair |
| ME full-mass histogram has only the forward direction | Direct fill exists for current Lambda × pooled nucleus, not the reverse | Fill both enabled directions exactly once |
| ME sideband centrality output is incomplete | Split Maker does not cover the original per-cent sideband family | Preserve the original family |
| Nuclear rigidity cut is lost in the split path | Original `StNuclearIdMaker` applies `maxPOverQ`; split `IdentifyNuclei` does not | Restore the original candidate selection |
| Pre-PID and final-candidate QA populations differ | Original nuclear kinematic QA precedes best-species/TOF selection | Keep stage-specific histogram semantics |
| Original merging QA confuses IDs and indices | `anaLambdaNuclearId.C` calls `picoDst->track(...)` with stored track IDs | Use verified Pico array indices |

ヒストグラム名だけでなく、型、軸、fill段階、重み、directory、SE/MEの意味を契約として固定する。KFとHelixで異なる再構成量を「同じ数値になるはず」とは扱わない。

## 3. Architecture

```text
run_anaFemtoLambda_<species>.C
  load/link libraries -> ACLiC anaFemtoLambda_<species>.C
    StChain
      StPicoDstMaker  [Event, Track, TrackCovMatrix, required PID branches]
      StFemtoMaker
        event policy + CentralityHelper (one decision per event)
        optional Lambda provider
          StPicoKFParticleInterface -> Finder / Topo -> shared final selector
        nuclear candidate builder (selected bachelor + legacy QA)
        same-event pairs / mixed-event pool
        canonical Femto histograms + legacy-compatible output
```

### 3.1 Optional Lambda provider

`StFemtoMaker` に、小さな **任意注入の候補provider interface** を追加する。汎用plugin registryや新しいevent-loop frameworkは作らない。

- interfaceはplain C++ / ROOT基本型だけを使用し、KFParticle/SIMD型を共通Femto headerに出さない。
- 具体的実装は `StMaker/kfparticle/` に置き、既存 `libStKfParticleCommon.so` に組み込む。factoryをLambda用ACLiCマクロから呼び、`StFemtoMaker`へ渡す。
- provider実装から `StFemtoMaker` の非inline関数を呼ばず、helper → Makerの逆向きlink依存を作らない。診断値を返してMakerがhistogramへ記録する。
- providerの責務は設定検証、mode-aware event policyの提供、イベントごとのfull KF再構成、最終cut、候補/診断値/エラーの返却。histogram ownershipや独自event loopは持たない。
- event policyは `StFemtoMaker` の対応箇所で一度だけ適用し、bad-run / pileup / centralityも一度だけ処理する。provider側で二重にcentrality選別しない。
- `Init`、`Clear`、破棄順、ownershipを明確にする。providerをMakerが所有する方式を基本とし、バッファには候補の値コピーだけを格納する。
- `species_lambda_builderType: resonance`、`species_lambda_particleKey: lambda_kf` を新設する。要求されたのにproviderや `kf:` がない場合は初期化失敗。Helixやdefault KFへの自動fallbackは行わない。
- species/channel数、対応するbuilder key、channelの両species、明示mass windowを初期化時に検証する。現在のPhi/proton向けの欠落key defaultはLambda設定に使用しない。現行の `SpeciesDef.cutsRef` は未使用なので、これでKF YAMLが接続されたとは扱わない。
- Phi/Kaon解析ではproviderを生成せず、従来経路を使う。

この分離が必要な理由はビルド依存にある。現在 `StFemtoMaker` は `make core` 対象であり、KF特有のcompile/link設定は `*KFParticleMaker` に限定されている。PIMPLだけではELFの共有ライブラリ依存は消えないため、`StFemtoMaker` から具体的KFクラスを直接参照して、既存Phi解析までKFライブラリ必須にする構成は採らない。

### 3.2 Shared KF selection

既存 `StPicoKFParticleInterface::ProcessEvent()` をfull PicoDstに対して使用する。Femtoのprimary/bachelor track loopを通過したtrackだけをKFへ渡してはならない。Λ娘はdisplaced trackであり、その前処理により選別が変わるためである。

`StLambdaKFParticleMaker::PassCandidateCuts()` を共通selectorへ抽出し、単独Λ MakerとFemto providerの双方が同じ関数を呼ぶ。これはstandalone KF Makerに対する唯一の必須選別リファクタリングとし、cut値、境界の不等号、disabled値、Imp5 path guardを維持する。上流Finder/Topo cutと、この最終selectorの**両方**を適用する。

`KfLambdaCandidate` には元のtrack IDとPico array indexの両方がある。Femtoへの移送は `protonIndex` / `pionIndex` を用い、各娘を `(eventIndex, trackIndex)` で記録する。既存standalone MakerのID getterを、そのままPico indexとして流用しない。

### 3.3 Momentum and mass conventions

- **Mass axis:** mass-unconstrained KF candidateの再構成質量を保存する。Λの質量制約でピークを人工的に狭くしない。
- **Lambda momentum:** 同じmass-unconstrained candidateの3-momentumを使う。PV-constrained copyから得る距離/有意度と混同しない。
- **k\*:** 旧解析と同じくΛの固定PDG質量と核種別質量で四元運動量を構成し、pair rest frameでの片粒子の運動量を使う。`q_lab = |p_Lambda - p_nucleus|` も別に保存する。
- Λの `reso.invMass` は測定質量、`P4()` のenergyは明示的な `lambdaPairMassMode: fixed` / `lambdaPairMass` に基づく。数値はYAMLへ置き、既存Phiのmass conventionは変えない。
- ³He/⁴HeはTPC rigidityから真の運動量へのZ=2補正を**一度だけ**行う。既に補正した運動量を再度 `NuclearP4` に渡して二重補正しない。
- signal/SB判定、k*、mass-vs-k*はSEとMEで同じ定義を使う。

## 4. Configuration and preservation

### 4.1 Adopted Lambda cuts

新しいmainconfは、次の保存済みファイルを `kf:` で直接参照する。原子核ごとに38個のKF parameterを複製しない。今回これらの値は変更しない。

| Dataset | Preset |
|---|---|
| AuAu13p5 FXT | [kf_auau13p5_anaLambda_KFParticle_highpurity.yaml](../../../config/cuts/kf/kf_auau13p5_anaLambda_KFParticle_highpurity.yaml) |
| AuAu3p85 FXT | [kf_auau3p85_anaLambda_KFParticle_highpurity.yaml](../../../config/cuts/kf/kf_auau3p85_anaLambda_KFParticle_highpurity.yaml) |

共通の `finder_topo / lambda_imp5 / pico_nsigma`、TOFなし、anti-Lambdaなし、保存質量範囲1.05–1.25 GeV/c²を維持する。設定が不足・不正なら停止する。

以前の独立100kで得たpurityは13p5で98.67%、3p85で97.73%だったが、これはinclusive Λのsideband推定である。Λ–原子核pair、低k*、centrality別や別productionで同じpurityを保証しない。Femto QAで改めて確認する。

### 4.2 Event selection must be explicit

現行の旧Λ解析および移植途中の専用Makerは、明示的event cutとして主に `maxNTr` を使用している。`StFemtoMaker` の通常経路はVz/Vr/RefMult/VPDも選別するため、無条件に置き換えるとΛカット以外で収量が変わる。

初回移植では **`lambdaEventPolicy: lambda_imp5_compat`** を明記し、保存済みKF studyと旧経路に対応する finite-PV安全確認 + `maxNTr` + centrality選別を使用する。娘粒子の `selectionProfile` からevent policyを暗黙決定しない。

- `analysis.mode` の `fxtmult/refmult` と `centrality.mode` の一致を検証する。
- `qaVertexByMode` の中心はQA座標であり、半径カットではない。稼働していないVz/Vr等を「有効なcut」として新設定に残さない。
- 物理的なmode別vertex cutを使う場合は、明示的な別policy `mode_vertex` と `vertexByMode` を使用する。元のcutとして無断で追加しない。
- このLambda policyを使う間はgeneric Femtoの `PassEventCuts` / `PhiCutConfig.maxNTr` を重ねない。VR QAも選択されたpolicyの中心で計算する。
- mixingのVz区間はevent acceptanceとは別の役割で、`mixing:` に明示する。Imp5ではvertex acceptanceが無いので、範囲外の扱いも設定・counterに残す。初期互換設定は旧clamp挙動を明記し、黙ってイベントを落とさない。
- cent9の意味を変更しない。pull後の `StRefMultCorr` のrun coverage / parameter indexを記録し、過去studyのcentrality上の注意点を引き継いで監査する。較正変更は今回に混ぜない。

### 4.3 Nuclear selection

核種別entryには `nuclearSelectionProfile: legacy_nuclearid` を明記し、元の `StNuclearIdMaker` の実効選別を基準にする。generic Phi-bachelorの追加pT/eta/DCA/TOF条件を自動適用しない。

- 核種の判定は `StNuclearIdHelper` / 既存PID計算を再利用する。ただし旧predicateとの一致をテストしてから採用する。
- 旧実効条件の `nHitsDedx >= 15`、`gPt >= 0.1`、正のdE/dx、TPC nSigma、best-species、`maxPOverQ`、任意の `m2_selection` を、新経路のYAML上で一意に定義する。hardcoded条件は値を維持して設定化する。
- 原子核候補用loopは既存generic `PassTrackCuts` の手前で独立させ、Λ娘と原子核で異なる選別を維持する。
- inclusive/TPC-hypothesis QAは最終best-species/TOF cutの前、candidate QAは後に埋める。QAのための仮説評価と、pairに使う対象核種の選別を分ける。
- 符号、TOF未検出時の扱い、PID同点時の選択、rigidity境界、Heの補正は回帰テストで固定する。新しい核種purity最適化は行わない。
- helperには、元の最終NuclearId選別にはない `pMag >= 0.1` やTOF選別での `minP_M2cut` 条件がある。この違いを無視して直接置換しない。必要な互換predicateはLambda用profile内へ分離し、既存Phi用helperの挙動は変更しない。

### 4.4 Config layout and rollback

新規production候補は、例として `main_auau3p85_anaFemtoLambda_4He_KFParticle_highpurity.yaml` のように区別し、energy × nucleusの8セットを作成する。macro名は既存の4つを維持し、`anaName` / ROOT出力先にはKF高purity条件を識別できるsuffixを付ける。

| Mainconf key | Responsibility |
|---|---|
| `analysis` | Dataset, mode, macro names, STAR release, output identity |
| `event` | Explicit Lambda event policy inputs / mode profiles |
| `centrality` | Bad runs, pileup, calibration, accepted bins, weights |
| `kf` | Saved Lambda high-purity preset |
| `nuclearid` | Nuclear PID and legacy quality requirements |
| `maker` | Species, channels, mass windows, pair conventions, QA switches |
| `mixing` | Pool bins/range, depth, direction, sampling policy |
| `femtoHist` | Histograms for `StFemtoMaker("femto", ...)` |

新経路で使わない `lambda:` / `v0:` / 重複PID値は新mainconfへコピーしない。既存Phiが共有する設定loaderのdefaultを変更せず、新Lambda経路では必要keyと実際の消費先を検証する。保存済みKF YAML、旧解析mainconf、旧出力を上書きしない。

4つの既存FemtoLambdaマクロは今回の移植対象であるため、その動作は変更される。旧per-species設定は新しいFemto設定と互換ではない点を明記し、新マクロで読んだ場合は黙って別解析をせず説明付きで停止させる。実装前に**旧マクロ/runner/設定/未追跡専用Makerを含む完全なsnapshotとhash**を保存し、旧状態に戻す手順を記録する。`git HEAD`だけでは未追跡ファイルを復元できない。元の `anaLambdaNuclearId.C` 経路はそのまま残す。

## 5. Histogram compatibility contract

### 5.1 Legacy output

`S = d, t, 3He, 4He`、`i = 0..8`。各原子核別ファイルに対象 `S` のpair histogramと共通QAを出力し、4種類を合わせて旧解析の全核種をカバーする。以下は代表ではなく、**各familyの全suffixを検査対象とする**。実装開始時に旧YAMLと全Fill箇所から全keyの機械可読manifestを作る。

| Family | Directory / keys | Fill population |
|---|---|---|
| Lambda QA | `/hLambda_InvMass`, kinematics, topology, daughter PID, mass correlations, centrality QA | Corresponding reconstruction/selection stage |
| Nuclear dE/dx QA | `true/hDedxP`, `_cut`, `_e`, `_pi`, `_K`, `_p`, `_else`, `_d`, `_t`, `_3He`, `_4He`, existing cent suffixes | Inclusive / TPC-hypothesis / post-cut stages as in legacy |
| Nuclear kinematic QA | `true/hPvsEta_S`, `hPvsY_S`, `hPvsPt_S`, existing `_CentBin{i}` variants | TPC-selected hypothesis before final PID |
| Nuclear TOF QA | `true/hM2P`, `true/hDedxP_d_m2`, `true/hDedxP_t_m2`, `true/hDedxP_3He_m2` | Raw-rigidity-derived mass squared and legacy TOF QA selections |
| SE signal | `true/hKstar_S`, `true/hQlab_S`, `true/hKstar_S_CentBin{i}` | Signal mass window |
| SE sidebands | `true/hKstar_S_SBPos`, `_SBNeg`, corresponding `_CentBin{i}`; `true/hQlab_S_SBPos`, `_SBNeg` | Right / left sideband |
| ME signal | `mix/hKstar_Mixed_S`, `mix/hQlab_Mixed_S`, `mix/hKstar_Mixed_S_CentBin{i}` | Signal window, cross-event only |
| ME sidebands | `mix/hKstar_Mixed_S_SBPos`, `_SBNeg`, corresponding `_CentBin{i}`; `mix/hQlab_Mixed_S_SBPos`, `_SBNeg` | Right / left sideband, cross-event only |
| SE full mass | `true/hKstarMass_S_CentBin{i}` | All selected Lambda candidates before signal/SB classification |
| ME full mass | `mix/hKstarMass_Mixed_S_CentBin{i}` | Same full-mass population, both enabled ME directions |
| Track-merging QA | `true/hDphiDeta_proton_S` | Signal SE; Lambda proton daughter versus nucleus |

現行13p5/3p85のNuclearId hist YAMLには `hDedxP_4He_m2` はない。追加する場合は新QAと明記する。また `hDedxP[_S]_CentBin{i}`、`hPvs{Eta,Y,Pt}_S_CentBin{i}`、`hVz_S`、`hMult_S[_CentBin{i}]` は元Makerにoptional fillがあるが、現在の該当hist YAMLでは定義されていない。**実際に設定されている出力と、コード上だけのoptional fillをmanifestで区別する。** 対象外核種のpairを空histだけで「解析済み」に見せない。共通のPID QAは、元の複数核種仮説を保持する。

`hDedxP_cut` は単なるno-TOF sampleではない。未match/無効TOFに加え、広いm²・p条件を満たさないmatched trackも含むので、元のfill条件を保持する。

トップ階層のΛ QAは次の全keyを対象とする。列挙は新しいfit精度の保証ではなく、保存内容の契約である。

- `hLambda_InvMass`、`hLambda_InvMass_CentBin{i}`、`hLambda_Pt`、`hLambda_Eta`、`hLambda_Phi`。
- `hDCA12`、`hDCAV0`、`hCosPointing`、`hNSigmaProton`、`hNSigmaPion`。
- `hLambda_InvMass_vs_{Pt,DecayLength,Y,TransDecayLength}` は X = その物理量、Y = mass。`hDCAV0_vs_InvMass` / `hCosPointing_vs_InvMass` は逆に X = mass、Y = geometryであり、転置しない。
- `hLambda_InvMass_vs_Cent9` / `hLambda_InvMass_vs_RefMultCorr` は X = cent9/refMultCorr、Y = mass。
- `hVz`、`hRefMult`、`hCentrality`、`hCentralityRaw`、`hCentrality16`、`hRefMultCorr`、`hRefMultWeight`、`hRawMult`、`hRefMultVsNTOFMatch`、`hRefMultVsNTOFMatchAfter`、`hCentralityVsVz`、`hRawMult_vs_RefMultCorr`、`hN`。
- `h{RawMult,RefMultCorr,NTracks,TofMatchMult,NProtonCand,NPionCand,NLambdaPairs}_vs_Cent9`。

| Observable | Legacy axis contract |
|---|---|
| `hKstarMass*` | TH2; X = k*, 200 bins, 0–1 GeV/c; Y = M(p pi), 200 bins, 1.05–1.25 GeV/c2 |
| `hKstar*` | 200 bins, 0–1 GeV/c |
| `hQlab*` | 200 bins, 0–2 GeV/c |
| `hDphiDeta_proton_S` | X = DeltaPhiStar, Y = DeltaEta; 200 bins each, -0.25–0.25 |

**全質量のTH2は狭いsignal windowで絞る前に埋める。** sideband候補をprovider段階で削除しない。既存YAMLが未fillのkeyを持つ場合も一覧化し、意味のあるfillを実装するか、未使用と判定した理由を個別に記録する。単に空の同名histを作って完了とはしない。

Λ topology QAはKF値と旧Helix値の定義差を記録する。例えば decay vertex–PV の幾何学距離とPV-constrained decay lengthは区別する。旧QAと同じ定義で計算できる量は保持し、KF固有量は別名で追加する。旧Helix fitをもう一度実行して選別を混ぜない。

### 5.2 Canonical Femto output and single-fill ownership

新しい共通形式として `hKstarSE_<channel>` / `hKstarME_<channel>` / `hKstarSEVsCent_<channel>` / `hKstarMEVsCent_<channel>` も用意する。channelにはfull-mass用base、`_signal`、`_leftSB`、`_rightSB`を持たせる。

- pair計算・mass分類・cut判定は一箇所に集約する。
- legacy出力は同じaccepted-pair結果から作る。signal/SB各channelの処理ごとにfull-mass TH2を再fillしない。
- full-mass TH3を追加する場合でも、要求されたlegacy TH2を省略しない。Phiの `hPhiMKK_vs_Kstar*` へΛを詰めない。
- directory書込みは同名histのROOT cycleによる上書きを防ぎ、SEを空ME histで置き換えない。ROOT所有権もテストする。
- 旧 `anaLambdaNuclearId` 同様、`mix/` にはMEを出力する。移植途中MakerのようにMixed histを `true/` にも重複保存しない。核種名の残った誤ったtitle（4HeなのにLambda–d等）も新YAMLで修正する。
- 標準のhist class、bin edges、under/overflow、errors/Sumw2、weight policyをmanifestで比較する。

### 5.3 Mass windows and merging QA

初期signal/SBは元解析の `MeanLambda`、`sigmaLambda`、`nsigmaLambda` と外側係数4を一つのYAML定義から生成する。signal境界の包含とSB境界の除外も保持する。以前のpurity studyの固定窓 `[1.110,1.122)` とは異なるため、自動置換せず、QAに両定義を明示する。

track-merging QAは、旧 `picoDst->track(id)` の誤参照を修正したindex-basedな娘対核種で作る。旧固定磁場0.5 T、評価半径1.4 mはコードに埋め込まず、互換QA設定として記録する。実イベント磁場を使うcorrected QAは別名で併記し、kG/T換算とHeのrigidity/chargeを検証する。close-pair vetoを今回黙って追加しない。

## 6. Event mixing and CF

- `StFemtoMaker` のpoolを使用し、同じVz bin × cent9 × 設定されたEP binだけを混ぜる。混合後に現在イベントをpoolへ追加する。
- 初期Lambda profileでは `nCentralityBins: 9` を必須とし、cent9をclampして異なるcentrality binを混ぜない。
- 初期比較では `mixingMode: bufferAll` / `mixBothDirections: true` を明示し、上限samplingに由来する分布差を避ける。pool depthと旧 `maxMixEvents` の実効差も記録し、互換検証は同じpool条件で行う。
- 初期Lambda profileは `bufferAll` に限定する。将来 `randomSample` を対応させる場合はbase channel単位で一度sampleし、同じpairをfull/signal/SB/legacy出力へ振り分ける。mass channelごとの独立samplingで出力間の整合性を崩さない。
- current Λ × pooled nucleusとpooled Λ × current nucleusを独立に処理する。希少核種では、現在イベントに両speciesが揃っていなくても有効な片方向を捨てない。
- `(eventIndex, trackIndex)` が一致するtrack sharingはSEで除去する。別イベントの同じlocal indexを誤除去しない。
- mixing counterでforward/reverse、eligible、filled、rejectedを確認する。mass TH2とsignal/SB k*は同じME母集団から得る。
- MakerはSE/MEを保存し、CFはmerged histogramを読むQA macroで算出する。centrality結合はSE/MEを足してから正規化し、完成したCFを平均しない。
- 既存のΛ QA/CFツールが `true/` / `mix/` の対象histを読めることを確認する。最低限、核種別のmass、mass-vs-k*、SE/ME、raw CFを確認できるcheckHist入口を用意する。Phi専用のMKK fittingやbackground modelをΛへ無検証で転用しない。

## 7. Implementation sequence

| Step | Work | Acceptance gate |
|---|---|---|
| 0 | Snapshot dirty sources/configs and current output contract | Recoverable snapshot; exact histogram manifest; no unrelated changes |
| 1 | Extract shared final KF selector; introduce neutral provider API and concrete KF provider | Standalone candidate/output closure; Phi core has no KF dependency |
| 2 | Add Lambda dispatch and explicit event/nuclear policies to StFemtoMaker | Synthetic candidate/event tests; no implicit Phi cuts |
| 3 | Add legacy-compatible histograms and single-fill pair output | Key/axis/stage coverage; SE/ME counting closure |
| 4 | Migrate four macro/runner pairs; add energy/species mainconf sets and wrappers | Correct macro names; explicit config; TrackCovMatrix checks |
| 5 | Build all affected libraries and run tests in STAR SL24y / ROOT5 container | `make core` and `make all`; runner loading and ACLiC successful |
| 6 | Real-input validation for both energies and four nuclei; QA | Event/candidate ledgers, histogram/kinematic checks, no missing input |
| 7 | Document config usage, differences, QA and rollback | Reviewed result report; no automatic farm submission |

主な変更予定は `StMaker/StFemtoMaker/`、新provider interface、`StMaker/kfparticle/`、`StLambdaKFParticleMaker` のselector委譲箇所、`FemtoConfig`等の必要な設定型、4組のanalysis/runner/script、新規mainconf/cut/hist/QA/testである。`StLambdaMaker`、旧NuclearId/Mix Makerや既存Phiのproduction YAMLへ無関係な変更を加えない。

新KF helperは既存の自動検出を使う。新しいMakerライブラリを核種ごとに作らない。Makefile調整が必要ならheader/test依存等に限定し、`-msse4.1`、`HomogeneousField`、KF ABI flagsをcoreへ広げない。

runnerは少なくともconfig → RefMultCorr → KF core → KF helper → common → StFemtoMakerの依存を満たしてロード/linkし、Lambda用compiled macroをACLiCで実行する。Phi runnerはKF依存なしで動くことを別途検証する。

mainconf・jobid・入力・出力は引数で渡す。`.current_mainconf`やhardcoded mainconfへのfallbackを新経路に持ち込まない。解決後pathと出所を必ずログに出す。必要設定/hist/providerの欠落、初期化エラー、入力枝不足は非zero終了にし、部分出力を成功扱いしない。

## 8. Validation details

### 8.1 Build and regression

`./script/singularity_make.sh <new-mainconf>` でbatch-matched buildを行う。pull後のsourceと既存`.so`が一致するとは仮定しない。config structを変更した場合は `libStarAnaConfig.so` を含む全依存先とACLiCキャッシュを再構築する。

- `make core` / `make all` と既存KF full-chain / Pico adapter試験。
- providerなしの既存Phi/Kaon runnerのロード試験、短い固定入力でのhistogram nonregression。
- standalone `anaLambda_KFParticle` の同じ入力・同じ設定でのrefactor前後比較。
- 既存 `anaLambda` / `anaLambdaNuclearId` が削除・置換されておらず、従来構成でロード/実行できること。

### 8.2 Deterministic tests

- KF全最終cutの境界、disabled設定、path guard、invalid covariance、ID/index mapping。
- 同じrun/eventとKF設定で、standaloneとFemto providerのselected候補多重集合、娘index、mass/momentum、診断値が一致すること。
- KF娘がgeneric primary/bachelor cutに巻き込まれないこと。event policyは同じ条件で照合し、異なるpolicyでの差を候補再構成差と混同しない。
- 核種PID、rigidity上限、Z=2補正、TOF有無、QA前後段階の旧コードとの一致。
- raw massを変えても固定PDG質量のk*が変わらないこと、PRF計算が独立式と一致すること。
- SEのshared-daughter除去、MEの異event同index許可、片speciesだけのイベント、両ME方向、pool境界。
- pair全質量histが一回だけ埋まること。signal/SB/full-rangeの整合性は境界・flowを考慮して照合する。
- 2ファイル目だけTrackCovMatrix欠落、Track/Cov数不一致、不正設定、既存出力path、途中I/Oエラーに対する失敗試験。

### 8.3 Real PicoDst validation

最初は各energyで同じ先頭10,000入力イベントを4核種に使用する。入力リスト順、実際に読めたfile/entry、run/event、accepted/reconstructedイベント、候補数を記録する。入力不足を10,000イベント成功として扱わない。

直近studyの3p85は2019 production、現行per-species analysis_infoは2021 productionである。この差を隠さず、まず検証済み入力で移植closureを行い、2021へ適用する場合はdataset/run情報を正しくそろえた別検証とする。analysis_infoの古いrunRangeやmacro名をそのまま新設定へコピーしない。

⁴He/³Heでは10kでpairがほとんど出ない可能性がある。0-entryを不具合とも成功証明とも断定しない。合成候補で必須fill経路を検証し、実データが不足する場合は明記する。必要な追加100k局所検証は実装時に所要時間を報告して進め、farm jobは別途確認なしに投げない。

QAはinclusive Λ massに加え、**pairに入ったΛのmass / k*別mass / 核種別・cent別のpurityと候補数**を確認する。高purity cutによるpair収量・受容の変化はinclusive signal retentionや真の検出効率と区別する。

## 9. Completion criteria and non-goals

完了条件は、4種類の入口が実際に `StFemtoMaker` を使い、保存済み高purity KF設定を読み、元解析の対象histogram一式と全質量TH2をSE/MEともに正しく出し、既存解析の回帰試験を通ること。

ROOT6導入、KFParticle upstreamの再移植、cut再最適化、centrality較正変更、close-pair vetoの新規有効化、旧Makerの削除、既存ROOTの上書き、farm大量投入、commit/pushは本計画の自動実行対象ではない。

計画・実装ではリポジトリの `femto-species-naming`、`add-new-analysis`、`stmaker-add-histograms` 手順を適用し、channel命名、histogram lifecycle、StChain/ACLiC構成、既存解析非退行を優先する。centralityとビルドの確認には対応する手順を使用する。

## 10. Sources inspected

- [Original analysis](../../../analysis/anaLambdaNuclearId.C), [current 4He entry](../../../analysis/anaFemtoLambda_4He.C), [Phi–deuteron reference](../../../analysis/anaFemtoPhiDeuteron.C)
- [StFemtoMaker implementation](../../../StMaker/StFemtoMaker/StFemtoMaker.cxx), [naming and mixing rules](../../../StMaker/StFemtoMaker/README.md), [candidate provenance](../../../include/FemtoCandidate.h)
- [Intermediate split Maker](../../../StMaker/StLambdaNuclearMaker/StLambdaNuclearMaker.cxx), [original NuclearId Maker](../../../StMaker/StNuclearIdMaker/StNuclearIdMaker.cxx), [original ME Maker](../../../StMaker/StLambdaNuclearMixMaker/StLambdaNuclearMixMaker.cxx)
- [Legacy nuclear histogram definitions](../../../config/hist/hist_auau3p85_anaLambdaNuclearId.yaml), [legacy Lambda histogram definitions](../../../config/hist/hist_auau3p85_anaLambda.yaml)
- [KF Pico interface](../../../StMaker/kfparticle/StPicoKFParticleInterface.h), [KF candidate schema](../../../StMaker/kfparticle/KfParticleHelper.h), [KF event policy](../../../StMaker/kfparticle/KfEventSelection.cxx), [standalone final selector](../../../StMaker/StLambdaKFParticleMaker/StLambdaKFParticleMaker.cxx)
- [Validated purity study](../../kfparticle/comparisons/lambda_kf_upstream_purity_20260928.md), [centrality conventions](../../../StRoot/StRefMultCorr/README.md), [build dependencies](../../../Makefile)

## 11. Implementation status (2026-09-30)

本計画に沿って4核種の入口を StFemtoMaker + full KF provider へ移植した。**比較対象は元の `analysis/anaLambdaNuclearId.C` とし、未デバッグだった旧 `anaFemtoLambda_4He.C` は収量基準に使用していない。**

実装内容、保存済みカットとの接続、ビルド・試験結果、10,000イベントの比較、QA、復元方法は [実装・検証記録](../implementation/femto_lambda_nuclei_kfparticle_20260930.md) に集約する。冒頭〜10節は実装前の設計履歴として残す。

新4核種×両エネルギーの10k実行、単独KF候補照合、元の anaLambdaNuclearId との原子核QA・SE/ME全出力監査がすべて成功した。13p5の旧Helix参照は最初の実行がsignal終了したため、元ログを残して単独再実行し、10,000イベント・exit 0で完了した。

ROOT5/SL24yの全ビルド、KF既存回帰試験、合成試験、入力・不正設定試験も成功。旧設定・保護対象289ファイルは実装前と同一。既存Phi/Kaonはcompile/loadとKF非依存を確認したが、該当mainconfがないため実イベント収量の非退行は未検証。核種・k*・cent別purityの数値fitは今回行っておらず、低統計CFと合わせ物理検証の残項目として報告書に明記した。farm投入・commit/pushは行っていない。

