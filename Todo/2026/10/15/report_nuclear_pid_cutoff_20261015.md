# 原子核PID cutoff研究：実施状況（2026-10-15）

## 現時点の状態

**2026-10-17（日本時間）の再開分は完了しました。** 前回の最終80条件集計を使用する[独立mainconf](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_purity50_20261015.yaml)と[核PID YAML](../../../../config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_20261015.yaml)を作成し、同条件のbatch一致build、20イベントの接続試験、閉ROOT検査、回帰試験、最終PDFを確認しました。**新条件は明示的に選べる独立プロファイルとして提供し、既存本番mainconfは切り替えていません。**

以下の段階別記録は途中経過を含みます。最終数値は「2026-10-15終了時点の集計・再開時の引継ぎ」、今回の検証と成果物は末尾の「2026-10-17（日本時間）の再開・完了結果」、過去の日次ログは[2026-10-15](../../../../analysisnote/20261015/summary20261015.md)を参照してください。

## 確認・実施した内容

- [依頼と計画](plan_nuclear_pid_cutoff_20261015.md)を整理しました。
- 前日のALL purity研究は準備ファイルまでで、最終計算は未実行でした。旧⁴Heの参考purityは、背景固定・核種帰属・誤差評価の制限があるため、今回の下限値として流用しません。
- d、t、³HeのTPC-only ALL purityは新たに評価する必要があります。
- 全4核種の高R TOFと核種別下限を扱う後方互換のコード、および専用テストを追加し、batch一致環境でのbuild・検証を進めています。旧設定ファイルの内容は保持しています。

## 入力と限界

- 対象ROOTは `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root`。既存の26956/26994 subjobの部分統合であり、38件を含む完全統計ではありません。
- 解析用cacheは `pdf/nuclear_pid_tof_20261014/inputs_001/cache.root`。ALLとTOF-tagのnative二次元分布を使用します。保存済みhistogramからの研究であり、今回新たにPicoDst全量を読み直すものではありません。
- 高Rの新TOF条件を適用した候補のpurityを、TPC-only ALLの二次元分布だけから直接決定することはできません。低R cutoff評価と高R TOF実装は区別して報告します。
- TOF-tagからALLへの応答転用には仮定があるため、モデル差と選択窓内残差を必ず記録します。

## 残作業

前回残した独立YAML・mainconf、同条件のbuild・先頭20イベント、最終PDFと索引の更新は完了しました。本番再投入、centrality／CFペア条件付き核purity、高R TOF選択後のpurity・収量の検証は今回の接続試験の範囲外であり、実施していません。

## ALL fitの第1段階：広い共通窓の診断

Gaussian／log-normal × positive linear／exponential背景について、R=0.4〜2.0、幅0.01の計640スペクトルを実行しました。TOF-tagで求めた同じnative R binの応答をALLへ固定転用し、ALLの成分収量・背景をfitする条件付き計算です。PID窓は凍結した現行設定のままです。

| ALL応答 | 背景 | 実行区間 | 数値条件を通過 |
|---|---|---:|---:|
| Gaussian | linear | 160 | 160 |
| Gaussian | exponential | 160 | 83 |
| log-normal | linear | 160 | 160 |
| log-normal | exponential | 160 | 81 |

[初回比較表](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/summary_nominal_001/summary.tsv)と[判定・詳細JSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/summary_nominal_001/summary.json)を保存しました。0.05幅でのpurityはnativeごとの信号数・総数を合計してから計算しています。

**この段階では、全核種について単一の推奨下限は確定できません。** ³Heには0.95〜1.30の支持区間がありましたが、2.0まで連続して支持されているわけではありません。d・t・⁴Heにもpurity推定値が高い区間はありますが、全比較モデルでの支持を確認できず「未評価」としています。

診断で次を区別しました。

- 高統計のχ²だけを理由に除外した結果ではありません。
- exponential背景の数値上限は、元のTAG局所窓用の傾き上限8を広いALL共通窓へそのまま適用しており、境界に達する例があります。数値guardの不足と背景関数自体の不適合を切り分ける必要があります。
- d〜⁴Heを跨ぐ広い共通窓で一つの背景形を当てると、稀な⁴Heがbackgroundへ吸収される例があります。例えばR=1.505で名目⁴He収量が0になるモデルがあり、数値収束だけでは適切なモデルとはいえません。
- 過去TAG exponential fitにも全体境界の例があり、中心・幅の支持判定への影響を調べています。

採用条件を緩めて未評価点を通すのではなく、核種ごとの局所fit窓による背景モデルの検証を次に行います。旧広窓の結果・失敗点は診断結果として保持します。まだ中心・幅の係数や最終cutoffを置換しません。

## ALL fitの第2段階：核種ごとの局所窓

全4核種それぞれについて4モデルを実行し、計2560区間の処理と閉ROOTの加算・読み戻し検査が完了しました。元の広窓出力は保持しています。局所窓から数学的に無視できる核種tailを除く場合も、その最大混入上限を計算し、purityの低側へ含めています。

補正済みの集計は [summary_local_nominal_002](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/summary_local_nominal_002/summary.json) です。`001`のpurity値そのものは同じですが、4核種×4モデルの比較数を核種ごとの4モデルと混同した集計バグで全点を未評価にしていたため、**支持判定は002を参照**してください。旧001は保持し、target group分離と穴の検出を回帰試験へ追加しました。

| 核種 | 固定応答の連続支持区間（GeV/c） | 単一下限だけで2.0まで支持されるか |
|---|---|---|
| d | 0.80〜1.70 | いいえ |
| t | 1.60〜2.00 | この段階でははい。ただし応答転用による未評価領域の検証を継続 |
| ³He | 0.95〜1.75 | いいえ |
| ⁴He | 0.85〜1.85 | いいえ |

これは**最終推奨cutではありません**。0.05幅で3区間以上連続する支持域を列挙しています。purityが50%未満の領域と、固定応答では残差・支持条件を確認できない領域を区別します。

例えば⁴HeのR=1.005では4モデルのpurityは約0.910〜0.953、R=1.505では約0.767〜0.882となり、広窓linear fitで収量が0になっていた問題は解消しました。一方、R=0.805では選択窓内の形の残差が約19〜22%、R=1.995のGaussian＋linearでは約16%残るため、数値だけで採用しません。

次段階はTOF-tag形状の固定転用ではなく、**ALL自身で目的核種の中心・幅をfree fit**する比較です。Gaussian／log-normalと2背景を維持し、他核種はnuisanceとして扱います。TOF-tagは初期値・帰属の参照に使い、現行のTPC選択窓はまだ動かしません。fit結果の中心移動・幅比・核種bandの連続性を検査してから、下限と校正更新の必要性を判断します。

raw countは保存されたR／dE/dx binから選択境界のbin内一様近似で計算した値です。bin全体が入る場合／境界に触れる場合の上下限も出力し、真の再構成効率やイベント数と同一視しません。

## ALL fitの第3段階：中心・幅をALLから求める比較

[free_local_settings.yaml](../../../../tools/nuclear_pid_cutoff/free_local_settings.yaml)を保存し、全4核種×4モデルの実行を開始しました。出力は `all_<species>_<response>_<background>_nominal_free001/` に分離します。

- 目的核種の中心・幅はALLで自由にfitし、TPC選択窓は凍結したままです。近傍の他核種はnuisanceとして残します。
- 初期値・帰属の参照は、各応答familyで支持されたTOF-tagのlinear背景fitに統一します。元のTAG exponential fitは診断用に保存し、そこの背景境界だけを理由に新しいALL fitを不採用にはしません。
- ALL側ではlinear／exponential背景とGaussian／log-normalの4モデルを比較します。参照TAGの背景選択と、ALLの背景モデルは別の要素です。
- 目的核種の中心移動1参照幅以内・幅比0.5〜2などはfit後の帰属確認であり、最小化へかける制約ではありません。数値共分散・支持域・残差も検査します。
- free fitのパラメータをnative binごとにそのまま本番PIDへ貼り付ける処理は行いません。確定支持域が得られた後に滑らかな関数化の必要性を評価します。

ROOT5数値自己テストは広窓100123、局所固定100127、局所free100154項目がPASSしました。ログ・各段階のソースは `validation_fit_001/` に段階別で保存しています。新しいfree版でも旧固定modeを保持し、旧結果のsnapshotは上書きしていません。

候補数はさらに[独立した640区間の抽出](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/independent_gate_counts_001/gate_counts.tsv)で確認しました。fitが実行不能だった区間では従来の計算経路が候補数の測定まで到達せず、初期値0が残り得るためです。free版の最終集計はこの独立countを必須とし、「purity未評価だが実際の候補数は存在する」状態を保持します。実測に到達したfitとのcount一致、全4種×160区間のcoverage、元cache・選択窓のSHAを検査しています。

## 検証済み区間を使う方針（ユーザー確認済み）

単一の下限値だけではなく、**未検証の小区間を明示的に除外し、検証済み区間を組み合わせて使う**方針をdについて確認しました。続いてt・⁴Heについても同じ方針を確認しました。

dの高R側では、共通fit下端が3.5→3.0 keV/cmへ切り替わるR=1.765→1.775で、linear背景だけがピーク幅を約半分へ縮め、purity推定0.935→0.228、gate残差2.4→28.1%と急変しました。同点のexponential背景はpurity約0.999、残差約3.6%を保ちます。**物理的なpurity急変と断定せず、背景形と窓の頑健性が未確認として扱います。** 1.75〜2.0を除外し、R≥2ではTOF必須の別選択へ戻す設定を可能にします。

t・⁴HeにはGaussian＋exponentialの各1 native点の数値未収束があります。候補数を残すため、下限を大きく引き上げるのではなく、確認できない小区間を除外します。暫定nominal比較でのTPC候補数保持は、単一下限ではt約23%／⁴He約31%、検証済み区間の組合せでは約64%／約52%でした。これは信号効率・CFペア収量ではなく、R=0.4〜2.0の従来TPC候補数に対する比です。

数値精度だけを厳しくした[再試験設定](../../../../tools/nuclear_pid_cutoff/numerical_recovery_settings.yaml)でt・⁴HeのGaussian＋exponential全160区間をそれぞれ再実行しましたが、未収束点は解消せず、他の未収束点が増えたため採用しません。双方が成功した点のpurity差は約10⁻⁶以下でした。都合のよいnative点だけを旧・新結果から抜き出す置換は行っていません。

`speciesPOverQVeto` に全4核種の半開区間 `[min,max)` を明記できるコードを追加し、batch一致build・境界試験を進めています。最良TPC核種の決定後に適用するため、棄却trackの別核種への繰上げは起こりません。新しいブロックがない旧設定は変更しません。最終区間は、実行中の64感度ケースも含めて決定します。

## 高R TOF実装と検証

`FemtoLambdaLegacy`に核種別下限 `speciesMinPOverQ` と全核種の `highRigidityTof` を追加しました。ブロックが存在しない旧設定は従来動作のままです。設定の未知キー・欠落・非有限値・両TOFモードの重複は初期化で拒否します。

新しい[高R TOF検証用mainconf](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_allNucleiTof20261015.yaml)と[核PID設定](../../../../config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_allNucleiTof20261015.yaml)を保存しました。これは**下限未確定の高R条件だけの検証用**であり、最終cutoff採用設定ではありません。既存mainconfの参照先は変更していません。

| 核種 | 高Rで必要なTOF m²窓 | 由来 |
|---|---:|---|
| d | 3.268773 < m² < 3.693147 | 既存中心±1.5σ |
| t | 7.2564275 < m² < 8.2816925 | 既存中心±1.5σ |
| ³He | 1.78285 < m² < 2.06485 | 既存中心±1.5σ |
| ⁴He | 3.0 < m² < 4.0 | 既存の⁴He高R TOF窓 |

m²はraw rigidityから計算した `(m/q)²`、単位は `(GeV/c²)²` です。閾値は全核種 `R ≥ 2.0 GeV/c`。この窓自体の再最適化や高R purityの測定は今回の低R cutoff評価とは別です。d・t・³Heは従来の `R < 2.5` を維持し、⁴Heのみ明示的に有効化していた上限拡張を維持します。新条件を全核種の無制限な高R延長と解釈しません。

SL24y／GCC4.8.5／ROOT5.34/38のbatch一致Singularityで次を確認しました。

| 検証 | 結果 |
|---|---|
| `singularity_make.sh ..._tofExtended.yaml --no-clean -j4` | 成功。common helperとStFemtoMakerを再コンパイル |
| 新しいPID domainテスト | 840項目PASS |
| 既存⁴He高R TOFテスト | 845項目PASS |
| 既存校正テスト（校正＋TOF拡張mainconf） | 8840項目PASS |

新テストは境界、全4核種のTOF、欠損TOF／NaN、最良核種からの再割当て禁止、旧動作、異常設定、実際の`CollectNuclei`、cut前QA／cut後候補、Z=2補正が一度だけであること、ROOT metadataを検査しました。人工データのテスト成功はpurity値の証拠ではありません。

校正テストの初回は、⁴He高R拡張を前提とするテストへ非拡張mainconfを渡したため失敗しました。原因を引数の不一致として記録し、適切なmainconfで再実行してPASSしました。ソースや期待値を変えて隠していません。[検証ログ・fixture](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/validation_domain_001/)へ失敗ログも含めて保存しました。

新しい `..._allNucleiTof20261015.yaml` でもbatch一致buildを実行し、2021 P24iyのローカル復旧cacheの**先頭20イベント**を解析しました。20件すべて再構成、選択Λは1件、`status=0`、STARの警告・エラー・fatalは0でした。ログに全4核種の実際のTOF窓と閾値2.0が出力されています。これは設定・chain・出力の接続確認であり、高R核種のpurityやCF統計を検証する統計ではありません。

- [20イベントROOT](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/allNucleiTof20261015/smoke_20events.root)
- [実行ログ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/allNucleiTof20261015/smoke_20events.log)
- [ROOT読み戻し検査](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/allNucleiTof20261015/readback.log)：zombie／recoveredでなく、raw QA 1648 tracksと全核種TOFのmetadataを確認。

## 2026-10-15終了時点の集計・再開時の引継ぎ

### 完了した評価

Gaussian／log-normal × linear／exponential背景の4モデルを全4核種に適用し、nominal 16実行＋nuisance感度64実行、計80実行が完了しました。感度比較は近傍核種の中心±0.1σ・幅±10%で、目的核種の中心・幅はALLから再fitしています。各核種20条件の比較です。

最終集約は [summary_free_sensitivity_001/summary.json](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/summary_free_sensitivity_001/summary.json)、使用候補範囲と保持率は [adopted_ranges_001/adopted_ranges.json](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/adopted_ranges_001/adopted_ranges.json) に保存しました。`adopted`は研究上の範囲の書き出し名であり、`productionEnabled=false`です。

| 核種 | 下限R（GeV/c） | 除外区間R（GeV/c、半開区間） | 従来TPC候補数に対する保持率 | 使用区間の最小モデル低側purity |
|---|---:|---|---:|---:|
| d | 0.65 | [1.45, 1.50), [1.70, 2.00) | 69.59% | 88.48% |
| t | 0.80 | [1.35, 1.40) | 64.41% | 81.09% |
| ³He | 0.75 | なし | 61.32% | 51.95% |
| ⁴He | 0.70 | [0.85, 0.95), [1.75, 1.80) | 41.95% | 56.73% |

- 数値fitはnative幅0.01 GeV/c、採用判定は幅0.05の区間で信号数・総数を加算してから行いました。各無限小R・centrality・CFペアのpurityを保証するものではありません。
- 保持率の分母は従来TPC選択のR=0.4〜2.0の候補track数です。真の信号効率、イベント数、CFペア収量、高R TOFの保持率ではありません。
- モデル低側50%以上は、今回比較したモデル・感度範囲での条件です。全系統誤差・信頼水準の検証完了とはしていません。
- 感度比較を加えたため、先のnominal暫定値からdの高R除外と⁴Heの除外域・保持率が変わっています。例えば⁴Heの保持率は約52%ではなく、最終比較で41.95%です。
- 数値共分散・応答支持などを確認できない除外区間を「低purityが測定された区間」とは扱いません。⁴HeのR=0.65〜0.70では一方でモデル低側が約41.4%となり、50%基準を満たしません。
- PIDの中心・幅の関数係数は更新していません。ALL fitとのずれは残りますが、固定した実TPC選択窓を評価した研究なので、関数を変える場合は選択窓・best species・purityの再評価が必要です。

### 完了した実装検証

`speciesPOverQVeto`まで含むbatch一致buildと拡張PID domain試験は、絶対出力パス・相対出力パスでそれぞれ**1359 checks PASS**でした。[validation_domain_veto_002](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/validation_domain_veto_002/run.log)に保存しています。初回001の相対出力パス二重解決によるテストharness失敗は保持し、harnessの修正後に再検証しました。Python回帰試験は17件PASSです。これらは実装・集約の検証であり、物理purityの独立検証ではありません。

### 中断時に残した作業（以下は2026-10-15時点の記録）

1. 上のJSONを正として、核種別下限・除外域・全4核種R≥2のTOF条件を組み合わせた**新しい独立YAML・mainconf**を作成します。既存本番mainconf・校正係数を上書きしません。
2. 最終mainconfのanalysis情報は既存の2021 P24iy・canonical出力先を参照し、試験出力はその下の`tests/purity50_20261015/`などへ分離します。2019年は使用しません。
3. 最終mainconfで`singularity_make.sh ... --no-clean -j2`を実行し、既存ローカル2021 cacheの先頭20イベントで新条件・ROOT metadata・raw QAを確認します。高R TOFだけの20イベント試験と混同しません。
4. 保存済み集計から最終purity／COLZ／保持数／応答比較PDFを作成します。予定の`views_final_001/`は終了時点で未作成です。新しいfitは原則不要です。
5. 結果図を確認してTodo・出力README・運用索引を更新し、最終プロファイルの採用状態を明記します。全4核種TOFのみの既存新プロファイルを、低R変更の切り戻し先として残します。

ユーザーの終了指示後は新規fit・描画・イベント解析を開始していません。元ROOT・旧成果物・既存設定・未commitの変更は保持し、job投入、hadd、メール、Git commit/pushは行っていません。

## 2026-10-17（日本時間）の再開・完了結果

### 最終設定と切り戻し

[最終核PID YAML](../../../../config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_20261015.yaml)へ、上の最終表の核種別下限と除外域をそのまま設定しました。R<2 GeV/cはTOFを要求せず、R≥2では全4核種のTOFを要求します。下限はinclusive、vetoは`[min,max)`、TOF窓と従来の最大Rはstrictです。PID中心・幅・nσ・best-species選択は変更していません。

[最終mainconf](../../../../config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle_purity50_20261015.yaml)はcanonicalな2021 P24iyのanalysis情報を参照します。出力先は従来の`rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/`で、試験出力だけを`tests/purity50_20261015/`へ分離しました。ファイル名の`20261015`は元studyの日付で、作成・接続検証日は2026-10-17 JSTです。

下限・vetoを外して全核種の高R TOFだけを残すには、既存の`..._allNucleiTof20261015.yaml`を指定できます。その検証用mainconfは別のanalysis名・出力先なので確認が必要です。旧本番設定、旧係数、元ROOTは保持しています。

### 同一mainconfでの検証

| 検証 | 結果 |
|---|---|
| `singularity_make.sh 最終mainconf --no-clean -j2` | SL24y／sl73_x8664_gcc485／ROOT5.34/38、exit 0 |
| 2021ローカルcacheの先頭20イベント | 入力20／再構成20／選択Λ1、status 0 |
| STAR警告／エラー／fatal | すべて0 |
| 閉ROOT読み戻し | 正常、raw QA 1648 tracks、最終下限・veto・TOFのmetadataを確認 |
| 同一mainconfを入口にしたdomain回帰 | 1359 checks PASS |
| Python集約・範囲出力回帰 | 17 tests PASS |

[検証README](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/purity50_20261015/README.md)から全ログ・ROOT・source snapshot・libraryを含むSHAへ辿れます。これは接続・実装検証であり、20イベントから核purityやCF効率を測定したものではありません。

独立監査では80実行のcoverage、1062件のSHA、範囲・保持率・モデル低側purityを照合し、全て一致しました。旧本番mainconf・旧核PID・係数header・校正YAMLのbaseline不変も確認しています。[監査READMEと再現スクリプト](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/validation_resume_audit_20261017_001/README.md)、[監査JSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/validation_resume_audit_20261017_001/report.json)を保存しました。元10GB統合ROOTの再ハッシュは行っていません。

### 最終PDF

| PDF | 内容 |
|---|---|
| [purity・cutoffの全体図](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_final_001/purity_cutoff_overview.pdf) | 3ページ。モデル包絡と使用／除外範囲、広いALL COLZへの現行PID窓の重ね、候補数。 |
| [中心・幅の比較](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_final_001/response_parameters.pdf) | 5ページ。現行関数とnominal ALL free fit、COLZへの重ね。 |
| [dの代表fit](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_final_001/fit_spectra_d.pdf) | 12ページ。下限近傍と高R未評価域。 |
| [tの代表fit](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_final_001/fit_spectra_t.pdf) | 12ページ。下限近傍と除外区間。 |
| [³Heの代表fit](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_final_001/fit_spectra_he3.pdf) | 12ページ。下限近傍とTPC-only上端。 |
| [⁴Heの代表fit](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_final_001/fit_spectra_he4.pdf) | 16ページ。下限の前後、除外区間と高R端。 |

代表fitは全dE/dx=0〜120 keV/cmと局所fit窓、linear／log、4モデルを保存済みcountのまま表示しました。fitを再実行したり、表示の規格化を変えたりしていません。全6 PDFのページ数・全頁読み戻しを検査し、主要図を目視確認しました。[図の説明・再現手順](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_final_001/README.md)と入力／出力SHAを保存しています。

中心・幅の関数は今回更新していません。例えばtでは現行幅がALL Gaussian fitより約14%狭いなどの差はありますが、評価した実選択窓自体を変更すると今回のpurityをそのまま適用できません。現行関数が全域で完全に正しいと結論せず、今回の採用候補範囲は現行窓のままとします。モデル包絡を信頼区間・全系統誤差の完了とも扱いません。

新規batch投入、hadd、メール、Git commit/pushは行っていません。旧結果の削除・上書きも行っていません。
