# 高R TOF選択後purity・全p/q fit集：実施記録

## 現時点の状態

2026-10-17（日本時間）の依頼です。高R評価・低R統合PDF・全R fit集・新旧比較、本番5入口のcutoffと検証を完了しました。追加指令のk*–核PID TH2／PID前対照／8D sparseも実装し、実イベント・合成入力・mixing不変・容量・toy haddの検証まで全てPASSしました。

## 入力監査で確認した点

### 追加タスク7：k*別の原子核PID診断

実行入口は anaFemtoLambda_4He.C＋KFParticle providerで、物理処理は共通Makerへ置きます。通常選択後の⁴He／³He、および既存³Heを⁴He質量へ写す対照ペアについて、k*–nσとk*–raw p/qを追加します。SE/ME・cent9・Λのsignal/SBPos/SBNeg・signal低側/高側を既存と揃えます。

任意のnuclearPidPairQaHist参照と専用HistManager、THnSparseDを実装しました。通常4He入口と4He highpurityだけで有効化し、旧ヒスト・cut・物理候補・mixing保存判定は変えません。binningはYAMLへ明示。raw p/qはZ=2補正済み物理運動量を2で割り、³HeAs⁴Heも元³Heのpullを保持します。

選択後nσ分布は±2σとbest-speciesですでに切断され、fit外側の背景を見られません。ユーザー回答により、通常選択後＋PID前の両方と、(k*,p/q,nσ)の疎な同時分布保存が正式承認されました。物理候補storeへQAを追加するとmixing poolのイベント採否が変わり得るため、別QA store/payloadで分離します。sparseにはcent9・SE/ME・元核種・排他的Λ質量区分・PID段階も保持し、全signalは低側＋高側の和で重複なく取り出せるようにします。実装・容量・時間・回帰試験の完了まで、承認と完成を区別して記録します。

既存10GB統合ROOTには今回の新ヒストはなく、今後再解析した出力から利用できます。追加後はtoyで射影和・species pull・Z補正・mass halves、QA off/onで既存全ヒスト不変を検証し、batch一致buildと短い実イベント読戻しを行います。大規模job投入や既存ROOT上書きは行いません。

実装後の[検証索引](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/nuclear_pid_pair_qa_20261017_001/README.md)：ROOT5 build成功、同一20eventのQA off/onで旧690hist・2tree完全一致、新1080 TH2と8D joint 87 Entries保存、RSS約314→522 MiB。人工100万疎binの保存・読戻しと、toy-only ROOT5 hadd（全内容・誤差²の加算）もPASS。合成Sparse/pair/poolはtoys_005で245,073,676項目PASS、pool_nonzero_001で通常SE=1・ME=3という非ゼロCFも確認しました。87 Entriesはイベント数や独立核数でなく、stage・核種/質量仮説別の診断pair数です。新histの[名前と意味](../../../../mdfiles/femto/implementation/femto_lambda_histograms_20260930.md)も更新済みです。

今回のhaddは意図的な同一入力×2の加算試験だけです。本番hadd・job投入はしていません。既存自動merge監査はTH1系に限定されているため、新Sparseの監査追加と、本番入力1〜数ファイルでの容量・workerメモリ確認を次回Todoへ残します。

追加の完成品：[高R fit集24冊840頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/highrigidity_tof_20261017_001/fit_atlas_001/README.md)、[新旧PID比較23頁](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/pid_parameter_comparison_final_20261017_001/raster_report_001/new_old_pid_final.pdf)。高Rはd/t/³He各50 native区間と⁴He全60 coarse区間、全4モデル・未支持点も掲載し、全冊読戻し・SHAを確認しました。新旧比較は設定と実fitを区別し、全低Rと採用低Rを別表示。⁴Heの旧高R fallback幅はR≈4.5でfit幅の約1.8〜2倍ですが、有限R投影幅なので係数を自動置換していません。

- 2021 P24iy部分統合ROOT由来のnative cacheにはALLとd／t／³HeのTOFタグがあります。R=0〜6 GeV/c・幅0.01、dE/dx=0〜120 keV/cm・幅0.1です。
- d・t・³Heのタグ窓は現行高R TOF窓と一致します。TPC・best-species選択前の分布から、現行TPC窓内のpurityを評価します。
- 元ROOTの282ヒストグラムと提出時runtimeを確認した結果、⁴Heの3<m²<4に対応する専用分布やdE/dx×m²の同時分布はありません。dタグの3.268773<m²<3.693147とは異なるため、同じ条件の本番⁴He purityとして代用できません。
- 最初は狭窓の参考値として扱う方針でしたが、ユーザーの追加指示により、**今回のstudyでは⁴Heもdと同じ3.268773<m²<3.693147を選択条件とします。** dと⁴Heは(m/q)²で重なるため、TOF窓だけで区別せず、窓内dE/dxで分離してpurityを求めます。fit中は旧設定を保持し、後続の本番適用指示に従って新profileにもd同一窓を採用しました。旧3<m²<4のpurityは未評価のままです。
- 高Rの選択順序は元TPC best species→核種別下限・veto→その核種のTOF→最大Rです。TOF不合格を理由に次点へ繰り上げません。
- d・t・³HeはR<2.5、⁴Heだけ明示的高R拡張です。保存軸が広いことと、既存fitが支持された範囲を区別します。
- 既存高R TOFタグfitは2〜2.5のnative全50点で数値収束していますが、旧χ²基準では形状支持を得ていません。統計不足やfit不可能と同一視せず、局所形状・背景・実TPC窓の残差を検証します。旧結果を無条件に確定purityへ流用しません。

## 全p/q fit集

概要PDFは0.05 GeV/c幅のpurityですが、元fitは0.01幅です。保存済みnative fitを全区間表示するPDFを核種別・R範囲で分冊します。未収束・不採用点も診断として掲載し、曲線の存在を支持済みと同一視しません。全dE/dx=0〜120と局所fit窓、linear／log、4モデルを保存countのまま表示します。

## 完了した最終数値と統合図

高R native 0.01幅の1800fitに加え、⁴HeはR=2〜5を全域0.05幅で240fitしました。いずれもGaussian／lognormal×linear／exponentialの4モデルで、都合のよいnative点の抜粋やpurityの単純平均はしていません。d/t/³HeはnativeのSとTを加算して0.05へ集約し、⁴Heは全域coarseの独立fitを採用しました。coarseの幅は有限R投影の幅で、検出器固有の校正幅とは区別します。

| 核種 | 検証済み高R範囲（GeV/c） | 範囲内の最小モデル低側purity | TOF m²窓（strict） |
|---|---|---:|---|
| d | 2.0 ≤ R < 2.5 | 99.43% | 3.268773 < m² < 3.693147 |
| t | 2.0 ≤ R < 2.5 | 99.51% | 7.2564275 < m² < 8.2816925 |
| ³He | 2.0 ≤ R < 2.5 | 98.68% | 1.78285 < m² < 2.06485 |
| ⁴He | 2.0 ≤ R < 3.2 | 98.52% | dと同じ |

⁴Heの3.2≤R<5は低purityが確定したのではなく、無視できない寄与を持つnuisance（他核種・背景）成分の応答が未検証であるためです。さらに高Rでは低統計・幅帰属の問題もあります。0.05へ再fitすると統計残差は改善しましたが、未検証nuisanceの基準は緩めていません。支持された高R範囲では全核種ともモデル低側が50%を超えています。

⁴Heの高R条件後候補は、2≤R<5の3,562,460.13から支持域の3,376,788.51へ、94.79%を残します。これはinclusive QAの候補countで、真の信号効率でも低Rを含む全収量の保持率でもありません。小数は保存native-binの選択窓重なりによるcountです。

- [統合PDF（6ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/views_tpc_tof_final_20261017_001/purity_tpc_tof_combined.pdf)：低R＋高R、TOF後高R拡大、条件説明、元の低R3ページ。緑線は本番採用区間、斜線は除外・未評価、赤白抜きは未支持fitです。
- [最終数値・quality・残差](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/highrigidity_tof_20261017_001/summary_final_001/summary.json)：全モデルの検査理由を保持。
- [全binの採用判定](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/production_ranges_20261017_001/bin_decisions.tsv)と[本番採用区間](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/production_ranges_20261017_001/README.md)。

低Rは4モデル×5応答変動の20条件、高Rは4 nominalモデルです。帯はモデル感度で信頼区間ではありません。全て2021 P24iy部分統合ROOT由来のtrack QAで、全23Bイベントを確認したものでも、centrality・mT・CFペア条件付きpurityでもありません。

## 完了：低R全p/qのfit集

[正式索引](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/fit_atlas_native_20261017_003/README.md)に全32冊・2560ページを保存しました。4核種それぞれ160 native区間×4モデルを全て表示し、各核種8冊・640ページ、PDF合計619,296,231 bytesです。数値fit・応答検査の状態を各パネルへ明記し、未収束・不採用区間も掲載しています。全冊の閉PDF読み戻し・頁数・SHA・保存fit入力の不変性を検査しました。新fitやPID変更はしていません。

## 追加指令：新旧選択パラメータ比較と本番cutoff

新旧dE/dxの中心・幅・選択帯と、TOF m²窓を比較表示します。TOF窓は今回中心・幅を再fitした値ではなく、設定の比較として明示します。高RのdE/dx応答は実TOF選択後の分布でfitし直し、校正係数へそのまま貼り付けたものとは扱いません。

本番への適用範囲は追加確認済みで、全4核種の入口を共通校正・cutoffへ揃えます。これまではcutoff機能と独立設定までの実装だったため、本番入口の接続が必要です。旧d・t・³He専用入口はlegacy校正で、共通校正へ移すとTPC最良核種の競合判定も変わります。この差を比較・検証します。fit中は基準となる旧設定を凍結して保護し、旧設定への切り戻しも残します。本番jobは投入しません。

## 本番5入口の採用・検証（後続の最終状態）

canonical⁴He／⁴He highpurity／d・t・³He highpurityの5入口を[共通本番YAML](../../../../config/cuts/nuclearid/nuclearid_femtoLambda_auau3p85_purity50_production_20261017.yaml)へ接続しました。cutoff切替時の他参照は不変で、後続TH2タスクで⁴Heの2入口だけQA参照を追加しました。対応するmainconf名末尾の _prePurity20261017.yaml で旧設定へ戻せ、参照が切替直前と同値であることを独立確認済みです。旧係数・YAML・study profileは保持しています。

d/t/³He専用入口はlegacyから共通校正へ変更しました。d/t自身の応答が同じでも競合³He/⁴Heのpullが変わり、候補選択が変わり得ます。再fitの投影幅は校正係数へ直接採用せず、校正有効域外は既存fallbackを維持します。

TPC winnerを先に決め、下限・veto・TOFを適用し、落ちても次点核種へ繰り上げません。全体上限R<5、d/t/³Heの[2.5,5)と⁴Heの[3.2,5)をvetoし、⁴Heの特別な上限超え許可はfalseにしました。

canonical入口のbatch一致buildはexit 0。旧fixture1359＋実5入口各69702、計349869の境界試験がPASSしました。[境界試験](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/production_purity50_20261017_001/domains_001/completed.json)と[独立設定照合](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/production_purity50_20261017_001/profile_audit/completed.json)を参照。productionEnabled:falseはテスト自体の設定変更をしない意味で、実mainconfが未接続という意味ではありません。最新cutoffツールのPython36試験もPASS。

全5入口で同じ2021先頭20入力／20再構成／選択Λ=1、閉ROOTのraw核QA Entries=1648とdomain metadata一致を確認しました。[全実行・読戻しログ](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/production_purity50_20261017_001/smokes_001/README.md)を参照。接続検証であり、purity・信号効率の再測定ではありません。元ROOTは不変更で、job投入・hadd・メール・git commit/pushは実施していません。

高R専用adapterは既存のlocal数学を再利用し、低R80実行は変更していません。ROOT5/GCC4.8.5でのcompile・数値自己試験100009項目と、高R Python7試験がPASSしました。最終stageでlow側31試験・high側7試験と元study設定／数学ソースSHA不変を確認しています。予備fitだけではなく、上記の全fit完了を最終判定に使いました。
