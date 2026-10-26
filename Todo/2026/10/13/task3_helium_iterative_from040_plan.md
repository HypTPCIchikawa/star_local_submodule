# Task 3方針変更: 0.40開始の区間別fitから反復して共通関数を得る

## 今回の依頼と確認

`refined_lowrigidity_from040/views/he4_low040_components.pdf`は現存する。
旧`refined_lowrigidity`は低剛性の探索を0.50から始める設定であり、0.40〜0.50の出発点にはしない。
ユーザーの依頼に従い、0.40開始の区間別fitを初期値にする。
従来式から直接始めた未収束global prototypeは、今回の初期値・最終結果には使わない。

保存された0.40開始fitにも旧nσ相対の探索範囲・幅・近傍核種の制約がある。
「完全に無制約だった」とは記録せず、今回の最初のfree phaseではそれらのlegacy相対制約を外す。
Gaussianの幅正値・非負収量・有限値と、宣言した検出器軸内の数値supportは維持する。

対象は3.85 GeV・2021 P24iy。元ROOT・旧fit/PDF・本番PID・mainconf・共有ヘッダを保持し、
`tools/he4_purity/iterative_logpoly/`と新しいstudyへ追加する。本投入・git操作は行わない。

## 固定する入力とモデル

- 保存済み`from040/he4_fits.root`のraw native射影とJSON/TF1をSHA照合して使う。
  R幅0.01、0.40〜2.50の210区間。元dE/dx幅0.1をfit用0.5 keV/cmにまとめ、counts/error/flowを保持する。
- 初期Gaussian中心・幅・全収量は旧TF1からraw keV/cmへ厳密変換する。
  不採用sliceの旧modelを初期値に使う場合も、不採用だったことを記録する。
- T・³He・⁴He Gaussianと正値smooth背景を用いる。Tの窓外shapeなど情報のない自由度は識別不能と記録し、
  fitted/固定/非使用を明示する。無情報のHesseを「良いfit」と扱わない。
- fit窓は全phaseで同じ。T・³Heのピークが見える広さを設定し、保存軸0〜120 keV/cmでclipする。
  0〜120全域には他の核種もあるため、表示全域と実fit範囲を区別する。
- GaussianはdE/dx bin積分で予測する。最後の共通関数fitではR-bin内変化も平均し、
  区間別fitの投影幅と関数のintrinsic幅を混同・二重加算しない。
- Gaussian templateは固定fit窓で正規化し、非負の窓内収量をfitする。
  旧fitの全Gaussian収量capは引き継がない。窓外への全収量換算を物理収量とは扱わない。
  各sliceのlinear／exponential背景の種別は旧採用modelから引き継ぎ、全phaseで固定する。

## 反復の手順

1. **区間別free phase**: 旧0.40開始fitを初期値として、Rごとの中心・幅を自由にfitする。
   legacy中心への固定・legacy±nσのhard bounds・隣接Rへの拘束を最初からかけない。
   Gaussian成分の入れ替わり、zero-yield、sub-bin幅、窓外成分の非識別を記録する。
2. **関数の初期化**: free phaseから、³He/⁴Heの中心を`log10(R)`の3次式、幅を2次式へfitする。
   R中心化・正規化した基底で数値条件を整え、最終係数は解析コードの原形式へ変換する。
   正値幅は区間端・二次式の頂点で確認する。非負Bernstein係数に限定して関数familyを変えない。
3. **反復phase**: 関数を初期値に区間別fitを再実行し、関数からのずれへsoftな平滑化penaltyを段階的に加える。
   hardなlegacy相対制約へは戻さない。penalty強度・情報量の定義・減衰更新はYAMLへ明示する。
   各stage内でraw data devianceとpenaltyを分け、非改善stepは保持・減衰する。
   局所phaseはR中心で評価し、この段階の幅は投影幅である。penalty情報量と中心・幅のscaleは
   free phaseから固定し、目的関数評価中には更新しない。中心と幅の差の二乗を用い、同じ重みのQRで関数を更新する。
4. **最終共通関数**: 平滑化された係数を初期値に、R-bin平均を含むraw分布の共通関数fitを行い、幅も再最適化する。
   規格化・背景・Tは適切なnuisanceとして扱う。各Rの独立free中心・幅を最終PID式として保存しない。
   global14係数をcold startで一度に探すのではなく、区間別fitを追ってつなぐ。

平滑化に同じデータを再利用するため、penaltyやfit共分散を独立な測定精度・核種の証明としない。
failed binで突然legacy値へ戻す補間は作らず、未支持区間はreferenceとして表示する。

## 出力・図

- 各phaseの入力・係数・局所model・penalty/raw deviance・収束・成分窓内収量をROOT/JSON/TSVへ保存する。
- 全210区間のnσ linearと**元dE/dx linear**を、それぞれ4区間/ページで保存する。
- 低Rではnσ log/linearとdE/dx log/linearの4面を同じ範囲で表示し、T・³He・⁴He・背景を重ねる。
- nσ軸は元のlegacy anchorで固定し、同じbin count/errorを表示する。座標変換でy値のJacobianを加えない。
- free→反復→最終の曲線・幅・残差・収束履歴、生p/q×dE/dx COLZへの重ね図をまとめる。
  model幅と係数の不確かさを混同しない。poor fit・境界・局所ピークのない肩を明示する。
- 関数・domain・来歴は新しいstudy-only YAMLに保存する。本番PIDへの自動適用は行わない。

## 検証・採用判断

元のcounts/error/flow・原本hash不変、Gaussian bin積分・成分和、raw/nσ表示等価性、
係数基底の往復・幅正値、zero-yield/label交換・反復更新を数値テストする。
ROOT処理はSL24y／ROOT5／GCC4.8.5で行う。

最終の滑らかさだけではなく、同じfit範囲のfree phaseに対するraw devianceの変化、
各Rの残差・未識別成分、初期値・平滑化強度・背景・R平均の感度を確認する。
収束しない場合やfit劣化が大きい場合は理由を示し、良好な関数が得られたとは主張しない。
Todoの報告・索引とtools/REFERENCEを途中から更新する。

## 全域fit後の限定した追加診断

全域の反復で低Rの成分が崩れた場合、数値精度・背景初期値・T形状固定・係数初期値を切り分ける。
これらで解決しない場合の最後のdomain診断は、同じfree入力のR=0.40〜0.70の30区間だけで行う。
関数family・fit窓・観測count・背景modelは維持し、同じ30区間のfree Dと比べる。
T可変とfree T形状固定を別出力へ保存する。限定domainを本番の範囲cutやpiecewise式と扱わない。
保存済みTOF³He／DタグもTPCとは独立した対照として用いるが、Dタグを純粋な⁴Heとは呼ばない。
有効な共通式が得られなければ、本番を変更せず、model・関数形を変更する前に結果と未解決点を説明する。

## 追加依頼: 表示全域と生COLZの確認

元dE/dxを高い側だけに制限せず、保存軸0〜120 keV/cmの全normal binを表示する。
対応するnσも同じnative bin・count・errorを用いる。下側の大きな通常粒子ピークでlinear表示の核種が
埋もれる場合に備え、全域図を保持した上で低R10区間の45〜120 keV/cm補助図を併記する。
表示範囲を広げてもfit窓は変えず、窓外は観測データだけを表示する。

保存済み`true/hDedxP`のCOLZへ、最新の共通式候補の³He／⁴He中心と±1σ・±2σを
TF1の`Draw("same")`で重ねる。係数を同じ原式・raw keV/cmで評価し、支持domain0.40〜0.70の外へは外挿しない。
R-bin幅を含む区間別free幅と、最終共通fitのintrinsic幅を混同しない。
元ROOT・旧図・fit結果・本番PIDを変更せず、新しい表示専用出力へPDF/PNGと数値監査を保存する。
未収束・model依存・未採用を図中にも記載する。視覚的な重なりだけで核種同定やpurityを確定しない。
