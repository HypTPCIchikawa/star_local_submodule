# TOFタグを手掛かりにしたd・t・³He・⁴HeのPID再校正計画

## 依頼と対象

保存し直された[Todo](STAR_Todo_20261014.txt)の最新版に従う。
無選別dE/dxだけでの4He分解を出発点とせず、TOF窓で選んだ分布から各核種の帯を調べる。
その中心・幅・支持域を使って未選別分布を再fitし、PID窓のpurityを評価する。
不明点は対話で相談する。各段階の結果と図の場所を[報告](report_tof_assisted_nuclear_pid.md)へ更新する。

- 入力: `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root`。
- 実行環境: `config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml`、SL24y／ROOT5。
- 未選別QA: `true/hDedxP`。
- TOFタグ: `true/hDedxP_d_m2`、`true/hDedxP_t_m2`、`true/hDedxP_3He_m2`。
- 3.85 GeVは2021 P24iyのみ。2019出力・13.5 GeVデータは今回の入力としない。
- 元ROOT・既存fit・旧PDF・現在の解析PIDは保護する。本投入・Git整理・pushはしない。

## 入力の意味と監査

剛性は`R=|p/q|`。TOF量は`R²(1/β²−1)=(m/|q|)²`で、単位電荷を仮定した物理m²ではない。
dと⁴Heのm/|q|は近いため、dタグはd／⁴He互換のタグである。
同じdタグの低いdE/dxのd帯と高いdE/dxの⁴He候補を分けて調べる。
タグ名だけで核種の純度やピークの帰属を確定しない。

元ヒストは600×1200、R=0〜6で幅0.01 GeV/c、dE/dx=0〜120で幅0.1 keV/cm。
実ROOTの軸・有限非負count・flow・Entriesとstored countの差・Sumw2の有無を記録する。
Sumw2がないため、保存countから得たPoisson型誤差を使うことを明記する。
TH2Fの丸めを理由に元のcountを整数へ書き換えない。

TOF3タグとd・tタグは互いに非重複のm²/q²窓。全binで`TAGd+TAGt+TAG3 <= ALL`を照合する。
TAGはALLの部分集合なので、ALLとTAGを独立データとして尤度に重複加算しない。
基本はTAGでshapeを推定→ALLで条件付fitする2段階とし、TAGの不確かさと転用biasを分ける。
必要ならALLから3タグを差し引いた非重複remaining分布で転用の整合を確認する。

## 段階1: 可用統計とピークの位置

1. 新しい入力cache・来歴を作り、4本の未切断ヒストとR射影を保存する。
2. 全native binを使ったCOLZと、RごとのdE/dxのlog／linearを出す。
   低dE/dx側を隠さず、全保存軸0〜120と核種の補助図を併記する。
3. R分割はnative bin端と一致した非重複区間だけとする。0.005刻みの独立fitは作らない。
   粗い0.05／0.10で全域を見渡した後、十分な統計がある場合は0.01まで細かくする。
   低統計の隣接bin合併は明示し、重複windowを独立測定点にしない。
4. R≈0.8〜0.9など帯が分離する領域を初期確認し、低Rへ順次調べる。
   最強の全域ピークを原子核ピークとは呼ばない。

研究範囲は保存軸と元の関数の範囲を調べてYAMLへ明記する。
低Rのfit可能下限に加え、高Rの識別・混合による支持域上限も記録する。

## 段階2: TOFタグでのshape fitと関数化

- Gaussian信号とpositiveなpol1／exp背景を基本に、必要な近傍核種成分を比較する。
  旧関数は初期値と照合用で、旧中心・幅の近傍へ強く固定して「検証」としない。
- d／t／³HeのTAG分布で明瞭な帯を先に調べ、そのshapeをdタグの高い側の分解にも活用する。
  ⁴Heが未識別のbinに旧ALL free fitの⁴Heラベルをtruthとして押し込まない。
- bin積分Gaussian、固定したmodel比較窓、非負yield・正値幅・有限値を使う。
  native R幅／合併R幅による投影幅とintrinsic幅を区別し、最終関数のR-bin平均で確認する。
- まず現行と同じ`log10(R)`の中心3次式・幅2次式を比較する。
  ユーザーは必要に応じた低R／高Rの連続接続関数の比較を許可した。
  正値幅・連続性・支持域・残差を監査し、細かいfree点のガタガタをそのまま本番式にしない。
- 収束、成分yieldと不確かさ、境界・zero-yield、残差、背景／初期値／binning感度を保存する。
  支持域の判定値・model・R分割・計算数はYAMLに置く。
- 下限は「式を描ける最小R」ではなく、ピーク／成分が識別できてmodel変更に安定する最低R。
  軸0〜120からピークが欠ける・統計不足の区間は未支持であり、核種の不存在とはしない。

## 段階3: 未選別ALLの再fitとpurity

TAGで得た関数をALLの固定shapeまたは不確かさつき候補として用いる。
TAG・ALLの受理の違いがshapeへ影響する可能性を、R域別の残差・TAG間／remainingとの比較で確認する。
pol1／expと必要な近傍Gaussianを比較し、一般粒子の低dE/dxピークを少数の原子核Gaussianへ押し込まない。
全軸表示と、各成分を推定するfit窓を区別する。

各核種について`P=S/(S+B)`を同じ選択窓で積分する。
現在の±2σとbest-species選択を基準にし、候補の単独Gaussian±2σ窓の純度は別列にする。
TOF併用後の最終候補やペアCFのpurityと混同しない。
TAG shapeの不確かさ・背景model差・TAG転用の仮定を保存し、条件付きfit誤差だけを最終誤差にしない。
未支持binの値を確定purityとしてつなげない。

## 出力とコード更新

- 支援コード: `tools/nuclear_pid_tof/`（旧`tools/he4_purity/`は保持）。
- 出力: `rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/<study>/`。
- 入力cache・TAG射影fit・中心／幅／COLZ重ね・支持域表・ALL fit・R vs purityを整理して保存する。
- 条件・入力fingerprint・source snapshot・設定・ROOT/JSON/CSV・検証記録を保存する。
- 検証済みパラメータは新しい校正YAMLと設定profileへ保存し、元へ戻せる構成にする。
  データが支持しないRへの外挿や旧式への急な切替を検証済み校正と呼ばない。
  支持域外の本番扱いは、結果を示してユーザーと相談する。検証前に現在のPIDを置換しない。
- C++やconfig構造体を変更する場合はSL24yのSingularity環境で再buildし、旧解析の回帰試験を行う。
- 新しい入口・設定・図のリンクをREADME／REFERENCEと今回Todoへ更新する。
