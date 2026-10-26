# ⁴Heと³HeのPID校正と背景CF入力の実装計画

3.85 GeVの2021 P24iyで得た0.40開始のfitを、Λ–⁴He解析のヘリウムPIDへ明示的に導入する。新しい選択窓に対応するPurity対p/qを作り、同じ出力ROOTへΛ–³HeのSE/MEスペクトルを追加する。既存ROOTと以前の設定は保持し、13.5 GeV・Helix解析・Phi解析の共通校正は変更しない。

## PID校正の導入

fitで得た⁴He信号と³He背景Gaussianの中心・幅を、エネルギー別の校正YAMLへ保存する。大域的な多項式候補を無条件に適用せず、利用できるfit点の間を補間する。校正domain、端点の扱い、範囲外で従来校正へ戻す方針、採用・参考の区別を設定と来歴に明記する。

最初の0.40–0.41の⁴Heは平均境界のため採用しない。³Heもその成分自身の境界・covariance・収束を確認して点を選ぶ。多くの区間はモデル品質外の参考fitであり、作業用校正としての限界を残す。今回導入するのは中心とGaussian幅で、TPCの許容nSigma閾値は別の既存設定を維持する。

新しい校正providerは`FemtoLambdaLegacy`が所有し、対象のnuclearid YAMLが明示的に有効化した場合だけ使用する。d・tは従来の中心・幅を使い、³He・⁴Heの新しいpullを同じbest-species競合とPID QAへ渡す。共有の`NuclearIdDeDxVsMom.h`は書き換えず、他の解析への暗黙の波及を避ける。ROOTへ校正設定・出典・補間方針も記録する。

3p85の現行Λ–⁴He mainconfには専用のnuclearid・maker・hist設定を参照させる。旧参照のmainconfを別名で保持し、従来条件へ戻せるようにする。提出済みjobのsnapshotや過去のpilot設定は更新しない。

## Purity対p/qの更新

元校正±2窓と新しい中心・幅の±2窓を別曲線で比較する。可能な範囲で新しいbest-species選択に対応するTPCモデル推定も示す。範囲端・不採用fit・背景関数依存を表示し、旧gateの純度を新しいgateの純度として流用しない。

保存ヒストはcharge-inclusive・中央度全体のpre-PID分布である。TOF後、Λとのペア、k*やmT別の純度を確定したことにはならず、CFの差し引き係数へ自動適用しない。以前のPurity PDFは保持して別のサブフォルダへ出力する。

## ³Heの背景CF入力

明示的な追加flagを持つΛ–⁴He profileで、`lambda,he4,he3`の3 speciesと核種ごとの4 mass channelを許可する。新flagがない既存profileの単一核種・4 channelの仕様は維持する。中央度、Λの±3σとsideband、低側／高側、SE/ME、k*–mT、質量–k*、アクセプタンスを⁴Heと同じ方法で³Heにも記録する。

³Heを⁴Heと誤識別したペアは、³Heの実質量で求めたk*とは異なる。そのため、同じ³He採用ペアについて通常の`lambda_he3`と、エネルギーだけを⁴He質量で再計算した`lambda_he3_as_he4`の2 namespaceを保存する。運動量のZ=2補正は一度だけとし、mapped controlのためにmixingを再実行したりペアを追加抽選したりしない。

MakerはSE/MEを保存するだけで、最終CFは統合後に作る。³He core選択と⁴He窓へ漏れる³He tailは母集団が同一とは限らないため、control CFを無補正で背景CFと断定しない。今回は差し引きに必要な入力を整備し、係数や転送補正は別の検証対象として記録する。

## 検証と記録

設定・補間・範囲外policy・旧条件への復帰を人工例で検証する。3核種を含む出力の全hist名・軸・SE/ME・centrality・signalLow+High一致・mT射影・mapped質量仮説と運動量を確認する。StMaker変更後はSingularityのSL7/ROOT5環境で再buildし、2021 P24iyの小さい実イベント試験を行う。

Task 3の計画・最新結果・変更箇所・検証・残るCF補正課題を本TodoのMDへ継続して更新する。ヒストグラム一覧、Femto README、実行手順も同期する。job投入・git pull/push・既存ROOTの書換えは今回行わない。
