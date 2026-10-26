# 新しいヘリウム校正のTPC選択窓QA

[Purity比較PDF](purity_calibrated_windows.pdf)は元校正±2、新校正±2、新校正best-speciesの3選択を比較します。背景関数は保存済みlinear/exponentialを同じ選択窓で積分しています。

205のchosen component-evidence knotsを両核種へ適用し、6品質通過／199参考fitを区別します。範囲端は支持half-bin内のみclamp、未支持gapや範囲外はlegacy校正です。0.405は採用していません。

charge-inclusive・中央度積分のpre-PID分布のTPCモデル推定で、TOF後やペア別の確定純度ではありません。保存されていない共分散から新しい統計誤差を作らず、CFへの純度補正も自動適用しません。

`purity_windows.json/csv`にmodel別の純度、信号数、元gate比、窓の区間・支持範囲を保存します。
