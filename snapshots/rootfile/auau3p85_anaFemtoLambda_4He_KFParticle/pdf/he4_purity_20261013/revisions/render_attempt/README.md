# ⁴HeのTPC校正窓内の純度推定

[剛性対純度](he4_purity.pdf) / [背景モデル比較](he4_background_models.pdf) / [全射影fit](he4_fits.pdf) / [数値表](purity.csv)

対象: Au+Au FXT: 3.85 GeV beam energy | 2021 P24iy | inclusive accepted events。

- `true/hDedxP` は核種PID・TOF前のtrack分布です。切断済み `hDedxP_4He` はfitには使いません。
- R=|p/q|を保存TH2の0.01 GeV/c刻みで射影し、元校正のnSigmaに軸を変換しています。
- 信号Gaussian＋正値linear/exponential背景＋範囲内の他核種Gaussianを比較し、最小Poisson AICを選びます。
- 積分窓は校正nSigma±2です。fit Gaussianの±2σに変更していません。
- 誤差棒はfit covarianceによる誤差です。CSVのmodel_differenceは背景関数間の差で、誤差棒に合算していません。
- 赤の白抜き点はfitが収束してもχ²/ndf等のモデル品質基準を満たさない参考推定です。信頼できる純度の確定値とは扱えません。
- 正負電荷・全採用中央度を合算しています。最終best-species選択後、新しいTOF併用後の純度ではありません。
- 校正窓は剛性bin中央で評価しています。bin内の連続track位置は復元できず、端での窓変化をJSON/CSVに記録しています。

fit品質: {"quality": 0, "exploratory": 175, "unavailable": 35, "linear": 40, "exponential": 135}。

`context.json` に入力fingerprint、mainconf、元summary、設定・ソースのSHA256、`purity.json` に全モデルfitの状態・幅・積分・残差を保存しています。入力ROOTは読み取り専用、変更なしを確認しました。
