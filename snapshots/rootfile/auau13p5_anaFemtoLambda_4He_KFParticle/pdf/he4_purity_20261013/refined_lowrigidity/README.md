# ⁴HeのTPC校正窓内の純度推定

[剛性対純度](he4_purity.pdf) / [背景モデル比較](he4_background_models.pdf) / [全射影fit・log](he4_fits.pdf) / [全射影fit・linear](he4_fits_linear.pdf) / [数値表](purity.csv)

対象: Au+Au FXT: 13.5 GeV beam energy | 2020 P24iy | inclusive accepted events。

- `true/hDedxP` は核種PID・TOF前のtrack分布です。切断済み `hDedxP_4He` はfitには使いません。
- R=|p/q|を保存TH2の0.01 GeV/c刻みで射影し、元校正のnSigmaに軸を変換しています。
- dE/dxの元bin幅は 0.1 keV/cm、fit bin幅は 0.5 keV/cmです。剛性bin幅・総count・元窓内countは変更していません。
- 信号Gaussian＋正値linear/exponential背景＋近傍核種Gaussianを比較し、最小Poisson AICを選びます。追補設定では平均が範囲外でもtailがfit範囲に重なる近傍成分を残します。
- 積分窓は校正nSigma±2です。fit Gaussianの±2σに変更していません。
- 誤差棒はfit covarianceによる誤差です。CSVのmodel_differenceは背景関数間の差で、誤差棒に合算していません。
- 赤の白抜き点はfitが収束してもχ²/ndf等のモデル品質基準を満たさない参考推定です。信頼できる純度の確定値とは扱えません。
- 正負電荷・全採用中央度を合算しています。最終best-species選択後、新しいTOF併用後の純度ではありません。
- 校正窓は剛性bin中央で評価しています。bin内の連続track位置は復元できず、端での窓変化をJSON/CSVに記録しています。

- `raw_native_*`（rebin時のみ）は元binの射影、`raw_*`・`projection_*`はfit binです。fit曲線はbin平均を用いて最適化しています。
- 低剛性の局所最大・肩、信号なしモデルとの比較、fit品質は別項目です。fit不採用はピーク不存在を意味しません。

fit品質: {"quality": 169, "exploratory": 30, "unavailable": 11, "linear": 71, "exponential": 128}。

`context.json` に入力fingerprint、mainconf、元summary、設定・ソースのSHA256、`purity.json` に全モデルfitの状態・幅・積分・残差を保存しています。入力ROOTは読み取り専用、変更なしを確認しました。


[校正候補・元窓／fit中心窓の比較表](he4_locus_fit_candidates.csv) / [study用PID曲線候補](he4_pid_curve_candidate_13p5.yaml)。数値的な成分fitからの候補であり、モデル依存性・品質判定を含めて確認してください。productionのPID設定は変更していません。
