# 保存済み⁴He fitの追加表示

[全剛性bin・linear scale](he4_fits_linear.pdf) / [低剛性dE/dx・log/linear比較](he4_low_rigidity_dedx.pdf)

- 元study: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013`。追加のfitは行っていません。
- 全binのlinear図は元の校正nSigma軸。低剛性比較は物理dE/dx軸（keV/cm）です。
- 黒点は保存射影、赤/青はlinear/exponential背景を使った全fit曲線、緑は採用モデルの⁴He Gaussian、紫はその背景です。
- 校正窓±2σを灰点線、採用Gaussian中心を緑破線で表示。fit曲線は座標だけを変換し、元保存binあたりの高さを保っています。
- `No accepted fit` はfitの品質・境界条件で不採用になったことを意味し、ピークが存在しないとの判定ではありません。
- TPC前PID、正負電荷・全採用中央度の合算。最終⁴He候補やTOF併用後の純度ではありません。
- `completed.json`のハッシュと入力fingerprintを照合し、元ROOTと旧PDFの変更なしを確認しました。
