# 保存済み⁴He fitの追加表示

[全剛性bin・linear scale](he4_fits_linear.pdf) / [低剛性・log/linear比較](he4_low040_components_full_dedx.pdf)

- 元study: `/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer/rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040`。追加のfitは行っていません。
- 全binのlinear図は元の校正nSigma軸。低剛性比較は物理dE/dx軸（keV/cm）です。
- 低剛性PDFは各binを4面で表示：校正nSigma軸のlog/linearと、物理dE/dx軸のlog/linear。
- 上下とも校正nSigma -8.0〜4.0 の同じ物理範囲です。保存軸0〜120 keV/cmに達する場合は、両軸とも同じbinまでclipします（拡大表示ではありません）。
- 現在mainconfの差異を明示的に許可したruntime専用実行です。旧/現在SHAと設定全文を保存し、STAR SL24y・2021 P24iy・3.85 GeV FXTを照合しました。現在PID/cutは元fitへ再適用していません。
- 黒点は保存射影、赤/青はlinear/exponential背景を使った全fit曲線、緑は採用モデルの⁴He Gaussian、紫はその背景です。
- T・³HeのGaussian背景を個別表示し、D tailが含まれる場合も区別します。不採用fitしかないbinでは最低AICモデルの成分を診断表示し、`unaccepted`と明記します。
- 校正窓±2σを灰点線、採用Gaussian中心を緑破線で表示。fit曲線は座標だけを変換し、保存済みfit用binあたりの高さを保っています。bin幅は各図の縦軸・脚注と`render_audit.json`に記録しています。
- `No accepted fit` はfitの品質・境界条件で不採用になったことを意味し、ピークが存在しないとの判定ではありません。
- TPC前PID、正負電荷・全採用中央度の合算。最終⁴He候補やTOF併用後の純度ではありません。
- y軸は元fit用binあたりの計数です。座標変換によるJacobian倍率は掛けず、関数値の一致を数値検証しました。
- `completed.json`のハッシュと入力fingerprintを照合し、元ROOT・fitソースsnapshot・旧PDFの変更なしを確認しました。
