# Gaussian／log-normalの保存fit比較

既存のfit結果を同じnative入力・fit窓で比較しました。再fit、PID変更、本番採用は行っていません。

- 中心はGaussianの平均とlog-normalのmodeです。幅はGaussianの標準偏差とlog-normalの局所曲率幅で、同じ定義ではありません。
- 赤点は核種同定未検証、または数値失敗・zero yield・軸外・境界などの注意がある点です。未支持の点は線で接続していません。
- dと4Heのdevianceは同じTAG_D fitの値です。独立データとして二重計数しません。
- Δdevianceは非nestedな関数モデルの比較指標であり、有意度ではありません。

| Background | TAG sample | Native slices | Gaussian D | Log-normal D | Log/Gauss D | Log-normal lower D |
|---|---|---:|---:|---:|---:|---:|
| linear | TAG_D | 220 | 2.61287e+07 | 4.60636e+07 | 1.763 | 99 |
| linear | TAG_T | 220 | 2.92602e+06 | 6.21028e+06 | 2.122 | 58 |
| linear | TAG_HE3 | 220 | 3.4017e+06 | 5.80012e+06 | 1.705 | 39 |
| exponential | TAG_D | 220 | 2.3064e+07 | 4.33672e+07 | 1.88 | 96 |
| exponential | TAG_T | 220 | 2.96563e+06 | 6.24871e+06 | 2.107 | 59 |
| exponential | TAG_HE3 | 220 | 3.46273e+06 | 5.93755e+06 | 1.715 | 36 |

## 出力

- [中心・幅・fit-window devianceの比較](response_parameters.pdf)
- [保存スペクトル：全dE/dx軸とpeak ROI、linear/log](saved_spectra_overlays.pdf)
- [全native区間の数値・警告](comparison.tsv)
- [同じnative count/error/flow・軸・fit maskの照合](native_comparison_audit.json)
- comparison_display_inputs.rootには未変更のnative観測分布と表示用に加算した分布を別名で保存しています。

数値収束・モデル比較の改善だけでは核種同定、TOF→ALLの移植、intrinsic幅、purityの検証を意味しません。
