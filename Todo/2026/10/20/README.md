# 2026-10-20 核PIDをGitの元関数から再構築する

**最新：[共通多項式とALL純度の評価結果](summary_polynomial_all_purity_20261020.md)。**
候補係数、旧新COLZ、ALL低R＋TOF高Rの全352点×3条件、Purity vs Rを保存した。
同じ次数の制約付きGLSも比較したが、低Rの良好fitが足りず、新cutoffは未認定。
新しい応答ファミリーは増やしておらず、本番6入口はGitの元PIDを保持する。

**前段階の区切り： [現在の応答・fit範囲比較のまとめ](summary_current_responses_20261020.md)。**
応答をさらに増やさず、実行済みの380点surveyと代表8点の比較・品質検査を整理した。
新係数・純度・cutoffはまだ本番採用していない。

[STAR_Todo20261020.txt](STAR_Todo20261020.txt)に、2026-10-19終了時の新しい依頼を保存しました。
[10/19の日報](../../../../analysisnote/20261019/summary20261019.md)から再開できます。

今回の方針は、以前の補間校正・fitの採否・cutoffを前提にせず仕切り直すことです。
過去のファイルを削除するのではなく、過去の結論を新しい解析の既定値へ持ち込みません。

1. **最初に本番PIDをGitの元式へ戻す。** 元headerは不変なので、現在有効な³He・⁴Heの校正上書きを解除し、実行時の関数値も確認する。
2. m²選択したd・t・³He・⁴HeのdE/dxをRでスライスしてfitし、元と同じ多項式の係数だけを必要に応じて更新する。⁴Heのm²窓はdと同じ。
3. m²なしALLへ更新関数を適用し、低Rを改善する。信頼できる純度S/(S+B)が50%未満なら核種別cutoffを設ける。
4. R≥2では全核種にTOF m²を要求し、統計から高R上限を決める。低Rと合わせてPurity vs p/qを作る。

中心は`log10(R)`の3次式、幅は2次式を維持します。生COLZへの関数重ね描きと、各Rの1D fit図を保存します。
未評価・不良fitを低純度と誤認せず、途中だけの除外は前後のfitを示して相談します。

**2026-10-20再開時点：6本番入口をGitの元PIDへ復帰し、ROOT5で864行・2592項目を照合しました。**
元の式をraw ALL／TOF分布へ重ねた8ページPDFと、8スライスのGaussian試験fitを保存しました。
380スライスの診断fitと190ページPDF、fit範囲を狭めた8点比較も保存・検査しました。
ユーザー承認のGaussian＋左右指数tailも代表8点で比較し、[8ページPDF](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/response_comparison_003/response_comparison.pdf)へ保存しました。
さらにFloat16保存の量子化をROOT5で57,075項目検証し、同じ8点の予測積分へ取り入れました。
[保存精度の診断PDF（4ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/storage_default12_smoke_001/views_001/tag_peak_fit_diagnostics.pdf)も保存済みです。
細かい周期構造の一部は説明できますが高χ²は残るため、続きはピーク応答とfit範囲を検証するところからです。
新しい多項式係数・純度・cutoffは未認定です。
既存のcutoff／除外区間／TOF条件は今回の復帰では変更しておらず、再評価の結論とは区別します。
経過・出力リンクは[再開作業レポート](report_pid_restart_20261020.md)に記録します。

10/13のTodoは[../13/](../13/README.md)へ移動済みです。過去のROOT／PDF／job snapshotは変更していません。
