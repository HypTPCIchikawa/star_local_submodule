# 共通PID多項式とALL純度の再評価（2026-10-20）

## 結論

候補多項式の作成、旧新COLZ重ね描き、ALL低R＋TOF高Rの核種分離、Purity vs Rの出力まで実行した。
ただし、**低Rのfit品質を満たせず、新しい低R cutoffは確定できていない**。
未評価を「純度50%未満」へ読み替えず、cutoff・孤立した除外区間・高R最終上限は本番へ反映していない。
本番6入口はGitの元PIDを保持する。[診断後の保持監査](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/pid_git_restore_20261020_001/restore_audit_after_shared_study.json)も成功した。

これは多項式という関数族が原理的に不可能という証明ではない。
**今回の係数推定と固定応答では、低Rの純度からcutoffを認定できない**という結果である。

## 主な出力

| File | Content |
|---|---|
| [Purity vs R（基準候補、2ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/shared_purity_summary_002/purity_vs_r_quality_masked.pdf) | 全Rと低R拡大。良好なモデルだけ表示、未評価は灰色× |
| [全352fit、linear/log（176ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/shared_curve_linear_001/views_001/shared_curve_fit_diagnostics.pdf) | 共通多項式固定の4核種成分・滑らかな背景・合計 |
| [旧新COLZ、基準候補（8ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/polynomial_calibration_001/colz_002/pid_polynomial_colz_comparison.pdf) | ALL、d/4He共通TOF、t TOF、3He TOF、各低R/高R |
| [多項式圧縮／ALL転用の比較表](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/transfer_comparison_001/transfer_comparison.md) | 同じTOF・同じROIで局所fitと固定多項式を比較 |
| [基準候補の係数YAML](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/polynomial_calibration_002/baseline/coefficients.yaml) | 元と同じ次数、支持点・曲線の定義域・不自然な極値の監査 |
| [cutoff採否YAML](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/shared_purity_summary_002/cutoff_candidates.yaml) | 未認定の下限／上限はnull。productionEnabled=false |
| [純度・品質・統計の全点JSON](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/shared_purity_summary_002/purity_summary.json) | 値・条件付き誤差・有効モデル・棄却理由・統計閾値感度 |

### 同じ関数形の制約付きGLS比較

独立な重み付き回帰だけに結論を依存させないため、保存済みの中心・幅・非対称性の全共分散を使うGLSも試した。
無制約GLSでは非対称性が既定の0.25〜4を外れたので採用せず、
同じ10係数のまま正の中心・幅、幅≥0.05、非対称性範囲を課す制約付き解も検証した。
関数族やピーク応答は追加していない。

| File | Content |
|---|---|
| [制約付きGLS係数](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/polynomial_calibration_004/joint_gls_bounded/coefficients.yaml) | 制約・収束・元係数との差を保存。本番未採用 |
| [制約付き曲線のCOLZ（8ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/polynomial_calibration_004/colz_bounded_001/pid_polynomial_colz_comparison.pdf) | 同じraw分布へ重ね描き |
| [制約付き候補のPurity vs R（2ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/shared_bounded_purity_summary_001/purity_vs_r_quality_masked.pdf) | 改善とは認定しない。低Rは引き続き未評価 |
| [制約付き候補の全352fit（176ページ）](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/shared_bounded_linear_001/views_001/shared_curve_fit_diagnostics.pdf) | 局所linear/log、4核種成分を全て保存 |
| [制約付き候補のTOF比較表](../../../../rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_restart_20261020/transfer_comparison_bounded_001/transfer_comparison.md) | 制約付きでも全Rを十分再現したとは言えない |

## 入力・選択・関数の定義

- 3.85 GeV、2021 P24iyの既存部分統合ROOTに由来する、PID前のnativeヒストcacheを使用した。
  26,956/26,994 subjob、欠損38件の部分統合であり、23B全統計を完全網羅したという意味ではない。
- R=|p/q|、x=log10(R)。中心はxの3次、幅はxの2次のまま。
- 確定した応答はsplit-normal＋Float16保存精度の積分。
  中心はmode、幅はsqrt(σ_left σ_right)で、Gaussian標準偏差や全RMSではない。
  log(σ_right/σ_left)の2次式は応答の補助量で、PIDの中心・幅を別関数族へ変えたものではない。
- 選択は本番の核種ゲートに合わせ、|nσ|<2、かつd/t/3He/4Heの最小|nσ|を持つ核種。
  同値の場合は本番と同じd、t、3He、4Heの順。native R-bin内のGL3各点でゲートを計算する。
- 低RはALL（m²要求なし）、R≥2はTOFタグ。
  dと4Heは同じ **3.268773<m²<3.693147**、tは7.2564275<m²<8.2816925、3Heは1.78285<m²<2.06485。
- cacheはcharge-inclusive・centrality積分のtrack QAで、Λとのペアやk*で重み付けしたpurityではない。
  TAGはALLの部分集合なので、独立データとして加算していない。
- 元ヒストのcount/error/flow・native軸は変更していない。Float16は予測積分の側で扱った。

## 実行結果と採否

全380 TOF局所fitを実ピーク中心±1.5元σで再実行し、支持点はd=57、t=86、3He=86、4He=90を得た。
ゼロ信号、未収束、境界解、無効な共分散、ROI切れ、未確認の検索境界を支持点と分離した。
旧Gitの他核種曲線が近いことだけでは、正しいTOF支持点を棄却しない。
支持範囲の共通部分のみで352fitを実行し、共通範囲外28点を「低純度」にしていない。

各係数候補で、線形BG、一定BG、fit幅3σ→3.5σの3条件を比較した。
4Heでは3Heの頂上と左肩を含める。背景核種の中心・幅を局所スライスだけ独立freeにしていない。

| Coefficient candidate | Linear BG accepted | Constant BG accepted | Wider ROI accepted |
|---|---:|---:|---:|
| Independent weighted polynomial | 94 / 352 | 73 / 352 | 89 / 352 |
| Bounded joint GLS polynomial | 73 / 352 | 44 / 352 | 73 / 352 |

上表は今回の**暫定機械基準**による採否で、曲線や物理的粒子帰属の本番認定ではない。
局所χ²/NDF≤25、条件付き共分散有効、target有効count≥100、target収量相対誤差≤20%、
targetを除いたfitからのΔdeviance≥9、情報行列の条件数≤10^8等を用いた。
係数用支持点にはχ²閾値5/10/25の感度比較を保存した。
χ²の悪い点を採るために閾値を追加緩和してはいない。

Purity図は基準fitと、少なくとも他1つの良好モデルがある点のみを表示する。
モデル差の低い側は**良好なfitだけ**の最小purityで、不良fitを系統誤差へ混入させない。

| Species | Low R assessed, weighted | High R assessed, weighted | Low R assessed, bounded GLS | High R assessed, bounded GLS |
|---|---:|---:|---:|---:|
| d | 0 | 2 | 0 | 0 |
| t | 0 | 4 | 0 | 4 |
| 3He | 0 | 38 | 0 | 24 |
| 4He | 0 | 47 | 0 | 43 |

高Rの3He/4Heには評価できる点が残る。一方、低Rを50%で打ち切る判断に必要な良好fitは得られていない。
未評価の空白を跨ぐ単一cutoffは作らず、孤立vetoも追加していない。
高Rの評価可能範囲の端は、統計不足・fit品質・共通曲線の支持端を区別する必要があるため、最終上限と呼ばない。

## 不一致の切り分け

保存した局所TOF応答を新しい分離エンジンへそのまま渡すと、devianceが再現された。
したがって、少なくとも代表8点ではエンジン変更だけが不一致を作ったとは考えにくい。
同じTOF標本・同じROIでも、局所fitを1本の多項式に固定するとfit品質が悪化した。

例：R=0.875の3Heは局所χ²/NDF=1.94に対し、固定多項式のTOF fitは484.94。
中心の差は局所幅の約−0.34倍。dでは多項式幅が局所幅の約0.67倍になった。
単に「高統計なので小さな差が見える」だけではなく、幅や中心の圧縮誤差を確認すべき結果である。
ALLでは標本の違い・隣接峰・背景の形も加わり得るため、TOF圧縮だけを唯一の原因とはしない。

## 次に必要な作業

1. ALLの各核種ピークから支持点を取り直し、**背景核種は共通曲線固定**のまま対象核種を確認する。
2. その支持点を同じ次数の共通多項式へ戻し、全核種を同時に再評価する反復を試す。
3. それでも不一致が残る場合、PIDゲートの多項式とpurity推定用の応答校正を分けるか、
   同じ関数形を維持する制約との優先順位を相談する。勝手にスプラインや独立free背景peakは導入しない。

今回の結果から本番PIDやcutoffを更新する段階ではない。
閾値を緩めて「新cutoffが完成した」とは報告しない。

## 検証・保持

- Python回帰165件（資料リンク検査を含む）、ROOT5人工試験4,041項目成功。
- 各fitのraw native count/error/flow、予測成分の加算、effective-bin対応を再読込検証。
- 出力PDFはGhostscriptで全ページを読み戻した。代表ページを画像化して表示範囲も確認。
- 全出力を新規フォルダへ分離。初期のコンパイル失敗・未完了フォルダも履歴として保持した。
- 本番6入口のGit PID、元ROOT、過去の解析結果は保持。job投入、hadd、git commit/pushはしていない。
- 子タスクの承認待ちが長時間不可視に滞留したため、実行と保存を親へ集約。
  必要な主要コマンドの手動承認後、自動承認へ戻せる旨をユーザーへ通知した。
