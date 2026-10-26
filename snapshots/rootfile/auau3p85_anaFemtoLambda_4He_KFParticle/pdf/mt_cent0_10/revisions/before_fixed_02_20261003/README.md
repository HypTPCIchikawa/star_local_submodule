# Λ–⁴He CFのm_T 4分割（3.85 GeV、Centrality 0–10%）

[まとめPDF](cf_mt_cent0_10_constant_vs_method5.pdf)

- Centrality: cent9=7・8をSE/MEそれぞれ加算。signalは±3σ。
- CF: 20 MeV/c幅、表示0–0.3 GeV/c、規格化0.30–0.40 GeV/c。
- m_T: sqrt((|pTΛ+pTA|/2)² + ((mΛ+mA)/2)²)。
- 境界はk*＝0–0.3 GeV/cの信号SE数の累積25/50/75%に近い保存bin境界を選択。CF形状は最適化していません。

| Slice | m_T [GeV/c²] | SE (0–0.3) | ME (0–0.3) | SE norm | ME norm | SE leftSB (0–0.3) | SE rightSB (0–0.3) |
|---|---|---:|---:|---:|---:|---:|---:|
| 1 | 2.40–2.60 | 12175 | 436803 | 20366 | 759268 | 1202 | 2973 |
| 2 | 2.60–2.70 | 17473 | 625983 | 23059 | 843567 | 672 | 1694 |
| 3 | 2.70–2.80 | 14273 | 515692 | 16562 | 611419 | 345 | 929 |
| 4 | 2.80–3.90 | 14993 | 541619 | 16288 | 590353 | 191 | 479 |

各区間は[下限,上限)で、同じ境界を両方法・SE/ME・sidebandに使用します。

## 注意点

- Method5は0–10%・m_T積分のP_SE(k*)とP_ME(k*)を4区間で共通使用する参考比較です。m_T別に純度をfitした結果ではありません。
- 保存済みROOTにk*×m_T×Λ質量の同時分布がないため、純度のm_T依存は検証できません。
- constant purityは前回と同じinclusive質量分布から求めた一定値P=0.96721065を全m_T・k*に使用します。各m_Tのsideband CFを個別に規格化して補正します。
- 両方式とも各m_T区間のSE・MEから規格化します。完成したCF同士を加算して中央度を統合しません。
- Method5は純度50 MeV/c刻み・5bin（内部250 MeV/c）窓を維持。元5 MeV/cのカウントに共通Pを掛けて20 MeV/cへ合算します。
- 同一純度Pの再利用と規格化に由来する相関は計算に含めます。異なるfit間の共分散、fit背景とカウントの相関、m_Tカウントとm_T積分純度の相関、イベント・混合ペア再利用の相関は未評価です。
- constant purityのPは固定値として誤差計算します。誤差は近似で、最小k*でSEが0または少数の点はPoisson区間ではありません。
- Method5純度は質量bin全体、m_T分布は候補への厳密な±3σゲートから作られます。この境界差も残るため、m_T別独立fitの代替や精密な物理結論にはしません。
- 対象は2021 P24iyの26,956 / 26,994 subjob部分統合。元ROOT、選択cut、mixingは変更していません。
- 表示内の1 binあたり最少signal SEペア数はm_T区間順に1 / 2 / 5 / 3です。全区間の総数が多くても、最小k*付近は低統計です。
- flow監査のwarningは比較・規格化の範囲外のカウントに関するものです。今回のm_T underflow/overflowは0で、k*≥1 GeV/cのoverflowはCFに使用していません。

## ファイル

- `mt_distribution.pdf`: 選択した境界とSE/MEのk*–m_T分布。
- `mt_cf_comparison.pdf`: m_T区間ごとの2方式の比較。
- `mt_cf_overlays.pdf`: 方式ごとの4区間の重ね描き。
- `mt_selection.json`: 区間選択、規格化・sideband統計、ゼロbin・flow検査。
- `mt_comparison.json`: CF・近似誤差・共通純度・同一k*のm_T区間間共分散。
- `mt_extract.json`: 入力ROOTのキー一覧と抽出ヒスト、射影検査。
