# Lambda Cut 条件の比較サマリー (gaus + pol2 fit)

不変質量分布を `gaus + pol2` でFitし、S/N比および規格化収量（Yield / Event）を抽出しました。

Fit結果をまとめたPDFは以下のパスに出力されています：
[lambda_cut_comparison.pdf](file:///star/u/yichikawa/pwg/git_STAR_3/star-analyzer/rootfile/auau13p5_anaLambda/lambda_cut_comparison.pdf)

参考Notion: [Lambda SN investigation](https://app.notion.com/p/2026-06-30-38f2284fda2880c5bfc4c52b2287161f?source=copy_link)

---

## 📊 Fit結果とCut条件の一覧表

| 変更フェーズ | S/N比 | 規格化収量<br>(×10⁻³) | σ<br>(MeV) | 解析Event数 | 主なCut条件の特徴 |
| :--- | :---: | :---: | :---: | :---: | :--- |
| **Default (main)** | 0.540 | 296.358 | 1.77 | 689,838 | 元の緩い設定。一次粒子由来のBGが大量に混入。 |
| **Improvement 1**<br>（厳格版/ichikawa_tmp） | 4.481 | 144.498 | 1.63 | 716,241 | DCA・nHits を全て厳格化。S/Nは劇的に良いが収量が半減。 |
| **Improvement 2**<br>（提案A: DCA緩和版） | 0.847 | 281.870 | 1.72 | 117,324 | DCAを緩めすぎたため、一次粒子のBGが戻った。 |
| **Improvement 3**<br>（提案D: 収量保護版） | 1.718 | 221.605 | 1.70 | 134,209 | 一次粒子を弾きつつ、Helixの分解能に配慮した最適値。 |
| **Improvement 4**<br>（Imp3 + Pointing 0.998） | 2.329 | 214.433 | 1.70 | 103,955 | Imp 3 に `minCosPointing: 0.998` を追加。収量ロス3%でS/N+35%。 |
| **⭐ Improvement 5**<br>（Imp1 + Pointing 0.998） | **5.217** | **142.829** | **1.63** | 716,241 | **Imp 1 の高S/Nに `minCosPointing: 0.998` を追加。収量ほぼ維持でS/N更に向上。** |

<br>

### 📝 各フェーズの具体的なCut条件パラメータ

| パラメータ | Default | Imp 1 | Imp 2 | Imp 3 | Imp 4 | ⭐ Imp 5 |
| :--- | :---: | :---: | :---: | :---: | :---: | :---: |
| **`nSigmaProton`** | 2.0 | 3.0 | 3.0 | 3.0 | 3.0 | 3.0 |
| **`nSigmaPion`** | 2.0 | 3.0 | 3.0 | 3.0 | 3.0 | 3.0 |
| **`minDCAProton`** | 0.5 | 0.7 | 0.5 | 0.7 | 0.7 | 0.7 |
| **`minDCAPion`**   | 0.8 | 1.0 | 0.8 | 1.0 | 1.0 | 1.0 |
| **`maxDaughterDCA`** | 1.0 | **0.5** | 0.8 | 0.8 | 0.8 | **0.5** |
| **`maxDCAV0`**     | 1.0 | **0.5** | 0.8 | 0.8 | 0.8 | **0.5** |
| **`minNHitsFit`**  | — | **15** | 12 | 12 | 12 | **15** |
| **`minNHitsRatio`**| — | 0.52 | 0.52 | 0.52 | 0.52 | 0.52 |
| **`minCosPointing`** | 0.995 | 0.995 | 0.995 | 0.995 | **0.998** | **0.998** |
| `maxPathLength` | 100.0 | 100.0 | 100.0 | 100.0 | 100.0 | 100.0 |

---

## 💡 考察

### Improvement 5 の評価 ⭐
- **Improvement 1**（厳格カット全部入り）に対して `minCosPointing: 0.995 → 0.998` を追加しただけ。
- 収量はほぼ変化なし（144.5 → 142.8、**約1%の減少のみ**）にもかかわらず、S/N比は **4.48 → 5.22（+16%向上）** となりました。
- これは「本物のΛは運動量ベクトルがほぼ完全に一次頂点を向く」という物理的性質を利用し、cosθ < 0.998 の偽ペアだけを非常に効率良く排除できたためです。
- デフォルトと比較すると、**収量をデフォルトの 48% に保ちながら S/N 比を約9.7倍に向上** させています。

### 今後の展望
- **Improvement 5** が現状で最もバランスの取れた最終候補です。
- これ以上の S/N 改善を目指す場合は、KFParticle（高精度頂点再構成）や機械学習（BDT/MVA）の導入が必要になるフェーズと考えられます。
