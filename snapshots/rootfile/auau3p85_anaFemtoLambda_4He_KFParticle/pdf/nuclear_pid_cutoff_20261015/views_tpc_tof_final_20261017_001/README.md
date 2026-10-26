# 低R TPC-only＋高R TPC/TOF purityの統合図

`purity_tpc_tof_combined.pdf`は計6ページです。新しい概要・高R拡大・条件説明の3ページに、元の`views_final_001/purity_cutoff_overview.pdf`3ページをそのまま連結しました。元ファイルは変更していません。

緑線は本番採用R域、斜線は除外域、赤白抜き点は未評価です。⁴HeのR≥3.2は数値purityが高く見えても応答支持を通らないため除外します。帯はモデル感度であり信頼区間ではありません。

高Rのd/t/³Heはnative幅0.01のS/T加算、⁴Heは全域幅0.05の独立fitです。⁴HeのTOFはユーザー指定のd同一窓です。元の広い3<m²<4窓のpurityは未評価です。元ROOTは2021 P24iyの部分統合で、inclusive track QA由来の結果です。

`completed.json`に閉PDF全頁読戻し、3／6頁確認、入力不変と出力SHAを記録しています。`source/`に実行時コードと設定を保存。各ページのPNGも併記しました。このREADMEは実行完了後の説明で、completedの出力SHA対象外です。
