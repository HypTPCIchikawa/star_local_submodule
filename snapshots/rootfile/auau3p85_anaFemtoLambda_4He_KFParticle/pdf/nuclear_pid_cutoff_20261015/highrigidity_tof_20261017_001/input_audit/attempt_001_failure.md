# 初回監査実行の停止記録

2026年10月17日。初回の`replay_audit.py`は、snapshot保存後、`merged_audit.tsv`の
区切り文字を文字列`\\t`と誤認したため、`Unexpected merged histogram count`で停止した。
元監査表は実タブ区切りで282ヒストグラムを記録している。元データの不一致ではない。

監査readerのみを実タブ／表示用escape双方に対応させた。初回作成snapshotは残し、
再実行では内容が完全一致する場合だけ再利用し、異なるファイルへの上書きは拒否する。
元ROOT、元監査表、過去runtime bundle、既存cacheやfitは変更していない。
