# Task 2 実装計画：高剛性⁴HeのTOF選別

## 現状と依頼の対応

- 対象は`anaFemtoLambda_4He.C`から呼ばれる`StFemtoMaker`のΛ経路、具体的には`FemtoLambdaLegacy::CollectNuclei()`。
- 現行`nuclearid_femtoLambda_legacy.yaml`は`maxPOverQ: 2.5`、`maxNSigmaNuclear: 2.0`、`m2_selection: false`である。Todoの上限2.0という記述とは一致しない。
- 既存のdE/dx較正は`NuclearIdDeDxVsMom.h`の関数を使う。⁴Heの有効な剛性範囲は0.4～5.0 GeV/cであり、この範囲外へ較正を外挿しない。
- TOFの既存計算は剛性R=|p/q|から`R²(1/β²−1)`を求める。これは`(m/q)²`であり、⁴Heの中心は`(3.72742/2)² ≈ 3.4734 (GeV/c²)²`。物理的なm²の約13.9と混同しない。

## 修正方針

1. 核PID YAMLに任意の`he4HighRigidityTof`ブロックを追加する。ブロックがない既存設定は従来の選別を維持する。有効化時は境界・m²窓を必須キーとして読み、誤記・不正値でInitを失敗させる。
2. 最初に従来通りd/t/³He/⁴HeのTPC pullからbest speciesを決める。TOFで他種を除外したために⁴Heが繰り上がることを防ぐ。
3. best speciesが⁴Heの場合だけ、R < 2.0 GeV/cは従来のTPC選別、R ≥ 2.0はTPCに加えて有効なTOFと設定された`(m/q)²`窓を要求する。境界R=2.0は高剛性側に含める。
4. この新しい⁴He分岐は従来の上限2.5を外す。他核種の選別・上限は維持する。β、剛性、m²の非有限値を有効な測定として扱わない。
5. 初期TOF窓は3.0～4.0とする。³Heの中心約1.97を除外するための初期設定であり、⁴HeのTOF分解能をfitして最適化した窓とは主張しない。後からYAMLだけで変更できる。
6. 新機能専用の核PID YAMLと13.5/3.85 GeV用mainconfを作る。従来mainconfを指定すれば従来選別に戻せる。
7. 高剛性TPC候補のTOF分布、最終⁴HeのdE/dx/TOF分布と、TOF選別段階ごとの収量を追加する。初期段階の保存済みROOTからはTPCとTOFの同時選別純度を復元できないため、今後の検証用に保存する。

## 検証

- 境界の上下、TOFなし、m²窓内外、非有限値、best-species競合、他核種の上限保持、旧設定の回帰を人工trackテストで確認する。
- `script/singularity_make.sh <新mainconf> --no-clean`でSTAR SL24y/ROOT5と一致する環境で全ライブラリを再ビルドする。
- 既存Femto Λテストを実行し、対応した入力があれば2021 P24iyの短いローカルイベントループでInit・出力を確認する。
- 本番ジョブの投入は今回のTodoに含まれていない。検証結果と実行コマンドを日本語の報告にまとめる。
