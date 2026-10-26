# AuAu13p5 KFParticle ROOT出力一覧

整理日: 2026-09-08。今回の作業はファイルの移動・索引作成のみで、再解析・カット変更・ROOT内容の書換えは行っていない。

## 2026-09-14 KF-specific study

通常Imp5の最終距離・pointing cutを固定した再処理と、上流KF条件の検証結果は [KF-specific study出力索引](kf_specific_scan_20260914.md) を参照。実行中・未確定の結果は同索引で `pending` と区別する。以下の2026-09-08整理記録と既存出力は保持している。

## 質量比較に使うファイル

- **Imp5に共通カットを合わせたKF（10,000イベント）**: [Imp5/local_compare_10000_20260908.root](../Imp5/local_compare_10000_20260908.root)
- 標準KFの最初の比較（10,000イベント）: [local_compare_10000_20260907.root](../local_compare_10000_20260907.root)
- 比較相手の旧Helix（10,000イベント、移動していない）: [auau13p5_anaLambda/local_compare_10000_20260907.root](../../auau13p5_anaLambda/local_compare_10000_20260907.root)

Λのみの比較は、KFの `hKfLambdaMassSelected` と旧Helixの `hLambda_InvMass` を使う。標準KFの互換 `hLambda_InvMass` はanti-Λも含むため、そのまま旧Helixと重ねない。候補数は背景差引き後の信号収量ではない。

## 配置

以下は `rootfile/auau13p5_anaLambda_KFParticle/` からの相対パス。

```text
auau13p5_anaLambda_KFParticle/
├── local_compare_10000_20260907.root     標準KFの10,000イベント比較
├── Imp5/                               Imp5共通カットでの比較
├── validation_20260907/                 初期動作・vertex mode対応試験
├── validation_20260908_imp5/            Imp5実装後の標準KF非退行試験
├── validation_20260908_config_cleanup/  未使用設定整理後の非退行試験
└── mdfiles/                            この一覧・移動対応表・ハッシュ・初期失敗ログ
```

合計14 ROOTファイル（KF関係11、検証用の旧Helix対照3）。検証セット内の対照出力はセットと一緒に移動した。独立した `rootfile/auau13p5_anaLambda/` と `rootfile/auau13p5_anaLambda_Imp5_validation_20260908/` は変更していない。

## 条件の違い

成功したKF解析はROOT 5.34/38 / STAR SL24yのSingularity環境でローカル実行したもの。いずれもfull Finder/Topo経路で、親Λのmass constraintはかけていない。farm jobではない。

| 条件 | 標準KF（reference） | Imp5に合わせたKF |
|---|---|---|
| mainconf | `main_auau13p5_anaLambda_KFParticle.yaml` | `main_auau13p5_anaLambda_KFParticle_Imp5.yaml` |
| PID | `dedx_pull`、p/πは3、Kは2、TOF併用 | 保存済み `nSigmaProton/Pion` 各3、TOF・K仮説なし |
| Daughter quality | nHitsFit≥15、比≥0.52、nHitsDedx≥5、pT 0.15–10 GeV/c、lab η −2.4–0、dE/dx error 0.04–0.12 | nHitsFit≥15、比≥0.52。旧Makerで使われないpT/η/nHitsDedx/dE/dx error選別は追加しない |
| 元trackのPV距離・経路長 | Imp5用追加条件なし | proton gDCA≥0.7 cm、pion≥1.0 cm、元helixの絶対経路長≤100 cm |
| KF最終topology | Imp5用追加条件なし | daughter間距離≤0.5 cm、parent/PV距離≤0.5 cm、cos pointing≥0.998 |
| Event | mode対応後はFXT: Vz [198,202] cm、XY中心(−0.4,−2.0) cmから半径2 cm。RefMult/VPD条件も適用 | 旧Helix同様maxNTr＋centrality。Vz/Vr/RefMult/VPDの追加選別なし。XY中心はQA専用 |
| anti-Λ | 有効 | 無効 |
| 最終mass | [1.05,1.25] GeV/c² | 同左 |

共通で共分散の妥当性、Interface/Finderの条件を維持している（interfaceChiPrimaryCut=3、Finder daughter距離1.5 cm、L=1 cm、chiPrimary2D=3、chi2/NDF=10、L/dL=3、upstream Λ PV χ²/NDF<3）。標準KFのTOF係数はxwu2参照値で、AuAu13p5で較正済みという意味ではない。

古いファイルに現在のYAMLを遡って適用したものではない。正確な生成時の条件はROOT内の `KFParticleEffectiveConfiguration`、対応する実行ログ・設定snapshotを優先する。初期ファイルには現在のmetadataがないものもある。

## 各ROOTファイル

「読取」は指定件数ではなく実際の読取数。「再構成」はevent/centrality選別後に再構成処理へ進んだイベント数。「Λ / anti-Λ」は最終選別候補数。

### 比較用・最新の非退行検証

| ROOTファイル | 条件・目的 | 読取 / 再構成 | Λ / anti-Λ |
|---|---|---:|---:|
| [local_compare_10000_20260907.root](../local_compare_10000_20260907.root) | 標準KF。初回Helix比較 | 10000 / 9728 | 31248 / 4731 |
| [Imp5/imp5_smoke100_20260908.root](../Imp5/imp5_smoke100_20260908.root) | Imp5実装の短時間確認、完了・QA成功 | 100 / 93 | 30 / 0 |
| [Imp5/local_compare_10000_20260908.root](../Imp5/local_compare_10000_20260908.root) | Imp5共通カットのHelix比較、完了・QA成功 | 10000 / 9744 | 3135 / 0 |
| [validation_20260908_imp5/default_kf1000_after.root](../validation_20260908_imp5/default_kf1000_after.root) | **標準KF**。Imp5機能追加後も標準条件が変わらないことを確認 | 1000 / 957 | 3189 / 463 |
| [validation_20260908_config_cleanup/standard_kf1000_after_cleanup.root](../validation_20260908_config_cleanup/standard_kf1000_after_cleanup.root) | 標準KF。未使用設定整理後 | 1000 / 957 | 3189 / 463 |
| [validation_20260908_config_cleanup/imp5_10000_after_cleanup.root](../validation_20260908_config_cleanup/imp5_10000_after_cleanup.root) | Imp5 KF。未使用設定整理後 | 10000 / 9744 | 3135 / 0 |

`validation_20260908_imp5` というフォルダ名は「Imp5実装時の検証」を意味し、中の `default_kf1000_after.root` はImp5カットの出力ではない。

未使用設定の整理前後で標準KF・Imp5の全68ヒストグラム（軸・entries・全bin内容/誤差・flow bins）が一致した。ROOT metadataなどが異なるため、ファイル全体が同じハッシュという意味ではない。10,000イベントの比較図は変更していない。

### 初期KF動作・vertex mode対応試験（2026-09-07）

| ROOTファイル | 条件・状態 | 読取 / 再構成 | Λ / anti-Λ |
|---|---|---:|---:|
| [validation_20260907/smoke_1000.root](../validation_20260907/smoke_1000.root) | **失敗記録**。1000指定、ROOT5のNetx plugin不在でremote入力を開けず、status=2 | 0 / 0 | 0 / 0 |
| [validation_20260907/local_smoke_1000.root](../validation_20260907/local_smoke_1000.root) | **途中失敗記録**。local化後、CINT/StdEdxModel初期化でsegfault。shell exit=0でも正常完了ではない | 17 / 0 | 0 / 0 |
| [validation_20260907/local_smoke_1000_context.root](../validation_20260907/local_smoke_1000_context.root) | CINT問題修正後。標準KFだが**原点中心**Vr≤2 cmの旧vertex条件。982イベントがVrで除外 | 1000 / 12 | 1 / 0 |
| [validation_20260907/local_mode_1000.root](../validation_20260907/local_mode_1000.root) | 標準KF、FXT中心(−0.4,−2.0) cmへmode対応後。完了・QA成功 | 1000 / 957 | 3189 / 463 |
| [validation_20260907/local_mode_1479.root](../validation_20260907/local_mode_1479.root) | 直上と同条件、先頭PicoDst全イベント。完了・QA成功 | 1479 / 1422 | 4701 / 716 |

最初の2ファイルは物理比較には使わない。`local_smoke_1000_context.root` も現在のFXT頂点条件とは異なる過去の診断用出力で、mode metadata導入前のため現在の出力QAでは拒否される。失敗の根拠ログは [smoke_1000.log](logs/smoke_1000.log) / [local_smoke_1000.log](logs/local_smoke_1000.log) に保存した。

### 検証セットに含まれる旧Helix対照（KF出力ではない）

| ROOTファイル | 目的 | 読取 / centrality通過 |
|---|---|---:|
| [validation_20260907/helix_control_20.root](../validation_20260907/helix_control_20.root) | full-KF導入時の旧解析対照 | 20 / 17 |
| [validation_20260907/helix_mode_control_20.root](../validation_20260907/helix_mode_control_20.root) | mode対応後の旧解析対照。直上と全46ヒストグラム一致 | 20 / 17 |
| [validation_20260908_config_cleanup/helix20_after_cleanup.root](../validation_20260908_config_cleanup/helix20_after_cleanup.root) | 未使用設定整理後。元の旧Helix検証出力と全46ヒストグラム一致 | 20 / 17 |

旧 `StLambdaMaker` と旧mainconfで生成された対照出力であり、KFの候補数に合算しない。

## 入力と記録の参照先

元入力は [config/picoDstList/auau13p5GeV.list](../../../config/picoDstList/auau13p5GeV.list)。初期試験と20/100/1000/1479イベント試験は先頭ファイル `st_physics_adc_21033026_raw_6500002.picoDst.root`（全1479イベント）の先頭から。10,000イベント試験は同リスト先頭7ファイル（最初の6ファイル9293件＋7番目707件）を順に読んだ。

ローカル入力は `/tmp/star-kf-validation-20260907.Ho7eLz/auau13p5_first.picoDst.root` および `/tmp/star-kf-compare-input-20260907.c7jFYa/local.list`。/tmpは永続保存ではない。元URI、入力順、イベントID監査と設定snapshotは以下の記録・archiveで確認できる。

- [初期実装・mode対応の記録](../../../mdfiles/kfparticle/implementation/kfparticle_full_implementation_20260907.md)、[2026-09-07作業ログ](../../../analysisnote/20260907/summary20260907.md)
- [標準KF 10,000イベント比較](../../../mdfiles/kfparticle/comparisons/lambda_helix_kf_comparison_10000_20260907.md)、[設定・入力監査・ログarchive](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle/comparison_provenance_20260907.tar.gz)
- [Imp5 10,000イベント比較](../../../mdfiles/kfparticle/comparisons/lambda_helix_kf_imp5_comparison_10000_20260908.md)、[設定・入力監査・ログarchive](../../../share/figure/auau13p5_Lambda_Helix_vs_KFParticle_Imp5/comparison_provenance_20260908.tar.gz)
- [未使用設定の整理と検証](../../../mdfiles/configuration/unused_lambda_cut_cleanup_20260908.md)、[整理後設定・ソース・検証ログarchive](../../../share/figure/auau13p5_config_cleanup_20260908/cleanup_validation_20260908.tar.gz)

## 旧パスからの対応と移動検証

| 整理前（rootfile/直下） | 整理後（このKFParticleディレクトリ内） |
|---|---|
| `auau13p5_anaLambda_KFParticle_Imp5/` | `Imp5/` |
| `auau13p5_anaLambda_KFParticle_validation_20260907/` | `validation_20260907/` |
| `auau13p5_anaLambda_KFParticle_Imp5_validation_20260908/` | `validation_20260908_imp5/` |
| `auau13p5_config_cleanup_validation_20260908/` | `validation_20260908_config_cleanup/` |

13ファイルをフォルダごと移動し、元から指定先にあった標準KF比較1ファイルも含め全14ファイルのSHA256不変を確認した。削除・上書きはしていない。旧位置にリンクや空フォルダは残していない。

- [全ファイルの旧→新パスとSHA256](rootfile_move_map_20260908.tsv)
- [移動後のSHA256一覧](rootfiles_SHA256SUMS_20260908)

過去の実行ログ・再現用archive・ROOT内部metadataは履歴として変更していない。そこに残る旧パス（古いMarkdownのリンク・実行例を含む）は上表で読み替える。今回は保管済み出力の整理のみで、解析コード・mainconf・analysis_info・ジョブ生成設定は変更していない。以後ローカル再試験を行う際も、このKFParticleディレクトリ内の未使用出力名をrunnerに指定する。

プロジェクトrootで内容不変を再確認するコマンド:

```bash
sha256sum -c rootfile/auau13p5_anaLambda_KFParticle/mdfiles/rootfiles_SHA256SUMS_20260908
```
