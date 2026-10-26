# KFParticle Λ–4He summary 用の後処理

既存の `tools/summary_4He{,_3p85}.tex` と描画マクロを参考に、KFParticle で作成済みの部分統合 ROOT を読み取り専用で扱います。元の解析マクロ・旧 PDF・ROOT を上書きしません。

## 一括実行

プロジェクト直下で実行します。`mainconf` は ROOT5 / STAR の実行環境と Maker の質量ゲートの入口です。
`summary_settings.yaml` は描画の rebin・規格化区間・質量フィット・LLフィット・centrality 集約の設定で、
production の選択cutとは別です。設定が読めなければ停止します。

現在の共有後処理設定は **3.85 GeV の再描画用**です。CF bin は 20 MeV/c、表示は 0–0.3 GeV/c、
規格化は **0.3–0.4 GeV/c**（20 MeV/c幅で5 bin）です。Method 5 の純度推定は別の 50 MeV/c 格子と
5 bin（内部で 250 MeV/c、端では短縮）の fit 窓を維持します。元の 5 MeV/c bin に所属する
純度を掛けてから 20 MeV/c へ合算するため、50/20 MeV の境界が揃わない部分も保持されます。
同じ出力binで共有する純度誤差は合算して伝播しますが、隣接fit間や観測数–fit間の共分散は未評価です。
LL fit の計算区間は描画軸と独立で、0–0.3 への表示変更だけで fit 範囲は切り詰めません。
LL fit は従来の0–0.4 GeV/cを維持し、新しい規格化区間と一部重なります。
従来のLL fitは規格化によるbin間相関を完全には扱わないため、fit誤差・χ²は参考値です。
規格化0.4–0.6版は3.85 GeVの `pdf/revisions/norm04_06_before_norm03_04_20261003/`、
質量半窓studyの旧版は `pdf/revisions/mass_halves_norm04_06_before_norm03_04_20261003/` に保存しています。

**13.5 GeV の現在の図・レポートは旧50 MeV/c版のままです。** 再描画を実行する指示は取り消され、
元の成果物をhash一致で復帰しました。元の設定・生成コードはその `pdf/revisions/50MeV_0p6_before_20MeV_20261003/source_snapshots/`
に保存されています。13.5 GeVを再現する場合はこの保存設定を明示し、新しい出力先へ生成してください。

```bash
python3 tools/kfparticle_summary/run_summary.py \
  --energy 13p5 \
  --mainconf config/mainconf/main_auau13p5_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau13p5_anaFemtoLambda_4He_KFParticle/auau13p5_anaFemtoLambda_4He_KFParticle_2020P24iy_partial_merge.root

python3 tools/kfparticle_summary/run_summary.py \
  --energy 3p85 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root
```

既定の保存先は入力ROOTと同じ階層の `pdf/`。同じ入力・設定の中断再開には `--resume` を付けます。
入力・設定を変えた再計算は、`--output-dir <新しい比較用ディレクトリ>/pdf` へ分けてください。
`pdflatex`・Ghostscript (`gs`)・Python3/PyYAML と既存の STAR Singularity 環境が必要です。
`script/singularity_root_batch.sh` 経由で実行するため、ホストの ROOT を直接起動しません。
SSH切断後も継続させたい場合は、上記コマンドを `tmux` 内で実行します。

`--workers 2` は ROOT マクロを独立した2プロセスで実行します。stateの書込みは直列化されます。
既存完成版を残して同じ出力先へ再描画する場合は、`--archive-previous <保存名>` を指定すると、
既存ファイルのhash検証後に `pdf/revisions/<保存名>/` へ移動してから開始します。
設定を変更する前に退避だけを行う場合は、同じ入力・mainconfを指定し、
`--archive-previous <保存名> --archive-only` を併用します。検証・移動が完了したところで終了し、
新しい描画は開始しません。その後に設定を編集し、通常の実行コマンドで新規生成します。
変更前に保存した設定・生成コードを検証に使う場合は `--previous-tool-snapshot <保存ディレクトリ>` を明示します。
旧ROOTは変更せず、`mass_halves/` や `mass_halves20MeV/` の独立した比較結果は移動しません。
移動ファイルのhash・保存先は `revision_manifest.json`、旧設定は `source_snapshots/` に記録されます。
保存済みの過去設定・生成コードに対する検証は `validate_outputs.py --provenance-snapshot <source_snapshots>` で行います。

処理順は、質量・ヒストの読取り検査 → 旧マクロのKF用コピー生成 → 13本の描画 → PDFページ検査 → summary組版です。
`macros/` に実際のコピー、`logs/` に実行ログとフィット監査TSV、`macro_manifest.json` に原本と生成物のSHA256、
`run_state.json` に入力ROOTのfingerprint・設定hash・各段階の状態を保存します。
`build_report.py` は図を読み直してsummaryだけを再作成でき、ROOTの再処理はしません。
`validate_outputs.py --output-dir <ROOTのディレクトリ>/pdf` で、完成状態、設定・マクロ・図・レポートのhash、
元ROOTのsize/inode/mtimeが変わっていないことを再確認できます。PDFページ数は組版時にも検査しています。

`prepare_macros.py --refresh-generated` は、旧manifestとhashが一致する、このツールが生成したファイルだけを
`macros/archive_<日時>/` に退避して更新します。旧マクロ原本や未知の編集を上書きする機能ではありません。
LLフィットが不成立・不安定な図では数値欄を警告表示に置き換えます。数値収束した質量fitもモデル妥当性の保証ではありません。

summaryは今回の2つのjobの来歴を固定しているため、別のjobに流用する際は `build_report.py` の sample 定義も照合が必要です。
3.85 GeVでは2021 P24iyの統合ROOTのみを使用し、2019のarchiveは入力にしません。

## `InspectKFSummary.C`

ROOT 5 / ACLiC 用です。入力 ROOT、出力先、質量窓を引数で渡します。運用では Maker YAML / 投入時スナップショットの値を明示してください。

```text
InspectKFSummary(input_root, output_dir, mean, sigma, nsigma, outer,
                 fit_min=1.08, fit_max=1.15)
```

出力は `inspect.json`、`Lambda_InvMass.pdf`、`Lambda_InvMass.png`。同名出力が存在する場合は上書きせず終了します。

- `inspect.json` は inputEvents などの TParameter、主要ヒストの有無・型・Entries・軸、centrality 等の一次元 bin 内容を記録します。
- 純度の主推定は、左右 sideband の合計を **実際に積算した質量 bin 幅の比**で補正した `P = 1 − α(N_left+N_right)/N_signal` です。旧解析の固定 S/N は使いません。
- Maker の信号窓は未ビン化した候補への厳密なカットです。一方、この事後計算は bin 中心が nominal `[lower,upper)` に入る **bin 全体**を積算します。signal/SB の重複はなく、実際の bin 境界・bin 番号・幅比を JSON に保存します。保存 bin 内の候補分布は復元できないため、両者は完全同一ではありません。
- nominal の ±3σ と左右 3–12σ では合計 sideband の幅比は 1/3 ですが、実 bin 境界では微調整されます。背景を sideband 間で平坦とみなす近似であり、そのモデル系統誤差は統計誤差に含めません。
- `gaus+pol2` fit を比較として保存します。fit status、covariance status、χ²/NDF、その確率、背景積分を記録します。`numericallyValid` は数値収束・有限な積分値の判定だけであり、モデル適合の良さを保証しません。隠れた goodness 閾値は設けず、`usable` / `primaryEstimator` は常に false（CF 補正には採用しない）です。fit で production の質量ゲートは変更しません。
- inclusive purity は **ペア・k*・centrality ごとの純度ではない proxy**です。CFsub に使う場合はこの近似を明記し、`usableForCFsub=false` の値で補正してはいけません。

## 旧レポートからの注意点

- `tools/summary_4He_3p85.tex` は S/N=5.217 と記載しますが、現行の 3p85 描画 CFsub マクロは 2.172、LL フィット CFsub マクロは 5.217 で不一致です。新しい ROOT に転用できる定数ではありません。
- 通常 SB-sub / CFsub の一部 group は 20–60%、Method 3/5 の対応 group は 20–80% です。方法比較では centrality を明示・統一する必要があります。
- Method 5 の滑動窓は隣接純度推定に相関を導入します。従来マクロの独立 bin の χ² や誤差のみで確定的な物理解釈をしないでください。

## Λ質量の低側・高側のCF比較（3.85 GeV）

今回のstudyは **3.85 GeV・2021 P24iyのみ**です。13.5 GeVの図・summaryは変更前のものを維持し、質量半窓のstudy対象にも含めません。

```bash
python3 tools/kfparticle_summary/compare_mass_halves.py \
  --energy 3p85 \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root
```

出力は同じROOTディレクトリの `pdf/mass_halves/` です。同じ入力・設定での再開は `--resume`。
日本語の `README.md`、図をまとめた `mass_halves_comparison.pdf`、方法別PDF/PNG、
計算値 `comparison.json`、元ヒストとfitの記録 `mass_half_extract.json`、描画canvasを保存します。
full（±3σ）、low（−3σ～0）、high（0～+3σ）をRaw、SBsub、CFsub、Method 3、Method 5で比較します。

- CFのbin幅は20 MeV/c、表示範囲は0～0.3 GeV/c。規格化は0.3～0.4 GeV/cの5 binです。
- Method 5の純度fitは50 MeV/c刻み、隣接5binの250 MeV/c窓を維持します。元の5 MeV/cカウントへ所属する純度binのPを適用してから20 MeV/cへ合算し、50と20の境界不一致を扱います。
- full/half、共通SB、規格化、同じcoarse純度の再利用による相関を誤差に反映します。ペア再利用、異なるfitの重複窓、カウントとfit背景の相関など未評価の項目はレポートに明記します。
- 失敗fit・規格化不能はNAです。有効なbinだけを使って比較区間を黙って変更することはありません。
- 検証は `python3 -m unittest discover -s tests -p 'test_mass_half*.py'`。ROOT/PicoDstの解析cutは変更しません。

## m_Tの4分割CF（3.85 GeV、Centrality 0–10%）

```bash
python3 tools/kfparticle_summary/run_mt_cf.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root \
  --purity-study rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mass_halves \
  --shared-method5-approved
```

出力は `pdf/mt_cent0_10/` の `cf_mt_cent0_10_constant_vs_method5.pdf` と日本語READMEです。
既存のsummary・質量窓別比較は変更しません。再開は同じ引数に `--resume` を追加します。
centralityと区間選択は `mt_settings.yaml`、CF・純度の刻みと規格化は `summary_settings.yaml` に記録します。

- cent9=7・8を加算。k*＝0～0.3 GeV/cの信号SEの累積25/50/75%に近いm_T bin境界を選び、規格化SE/MEも正であることを確認します。
- 今回はm_T = **2.40–2.60、2.60–2.70、2.70–2.80、2.80–3.90 GeV/c²**。すべて[下限,上限)です。
- constant purityは前回と同じ一定P、sidebandと規格化は各m_Tで個別計算します。
- **Method5はm_T積分のP_SE(k*)・P_ME(k*)を共通適用する、ユーザー承認済みの参考比較です。** m_T別の独立したfitではありません。この近似を明示するため `--shared-method5-approved` を必須にしています。
- 保存ROOTにはΛのtreeはありますが、原子核運動量・pair情報やk*×m_T×質量の共同分布はないため、m_T別の質量fitは復元できません。
- `mt_extract.json` に型・キー・軸・計数・射影検査、`mt_selection.json` に区間と統計、`mt_comparison.json` にCF・共有P由来の区間間共分散を保存します。未取得のfit–data相関などはPDF/READMEに記載します。
- 検証: `python3 -m unittest discover -s tests -p 'test_mt*.py'`。元ROOTやproductionの選択cutは変更しません。

### 固定幅m_Tの追加比較

`mt_fixed_02_settings.yaml` は m_T = **2.4–2.6、2.6–2.8、2.8–3.0、3.0–3.2 GeV/c²** を指定します。
分位点分割とは別に計算し、同じPDFの後半へ追記します。範囲外のtailは混ぜず監査に残します。
元の4ページとそのREADME・設定・生成コードは `revisions/before_fixed_02_20261003/` に保存済みです。

```bash
python3 tools/kfparticle_summary/run_mt_cf.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root \
  --purity-study rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mass_halves \
  --shared-method5-approved \
  --settings tools/kfparticle_summary/mt_fixed_02_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_10/fixed_02 \
  --reuse-extraction rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_10

python3 tools/kfparticle_summary/append_mt_report.py \
  --study-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_10 \
  --original-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_10/revisions/before_fixed_02_20261003 \
  --appendix-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_10/fixed_02
```

- `--reuse-extraction` は完成済み抽出の入力・centrality・checksum・抽出コードを検証して再利用し、ROOTの再読込みを省きます。
- `fixed_02/` の再開は上の計算コマンドに `--resume` を追加します。結合後のtop-levelを元の分位点コマンドで上書きすることは拒否します。
- 追加PDFは統計表・m_T分布・区間ごとの2方式比較・方式ごとの重ね描きの4ページです。元4ページは変更しません。
- m_T=3.0–3.2の最初のk* bin [0,0.02) はME sideband=0のためconstant purityのCFは欠測です。Method5は共通P(k*)の参考比較を維持します。
- 結合結果の来歴とhashは `combined_report.json` に保存します。各分割のJSON/TSV/図は別々に保持します。

### Centrality 0–20%の同じ比較

0–20%版は `mt_cent0_20_settings.yaml` と `mt_cent0_20_fixed_02_settings.yaml` を使用します。
正式な[cent9対応表](../../StRoot/StRefMultCorr/README.md#cent9-bin-index-convention-read-this-first)と
`StRefMultCorr.h` に従い **cent9=6～8** のSE/MEを合算してからCFを規格化します。
元の0–10%結果と解析cutは変更しません。

```bash
python3 tools/kfparticle_summary/run_mt_cf.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root \
  --purity-study rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mass_halves \
  --shared-method5-approved \
  --settings tools/kfparticle_summary/mt_cent0_20_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_20

python3 tools/kfparticle_summary/run_mt_cf.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root \
  --purity-study rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mass_halves \
  --shared-method5-approved \
  --settings tools/kfparticle_summary/mt_cent0_20_fixed_02_settings.yaml \
  --output-dir rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_20/fixed_02 \
  --reuse-extraction rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/mt_cent0_20
```

- 出力名は設定の `report_stem` で指定し、今回は `cf_mt_cent0_20_constant_vs_method5.pdf` です。
- Method5は保存済みの**0–20%質量分布のfit**を使います。0–10%のP(k*)を代用しません。ただしm_T積分P(k*)の4区間共通適用という近似は維持します。
- constant purityは比較条件を揃えるため、従来と同じinclusive一定Pを使用します。0–20%に独立fitした一定Pではありません。
- 分位点境界は0–20%のSE数で再選択します。固定幅境界は2.4～3.2 GeV/c²、0.2刻みです。
- 元の4ページを `mt_cent0_20/revisions/before_fixed_02_20261003/` に保存後、上の結合コマンドの対象を0–20%に変更し、`--report-stem cf_mt_cent0_20_constant_vs_method5` を追加して8ページにまとめます。
- bin幅20 MeV/c、表示0～0.3 GeV/c、規格化0.3～0.4 GeV/c、純度fitの50 MeV/c刻み・5bin窓はすべて維持します。

## Λ・娘p/π⁻・⁴HeのpT–Y_Lab分布

```bash
python3 tools/kfparticle_summary/run_acceptance.py \
  --mainconf config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml \
  --input-root rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/auau3p85_anaFemtoLambda_4He_KFParticle_2021P24iy_partial_merge.root \
  --source-summary rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf
```

- `acceptance_3p85_settings.yaml` にヒスト名・図のラベル・表示余白を保存します。
- `pdf/acceptance/acceptance_pt_ylab.pdf` は5ページ。4粒子の一覧と各粒子の拡大図です。PNG・個別PDF・canvas ROOT・監査JSON・日本語READMEも保存します。
- x=実験室系ラピディティY_Lab、y=pT [GeV/c]、z=元binの候補数（log）。擬ラピディティや重心系ラピディティに変換しません。
- Λと娘はKF高purity＋±3σのΛ候補ごと、⁴Heは最終選別trackごとに、組合せ前の1回を保存した分布です。⁴Heは同じイベントにΛがない場合も含みます。ペア重みはありません。
- 保存TH2には中央度軸がないため、採用イベント全体の分布を示します。前回の0–10%／0–20% CF選択は適用していません。
- これは候補の運動学的な占有範囲であり、生成粒子数で割ったアクセプタンス／効率ではありません。
- 共通Y_Lab範囲・粒子別pT範囲は全regular-bin候補を含むように自動設定します。flowは別計上し、rebin・規格化・cut変更はしません。
- 完成した元summaryの入力fingerprintとmainconf hashを照合します。既存出力は上書きせず、再作成する場合は新しい `--output-dir` を指定してください。
