#!/usr/bin/env python3
"""Draw the saved inclusive pT-YLab maps without changing the analysis."""
import argparse
from datetime import datetime, timezone, timedelta
import json
from pathlib import Path
import re
import shutil

import yaml

from compare_mass_halves import HERE, PROJECT, command, digest, fingerprint, save
from report_mass_halves import pdf_pages


def validate_settings(settings):
    if settings['energy'] != '3p85' or '2021 P24iy' not in settings['sample_label']:
        raise ValueError('This preset requires the established 2021 P24iy 3p85 sample.')
    padding = settings['padding_bins']
    if isinstance(padding, bool) or not isinstance(padding, int) or padding < 0:
        raise ValueError('padding_bins must be a nonnegative integer.')
    rows = settings['histograms']
    if len(rows) != 4 or len({r['id'] for r in rows}) != 4:
        raise ValueError('Exactly four distinct particle rows are required.')
    for row in rows:
        if not re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', row['id']):
            raise ValueError('Unsafe particle ID.')
        for key in ('id', 'histogram', 'title', 'selection'):
            if not isinstance(row[key], str) or not row[key] or any(c in row[key] for c in '\t\r\n'):
                raise ValueError('Manifest fields must be nonempty single-line text.')
    return rows


def write_readme(out, settings, audit):
    lines = ['# Λ・娘p/π⁻・⁴HeのpT–実験室系ラピディティ分布', '',
             '[まとめPDF（5ページ）](acceptance_pt_ylab.pdf)', '',
             '- 対象は3.85 GeV FXT、2021 P24iyの26,956/26,994 subjob部分統合。38件は欠損のままです。',
             '- 1ページ目は4粒子の一覧、2～5ページはΛ、p、π⁻、⁴Heの個別拡大図です。',
             '- 横軸は実験室系ラピディティY_Lab、縦軸はpT [GeV/c]、色は元bin内の候補数（対数表示）です。',
             '- Y_Lab = 0.5 ln[(E+pz)/(E−pz)]。ηではなく、重心系への移動・ビームラピディティの差引きはありません。',
             '- 保存済みヒストは採用イベント全体を合算しています。前回CFの0–10%・0–20%限定ではありません。中央度軸がなく、このTH2だけから中央度別の分布は復元できません。',
             '- 選択された候補の運動学的範囲を示す図であり、生成粒子数で割ったアクセプタンス／再構成効率ではありません。背景差引き・効率補正・面積規格化はしていません。',
             '- Λと娘粒子は高purity KF選別＋±3σ窓を通ったΛ候補ごとに1回、原子核との組合せ前に記録します。娘trackが異なるΛ候補で共有される場合は候補ごとに入ります。',
             '- Λの運動量はKF値、ラピディティ計算の質量は固定値1.115683 GeV/c²。娘はgMom()とp質量0.9382720813、π質量0.13957039 GeV/c²を使用します。娘の運動量をKF崩壊点へ輸送した図ではありません。',
             '- ⁴Heは原子核PID・最終選別後の1trackにつき1回。同じイベントにΛがなくても記録します。Heの運動量はZ=2補正済みの物理運動量で、p/qではありません。',
             '- ペア数・混合回数による重み付けはしていません。4粒子の母集団と色の範囲は異なるため、色だけで収量を比較しないでください。',
             '- 表示は4粒子共通のY_Lab範囲と粒子ごとのpT範囲を、非ゼロの保存binを全て含むように選びます。rebinはせず、空の外側だけを表示範囲から外します。flowは移し替えず、下表に記録します。', '',
             '| 粒子ID | 保存ヒスト | Entries | regular-bin合計 | under/overflow合計 | 表示Y_Lab | 表示pT [GeV/c] |',
             '|---|---|---:|---:|---:|---|---|']
    for h in audit['histograms']:
        lines.append('| %s | `%s` | %.0f | %.0f | %.0f | %.2f～%.2f | %.2f～%.2f |' %
                     (h['id'], h['histogram'], h['entries'], h['regular_count'], h['flow_count'],
                      *h['display_x'], *h['display_y']))
    lines += ['', '## 再現情報', '',
              '- `context.json`: 入力ROOTのfingerprint、mainconf、設定、元summaryの来歴。',
              '- `settings.yaml` / `histograms.tsv` / `macros/DrawAcceptance.C`: この描画に使用した設定・マクロ。',
              '- `acceptance_audit.json`: 軸・収量・flow・表示範囲の検査。',
              '- `completed.json`: PDFページ数、入力不変性、出力のSHA256。',
              '- `logs/plot.log`: mainconf、STAR/ROOT環境、実行コマンドと描画ログ。', '',
              '元ROOT・production設定・既存summary/CF PDFは変更していません。']
    (out/'README.md').write_text('\n'.join(lines)+'\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--input-root', type=Path, required=True)
    parser.add_argument('--source-summary', type=Path, required=True)
    parser.add_argument('--settings', type=Path, default=HERE/'acceptance_3p85_settings.yaml')
    parser.add_argument('--output-dir', type=Path)
    args = parser.parse_args()
    source = args.input_root.resolve(strict=True)
    mainconf = args.mainconf.resolve(strict=True)
    settings_path = args.settings.resolve(strict=True)
    settings = yaml.safe_load(settings_path.read_text())
    rows = validate_settings(settings)
    old_state = args.source_summary.resolve(strict=True)/'run_state.json'
    state = json.loads(old_state.read_text())
    prior = state['provenance']
    initial = fingerprint(source)
    if state['phase'] != 'complete' or prior['energy'] != settings['energy']:
        raise ValueError('Source summary must be the completed study of this energy.')
    if any(prior['input'].get(k) != v for k, v in initial.items()):
        raise ValueError('Input ROOT differs from the verified source summary.')
    if prior['configs'].get(str(mainconf)) != digest(mainconf):
        raise ValueError('Mainconf differs from the source summary.')
    out = (args.output_dir or source.parent/'pdf'/'acceptance').resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError('Refusing to overwrite existing acceptance output; use a new output directory.')
    out.mkdir(parents=True, exist_ok=True)
    (out/'macros').mkdir(exist_ok=True)
    (out/'logs').mkdir(exist_ok=True)
    macro = out/'macros/DrawAcceptance.C'
    shutil.copy2(str(HERE/macro.name), str(macro))
    shutil.copy2(str(settings_path), str(out/'settings.yaml'))
    manifest = out/'histograms.tsv'
    manifest.write_text('id\thistogram\ttitle\tselection\n'+''.join(
        '\t'.join(r[k] for k in ('id','histogram','title','selection'))+'\n' for r in rows))
    context = dict(input=initial, mainconf=str(mainconf), mainconf_sha256=digest(mainconf),
                   settings=settings, settings_sha256=digest(settings_path),
                   source_summary=str(old_state), source_summary_sha256=digest(old_state),
                   macro_sha256=digest(macro), centrality='inclusive stored accepted events',
                   created_jst=datetime.now(timezone(timedelta(hours=9))).isoformat())
    save(out/'context.json', context)
    call = str(macro)+'+('+','.join(json.dumps(str(s)) for s in
        [source,out,manifest,settings['energy_label'],settings['sample_label']])+','+str(settings['padding_bins'])+')'
    command(['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),
             'root4star','-l','-b','-q',call],out/'logs/plot.log')
    audit = json.loads((out/'acceptance_audit.json').read_text())
    if audit['inputRoot'] != str(source) or len(audit['histograms']) != len(rows):
        raise ValueError('Unexpected audit input or histogram count.')
    for expected,h in zip(rows,audit['histograms']):
        if h['histogram'] != expected['histogram'] or h['display_count'] != h['regular_count']:
            raise ValueError('Histogram mismatch or populated regular bins were hidden.')
    if fingerprint(source) != initial:
        raise ValueError('Input ROOT changed during rendering.')
    if pdf_pages(out/'acceptance_pt_ylab.pdf') != 5:
        raise ValueError('Expected five acceptance PDF pages.')
    write_readme(out,settings,audit)
    save(out/'completed.json',dict(phase='complete',pages=5,input_unchanged=True,
        sources={str(p):digest(p) for p in (Path(__file__).resolve(),HERE/'DrawAcceptance.C')},
        products={p.name:digest(p) for p in out.iterdir() if p.is_file() and p.name!='completed.json'}))
    print('COMPLETE',out/'acceptance_pt_ylab.pdf',flush=True)


if __name__=='__main__':
    main()
