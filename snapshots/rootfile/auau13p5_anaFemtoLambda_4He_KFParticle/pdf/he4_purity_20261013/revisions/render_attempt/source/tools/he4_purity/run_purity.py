#!/usr/bin/env python3
"""Fit existing pre-PID dE/dx TH2s using an explicit mainconf and YAML preset."""
import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(path):
    st = path.stat()
    return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns,
                device=st.st_dev, inode=st.st_ino)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False, allow_nan=False) + '\n')


def page_count(path):
    code = '(%s) (r) file runpdfbegin pdfpagecount = quit' % path
    result = subprocess.run(['gs', '-q', '-dNOSAFER', '-dNODISPLAY', '-c', code],
                            capture_output=True, text=True, check=True)
    return int(result.stdout.strip())


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--settings', type=Path, default=HERE/'purity_settings.yaml')
    p.add_argument('--sample', choices=['3p85', '13p5'], required=True)
    p.add_argument('--mainconf', type=Path, required=True)
    p.add_argument('--output-dir', type=Path)
    p.add_argument('--archive-attempt', type=Path, action='append', default=[],
                   help='Preserve a prior attempt from this output parent under revisions/.')
    args = p.parse_args()
    settings_path = args.settings.resolve(strict=True)
    cfg = yaml.safe_load(settings_path.read_text())
    sample = cfg['samples'][args.sample]
    mainconf = args.mainconf.resolve(strict=True)
    if mainconf != (PROJECT/sample['mainconf']).resolve(strict=True):
        raise ValueError('Explicit mainconf does not match the selected preset.')
    if cfg['models'] != ['linear', 'exponential'] or cfg['fitPanelsPerPage'] != 4:
        raise ValueError('This implementation requires linear/exponential models and four panels/page.')
    if not (0.4 <= cfg['rigidityMin'] < cfg['rigidityMax'] <= 5 and cfg['rigidityStep'] > 0
            and cfg['gateNSigma'] > 0 and cfg['signalSigmaMin'] > 0
            and cfg['signalSigmaMax'] > cfg['signalSigmaMin']):
        raise ValueError('Invalid calibration range, step or width constraints.')
    if args.sample == '3p85' and '2021P24iy' not in sample['input']:
        raise ValueError('3.85 GeV must use 2021 P24iy only.')
    source = (PROJECT/sample['input']).resolve(strict=True)
    old_state_path = (PROJECT/sample['sourceSummary']).resolve(strict=True)
    state = json.loads(old_state_path.read_text())
    initial = fingerprint(source)
    if state['phase'] != 'complete' or state['provenance']['input'] != initial:
        raise ValueError('Input differs from the verified completed source summary.')
    output = (args.output_dir or PROJECT/sample['output']).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError('Output already exists; use a new directory to preserve existing products.')
    output.mkdir(parents=True, exist_ok=True)
    (output/'logs').mkdir()
    macro = output/'source/tools/he4_purity/FitHe4Purity.C'
    calibration = output/'source/StMaker/common/NuclearIdDeDxVsMom.h'
    for src, dst in [(HERE/'FitHe4Purity.C', macro),
                     (PROJECT/'StMaker/common/NuclearIdDeDxVsMom.h', calibration),
                     (settings_path, output/'settings.yaml'),
                     (Path(__file__).resolve(), output/'source/tools/he4_purity/run_purity.py')]:
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(src), str(dst))
    numeric = {key: value for key, value in cfg.items()
               if key not in ('models', 'samples')}
    manifest = output/'settings.tsv'
    manifest.write_text(''.join('%s\t%s\n' % (k, int(v) if isinstance(v, bool) else v)
                                for k, v in numeric.items()))
    save(output/'context.json', dict(input=initial, settings=cfg, sample=args.sample,
        mainconf=str(mainconf), mainconf_sha256=digest(mainconf), source_summary=str(old_state_path),
        source_summary_sha256=digest(old_state_path), sources={str(f):digest(f) for f in (macro, calibration)},
        created_jst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        note='TPC-only charge-inclusive, inclusive centrality; existing ROOT has no joint TOF/best-species axes.'))
    call = str(macro)+'+('+','.join(json.dumps(str(x)) for x in (source, output, manifest, sample['label']))+')'
    command = ['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
               'root4star', '-l', '-b', '-q', call]
    log = output/'logs/fit.log'
    print('RUN', args.sample, output, flush=True)
    with log.open('w') as handle:
        handle.write('COMMAND: '+json.dumps(command)+'\n'); handle.flush()
        result = subprocess.run(command, cwd=str(PROJECT), stdout=handle,
                                stderr=subprocess.STDOUT, timeout=3600)
    if result.returncode:
        raise RuntimeError('ROOT fit failed; see '+str(log))
    content = log.read_text(errors='replace')
    if 'HE4_PURITY_COMPLETE' not in content or any(t in content for t in
       ('Error in <ACLiC>', 'segmentation violation', 'fatal error:', 'Error: Function')):
        raise RuntimeError('ROOT fit did not complete cleanly; see '+str(log))
    audit = json.loads((output/'purity.json').read_text(), parse_constant=lambda x: (_ for _ in ()).throw(ValueError(x)))
    expected_bins = round((cfg['rigidityMax']-cfg['rigidityMin'])/cfg['rigidityStep'])
    if audit['inputRoot'] != str(source) or len(audit['slices']) != expected_bins:
        raise ValueError('Unexpected input or rigidity-bin count in fit audit.')
    counts = dict(quality=0, exploratory=0, unavailable=0, linear=0, exponential=0)
    for item in audit['slices']:
        if item['chosen'] < 0:
            counts['unavailable'] += 1
            continue
        fit = item['fits'][item['chosen']]
        if not fit['valid'] or not 0 < fit['purity'] <= 1:
            raise ValueError('Chosen fit fails physical/fit validity checks.')
        counts['quality' if fit['quality'] else 'exploratory'] += 1
        counts[fit['model']] += 1
    if fingerprint(source) != initial:
        raise ValueError('Original input changed during purity estimation.')
    pages = dict(purity=page_count(output/'he4_purity.pdf'), fits=page_count(output/'he4_fits.pdf'),
                 background_models=page_count(output/'he4_background_models.pdf'))
    if pages != dict(purity=1, fits=(expected_bins+3)//4, background_models=2):
        raise ValueError('Unexpected PDF page count.')
    (output/'README.md').write_text(
        '# ⁴HeのTPC校正窓内の純度推定\n\n'
        '[剛性対純度](he4_purity.pdf) / [背景モデル比較](he4_background_models.pdf) / '
        '[全射影fit](he4_fits.pdf) / [数値表](purity.csv)\n\n'
        '対象: '+sample['label']+'。\n\n'
        '- `true/hDedxP` は核種PID・TOF前のtrack分布です。切断済み `hDedxP_4He` はfitには使いません。\n'
        '- R=|p/q|を保存TH2の0.01 GeV/c刻みで射影し、元校正のnSigmaに軸を変換しています。\n'
        '- 信号Gaussian＋正値linear/exponential背景＋範囲内の他核種Gaussianを比較し、最小Poisson AICを選びます。\n'
        '- 積分窓は校正nSigma±2です。fit Gaussianの±2σに変更していません。\n'
        '- 誤差棒はfit covarianceによる誤差です。CSVのmodel_differenceは背景関数間の差で、誤差棒に合算していません。\n'
        '- 赤の白抜き点はfitが収束してもχ²/ndf等のモデル品質基準を満たさない参考推定です。信頼できる純度の確定値とは扱えません。\n'
        '- 正負電荷・全採用中央度を合算しています。最終best-species選択後、新しいTOF併用後の純度ではありません。\n'
        '- 校正窓は剛性bin中央で評価しています。bin内の連続track位置は復元できず、端での窓変化をJSON/CSVに記録しています。\n\n'
        'fit品質: '+json.dumps(counts, ensure_ascii=False)+'。\n\n'
        '`context.json` に入力fingerprint、mainconf、元summary、設定・ソースのSHA256、'
        '`purity.json` に全モデルfitの状態・幅・積分・残差を保存しています。'
        '入力ROOTは読み取り専用、変更なしを確認しました。\n')
    archives = []
    for attempt_arg in args.archive_attempt:
        attempt = attempt_arg.resolve(strict=True)
        if attempt.parent != output.parent or not attempt.name.startswith(output.name+'_'):
            raise ValueError('Only explicitly named prior attempts for this output are archived.')
        prior_context = json.loads((attempt/'context.json').read_text())
        if prior_context['input'] != initial:
            raise ValueError('Prior attempt belongs to another input ROOT.')
        target = output/'revisions'/attempt.name[len(output.name)+1:]
        if target.exists() or attempt.is_symlink():
            raise ValueError('Archive target exists or prior attempt is a symlink.')
        records = {str(f.relative_to(attempt)):digest(f) for f in attempt.rglob('*') if f.is_file()}
        target.parent.mkdir(exist_ok=True)
        attempt.rename(target)
        if any(digest(target/name) != checksum for name, checksum in records.items()):
            raise ValueError('Preserved attempt hash changed during archival.')
        archives.append(dict(original=str(attempt), current=str(target), hashes=records,
            note='Historical absolute paths are retained; diagnostic attempt, not the final physics result.'))
    if archives:
        save(output/'revisions/relocation_manifest.json', dict(schemaVersion=1, attempts=archives))
    save(output/'completed.json', dict(phase='complete', input_unchanged=True, pages=pages, fits=counts,
        products={str(f.relative_to(output)):digest(f) for f in output.rglob('*')
                  if f.is_file() and f.name != 'completed.json' and 'revisions' not in f.relative_to(output).parts
                  and not f.suffix in ('.so', '.d', '.o')}, archives=archives))
    print('COMPLETE', args.sample, counts, flush=True)


if __name__ == '__main__':
    main()
