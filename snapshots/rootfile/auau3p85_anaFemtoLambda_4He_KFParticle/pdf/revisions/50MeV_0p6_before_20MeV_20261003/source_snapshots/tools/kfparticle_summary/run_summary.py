#!/usr/bin/env python3
"""Reproduce KFParticle Lambda--4He QA figures without modifying input ROOTs."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(path):
    st = path.stat()
    return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns,
                device=st.st_dev, inode=st.st_ino)


def timestamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def save(path, data):
    tmp = path.with_suffix(path.suffix + '.tmp')
    tmp.write_text(json.dumps(data, indent=2, ensure_ascii=False) + '\n')
    tmp.replace(path)


def run(command, log, timeout=1800):
    print(timestamp(), 'RUN', log.name, flush=True)
    with log.open('w') as handle:
        handle.write('COMMAND: ' + json.dumps(command) + '\n')
        handle.flush()
        result = subprocess.run(command, cwd=str(PROJECT), stdout=handle,
                                stderr=subprocess.STDOUT, timeout=timeout)
    if result.returncode:
        raise RuntimeError('Command failed (%s): %s' % (result.returncode, log))
    content = log.read_text(errors='replace')
    for token in ('Error in <ACLiC>', 'fatal error:', 'segmentation violation',
                  'Error in <TFile::TFile>', 'Error: Function ', 'Syntax Error'):
        if token in content:
            raise RuntimeError('ROOT diagnostic %r; see %s' % (token, log))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energy', required=True, choices=['13p5', '3p85'])
    parser.add_argument('--mainconf', required=True, type=Path)
    parser.add_argument('--input-root', required=True, type=Path)
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--settings', type=Path, default=HERE / 'summary_settings.yaml',
                        help='Explicit downstream plot/fit YAML; never changes production cuts.')
    parser.add_argument('--resume', action='store_true',
                        help='Resume only this same input/config; reuse successful stages.')
    parser.add_argument('--only-macros', nargs='+',
                        help='With --resume, redraw only named stages; all report inputs must still exist.')
    args = parser.parse_args()
    if args.only_macros and not args.resume:
        parser.error('--only-macros requires --resume')
    source = args.input_root.resolve(strict=True)
    mainconf = args.mainconf.resolve(strict=True)
    summary_settings = args.settings.resolve(strict=True)
    plotting = yaml.safe_load(summary_settings.read_text())
    out = (args.output_dir or source.parent / 'pdf').resolve()
    if out == source.parent or source in out.parents:
        raise ValueError('Output must be a separate directory, normally ROOT_DIR/pdf.')
    mainconfig = yaml.safe_load(mainconf.read_text())
    maker = (PROJECT / 'config' / mainconfig['maker']).resolve(strict=True)
    settings = yaml.safe_load(maker.read_text())
    if 'maker' in settings:
        settings = settings['maker']
    gate = [float(settings[key]) for key in ('lambdaSignalMean', 'lambdaSignalSigma',
                                             'lambdaSignalNSigma', 'lambdaSidebandOuterFactor')]
    config_hashes = {str(mainconf): digest(mainconf)}
    config_hashes[str(summary_settings)] = digest(summary_settings)
    for rel in mainconfig.values():
        if isinstance(rel, str) and rel.endswith(('.yaml', '.yml')):
            cfg = (PROJECT / 'config' / rel).resolve(strict=True)
            config_hashes[str(cfg)] = digest(cfg)
    provenance = dict(energy=args.energy, input=fingerprint(source),
                      configs=config_hashes, gate=gate)
    state_file = out / 'run_state.json'
    if state_file.exists():
        if not args.resume:
            raise RuntimeError('Existing run; use --resume after inspecting %s' % state_file)
        state = json.loads(state_file.read_text())
        if state['provenance'] != provenance:
            raise RuntimeError('Input/config changed; choose a NEW output directory.')
    else:
        if out.exists() and any(out.iterdir()):
            raise RuntimeError('Nonempty output without run_state; choose a new directory.')
        out.mkdir(parents=True, exist_ok=True)
        state = dict(provenance=provenance, started=timestamp(), stages={})
    (out / 'logs').mkdir(exist_ok=True)
    (out / 'macros').mkdir(exist_ok=True)
    state.update(phase='running', updated=timestamp())
    save(state_file, state)

    def stage(name, signature, command, products):
        previous = state['stages'].get(name, {})
        if (previous.get('signature') == signature and previous.get('phase') == 'complete'
                and all(p.is_file() and p.stat().st_size for p in products)):
            print(timestamp(), 'REUSE', name, flush=True)
            return
        state['stages'][name] = dict(phase='running', signature=signature, started=timestamp())
        save(state_file, state)
        run(command, out / 'logs' / (name + '.log'))
        missing = [str(p) for p in products if not p.is_file() or not p.stat().st_size]
        if missing:
            raise RuntimeError('Missing outputs: ' + ', '.join(missing))
        state['stages'][name].update(phase='complete', finished=timestamp(),
                                   products=[str(p) for p in products])
        save(state_file, state)

    try:
        inspect_source = HERE / 'InspectKFSummary.C'
        inspect_copy = out / 'macros' / inspect_source.name
        shutil.copy2(str(inspect_source), str(inspect_copy))
        root_prefix = ['bash', str(PROJECT / 'script/singularity_root_batch.sh'),
                       str(mainconf), 'root4star', '-l', '-b', '-q']
        call = str(inspect_copy) + '+(' + ','.join(
            [json.dumps(str(source)), json.dumps(str(out))] +
            [repr(x) for x in gate + plotting['mass_fit']['range']]) + ')'
        stage('inspect', digest(inspect_source), root_prefix + [call], [out / 'inspect.json'])
        info = json.loads((out / 'inspect.json').read_text())
        purity = info['purity_for_cfsub']
        if not isinstance(purity, (float, int)) or not 0 < purity < 1:
            raise RuntimeError('No physical inclusive-purity proxy; inspect diagnostics before CFsub.')
        prepare = [sys.executable, str(HERE / 'prepare_macros.py'), '--energy', args.energy,
                   '--input-root', str(source), '--output-dir', str(out), '--purity', str(purity),
                   '--maker-yaml', str(maker), '--settings', str(summary_settings)]
        run(prepare, out / 'logs' / 'prepare.log')
        manifest = json.loads((out / 'macro_manifest.json').read_text())
        selected = set(args.only_macros or [])
        available = {Path(item['generated']).stem for item in manifest['macros']}
        if selected - available:
            raise ValueError('Unknown macro stages: ' + ', '.join(sorted(selected - available)))
        for item in manifest['macros']:
            macro = Path(item['generated'])
            if selected and macro.stem not in selected:
                if state['stages'].get(macro.stem, {}).get('phase') != 'complete':
                    raise RuntimeError('Cannot skip an incomplete stage: ' + macro.stem)
                continue
            products = [out / name for name in item['expected_pdfs']]
            signature = digest(macro) + ':' + digest(out / 'macros' / 'KFSummarySupport.h')
            stage(macro.stem, signature, root_prefix + [str(macro) + '+()'], products)
        report = HERE / 'build_report.py'
        run([sys.executable, str(report), '--output-dir', str(out), '--energy', args.energy],
            out / 'logs' / 'build_report.log')
        if fingerprint(source) != provenance['input']:
            raise RuntimeError('Input ROOT fingerprint changed during reporting.')
        state.pop('error', None)
        state.update(phase='complete', finished=timestamp())
    except Exception as exc:
        state.update(phase='failed', error=str(exc), updated=timestamp())
        raise
    finally:
        save(state_file, state)
    print(timestamp(), 'COMPLETE', out, flush=True)


if __name__ == '__main__':
    main()
