#!/usr/bin/env python3
"""Reproduce KFParticle Lambda--4He QA figures without modifying input ROOTs."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import datetime
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys
import threading

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


def independent_subdirectories(out):
    """Only explicitly named, untouched sibling studies / history are allowed."""
    names = {'revisions', 'mass_halves', 'mass_halves20MeV'}
    return [p for p in out.iterdir() if p.name not in names or not p.is_dir() or p.is_symlink()]


def archive_previous(out, label, source, snapshot=None):
    """Move only verified prior summary products; retain inputs and other studies.

    Historical JSON paths intentionally remain as originally recorded. The
    relocation manifest records current paths and hashes, and config snapshots
    preserve the settings used by the previous report.
    """
    if not re.fullmatch(r'[A-Za-z0-9][A-Za-z0-9_.-]*', label):
        raise ValueError('Revision label must be a single safe directory name.')
    if out.name != 'pdf':
        raise ValueError('Revision archival is restricted to a dedicated pdf directory.')
    state = json.loads((out/'run_state.json').read_text())
    manifest = json.loads((out/'macro_manifest.json').read_text())
    report = json.loads((out/'report_metadata.json').read_text())
    if state['phase'] != 'complete' or report['phase'] != 'complete':
        raise RuntimeError('Only a completed prior report may be archived.')
    if state['provenance']['input'] != fingerprint(source) or report['inputRoot'] != str(source):
        raise RuntimeError('Prior report input does not match this unchanged ROOT.')
    snapshots = {}
    for name, expected in state['provenance']['configs'].items():
        path = Path(name)
        candidate = path
        if digest(path) != expected:
            if snapshot is None or path.parent != HERE:
                raise RuntimeError('Prior config changed without its verified snapshot: '+name)
            candidate = snapshot/path.name
        if not candidate.is_file() or digest(candidate) != expected:
            raise RuntimeError('Prior config snapshot hash mismatch: '+name)
        snapshots[path.relative_to(PROJECT)] = candidate
    checks = {out/'inspect.json': report['inspect_sha256'],
              out/'macro_manifest.json': report['manifest_sha256'],
              out/report['report']['name']: report['report']['sha256'],
              out/Path(report['report']['name']).with_suffix('.tex'): report['tex_sha256']}
    for name, record in report['figures'].items():
        checks[out/name] = record['sha256']
    for record in report['fit_audits']:
        checks[Path(record['path'])] = record['sha256']
    for macro in manifest['macros']:
        checks[Path(macro['source'])] = macro['source_sha256']
        checks[Path(macro['generated'])] = macro['generated_sha256']
    checks[out/'macros/KFSummarySupport.h'] = manifest['support']['sha256']
    generator = HERE/'build_report.py'
    if digest(generator) != report['generator_sha256']:
        generator = snapshot/'build_report.py' if snapshot is not None else generator
    if not generator.is_file() or digest(generator) != report['generator_sha256']:
        raise RuntimeError('Prior report generator has no matching saved snapshot.')
    snapshots[HERE.relative_to(PROJECT)/'build_report.py'] = generator
    for path, expected in checks.items():
        if not path.is_file() or path.is_symlink() or digest(path) != expected:
            raise RuntimeError('Prior report/source hash mismatch: '+str(path))
    names = {'logs','macros','run_state.json','inspect.json','Lambda_InvMass.png',
             'macro_manifest.json','report_metadata.json','README.md'} | set(report['figures'])
    stem = Path(report['report']['name']).stem
    names.update(stem+suffix for suffix in ('.pdf','.tex','.log','.aux','.out'))
    items = [out/name for name in sorted(names) if (out/name).exists()]
    unowned = [p.name for p in independent_subdirectories(out) if p.name not in names]
    if unowned:
        raise RuntimeError('Unowned top-level products prevent archival: '+', '.join(unowned))
    records = {}
    for item in items:
        for path in ([item] if item.is_file() else sorted(item.rglob('*'))):
            if path.is_symlink():
                raise RuntimeError('Refusing symlink in archived report: '+str(path))
            if path.is_file():
                records[str(path.relative_to(out))] = dict(sha256=digest(path), bytes=path.stat().st_size)
    archive = out/'revisions'/label
    archive.mkdir(parents=True, exist_ok=False)
    relocation = dict(schemaVersion=1, phase='prepared', created=timestamp(),
        source_directory=str(out), archive_directory=str(archive), input=fingerprint(source),
        old_settings=manifest['settings'], files=records,
        untouched_subdirectories=['mass_halves','mass_halves20MeV'],
        note='Historical JSON paths retained; use this relocation record for archived paths.')
    save(archive/'revision_manifest.json', relocation)
    for relative, old in snapshots.items():
        target = archive/'source_snapshots'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(old),str(target))
    if snapshot is not None:
        for name in ('prepare_macros.py','run_summary.py','validate_outputs.py'):
            old = snapshot/name
            if old.is_file():
                target = archive/'source_snapshots'/HERE.relative_to(PROJECT)/name
                target.parent.mkdir(parents=True, exist_ok=True)
                shutil.copy2(str(old),str(target))
    for item in items:
        item.rename(archive/item.name)
    for relative, expected in records.items():
        path = archive/relative
        if not path.is_file() or path.stat().st_size != expected['bytes'] or digest(path) != expected['sha256']:
            raise RuntimeError('Archived product failed verification: '+str(path))
    if fingerprint(source) != relocation['input']:
        raise RuntimeError('Input ROOT changed during report archival.')
    relocation.update(phase='complete', finished=timestamp(), moved_files=len(records))
    save(archive/'revision_manifest.json', relocation)
    (archive/'REVISION_README.md').write_text(
        '# 旧summaryの保存版\n\n'
        'summaryの切替前の図・レポート・ログを、削除せず退避しました。保存版のrebin・表示条件は '
        '`revision_manifest.json` の `old_settings` を参照してください。\n'
        '旧設定と生成コードは `source_snapshots/`、全移動ファイルのhashと移動対応は '
        '`revision_manifest.json` に保存しています。\n\n'
        '旧JSONやマクロ内の絶対パスは歴史的記録としてそのままです。旧マクロをここから直接実行すると '
        '現行出力先を参照するため、再計算時は別の出力ディレクトリへ新しく生成してください。\n'
        '元ROOTと新しい質量low/high比較のサブディレクトリは変更していません。\n')
    print(timestamp(), 'ARCHIVED', archive, len(records), 'files', flush=True)
    return archive


def restore_revision(out, label, source):
    """Restore verified archived products while keeping the archive unchanged."""
    archive = (out/'revisions'/label).resolve(strict=True)
    if archive.parent != (out/'revisions').resolve(strict=True):
        raise ValueError('Revision must be directly under this report revisions directory.')
    record = json.loads((archive/'revision_manifest.json').read_text())
    if record['phase'] != 'complete' or record['input'] != fingerprint(source):
        raise RuntimeError('Completed revision / unchanged input required for restoration.')
    if independent_subdirectories(out):
        raise RuntimeError('Archive current products before restoring an older revision.')
    top_names = set()
    for relative, expected in record['files'].items():
        path = archive/relative
        if (path.is_symlink() or not path.is_file() or path.stat().st_size != expected['bytes']
                or digest(path) != expected['sha256']):
            raise RuntimeError('Archive hash mismatch before restoration: '+str(path))
        top_names.add(Path(relative).parts[0])
    for name in sorted(top_names):
        old, new = archive/name, out/name
        if new.exists():
            raise RuntimeError('Refusing to replace existing product: '+str(new))
        if old.is_dir():
            shutil.copytree(str(old),str(new),copy_function=shutil.copy2)
        else:
            shutil.copy2(str(old),str(new))
    for relative, expected in record['files'].items():
        if digest(out/relative) != expected['sha256']:
            raise RuntimeError('Restored file hash mismatch: '+relative)
    if fingerprint(source) != record['input']:
        raise RuntimeError('Input ROOT changed during restoration.')
    saved = dict(phase='restored', restored=timestamp(), output_directory=str(out),
                 input_unchanged=True, files_verified=len(record['files']),
                 archived_originals_retained=True,
                 note='Validate historical config/generator hashes using source_snapshots.')
    save(archive/'restoration.json',saved)
    print(timestamp(),'RESTORED',out,len(record['files']),'files',flush=True)
    return saved


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
    parser.add_argument('--archive-previous', metavar='LABEL',
                        help='Verify and retain a completed prior summary under pdf/revisions/LABEL before a fresh run.')
    parser.add_argument('--archive-only',action='store_true',
                        help='With --archive-previous, stop after verified archival, before new plotting.')
    parser.add_argument('--previous-tool-snapshot', type=Path,
                        help='Verified pre-edit settings/generator copies when archiving a revised report.')
    parser.add_argument('--workers',type=int,choices=[1,2],default=1,
                        help='Independent ROOT macro processes; state writes stay serialized.')
    args = parser.parse_args()
    if args.only_macros and not args.resume:
        parser.error('--only-macros requires --resume')
    if args.archive_previous and args.resume:
        parser.error('--archive-previous starts a new run; it cannot be combined with --resume')
    if args.archive_only and not args.archive_previous:
        parser.error('--archive-only requires --archive-previous')
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
    previous_revision = None
    if args.archive_previous:
        previous_revision = archive_previous(out,args.archive_previous,source,
            args.previous_tool_snapshot.resolve(strict=True) if args.previous_tool_snapshot else None)
        if args.archive_only:
            return
    state_file = out / 'run_state.json'
    if state_file.exists():
        if not args.resume:
            raise RuntimeError('Existing run; use --resume after inspecting %s' % state_file)
        state = json.loads(state_file.read_text())
        if state['provenance'] != provenance:
            raise RuntimeError('Input/config changed; choose a NEW output directory.')
    else:
        if out.exists() and independent_subdirectories(out):
            raise RuntimeError('Nonempty output without run_state; choose a new directory.')
        out.mkdir(parents=True, exist_ok=True)
        state = dict(provenance=provenance, started=timestamp(), stages={})
        if previous_revision:
            state['previous_revision'] = str(previous_revision)
    (out / 'logs').mkdir(exist_ok=True)
    (out / 'macros').mkdir(exist_ok=True)
    state.update(phase='running', updated=timestamp())
    save(state_file, state)
    state_lock=threading.RLock()

    def stage(name, signature, command, products):
        with state_lock:
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
        with state_lock:
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
        def macro_stage(item):
            macro = Path(item['generated'])
            if selected and macro.stem not in selected:
                if state['stages'].get(macro.stem, {}).get('phase') != 'complete':
                    raise RuntimeError('Cannot skip an incomplete stage: ' + macro.stem)
                return
            products = [out / name for name in item['expected_pdfs']]
            signature = digest(macro) + ':' + digest(out / 'macros' / 'KFSummarySupport.h')
            stage(macro.stem, signature, root_prefix + [str(macro) + '+()'], products)
        with ThreadPoolExecutor(max_workers=args.workers) as pool:
            list(pool.map(macro_stage,manifest['macros']))
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
