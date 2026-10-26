#!/usr/bin/env python3
"""Run explicit, local QA-off/on analyses and compare closed ROOT outputs."""
import argparse
import hashlib
import json
import os
from pathlib import Path
import shutil
import subprocess
import yaml


def sha(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            h.update(block)
    return h.hexdigest()


def fingerprint(path):
    stat = path.stat()
    return dict(size=stat.st_size, mtimeNs=stat.st_mtime_ns, inode=stat.st_ino, device=stat.st_dev)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('main-on', 'main-off', 'input', 'output-dir'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--events', required=True, type=int)
    parser.add_argument('--expected-new-th2', required=True, type=int)
    parser.add_argument('--timeout', required=True, type=int)
    args = parser.parse_args()
    project = Path(__file__).resolve().parents[1]
    on_path, off_path, input_path = (p.resolve() for p in (args.main_on, args.main_off, args.input))
    on, off = (yaml.safe_load(p.read_text()) for p in (on_path, off_path))
    if 'nuclearPidPairQaHist' not in on or 'nuclearPidPairQaHist' in off:
        raise ValueError('Require explicit QA-on and QA-off mainconfs')
    if {k: v for k, v in on.items() if k != 'nuclearPidPairQaHist'} != off:
        raise ValueError('On/off differ in physics selection')
    if args.events <= 0 or args.expected_new_th2 <= 0 or args.timeout <= 0:
        raise ValueError('Explicit positive test parameters required')
    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=False)
    source = output / 'source'
    source.mkdir()
    paths = {on_path, off_path, Path(__file__).resolve(), project / 'tests/check_femto_nuclear_pair_pid_qa_output.C'}
    paths.update((project / 'config' / value).resolve() for value in on.values())
    nuclear = yaml.safe_load((project / 'config' / on['nuclearid']).read_text())
    calibration = nuclear.get('dedxCalibration', {})
    if calibration.get('enabled'):
        paths.add((project / 'config' / calibration['calibrationFile']).resolve())
    paths.update(project.glob('StMaker/common/FemtoNuclearPairPidSparse.*'))
    paths.update(project.glob('StMaker/common/FemtoLambdaLegacy.*'))
    paths.update(project.glob('StMaker/common/*.h'))
    paths.update(project.glob('StMaker/common/*.cxx'))
    paths.update(project.glob('StMaker/StFemtoMaker/StFemtoMaker.*'))
    paths.update(project.glob('lib/*.so'))
    paths.add(project / 'Makefile')
    paths.update(project / name for name in (
        'script/generate_femto_nuclear_pid_pair_qa_hist.py', 'tools/femto_nuclear_pid_pair_qa_settings.yaml',
        'include/FemtoCandidate.h', 'include/ConfigManager.h', 'include/cuts/NuclearIdCutConfig.h',
        'analysis/anaFemtoLambda_4He.C', 'analysis/run_anaFemtoLambda_4He.C', 'analysis/FemtoLambdaRunSupport.h',
        'script/singularity_run_anaFemtoLambda_4He.sh', 'script/run_femtoLambda_common.sh',
        'script/femtoLambda_environment.sh', 'script/singularity_root_batch.sh'))
    protected = {str(p): sha(p) for p in paths if p.is_file()}
    for path in paths:
        if path.is_file() and path.suffix != '.so':
            dest = source / path.relative_to(project)
            dest.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(path, dest)
    before_input = fingerprint(input_path)
    commands = []
    env = os.environ.copy()
    env['LC_ALL'] = 'C'
    for label, main_path in (('off', off_path), ('on', on_path)):
        root = output / (label + '.root')
        command = ['/usr/bin/time', '-v', '-o', str(output / (label + '.time.txt')),
                   'bash', str(project / 'script/singularity_run_anaFemtoLambda_4He.sh'),
                   str(main_path), str(input_path), str(root), 'nuclear_pair_pid_qa_smoke', str(args.events)]
        commands.append(command)
        with (output / (label + '.log')).open('w') as log:
            subprocess.run(command, cwd=project, stdout=log, stderr=subprocess.STDOUT,
                           timeout=args.timeout, env=env, check=True)
    macro = source / 'tests/check_femto_nuclear_pair_pid_qa_output.C'
    call = str(macro) + '+(' + ','.join(json.dumps(str(output / name)) for name in ('off.root', 'on.root', 'readback.json'))
    call += ',' + str(args.expected_new_th2) + ',' + str(args.events) + ')'
    command = ['bash', str(project / 'script/singularity_root_batch.sh'), str(on_path), 'root4star', '-l', '-b', '-q', call]
    commands.append(command)
    with (output / 'readback.log').open('w') as log:
        subprocess.run(command, cwd=project, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout, env=env, check=True)
    if any(sha(p) != value for p, value in protected.items()) or fingerprint(input_path) != before_input:
        raise ValueError('Protected source/library/input changed during test')
    report = dict(complete=True, physicsConfigIdentical=True, input=str(input_path), inputFingerprint=before_input,
                  requestedEventsPerRun=args.events, commands=commands, protected=protected,
                  readback=json.loads((output / 'readback.json').read_text()),
                  timingNote='One short run per mode; includes process setup and possible first-run ACLiC/cache effects, not a benchmark',
                  products={str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file()})
    (output / 'completed.json').write_text(json.dumps(report, indent=2) + '\n')
    print('NUCLEAR_PID_PAIR_SMOKES_COMPLETE', output)


if __name__ == '__main__':
    main()
