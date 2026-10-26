#!/usr/bin/env python3
"""Run fresh ROOT5 toy QA/CF/pool regressions with protected inputs and logs."""
import argparse
import hashlib
import json
from pathlib import Path
import shutil
import subprocess
import time
import yaml

PROJECT = Path(__file__).resolve().parents[1]
CASES = {
    'sparse': ('run_femto_nuclear_pair_pid_sparse.C', 'PAIR_PID_SPARSE_PASS'),
    'pair': ('run_femto_nuclear_pair_pid_qa.C', 'NUCLEAR_PAIR_PID_QA_PASS'),
    'pool': ('run_femto_nuclear_pair_pid_pool.C', 'NUCLEAR_PAIR_PID_POOL_PASS'),
}


def sha(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def write(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def config_path(ref):
    p = Path(ref)
    return p.resolve(strict=True) if p.is_absolute() else (PROJECT / 'config' / p).resolve(strict=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--on-mainconf', type=Path, required=True)
    p.add_argument('--off-mainconf', type=Path, required=True)
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--cases', default='sparse,pair,pool')
    args = p.parse_args()
    on = args.on_mainconf.resolve(strict=True)
    off = args.off_mainconf.resolve(strict=True)
    new = yaml.safe_load(on.read_text())
    old = yaml.safe_load(off.read_text())
    if 'nuclearPidPairQaHist' not in new or 'nuclearPidPairQaHist' in old:
        raise ValueError('Require explicit on/off QA references.')
    if {k: v for k, v in new.items() if k != 'nuclearPidPairQaHist'} != old:
        raise ValueError('Physical configuration differs between on/off.')
    cases = args.cases.split(',')
    if len(set(cases)) != len(cases) or any(c not in CASES for c in cases):
        raise ValueError('Choose distinct declared cases.')
    output = args.output_dir.resolve()
    allowed = PROJECT / 'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/nuclear_pid_pair_qa_20261017_001'
    if output.exists() or output.parent != allowed:
        raise ValueError('Require new child folder of the designated QA test study.')
    files = {Path(__file__).resolve(), on, off, PROJECT / 'tests/test_femto_lambda_entries.py', PROJECT / 'tests/test_femto_lambda_hist_contract.py', PROJECT / 'script/singularity_root_batch.sh'}
    files.update(PROJECT / 'tests' / f for f in (
        'femto_nuclear_pair_pid_sparse.C', 'run_femto_nuclear_pair_pid_sparse.C',
        'femto_nuclear_pair_pid_qa.C', 'run_femto_nuclear_pair_pid_qa.C',
        'femto_nuclear_pair_pid_pool.C', 'run_femto_nuclear_pair_pid_pool.C'))
    files.update(PROJECT / f for f in (
        'Makefile', 'include/FemtoNuclearPidPairQaCandidate.h',
        'StMaker/common/FemtoLambdaLegacy.h', 'StMaker/common/FemtoLambdaLegacy.cxx',
        'StMaker/common/FemtoNuclearPairPidSparse.h', 'StMaker/common/FemtoNuclearPairPidSparse.cxx',
        'StMaker/StFemtoMaker/StFemtoMaker.h', 'StMaker/StFemtoMaker/StFemtoMaker.cxx'))
    for mainconf in (new, old):
        files.update(config_path(v) for v in mainconf.values())
        cut = yaml.safe_load(config_path(mainconf['nuclearid']).read_text())
        if cut.get('dedxCalibration', {}).get('enabled'):
            files.add(config_path(cut['dedxCalibration']['calibrationFile']))
    libraries = [PROJECT / 'lib' / name for name in ('libStarAnaConfig.so', 'libStRefMultCorr.so', 'libStCommon.so', 'libStFemtoMaker.so')]
    protected = {str(v): sha(v) for v in sorted(files | set(libraries))}
    output.mkdir(parents=True)
    for path in sorted(files):
        dest = output / 'source' / path.relative_to(PROJECT)
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, dest)
    write(output / 'context.json', dict(onMainconf=str(on), offMainconf=str(off), cases=cases, protected=protected, noPicoInput=True, noNewFit=True, separateQaStore=True, sourceLibrariesFrozen=True))
    records = []
    for name, cmd in (
        ('static_entries', ['python3', '-B', '-m', 'unittest', 'discover', '-s', 'tests', '-p', 'test_femto_lambda_entries.py', '-v']),
        ('static_hist_contract', ['python3', '-B', 'tests/test_femto_lambda_hist_contract.py']),
    ):
        with (output / (name + '.log')).open('x') as log:
            r = subprocess.run(cmd, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, timeout=120)
        records.append(dict(case=name, exitCode=r.returncode))
        if r.returncode:
            raise RuntimeError(name + ' failed; see retained log.')
    for case in cases:
        runner, token = CASES[case]
        target = output / (case + '_001')
        parameters = [str(on), str(target)] if case == 'sparse' else [str(on), str(off), str(target)]
        call = 'tests/' + runner + '(' + ','.join(json.dumps(v) for v in parameters) + ')'
        cmd = ['bash', str(PROJECT / 'script/singularity_root_batch.sh'), str(on), 'root4star', '-l', '-b', '-q', call]
        write(output / (case + '_command.json'), cmd)
        started = time.monotonic()
        with (output / (case + '.log')).open('x') as log:
            result = subprocess.run(cmd, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, timeout=1200)
        record = dict(case=case, exitCode=result.returncode, seconds=time.monotonic()-started)
        write(output / (case + '_result.json'), record)
        records.append(record)
        if result.returncode or token not in (output / (case + '.log')).read_text():
            raise RuntimeError(case + ' failed; see retained log and new attempt folder.')
        record['summary'] = json.loads((target / 'summary.json').read_text())
    changed = [name for name, expected in protected.items() if sha(Path(name)) != expected]
    if changed:
        raise RuntimeError('Protected source/config/libraries changed: ' + ','.join(changed))
    write(output / 'completed.json', dict(executionComplete=True, protectedInputsVerified=True, noPicoInput=True, allDeclaredCasesPassed=True, runs=records, products={str(v.relative_to(output)): sha(v) for v in output.rglob('*') if v.is_file()}))
    print('NUCLEAR_PID_PAIR_QA_TOYS_PASS', output)


if __name__ == '__main__':
    main()
