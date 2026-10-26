#!/usr/bin/env python3
"""同一QA ROOTを意図的に2入力とする、非本番hadd互換性試験。"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess


def fingerprint(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return dict(path=str(path), size=path.stat().st_size, sha256=digest.hexdigest())


def save(path, value):
    path.write_text(json.dumps(value, indent=2, sort_keys=True) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mainconf', required=True, type=Path)
    parser.add_argument('--input', required=True, type=Path)
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args()
    repo = Path(__file__).resolve().parents[1]
    mainconf, source, output = (p.resolve() for p in (args.mainconf, args.input, args.output))
    if not mainconf.is_file() or not source.is_file():
        raise RuntimeError('Existing explicit mainconf and input ROOT are required')
    output.mkdir(parents=True, exist_ok=False)
    checker = repo / 'tests/check_femto_nuclear_pair_pid_qa_hadd.C'
    protected_paths = [source, mainconf, checker, Path(__file__).resolve()]
    for relative in (
            'StMaker/common/FemtoLambdaLegacy.cxx', 'StMaker/common/FemtoLambdaLegacy.h',
            'StMaker/common/FemtoNuclearPairPidSparse.cxx', 'StMaker/common/FemtoNuclearPairPidSparse.h',
            'StMaker/StFemtoMaker/StFemtoMaker.cxx', 'StMaker/StFemtoMaker/StFemtoMaker.h',
            'include/FemtoNuclearPidPairQaCandidate.h', 'Makefile',
            'lib/libStCommon.so', 'lib/libStFemtoMaker.so',
            'lib/libKFParticle.so', 'lib/libStKfParticleCommon.so',
            'config/hist/hist_auau3p85_femtoLambda_he4_he3_pid_pair_qa.yaml'):
        protected_paths.append(repo / relative)
    before = [fingerprint(path) for path in protected_paths]
    snapshots = output / 'source'
    snapshots.mkdir()
    frozen_checker = snapshots / checker.name
    shutil.copy2(str(checker), str(frozen_checker))
    shutil.copy2(__file__, str(snapshots / Path(__file__).name))
    merged = output / 'intentionally_duplicated_input.root'
    wrapper = ['bash', str(repo / 'script/singularity_root_batch.sh'), str(mainconf)]
    commands = [
        wrapper + ['hadd', str(merged), str(source), str(source)],
        wrapper + ['root4star', '-l', '-b', '-q',
                   '{}+("{}","{}")'.format(frozen_checker, source, merged)],
    ]
    context = dict(
        intentionalDuplicateInput=True, productionHadd=False,
        inputMultiplicity=2, inputRoot=str(source), outputRoot=str(merged),
        commands=commands, protectedBefore=before,
        metadataCountersAreNotAdditive=True,
        existingMergeMonitorAuditsTH1Only=True,
        startedUTC=datetime.now(timezone.utc).isoformat())
    save(output / 'context.json', context)
    try:
        for command, basename in zip(commands, ('hadd.log', 'readback.log')):
            with (output / basename).open('w') as log:
                result = subprocess.run(command, cwd=str(repo), stdout=log, stderr=subprocess.STDOUT)
            if result.returncode:
                raise RuntimeError('{} exited {}'.format(basename, result.returncode))
        match = re.search(r'NUCLEAR_PID_PAIR_QA_HADD_PASS qaTH2=(\d+) sparseEntries=([0-9.]+) sparseFilledBins=(\d+) checks=(\d+)',
                          (output / 'readback.log').read_text())
        if not match or int(match.group(1)) != 1080:
            raise RuntimeError('Missing full dense/sparse hadd validation marker')
        after = [fingerprint(path) for path in protected_paths]
        if before != after:
            raise RuntimeError('Immutable input/source/library changed during test')
        result = dict(executionComplete=True, intentionalDuplicateInput=True, productionHadd=False,
                      qaTH2=int(match.group(1)), sparseEntries=float(match.group(2)),
                      sparseFilledBins=int(match.group(3)), checks=int(match.group(4)),
                      immutableInputsVerified=True, protectedInputs=after,
                      completedUTC=datetime.now(timezone.utc).isoformat(),
                      note='THnSparse content/variance checked by coordinates; TNamed counters are not summed.')
        result['products'] = [fingerprint(path) for path in (
            merged, output / 'hadd.log', output / 'readback.log',
            output / 'context.json', frozen_checker, snapshots / Path(__file__).name)]
        save(output / 'completed.json', result)
        print('NUCLEAR_PID_PAIR_QA_HADD_COMPLETE qaTH2={} sparseEntries={} checks={}'.format(
            result['qaTH2'], result['sparseEntries'], result['checks']))
    except Exception as error:
        save(output / 'failed.json', dict(executionComplete=False, error=str(error),
                                         protectedAfter=[fingerprint(path) for path in protected_paths]))
        raise


if __name__ == '__main__':
    main()
