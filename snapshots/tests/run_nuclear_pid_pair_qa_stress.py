#!/usr/bin/env python3
"""Measure unique-cell capacity only; no physics/Pico input or production edits."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import yaml

PROJECT = Path(__file__).resolve().parents[1]


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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--unique-counts', type=int, nargs='+', required=True)
    parser.add_argument('--timeout-seconds', type=int, required=True)
    args = parser.parse_args()
    mainconf = args.mainconf.resolve(strict=True)
    config = yaml.safe_load(mainconf.read_text())
    ref = Path(config['nuclearPidPairQaHist'])
    hist = ref.resolve(strict=True) if ref.is_absolute() else (PROJECT / 'config' / ref).resolve(strict=True)
    if len(set(args.unique_counts)) != len(args.unique_counts) or any(n < 0 for n in args.unique_counts):
        raise ValueError('Require distinct nonnegative explicit cell counts.')
    output = args.output_dir.resolve()
    allowed = PROJECT / 'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/tests/nuclear_pid_pair_qa_20261017_001'
    if output.exists() or output.parent != allowed:
        raise ValueError('Require a new direct child of the designated study folder.')
    files = {Path(__file__).resolve(), mainconf, hist, PROJECT / 'script/singularity_root_batch.sh'}
    files.update(PROJECT / name for name in (
        'tests/stress_femto_nuclear_pid_pair_sparse.C', 'tests/run_stress_femto_nuclear_pid_pair_sparse.C',
        'StMaker/common/FemtoNuclearPairPidSparse.h', 'StMaker/common/FemtoNuclearPairPidSparse.cxx'))
    libraries = {PROJECT / 'lib' / name for name in ('libStarAnaConfig.so', 'libStRefMultCorr.so', 'libStCommon.so')}
    protected = {str(path): sha(path) for path in sorted(files | libraries)}
    output.mkdir(parents=True)
    for path in sorted(files):
        target = output / 'source' / path.relative_to(PROJECT)
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
    settings = dict(mainconf=str(mainconf), histogramYaml=str(hist), uniqueRegularCellCounts=args.unique_counts,
                    category='prePID/he4/SE/cent0/signalLow', additionalDistinctFlowCellsPerNonemptySample=6,
                    noPicoInput=True, noNewFit=True, physicalDistribution=False,
                    timingNote='Each fresh ROOT process includes runtime/library startup and full audits, but not ACLiC compilation.',
                    extrapolationNote='Uniform unique-cell stress is not a farm-distribution prediction or a guarantee of merged ROOT/hadd capacity.',
                    protected=protected, timeoutSeconds=args.timeout_seconds)
    write(output / 'context.json', settings)
    build = output / 'aclic'
    source = output / 'source/tests/stress_femto_nuclear_pid_pair_sparse.C'
    runner = output / 'source/tests/run_stress_femto_nuclear_pid_pair_sparse.C'

    def command(count, target, library):
        arguments = [str(mainconf), str(target), count, str(build), str(source), str(library)]
        call = str(runner) + '(' + ','.join(json.dumps(v) for v in arguments) + ')'
        return ['bash', str(PROJECT / 'script/singularity_root_batch.sh'), str(mainconf), 'root4star', '-l', '-b', '-q', call]

    compile_command = command(-1, output / 'compile_only', '')
    write(output / 'compile_command.json', compile_command)
    with (output / 'compile.log').open('x') as log:
        result = subprocess.run(compile_command, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout_seconds)
    if result.returncode:
        raise RuntimeError('Compile failed; retained compile.log.')
    candidates = list(build.rglob('stress_femto_nuclear_pid_pair_sparse_C.so'))
    if len(candidates) != 1:
        raise RuntimeError('Cannot uniquely resolve the once-compiled macro.')
    library = candidates[0]
    library_sha = sha(library)
    records = []
    for count in args.unique_counts:
        target = output / ('cells_%d' % count)
        time_file = output / ('cells_%d.time.txt' % count)
        cmd = ['/usr/bin/time', '-v', '-o', str(time_file)] + command(count, target, library)
        write(output / ('cells_%d_command.json' % count), cmd)
        started = time.monotonic()
        with (output / ('cells_%d.log' % count)).open('x') as log:
            result = subprocess.run(cmd, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, timeout=args.timeout_seconds)
        elapsed = time.monotonic() - started
        if result.returncode or 'PAIR_PID_SPARSE_STRESS_PASS' not in (output / ('cells_%d.log' % count)).read_text():
            raise RuntimeError('Stress sample failed; retained individual log.')
        report = json.loads((target / 'summary.json').read_text())
        peak = re.search(r'Maximum resident set size \(kbytes\):\s*(\d+)', time_file.read_text())
        if not peak or not report.get('passed'):
            raise RuntimeError('Missing sample readback or RSS report.')
        report.update(rootBytes=(target / 'capacity.root').stat().st_size, maxRssKiB=int(peak.group(1)), wallSeconds=elapsed)
        records.append(report)
        if sha(library) != library_sha:
            raise RuntimeError('Compiled stress macro changed during samples.')
    changed = [name for name, expected in protected.items() if sha(Path(name)) != expected]
    if changed:
        raise RuntimeError('Protected source/config/libraries changed: ' + ','.join(changed))
    write(output / 'completed.json', dict(executionComplete=True, protectedPostRunVerified=True,
          noPicoInput=True, physicalDistribution=False, compilationOnce=True, compiledMacroSha256=library_sha,
          samples=records, products={str(p.relative_to(output)): sha(p) for p in output.rglob('*') if p.is_file()}))
    print('NUCLEAR_PID_PAIR_QA_STRESS_PASS', output)


if __name__ == '__main__':
    main()
