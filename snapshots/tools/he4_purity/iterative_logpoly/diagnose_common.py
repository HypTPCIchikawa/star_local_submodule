#!/usr/bin/env python3
"""Probe a completed frozen common point; no shape fitting or PID publication."""
import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def tsv(path, rows):
    path.write_text(''.join('\t'.join(format(v, '.17g') if isinstance(v, float) else str(v)
                                    for v in row)+'\n' for row in rows))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, required=True)
    parser.add_argument('--phase', required=True)
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--variant', action='append', required=True,
                        help='Explicit name:tolerance:iterations, repeat for each diagnostic policy')
    parser.add_argument('--timeout-seconds', type=int, required=True)
    args = parser.parse_args()
    study = args.study.resolve(strict=True)
    phase = (study/args.phase).resolve(strict=True)
    if study not in phase.parents:
        raise ValueError('Phase must be inside the completed source study')
    completed = json.loads((study/'completed.json').read_text())
    context = json.loads((study/'context.json').read_text())
    result_path = phase/'fit_result.json'
    result = json.loads(result_path.read_text())
    if not completed.get('executionComplete') or result['phase'] != 'common':
        raise ValueError('Only a completed common-phase source is allowed')
    for path in (result_path, phase/'phase.root', phase/'settings.tsv'):
        if digest(path) != completed['products'][str(path.relative_to(study))]:
            raise ValueError('Source product hash mismatch '+str(path))
    cfg = context['settings']
    identity = cfg['sample']
    if (identity['collision'], identity['year'], identity['production']) != ('auau3p85', 2021, 'P24iy'):
        raise ValueError('Only 3.85 GeV / 2021 P24iy is authorized')
    mainconf = args.mainconf.resolve(strict=True)
    if mainconf != (PROJECT/identity['mainconf']).resolve(strict=True):
        raise ValueError('Explicit mainconf differs from frozen runtime')
    cache = Path(context['cache']).resolve(strict=True)
    output = args.output_dir.resolve()
    if output.exists() or output.parent != study.parent:
        raise ValueError('Choose a nonexistent sibling diagnostic directory')
    variants = []
    for spec in args.variant:
        name, tolerance, iterations = spec.split(':')
        tolerance, iterations = float(tolerance), int(iterations)
        if not name.isidentifier() or not 0 < tolerance < 1 or iterations < 1:
            raise ValueError('Invalid explicit diagnostic policy')
        variants.append([name, tolerance, iterations])
    if len({v[0] for v in variants}) != len(variants) or args.timeout_seconds <= 0:
        raise ValueError('Duplicate variant or invalid timeout')
    output.mkdir()
    protected = {str(path): digest(path) for path in
                 (cache, result_path, phase/'phase.root', phase/'settings.tsv',
                  study/'completed.json', study/'context.json', mainconf)}
    source_relative = ('tools/he4_purity/iterative_logpoly/IterativeHeliumFit.C',
                       'tools/he4_purity/global_logpoly/GlobalHeliumFit.C',
                       'StMaker/common/NuclearIdDeDxVsMom.h')
    snapshots = {}
    for relative in source_relative:
        src = study/'source'/relative
        expected = context['snapshots'][str(src)]
        if digest(src) != expected:
            raise ValueError('Frozen compiled source changed '+str(src))
        dst = output/'source'/relative
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        snapshots[str(dst)] = digest(dst)
        protected[str(src)] = expected
    for src in (HERE/'DiagnoseCommonPoint.C', Path(__file__).resolve()):
        dst = output/'source'/src.relative_to(PROJECT)
        dst.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dst)
        snapshots[str(dst)] = digest(dst)
        protected[str(src)] = digest(src)
    wrapper = PROJECT/'script/singularity_root_batch.sh'
    protected[str(wrapper)] = digest(wrapper)
    shutil.copy2(phase/'settings.tsv', output/'settings.tsv')
    tsv(output/'variants.tsv', variants)
    tsv(output/'targets.tsv', [['coefficients']+result['sharedCoefficients']])
    seeds, amps = [], []
    for row in result['slices']:
        seeds.append([row['index'], row['lo'], row['hi'], row['anchorMean'], row['anchorSigma'],
                      row['oldChosen'], int(row['oldValid']), int(row['oldQuality']),
                      row['backgroundModel'], int(row['t']['included'])]+
                     [row[key][field] for key in ('t', 'he3', 'he4')
                      for field in ('mean', 'sigma', 'windowYield')]+
                     [row['backgroundWindowYield'], row['backgroundSlope']])
        amps.append([row['index']]+[row[key]['windowYield'] for key in ('t', 'he3', 'he4')]+
                    [row['backgroundWindowYield']])
    tsv(output/'seeds.tsv', seeds)
    tsv(output/'amplitudes.tsv', amps)
    (output/'diagnostic_settings.yaml').write_text(yaml.safe_dump(
        dict(schemaVersion=1, sourceStudy=str(study), sourcePhase=args.phase,
             shapeRefit=False, probeStepFrom='commonCoefficientTolerance', variants=variants,
             timeoutSeconds=args.timeout_seconds), sort_keys=False))
    for name in ('settings.tsv', 'variants.tsv', 'targets.tsv', 'seeds.tsv', 'amplitudes.tsv',
                 'diagnostic_settings.yaml'):
        snapshots[str(output/name)] = digest(output/name)
    save(output/'context.json', dict(schemaVersion=1, sourceStudy=str(study), sourcePhase=args.phase,
         productionEnabled=False, protected=protected, consumedSnapshots=snapshots,
         createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))

    def preserve():
        for path, expected in {**protected, **snapshots}.items():
            if digest(Path(path)) != expected:
                raise ValueError('Protected source/input changed '+path)

    macro = output/'source/tools/he4_purity/iterative_logpoly/DiagnoseCommonPoint.C'
    call = str(macro)+'++('+','.join(json.dumps(str(p)) for p in
        (cache, phase/'phase.root', output, output/'settings.tsv', output/'seeds.tsv',
         output/'targets.tsv', output/'amplitudes.tsv', output/'variants.tsv'))+')'
    command = ['bash', str(wrapper), str(mainconf), 'root4star', '-l', '-b', '-q', call]
    save(output/'command.json', command)
    preserve()
    print('DIAGNOSTIC_RUN', output, flush=True)
    with (output/'diagnostic.log').open('w') as log:
        run = subprocess.run(command, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT,
                             timeout=args.timeout_seconds)
    if run.returncode or 'COMMON_POINT_DIAGNOSTIC_COMPLETE' not in (output/'diagnostic.log').read_text():
        raise RuntimeError('Diagnostic ROOT process failed; see '+str(output/'diagnostic.log'))
    data = json.loads((output/'diagnostic.json').read_text(),
                      parse_constant=lambda v: (_ for _ in ()).throw(ValueError(v)))
    if len(data['variants']) != len(variants) or any(len(v['cases']) != 29 for v in data['variants']):
        raise ValueError('Unexpected diagnostic coverage')
    preserve()
    save(output/'completed.json', dict(executionComplete=True, productionEnabled=False,
         physicsValidationComplete=False, products={name:digest(output/name)
         for name in ('diagnostic.json', 'diagnostic.log', 'diagnostic_settings.yaml')},
         completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))
    print('DIAGNOSTIC_COMPLETE', output, flush=True)


if __name__ == '__main__':
    main()
