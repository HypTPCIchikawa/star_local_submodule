#!/usr/bin/env python3
"""Snapshot and verify ROOT5 free-core/tail math only; no data fit or PID change."""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
OUTPUT_BASE = PROJECT / 'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_quality_20261019'


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(block)
    return digest.hexdigest()


def save(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def run(args):
    mainconf = args.mainconf.resolve(strict=True)
    if PROJECT not in mainconf.parents:
        raise ValueError('Require an explicit project-owned mainconf.')
    main = yaml.safe_load(mainconf.read_text())
    analysis = (PROJECT / 'config' / main['analysis']).resolve(strict=True)
    tags = yaml.safe_load(analysis.read_text())['starTag']
    if tuple(tags[key] for key in ('collision', 'year', 'productionTag', 'libraryTag')) != (
            'auau3p85', 'year2021', 'P24iy', 'SL24y'):
        raise ValueError('Use AuAu3p85 2021 P24iy / SL24y ROOT5 mainconf only.')
    output = args.output_dir.resolve()
    if output.exists() or output.parent != OUTPUT_BASE.resolve():
        raise ValueError('Require a NEW direct child of the quality output base; no overwrite.')
    sources = [Path(__file__).resolve(), HERE / 'SelfTestFreeCoreTail.C',
               HERE / 'FreeGaussianCoreTail.h', HERE / 'NuclearCoreTailFitCommon.h',
               HERE / 'free_coretail_math_settings.tsv',
               PROJECT / 'tools/nuclear_pid_tof/NuclearFitCommon.h',
               PROJECT / 'script/singularity_root_batch.sh',
               PROJECT / 'script/analysis_info_helper.py', mainconf, analysis]
    protected = {str(path.resolve(strict=True)): sha(path) for path in sources}
    output.mkdir(parents=True)
    source_root = output / 'source'
    source_root.mkdir()
    for path in sources:
        path = path.resolve(strict=True)
        relative = path.relative_to(PROJECT)
        target = source_root / relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(path, target)
        if sha(target) != protected[str(path)]:
            raise ValueError('Source changed while snapshotting: ' + str(path))
    snapshot_here = source_root / HERE.relative_to(PROJECT)
    expression = str(snapshot_here / 'SelfTestFreeCoreTail.C') + '+(' + ','.join(
        json.dumps(str(path)) for path in (
            snapshot_here / 'free_coretail_math_settings.tsv', output / 'selftest.json')) + ')'
    command = ['bash', str(PROJECT / 'script/singularity_root_batch.sh'), str(mainconf),
               'root4star', '-l', '-b', '-q', expression]
    save(output / 'context.json', dict(
        schemaVersion=1, createdUtc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        runtimeMainconf=str(mainconf), protectedSources=protected,
        sourceHierarchyPreserved=True, noRealDataInput=True, noFit=True,
        productionChanged=False,
        scope='Free Gaussian core/tail normalization, independent parameters, display folding, post-fit screens only.'))
    save(output / 'command.json', command)
    try:
        with (output / 'root.log').open('x') as log:
            result = subprocess.run(command, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT,
                                    timeout=300)
        log_text = (output / 'root.log').read_text()
        if result.returncode or 'NUCLEAR_FREE_CORE_TAIL_SELFTEST_PASS' not in log_text:
            raise RuntimeError('ROOT5 free-core/tail selftest did not complete; retained root.log.')
        report = json.loads((output / 'selftest.json').read_text())
        checks = report.get('checks')
        if (report.get('passed') is not True or not isinstance(checks, int) or isinstance(checks, bool)
                or checks <= 0 or report.get('fitInvoked') is not False
                or report.get('productionChanged') is not False
                or report.get('noRealDataInput') is not True):
            raise ValueError('Invalid/incomplete mathematical selftest report.')
        for path, expected in protected.items():
            if sha(Path(path)) != expected:
                raise ValueError('Protected source/config changed during test: ' + path)
            snapshot = source_root / Path(path).relative_to(PROJECT)
            if sha(snapshot) != expected:
                raise ValueError('Consumed source snapshot changed during test: ' + str(snapshot))
        products = {str(path.relative_to(output)): sha(path) for path in output.rglob('*')
                    if path.is_file() and path.suffix not in ('.so', '.d', '.pcm')}
        save(output / 'completed.json', dict(
            executionComplete=True, passed=True, checks=checks, noFit=True,
            noRealDataInput=True, productionChanged=False, protectedSourcesVerified=True,
            consumedSnapshotsVerified=True, sourceHierarchyPreserved=True,
            groupedPurityAvailable=False, physicsValidationComplete=False, products=products))
        print('FREE_CORE_TAIL_MATH_SELFTEST_COMPLETE', output, 'checks=' + str(checks), flush=True)
    except Exception as error:
        save(output / 'failed.json', dict(executionComplete=False, noFit=True, error=str(error)))
        raise


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mainconf', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    run(parser.parse_args())
