#!/usr/bin/env python3
"""Snapshot the actual TPC selector before an authorized production-profile change."""
import argparse
from datetime import datetime, timezone
import hashlib
import importlib.util
import json
from pathlib import Path
import shutil

PROJECT = Path(__file__).resolve().parents[2]


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2,
                                    allow_nan=False) + '\n')


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mainconf', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    mainconf = args.mainconf.resolve(strict=True)
    output = args.output_dir.resolve()
    if PROJECT not in mainconf.parents or PROJECT not in output.parents or output.exists():
        raise ValueError('Require a project-owned mainconf and a new project-owned output directory.')
    source = PROJECT/'tools/nuclear_pid_tof/runtime_selection.py'
    runtime = module('baseline_runtime', source)
    current = runtime.load_current_selection(PROJECT, mainconf)
    output.mkdir(parents=True)
    snapshot = output/'source'
    snapshot.mkdir()
    signatures = {}
    paths = [Path(path) for path in current['protectedFiles']]
    paths += [source, Path(__file__).resolve(), PROJECT/'StMaker/common/FemtoLambdaLegacy.cxx',
              PROJECT/'StMaker/common/FemtoNuclearDedxCalibration.cxx',
              PROJECT/'StMaker/common/FemtoNuclearDedxCalibration.h']
    for path in paths:
        target = snapshot/path.relative_to(PROJECT)
        target.parent.mkdir(parents=True, exist_ok=True)
        expected = digest(path)
        shutil.copy2(path, target)
        if digest(target) != expected or digest(path) != expected:
            raise ValueError('Baseline source changed during snapshot: '+str(path))
        signatures[str(path.relative_to(PROJECT))] = expected
    save(output/'baseline.json', dict(schemaVersion=1, productionEnabled=False,
         createdUtc=datetime.now(timezone.utc).isoformat(), project=str(PROJECT),
         runtimeMainconf=str(mainconf), current=current, sourceSha256=signatures,
         policy='Frozen TPC-only reference, not a claim that working helium calibration is validated. '
                'Authorized subsequent adapter/profile edits do not alter this reference snapshot. '
                'No active TOF cut is inferred from inclusive ALL data.'))
    print('NUCLEAR_CUTOFF_BASELINE_COMPLETE', output, flush=True)


if __name__ == '__main__':
    main()
