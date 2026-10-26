#!/usr/bin/env python3
"""Execute declared independent response/background comparisons, with limited workers."""
import argparse
from concurrent.futures import ThreadPoolExecutor
import subprocess
import sys
from pathlib import Path
import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('settings', 'baseline', 'mainconf'):
        parser.add_argument('--'+name, required=True, type=Path)
    parser.add_argument('--variation', action='append', required=True)
    parser.add_argument('--target', choices=('d','t','he3','he4'), action='append')
    parser.add_argument('--response', choices=('gaussian','lognormal'), action='append')
    parser.add_argument('--background', choices=('linear','exponential'), action='append')
    parser.add_argument('--suffix', required=True)
    parser.add_argument('--workers', type=int, required=True)
    args = parser.parse_args()
    if args.workers < 1 or args.workers > 4 or not args.suffix.isalnum():
        raise ValueError('Require 1–4 local workers and an alphanumeric new run suffix.')
    cfg = yaml.safe_load(args.settings.read_text())
    responses = args.response or cfg['responseModels']
    backgrounds = args.background or cfg['backgroundModels']
    if set(responses)-set(cfg['responseModels']) or set(backgrounds)-set(cfg['backgroundModels']):
        raise ValueError('Requested execution subset is not declared by the study settings.')
    base = PROJECT/cfg['outputBase']
    jobs = []
    for variation in args.variation:
        if variation not in cfg['shapeVariations']:
            raise ValueError('Undeclared sensitivity variation.')
        for target in args.target or [None]:
            for response in responses:
                for background in backgrounds:
                    out = base/('all_'+(target+'_' if target else '')+response+'_'+background+'_'+variation+'_'+args.suffix)
                    if out.exists():
                        raise ValueError('Never overwrite an existing run: '+str(out))
                    command = [sys.executable, '-B', str(HERE/'run_cutoff.py'),
                        '--settings', str(args.settings.resolve()), '--baseline', str(args.baseline.resolve()),
                        '--mainconf', str(args.mainconf.resolve()), '--response', response,
                        '--background', background, '--variation', variation, '--output-dir', str(out)]
                    if target:
                        command += ['--target', target]
                    jobs.append((out, command))
    def execute(job):
        out, command = job
        print('SUITE_START', out.name, flush=True)
        result = subprocess.run(command, cwd=PROJECT)
        print('SUITE_END', out.name, 'exit='+str(result.returncode), flush=True)
        return result.returncode
    with ThreadPoolExecutor(max_workers=args.workers) as pool:
        status = list(pool.map(execute, jobs))
    if any(status):
        raise SystemExit('Some runs failed; preserve each log/output and choose a new suffix after diagnosis.')
    print('NUCLEAR_CUTOFF_SUITE_COMPLETE', len(jobs), flush=True)


if __name__ == '__main__':
    main()
