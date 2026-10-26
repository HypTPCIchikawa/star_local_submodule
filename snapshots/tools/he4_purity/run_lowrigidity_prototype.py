#!/usr/bin/env python3
"""Refit selected cached rigidity slices before running the complete study."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--sample', choices=['3p85', '13p5'], required=True)
    p.add_argument('--mainconf', type=Path, required=True)
    p.add_argument('--settings', type=Path, default=HERE/'purity_lowrigidity_settings.yaml')
    p.add_argument('--output-dir', type=Path, required=True)
    p.add_argument('--cache', type=Path, help='Read an explicit cached he4_fits.root; default is the historical cache.')
    p.add_argument('--indices', default='', help='Comma-separated 0-based rigidity slice indices; omitted preserves historical points.')
    args = p.parse_args()
    cfg = yaml.safe_load(args.settings.read_text())
    sample = cfg['samples'][args.sample]
    mainconf = args.mainconf.resolve(strict=True)
    if mainconf != (PROJECT/sample['mainconf']).resolve(strict=True):
        raise ValueError('Explicit mainconf differs from the selected preset.')
    cache = args.cache if args.cache is not None else (PROJECT/sample['input']).parent/'pdf/he4_purity_20261013/he4_fits.root'
    cache = cache.resolve(strict=True)
    indices = [] if not args.indices else [int(item) for item in args.indices.split(',')]
    if any(index < 0 or index >= 210 for index in indices) or len(indices) != len(set(indices)):
        raise ValueError('Indices must be unique integers in [0, 209].')
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=False)
    manifest = out/'settings.tsv'
    manifest.write_text(''.join('%s\t%s\n' % (k, int(v) if isinstance(v, bool) else v)
                               for k, v in cfg.items() if k not in ('models', 'samples')))
    sources = {str(p):hashlib.sha256(p.read_bytes()).hexdigest()
               for p in (HERE/'FitHe4Purity.C', HERE/'CheckLowRigidity.C', args.settings.resolve())}
    (out/'context.json').write_text(json.dumps(dict(cache=str(cache), cache_sha256=hashlib.sha256(cache.read_bytes()).hexdigest(),
        sample=args.sample, settings=cfg, slice_indices=indices if indices else [10,15,20,25,30,35,45,80,120], sources=sources), indent=2)+'\n')
    call = str(HERE/'CheckLowRigidity.C')+'++('+','.join(json.dumps(str(x)) for x in (cache, manifest, out, args.indices))+')'
    with (out/'prototype.log').open('w') as handle:
        result = subprocess.run(['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
                                 'root4star', '-l', '-b', '-q', call], cwd=str(PROJECT),
                                stdout=handle, stderr=subprocess.STDOUT, timeout=600)
    if result.returncode or 'LOW_R_PROTOTYPE_COMPLETE' not in (out/'prototype.log').read_text():
        raise RuntimeError('Prototype failed; see '+str(out/'prototype.log'))
    print('COMPLETE', args.sample, out, flush=True)


if __name__ == '__main__':
    main()
