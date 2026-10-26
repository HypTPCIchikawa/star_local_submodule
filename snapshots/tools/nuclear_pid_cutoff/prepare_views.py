#!/usr/bin/env python3
"""Prepare explicit model/slice tables for the ROOT 5 spectrum renderer."""
import argparse
from pathlib import Path
import yaml
from run_cutoff import read, resolve, save, tsv


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings', required=True, type=Path)
    parser.add_argument('--runs', required=True, nargs='+', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    cfg = yaml.safe_load(resolve(args.settings).read_text())
    output = args.output_dir.resolve()
    if output.exists() or output.parent != resolve(cfg['outputBase']):
        raise ValueError('Choose a new views sibling below outputBase.')
    models, reference = [], None
    for folder in map(resolve, args.runs):
        if not read(folder/'completed.json')['executionComplete']:
            raise ValueError('Incomplete spectra source.')
        c = read(folder/'context.json')
        models.append(['_'.join(c[k] for k in ('response','background','variation')),
                       str(folder/'cutoff_fits.root')])
        if reference is None:
            reference = read(folder/'cutoff_fits.json')['slices']
    output.mkdir()
    tsv(output/'models.tsv', 'label path', models)
    slices = []
    for r in cfg['plots']['representativeRigidity']:
        row = min(reference, key=lambda x: abs(x['r']-r))
        if abs(row['r']-r) > cfg['numeric']['axisTolerance']:
            raise ValueError('Requested representative R is not a saved native-bin center.')
        slices.append([row['index'], row['r']])
    tsv(output/'slices.tsv', 'index rigidity', slices)
    tsv(output/'settings.tsv', '', [[k,cfg['plots'][k]] for k in
        ('dedxMin','dedxMax','displayDedxRebinFactor','logMinimum')])
    save(output/'context.json', dict(productionEnabled=False, inputModels=models,
                                   slices=slices, plots=cfg['plots']))
    print('NUCLEAR_CUTOFF_VIEWS_PREPARED', output, flush=True)


if __name__ == '__main__':
    main()
