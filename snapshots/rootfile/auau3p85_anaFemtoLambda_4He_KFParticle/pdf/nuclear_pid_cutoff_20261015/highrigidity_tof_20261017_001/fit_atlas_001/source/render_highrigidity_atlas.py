#!/usr/bin/env python3
"""高R保存TAG fitの全R分冊表示。native/coarseを明示し新fitなし。"""
import argparse
import concurrent.futures
from pathlib import Path
import shutil
import sys
import yaml
HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import render_native_atlas as common


def flag(value):
    common.require(str(value) in ('True', 'False', '1', '0'), 'Invalid saved flag')
    return str(value) in ('True', '1')


def prepare(settings, output):
    cfg = yaml.safe_load(settings.read_text())
    common.require(cfg['schemaVersion'] == 1 and not output.exists() and
                   output.parent == common.resolve(cfg['outputBase']), 'Choose a new high-R atlas output')
    param = common.resolve(cfg['parameterTable'])
    parameters = {(r['species'], r['label'], round(float(r['r']), 9)): r for r in common.rows(param)}
    done = common.read(param.parent / 'completed.json')
    common.require(done['executionComplete'] and common.sha(param) == done['products'][param.name],
                   'Final parameter table changed/incomplete')
    model_inputs = {str(param.relative_to(common.PROJECT)): {'bytes': param.stat().st_size, 'sha256': common.sha(param)}}
    pending, coverage = [], []
    tol = cfg['axisTolerance']
    for species, target in cfg['targets'].items():
        models, reference, quality, contexts = [], None, [], []
        expected = round((target['rigidityMax'] - target['rigidityMin']) / target['fitWidth'])
        for response in cfg['responseModels']:
            for background in cfg['backgroundModels']:
                folder = common.resolve(cfg['outputBase']) / ('%s_%s_%s_%s_001' %
                         (target['prefix'], species, response, background))
                done = common.read(folder / 'completed.json')
                common.require(done['executionComplete'] and done['protectedSourcesVerified'] and
                               not done.get('selfTestOnly'), 'Incomplete saved high-R fit')
                for name in ('cutoff_fits.root', 'cutoff_fits.json', 'context.json', 'native_summary.tsv'):
                    path = folder / name
                    common.require(common.sha(path) == done['products'][name], 'Changed saved high-R artifact')
                    model_inputs[str(path.relative_to(common.PROJECT))] = {'bytes': path.stat().st_size, 'sha256': common.sha(path)}
                path = folder / 'completed.json'
                model_inputs[str(path.relative_to(common.PROJECT))] = {'bytes': path.stat().st_size, 'sha256': common.sha(path)}
                context, fit = common.read(folder / 'context.json'), common.read(folder / 'cutoff_fits.json')
                common.require(context['targetSpecies'] == species and context['sampleHistogram'] == fit['sampleHistogram'] == target['sample'] and
                    context['mass2Window'] == fit['mass2Window'] and context['exactProductionTofWindow'] == target['exactProductionTofWindow'],
                    'Contradictory TAG/window provenance')
                rows = sorted(fit['slices'], key=lambda r: r['r'])
                metrics = {int(r['rBin']): r for r in common.rows(folder / 'native_summary.tsv')}
                common.require(len(rows) == expected and len(metrics) == expected, 'Incomplete high-R coverage')
                label = ('Gaussian' if response == 'gaussian' else 'Log-normal') + ' + ' + background
                models.append((label, str((folder / 'cutoff_fits.root').relative_to(common.PROJECT))))
                contexts.append((response + '_' + background, folder / 'context.json'))
                for i, row in enumerate(rows):
                    common.require(abs(row['lo'] - (target['rigidityMin'] + i * target['fitWidth'])) < tol and
                        abs(row['hi'] - (target['rigidityMin'] + (i + 1) * target['fitWidth'])) < tol, 'Saved R coverage has gaps/overlap')
                    if reference is not None:
                        common.require(row['index'] == reference[i]['index'] and row['r'] == reference[i]['r'], 'Model index/R mismatch')
                    response_quality = parameters[(species, response + '_' + background + '_nominal', round(row['r'], 9))]
                    metric = metrics[row['rBin']]
                    common.require(metric['sample'] == target['sample'] and abs(float(metric['rLo']) - row['lo']) < tol and
                        abs(float(response_quality['fitRigidityWidth']) - target['fitWidth']) < tol, 'Quality sample/R-width mismatch')
                    item = [species, label, row['index'], row['lo'], row['hi'], row['r'], target['fitWidth'], target['sample'],
                            *context['mass2Window'], row['numericalValid'], response_quality['fitResponseQuality'], flag(metric['quality']),
                            context['exactProductionTofWindow'], target['prefix'] == 'coarse_native', row['fitFailure'], metric['reasons'],
                            'included_regardless_of_support']
                    quality.append(item)
                    coverage.append(item)
                if reference is None:
                    reference = rows
        common.require(len(models) == cfg['modelsPerSpecies'], 'Exactly four saved models required')
        pending.append((species, target, models, reference, quality, contexts))
    output.mkdir()
    (output / 'source').mkdir()
    (output / 'inputs').mkdir()
    for source in (settings, Path(__file__).resolve(), common.resolve(cfg['renderer']), common.resolve(cfg['rendererDependency']),
                   common.resolve(cfg['commonDriver']), common.resolve(cfg['runtimeWrapper']), common.resolve(cfg['runtimeMainconf'])):
        shutil.copy2(str(source), str(output / 'source' / source.name))
    shutil.copy2(str(param), str(output / 'inputs/native_parameters.tsv'))
    header = ('species', 'modelLabel', 'index', 'rLo', 'rHi', 'r', 'fitWidth', 'sample', 'mass2Min', 'mass2Max',
              'numericalValid', 'fitResponseQuality', 'quality', 'exactProductionTofWindow', 'coarse', 'fitFailure',
              'savedQualityReasons', 'displayPolicy')
    books = []
    for species, target, models, reference, quality, contexts in pending:
        common.table(output / ('models_' + species + '.tsv'), ('label', 'path'), models)
        common.table(output / ('quality_' + species + '.tsv'), header, quality)
        for label, context in contexts:
            shutil.copy2(str(context), str(output / 'inputs' / (species + '_' + label + '_context.json')))
        for start in range(0, len(reference), round(cfg['bookRigidityWidth'] / target['fitWidth'])):
            selected = reference[start:start + round(cfg['bookRigidityWidth'] / target['fitWidth'])]
            lo, hi = selected[0]['lo'], selected[-1]['hi']
            identifier = '%s_R%03d_%03d' % (species, round(lo * 100), round(hi * 100))
            common.table(output / (identifier + '.tsv'), ('index', 'rigidity'), [(r['index'], r['r']) for r in selected])
            books.append(dict(id=identifier, species=species, rLo=lo, rHi=hi, fitWidth=target['fitWidth'], sample=target['sample'],
                              observedBins=len(selected), pages=len(selected) * cfg['pagesPerBin']))
    with (output / 'settings.tsv').open('x') as stream:
        for key, value in cfg['plots'].items():
            stream.write(str(key) + '\t' + str(value) + '\n')
    common.table(output / 'coverage.tsv', header, coverage)
    common.save(output / 'plan.json', dict(schemaVersion=1, createdJST=common.stamp(), presentationOnly=True, refitted=False,
        productionEnabled=False, config=cfg, books=books, expectedPages=sum(r['pages'] for r in books), coverageRows=len(coverage),
        modelInputs=model_inputs, failedFitsOmitted=False,
        warning='He4 saved narrow D window differs from pre-study baseline 3<m2<4; provenance does not encode later adoption.',
        definitions=['d/t/He3 use native 0.01-R spectra; He4 is separate 0.05-R directly fitted projected data.',
                     'Response is the saved response screen, runQ is per-model saved quality, not full uncertainty.']))
    print('HIGH_R_ATLAS_PREPARED', output, len(books), sum(r['pages'] for r in books), flush=True)


def run(output, book, workers):
    plan = common.read(output / 'plan.json')
    common.require(1 <= workers <= plan['config']['maxRenderWorkers'], 'Too many render workers')
    common.require(common.sha(Path(__file__).resolve()) == common.sha(output / 'source' / Path(__file__).name),
                   'Driver differs from prepared snapshot')
    common.require(common.sha(Path(common.__file__).resolve()) == common.sha(output / 'source' / Path(common.__file__).name),
                   'Shared rendering helper differs from prepared snapshot')
    selected = [r for r in plan['books'] if book == 'all' or r['id'] == book]
    common.require(bool(selected), 'Unknown explicit book')
    if book == 'all':
        common.save(output / 'parallel_launch.json', dict(workers=workers, argv=sys.argv, createdJST=common.stamp(), presentationOnly=True))
    with concurrent.futures.ThreadPoolExecutor(max_workers=workers) as executor:
        list(executor.map(lambda row: common.run(output, row['id']), selected))
    if book == 'all':
        common.check(output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'check'))
    parser.add_argument('--settings', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--book', default='all')
    parser.add_argument('--workers', type=int, default=1)
    args = parser.parse_args()
    output = common.resolve(args.output_dir)
    if args.action == 'prepare':
        common.require(args.settings is not None, 'Explicit display settings required')
        prepare(common.resolve(args.settings), output)
    elif args.action == 'run':
        run(output, args.book, args.workers)
    else:
        common.check(output)


if __name__ == '__main__':
    main()
