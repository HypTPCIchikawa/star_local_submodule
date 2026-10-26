#!/usr/bin/env python3
"""保存済みnative ALL fitを全R分冊表示する。fit・正規化・選択変更なし。"""
import argparse
import csv
import datetime
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time

import yaml


PROJECT = Path(__file__).resolve().parents[2]


def resolve(path):
    p = Path(path)
    return p.resolve() if p.is_absolute() else (PROJECT / p).resolve()


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            digest.update(chunk)
    return digest.hexdigest()


def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def table(path, header, data):
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(data)


def require(ok, message):
    if not ok:
        raise ValueError(message)


def prepare(settings, output):
    cfg = yaml.safe_load(settings.read_text())
    require(cfg['schemaVersion'] == 1, 'Unsupported settings schema')
    require(not output.exists() and output.parent == resolve(cfg['outputBase']),
            'Choose a new atlas directory directly below outputBase')
    width = cfg['nativeRigidityWidth']
    count = round((cfg['rigidityMax'] - cfg['rigidityMin']) / width)
    per_book = round(cfg['bookRigidityWidth'] / width)
    require(count > 0 and per_book > 0 and count % per_book == 0,
            'Book partition must exactly cover the native range')
    tolerance = cfg['axisTolerance']
    quality = {(r['label'], r['species'], round(float(r['r']), 9)): r
               for r in rows(resolve(cfg['parameterTable'])) if r['label'].endswith('_nominal')}
    plans, coverage, inputs = [], [], {}
    pending = []
    for species in cfg['species']:
        model_table = resolve(cfg['modelTablesDirectory']) / ('models_' + species + '.tsv')
        models = rows(model_table)
        require(len(models) == cfg['modelsPerSpecies'], 'Wrong number of saved nominal models')
        reference = None
        for model in models:
            root = resolve(model['path'])
            folder = root.parent
            done = read(folder / 'completed.json')
            require(done['executionComplete'], 'Incomplete saved fit: ' + str(folder))
            fit = read(folder / 'cutoff_fits.json')
            require(fit['targetSpecies'] == species and fit['sampleHistogram'] == 'ALL',
                    'Atlas must contain the requested ALL target fits')
            current = sorted(fit['slices'], key=lambda x: x['r'])
            require(len(current) == count, 'Missing native bins: ' + str(folder))
            for i, item in enumerate(current):
                require(abs(item['lo'] - (cfg['rigidityMin'] + i * width)) < tolerance
                        and abs(item['hi'] - (cfg['rigidityMin'] + (i + 1) * width)) < tolerance,
                        'Wrong native-bin coverage')
                if reference is not None:
                    require(item['index'] == reference[i]['index'] and item['r'] == reference[i]['r'],
                            'Saved model index/R mismatch')
                label = fit['responseModel'] + '_' + fit['backgroundModel'] + '_nominal'
                audit = quality[(label, species, round(item['r'], 9))]
                coverage.append([species, label, item['index'], item['lo'], item['hi'], item['r'],
                                 item['numericalValid'], item['fitFailure'], item['atBoundary'],
                                 audit['fitResponseQuality'], audit['targetFreeSupported'],
                                 'included_regardless_of_support'])
            if reference is None:
                reference = current
            for path in (root, folder / 'cutoff_fits.json', folder / 'completed.json'):
                inputs[str(path.relative_to(PROJECT))] = {'bytes': path.stat().st_size, 'sha256': sha(path)}
        pending.append((species, model_table, reference))
    output.mkdir()
    (output / 'source').mkdir()
    (output / 'inputs').mkdir()
    for source in (settings, Path(__file__).resolve(), resolve(cfg['renderer']), resolve(cfg['rendererDependency']), resolve(cfg['runtimeWrapper']),
                   resolve(cfg['runtimeMainconf'])):
        shutil.copy2(str(source), str(output / 'source' / source.name))
    shutil.copy2(str(resolve(cfg['parameterTable'])), str(output / 'inputs/native_parameters.tsv'))
    # The existing ROOT renderer requires a headerless two-column settings table.
    with (output / 'settings.tsv').open('x') as stream:
        for key, value in cfg['plots'].items():
            stream.write(str(key) + '\t' + str(value) + '\n')
    for species, model_table, reference in pending:
        shutil.copy2(str(model_table), str(output / ('models_' + species + '.tsv')))
        table(output / ('quality_' + species + '.tsv'), ('species', 'model', 'index', 'rLo', 'rHi', 'r',
              'numericalValid', 'fitFailure', 'atBoundary', 'fitResponseQuality',
              'targetFreeSupported', 'displayPolicy'), [r for r in coverage if r[0] == species])
        for start in range(0, count, per_book):
            selected = reference[start:start + per_book]
            low, high = selected[0]['lo'], selected[-1]['hi']
            identifier = '%s_R%03d_%03d' % (species, round(low * 100), round(high * 100))
            table(output / (identifier + '.tsv'), ('index', 'rigidity'),
                  [(x['index'], x['r']) for x in selected])
            plans.append(dict(id=identifier, species=species, rLo=low, rHi=high,
                              nativeBins=len(selected), pages=len(selected) * cfg['pagesPerNativeBin']))
    table(output / 'coverage.tsv', ('species', 'model', 'index', 'rLo', 'rHi', 'r',
          'numericalValid', 'fitFailure', 'atBoundary', 'fitResponseQuality',
          'targetFreeSupported', 'displayPolicy'), coverage)
    save(output / 'plan.json', dict(schemaVersion=1, createdJST=stamp(), presentationOnly=True,
         refitted=False, productionEnabled=False, config=cfg, books=plans,
         expectedPages=sum(x['pages'] for x in plans), coverageRows=len(coverage),
         modelInputs=inputs, failedFitsOmitted=False,
         limitations=['Saved unsupported fits are shown as diagnostics, not validated physics.',
                     'Each .01 rigidity interval is drawn separately; display does not refit or normalize.',
                     'Each native bin has 4 pages: 2 model pairs times full/fit-range views; linear and log panels.']))
    print('NATIVE_ATLAS_PREPARED', output, len(plans), sum(x['pages'] for x in plans), flush=True)


def run(output, book_id):
    plan = read(output / 'plan.json')
    cfg = plan['config']
    selected = [x for x in plan['books'] if book_id == 'all' or x['id'] == book_id]
    require(bool(selected), 'Unknown explicit book ID')
    for book in selected:
        identifier = book['id']
        pdf = output / (identifier + '.pdf')
        done = output / (identifier + '.completed.json')
        if done.exists():
            require(sha(pdf) == read(done)['sha256'], 'Previously completed PDF changed')
            continue
        log = output / (identifier + '.log')
        require(not pdf.exists() and not log.exists(), 'Refusing to overwrite unfinished/previous rendering')
        def quoted(path):
            value = str(path)
            require('"' not in value and '\\' not in value, 'Unsupported path quoting')
            return '"' + value + '"'
        invocation = str(output / 'source' / Path(cfg['renderer']).name) + '+(' + ','.join(
            quoted(x) for x in (output / ('models_' + book['species'] + '.tsv'),
                               output / (identifier + '.tsv'), output / 'settings.tsv',
                               output / ('quality_' + book['species'] + '.tsv'), pdf)) + ')'
        wrapper = resolve(cfg['runtimeWrapper'])
        require(sha(wrapper) == sha(output / 'source' / wrapper.name), 'Runtime wrapper changed since preparation')
        command = ['bash', str(wrapper),
                   str(resolve(cfg['runtimeMainconf'])), 'root4star', '-l', '-b', '-q', invocation]
        save(output / (identifier + '.command.json'), dict(cwd=str(PROJECT), argv=command))
        started = time.monotonic()
        with log.open('x') as stream:
            subprocess.run(command, cwd=str(PROJECT), stdout=stream, stderr=subprocess.STDOUT,
                           check=True, timeout=cfg['executionTimeoutSeconds'])
        expected = 'NUCLEAR_CUTOFF_ATLAS_COMPLETE pages=' + str(book['pages'])
        require(expected in log.read_text(), 'ROOT did not confirm the expected closed book')
        result = subprocess.run(['gs', '-dBATCH', '-dNOPAUSE', '-sDEVICE=nullpage', str(pdf)],
                                stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                                universal_newlines=True, check=True)
        with (output / (identifier + '.readback.log')).open('x') as stream:
            stream.write(result.stdout)
        match = re.search(r'Processing pages 1 through ([0-9]+)\.', result.stdout)
        require(match and int(match.group(1)) == book['pages'] and 'Error:' not in result.stdout,
                'Closed PDF readback/page count failed')
        info = dict(id=identifier, pages=book['pages'], bytes=pdf.stat().st_size,
                    elapsedSeconds=time.monotonic() - started, sha256=sha(pdf),
                    completedJST=stamp(), closedPdfReadback=True, refitted=False)
        save(done, info)
        print('NATIVE_ATLAS_BOOK_COMPLETE', json.dumps(info, sort_keys=True), flush=True)


def check(output):
    plan = read(output / 'plan.json')
    books = []
    for item in plan['books']:
        done = read(output / (item['id'] + '.completed.json'))
        pdf = output / (item['id'] + '.pdf')
        require(done['pages'] == item['pages'] and done['sha256'] == sha(pdf),
                'Book coverage or closed-file digest changed')
        books.append(done)
    for relative, metadata in plan['modelInputs'].items():
        require(sha(PROJECT / relative) == metadata['sha256'], 'Immutable saved fit input changed')
    sources = {str(p.relative_to(output)): sha(p) for p in sorted(output.rglob('*'))
               if p.is_file() and p.suffix in ('.py', '.C', '.sh', '.yaml', '.tsv', '.log')}
    save(output / 'completed.json', dict(status='completed', completedJST=stamp(),
         presentationOnly=True, refitted=False, productionEnabled=False, books=books,
         totalPages=sum(x['pages'] for x in books), totalBytes=sum(x['bytes'] for x in books),
         nativeCoverageRows=plan['coverageRows'], failedFitsOmitted=False,
         immutableInputChecks=True, sourceAndTableHashes=sources))
    print('NATIVE_ATLAS_COMPLETE', len(books), sum(x['pages'] for x in books), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'check'))
    parser.add_argument('--settings', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--book', default='all')
    args = parser.parse_args()
    output = resolve(args.output_dir)
    if args.action == 'prepare':
        require(args.settings is not None, 'prepare requires explicit settings YAML')
        prepare(resolve(args.settings), output)
    elif args.action == 'run':
        run(output, args.book)
    else:
        check(output)


if __name__ == '__main__':
    main()
