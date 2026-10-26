#!/usr/bin/env python3
"""保存fitだけから旧holeと隣接Rのwide/local・linear/log診断図を新規保存。"""
import argparse
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess
import yaml

PROJECT = Path(__file__).resolve().parents[3]


def resolve(value):
    path = Path(value)
    return path.resolve() if path.is_absolute() else (PROJECT / path).resolve()


def require(ok, message):
    if not ok:
        raise ValueError(message)


def sha(path):
    digest = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            digest.update(chunk)
    return digest.hexdigest()


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def table(path, header, data):
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        if header:
            writer.writerow(header)
        writer.writerows(data)


def rows(path):
    with path.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def flag(value):
    require(str(value).lower() in ('true', 'false', '0', '1'), 'Invalid audit flag')
    return str(value).lower() in ('true', '1')


def quoted(path):
    value = str(path)
    require('"' not in value and '\\' not in value, 'Unsafe ROOT path quoting')
    return '"' + value + '"'


def display_number(value):
    if value is None or str(value).strip().lower() in ('', 'na', 'none', 'null'):
        return 'NA'
    number = float(value)
    return '%.4g' % number if math.isfinite(number) else 'NA'


def focus_extent(item, species, cfg):
    intervals = []
    for name in cfg['focusSpecies'][species]:
        component = item['components'][name]
        mean, width = component['initialMean'], component['initialSigma']
        require(math.isfinite(mean) and math.isfinite(width) and width > 0,
                'Invalid saved initial focus shape')
        intervals.append((mean-cfg['focusSigma']*width, mean+cfg['focusSigma']*width))
        if item['numericalValid']:
            mean, width = component['mean'], component['sigma']
            require(math.isfinite(mean) and math.isfinite(width) and width > 0,
                    'Invalid saved fitted focus shape')
            intervals.append((mean-cfg['focusSigma']*width, mean+cfg['focusSigma']*width))
    low, high = min(x[0] for x in intervals), max(x[1] for x in intervals)
    if high-low < cfg['focusMinimumWidth']:
        midpoint = .5*(low+high)
        low, high = midpoint-.5*cfg['focusMinimumWidth'], midpoint+.5*cfg['focusMinimumWidth']
    return max(cfg['plots']['dedxMin'], low), min(cfg['plots']['dedxMax'], high)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings', required=True)
    parser.add_argument('--audit-table', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    settings, audit_path, output = resolve(args.settings), resolve(args.audit_table), resolve(args.output_dir)
    cfg = yaml.safe_load(settings.read_text())
    require(cfg['schemaVersion'] == 1, 'Unsupported renderer schema')
    require(not output.exists() and output.parent == resolve(cfg['outputBase']),
            'Use a new output sibling under outputBase')
    audit = {}
    for row in rows(audit_path):
        if row['stage'] != 'low_ALL' or row['variation'] != 'nominal':
            continue
        key = (row['target'], row['response'], row['background'], int(row['rowIndex']))
        require(key not in audit, 'Duplicate nominal audit identity')
        audit[key] = row
    require(bool(audit), 'Missing nominal low-R audit')
    immutable = {str(audit_path): sha(audit_path)}
    books = []
    output.mkdir()
    source = output / 'source'
    source.mkdir()
    for path in (settings, Path(__file__).resolve(), resolve(cfg['renderer']),
                 resolve(cfg['rendererDependency']), resolve(cfg['runtimeWrapper']),
                 resolve(cfg['runtimeMainconf'])):
        shutil.copy2(str(path), str(source / path.name))
    table(output/'settings.tsv', (), cfg['plots'].items())
    for species, requested in cfg['representativeRigidity'].items():
        model_table = resolve(cfg['modelTablesDirectory']) / ('models_'+species+'.tsv')
        models = rows(model_table)
        require(len(models) == cfg['modelsPerSpecies'], 'Expected all four nominal models')
        immutable[str(model_table)] = sha(model_table)
        model_rows, selected_rows, reference = [], [], None
        for model in models:
            root = resolve(model['path'])
            fit_path, complete = root.parent/'cutoff_fits.json', root.parent/'completed.json'
            saved = read(complete)
            require(saved['executionComplete'], 'Incomplete saved fit')
            for path in (root, fit_path):
                require(sha(path) == saved['products'][path.name], 'Saved fit product changed')
            fit = read(fit_path)
            require(fit['targetSpecies'] == species and fit['sampleHistogram'] == 'ALL',
                    'Wrong saved target/sample')
            immutable.update({str(p): sha(p) for p in (root, fit_path, complete)})
            model_rows.append([model['label'], str(root)])
            selected = []
            for rigidity in requested:
                item = min(fit['slices'], key=lambda x: abs(x['r']-rigidity))
                require(abs(item['r']-rigidity) < cfg['axisTolerance'],
                        'Requested R is not a native center')
                selected.append([item['index'], item['r']])
                key = (species, fit['responseModel'], fit['backgroundModel'], item['index'])
                require(key in audit, 'Missing nominal quality row')
                a = audit[key]
                low, high = focus_extent(item, species, cfg)
                reason = '; '.join(json.loads(a['reasons'])).replace('_', ' ')
                selected_rows.append([model['label'], item['index'], species,
                    display_number(a['chi2']), display_number(a['ndf']), display_number(a['chi2Ndf']),
                    int(flag(a['numericalValid'])), int(flag(a['supportedWithoutChi2'])), reason,
                    item['fitLo'], item['fitHi'], low, high])
            require(reference is None or reference == selected, 'Saved model/R indices disagree')
            reference = selected
        table(output/('models_'+species+'.tsv'), ('label', 'path'), model_rows)
        table(output/('slices_'+species+'.tsv'), ('index', 'rigidity'), reference)
        table(output/('quality_'+species+'.tsv'), ('label', 'index', 'target', 'chi2', 'ndf',
            'chi2Ndf', 'numerical', 'supported', 'reason', 'fitLo', 'fitHi', 'focusLo', 'focusHi'),
            selected_rows)
        books.append(dict(species=species, pages=len(reference)*cfg['pagesPerRigidity'], rigidity=requested))
    save(output/'plan.json', dict(config=cfg, books=books, immutableInputs=immutable,
        refitted=False, normalized=False, productionEnabled=False, qualityThresholdAdopted=False))
    results = []
    for book in books:
        species = book['species']
        pdf = output/('fit_diagnostics_'+species+'.pdf')
        invocation = str(source/Path(cfg['renderer']).name)+'+('+','.join(quoted(p) for p in (
            output/('models_'+species+'.tsv'), output/('slices_'+species+'.tsv'),
            output/'settings.tsv', output/('quality_'+species+'.tsv'), pdf))+')'
        command = ['bash', str(resolve(cfg['runtimeWrapper'])), str(resolve(cfg['runtimeMainconf'])),
                   'root4star', '-l', '-b', '-q', invocation]
        save(output/('command_'+species+'.json'), dict(cwd=str(PROJECT), argv=command))
        with (output/('render_'+species+'.log')).open('x') as log:
            subprocess.run(command, cwd=str(PROJECT), stdout=log, stderr=subprocess.STDOUT,
                           check=True, timeout=cfg['executionTimeoutSeconds'])
        log_text = (output/('render_'+species+'.log')).read_text()
        require('QUALITY_DIAGNOSTICS_COMPLETE pages='+str(book['pages']) in log_text,
                'Missing renderer completion')
        require(log_text.count('SPECTRUM_VISIBLE_AXES ') == book['pages']*4,
                'Missing per-panel axis audit')
        check = subprocess.run(['gs', '-dBATCH', '-dNOPAUSE', '-sDEVICE=nullpage', str(pdf)],
            stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=True)
        with (output/('readback_'+species+'.log')).open('x') as log:
            log.write(check.stdout)
        match = re.search(r'Processing pages 1 through ([0-9]+)\.', check.stdout)
        require(match and int(match.group(1)) == book['pages'] and 'Error:' not in check.stdout,
                'Closed-PDF readback/page check failed')
        results.append(dict(species=species, pages=book['pages'], sha256=sha(pdf),
                            bytes=pdf.stat().st_size, closedReadback=True))
    for path, digest in immutable.items():
        require(sha(Path(path)) == digest, 'Immutable fit/count input changed')
    save(output/'completed.json', dict(status='completed',
        completedJST=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
        refitted=False, normalized=False, productionEnabled=False, qualityThresholdAdopted=False,
        failedFitsOmitted=False, results=results, totalPages=sum(x['pages'] for x in results),
        immutableInputsUnchanged=True, outputHashes={str(p.relative_to(output)): sha(p)
                                                    for p in output.rglob('*') if p.is_file()}))
    print('QUALITY_DIAGNOSTICS_COMPLETE', output, sum(x['pages'] for x in results), flush=True)


if __name__ == '__main__':
    main()
