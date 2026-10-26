#!/usr/bin/env python3
"""Render saved fixed-polynomial peak diagnostics; never refit or certify purity."""
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
SPECIES = ('d', 't', 'he3', 'he4')
REGIONS = ('helium', 'tControl')
METRIC_HEADER = ('label', 'index', 'response', 'numerical', 'fitLo', 'fitHi', 'heLo', 'heHi',
                 'tLo', 'tHi', 'muT', 'sigmaT', 'muHe3', 'sigmaHe3', 'muHe4', 'sigmaHe4',
                 'jointChi2', 'jointNdf', 'jointRatio', 'heChi2', 'heBins', 'hePull',
                 'tChi2', 'tBins', 'tPull', 'reason')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def owned(value, exists=True):
    path = Path(value)
    path = (path if path.is_absolute() else PROJECT/path).resolve()
    require(PROJECT in path.parents, 'Require an explicit project-owned path')
    if exists:
        require(path.exists(), 'Missing explicit input: '+str(path))
    return path


def sha(path):
    digest = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            digest.update(block)
    return digest.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def table(path, header, rows):
    with Path(path).open('x', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        if header:
            writer.writerow(header)
        writer.writerows(rows)


def quoted(value):
    value = str(value)
    require(not any(character in value for character in ('"', '\\', '\n')), 'Unsafe ROOT quoting')
    return '"'+value+'"'


def finite(value):
    return type(value) in (int, float) and math.isfinite(value)


def display(value):
    return '%.5g' % value if finite(value) else 'NA'


def fixed_components(fit, row):
    """The physical nuclear centers/widths must be the original PID curves, not free BG."""
    require(fit.get('shapeSource') == 'legacy_log10_polynomial',
            'Require explicitly declared original log10 polynomial response shapes')
    require(fit.get('backgroundPeakShapesFree') is False,
            'Background nuclear peak centers/widths may not be free parameters')
    components = row['components']
    require(set(components) == set(SPECIES), 'Require all four named nuclear components')
    for name in SPECIES:
        component = components[name]
        require(component.get('freeMean') is False and component.get('freeSigma') is False,
                'Physical nuclear peak parameters must stay fixed: '+name)
        for key, initial in (('mean', 'initialMean'), ('sigma', 'initialSigma')):
            require(finite(component.get(key)) and finite(component.get(initial))
                    and component[key] == component[initial],
                    'Saved nuclear shape differs from its original PID value: '+name+'/'+key)
        require(component['sigma'] > 0, 'Nonpositive fixed PID width')
    return components


def metric_row(job, fit, row, plots):
    components = fixed_components(fit, row)
    require(type(row['numericalValid']) is bool, 'Require saved numerical bool, never infer convergence')
    joint, helium, triton = (row['regions'][name] for name in ('joint', 'helium', 'tControl'))
    require(joint['lo'] == row['fitLo'] and joint['hi'] == row['fitHi'], 'Joint region differs from actual fit window')
    for region in (joint, helium, triton):
        require(finite(region['lo']) and finite(region['hi']) and region['hi'] > region['lo'],
                'Invalid saved diagnostic region')
        require(joint['lo'] <= region['lo'] < region['hi'] <= joint['hi'],
                'Diagnostic ROI exceeds saved joint fit window')
        require(type(region['binCount']) is int and region['binCount'] > 0, 'Require actual regional bin count')
    for region in (helium, triton):
        require(plots['dedxMin'] <= region['lo'] < region['hi'] <= plots['dedxMax'],
                'Focused region exceeds explicitly saved display axis')
        require((region['lo'], region['hi']) != (plots['dedxMin'], plots['dedxMax']),
                'Linear nuclear diagnostic may not use the full dE/dx axis')
    chi2, ndf = row.get('pearson'), row.get('conditionalNdf')
    ratio = chi2/ndf if finite(chi2) and finite(ndf) and ndf > 0 else None
    values = [job['id'], row['index'], job['response'], int(row['numericalValid']),
              row['fitLo'], row['fitHi'], helium['lo'], helium['hi'], triton['lo'], triton['hi']]
    for name in ('t', 'he3', 'he4'):
        values.extend((components[name]['mean'], components[name]['sigma']))
    values.extend((display(chi2), display(ndf), display(ratio), display(helium.get('pearson')),
                   helium['binCount'], display(helium.get('maxAbsPearsonPull')),
                   display(triton.get('pearson')), triton['binCount'],
                   display(triton.get('maxAbsPearsonPull')), str(row.get('fitFailure', '')).replace('\t', ' ').replace('\n', ' ')))
    return values


def verify_map(mapping):
    for path, expected in mapping.items():
        require(sha(path) == expected, 'Protected saved input changed: '+str(path))


def native_selection(saved, planned, tolerance):
    """Match actual native-bin identities; tolerate only requested-center roundoff."""
    require(finite(tolerance) and tolerance > 0 and len(saved) == len(planned) and bool(saved),
            'Invalid native selection tolerance or incomplete saved selection')
    result = []
    for row, expected in zip(saved, planned):
        require(all(row[key] == expected[key] for key in ('index', 'rBin', 'lo', 'hi')),
                'Saved native index/bin/edges differ from prepared plan')
        require(type(row['index']) is int and row['index'] >= 0
                and type(row['rBin']) is int and row['rBin'] > 0
                and finite(row['lo']) and finite(row['hi']) and row['hi'] > row['lo'],
                'Invalid saved native bin identity')
        midpoint = .5*(row['lo']+row['hi'])
        require(finite(row['r']) and row['r'] == midpoint,
                'Saved rigidity is not the actual native-bin midpoint')
        require(finite(expected['r']) and abs(midpoint-expected['r']) <= tolerance,
                'Requested rigidity materially differs from the actual native-bin midpoint')
        result.append((row['index'], midpoint))
    return result


def verify_panels(text, selected, models, expected_pages):
    lines = [line for line in text.splitlines() if line.startswith('PEAK_ANCHORED_PANEL ')]
    expected_panels = len(selected)*len(models)*len(REGIONS)*2
    require(len(lines) == expected_panels == expected_pages*4, 'Incomplete saved peak diagnostic panels')
    rows = [dict(item.split('=', 1) for item in line.split()[1:]) for line in lines]
    require([int(row['serial']) for row in rows] == list(range(expected_panels)), 'Repeated/out-of-order peak panel')
    for index, _ in selected:
        for region in REGIONS:
            for logarithmic in (0, 1):
                matching = [row for row in rows if int(row['index']) == index
                            and row['role'] == region and int(row['log']) == logarithmic]
                require(len(matching) == len(models)
                        and {row['label'] for row in matching} == {model[0] for model in models},
                        'A model/R/ROI/scale combination was omitted')
    require(text.count('SPECTRUM_VISIBLE_AXES ') == expected_panels, 'Missing independent visible-axis audit')
    return rows


def prepare(study, mainconf, output):
    require(not output.exists() and study not in output.parents, 'Use a fresh sibling output, never overwrite a study')
    plan_path, done_path = study/'plan.json', study/'completed.json'
    plan, completed = read(plan_path), read(done_path)
    cfg = plan['config']
    require(plan.get('productionEnabled') is False and completed.get('productionEnabled') is False
            and completed.get('status') == 'completed' and completed.get('immutableInputsUnchanged') is True,
            'Require a fully completed immutable diagnostic study')
    require(mainconf == owned(plan['explicitMainconf']) == owned(cfg['runtimeMainconf']),
            'Explicit rendering mainconf differs from study runtime')
    jobs = sorted(plan['jobs'], key=lambda job: job['id'])
    require(len(jobs) == 4 and {(job['response'], job['background']) for job in jobs}
            == {(response, background) for response in ('gaussian', 'lognormal')
                for background in ('linear', 'exponential')}, 'Require all four declared response/background models')
    plots = cfg['plots']
    require(plots['displayDedxRebinFactor'] >= 1 and plots['displayDedxRebinFactor'] == int(plots['displayDedxRebinFactor'])
            and plots['logMinimum'] > 0 and plots['dedxMax'] > plots['dedxMin'], 'Invalid explicit display settings')
    immutable = {str(plan_path): sha(plan_path), str(done_path): sha(done_path)}
    immutable.update(plan['immutableInputs'])
    models, metrics, selected = [], [], None
    for job in jobs:
        folder = owned(job['folder'])
        job_done_path = folder/'completed.json'
        job_done = read(job_done_path)
        require(job_done.get('executionComplete') is True and job_done.get('rawReadback') is True,
                'Incomplete saved peak-anchored job')
        immutable[str(job_done_path)] = sha(job_done_path)
        immutable.update({str(folder/name): digest for name, digest in job_done['products'].items()})
        fit = read(folder/'cutoff_fits.json')
        require(fit['responseModel'] == job['response'] and fit['backgroundModel'] == job['background']
                and fit['sampleHistogram'] == 'ALL' and fit['productionEnabled'] is False
                and fit.get('responseStudyMode') == 'peak_anchored_fixed_pid',
                'Wrong saved response/background/sample identity')
        current = native_selection(fit['slices'], job['slices'], cfg['numeric']['axisTolerance'])
        require(selected is None or selected == current, 'Compared models do not use the same R slices')
        selected = current
        models.append([job['id'], str(folder/'cutoff_fits.root')])
        metrics.extend(metric_row(job, fit, row, plots) for row in fit['slices'])
    require(completed.get('totalJobs') == len(jobs) and completed.get('totalFits') == len(selected)*len(jobs),
            'Completed study omitted a fit/model')
    wrapper = owned(cfg['runtimeWrapper'])
    renderer = Path(__file__).resolve().with_name('PlotPeakAnchored.C')
    dependency = owned('tools/nuclear_pid_cutoff/PlotCutoffSpectra.C')
    sources = {Path(__file__).resolve(), renderer, dependency, wrapper, mainconf}
    main = yaml.safe_load(mainconf.read_text())
    sources.update(owned(Path('config')/reference) for reference in main.values())
    immutable.update({str(path): sha(path) for path in sources})
    verify_map(immutable)
    output.mkdir(parents=True)
    snapshot = output/'source'
    snapshot.mkdir()
    for source in sources:
        destination = snapshot/source.relative_to(PROJECT) if source.suffix == '.yaml' else snapshot/source.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(source), str(destination))
        require(sha(destination) == immutable[str(source)], 'Renderer snapshot differs from source')
    table(output/'models.tsv', ('label', 'path'), models)
    table(output/'slices.tsv', ('index', 'rigidity'), selected)
    table(output/'settings.tsv', (), plots.items())
    table(output/'metrics.tsv', METRIC_HEADER, metrics)
    pages = len(REGIONS)*((len(models)+1)//2)*len(selected)
    rendered_plan = dict(schemaVersion=1, study=str(study), mainconf=str(mainconf),
                         wrapper=str(wrapper), renderer=str(snapshot/renderer.name),
                         immutableInputs=immutable, selected=selected, models=models, totalPages=pages,
                         plots=plots, executionTimeoutSeconds=cfg['executionTimeoutSeconds'],
                         shapeSource='legacy_log10_polynomial', allNuclearShapesFixed=True,
                         refitted=False, normalized=False, productionEnabled=False, purityCertified=False)
    save(output/'plan.json', rendered_plan)
    return rendered_plan


def render(output, plan):
    pdf = output/'peak_anchored_diagnostics.pdf'
    invocation = plan['renderer']+'+('+','.join(quoted(path) for path in
        (output/'models.tsv', output/'slices.tsv', output/'settings.tsv', output/'metrics.tsv', pdf))+')'
    command = ['bash', plan['wrapper'], plan['mainconf'], 'root4star', '-l', '-b', '-q', invocation]
    save(output/'command.json', dict(cwd=str(PROJECT), argv=command))
    with (output/'render.log').open('x') as stream:
        subprocess.run(command, cwd=str(PROJECT), stdout=stream, stderr=subprocess.STDOUT,
                       check=True, timeout=plan['executionTimeoutSeconds'])
    text = (output/'render.log').read_text()
    require('PEAK_ANCHORED_VIEWS_COMPLETE pages='+str(plan['totalPages']) in text,
            'Missing closed peak PDF completion marker')
    verify_panels(text, plan['selected'], plan['models'], plan['totalPages'])
    check = subprocess.run(['gs', '-dBATCH', '-dNOPAUSE', '-dPDFSTOPONERROR', '-sDEVICE=nullpage', str(pdf)],
                           stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True, check=True)
    with (output/'readback.log').open('x') as stream:
        stream.write(check.stdout)
    match = re.search(r'Processing pages 1 through ([0-9]+)\.', check.stdout)
    require(match and int(match.group(1)) == plan['totalPages'] and 'Error:' not in check.stdout,
            'Closed peak PDF readback/page count failed')
    verify_map(plan['immutableInputs'])
    save(output/'completed.json', dict(status='completed',
         completedJST=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
         pages=plan['totalPages'], bytes=pdf.stat().st_size, sha256=sha(pdf), closedReadback=True,
         immutableInputsUnchanged=True, refitted=False, normalized=False, productionEnabled=False,
         allNuclearShapesFixed=True, shapeSource='legacy_log10_polynomial', purityCertified=False,
         failedFitsOmitted=False,
         outputHashes={str(path.relative_to(output)): sha(path) for path in output.rglob('*') if path.is_file()}))
    print('PEAK_ANCHORED_PDFS_COMPLETE', plan['totalPages'], output, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('study-dir', 'mainconf', 'output-dir'):
        parser.add_argument('--'+name, required=True)
    args = parser.parse_args()
    study, mainconf, output = (owned(value, exists=name != 'output') for name, value in
                               (('study', args.study_dir), ('mainconf', args.mainconf), ('output', args.output_dir)))
    render(output, prepare(study, mainconf, output))


if __name__ == '__main__':
    main()
