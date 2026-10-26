#!/usr/bin/env python3
"""Render completed pilot counts without refitting or certifying PID/purity."""
import argparse
import csv
import datetime
import hashlib
import importlib.util
import itertools
import json
import math
from pathlib import Path
import shutil
import subprocess

import yaml


PROJECT = Path(__file__).resolve().parents[3]
QUALITY = PROJECT/'tools/nuclear_pid_cutoff/quality_20261019'
PROFILE_ORDER = ('gaussian_control', 'symmetric_tail', 'positive_tail', 'negative_tail')
FREE_PROFILE_ORDER = ('gaussian_control', 'free_symmetric_seed')
FOCUS = {'d': ['d', 't'], 't': ['t'], 'he4': ['he3', 'he4']}
PENDING = 'Diagnostic pilot only: identity/calibration/purity certification pending'
QUALITY_HEADER = ('label', 'index', 'target', 'chi2', 'ndf', 'chi2Ndf', 'numerical',
                  'supported', 'supportKnown', 'reason', 'fitLo', 'fitHi', 'focusLo', 'focusHi')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def sha(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            value.update(block)
    return value.hexdigest()


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


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, str(path))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def owned(value):
    path = Path(value)
    path = (path if path.is_absolute() else PROJECT/path).resolve()
    require(PROJECT in path.parents, 'Require a project-owned path')
    return path


def mode_contract(cfg, mode):
    """Require an explicit free-model declaration; never infer the model from counts."""
    require(mode in ('fixed', 'free_coretail'), 'Unknown explicit renderer mode')
    if mode == 'free_coretail':
        require(cfg.get('responseStudyMode') == 'free_coretail'
                and cfg.get('responseModel') == 'gaussian_coretail_free_components'
                and Path(cfg.get('adapter', '')).name == 'FitFreeCoreTailPilot.C'
                and Path(cfg.get('engine', '')).name == 'FreeGaussianCoreTail.h',
                'Free mode requires the declared free adapter, engine and response model')
        return dict(profiles=FREE_PROFILE_ORDER, jobs=24, fits=32,
                    responseModel='gaussian_coretail_free_components',
                    parameterMeaning='Seed values only; Gaussian centers, widths and yields are fitted')
    require(cfg.get('responseStudyMode', 'fixed') == 'fixed',
            'Free-model plan cannot be interpreted as a fixed-profile comparison')
    return dict(profiles=PROFILE_ORDER, jobs=48, fits=64,
                responseModel='gaussian_coretail_fixed_profile',
                parameterMeaning='Fixed tail profiles; only declared core/nuisance parameters are fitted')


def grouped_jobs(plan, mode='fixed'):
    """Require the full declared 3-target x 2-background x 2-policy comparison."""
    cfg = plan['config']
    contract = mode_contract(cfg, mode)
    profiles = contract['profiles']
    require(set(cfg['representativeRigidity']) == set(FOCUS), 'Unexpected pilot targets')
    require(set(cfg['backgrounds']) == {'linear', 'exponential'}, 'Compare both backgrounds')
    require(set(cfg['freePolicies']) == {'target_only', 'local_neighbors'}, 'Compare both free policies')
    require(set(cfg['profiles']) == set(profiles), 'Compare every declared model/profile')
    expected = set(itertools.product(FOCUS, cfg['backgrounds'], cfg['freePolicies']))
    groups, identifiers = {}, set()
    for job in plan['jobs']:
        require(job['id'] not in identifiers, 'Duplicate pilot job id')
        identifiers.add(job['id'])
        key = job['target'], job['background'], job['freePolicy']
        require(key in expected and job['profile'] in profiles, 'Unknown pilot group/profile')
        members = groups.setdefault(key, {})
        require(job['profile'] not in members, 'Repeated profile within a pilot group')
        members[job['profile']] = job
    require(set(groups) == expected and len(plan['jobs']) == contract['jobs'],
            'Incomplete planned pilot groups')
    result = []
    for key in sorted(groups):
        members = groups[key]
        require(set(members) == set(profiles), 'Missing a profile; do not omit failed fits')
        jobs = [members[name] for name in profiles]
        indices = [(row['index'], row['r'], row['fitLo'], row['fitHi']) for row in jobs[0]['slices']]
        require(bool(indices) and len(set(row[0] for row in indices)) == len(indices), 'Invalid pilot slices')
        for job in jobs:
            require([(row['index'], row['r'], row['fitLo'], row['fitHi']) for row in job['slices']]
                    == indices and job['rows'] == len(indices), 'Profiles changed saved R/window selection')
        expected_rows = 2 if key[0] == 'he4' else 1
        require(len(indices) == expected_rows, 'Wrong representative-R coverage')
        result.append(dict(target=key[0], background=key[1], freePolicy=key[2],
                           name='_'.join(key), jobs=jobs,
                           pages=2*((len(profiles)+1)//2)*len(indices)))
    require(sum(group['pages'] for group in result) == contract['fits'],
            'Wrong total page count for the explicitly declared mode')
    return result


def quality_row(job, item, view_cfg, helper):
    require(type(item['numericalValid']) is bool, 'Saved numerical flag must be explicit bool')
    chi2, ndf = item.get('pearson'), item.get('conditionalNdf')
    finite_chi2 = isinstance(chi2, (int, float)) and math.isfinite(chi2)
    valid_ndf = isinstance(ndf, (int, float)) and math.isfinite(ndf) and ndf > 0
    ratio = chi2/ndf if finite_chi2 and valid_ndf else None
    low, high = helper.focus_extent(item, job['target'], view_cfg)
    require(high > low and item['fitHi'] > item['fitLo'], 'Empty saved display/fit extent')
    return [job['profile'], item['index'], job['target'], helper.display_number(chi2),
            helper.display_number(ndf), helper.display_number(ratio), int(item['numericalValid']),
            0, 0, PENDING, item['fitLo'], item['fitHi'], low, high]


def verify_map(signatures):
    for path, expected in signatures.items():
        require(sha(Path(path)) == expected, 'Immutable pilot input changed: ' + str(path))


def prepare(pilot, mainconf, output, mode='fixed'):
    require(not output.exists() and pilot not in output.parents,
            'Use a new output folder, never inside the completed pilot')
    plan_path, done_path = pilot/'plan.json', pilot/'completed.json'
    plan, completed = read(plan_path), read(done_path)
    contract = mode_contract(plan['config'], mode)
    require(plan.get('productionEnabled') is False
            and completed.get('status') == 'completed'
            and completed.get('productionEnabled') is False
            and completed.get('immutableInputsUnchanged') is True
            and completed.get('oldWindowsAndGatesPreserved') is True
            and completed.get('totalJobs') == contract['jobs']
            and completed.get('totalFits') == contract['fits'],
            'Require the completed full diagnostic pilot, not pending/partial jobs')
    require(mainconf == owned(plan['explicitMainconf']), 'Explicit mainconf must match pilot runtime')
    groups = grouped_jobs(plan, mode)
    helper_path = QUALITY/'render_quality_diagnostics.py'
    helper = module('coretail_quality_view_helper', helper_path)
    view_settings = QUALITY/'render_quality_settings.yaml'
    view_cfg = yaml.safe_load(view_settings.read_text())
    require(view_cfg['focusSpecies']['d'] == FOCUS['d']
            and view_cfg['focusSpecies']['t'] == FOCUS['t']
            and view_cfg['focusSpecies']['he4'] == FOCUS['he4'], 'Unexpected nuclear-neighborhood focus')
    renderer = owned(view_cfg['renderer'])
    dependency = owned(view_cfg['rendererDependency'])
    require('supportKnown' in renderer.read_text(),
            'Renderer must display unknown support as not evaluated, never as FAIL')
    wrapper = owned(plan['config']['runtimeWrapper'])
    immutable = {str(plan_path): sha(plan_path), str(done_path): sha(done_path)}
    immutable.update(plan['immutableInputs'])
    immutable.update({str(pilot/name): expected
                      for name, expected in completed.get('sourceAndCompiledArtifacts', {}).items()})
    sources = {Path(__file__).resolve(), Path(__file__).resolve().with_name('test_render_pilot.py'),
               helper_path, view_settings, renderer, dependency, wrapper, mainconf}
    main_cfg = yaml.safe_load(mainconf.read_text())
    sources.update(owned(Path('config')/ref) for ref in main_cfg.values())
    immutable.update({str(path): sha(path) for path in sources})
    prepared_books = []
    for group in groups:
        model_rows, slices, quality_rows = [], None, []
        for job in group['jobs']:
            folder = owned(job['folder'])
            job_done_path = folder/'completed.json'
            job_done = read(job_done_path)
            require(job_done.get('executionComplete') and job_done.get('rawReadback'), 'Incomplete pilot subjob')
            for name, expected in job_done['products'].items():
                immutable[str(folder/name)] = expected
            immutable[str(job_done_path)] = sha(job_done_path)
            fit = read(folder/'cutoff_fits.json')
            require(fit['targetSpecies'] == group['target']
                    and fit['backgroundModel'] == group['background']
                    and fit['responseModel'] == contract['responseModel']
                    and fit['sampleHistogram'] == 'ALL'
                    and fit['profileName'] == job['profile']
                    and fit['productionEnabled'] is False, 'Wrong saved pilot fit identity')
            selected = [(item['index'], item['r']) for item in fit['slices']]
            expected_slices = [(item['index'], item['r']) for item in job['slices']]
            require(selected == expected_slices, 'Saved pilot slices differ from plan')
            require(slices is None or slices == selected, 'Profiles have different native indices')
            slices = selected
            model_rows.append([job['profile'], str(folder/'cutoff_fits.root')])
            planned = {item['index']: item for item in job['slices']}
            for item in fit['slices']:
                require((item['fitLo'], item['fitHi']) ==
                        (planned[item['index']]['fitLo'], planned[item['index']]['fitHi']),
                        'Pilot changed old saved fit window')
                quality_rows.append(quality_row(job, item, view_cfg, helper))
        prepared_books.append((group, model_rows, slices, quality_rows))
    verify_map(immutable)
    output.mkdir(parents=True)
    source = output/'source'
    source.mkdir()
    # Both C++ files must be adjacent: the diagnostic macro includes the common painter.
    for original in sources:
        destination = source/original.relative_to(PROJECT) if original.suffix == '.yaml' else source/original.name
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(original), str(destination))
        require(sha(destination) == immutable[str(original)], 'Renderer snapshot differs')
    table(output/'settings.tsv', (), view_cfg['plots'].items())
    books = []
    for group, models, slices, quality in prepared_books:
        name = group['name']
        table(output/(name+'_models.tsv'), ('label', 'path'), models)
        table(output/(name+'_slices.tsv'), ('index', 'rigidity'), slices)
        table(output/(name+'_quality.tsv'), QUALITY_HEADER, quality)
        books.append({key: group[key] for key in ('name', 'target', 'background', 'freePolicy', 'pages')})
    save(output/'plan.json', dict(schemaVersion=1, pilot=str(pilot), mainconf=str(mainconf),
         mode=mode, books=books, totalPages=contract['fits'],
         immutableInputs=immutable, plots=view_cfg['plots'],
         runtimeWrapper=str(wrapper), renderer=str(source/renderer.name),
         executionTimeoutSeconds=view_cfg['executionTimeoutSeconds'],
         refitted=False, normalized=False, productionEnabled=False, supportCertified=False,
         supportPolicy=PENDING, oldSavedWindowsAndGatesPreserved=True,
         profileParameters=plan['config']['profiles'],
         parameterInterpretation=contract['parameterMeaning']))
    return read(output/'plan.json'), helper


def render(output, plan, helper):
    results = []
    for book in plan['books']:
        name = book['name']
        pdf = output/(name+'.pdf')
        invocation = plan['renderer']+'+('+','.join(helper.quoted(path) for path in
            (output/(name+'_models.tsv'), output/(name+'_slices.tsv'), output/'settings.tsv',
             output/(name+'_quality.tsv'), pdf))+')'
        command = ['bash', plan['runtimeWrapper'], plan['mainconf'],
                   'root4star', '-l', '-b', '-q', invocation]
        save(output/(name+'_command.json'), dict(cwd=str(PROJECT), argv=command))
        with (output/(name+'_render.log')).open('x') as stream:
            subprocess.run(command, cwd=str(PROJECT), stdout=stream, stderr=subprocess.STDOUT,
                           check=True, timeout=plan['executionTimeoutSeconds'])
        text = (output/(name+'_render.log')).read_text()
        require('QUALITY_DIAGNOSTICS_COMPLETE pages='+str(book['pages']) in text,
                'Missing closed renderer completion marker')
        require(text.count('SPECTRUM_VISIBLE_AXES ') == book['pages']*4,
                'Missing a per-panel axis-range audit')
        helper.verify_panel_policy(text, book['pages'], plan['plots'])
        command = ['gs', '-dBATCH', '-dNOPAUSE', '-dPDFSTOPONERROR', '-sDEVICE=nullpage', str(pdf)]
        check = subprocess.run(command, stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                               universal_newlines=True, check=True)
        with (output/(name+'_readback.log')).open('x') as stream:
            stream.write(check.stdout)
        import re
        match = re.search(r'Processing pages 1 through ([0-9]+)\.', check.stdout)
        require(match and int(match.group(1)) == book['pages'] and 'Error:' not in check.stdout,
                'Closed PDF readback/page count failed')
        results.append(dict(name=name, pages=book['pages'], bytes=pdf.stat().st_size,
                            sha256=sha(pdf), closedReadback=True))
    verify_map(plan['immutableInputs'])
    save(output/'completed.json', dict(status='completed',
         completedJST=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
         mode=plan['mode'], totalBooks=len(results), totalPages=sum(item['pages'] for item in results),
         results=results, immutableInputsUnchanged=True, refitted=False, normalized=False,
         productionEnabled=False, supportCertified=False, failedFitsOmitted=False,
         supportPolicy=PENDING,
         outputHashes={str(path.relative_to(output)): sha(path)
                       for path in output.rglob('*') if path.is_file()}))
    print('CORETAIL_PILOT_PDFS_COMPLETE', len(results), sum(item['pages'] for item in results), output)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('pilot-dir', 'mainconf', 'output-dir'):
        parser.add_argument('--'+name, required=True)
    parser.add_argument('--mode', choices=('fixed', 'free_coretail'), default='fixed',
                        help='Free mode requires matching explicit adapter/engine/model declarations')
    args = parser.parse_args()
    pilot, mainconf, output = (owned(path) for path in (args.pilot_dir, args.mainconf, args.output_dir))
    plan, helper = prepare(pilot, mainconf, output, args.mode)
    render(output, plan, helper)


if __name__ == '__main__':
    main()
