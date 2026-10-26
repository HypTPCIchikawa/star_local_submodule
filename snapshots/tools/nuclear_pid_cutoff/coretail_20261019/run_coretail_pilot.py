#!/usr/bin/env python3
"""Small immutable core+tail pilot: prepare, run explicit ROOT5 jobs, then audit."""
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


def require(ok, message):
    if not ok:
        raise ValueError(message)


def path(value, exists=True):
    candidate = Path(value)
    candidate = (candidate if candidate.is_absolute() else PROJECT / candidate).resolve()
    require(PROJECT in candidate.parents, 'Expected a project-owned pilot path')
    if exists:
        require(candidate.exists(), 'Missing explicit input: ' + str(candidate))
    return candidate


def sha(source):
    result = hashlib.sha256()
    with source.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1048576), b''):
            result.update(chunk)
    return result.hexdigest()


def read(source):
    return json.loads(source.read_text())


def save(destination, value):
    with destination.open('x') as stream:
        json.dump(value, stream, indent=2, sort_keys=True, allow_nan=False)
        stream.write('\n')


def rows(source):
    with source.open() as stream:
        return list(csv.DictReader(stream, delimiter='\t'))


def table(destination, header, values):
    with destination.open('x', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        if header:
            writer.writerow(header)
        writer.writerows(values)


def stamp():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def settings_tsv(source):
    values = {}
    for line in source.read_text().splitlines():
        key, value = line.split('\t')
        require(key not in values, 'Duplicate old math setting')
        values[key] = value
    return values


def quote(value):
    text = str(value)
    require('"' not in text and '\\' not in text and '\n' not in text, 'Unsupported ROOT path quoting')
    return '"' + text + '"'


def identifier(*parts):
    require(parts and all(re.fullmatch(r'[A-Za-z][A-Za-z0-9_]*', value) for value in parts),
            'Pilot identifiers must be simple labels, never paths')
    return '_'.join(parts)


def select_slices(slices, rigidity, tolerance):
    require(slices and rigidity and tolerance > 0, 'Missing explicit native slices')
    selected = []
    for requested in rigidity:
        require(math.isfinite(requested) and requested > 0, 'Invalid representative rigidity')
        item = min(slices, key=lambda item: abs(item['r'] - requested))
        require(abs(item['r'] - requested) < tolerance, 'Requested R is not a saved native bin')
        selected.append(item)
    require(len({item['index'] for item in selected}) == len(selected), 'Repeated native pilot slice')
    return selected


def prepare(settings, output, explicit_mainconf):
    cfg = yaml.safe_load(settings.read_text())
    require(cfg['schemaVersion'] == 1, 'Unsupported pilot schema')
    main = path(explicit_mainconf)
    require(main == path(cfg['runtimeMainconf']), 'Explicit mainconf differs from settings')
    require(not output.exists() and output.parent == path(cfg['outputBase']), 'Use a new pilot sibling output')
    require(set(cfg['backgrounds']) == {'linear', 'exponential'}, 'Pilot must compare both saved backgrounds')
    profiles = cfg['profiles']
    adapter_function = cfg.get('adapterFunction', Path(cfg['adapter']).stem)
    identifier(adapter_function)
    require(adapter_function == Path(cfg['adapter']).stem, 'Adapter entry point must match the explicit macro basename')
    if cfg.get('responseStudyMode') == 'free_coretail':
        require(cfg['responseModel'] == 'gaussian_coretail_free_components'
                and adapter_function == 'FitFreeCoreTailPilot'
                and Path(cfg['engine']).name == 'FreeGaussianCoreTail.h', 'Wrong explicit free-coretail adapter/engine/response')
    require(profiles['gaussian_control']['fraction'] == 0, 'Gaussian control must have zero tail')
    sources = [settings, Path(__file__).resolve(), path(cfg['adapter']), path(cfg['adapterDependency']),
               path(cfg['engine']), path(cfg['runtimeWrapper']), main,
               PROJECT / 'script/analysis_info_helper.py']
    sources += [path(value) for value in cfg['engineDependencies']]
    main_cfg = yaml.safe_load(main.read_text())
    sources.append(path(Path('config') / main_cfg['analysis']))
    source_meta = {str(source): sha(source) for source in sources}
    cache = path(cfg['cache'])
    cache_done = cache.parent / 'completed.json'
    cache_digest = sha(cache)
    require(read(cache_done)['products'][cache.name] == cache_digest, 'Native ALL cache differs from its completed extraction')
    immutable = {str(cache): cache_digest, str(cache_done): sha(cache_done)}
    pending = []
    for species, rigidity in cfg['representativeRigidity'].items():
        require(species in SPECIES and rigidity, 'Invalid pilot target')
        for background in cfg['backgrounds']:
            old = path(Path(cfg['oldRunsBase']) / cfg['oldRunPattern'].format(species=species, background=background))
            old_done = read(old / 'completed.json')
            require(old_done['executionComplete'] and old_done['productionEnabled'] is False, 'Incomplete old Gaussian reference')
            fit = read(old / 'cutoff_fits.json')
            require(fit['targetSpecies'] == species and fit['responseModel'] == 'gaussian'
                    and fit['backgroundModel'] == background and fit['sampleHistogram'] == 'ALL', 'Wrong old reference identity')
            selected = select_slices(fit['slices'], rigidity, cfg['axisTolerance'])
            for name in ('cutoff_fits.json', 'cutoff_fits.root', 'completed.json', 'settings.tsv', 'shapes.tsv', 'gates.tsv'):
                actual = sha(old / name)
                if name != 'completed.json':
                    require(actual == old_done['products'][name], 'Changed saved Gaussian reference: ' + str(old / name))
                immutable[str(old / name)] = actual
            shape_rows, gate_rows = rows(old / 'shapes.tsv'), rows(old / 'gates.tsv')
            indices = {item['index'] for item in selected}
            shape_rows = [row for row in shape_rows if int(row['index']) in indices]
            gate_rows = [row for row in gate_rows if int(row['index']) in indices]
            require(len(shape_rows) == 12 * len(indices) and len(gate_rows) == 12 * len(indices), 'Incomplete old GL3 shape/gate subset')
            for policy, policies in cfg['freePolicies'].items():
                free = policies[species]
                require(species in free and len(set(free)) == len(free) and set(free) <= set(SPECIES), 'Invalid explicit free-species policy')
                tail = cfg['tailPolicies'][policy][species] if 'tailPolicies' in cfg else cfg['tailSpecies'][species]
                require(tail and len(set(tail)) == len(tail) and set(tail) <= set(SPECIES), 'Invalid explicit tail species')
                if cfg.get('responseStudyMode') == 'free_coretail':
                    require(free == tail, 'Free-coretail hypotheses must match the explicitly free core hypotheses')
                for profile, shape in profiles.items():
                    require(0 <= shape['fraction'] < .5 and shape['widthRatio'] >= 1
                            and (shape['fraction'] == 0 or shape['widthRatio'] > 1)
                            and all(math.isfinite(value) for value in shape.values()), 'Invalid fixed profile')
                    job_id = identifier(species, background, policy, profile)
                    numeric = settings_tsv(old / 'settings.tsv')
                    numeric.update(responseModel=cfg.get('responseModel', 'gaussian_coretail_fixed_profile'), freeSpecies=','.join(free),
                                   tailSpecies=','.join(tail), profileName=profile,
                                   tailFraction=str(shape['fraction']), tailOffsetSeedSigma=str(shape['offsetSeedSigma']),
                                   tailWidthRatio=str(shape['widthRatio']))
                    if 'responseStudyMode' in cfg:
                        numeric['responseStudyMode'] = cfg['responseStudyMode']
                    pending.append((job_id, species, background, policy, profile, old, selected,
                                    numeric, shape_rows, gate_rows))
    output.mkdir()
    source_root = output / 'source'
    for original, digest in source_meta.items():
        original = Path(original)
        destination = source_root / original.relative_to(PROJECT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(original), str(destination))
        require(sha(original) == digest and sha(destination) == digest, 'Source changed during snapshot')
    jobs = []
    for job_id, species, background, policy, profile, old, selected, numeric, shapes, gates in pending:
        folder = output / job_id
        folder.mkdir()
        table(folder / 'settings.tsv', (), numeric.items())
        table(folder / 'shapes.tsv', tuple(shapes[0]), [tuple(row.values()) for row in shapes])
        table(folder / 'gates.tsv', tuple(gates[0]), [tuple(row.values()) for row in gates])
        copied_macro = source_root / Path(cfg['adapter'])
        includes = source_root / 'tools/nuclear_pid_tof'
        invocation = adapter_function + '(' + ','.join(quote(value) for value in
            (path(cfg['cache']), folder, folder / 'settings.tsv', folder / 'shapes.tsv', folder / 'gates.tsv')) + ');'
        # The CINT entry controls ACLiC only; heavy math lives in the compiled adapter.
        bootstrap = 'void RunCoreTailPilot(){\n'
        bootstrap += 'gSystem->AddIncludePath(' + quote('-I' + str(includes)) + ');\n'
        bootstrap += 'Int_t error=0;gROOT->ProcessLine(' + quote('.L ' + str(copied_macro) + '+') + ',&error);if(error){gSystem->Exit(2);return;}\n'
        # Escape the inner C++ call as one CINT string, not shell interpolation.
        bootstrap += 'gROOT->ProcessLine(' + json.dumps(invocation) + ',&error);if(error)gSystem->Exit(2);\n}\n'
        with (folder / 'RunCoreTailPilot.C').open('x') as stream:
            stream.write(bootstrap)
        for name in ('settings.tsv', 'shapes.tsv', 'gates.tsv', 'RunCoreTailPilot.C'):
            immutable[str(folder / name)] = sha(folder / name)
        jobs.append(dict(id=job_id, target=species, background=background, freePolicy=policy,
                         profile=profile, oldRun=str(old), slices=selected, rows=len(selected),
                         folder=str(folder), referenceParityRequired=policy == 'target_only' and profile == 'gaussian_control'))
    snapshots = {str(p): sha(p) for p in source_root.rglob('*') if p.is_file()}
    immutable.update(snapshots)
    save(output / 'plan.json', dict(schemaVersion=1, preparedJST=stamp(), config=cfg,
         explicitMainconf=str(main), jobs=jobs, immutableInputs=immutable, sourceInputs=source_meta,
         productionEnabled=False, oldWindowsAndGatesPreserved=True,
         profileScanIsNotIndependentGoFPvalue=True))
    print('CORETAIL_PILOT_PREPARED', len(jobs), sum(job['rows'] for job in jobs), output, flush=True)


def unchanged(plan):
    for original, digest in plan['immutableInputs'].items():
        require(sha(Path(original)) == digest, 'Immutable pilot input changed: ' + original)


def validate_result_identity(result, cfg, job):
    require(result['responseModel'] == cfg.get('responseModel', 'gaussian_coretail_fixed_profile')
            and result['sampleHistogram'] == 'ALL' and result['targetSpecies'] == job['target']
            and result['backgroundModel'] == job['background'] and result['profileName'] == job['profile'],
            'Pilot output has the wrong explicit response/target/background/profile')
    if cfg.get('responseStudyMode') == 'free_coretail':
        require(result['responseStudyMode'] == 'free_coretail'
                and result['fullUncertaintyValidated'] is False
                and result['purityPolicy'] == 'Not estimated in this response-quality pilot; null is not zero',
                'Free response-quality pilot must not claim validated grouped purity')
        for row in result['slices']:
            require(set(row['components']) == set(SPECIES), 'Free-pilot folded isotope schema differs')
            for component in row['components'].values():
                selection = component['selection']
                require(selection['available'] is False
                        and selection['reason'] == 'group_purity_not_evaluated_in_response_quality_pilot'
                        and all(selection[key] is None for key in
                                ('purity', 'signalSelected', 'totalSelected', 'backgroundSelected',
                                 'conditionalError', 'purityLowerIncludingOmittedBound',
                                 'selectedColumnFractions', 'selectedColumnCounts')),
                        'Unavailable grouped purity must be null, not a fabricated number')


def run(output, selection):
    require(not (output / 'superseded.json').exists(),
            'This pilot was superseded by the user before fitting; consult its superseded.json and current Todo')
    plan = read(output / 'plan.json')
    unchanged(plan)
    if selection == 'pending':
        jobs = []
        for job in plan['jobs']:
            folder = Path(job['folder'])
            if (folder / 'completed.json').exists():
                done = read(folder / 'completed.json')
                require(done['executionComplete'], 'Invalid previous pilot completion')
                for name, digest in done['products'].items():
                    require(sha(folder / name) == digest, 'Previous completed pilot output changed')
            else:
                require(not (folder / 'run.log').exists(), 'Failed/attempted pilot is preserved; choose a new output instead of retrying it')
                jobs.append(job)
        if not jobs:
            print('CORETAIL_PILOT_NO_PENDING_JOBS', flush=True)
            return
    else:
        jobs = plan['jobs'] if selection == 'all' else [job for job in plan['jobs'] if job['id'] == selection]
        require(jobs, 'Unknown explicit pilot job')
    for job in jobs:
        folder = Path(job['folder'])
        require(not (folder / 'completed.json').exists() and not (folder / 'run.log').exists(), 'Pilot job already attempted; preserve it and prepare a new output')
        cfg = plan['config']
        wrapper = path(cfg['runtimeWrapper'])
        require(sha(wrapper) == plan['sourceInputs'][str(wrapper)], 'Runtime wrapper changed since snapshot')
        main = path(plan['explicitMainconf'])
        require(sha(main) == plan['sourceInputs'][str(main)], 'Explicit runtime mainconf changed since snapshot')
        command = ['bash', str(wrapper), str(main), 'root4star', '-l', '-b', '-q', str(folder / 'RunCoreTailPilot.C')]
        save(folder / 'command.json', dict(cwd=str(PROJECT), argv=command))
        with (folder / 'run.log').open('x') as stream:
            subprocess.run(command, cwd=str(PROJECT), stdout=stream, stderr=subprocess.STDOUT,
                           check=True, timeout=cfg['executionTimeoutSeconds'])
        text = (folder / 'run.log').read_text()
        require('CORETAIL_PILOT_COMPLETE rows=' + str(job['rows']) in text
                and 'NUCLEAR_CUTOFF_READBACK_PASS rows=' + str(job['rows']) in text, 'Missing closed-ROOT readback/completion marker')
        result = read(folder / 'cutoff_fits.json')
        require(len(result['slices']) == job['rows'] and result['productionEnabled'] is False, 'Unexpected pilot result schema')
        validate_result_identity(result, cfg, job)
        if cfg.get('responseStudyMode') == 'free_coretail':
            require('FREE_CORETAIL_DYNAMIC_READBACK_PASS rows=' + str(job['rows']) in text,
                    'Missing closed-ROOT dynamic-column readback marker')
        unchanged(plan)
        save(folder / 'completed.json', dict(executionComplete=True, productionEnabled=False,
             completedJST=stamp(), referenceParityRequired=job['referenceParityRequired'],
             protectedInputsUnchanged=True, rawReadback=True,
             products={p.name: sha(p) for p in folder.iterdir() if p.is_file()}))
        print('CORETAIL_PILOT_JOB_COMPLETE', job['id'], flush=True)


def close(a, b, tolerance):
    return math.isfinite(a) and math.isfinite(b) and abs(a-b) <= tolerance * max(1., abs(a), abs(b))


def check(output):
    plan = read(output / 'plan.json')
    unchanged(plan)
    comparisons = []
    for job in plan['jobs']:
        folder = Path(job['folder'])
        done = read(folder / 'completed.json')
        require(done['executionComplete'] and done['rawReadback'], 'Unfinished pilot job')
        for name, digest in done['products'].items():
            require(sha(folder / name) == digest, 'Pilot output changed')
        new = read(folder / 'cutoff_fits.json')
        validate_result_identity(new, plan['config'], job)
        old = {row['index']: row for row in read(Path(job['oldRun']) / 'cutoff_fits.json')['slices']}
        for row in new['slices']:
            reference = old[row['index']]
            require(row['lo'] == reference['lo'] and row['hi'] == reference['hi']
                    and row['fitLo'] == reference['fitLo'] and row['fitHi'] == reference['fitHi'], 'Pilot changed old bin/window')
            item = dict(job=job['id'], r=row['r'], numericalValid=row['numericalValid'],
                        chi2=row['pearson'], ndf=row['conditionalNdf'], oldChi2=reference['pearson'],
                        oldNumericalValid=reference['numericalValid'], controlParityRequired=job['referenceParityRequired'])
            if job['referenceParityRequired']:
                tolerance = plan['config']['parityRelativeTolerance']
                item['controlParity'] = all(close(row[key], reference[key], tolerance)
                    for key in ('pearson', 'dataDeviance', 'conditionalNdf')) and row['status'] == reference['status']
                require(item['controlParity'], 'Zero-tail target-only control does not reproduce old Gaussian; do not claim isolated response improvement')
            comparisons.append(item)
    save(output / 'completed.json', dict(status='completed', completedJST=stamp(),
         productionEnabled=False, oldWindowsAndGatesPreserved=True, immutableInputsUnchanged=True,
         totalJobs=len(plan['jobs']), totalFits=len(comparisons), comparisons=comparisons,
         profileScanIsNotIndependentGoFPvalue=True,
         sourceAndCompiledArtifacts={str(p.relative_to(output)): sha(p)
             for p in (output / 'source').rglob('*') if p.is_file()},
         limitations=plan['config'].get('studyLimitations',
                     ['Frozen-old-gate purity is a reference estimate, not new-polynomial purity.',
                      'Fixed tail-profile scans do not fit or validate a smooth detector calibration.',
                      'Core width is not PID sigma or mixture RMS; no production coefficient adopted.'])))
    print('CORETAIL_PILOT_ALL_COMPLETE', len(comparisons), flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'check'))
    parser.add_argument('--settings')
    parser.add_argument('--mainconf')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--job', default='all')
    args = parser.parse_args()
    output = path(args.output_dir, exists=args.action != 'prepare')
    if args.action == 'prepare':
        require(args.settings and args.mainconf, 'prepare requires explicit --settings and --mainconf')
        prepare(path(args.settings), output, args.mainconf)
    elif args.action == 'run':
        run(output, args.job)
    else:
        check(output)


if __name__ == '__main__':
    main()
