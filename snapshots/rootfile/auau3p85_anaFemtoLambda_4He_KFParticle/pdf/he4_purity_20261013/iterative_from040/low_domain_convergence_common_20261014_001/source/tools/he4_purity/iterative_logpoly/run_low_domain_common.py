#!/usr/bin/env python3
"""Study-only R=.4-.7 COMMON diagnostic; never change the completed full study."""
import argparse
import copy
from datetime import datetime, timedelta, timezone
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
INDICES = list(range(30))
DOMAIN = [0.4, 0.7]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study-folder', type=Path, required=True)
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--domain-source', type=Path, required=True)
    parser.add_argument('--settings', type=Path, default=HERE/'low_domain_common_settings.yaml')
    parser.add_argument('--fix-triton-shape', action='store_true')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--no-render', action='store_true')
    args = parser.parse_args()
    # This small dependency provides strict read/hash/seed-floor guards only;
    # the frozen COMMON helper itself is never executed or modified.
    import run_common_refit as guards
    digest, read, save, inside = guards.digest, guards.read_json, guards.save, guards.inside
    study = args.study_folder.resolve(strict=True)
    context_path, completed_path = study/'context.json', study/'completed.json'
    context, completed = read(context_path), read(completed_path)
    if (context['schemaVersion'] != 1 or context['indices'] != list(range(210)) or
            context['productionEnabled'] or not completed['executionComplete'] or completed['productionEnabled']):
        raise ValueError('A completed study-only full 210-slice study is required.')
    protected = {}
    for encoded, expected in context['protected'].items():
        path = Path(encoded)
        path = inside(path if path.is_absolute() else PROJECT/path, PROJECT)
        if str(path) in protected and protected[str(path)] != expected:
            raise ValueError('Conflicting protected source path')
        protected[str(path)] = expected
    protected[str(context_path)], protected[str(completed_path)] = digest(context_path), digest(completed_path)
    for relative, expected in completed['products'].items():
        path = guards.checked_product(study, relative)
        if digest(path) != expected:
            raise ValueError('Completed-study product changed: '+relative)
        protected[str(path)] = expected
    for encoded, expected in context['snapshots'].items():
        path = inside(Path(encoded), study/'source')
        if digest(path) != expected:
            raise ValueError('Frozen source changed: '+str(path))
        protected[str(path)] = expected
    for path, expected in protected.items():
        if digest(Path(path)) != expected:
            raise ValueError('Protected source changed: '+path)
    helper, smoother = HERE/'run_iterative_fit.py', HERE/'smoothing.py'
    for path in (helper, smoother):
        if protected.get(str(path)) != digest(path):
            raise ValueError('Trusted helper differs from completed study: '+str(path))
    sys.path.insert(0, str(HERE))
    spec = importlib.util.spec_from_file_location('low_domain_frozen_helpers', helper)
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    cfg = copy.deepcopy(context['settings'])
    helpers.validate(cfg)  # Validate the ORIGINAL full-domain contract first.
    sample = cfg['sample']
    if (sample['collision'], sample['year'], sample['production']) != ('auau3p85', 2021, 'P24iy'):
        raise ValueError('Only auau3p85 2021 P24iy is authorized.')
    allowed = (PROJECT/sample['outputBase']).resolve()
    inside(study, allowed)
    mainconf = args.mainconf.resolve(strict=True)
    if mainconf != (PROJECT/sample['mainconf']).resolve(strict=True):
        raise ValueError('Explicit mainconf differs from the source study.')
    runtime = yaml.safe_load(mainconf.read_text())
    analysis_path = (PROJECT/'config'/runtime['analysis']).resolve(strict=True)
    analysis = yaml.safe_load(analysis_path.read_text())
    tags = analysis.get('starTag', {})
    if (tags.get('libraryTag'), tags.get('year'), tags.get('productionTag'), tags.get('collision'),
            analysis.get('analysis', {}).get('mode')) != ('SL24y', 'year2021', 'P24iy', 'auau3p85', 'fxtmult'):
        raise ValueError('Runtime must be SL24y / 2021 P24iy / auau3p85 / FXT.')
    source, cache = (PROJECT/sample['input']).resolve(strict=True), Path(context['cache']).resolve(strict=True)
    initial = helpers.fingerprint(source)
    if initial != context['source']:
        raise ValueError('Source ROOT fingerprint changed.')
    source_summary = (PROJECT/sample['sourceSummary']).resolve(strict=True)
    state = read(source_summary)
    if state['phase'] != 'complete' or state['provenance']['input'] != initial:
        raise ValueError('Verified source run state differs.')
    settings_path, domain_path = args.settings.resolve(strict=True), args.domain_source.resolve(strict=True)
    settings = yaml.safe_load(settings_path.read_text())
    if (set(settings) != {'schemaVersion', 'domain', 'engine'} or settings['schemaVersion'] != 1 or
            settings['domain'] != DOMAIN or not isinstance(settings['engine'], dict)):
        raise ValueError('Settings must declare exactly the authorized R=.4-.7 diagnostic domain.')
    override = {}
    for key, value in settings['engine'].items():
        if (key not in guards.NUMERICAL_KEYS or isinstance(value, bool) or not isinstance(value, (int, float)) or
                not math.isfinite(value) or value <= 0 or
                ((key.endswith('Calls') or key.endswith('Iterations')) and not isinstance(value, int))):
            raise ValueError('Invalid numerical-only override: '+key)
        override[key] = dict(original=cfg['engine'][key], override=value)
        cfg['engine'][key] = value
    cfg['engine']['rigidityMin'], cfg['engine']['rigidityMax'] = DOMAIN
    diagnostic = read(domain_path)
    if diagnostic['domain'] != DOMAIN:
        raise ValueError('Coefficient diagnostic domain differs.')
    coefficients = []
    for species in ('he3', 'he4'):
        item = diagnostic['coefficients'][species]
        if item['nSlices'] != len(INDICES):
            raise ValueError('Coefficient source does not use the same 30 FREE slices.')
        coefficients += item['meanCoefficients']+item['sigmaCoefficients']
    minima = guards.require_shared_seed(coefficients, cfg['engine'])
    free_path = study/'00_local_free/fit_result.json'
    free_all = read(free_path)
    if [row['index'] for row in free_all['slices']] != list(range(210)):
        raise ValueError('FREE baseline must cover the full ordered input.')
    free = copy.deepcopy(free_all)
    free['slices'] = free['slices'][:len(INDICES)]
    for index, row in enumerate(free['slices']):
        if (row['index'] != index or abs(row['lo']-(.4+.01*index)) > 1.e-9 or
                abs(row['hi']-(.41+.01*index)) > 1.e-9):
            raise ValueError('FREE subset rigidity edges differ.')
    for key in ('dataDeviance', 'penaltyLoss', 'totalLoss'):
        free[key] = math.fsum(row[key] for row in free['slices'])
    old = (PROJECT/sample['oldStudy']).resolve(strict=True)
    seed_rows, seed_details = helpers.old_seeds(read(old/'purity.json'), cfg['engine'], INDICES)
    seed_rows = helpers.seed_update(free, seed_rows)
    for seed, row in zip(seed_rows, free['slices']):
        seed[9] = int(row['t']['included'])
        low, high = helpers.used_window(seed[3], seed[4], cfg['engine'], .5)
        if abs(low-row['fitLo']) > 1.e-9 or abs(high-row['fitHi']) > 1.e-9:
            raise ValueError('Low-domain diagnostic would change the FREE fit window.')
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError('Existing output will not be overwritten.')
    inside(output, allowed)
    if study in output.parents:
        raise ValueError('Output must be a new sibling study.')
    output.mkdir(parents=True)
    (output/'logs').mkdir()
    snapshots = {}
    dependencies = ['tools/he4_purity/iterative_logpoly/IterativeHeliumFit.C',
        'tools/he4_purity/iterative_logpoly/RenderIterativeHelium.C',
        'tools/he4_purity/global_logpoly/GlobalHeliumFit.C', 'StMaker/common/NuclearIdDeDxVsMom.h',
        'tools/he4_purity/iterative_logpoly/run_iterative_fit.py',
        'tools/he4_purity/iterative_logpoly/smoothing.py', 'script/singularity_root_batch.sh']
    for relative in dependencies:
        original = study/'source'/relative
        expected = context['snapshots'].get(str(original))
        if expected is None or digest(original) != expected:
            raise ValueError('Frozen dependency changed: '+relative)
        target = output/'source'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, target)
        snapshots[str(target)] = expected
    additional = [Path(__file__).resolve(), HERE/'run_common_refit.py', HERE/'CommonHeliumFixedT.C',
                  settings_path, domain_path, mainconf, analysis_path, source_summary]
    for original in additional:
        expected = digest(original)
        protected[str(original)] = expected
        relative = original.relative_to(PROJECT) if PROJECT in original.parents else Path('external')/original.name
        target = output/'source'/relative
        if target.exists():
            if digest(target) != expected:
                raise ValueError('Snapshot path collision')
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, target)
        snapshots[str(target)] = expected
    engine = dict(cfg['engine'], phase='common', penaltyStrength=0.)
    helpers.tsv(output/'settings.tsv', engine.items())
    helpers.tsv(output/'seeds.tsv', seed_rows)
    helpers.tsv(output/'targets.tsv', [['coefficients']+coefficients])
    (output/'effective_settings.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False))
    t_policy = ('Fixed FREE empirical T projected mean/width; T yield and BG free; conditional diagnostic.'
                if args.fix_triton_shape else 'FREE nuisance start; original numerical T-information policy.')
    save(output/'context.json', dict(schemaVersion=1, indices=INDICES, studyDomain=DOMAIN,
        dataScope='Only FREE input rows 0-29, R=.4-.7; diagnostic domain, NOT an analysis selection or package change.',
        sourceStudy=str(study), sourceContextSha256=digest(context_path),
        sourceCompletionSha256=digest(completed_path), input=initial, cache=str(cache), settings=cfg,
        numericalOverrides=override, nuisanceSource=str(free_path), nuisanceSourceSha256=digest(free_path),
        coefficientSource=str(domain_path), coefficientSourceSha256=digest(domain_path),
        initialCoefficients=coefficients, initialSigmaMinima=minima, fixTritonShape=args.fix_triton_shape,
        tConstraintPolicy=t_policy, protected=protected, snapshots=snapshots, productionEnabled=False,
        physicsValidationComplete=False, createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))
    save(output/'seed_audit.json', dict(originalFrom040=seed_details, freeSubsetSummary=helpers.result_summary(free),
        dataScope=DOMAIN, tConstraintPolicy=t_policy))

    def preserve():
        if helpers.fingerprint(source) != initial:
            raise ValueError('Original merged ROOT changed.')
        for path, expected in dict(protected, **snapshots).items():
            if digest(Path(path)) != expected:
                raise ValueError('Protected source/product changed: '+path)

    macro = output/'source/tools/he4_purity/iterative_logpoly'/(
        'CommonHeliumFixedT.C' if args.fix_triton_shape else 'IterativeHeliumFit.C')
    call = str(macro)+'+('+','.join(json.dumps(str(p)) for p in
        (cache, output, output/'settings.tsv', output/'seeds.tsv', output/'targets.tsv'))+')'
    command = ['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
               'root4star', '-l', '-b', '-q', call]
    save(output/'command.json', command)
    preserve()
    print('LOW_DOMAIN_COMMON_PREPARED' if args.prepare_only else 'LOW_DOMAIN_COMMON_START', output, flush=True)
    if args.prepare_only:
        return
    helpers.run_command(command, output/'logs/fit.log', cfg['executionTimeoutSeconds'], 'ITERATIVE_HELIUM_COMPLETE')
    result = read(output/'fit_result.json')
    if result['phase'] != 'common' or [row['index'] for row in result['slices']] != INDICES:
        raise ValueError('Returned COMMON diagnostic scope differs.')
    for row, baseline in zip(result['slices'], free['slices']):
        for key in ('lo', 'hi', 'fitLo', 'fitHi', 'firstFitBin', 'lastFitBin', 'fitCount'):
            if row[key] != baseline[key]:
                raise ValueError('COMMON diagnostic changed FREE data/bin selection: '+key)
    summary = helpers.result_summary(result)
    save(output/'summary.json', dict(summary=summary, freeSubsetSummary=helpers.result_summary(free),
        studyDomain=DOMAIN, tConstraintPolicy=t_policy, productionEnabled=False, physicsValidationComplete=False))
    candidate = helpers.candidate_yaml(result['sharedCoefficients'], cfg, 'common')
    candidate.update(tConstraintPolicy=t_policy, diagnosticDomainOnly=True)
    (output/'analytic_candidate.yaml').write_text(yaml.safe_dump(candidate, sort_keys=False))
    if not args.no_render:
        views = output/'views'
        views.mkdir()
        helpers.renderer_manifest(result, views/'slices.tsv', cfg['display'])
        helpers.tsv(views/'settings.tsv', dict(cfg['display'], expectedSlices=30, allowSparseSlices=0).items())
        renderer = output/'source/tools/he4_purity/iterative_logpoly/RenderIterativeHelium.C'
        title = 'LOW_DOMAIN_COMMON_FREE'+('_T_shape_fixed' if args.fix_triton_shape else '')
        call = str(renderer)+'+('+','.join(json.dumps(str(p)) for p in
            (output/'phase.root', views/'slices.tsv', views/'settings.tsv', views, title))+')'
        helpers.run_command(['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
            'root4star', '-l', '-b', '-q', call], output/'logs/render.log', cfg['executionTimeoutSeconds'],
            'ITERATIVE_HELIUM_RENDER_COMPLETE')
    preserve()
    products = {str(p.relative_to(output)):digest(p) for p in output.rglob('*')
                if p.is_file() and 'source' not in p.relative_to(output).parts and p.name != 'context.json'}
    save(output/'completed.json', dict(executionComplete=True, productionEnabled=False,
        physicsValidationComplete=False, sourceInputPreserved=True, products=products,
        completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        note='Diagnostic completion is NOT fit convergence, isotope identification or validated PID.'))
    print('LOW_DOMAIN_COMMON_COMPLETE', output, summary, flush=True)


if __name__ == '__main__':
    main()
