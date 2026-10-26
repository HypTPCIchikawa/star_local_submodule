#!/usr/bin/env python3
"""Isolated COMMON-only sensitivity; preserve the completed iterative study."""
import argparse
import copy
from datetime import datetime, timedelta, timezone
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import sys

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
NUMERICAL_KEYS = {
    'localMaxCalls', 'localTolerance', 'amplitudeMaxIterations',
    'amplitudeGradientTolerance', 'commonMaxCalls', 'commonTolerance',
    'maximumCommonOuterIterations', 'commonRelativeTolerance',
    'commonCoefficientTolerance', 'meanStep', 'sigmaLogStep',
    'backgroundSlopeStep', 'commonMeanStep', 'commonWidthShapeStep',
    'commonWidthLogStep', 'minimumGaussianAcceptance',
}


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_json(path):
    return json.loads(path.read_text(), parse_constant=lambda v:
                      (_ for _ in ()).throw(ValueError('Nonfinite JSON '+v)))


def save(path, obj):
    if path.exists():
        raise ValueError('Existing output will not be replaced: '+str(path))
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def inside(path, directory):
    path = path.resolve()
    if directory.resolve() not in path.parents:
        raise ValueError('Path is outside the declared directory: '+str(path))
    return path


def checked_product(study, relative):
    path = inside(study/relative, study)
    if not path.is_file():
        raise ValueError('Missing completed-study product '+str(path))
    return path


def apply_override(cfg, path):
    result = copy.deepcopy(cfg)
    audit = {}
    if path is None:
        return result, audit
    doc = yaml.safe_load(path.read_text())
    if (not isinstance(doc, dict) or set(doc) != {'schemaVersion', 'engine'} or
            doc['schemaVersion'] != 1 or not isinstance(doc['engine'], dict) or
            not doc['engine']):
        raise ValueError('Numerical override requires schemaVersion: 1 and a nonempty engine mapping.')
    for key, value in doc['engine'].items():
        if key not in NUMERICAL_KEYS:
            raise ValueError('Only declared numerical settings may be overridden: '+key)
        if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value) or value <= 0:
            raise ValueError('Invalid positive numerical override '+key)
        if key.endswith('Calls') or key.endswith('Iterations'):
            if not isinstance(value, int):
                raise ValueError('Iteration/call limits must be integers: '+key)
        audit[key] = dict(original=result['engine'][key], override=value)
        result['engine'][key] = value
    return result, audit


def require_shared_seed(coefficients, engine):
    if len(coefficients) != 14 or not all(math.isfinite(v) for v in coefficients):
        raise ValueError('Exactly 14 finite shared coefficients are required.')
    lo, hi = math.log10(engine['rigidityMin']), math.log10(engine['rigidityMax'])
    minima = []
    for offset in (4, 11):
        c0, c1, c2 = coefficients[offset:offset+3]
        points = [lo, hi]
        if c2 > 0 and lo < -c1/(2*c2) < hi:
            points.append(-c1/(2*c2))
        minimum = min(c0+x*(c1+x*c2) for x in points)
        if minimum <= engine['minimumSigma']:
            raise ValueError('Shared sigma seed touches the strict optimizer floor; no silent shift is allowed. '
                             'Use initial_qr or explicitly study a documented coefficient-margin variant.')
        minima.append(minimum)
    return minima


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study-folder', type=Path, required=True,
                        help='Completed full 210-slice iterative study, not a running study.')
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--nuisance-start', choices=['free', 'last_regularized'], default='free')
    parser.add_argument('--coefficient-start', choices=['final', 'initial_qr'], default='final')
    parser.add_argument('--settings', type=Path,
                        help='Explicit numerical-only override YAML; schemaVersion: 1, engine: {keys}.')
    parser.add_argument('--fix-triton-shape', action='store_true',
                        help='Conditional sensitivity: fix T mean/width to FREE empirical values; yields/BG remain free.')
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--no-render', action='store_true')
    args = parser.parse_args()
    if args.fix_triton_shape and args.nuisance_start != 'free':
        raise ValueError('Fixed empirical T uses the FREE nuisance start only.')
    study = args.study_folder.resolve(strict=True)
    context_path, completion_path = study/'context.json', study/'completed.json'
    if not completion_path.is_file():
        raise ValueError('Wait for the original full study to complete; it will not be altered or resumed here.')
    context, completed = read_json(context_path), read_json(completion_path)
    if (context['schemaVersion'] != 1 or context['indices'] != list(range(210)) or
            context['productionEnabled'] or not completed['executionComplete'] or completed['productionEnabled']):
        raise ValueError('A completed study-only full 210-slice input is required.')
    original_protected = {}
    for encoded, expected in context['protected'].items():
        path = Path(encoded)
        path = inside(path if path.is_absolute() else PROJECT/path, PROJECT)
        if str(path) in original_protected and original_protected[str(path)] != expected:
            raise ValueError('Conflicting canonical protected path '+str(path))
        original_protected[str(path)] = expected
    protected = dict(original_protected)
    protected[str(context_path)] = digest(context_path)
    protected[str(completion_path)] = digest(completion_path)
    for relative, expected in completed['products'].items():
        path = checked_product(study, relative)
        if digest(path) != expected:
            raise ValueError('Completed product changed: '+relative)
        protected[str(path)] = expected
    for path, expected in context['snapshots'].items():
        checked = inside(Path(path), study/'source')
        if digest(checked) != expected:
            raise ValueError('Frozen study source changed: '+str(checked))
        protected[str(checked)] = expected
    for path, expected in original_protected.items():
        checked = Path(path)
        if digest(checked) != expected:
            raise ValueError('Original protected source/input changed: '+str(checked))
    # Import trusted repo helpers only after proving exact equality to the
    # completed study's frozen copies; no execution from an arbitrary study.
    helper = HERE/'run_iterative_fit.py'
    smoother = HERE/'smoothing.py'
    for path in (helper, smoother):
        if str(path) not in original_protected or digest(path) != original_protected[str(path)]:
            raise ValueError('Helper no longer matches the completed study: '+str(path))
    sys.path.insert(0, str(HERE))
    spec = importlib.util.spec_from_file_location('common_refit_frozen_helpers', helper)
    helpers = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helpers)
    cfg = copy.deepcopy(context['settings'])
    helpers.validate(cfg)
    sample = cfg['sample']
    if (sample['collision'], sample['year'], sample['production']) != ('auau3p85', 2021, 'P24iy'):
        raise ValueError('Only 3.85 GeV 2021 P24iy is authorized.')
    allowed = (PROJECT/sample['outputBase']).resolve()
    inside(study, allowed)
    mainconf = args.mainconf.resolve(strict=True)
    if mainconf != (PROJECT/sample['mainconf']).resolve(strict=True):
        raise ValueError('Explicit mainconf differs from the frozen study.')
    runtime = yaml.safe_load(mainconf.read_text())
    analysis_path = (PROJECT/'config'/runtime['analysis']).resolve(strict=True)
    analysis = yaml.safe_load(analysis_path.read_text())
    tags = analysis.get('starTag', {})
    if (tags.get('libraryTag'), tags.get('year'), tags.get('productionTag'),
            tags.get('collision'), analysis.get('analysis', {}).get('mode')) != (
            'SL24y', 'year2021', 'P24iy', 'auau3p85', 'fxtmult'):
        raise ValueError('Explicit runtime must be SL24y / 2021 P24iy / 3.85 FXT.')
    source = (PROJECT/sample['input']).resolve(strict=True)
    cache = Path(context['cache']).resolve(strict=True)
    initial = helpers.fingerprint(source)
    if initial != context['source']:
        raise ValueError('Merged ROOT input fingerprint changed.')
    source_summary = (PROJECT/sample['sourceSummary']).resolve(strict=True)
    state = read_json(source_summary)
    if state['phase'] != 'complete' or state['provenance']['input'] != initial:
        raise ValueError('Source ROOT no longer matches verified original run state.')
    override_path = args.settings.resolve(strict=True) if args.settings else None
    cfg, override_audit = apply_override(cfg, override_path)
    helpers.validate(cfg)
    history = read_json(study/'history.json')
    if not history or history[-1]['mode'] != 'common':
        raise ValueError('Completed input must end with the raw COMMON refit.')
    common_file = study/history[-1]['phase']/'fit_result.json'
    common_result = read_json(common_file)
    if args.coefficient_start == 'final':
        coefficients = list(common_result['sharedCoefficients'])
        coefficient_file = common_file
    else:
        coefficient_file = study/'initial_common.json'
        coefficients = list(read_json(coefficient_file)['coefficients'])
    minima = require_shared_seed(coefficients, cfg['engine'])
    if args.nuisance_start == 'free':
        nuisance_file = study/'00_local_free/fit_result.json'
    else:
        candidates = [h for h in history if h['mode'] == 'regularized']
        if not candidates:
            raise ValueError('No last_regularized phase exists.')
        nuisance_file = study/candidates[-1]['phase']/'fit_result.json'
    nuisance = read_json(nuisance_file)
    if [r['index'] for r in nuisance['slices']] != list(range(210)):
        raise ValueError('Nuisance input does not cover the same ordered 210 slices.')
    old = (PROJECT/sample['oldStudy']).resolve(strict=True)
    old_audit = read_json(old/'purity.json')
    seed_rows, seed_details = helpers.old_seeds(old_audit, cfg['engine'], list(range(210)))
    seed_rows = helpers.seed_update(nuisance, seed_rows)
    # Preserve the selected data model's absent T components even if an
    # acceptance precision threshold is relaxed. Do not introduce new tails.
    by_index = {r['index']: r for r in nuisance['slices']}
    for row in seed_rows:
        row[9] = int(by_index[row[0]]['t']['included'])
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError('Choose a nonexistent output directory; no study is overwritten.')
    inside(output, allowed)
    if study in output.parents:
        raise ValueError('Output must be a new sibling study, not inside the completed input.')
    output.mkdir(parents=True)
    (output/'logs').mkdir()
    sources = {}
    relative_sources = [
        'tools/he4_purity/iterative_logpoly/IterativeHeliumFit.C',
        'tools/he4_purity/iterative_logpoly/RenderIterativeHelium.C',
        'tools/he4_purity/global_logpoly/GlobalHeliumFit.C',
        'StMaker/common/NuclearIdDeDxVsMom.h',
        'script/singularity_root_batch.sh',
        'tools/he4_purity/iterative_logpoly/run_iterative_fit.py',
        'tools/he4_purity/iterative_logpoly/smoothing.py',
    ]
    for relative in relative_sources:
        original = study/'source'/relative
        expected = context['snapshots'].get(str(original))
        if expected is None or digest(original) != expected:
            raise ValueError('Missing or changed frozen dependency '+relative)
        target = output/'source'/relative
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(original, target)
        if digest(target) != expected:
            raise ValueError('Frozen source copy differs '+relative)
        sources[str(target)] = expected
    added = [Path(__file__).resolve(), HERE/'CommonHeliumFixedT.C', mainconf, analysis_path]
    if override_path:
        added.append(override_path)
    for original in added:
        protected[str(original)] = digest(original)
        relative = original.relative_to(PROJECT) if PROJECT in original.parents else Path('external')/original.name
        target = output/'source'/relative
        if target.exists():
            if digest(target) != digest(original):
                raise ValueError('Source snapshot collision '+str(target))
        else:
            target.parent.mkdir(parents=True, exist_ok=True)
            shutil.copy2(original, target)
        sources[str(target)] = digest(original)
    protected[str(source_summary)] = digest(source_summary)
    engine = dict(cfg['engine'], phase='common', penaltyStrength=0.)
    helpers.tsv(output/'settings.tsv', engine.items())
    helpers.tsv(output/'seeds.tsv', seed_rows)
    helpers.tsv(output/'targets.tsv', [['coefficients']+coefficients])
    (output/'effective_settings.yaml').write_text(yaml.safe_dump(cfg, sort_keys=False))
    t_policy = ('Fixed T mean/width from FREE empirical projected fits, yields and BG free; '
                'conditional sensitivity, NOT a legacy or independent prior.' if args.fix_triton_shape else
                'T and BG start from '+args.nuisance_start+'; T shape freedom follows the same numerical information guard.')
    new_context = dict(schemaVersion=1, input=initial, cache=str(cache),
        sourceStudy=str(study), sourceContextSha256=digest(context_path),
        sourceCompletionSha256=digest(completion_path), settings=cfg,
        numericalOverrides=override_audit, nuisanceStart=args.nuisance_start,
        nuisanceSource=str(nuisance_file), nuisanceSourceSha256=digest(nuisance_file),
        coefficientStart=args.coefficient_start, coefficientSource=str(coefficient_file),
        coefficientSourceSha256=digest(coefficient_file), initialCoefficients=coefficients,
        initialSigmaMinima=minima, tConstraintPolicy=t_policy,
        fixTritonShape=args.fix_triton_shape, protected=protected, snapshots=sources,
        productionEnabled=False, physicsValidationComplete=False,
        createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat())
    save(output/'context.json', new_context)
    save(output/'seed_audit.json', dict(originalFrom040=seed_details,
                                      start=args.nuisance_start, tConstraintPolicy=t_policy))
    macro_name = 'CommonHeliumFixedT.C' if args.fix_triton_shape else 'IterativeHeliumFit.C'
    macro = output/'source/tools/he4_purity/iterative_logpoly'/macro_name
    call = str(macro)+'+('+','.join(json.dumps(str(p)) for p in
                    (cache, output, output/'settings.tsv', output/'seeds.tsv', output/'targets.tsv'))+')'
    command = ['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
               'root4star', '-l', '-b', '-q', call]
    save(output/'command.json', command)

    def preserve():
        if helpers.fingerprint(source) != initial:
            raise ValueError('Original merged ROOT changed during COMMON sensitivity.')
        for path, expected in protected.items():
            if digest(Path(path)) != expected:
                raise ValueError('Original input/source changed during COMMON sensitivity: '+path)
        for path, expected in sources.items():
            if digest(Path(path)) != expected:
                raise ValueError('Consumed source snapshot changed: '+path)

    preserve()
    print('COMMON_REFIT_PREPARED' if args.prepare_only else 'COMMON_REFIT_START', output, flush=True)
    if args.prepare_only:
        return
    helpers.run_command(command, output/'logs/fit.log', cfg['executionTimeoutSeconds'],
                        'ITERATIVE_HELIUM_COMPLETE')
    result = read_json(output/'fit_result.json')
    if result['phase'] != 'common' or [r['index'] for r in result['slices']] != list(range(210)):
        raise ValueError('COMMON refit exported unexpected phase or slices.')
    fitted = result['sharedCoefficients']
    if len(fitted) != 14 or not all(math.isfinite(v) for v in fitted):
        raise ValueError('COMMON refit exported invalid shared coefficients.')
    summary = helpers.result_summary(result)
    save(output/'summary.json', dict(summary=summary, nuisanceStart=args.nuisance_start,
        coefficientStart=args.coefficient_start, numericalOverrides=override_audit,
        tConstraintPolicy=t_policy, freeSummary=read_json(study/'00_local_free/summary.json'),
        originalCommonSummary=helpers.result_summary(common_result),
        productionEnabled=False, physicsValidationComplete=False))
    candidate = helpers.candidate_yaml(fitted, cfg, 'common')
    candidate['tConstraintPolicy'] = t_policy
    (output/'analytic_candidate.yaml').write_text(yaml.safe_dump(candidate, sort_keys=False))
    if not args.no_render:
        views = output/'views'
        views.mkdir()
        helpers.renderer_manifest(result, views/'slices.tsv', cfg['display'])
        helpers.tsv(views/'settings.tsv', dict(cfg['display'], expectedSlices=210, allowSparseSlices=0).items())
        renderer = output/'source/tools/he4_purity/iterative_logpoly/RenderIterativeHelium.C'
        title = 'COMMON_'+args.nuisance_start+('_FREE_T_shape_fixed' if args.fix_triton_shape else '')
        call = str(renderer)+'+('+','.join(json.dumps(str(p)) for p in
                  (output/'phase.root', views/'slices.tsv', views/'settings.tsv', views, title))+')'
        helpers.run_command(['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
                             'root4star', '-l', '-b', '-q', call],
                            output/'logs/render.log', cfg['executionTimeoutSeconds'],
                            'ITERATIVE_HELIUM_RENDER_COMPLETE')
    preserve()
    products = {str(p.relative_to(output)): digest(p) for p in output.rglob('*')
                if p.is_file() and 'source' not in p.relative_to(output).parts and p.name != 'context.json'}
    save(output/'completed.json', dict(executionComplete=True, productionEnabled=False,
        physicsValidationComplete=False, sourceInputPreserved=True, products=products,
        completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        note='Execution completion is NOT fit convergence, isotope identification or validated PID.'))
    print('COMMON_REFIT_COMPLETE', output, summary, flush=True)


if __name__ == '__main__':
    main()
