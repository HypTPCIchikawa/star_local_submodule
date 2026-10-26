#!/usr/bin/env python3
"""Display-only raw COLZ + saved low-domain polynomial fits; never refit/adopt PID."""
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


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--latest-study', type=Path, required=True)
    parser.add_argument('--reference-study', type=Path, required=True)
    parser.add_argument('--qr-source', type=Path, required=True)
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--settings', type=Path, default=HERE/'low_domain_overlay_settings.yaml')
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    import run_common_refit as guards
    digest, read, save, inside = guards.digest, guards.read_json, guards.save, guards.inside
    latest, reference = args.latest_study.resolve(strict=True), args.reference_study.resolve(strict=True)
    protected = {}

    def protect(path, expected=None):
        actual = digest(path)
        if expected is not None and actual != expected:
            raise ValueError('Protected source/product differs: '+str(path))
        if str(path) in protected and protected[str(path)] != actual:
            raise ValueError('Conflicting protected provenance')
        protected[str(path)] = actual

    def checked_study(folder):
        context_path, complete_path = folder/'context.json', folder/'completed.json'
        if not complete_path.is_file():
            raise ValueError('Wait for source fit/render/hash completion before overlay; no running study is consumed.')
        context, complete = read(context_path), read(complete_path)
        if (context.get('indices') != list(range(30)) or context.get('studyDomain') != [.4,.7] or
                context['productionEnabled'] or not complete['executionComplete'] or complete['productionEnabled']):
            raise ValueError('A complete study-only 30-slice low-domain result is required.')
        protect(context_path);protect(complete_path)
        for encoded, expected in dict(context['protected'], **context['snapshots']).items():
            path = Path(encoded)
            path = path.resolve(strict=True) if path.is_absolute() else (PROJECT/path).resolve(strict=True)
            inside(path, PROJECT);protect(path, expected)
        for relative, expected in complete['products'].items():
            protect(guards.checked_product(folder, relative), expected)
        result = read(folder/'fit_result.json')
        if result['phase'] != 'common' or [row['index'] for row in result['slices']] != list(range(30)):
            raise ValueError('Unexpected result scope')
        guards.require_shared_seed(result['sharedCoefficients'], context['settings']['engine'])
        return context, result

    context, current = checked_study(latest)
    ref_context, ref_result = checked_study(reference)
    if context['input'] != ref_context['input'] or context['sourceStudy'] != ref_context['sourceStudy']:
        raise ValueError('Reference and latest use different source data.')
    cfg = copy.deepcopy(context['settings'])
    sample = cfg['sample']
    if (sample['collision'],sample['year'],sample['production']) != ('auau3p85',2021,'P24iy'):
        raise ValueError('Only 3.85 GeV 2021 P24iy is authorized.')
    allowed = (PROJECT/sample['outputBase']).resolve()
    inside(latest,allowed);inside(reference,allowed)
    helper = HERE/'run_iterative_fit.py'
    if protected.get(str(helper)) != digest(helper):
        raise ValueError('Runtime helper differs from the frozen study.')
    spec = importlib.util.spec_from_file_location('overlay_frozen_helpers',helper)
    helpers = importlib.util.module_from_spec(spec);spec.loader.exec_module(helpers)
    source = (PROJECT/sample['input']).resolve(strict=True)
    initial = helpers.fingerprint(source)
    if initial != context['input']:
        raise ValueError('Original input fingerprint changed')
    mainconf = args.mainconf.resolve(strict=True)
    if mainconf != (PROJECT/sample['mainconf']).resolve(strict=True):
        raise ValueError('Explicit mainconf differs from the study.')
    runtime = yaml.safe_load(mainconf.read_text())
    analysis_path = (PROJECT/'config'/runtime['analysis']).resolve(strict=True)
    info = yaml.safe_load(analysis_path.read_text());tag = info['starTag']
    if (tag['collision'],tag['year'],tag['productionTag'],tag['libraryTag'],info['analysis']['mode']) != (
            'auau3p85','year2021','P24iy','SL24y','fxtmult'):
        raise ValueError('Explicit runtime dataset/library mismatch')
    qr_path, settings_path = args.qr_source.resolve(strict=True), args.settings.resolve(strict=True)
    if (qr_path != Path(context['coefficientSource']).resolve(strict=True) or
            digest(qr_path) != context['coefficientSourceSha256']):
        raise ValueError('QR source differs from the actual saved fit initialization.')
    qr = read(qr_path)
    if qr['domain'] != [.4,.7]:
        raise ValueError('QR curve domain differs.')
    qr_coefficients = []
    for sp in ('he3','he4'):
        if qr['coefficients'][sp]['nSlices'] != 30:
            raise ValueError('QR source slice count differs')
        qr_coefficients += qr['coefficients'][sp]['meanCoefficients']+qr['coefficients'][sp]['sigmaCoefficients']
    guards.require_shared_seed(qr_coefficients,cfg['engine'])
    display = yaml.safe_load(settings_path.read_text())
    required = {'schemaVersion','rawHistogram','domain','displayRigidity','displayDedx','nuclearDedx',
        'dedxRebinFactor','curveStep','curveRelativeTolerance','axisTolerance','logZ','minimumZ','canvasWidth','canvasHeight'}
    if (set(display) != required or display['schemaVersion'] != 1 or display['rawHistogram'] != 'true/hDedxP' or
            display['domain'] != [.4,.7] or display['displayRigidity'] != [.35,.8] or
            display['displayDedx'] != [0.,120.] or display['nuclearDedx'] != [20.,120.]):
        raise ValueError('Display-only scope/axes must be explicitly declared.')
    if type(display['logZ']) is not bool:
        raise ValueError('Explicit logZ boolean required')
    for key in ('dedxRebinFactor','canvasWidth','canvasHeight'):
        if type(display[key]) is not int or display[key] <= 0:
            raise ValueError('Invalid positive integer display setting')
    for key in ('curveStep','curveRelativeTolerance','axisTolerance','minimumZ'):
        if type(display[key]) not in (int,float) or not math.isfinite(display[key]) or display[key] <= 0:
            raise ValueError('Invalid positive numeric display setting')
    if display['curveStep'] > .001 or display['curveStep'] < .00001:
        raise ValueError('Unreasonable display grid')
    curves = [current['sharedCoefficients'],ref_result['sharedCoefficients'],qr_coefficients]
    minima = [guards.require_shared_seed(c,cfg['engine']) for c in curves]
    grid = []
    count = math.ceil((.7-.4)/display['curveStep'])
    for mode,c in enumerate(curves):
        for sp in range(2):
            for i in range(count+1):
                r = .7 if i == count else .4+(.7-.4)*i/count
                z = math.log10(r)
                mu = math.fsum(c[7*sp+j]*z**j for j in range(4))
                sigma = math.fsum(c[7*sp+4+j]*z**j for j in range(3))
                if not math.isfinite(mu) or not math.isfinite(sigma) or sigma <= 0:
                    raise ValueError('Nonfinite/unsupported curve point; no fallback is permitted.')
                grid.append([mode,sp,r,mu,sigma])
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError('Existing output must not be overwritten.')
    inside(output,allowed)
    if latest in output.parents or reference in output.parents:
        raise ValueError('Overlay must be a new sibling, not inside a source study.')
    macro = HERE/'DrawLowDomainRawOverlay.C'
    wrapper = PROJECT/'script/singularity_root_batch.sh'
    for path in (Path(__file__).resolve(),HERE/'run_common_refit.py',macro,settings_path,mainconf,analysis_path,wrapper):
        protect(path)
    output.mkdir(parents=True)
    (output/'source').mkdir()
    snapshots = {}
    for path in (Path(__file__).resolve(),macro,settings_path,mainconf,analysis_path,wrapper):
        relative = path.relative_to(PROJECT)
        target = output/'source'/relative
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
        snapshots[str(target)] = digest(path)
    root_settings = {k:v for k,v in display.items() if not isinstance(v,list)}
    root_settings['logZ'] = int(display['logZ'])
    for src,prefix in [('domain','domain'),('displayRigidity','displayRigidity'),
                       ('displayDedx','displayDedx'),('nuclearDedx','nuclearDedx')]:
        root_settings[prefix+'Min'],root_settings[prefix+'Max'] = display[src]
    root_settings.update(latestConverged=int(current['numericalConverged']),latestGlobalStatus=current['common']['status'])
    helpers.tsv(output/'settings.tsv',root_settings.items())
    helpers.tsv(output/'coefficients.tsv',[[label]+c for label,c in zip(('latest','outer12','initialQR'),curves)])
    helpers.tsv(output/'curve_grid.tsv',grid)
    save(output/'context.json',dict(schemaVersion=1,displayOnly=True,fitPerformed=False,productionEnabled=False,
        sourceInput=initial,rawHistogram=display['rawHistogram'],sourceLatest=str(latest),
        sourceReference=str(reference),sourceQr=str(qr_path),settings=display,protected=protected,snapshots=snapshots,
        domain=[.4,.7],continuousSigmaMinima=minima,widthMeaning='Latest/outer12 intrinsic shared kernel; initial QR folded projection-width seed.',
        limitations=['Charge-inclusive pre-PID raw tracks, not independently identified nuclei.',
            'Assigned He3/He4 are model-dependent components; no identity/purity validation.',
            'Curves are not drawn outside .4-.7; no old knot/fallback/extrapolation.',
            'Execution completion and numerical convergence are separate; latest result not adopted.'],
        latestNumericalConverged=current['numericalConverged'],latestCommon=current['common'],
        createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))

    def preserve():
        if helpers.fingerprint(source) != initial:
            raise ValueError('Original ROOT fingerprint changed')
        for path,expected in dict(protected,**snapshots).items():
            if digest(Path(path)) != expected:
                raise ValueError('Original product/source changed during display: '+path)

    preserve()
    frozen_macro = output/'source'/macro.relative_to(PROJECT)
    call = str(frozen_macro)+'+('+','.join(json.dumps(str(p)) for p in (
        source,output,output/'settings.tsv',output/'coefficients.tsv',output/'curve_grid.tsv'))+')'
    command = ['bash',str(wrapper),str(mainconf),'root4star','-l','-b','-q',call]
    save(output/'command.json',command)
    print('LOW_DOMAIN_RAW_OVERLAY_PREPARED' if args.prepare_only else 'LOW_DOMAIN_RAW_OVERLAY_START',output,flush=True)
    if args.prepare_only:
        return
    helpers.run_command(command,output/'root.log',1200,'LOW_DOMAIN_RAW_OVERLAY_COMPLETE')
    audit = read(output/'overlay_audit.json')
    if (audit['pages'] != 4 or audit['rawCountsBeforeDisplayRebin'] != audit['rawCountsAfterDisplayRebin'] or
            audit['outsideDomainDrawn'] or audit['fitPerformed'] or audit['pidAdopted']):
        raise ValueError('Display-only count/domain audit failed')
    preserve()
    products = {str(p.relative_to(output)):digest(p) for p in output.rglob('*')
                if p.is_file() and 'source' not in p.relative_to(output).parts and p.name != 'context.json'}
    save(output/'completed.json',dict(executionComplete=True,displayOnly=True,productionEnabled=False,
        originalsPreserved=True,products=products,audit=audit,
        completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))
    print('LOW_DOMAIN_RAW_OVERLAY_COMPLETE',output,flush=True)


if __name__ == '__main__':
    main()
