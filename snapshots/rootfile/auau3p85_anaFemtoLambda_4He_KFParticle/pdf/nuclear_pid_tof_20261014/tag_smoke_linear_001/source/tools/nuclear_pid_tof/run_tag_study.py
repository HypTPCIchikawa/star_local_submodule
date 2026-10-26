#!/usr/bin/env python3
"""Run study-only TOF-tag shape fits from verified native, nonoverlapping bins."""
import argparse
from datetime import datetime, timedelta, timezone
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
MATH_KEYS = set('''amplitudeMaxIterations amplitudeGradientTolerance amplitudeHessianRidge
amplitudeLineSearchIterations shapeMaxCalls shapeTolerance minimumSigma maximumSigma meanAxisPadding
minimumGaussianAcceptance meanStep sigmaLogStep backgroundSlopeStep exponentialSlopeLimit boundaryTolerance
normalizeWarmStartToObserved freeBackgroundSlope displayDedxRebinFactor axisTolerance
countVarianceRelativeTolerance nativeRigidityBinWidth nativeDedxBinWidth'''.split())
TOP_KEYS = set('''schemaVersion runtimeMainconf inputStudy outputBase legacyHeader executionTimeoutSeconds
backgroundModels partitions tagComponents primaryComponents seed numeric support curves allFit'''.split())
SAMPLES = ('TAG_T', 'TAG_HE3', 'TAG_D')


def read_json(path):
    return json.loads(path.read_text(), parse_constant=lambda value:
                      (_ for _ in ()).throw(ValueError('Nonfinite JSON '+value)))


def digest(path):
    import hashlib
    h = hashlib.sha256()
    with path.open('rb') as source:
        for chunk in iter(lambda: source.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()


def fingerprint(path):
    value = path.stat()
    return dict(path=str(path), size=value.st_size, mtime_ns=value.st_mtime_ns,
                device=value.st_dev, inode=value.st_ino)


def save(path, obj):
    if path.exists():
        raise ValueError('Existing output will not be overwritten: '+str(path))
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def tsv(path, rows):
    if path.exists():
        raise ValueError('Existing TSV will not be overwritten.')
    def encode(value):
        return format(value, '.17g') if isinstance(value, float) else str(value)
    path.write_text(''.join('\t'.join(map(encode, row))+'\n' for row in rows))


def inside(path, folder):
    path = path.resolve()
    if folder.resolve() not in path.parents:
        raise ValueError('Path is outside allowed scope: '+str(path))
    return path


def probability(lo, hi, mean, sigma):
    a, b = (lo-mean)/(math.sqrt(2)*sigma), (hi-mean)/(math.sqrt(2)*sigma)
    if a >= 0:
        return .5*(math.erfc(a)-math.erfc(b))
    if b <= 0:
        return .5*(math.erfc(-b)-math.erfc(-a))
    return .5*(math.erf(b)-math.erf(a))


def nuisance_projection(candidate, curve_module, runtime_module, name, lo, hi, native_width):
    """Gaussian projection approximation, not a validated intrinsic response."""
    positions=[lo,hi,(lo+hi)/2]
    positions += [lo+(i+.5)*native_width for i in range(round((hi-lo)/native_width))]
    nodes=runtime_module.gl3_nodes(lo,hi)
    positions += [r for r,_ in nodes]
    if not all(curve_module.evaluate_curve(candidate,name,r)['supported'] for r in positions):
        return dict(supported=False,reason='Candidate lacks full selected native-bin interval support; legacy nuisance reference retained.')
    values=[(w,curve_module.evaluate_curve(candidate,name,r)) for r,w in nodes]
    mean=math.fsum(w*value['mean'] for w,value in values)
    intrinsic_variance=math.fsum(w*value['sigma']**2 for w,value in values)
    rigidity_variance=math.fsum(w*(value['mean']-mean)**2 for w,value in values)
    sigma=math.sqrt(intrinsic_variance+rigidity_variance)
    return dict(supported=True,mean=mean,sigma=sigma,intrinsicVariance=intrinsic_variance,
                rigidityMeanVariance=rigidity_variance,
                policy='TAG-derived curve GL3 Gaussian projection approximation; conditional nuisance reference, NOT truth or validation.')


def build_rows(cfg, partition, input_settings, tables, runtime_module, nuisance=None, curve_module=None):
    numeric, seed_cfg = cfg['numeric'], cfg['seed']
    spec = cfg['partitions'][partition]
    if set(spec) != {'min', 'max', 'step', 'useForIntrinsicCalibration'}:
        raise ValueError('Partition must explicitly declare min/max/step/calibration scope.')
    if type(spec['useForIntrinsicCalibration']) is not bool:
        raise ValueError('Partition calibration flag must be boolean.')
    lo, hi, step = (spec[key] for key in ('min', 'max', 'step'))
    if not all(type(v) in (int, float) and math.isfinite(v) for v in (lo, hi, step)) or not 0 < lo < hi or step <= 0:
        raise ValueError('Invalid positive R partition')
    axis_lo, axis_hi, native = (input_settings[key] for key in ('expectedXMin', 'expectedXMax', 'expectedXWidth'))
    tol = numeric['axisTolerance']
    if lo < axis_lo-tol or hi > axis_hi+tol or step < native-tol:
        raise ValueError('R partition lacks native support or is finer than native width.')
    for edge in (lo, hi):
        if abs((edge-axis_lo)/native-round((edge-axis_lo)/native))*native > tol:
            raise ValueError('R interval edges must equal native edges.')
    if abs(step/native-round(step/native))*native > tol:
        raise ValueError('R step must be an integer multiple of native bin width.')
    count = round((hi-lo)/step)
    if count <= 0 or abs(count*step-(hi-lo)) > tol:
        raise ValueError('R intervals must exactly partition the declared domain.')
    if spec['useForIntrinsicCalibration'] and abs(step-native) > tol:
        raise ValueError('Coarse projected-width fits cannot be advertised as native intrinsic calibration.')
    ylo, yhi, ywidth = (input_settings[key] for key in ('expectedYMin', 'expectedYMax', 'expectedYWidth'))
    align = seed_cfg['fitWindowAlignmentKeV']
    if abs(align/ywidth-round(align/ywidth))*ywidth > tol or align != ywidth*numeric['displayDedxRebinFactor']:
        raise ValueError('Fit-window alignment must equal exact native-to-display rebin width.')
    rows, seeds, details = [], [], []
    for sample in SAMPLES:
        components, primary = cfg['tagComponents'][sample], cfg['primaryComponents'][sample]
        if len(set(components)) != len(components) or not components or not set(primary) <= set(components):
            raise ValueError('Invalid or duplicate TAG component mapping')
        # Descending R is declared execution order only; no previous fit is
        # silently substituted as a seed in this first legacy-reference phase.
        for interval in reversed(range(count)):
            lower, upper = lo+interval*step, lo+(interval+1)*step
            center = (lower+upper)/2
            parameters = {}
            for name in components:
                if name not in runtime_module.SPECIES:
                    raise ValueError('Unknown nuclear component '+name)
                mu, sigma = runtime_module.legacy_parameters(tables, runtime_module.SPECIES.index(name), center,
                                                            ignore_domain=True)
                if not math.isfinite(mu) or not math.isfinite(sigma) or sigma <= 0:
                    raise ValueError('Legacy mathematical seed is nonfinite/nonpositive; no hidden replacement.')
                parameters[name] = (mu, sigma)
            half = seed_cfg['fitWindowSigmaHalfWidth']
            requested_lo = min(mu-half*sigma for mu,sigma in parameters.values())
            requested_hi = max(mu+half*sigma for mu,sigma in parameters.values())
            fit_lo = max(ylo, ylo+math.floor((requested_lo-ylo)/align)*align)
            fit_hi = min(yhi, ylo+math.ceil((requested_hi-ylo)/align)*align)
            if fit_hi <= fit_lo:
                raise ValueError('All seeded components lack a fit window inside the saved detector axis.')
            index = len(rows)
            rows.append([index, sample, lower, upper, fit_lo, fit_hi, ','.join(components)])
            component_details = {}
            for name in components:
                mu, sigma = parameters[name]
                replacement=None
                if nuisance is not None and sample=='TAG_D' and name not in primary:
                    replacement=nuisance_projection(nuisance,curve_module,runtime_module,name,lower,upper,native)
                    if replacement['supported']:
                        mu,sigma=replacement['mean'],replacement['sigma']
                acceptance = probability(fit_lo, fit_hi, mu, sigma)
                active = acceptance > numeric['minimumGaussianAcceptance']
                free = active and name in primary
                seeds.append([index, name, mu, sigma, seed_cfg['componentWindowWeight'],
                              int(active), int(free), int(free)])
                component_details[name] = dict(mean=mu, sigma=sigma, active=active, freeMean=free, freeSigma=free,
                    initialRelativeWindowWeight=seed_cfg['componentWindowWeight'], windowAcceptance=acceptance,
                    policy=('TAG_curve_fixed_unvalidated_projected_reference' if replacement and replacement['supported'] else
                            ('legacy_raw_log10_seed_ignore_domain' if free else 'legacy_fixed_unvalidated_nuisance_reference')),
                    nuisanceReplacement=replacement,legacyReferenceMean=parameters[name][0],legacyReferenceSigma=parameters[name][1],
                    reason=None if active else 'Gaussian seed acceptance below explicit numerical threshold',
                    identityValidated=False)
            seeds.append([index, 'BG', cfg['allFit']['backgroundStartSlope'], 1,
                          seed_cfg['backgroundWindowWeight'], 0, 0, 0])
            details.append(dict(index=index, sample=sample, lo=lower, hi=upper, r=center,
                firstNativeRBin=round((lower-axis_lo)/native)+1, lastNativeRBin=round((upper-axis_lo)/native),
                requestedFitWindow=[requested_lo,requested_hi], usedFitWindow=[fit_lo,fit_hi],
                components=component_details, widthMeaning='R-projection (folded), NOT intrinsic shared width',
                he4NuisanceValidationComplete=False if sample == 'TAG_D' else None))
    return rows, seeds, details


def validate_settings(cfg):
    if set(cfg) != TOP_KEYS or cfg['schemaVersion'] != 1:
        raise ValueError('Unknown/missing explicit TAG study settings.')
    if set(cfg['numeric']) != MATH_KEYS or set(cfg['tagComponents']) != set(SAMPLES) or set(cfg['primaryComponents']) != set(SAMPLES):
        raise ValueError('Unknown/missing TAG math or sample keys.')
    for key,value in cfg['numeric'].items():
        if isinstance(value,bool) or type(value) not in (int,float) or not math.isfinite(value):
            raise ValueError('Numeric setting must be explicit and finite: '+key)
        if key == 'meanAxisPadding':
            valid = value >= 0
        elif key in ('normalizeWarmStartToObserved','freeBackgroundSlope'):
            valid = value in (0,1)
        else:
            valid = value > 0
        if not valid:
            raise ValueError('Invalid numeric support: '+key)
        if (key.endswith('Iterations') or key.endswith('Calls') or key in ('displayDedxRebinFactor','normalizeWarmStartToObserved','freeBackgroundSlope')) and not isinstance(value,int):
            raise ValueError('Integer limit/flag required: '+key)
    seed = cfg['seed']
    if set(seed) != {'componentWindowWeight','backgroundWindowWeight','nuisanceMode','fitWindowSigmaHalfWidth','fitWindowAlignmentKeV'} or seed['nuisanceMode'] != 'legacy_fixed_reference':
        raise ValueError('First TAG phase requires the explicit unvalidated legacy nuisance reference policy.')
    for key in ('componentWindowWeight','backgroundWindowWeight','fitWindowSigmaHalfWidth','fitWindowAlignmentKeV'):
        if type(seed[key]) not in (int,float) or not math.isfinite(seed[key]) or seed[key] <= 0:
            raise ValueError('Bad warm-start/window scalar: '+key)
    if cfg['numeric']['normalizeWarmStartToObserved'] != 1:
        raise ValueError('Relative seed weights require normalization to observed fit counts.')
    if not isinstance(cfg['executionTimeoutSeconds'],int) or cfg['executionTimeoutSeconds'] <= 0:
        raise ValueError('Explicit positive timeout required.')
    if not (0 < cfg['numeric']['minimumSigma'] < cfg['numeric']['maximumSigma']):
        raise ValueError('Reversed sigma support.')
    slope=cfg['allFit']['backgroundStartSlope']
    if type(slope) not in (int,float) or not math.isfinite(slope) or abs(slope)>cfg['numeric']['exponentialSlopeLimit']:
        raise ValueError('Explicit background starting slope is outside mathematical support.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--project-dir', type=Path, required=True)
    parser.add_argument('--settings', type=Path, required=True)
    parser.add_argument('--partition', required=True, help='Explicit key in settings.partitions, e.g. smoke/pilot/native.')
    parser.add_argument('--background', choices=['linear','exponential'], required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--prepare-only', action='store_true')
    parser.add_argument('--nuisance-calibration', type=Path,
        help='Optional study-only TAG T/He3 curve JSON; replace only supported TAGD nuisance shapes as reference sensitivity.')
    args = parser.parse_args()
    project = args.project_dir.resolve(strict=True)
    if project != PROJECT.resolve(strict=True):
        raise ValueError('Explicit project-dir differs from this runner project.')
    settings_path = args.settings.resolve(strict=True)
    cfg = yaml.safe_load(settings_path.read_text());validate_settings(cfg)
    if args.partition not in cfg['partitions'] or args.background not in cfg['backgroundModels']:
        raise ValueError('Undeclared partition/background model.')
    def local(value):
        path=Path(value);path=(path if path.is_absolute() else project/path).resolve(strict=True)
        return inside(path,project)
    inputs = local(cfg['inputStudy'])
    context_path, complete_path = inputs/'context.json', inputs/'completed.json'
    if not complete_path.is_file():
        raise ValueError('Immutable INPUT stage must complete before TAG fits.')
    context, complete = read_json(context_path), read_json(complete_path)
    if (context['schemaVersion'] != 1 or context['productionEnabled'] or complete['productionEnabled'] or
            not complete['executionComplete'] or not complete['noFit'] or not complete['protectedPostRunVerified']):
        raise ValueError('Only verified immutable INPUT cache is authorized.')
    if (context['cacheKeys']!=['ALL','TAG_D','TAG_T','TAG_HE3','REMAINING'] or
            cfg['numeric']['nativeRigidityBinWidth']!=context['settings']['expectedXWidth'] or
            cfg['numeric']['nativeDedxBinWidth']!=context['settings']['expectedYWidth']):
        raise ValueError('Explicit native axes/cache labels differ from verified INPUT stage.')
    protected = {}
    def protect(path, expected=None):
        actual=digest(path)
        if expected is not None and actual!=expected:
            raise ValueError('Frozen input/source/product changed: '+str(path))
        if str(path) in protected and protected[str(path)]!=actual:
            raise ValueError('Conflicting protected provenance')
        protected[str(path)]=actual
    protect(context_path);protect(complete_path)
    for relative,expected in complete['products'].items():
        protect(inside(inputs/relative,inputs),expected)
    for relative,expected in complete['sourceSha256'].items():
        protect(inside(inputs/'source'/relative,inputs/'source'),expected)
    for path,expected in context['protectedPreRun'].items():
        protect(local(path),expected)
    original = local(context['inputFingerprint']['path']);initial=fingerprint(original)
    if initial != context['inputFingerprint']:
        raise ValueError('Merged source ROOT fingerprint changed.')
    cache = inside(inputs/'cache.root',inputs)
    mainconf = local(cfg['runtimeMainconf'])
    if mainconf != Path(context['runtimeMainconf']).resolve(strict=True):
        raise ValueError('Explicit runtime differs from immutable source stage.')
    runtime = yaml.safe_load(mainconf.read_text());info_path=local(Path('config')/runtime['analysis'])
    info=yaml.safe_load(info_path.read_text());tag=info['starTag']
    if (tag['collision'],tag['year'],tag['productionTag'],tag['libraryTag'],info['analysis']['mode']) != (
            'auau3p85','year2021','P24iy','SL24y','fxtmult'):
        raise ValueError('Only 3.85 GeV 2021 P24iy / ROOT5 SL24y is authorized.')
    runtime_module_path=HERE/'runtime_selection.py';protect(runtime_module_path)
    spec=importlib.util.spec_from_file_location('tag_runtime_selection',runtime_module_path)
    runtime_module=importlib.util.module_from_spec(spec);spec.loader.exec_module(runtime_module)
    legacy_path=local(cfg['legacyHeader']);tables=runtime_module.parse_legacy_tables(legacy_path)
    nuisance_path=None;nuisance=None;curve_module=None
    if args.nuisance_calibration:
        nuisance_path=local(args.nuisance_calibration);protect(nuisance_path);nuisance=read_json(nuisance_path)
        if nuisance.get('schemaVersion')!=1 or nuisance.get('productionEnabled') is not False or not isinstance(nuisance.get('species'),dict):
            raise ValueError('Nuisance curves must be explicit schema1 study-only candidates.')
        curve_path=HERE/'nuclear_curves.py';protect(curve_path);sys.path.insert(0,str(HERE))
        curve_spec=importlib.util.spec_from_file_location('tag_nuisance_curves',curve_path)
        curve_module=importlib.util.module_from_spec(curve_spec);curve_spec.loader.exec_module(curve_module)
    rows,seeds,seed_details=build_rows(cfg,args.partition,context['settings'],tables,runtime_module,nuisance,curve_module)
    output=args.output_dir.resolve();allowed=(project/cfg['outputBase']).resolve()
    authorized=(project/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014').resolve()
    if allowed!=authorized:
        raise ValueError('TAG outputs must remain in the authorized new nuclear_pid_tof study folder.')
    if output.exists():
        raise ValueError('Existing TAG result will not be overwritten.')
    inside(output,allowed)
    if inputs in output.parents:
        raise ValueError('TAG output must not be inside immutable inputs.')
    source_paths=[Path(__file__).resolve(),HERE/'NuclearFitCommon.h',HERE/'FitTofTagShapes.C',
                  runtime_module_path,settings_path,mainconf,info_path,legacy_path,project/'script/singularity_root_batch.sh']
    if nuisance_path:
        source_paths += [nuisance_path,HERE/'nuclear_curves.py']
    for path in source_paths:
        protect(path)
    output.mkdir(parents=True);(output/'logs').mkdir();snapshots={}
    for path in source_paths:
        relative=path.relative_to(project) if project in path.parents else Path('external')/path.name
        target=output/'source'/relative;target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
        snapshots[str(target)]=digest(path)
    engine=dict(cfg['numeric'],backgroundModel=args.background)
    tsv(output/'settings.tsv',engine.items());tsv(output/'r_slices.tsv',rows);tsv(output/'seeds.tsv',seeds)
    save(output/'seed_audit.json',dict(rows=seed_details,initializationPolicy='G/BG positive relative window weights normalized to observed count; no fit constraints or yield cap.',
        nuisancePolicy='TAGD T/He3 fixed reference only; optional supported TAG-derived shapes, otherwise legacy with explicit reason; neither validates isotope truth.',
        mathematicalLegacyTables=tables,legacyDomainIgnoredForSeedsOnly=True))
    save(output/'context.json',dict(schemaVersion=1,productionEnabled=False,physicsValidationComplete=False,
        inputFingerprint=initial,inputStudy=str(inputs),inputCompletionSha256=digest(complete_path),cacheRoot=str(cache),
        runtimeMainconf=str(mainconf),settings=cfg,selectedPartition=args.partition,background=args.background,
        fitOrder='Sample T/He3/D; descending R per sample; all initial shapes remain explicit legacy seeds.',
        seedPolicy='Raw log10 tables, PMomentumRange ignored ONLY for initial guesses; old working knots not used.',
        nuisancePolicy='TAGD T/He3 are fixed unvalidated conditional references; optional TAG curve replacement/fallback is audited per row; He4 scientific support is not established.',
        nuisanceCalibration=str(nuisance_path) if nuisance_path else None,
        nuisanceCalibrationSha256=digest(nuisance_path) if nuisance_path else None,
        foldedWidthMustNotBeAdopted=True,
        eligibleForSubsequentIntrinsicCalibration=cfg['partitions'][args.partition]['useForIntrinsicCalibration'],
        sourceNativeWidth=context['settings']['expectedXWidth'],rawDedxFitWidth=context['settings']['expectedYWidth'],
        fitUsesOnlySavedTagSpectraOnce=True,protected=protected,snapshots=snapshots,
        stageConsumedBlocks=['partitions','tagComponents','primaryComponents','seed','numeric','allFit.backgroundStartSlope'],
        deferredToLaterStages=['support','curves','allFit selector/purity controls'],
        createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))
    def preserve():
        if fingerprint(original)!=initial:
            raise ValueError('Original merged ROOT changed during TAG study.')
        for path,expected in dict(protected,**snapshots).items():
            if digest(Path(path))!=expected:
                raise ValueError('Protected input/source changed during TAG fit: '+path)
    preserve()
    macro=output/'source/tools/nuclear_pid_tof/FitTofTagShapes.C'
    call=str(macro)+'+('+','.join(json.dumps(str(p)) for p in
         (cache,output,output/'settings.tsv',output/'r_slices.tsv',output/'seeds.tsv'))+')'
    command=['bash',str(project/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',call]
    save(output/'command.json',command)
    print('NUCLEAR_TAG_PREPARED' if args.prepare_only else 'NUCLEAR_TAG_START',output,'rows='+str(len(rows)),flush=True)
    if args.prepare_only:
        return
    with (output/'logs/fit.log').open('w') as log:
        log.write(json.dumps(command)+'\n');log.flush()
        result=subprocess.run(command,cwd=project,stdout=log,stderr=subprocess.STDOUT,timeout=cfg['executionTimeoutSeconds'])
    content=(output/'logs/fit.log').read_text(errors='replace')
    if result.returncode or 'NUCLEAR_TAG_FIT_COMPLETE' not in content or any(marker in content for marker in
            ['Error in <ACLiC>','segmentation violation','fatal error:','Error: Function']):
        raise RuntimeError('TAG fit execution failed; preserve new log/output for diagnosis.')
    fitted=read_json(output/'tag_fits.json')
    if fitted['productionEnabled'] or fitted['fitStage']!='independent_tag_shapes' or len(fitted['slices'])!=len(rows):
        raise ValueError('Unexpected TAG result scope.')
    for observed,expected in zip(fitted['slices'],rows):
        for key,value in zip(('index','sample','lo','hi','fitLo','fitHi'),expected[:6]):
            if observed[key]!=value:
                raise ValueError('Fit result changed explicit data partition/window: '+key)
    preserve()
    summary=dict(executionComplete=True,productionEnabled=False,physicsValidationComplete=False,
        partition=args.partition,background=args.background,nSlices=len(rows),totalDataDeviance=fitted['totalDataDeviance'],
        numericalConverged=fitted['numericalConverged'],numericallyConvergedSlices=sum(r['numericalConverged'] for r in fitted['slices']),
        he4IndependentTagNuisanceValidationComplete=False,
        note='Execution completion is not convergence, isotope identification, supported calibration or validated purity.')
    save(output/'summary.json',summary)
    products={str(p.relative_to(output)):digest(p) for p in output.rglob('*')
              if p.is_file() and 'source' not in p.relative_to(output).parts and p.name!='context.json'}
    save(output/'completed.json',dict(summary,inputFingerprintPreserved=True,protectedSourcesPreserved=True,
        products=products,completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))
    print('NUCLEAR_TAG_COMPLETE',output,summary,flush=True)


if __name__=='__main__':
    main()
