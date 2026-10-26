#!/usr/bin/env python3
"""Immutable high-R TOF-conditional study; never modifies low-R products."""
import argparse
from datetime import datetime, timezone
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
import yaml

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import run_cutoff as common
from summarize_cutoff import aggregate_one
PROJECT = common.PROJECT
SPECIES = common.SPECIES
MAPPING = {'d': ('TAG_D', 'D'), 't': ('TAG_T', 'T'), 'he3': ('TAG_HE3', 'HE3'), 'he4': ('TAG_D', 'D')}


def save(path, value):
    if path.exists():
        raise ValueError('Refusing overwrite: '+str(path))
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def gate_parameters(runtime, provider, species, rigidity):
    mean, sigma, policy = provider.current_parameters(runtime, species, rigidity)
    # Competing hypotheses remain active before target-only TOF/domain rejection.
    return mean, sigma, bool(rigidity > 0 and mean > 0 and sigma > 0), policy


def require_exact_window(target, mapping, current, context):
    if target not in MAPPING or mapping['sample'] != MAPPING[target][0]:
        raise ValueError('Explicit saved sample/target mapping required.')
    high = current['cuts']['highRigidityTof']
    if high['enabled'] is not True or high['thresholdPOverQ'] != 2.0 or current['cuts']['m2_selection']:
        raise ValueError('Require explicit all-nuclei post-winner high-R profile.')
    suffix = MAPPING[target][1]
    expected = tuple(context['settings']['tag'+suffix+key] for key in ('Mass2Min', 'Mass2Max'))
    actual = high['species'][target]
    matches=all(abs(a-b)<=1.e-12 for a,b in zip(expected,(actual['mass2Min'],actual['mass2Max'])))
    if mapping['exactProductionTofWindow'] is not matches:
        raise ValueError('Saved and requested production TOF windows differ.')
    if target=='he4' and mapping.get('studyTofPolicy')!='same_saved_D_window':
        raise ValueError('He4 study requires authorized same-D window.')
    return expected


def selected_bins(cfg, stage):
    width=cfg['nativeWidth']
    first,last=round(cfg['rigidityMin']/width)+1,round(cfg['rigidityMax']/width)
    if first < 1 or last < first or any(abs(round(cfg[k]/width)*width-cfg[k])>1.e-9 for k in ('rigidityMin','rigidityMax')):
        raise ValueError('R edges must equal native edges.')
    if stage=='native': return list(range(first,last+1))
    result=[round(r/width-.5)+1 for r in cfg['pilotRigidity']]
    if len(result)!=len(set(result)) or any(not first<=b<=last or abs((b-.5)*width-r)>1.e-9 for b,r in zip(result,cfg['pilotRigidity'])):
        raise ValueError('Pilot must identify distinct native bin centers.')
    return result


def checked_settings(path):
    cfg=yaml.safe_load(path.read_text())
    required={'schemaVersion','runtimeMainconf','inputStudy','outputBase','executionTimeoutSeconds','nativeWidth','aggregateWidth',
        'rigidityMin','rigidityMax','pilotRigidity','targets','tagStudies','responseModels','backgroundModels','shapeBackgroundSource',
        'shapeVariations','meanShiftSigma','widthRelativeShift','numeric','allFit','support','cutoff'}
    if set(cfg)!=required or cfg['schemaVersion']!=1 or cfg['nativeWidth']!=.01:
        raise ValueError('Missing/unknown explicit high-R settings.')
    if set(cfg['targets'])!=set(MAPPING) or cfg['responseModels']!=['gaussian','lognormal'] or cfg['backgroundModels']!=['linear','exponential'] or cfg['shapeBackgroundSource']!='linear':
        raise ValueError('Explicit four-model linear-reference study required.')
    if cfg['allFit']['targetShapeMode']!='free_target_mean_sigma' or not cfg['allFit']['requireBestSpecies']:
        raise ValueError('Free target and original best-species gate required.')
    if cfg['rigidityMin']!=2.0 or cfg['rigidityMax']>5.0:
        raise ValueError('Study exceeds explicit response provider domain.')
    return cfg


def prepare(args):
    settings=common.resolve(args.settings);cfg=checked_settings(settings)
    if args.variation not in cfg['shapeVariations']: raise ValueError('Undeclared variation.')
    mainconf=common.resolve(cfg['runtimeMainconf'])
    provider=common.module('high_r_runtime',PROJECT/'tools/nuclear_pid_tof/runtime_selection.py')
    current=provider.load_current_selection(PROJECT,mainconf)
    main=yaml.safe_load(mainconf.read_text());info_path=PROJECT/'config'/main['analysis']
    info=yaml.safe_load(info_path.read_text())['starTag']
    if tuple(info[k] for k in ('collision','year','productionTag','libraryTag'))!=('auau3p85','year2021','P24iy','SL24y'):
        raise ValueError('Runtime must be 3.85 GeV 2021 P24iy SL24y.')
    if cfg['allFit']['maxNSigma']!=current['cuts']['maxNSigmaNuclear']:
        raise ValueError('Study and current TPC gate/domain differ.')
    cfg=dict(cfg,rigidityMax=min(cfg['rigidityMax'],5.0 if args.target=='he4' else current['cuts']['maxPOverQ']))
    if cfg['rigidityMin']<current['cuts']['speciesMinPOverQ'][args.target] or any(v['min']<cfg['rigidityMax'] and v['max']>cfg['rigidityMin'] for v in current['cuts']['speciesPOverQVeto'][args.target]):
        raise ValueError('Target min/veto overlaps evaluated domain.')
    inputs=common.resolve(cfg['inputStudy'])
    done,context,audit=(common.read(inputs/n) for n in ('completed.json','context.json','native_audit.json'))
    if not done.get('executionComplete') or not done.get('protectedPostRunVerified') or not all(audit.get(k) for k in ('axesExactlyEqual','cacheCountErrorFlowPreserved','writtenCacheReadbackExact')):
        raise ValueError('Input extraction/readback incomplete.')
    cache=inputs/'cache.root'
    if common.digest(cache)!=done['products']['cache.root']: raise ValueError('Native cache changed.')
    for name in ('context.json','native_audit.json'):
        if common.digest(inputs/name)!=done['products'][name]:
            raise ValueError('Changed original extraction metadata: '+name)
    original=Path(context['inputFingerprint']['path'])
    if '2021P24iy_partial_merge' not in original.name or common.fingerprint(original)!=context['inputFingerprint']:
        raise ValueError('Original 2021 P24iy ROOT changed.')
    window=require_exact_window(args.target,cfg['targets'][args.target],current,context)
    axis=audit['nativeAxis']
    if axis['rBinWidth']!=cfg['nativeWidth'] or axis['dedxBinWidth']!=cfg['numeric']['nativeDedxBinWidth']:
        raise ValueError('Native geometry differs.')
    refs,protected={},dict(current['protectedFiles'])
    for response in cfg['responseModels']:
        for background in cfg['backgroundModels']:
            rows,signatures=common.tag_index(common.resolve(cfg['tagStudies'][response][background]))
            refs[response,background]=rows;protected.update(signatures)
    for p in (settings,info_path,cache,inputs/'context.json',inputs/'completed.json',inputs/'native_audit.json',PROJECT/'StMaker/common/FemtoLambdaLegacy.cxx',PROJECT/'StMaker/common/FemtoLambdaLegacy.h'):
        protected[str(p)]=common.digest(p)
    shapes,gates,diagnostics,windows=[],[],[],[]
    for index,rbin in enumerate(selected_bins(cfg,args.stage)):
        lo,hi=(rbin-1)*cfg['nativeWidth'],rbin*cfg['nativeWidth']
        chosen,failures,source_rows,references={},{},{},[]
        for s,species in enumerate(SPECIES):
            sample=dict(d='TAG_D',t='TAG_T',he3='TAG_HE3',he4='TAG_D')[species]
            row=refs[args.response,'linear'].get((sample,rbin))
            if row is not None:
                component=row['components'][species]
                if not common.positive(component['mean']) or not common.positive(component['sigma']):
                    raise ValueError('No invented fallback for invalid response seed.')
                chosen[species]=component
                failures[species]=common.shape_failures(row,component,cfg['support'])
                bounds=common.named_component_boundaries(row,species,args.response,axis)
                if any(b['atBoundary'] for b in bounds): failures[species].append('named_response_math_boundary')
                source_rows[species]=dict(sample=sample,rowIndex=row['index'],numericalConverged=row['numericalConverged'],
                    pearson=row['chi2'],ndf=row['conditionalNdf'],parameterGuards=bounds,
                    referencePolicy='Conditional response seed; never an independent prior.')
            else:
                mean,sigma,_=provider.current_parameters(current,s,.5*(lo+hi))
                if not common.positive(mean) or not common.positive(sigma): raise ValueError('Provider seed outside physical domain.')
                chosen[species]=dict(mean=mean,sigma=sigma)
                failures[species]=[] if species==args.target else ['provider_only_nuisance_response_unvalidated']
                source_rows[species]=dict(sample=sample,rowIndex=None,numericalConverged=None,
                    referencePolicy='No TAG fit beyond R=2.5; provider initializes only, not a validated TAG extrapolation.')
            if species==args.target:
                for response in cfg['responseModels']:
                    other=refs[response,'linear'].get((sample,rbin))
                    c=other['components'][species] if other is not None else chosen[species]
                    if not common.positive(c['mean']) or not common.positive(c['sigma']): raise ValueError('Invalid target model seed.')
                    references.append((c['mean'],c['sigma']))
        for node,(r,weight) in enumerate(provider.gl3_nodes(lo,hi)):
            for s,species in enumerate(SPECIES):
                mean,sigma,active,_=gate_parameters(current,provider,s,r)
                gates.append([index,node,species,mean,sigma,int(active)])
                if species==args.target and active: references.append((mean,sigma))
        half,align=(cfg['allFit'][k] for k in ('fitWindowSigmaHalfWidth','fitWindowAlignmentKeV'))
        ylo,yhi=axis['dedxRange']
        fitlo=max(ylo,math.floor(min(m-half*s for m,s in references)/align)*align)
        fithi=min(yhi,math.ceil(max(m+half*s for m,s in references)/align)*align)
        if not fithi>fitlo: raise ValueError('Target lacks observable fit window.')
        windows.append(dict(index=index,rBin=rbin,lo=lo,hi=hi,fitLo=fitlo,fitHi=fithi))
        for node,(r,weight) in enumerate(provider.gl3_nodes(lo,hi)):
            for species in SPECIES:
                mean,sigma=common.varied_response(chosen[species],args.variation,cfg,species,args.target)
                shapes.append([index,rbin,lo,hi,node,weight,r,species,mean,sigma,int(not failures[species]),fitlo,fithi])
        diagnostics.append(dict(index=index,rBin=rbin,lo=lo,hi=hi,failures=failures,tagComponents=chosen,tagDiagnostics=source_rows))
    output=args.output_dir.resolve();allowed=(PROJECT/cfg['outputBase']).resolve()
    if output.exists() or output.parent!=allowed: raise ValueError('Require new direct child of outputBase.')
    output.mkdir(parents=True)
    source=output/'source'
    files=[Path(__file__).resolve(),HERE/'FitHighRigidity.C',HERE/'SelfTestHighRigidity.C',HERE/'test_highrigidity.py',
        settings,HERE.parent/'FitLocalNuclearCutoff.C',HERE.parent/'run_cutoff.py',HERE.parent/'summarize_cutoff.py',
        PROJECT/'tools/nuclear_pid_tof/runtime_selection.py',PROJECT/'tools/nuclear_pid_tof/NuclearFitCommon.h',
        PROJECT/'tools/nuclear_pid_tof/NuclearLognormalFitCommon.h',PROJECT/'script/singularity_root_batch.sh']
    for path in files:
        dest=source/path.relative_to(PROJECT);dest.parent.mkdir(parents=True,exist_ok=True)
        expected=common.digest(path);shutil.copy2(path,dest)
        if common.digest(dest)!=expected: raise ValueError('Source changed during snapshot.')
        protected[str(path)]=expected;protected[str(dest)]=expected
    for name in ('NuclearFitCommon.h','NuclearLognormalFitCommon.h'):
        dest=source/'tools/nuclear_pid_cutoff'/name;src=PROJECT/'tools/nuclear_pid_tof'/name
        shutil.copy2(src,dest);protected[str(dest)]=common.digest(src)
    snapshot=source/'tools/nuclear_pid_cutoff/highrigidity'
    common.tsv(output/'shapes.tsv','index rBin rLo rHi node weight rValue species mean sigma supported fitLo fitHi',shapes)
    common.tsv(output/'gates.tsv','index node species mean sigma active',gates)
    exact=cfg['targets'][args.target]['exactProductionTofWindow']
    numeric=dict(cfg['numeric'],responseModel=args.response,backgroundModel=args.background,
        sampleHistogram=cfg['targets'][args.target]['sample'],targetSpecies=args.target,
        mass2Min=window[0],mass2Max=window[1],exactProductionTofWindow=int(exact),
        rigidityMin=cfg['rigidityMin'],rigidityMax=cfg['rigidityMax'],tofThreshold=current['cuts']['highRigidityTof']['thresholdPOverQ'])
    for key,value in cfg['allFit'].items():
        if key not in ('fitWindowSigmaHalfWidth','fitWindowAlignmentKeV','requireBestSpecies'): numeric[key]=value
    common.tsv(output/'settings.tsv','',sorted(numeric.items()))
    save(output/'effective_settings.json',cfg);save(output/'fixed_windows.json',windows)
    save(output/'tag_shape_audit.json',diagnostics);save(output/'current_selection.json',current)
    save(output/'context.json',dict(schemaVersion=1,productionEnabled=False,createdUtc=datetime.now(timezone.utc).isoformat(),
        stage=args.stage,targetSpecies=args.target,sampleHistogram=numeric['sampleHistogram'],response=args.response,background=args.background,variation=args.variation,
        exactProductionTofWindow=exact,mass2Window=list(window),mass2Bounds='strict',runtimeMainconf=str(mainconf),
        sourceOriginalInput=context['inputFingerprint'],cache=str(cache),protectedSources=protected,
        historicalInputSelection='Saved pre-TPC charge-inclusive QA; independent of old final nuclear PID/domain selection.',
        evaluatedSelection='Original 4-way TPC winner, then target TOF and domain. No reassignment. He4 study uses saved D window, not production 3-4.',
        likelihoodPolicy='Selected TOF sample only; same-data TAG parameters are seeds/references, not independent likelihood or priors.',
        widthMeaning='Gaussian arithmetic mean/std versus lognormal mode/curvature; native projected response, not intrinsic.',
        limitations=['Inclusive track purity, not CF-pair-conditioned purity or true efficiency.',
            'He4 production 3<m2<4 remains unassessed; study uses saved D window.',
            'Sensitivity envelope is not a confidence interval or full systematic uncertainty.',
            'Unassessed points are never interpolated.']))
    return cfg,current,mainconf,cache,original,context,protected,output,snapshot,diagnostics


def summarize_run(cfg,output,target,audit):
    fit=common.read(output/'cutoff_fits.json');support_cfg=dict(cfg,fitWindowMode='target_local')
    audits={r['index']:r for r in audit};results=[]
    for row in fit['slices']:
        value=aggregate_one([row],target,audits,support_cfg)
        component=row['components'][target]
        if component.get('conditionalDropDeltaD') is None or component['conditionalDropDeltaD']<cfg['support']['minDropDeltaD']:
            value['reasons']=sorted(set(value['reasons']+['fitted_drop_delta_d_too_small']));value['supported']=False
        results.append(dict(species=target,rBin=row['rBin'],rLo=row['lo'],rHi=row['hi'],sample=fit['sampleHistogram'],
            signalSelected=value['signalSelected'],totalSelected=value['totalSelected'],rawGateCount=value['rawSelected'],
            purityNumerical=value['purity'],quality=value['supported'],reasons=value['reasons'],
            localResidual=value['localResidualFrac'],integralResidual=value['integralResidualFrac'],
            conditionalError=value['conditionalErrorWeightedBound'],numericalValid=row['numericalValid'],
            mean=component['mean'],sigma=component['sigma'],meanError=component['meanError'],sigmaError=component['sigmaError']))
    save(output/'native_summary.json',dict(schemaVersion=1,productionEnabled=False,sample=fit['sampleHistogram'],
        targetSpecies=target,response=fit['responseModel'],background=fit['backgroundModel'],nativeRows=results))
    keys=('species','rBin','rLo','rHi','sample','signalSelected','totalSelected','rawGateCount','purityNumerical','quality',
          'localResidual','integralResidual','conditionalError','mean','sigma','meanError','sigmaError')
    common.tsv(output/'native_summary.tsv',' '.join(keys)+' reasons',
        [[r[k] if r[k] is not None else 'nan' for k in keys]+[';'.join(r['reasons']) or 'supported'] for r in results])
    return results


def run(args):
    cfg,current,mainconf,cache,original,input_context,protected,output,snapshot,audit=prepare(args)
    expression=str(snapshot/('SelfTestHighRigidity.C' if args.self_test_only else 'FitHighRigidity.C'))+'+('
    if not args.self_test_only:
        expression+=','.join(json.dumps(str(p)) for p in (cache,output,output/'settings.tsv',output/'shapes.tsv',output/'gates.tsv'))
    expression+=')'
    command=['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',expression]
    save(output/'command.json',command)
    if args.prepare_only:
        save(output/'prepared.json',dict(executionComplete=False,prepared=True));print('NUCLEAR_HIGH_R_PREPARED',output,flush=True);return
    started=time.monotonic()
    try:
        with (output/'fit.log').open('w') as log:
            result=subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=cfg['executionTimeoutSeconds'])
        marker='NUCLEAR_HIGH_R_SELFTEST_PASS' if args.self_test_only else 'NUCLEAR_HIGH_R_COMPLETE'
        if result.returncode or marker not in (output/'fit.log').read_text():
            raise RuntimeError('ROOT failed/completion absent: '+str(output/'fit.log'))
        summary=[] if args.self_test_only else summarize_run(cfg,output,args.target,audit)
        for path,expected in protected.items():
            if common.digest(path)!=expected: raise ValueError('Protected source/cache changed during fit: '+path)
        if common.fingerprint(original)!=input_context['inputFingerprint']: raise ValueError('Original 2021 P24iy ROOT changed.')
        products={str(p.relative_to(output)):common.digest(p) for p in output.rglob('*') if p.is_file() and p.suffix not in ('.so','.d','.pcm')}
        save(output/'completed.json',dict(executionComplete=True,productionEnabled=False,fullUncertaintyValidated=False,
            protectedSourcesVerified=True,originalInputPreserved=True,runtimeSeconds=time.monotonic()-started,
            nativeRows=len(summary),supportedRows=sum(r['quality'] for r in summary),targetSpecies=args.target,
            stage=args.stage,response=args.response,background=args.background,variation=args.variation,
            selfTestOnly=args.self_test_only,products=products))
        print('NUCLEAR_HIGH_R_RUN_COMPLETE',output,'supported',sum(r['quality'] for r in summary),'/',len(summary),flush=True)
    except Exception as error:
        save(output/'failed.json',dict(executionComplete=False,error=str(error),runtimeSeconds=time.monotonic()-started));raise


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--target',required=True,choices=tuple(MAPPING))
    parser.add_argument('--stage',required=True,choices=('pilot','native'))
    parser.add_argument('--response',required=True,choices=('gaussian','lognormal'))
    parser.add_argument('--background',required=True,choices=('linear','exponential'))
    parser.add_argument('--variation',required=True)
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--self-test-only',action='store_true')
    run(parser.parse_args())


if __name__=='__main__': main()
