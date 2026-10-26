#!/usr/bin/env python3
"""Independent 0.05-R projected He4 fits; never averages/replaces native purity."""
import argparse
from datetime import datetime,timezone
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys
import time
import yaml
import run_highrigidity as high
common=high.common


def prepare(args):
    settings=common.resolve(args.settings);cfg=yaml.safe_load(settings.read_text())
    required={'schemaVersion','baseSettings','inputStudy','outputBase','targetSpecies','sampleHistogram','tofPolicy',
        'rigidityMin','rigidityMax','fitRigidityWidth','sourceNativeWidth','pilotRigidity','targetResponseProjection'}
    if set(cfg)!=required or cfg['schemaVersion']!=1 or (cfg['targetSpecies'],cfg['sampleHistogram'],cfg['tofPolicy'])!=('he4','TAG_D','same_saved_D_window'):
        raise ValueError('Explicit separate He4 same-D-window coarse settings required.')
    if cfg['targetResponseProjection']!='constant_projected_shape_free_mean_sigma' or cfg['fitRigidityWidth']!=.05 or cfg['sourceNativeWidth']!=.01:
        raise ValueError('Unsupported explicit coarse projection schema.')
    base_path=common.resolve(cfg['baseSettings']);base=high.checked_settings(base_path)
    mainconf=common.resolve(base['runtimeMainconf'])
    provider=common.module('coarse_runtime',common.PROJECT/'tools/nuclear_pid_tof/runtime_selection.py')
    current=provider.load_current_selection(common.PROJECT,mainconf)
    if base['allFit']['maxNSigma']!=current['cuts']['maxNSigmaNuclear'] or cfg['rigidityMin']<current['cuts']['speciesMinPOverQ']['he4']:
        raise ValueError('Configured TPC gate/minimum differs from the coarse study.')
    if any(v['min']<cfg['rigidityMax'] and v['max']>cfg['rigidityMin'] for v in current['cuts']['speciesPOverQVeto']['he4']):
        raise ValueError('Target post-winner veto overlaps the coarse study domain.')
    if not current['cuts']['highRigidityTof']['extendHe4BeyondMaxPOverQ'] or cfg['rigidityMax']>current['tables']['domain'][3][1]:
        raise ValueError('He4 study extension exceeds explicitly enabled provider domain.')
    inputs=common.resolve(cfg['inputStudy'])
    done,context,audit=(common.read(inputs/n) for n in ('completed.json','context.json','native_audit.json'))
    if not done.get('executionComplete') or not done.get('protectedPostRunVerified') or not all(audit.get(k) for k in
            ('exactFiveBinSumPassed','countErrorFlowPreservedBySummation','writtenCacheReadbackExact')):
        raise ValueError('Coarse input addition/readback audit incomplete.')
    protected=dict(current['protectedFiles'])
    for name in ('cache.root','context.json','native_audit.json'):
        if common.digest(inputs/name)!=done['products'][name]: raise ValueError('Changed coarse input artifact: '+name)
        protected[str(inputs/name)]=done['products'][name]
    axis=audit['coarseAxis']
    if abs(axis['rBinWidth']-cfg['fitRigidityWidth'])>base['numeric']['axisTolerance'] or axis['dedxBinWidth']!=base['numeric']['nativeDedxBinWidth']:
        raise ValueError('Coarse geometry differs from declared fit input.')
    original=Path(context['inputFingerprint']['path'])
    if '2021P24iy_partial_merge' not in original.name or common.fingerprint(original)!=context['inputFingerprint']:
        raise ValueError('Original 2021 P24iy ROOT changed.')
    native=Path(context['nativeCache']['path'])
    if common.digest(native)!=context['nativeCache']['sha256']: raise ValueError('Original native cache changed.')
    protected[str(native)]=context['nativeCache']['sha256']
    window=high.require_exact_window('he4',base['targets']['he4'],current,context)
    protected[str(inputs/'completed.json')]=common.digest(inputs/'completed.json')
    for p in (settings,base_path): protected[str(p)]=common.digest(p)
    references={}
    for response in base['responseModels']:
        rows,signatures=common.tag_index(common.resolve(base['tagStudies'][response]['linear']))
        references[response]=rows;protected.update(signatures)
    effective=dict(base,nativeWidth=cfg['fitRigidityWidth'],aggregateWidth=cfg['fitRigidityWidth'],
        rigidityMin=cfg['rigidityMin'],rigidityMax=cfg['rigidityMax'],pilotRigidity=cfg['pilotRigidity'])
    effective['numeric']=dict(base['numeric'],nativeRigidityBinWidth=cfg['fitRigidityWidth'])
    shapes,gates,diagnostics,windows=[],[],[],[]
    for index,rbin in enumerate(high.selected_bins(effective,args.stage)):
        lo,hi=(rbin-1)*cfg['fitRigidityWidth'],rbin*cfg['fitRigidityWidth'];mid=.5*(lo+hi)
        shape_nodes={};failures={s:[] for s in high.SPECIES};source_rows={};window_references=[]
        for q,(r,weight) in enumerate(provider.gl3_nodes(lo,hi)):
            for s,species in enumerate(high.SPECIES):
                mean,sigma,active,_=high.gate_parameters(current,provider,s,r)
                gates.append([index,q,species,mean,sigma,int(active)])
                sample=dict(d='TAG_D',t='TAG_T',he3='TAG_HE3',he4='TAG_D')[species]
                seed_r=mid if species=='he4' else r
                old=references[args.response].get((sample,int(math.floor(seed_r/cfg['sourceNativeWidth']))+1))
                if old is not None:
                    component=old['components'][species]
                    fail=common.shape_failures(old,component,base['support'])
                    if any(p['atBoundary'] for p in common.named_component_boundaries(old,species,args.response,axis)):
                        fail.append('named_response_math_boundary')
                    mu,width=component['mean'],component['sigma'];failures[species].extend(fail)
                    source_rows[species]=dict(policy='Native conditional reference, initialization for projected target; nuisance varies at GL3 nodes.',sample=sample)
                else:
                    mu,width,_=provider.current_parameters(current,s,seed_r)
                    if species!='he4': failures[species].append('provider_only_nuisance_response_unvalidated')
                    source_rows[species]=dict(policy='Provider initialization only; no validated TAG extrapolation beyond R=2.5.',sample=sample)
                if not common.positive(mu) or not common.positive(width): raise ValueError('Invalid coarse seed, no fallback.')
                shape_nodes[q,species]=(mu,width)
                if species=='he4':
                    window_references.append((mean,sigma))
                    for response in base['responseModels']:
                        ref=references[response].get((sample,int(math.floor(mid/cfg['sourceNativeWidth']))+1))
                        if ref is not None: window_references.append((ref['components'][species]['mean'],ref['components'][species]['sigma']))
                        else: window_references.append(provider.current_parameters(current,s,mid)[:2])
        half,align=(base['allFit'][k] for k in ('fitWindowSigmaHalfWidth','fitWindowAlignmentKeV'))
        fitlo=max(axis['dedxRange'][0],math.floor(min(m-half*s for m,s in window_references)/align)*align)
        fithi=min(axis['dedxRange'][1],math.ceil(max(m+half*s for m,s in window_references)/align)*align)
        for q,(r,weight) in enumerate(provider.gl3_nodes(lo,hi)):
            for species in high.SPECIES:
                mu,width=shape_nodes[q,species]
                mu,width=common.varied_response(dict(mean=mu,sigma=width),args.variation,base,species,'he4')
                shapes.append([index,rbin,lo,hi,q,weight,r,species,mu,width,int(not failures[species]),fitlo,fithi])
        windows.append(dict(index=index,rBin=rbin,lo=lo,hi=hi,fitLo=fitlo,fitHi=fithi))
        diagnostics.append(dict(index=index,rBin=rbin,lo=lo,hi=hi,failures=failures,tagDiagnostics=source_rows))
    output=args.output_dir.resolve();allowed=(common.PROJECT/cfg['outputBase']).resolve()
    if output.exists() or output.parent!=allowed: raise ValueError('Require new coarse-fit output child.')
    output.mkdir();source=output/'source'
    files=[Path(__file__).resolve(),settings,base_path,high.HERE/'FitHighRigidity.C',high.HERE/'run_highrigidity.py',
        high.HERE.parent/'FitLocalNuclearCutoff.C',high.HERE.parent/'run_cutoff.py',high.HERE.parent/'summarize_cutoff.py',
        common.PROJECT/'tools/nuclear_pid_tof/NuclearFitCommon.h',common.PROJECT/'tools/nuclear_pid_tof/NuclearLognormalFitCommon.h',
        common.PROJECT/'tools/nuclear_pid_tof/runtime_selection.py',common.PROJECT/'script/singularity_root_batch.sh']
    for p in files:
        dest=source/p.relative_to(common.PROJECT);dest.parent.mkdir(parents=True,exist_ok=True)
        expected=common.digest(p);shutil.copy2(p,dest)
        if common.digest(dest)!=expected: raise ValueError('Source changed during snapshot.')
        protected[str(p)]=expected;protected[str(dest)]=expected
    for name in ('NuclearFitCommon.h','NuclearLognormalFitCommon.h'):
        dest=source/'tools/nuclear_pid_cutoff'/name;src=common.PROJECT/'tools/nuclear_pid_tof'/name
        shutil.copy2(src,dest);protected[str(dest)]=common.digest(src)
    common.tsv(output/'shapes.tsv','index rBin rLo rHi node weight rValue species mean sigma supported fitLo fitHi',shapes)
    common.tsv(output/'gates.tsv','index node species mean sigma active',gates)
    numeric=dict(effective['numeric'],responseModel=args.response,backgroundModel=args.background,
        sampleHistogram='TAG_D',targetSpecies='he4',mass2Min=window[0],mass2Max=window[1],exactProductionTofWindow=0,
        rigidityMin=cfg['rigidityMin'],rigidityMax=cfg['rigidityMax'],tofThreshold=current['cuts']['highRigidityTof']['thresholdPOverQ'])
    for k,v in base['allFit'].items():
        if k not in ('fitWindowSigmaHalfWidth','fitWindowAlignmentKeV','requireBestSpecies'): numeric[k]=v
    common.tsv(output/'settings.tsv','',sorted(numeric.items()))
    high.save(output/'effective_settings.json',effective);high.save(output/'coarse_settings.json',cfg)
    high.save(output/'fixed_windows.json',windows);high.save(output/'tag_shape_audit.json',diagnostics)
    high.save(output/'current_selection.json',current)
    high.save(output/'context.json',dict(schemaVersion=1,productionEnabled=False,createdUtc=datetime.now(timezone.utc).isoformat(),
        stage='coarse_'+args.stage,targetSpecies='he4',sampleHistogram='TAG_D',response=args.response,background=args.background,
        variation=args.variation,exactProductionTofWindow=False,mass2Window=list(window),mass2Bounds='strict',
        runtimeMainconf=str(mainconf),sourceOriginalInput=context['inputFingerprint'],cache=str(inputs/'cache.root'),
        nativeCache=context['nativeCache'],fitRigidityWidth=cfg['fitRigidityWidth'],sourceNativeRigidityWidth=cfg['sourceNativeWidth'],
        protectedSources=protected,
        responsePolicy='Full 0.05-R projected target shape is free mean/width. Not an intrinsic calibration. Nuisance references vary at GL3 nodes.',
        projectionPolicy='Exact 5-native-R-bin count/error/flow addition; independent new fit, not an average/replacement of native purities.',
        selectionPolicy='Original 4-way TPC winner then user-requested same-D TOF gate; no reassignment; production He4 3-4 remains different.',
        likelihoodPolicy='Coarse TAG_D sample only; original same-data response fits are initialization, not independent priors.',
        fullUncertaintyValidated=False))
    return effective,mainconf,inputs/'cache.root',original,context,protected,output,source/'tools/nuclear_pid_cutoff/highrigidity',diagnostics


def run(args):
    cfg,mainconf,cache,original,context,protected,output,snapshot,diagnostics=prepare(args)
    expression=str(snapshot/'FitHighRigidity.C')+'+('+','.join(json.dumps(str(p)) for p in
        (cache,output,output/'settings.tsv',output/'shapes.tsv',output/'gates.tsv'))+')'
    command=['bash',str(common.PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',expression]
    high.save(output/'command.json',command);started=time.monotonic()
    try:
        with (output/'fit.log').open('w') as log:
            result=subprocess.run(command,cwd=common.PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=cfg['executionTimeoutSeconds'])
        if result.returncode or 'NUCLEAR_HIGH_R_COMPLETE' not in (output/'fit.log').read_text():
            raise RuntimeError('Coarse ROOT fit failed; retained fit.log.')
        rows=high.summarize_run(cfg,output,'he4',diagnostics)
        for p,expected in protected.items():
            if common.digest(p)!=expected: raise ValueError('Protected coarse source/input changed: '+p)
        if common.fingerprint(original)!=context['inputFingerprint']: raise ValueError('Original ROOT fingerprint changed.')
        high.save(output/'completed.json',dict(executionComplete=True,productionEnabled=False,fullUncertaintyValidated=False,
            protectedSourcesVerified=True,originalInputPreserved=True,runtimeSeconds=time.monotonic()-started,
            nativeRows=len(rows),supportedRows=sum(r['quality'] for r in rows),targetSpecies='he4',
            stage='coarse_'+args.stage,response=args.response,background=args.background,variation=args.variation,
            sourceNativeRigidityWidth=.01,fitRigidityWidth=.05,
            products={str(p.relative_to(output)):common.digest(p) for p in output.rglob('*') if p.is_file() and p.suffix not in ('.so','.d','.pcm')}))
        print('NUCLEAR_COARSE_HIGH_R_COMPLETE',output,'supported',sum(r['quality'] for r in rows),'/',len(rows),flush=True)
    except Exception as error:
        high.save(output/'failed.json',dict(executionComplete=False,error=str(error)));raise


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--stage',required=True,choices=('pilot','native'))
    parser.add_argument('--response',required=True,choices=('gaussian','lognormal'))
    parser.add_argument('--background',required=True,choices=('linear','exponential'))
    parser.add_argument('--variation',required=True)
    run(parser.parse_args())
