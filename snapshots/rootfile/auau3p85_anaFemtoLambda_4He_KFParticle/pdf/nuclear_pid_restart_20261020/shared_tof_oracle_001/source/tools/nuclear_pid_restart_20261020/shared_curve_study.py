#!/usr/bin/env python3
"""Fresh, frozen common-polynomial nuclear separation; never edits production."""
import argparse
import csv
import math
from pathlib import Path
import re
import shutil
import subprocess
import yaml
from prepare_tag_study import PROJECT,SPECIES,require,owned,sha,read,save,tsv,quote,gl3
from polynomial_calibration import polynomial,COEFFICIENT_HEADER

def coefficients(folder):
    meta=read(folder/'calibration.json')
    require(meta['productionEnabled'] is False,'Candidate calibration only')
    require(all(meta['species'][s]['curveCandidateValid'] for s in SPECIES),'All four common curves must be valid candidates')
    with (folder/'coefficients.tsv').open() as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        require(reader.fieldnames==COEFFICIENT_HEADER,'Unexpected coefficient schema')
        rows={}
        for q in reader:
            name=q.pop('species');require(name in SPECIES and name not in rows,'Duplicate/unknown species')
            rows[name]={k:float(v) for k,v in q.items()}
    require(set(rows)==set(SPECIES),'Missing common nuclear curve')
    return rows

def evaluate(table,r):
    require(table['rMin']-1.e-10<=r<=table['rMax']+1.e-10,'No unsupported coefficient extrapolation')
    x=math.log10(r)
    mu=polynomial([table['mean%d'%i] for i in range(4)],x)
    sigma=polynomial([table['sigma%d'%i] for i in range(3)],x)
    ratio=math.exp(polynomial([table['logRatio%d'%i] for i in range(3)],x))
    require(all(math.isfinite(v) and v>0 for v in (mu,sigma,ratio)) and .25<=ratio<=4.,'Invalid candidate response')
    return mu,sigma,ratio

def select_roi(shapes,target,cfg,variant,axis):
    center={s:tuple(sum(w*v[k] for w,v in shapes[s]) for k in range(3)) for s in SPECIES}
    mu,sigma,_=center[target];half=variant['targetFitHalfWidthSigma']
    lo,hi=mu-half*sigma,mu+half*sigma
    active={target}|set(cfg['requiredNeighbors'][target])
    for _ in range(6):
        old=set(active)
        for s in active:
            m,e,_=center[s];pad=half if s==target else cfg['neighborHalfWidthSigma']
            lo=min(lo,m-pad*e);hi=max(hi,m+pad*e)
        active.update(s for s in SPECIES if lo<=center[s][0]<=hi)
        if active==old: break
    step=axis['dedxBinWidth'];amin,amax=axis['dedxRange']
    clipped=lo<amin or hi>amax
    lo=max(amin,math.floor(lo/step)*step);hi=min(amax,math.ceil(hi/step)*step)
    require(lo<hi,'Empty local fit ROI')
    return lo,hi,active,center,clipped

def prepare(settings,coefficient_dir,output,variant_name,smoke,diagnostic_mode=None):
    cfg=yaml.safe_load(settings.read_text())
    require(cfg['schemaVersion']==1 and not output.exists() and output.parent==owned(cfg['outputBase'],False),'Fresh direct study output required')
    require(cfg['tofWindows']['d']==cfg['tofWindows']['he4'],'d and He4 must use identical TOF window')
    require(cfg['samplesAboveSwitch']==dict(d='TAG_D',t='TAG_T',he3='TAG_HE3',he4='TAG_D'),'Species TOF mapping mismatch')
    require(variant_name in cfg['variants'],'Unknown explicit model variant')
    variant=cfg['variants'][variant_name];tables=coefficients(coefficient_dir)
    source_rows,original_nodes={},{}
    if diagnostic_mode:
        source_dir=Path(read(coefficient_dir.parent/'completed.json')['sourceFitPath']).parent
        source_rows={q['id']:q for q in read(source_dir/'fits.json')['results']}
        with (source_dir/'nodes.tsv').open() as stream:
            for q in csv.DictReader(stream,delimiter='\t'):
                original_nodes.setdefault((q['id'],q['species']),[]).append(tuple(float(q[k]) for k in ('weight','mean','sigma')))
    inputs=owned(cfg['inputDirectory']);done,context,audit=(read(inputs/n) for n in ('completed.json','context.json','native_audit.json'))
    require(done['executionComplete'] and done['protectedPostRunVerified'] and audit['writtenCacheReadbackExact'],'Unverified native inputs')
    require('2021P24iy_partial_merge' in Path(context['inputFingerprint']['path']).name,'2021 P24iy input required')
    protected={str(inputs/n):sha(inputs/n) for n in ('cache.root','completed.json','context.json','native_audit.json')}
    for n in ('cache.root','context.json','native_audit.json'):
        require(protected[str(inputs/n)]==done['products'][n],'Input hash differs')
    for s,p in (('d','D'),('t','T'),('he3','HE3')):
        require(cfg['tofWindows'][s]==[context['settings']['tag'+p+'Mass2Min'],context['settings']['tag'+p+'Mass2Max']],'TOF windows differ from saved sample')
    axis=audit['nativeAxis'];bins=cfg['smokeNativeRBins'] if smoke else sorted(set(range(cfg['nativeRBinFirst'],cfg['nativeRBinLast']+1,cfg['nativeRBinStride']))|set(cfg['extraNativeRBins']))
    records,slices,nodes,skipped=[],[],[],[]
    for rb in bins:
        require(type(rb) is int and 1<=rb<=axis['rBins'],'Invalid R bin')
        lo=axis['rRange'][0]+(rb-1)*axis['rBinWidth'];hi=lo+axis['rBinWidth'];r=(lo+hi)/2.
        unsupported=[s for s in SPECIES if lo<tables[s]['rMin']-1.e-10 or hi>tables[s]['rMax']+1.e-10]
        if unsupported:
            skipped.extend(dict(id=s+'_b'+str(rb).zfill(4),target=s,r=r,rLo=lo,rHi=hi,reason='outside_common_TOF_calibration_support',unsupportedSpecies=unsupported,purity=None,lowPurityInferred=False) for s in SPECIES)
            continue
        shape={s:[(w,evaluate(tables[s],q)) for q,w in gl3(lo,hi)] for s in SPECIES}
        for target in SPECIES:
            identity=target+'_b'+str(rb).zfill(4)
            mode='ALL' if hi<=cfg['tofSwitchR']+1.e-10 else cfg['samplesAboveSwitch'][target]
            require(mode=='ALL' or lo>=cfg['tofSwitchR']-1.e-10,'R bin straddles transition')
            current_shape=dict(shape)
            if diagnostic_mode:
                source=source_rows[identity];require(source['fitCompleted'],'Incomplete oracle source')
                mode=cfg['samplesAboveSwitch'][target];fitlo,fithi=source['fitLo'],source['fitHi'];active={target};clipped=False
                if diagnostic_mode=='tof_oracle':
                    c=next(q for q in source['components'] if q['species']==target)
                    shift=c['mean']-c['initialMean'];scale=c['sigma']/c['initialSigma']
                    current_shape[target]=[(w,(m+shift,e*scale,source['responseWidthRatio'])) for w,m,e in original_nodes[(identity,target)]]
                centers={s:tuple(sum(w*v[k] for w,v in current_shape[s]) for k in range(3)) for s in SPECIES}
            else:
                fitlo,fithi,active,centers,clipped=select_roi(current_shape,target,cfg,variant,axis)
            mu,sigma,_=centers[target];n=cfg['gateNSigma']
            slices.append((identity,mode,target,rb,fitlo,fithi,mu-n*sigma,mu+n*sigma))
            for s in SPECIES:
                nodes.extend((identity,s,w,m,e,q,int(s in active)) for w,(m,e,q) in current_shape[s])
            records.append(dict(id=identity,target=target,mode=mode,rBin=rb,r=r,rLo=lo,rHi=hi,fitLo=fitlo,fitHi=fithi,mean=mu,sigma=sigma,activeSpecies=sorted(active),fitAxisClipped=clipped,withinTOFCalibrationSupport=True))
    require(slices,'No common-supported R slices')
    output.mkdir(parents=True)
    base=yaml.safe_load(owned(cfg['baseNumericSettings']).read_text())
    numeric=dict(base['numeric']);numeric.update(cfg['numeric']);numeric.update(backgroundModel=variant['backgroundModel'],backgroundSlope=0.,freeBackgroundSlope=0,gatePolicy=cfg['gatePolicy'],gateNSigma=cfg['gateNSigma'],tofTransitionR=cfg['tofSwitchR'])
    numeric['tofTransferDiagnostic']=int(bool(diagnostic_mode))
    numeric.update({'hist.'+k:k for k in ('ALL','TAG_D','TAG_T','TAG_HE3')})
    numeric.update(expectedRBins=axis['rBins'],expectedRMin=axis['rRange'][0],expectedRMax=axis['rRange'][1],expectedDedxBins=axis['dedxBins'],expectedDedxMin=axis['dedxRange'][0],expectedDedxMax=axis['dedxRange'][1])
    tsv(output/'settings.tsv',('key','value'),numeric.items())
    tsv(output/'slices.tsv',('id','mode','target','rBin','fitLo','fitHi','gateLo','gateHi'),slices)
    tsv(output/'nodes.tsv',('id','species','weight','mean','sigma','widthRatio','active'),nodes)
    sources={settings,Path(__file__).resolve(),owned(cfg['baseNumericSettings']),owned(cfg['fitSource']),owned(cfg['mainconf']),owned(cfg['wrapper'])}
    sources.update(owned(p) for p in cfg['extraSources'])
    sources.update(p.resolve() for p in Path(__file__).resolve().parent.glob('*.py'))
    for name in ('DrawPidRecalibration.C','DrawSharedPidFits.C','RunRestartMacro.C','TestAllFixedResponseFit.C'):
        p=Path(__file__).resolve().with_name(name)
        if p.exists(): sources.add(p)
    main=yaml.safe_load(owned(cfg['mainconf']).read_text());sources.update(owned('config/'+ref) for ref in main.values())
    nuclear=yaml.safe_load(owned('config/'+main['nuclearid']).read_text());require(nuclear['dedxCalibration']['enabled'] is False,'Production must remain original Git PID')
    for path in (coefficient_dir/'coefficients.tsv',coefficient_dir/'calibration.json',coefficient_dir/'coefficients.yaml'):
        protected[str(path)]=sha(path)
    if diagnostic_mode:
        for path in (source_dir/'fits.json',source_dir/'nodes.tsv'):protected[str(path)]=sha(path)
    for original in sources:
        snapshot=output/'source'/original.relative_to(PROJECT);snapshot.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(str(original),str(snapshot));protected[str(snapshot)]=sha(snapshot)
    for path in output.glob('*.tsv'):protected[str(path)]=sha(path)
    save(output/'plan.json',dict(schemaVersion=1,config=cfg,variant=variant_name,smoke=smoke,diagnosticMode=diagnostic_mode,coefficientDirectory=str(coefficient_dir),slices=records,skippedSlices=skipped,totalFits=len(records),inputs=protected,cache=str(inputs/'cache.root'),inputFingerprint=context['inputFingerprint'],productionEnabled=False,withinTOFCalibrationSupportOnly=True,limitations=['Centrality-integrated, charge-inclusive track QA, not Lambda-pair-weighted purity.','ALL and TAG overlap; no independent likelihood addition.','Uniform GL3 within each native R bin.','Fixed-curve conditional statistical errors; model differences must be assessed separately.']))
    print('SHARED_PID_PREPARED fits=%d unsupported=%d %s'%(len(records),len(skipped),output),flush=True)

def unchanged(plan):
    for name,digest in plan['inputs'].items():require(sha(Path(name))==digest,'Protected input changed: '+name)

def execute(output):
    plan=read(output/'plan.json');unchanged(plan);cfg=plan['config']
    script=output/'source'/cfg['fitSource'];args=[plan['cache'],output/'settings.tsv',output/'slices.tsv',output/'nodes.tsv',output]
    invocation=str(script)+'+('+','.join(quote(p) for p in args)+')'
    command=['bash',str(owned(cfg['wrapper'])),str(owned(cfg['mainconf'])),'root4star','-l','-b','-q',invocation]
    log=output/'fit.log';require(not log.exists(),'Preserve previous attempt; choose fresh output')
    save(output/'fit.command.json',dict(argv=command,cwd=str(PROJECT)))
    with log.open('x') as stream:
        subprocess.run(command,cwd=PROJECT,stdout=stream,stderr=subprocess.STDOUT,timeout=cfg['executionTimeoutSeconds'],check=True)
    require('ALL_RESTART_FITS_COMPLETE slices='+str(plan['totalFits']) in log.read_text(),'Missing completed fit/readback marker')
    result=read(output/'fits.json')
    require(result['status']=='completed' and result['rawReadback'] and result['predictionReadback'] and result['mappingReadback'] and result['rawCountsModified'] is False and result['totalSlices']==plan['totalFits'],'Incomplete or modified outputs')
    for q,want in zip(result['results'],plan['slices']):require(all(q[k]==want[k] for k in ('id','target','mode','rBin')),'Changed requested R slice')
    unchanged(plan);save(output/'completed.json',dict(status='completed',productionEnabled=False,inputsUnchanged=True,totalFits=result['totalSlices'],fitFailures=result['fitFailures'],qualityAcceptedSlices=result['qualityAcceptedSlices'],outputs={n:sha(output/n) for n in ('fits.json','native_spectra.root','native_effective_bins.tsv')}))
    print('SHARED_PID_FITS_COMPLETE accepted=%d/%d %s'%(result['qualityAcceptedSlices'],result['totalSlices'],output),flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('prepare','fit'));p.add_argument('--settings');p.add_argument('--coefficient-dir');p.add_argument('--variant',choices=('linear','constant','wider'));p.add_argument('--smoke',action='store_true');p.add_argument('--diagnostic-mode',choices=('tof_fixed','tof_oracle'));p.add_argument('--output-dir',required=True)
    a=p.parse_args();out=owned(a.output_dir,a.action!='prepare')
    if a.action=='prepare':
        require(a.settings and a.coefficient_dir and a.variant,'Explicit settings/coefficient/variant required');prepare(owned(a.settings),owned(a.coefficient_dir),out,a.variant,a.smoke,a.diagnostic_mode)
    else:
        require(not (a.settings or a.coefficient_dir or a.variant or a.smoke or a.diagnostic_mode),'Execution uses frozen plan');execute(out)
if __name__=='__main__': main()
