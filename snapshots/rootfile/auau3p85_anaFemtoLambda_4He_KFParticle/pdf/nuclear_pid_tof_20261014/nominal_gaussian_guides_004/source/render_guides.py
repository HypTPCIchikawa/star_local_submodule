#!/usr/bin/env python3
"""Display-only native TAG Gaussian points and explicitly unsupported curve gaps."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import yaml

from nuclear_curves import evaluate_curve

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]
SPECIES=('d','t','he3','he4')
PRIMARY={'d':'TAG_D','t':'TAG_T','he3':'TAG_HE3','he4':'TAG_D'}
CALIBRATION_HELPER=PROJECT/'tools/he4_purity/export_runtime_calibration.py'
spec=importlib.util.spec_from_file_location('frozen_legacy_reference',CALIBRATION_HELPER)
reference=importlib.util.module_from_spec(spec);spec.loader.exec_module(reference)

def digest(p):
    h=hashlib.sha256()
    with p.open('rb') as f:
        for b in iter(lambda:f.read(1048576),b''):h.update(b)
    return h.hexdigest()

def fingerprint(p):
    s=p.stat();return dict(path=str(p),size=s.st_size,mtime_ns=s.st_mtime_ns,device=s.st_dev,inode=s.st_ino)

def save(p,x):p.write_text(json.dumps(x,indent=2,allow_nan=False)+'\n')

def table(p,header,rows):
    with p.open('w',newline='') as f:
        w=csv.writer(f,delimiter='\t',lineterminator='\n');w.writerow(header);w.writerows(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    for k in ('input-study','tag-study','alternative-study','curve-study','mainconf','working-calibration','settings','output-dir'):
        p.add_argument('--'+k,type=Path,required=True)
    a=p.parse_args();cfg_path=a.settings.resolve(strict=True);cfg=yaml.safe_load(cfg_path.read_text())
    if set(cfg)!={'schemaVersion','executionTimeoutSeconds','display'} or cfg['schemaVersion']!=1:raise ValueError('Unknown overlay settings schema.')
    d=cfg['display']
    if not all(isinstance(v,(int,float)) and not isinstance(v,bool) and math.isfinite(v) for v in d.values()):raise ValueError('Display-only settings must be finite numeric.')
    if (d['schemaVersion'],d['expectedRBins'],d['expectedDedxBins'],d['dedxDisplayRebin'],d['guideNSigma'])!=(1,600,1200,5,2):raise ValueError('Fixed native cache and +/-2 projected width schema required.')
    if not (0<d['curveDisplayStep']<=.01 and 0<d['fullRMin']<d['fullRMax']<=6 and 0<d['lowRMin']<d['lowRMax']<=6):raise ValueError('Invalid display ranges/sampling.')
    if not isinstance(cfg['executionTimeoutSeconds'],int) or cfg['executionTimeoutSeconds']<=0:raise ValueError('Explicit positive timeout required.')
    mainconf=a.mainconf.resolve(strict=True);runtime=yaml.safe_load(mainconf.read_text());info_path=(PROJECT/'config'/runtime['analysis']).resolve(strict=True)
    info=yaml.safe_load(info_path.read_text());tag=info['starTag']
    if (tag['collision'],tag['year'],tag['productionTag'],tag['libraryTag'])!=('auau3p85','year2021','P24iy','SL24y'):raise ValueError('Explicit 3p85 2021P24iy ROOT5 runtime required.')
    inputs,tags,alternative,curve=(getattr(a,k).resolve(strict=True) for k in ('input_study','tag_study','alternative_study','curve_study'))
    calibration_path=a.working_calibration.resolve(strict=True);calibration=yaml.safe_load(calibration_path.read_text())
    own=[Path(__file__).resolve(),HERE/'PlotNuclearGuides.C',cfg_path,HERE/'nuclear_curves.py',HERE/'runtime_selection.py',CALIBRATION_HELPER,calibration_path,mainconf,info_path,PROJECT/'script/singularity_root_batch.sh',PROJECT/'StMaker/common/NuclearIdDeDxVsMom.h']
    protected={str(x):digest(x) for x in own};initial=None
    def add(x,h):
        if str(x) in protected and protected[str(x)]!=h:raise ValueError('Conflicting protected hash: '+str(x))
        protected[str(x)]=h
    for folder in (inputs,tags,alternative):
        cp,mp=folder/'context.json',folder/'completed.json';context=json.loads(cp.read_text());complete=json.loads(mp.read_text())
        if complete.get('executionComplete') is not True or complete.get('productionEnabled') is not False or context.get('productionEnabled') is not False:raise ValueError('Only authoritative complete unadopted studies supported.')
        if Path(context['runtimeMainconf']).resolve(strict=True)!=mainconf:raise ValueError('Study runtime mismatch.')
        fp=context['inputFingerprint']
        if initial is not None and initial!=fp:raise ValueError('Mixed original inputs forbidden.')
        initial=fp
        for field in ('protectedPreRun','protected','snapshots'):
            for encoded,h in context.get(field,{}).items():add(Path(encoded),h)
        for rel,h in complete['products'].items():
            target=(folder/rel).resolve(strict=True)
            if folder not in target.parents:raise ValueError('Study product escapes folder.')
            add(target,h)
        for rel,h in complete.get('sourceSha256',{}).items():add(folder/'source'/rel,h)
        add(cp,digest(cp));add(mp,digest(mp))
    curve_context=json.loads((curve/'context.json').read_text());curve_complete=json.loads((curve/'completed.json').read_text())
    if curve_complete.get('executionComplete') is not True or curve_complete.get('productionEnabled') is not False:raise ValueError('Curve study incomplete/adopted.')
    for rel,h in curve_complete['products'].items():add(curve/rel,h)
    for rel,h in curve_context['snapshots'].items():add(curve/'source'/rel,h)
    for encoded,h in curve_context['tagSources'].items():add(Path(encoded),h)
    add(Path(curve_context['settingsPath']),curve_context['settingsSha256'])
    for name in ('context.json','completed.json','support_audit.json'):add(curve/name,digest(curve/name))
    def preserve():
        if fingerprint(Path(initial['path']))!=initial:raise ValueError('Original ROOT fingerprint changed.')
        for encoded,h in protected.items():
            if digest(Path(encoded))!=h:raise ValueError('Protected source/input/product changed: '+encoded)
    preserve()
    primary=json.loads((tags/'tag_fits.json').read_text());alt=json.loads((alternative/'tag_fits.json').read_text())
    for x in (primary,alt):
        if x.get('schemaVersion')!=1 or x.get('fitStage')!='independent_tag_shapes' or x.get('responseModel','gaussian')!='gaussian' or x.get('productionEnabled') is not False:raise ValueError('This guide utility is explicit Gaussian only; log-normal guide semantics require a separate study.')
    key=lambda x:(x['sample'],round(x['lo'],10),round(x['hi'],10))
    keyed={key(x):x for x in alt['slices']}
    if len(keyed)!=len(alt['slices']) or set(keyed)!=set(map(key,primary['slices'])):raise ValueError('Alternative backgrounds must share exact native slices.')
    support=json.loads((curve/'support_audit.json').read_text())
    support_lookup={(x['species'],round(x['lo'],10),round(x['hi'],10)):x for x in support}
    if len(support_lookup)!=len(support):raise ValueError('Duplicate support intervals.')
    points=[];omitted=[]
    for s in primary['slices']:
        other=keyed[key(s)]
        if (s['fitLo'],s['fitHi'],s['firstNativeRBin'],s['lastNativeRBin'])!=(other['fitLo'],other['fitHi'],other['firstNativeRBin'],other['lastNativeRBin']):raise ValueError('Different actual input/fit window.')
        if s['firstNativeRBin']!=s['lastNativeRBin'] or abs(s['hi']-s['lo']-.01)>d['axisRelativeTolerance']:raise ValueError('One native R bin per point required.')
        for name in SPECIES:
            if s['sample']!=PRIMARY[name]:continue
            c=s.get('components',{}).get(name);record=support_lookup[(name,round(s['lo'],10),round(s['hi'],10))]
            if not c or not all(math.isfinite(c[k]) for k in ('mean','sigma')) or c['sigma']<=0:
                omitted.append(dict(species=name,index=s['index'],reason='No finite positive saved width; cannot approximate guide.'));continue
            if not c['freeMean'] or not c['freeSigma']:raise ValueError('Primary component is not a free center/width.')
            accepted=bool(record['acceptedNumericalShape'])
            points.append([name,s['sample'],s['index'],s['lo'],s['hi'],s['r'],c['mean'],c['sigma'],int(accepted),';'.join(record['failures']) or '-',int(c['zeroYield']),int(s['atBoundary']),int(s['numericalConverged'])])
    candidate=json.loads((curve/'pchip_logR.json').read_text())
    if candidate.get('productionEnabled') is not False or candidate.get('identityValidated') is not False:raise ValueError('Unadopted reference candidate required.')
    # He4 has insufficient contiguous accepted points: never fabricate its curve.
    if candidate['species']['he4']['segments']:raise ValueError('This nominal-He4-no-supported-curve display contract must be revised for a future candidate.')
    curves=[];counts={}
    lo=min(d['fullRMin'],d['lowRMin']);hi=max(d['fullRMax'],d['lowRMax']);step=d['curveDisplayStep']
    grids=[lo+j*step for j in range(int(math.ceil((hi-lo)/step))+1)];grids[-1]=hi
    for index,name in enumerate(SPECIES):
        # Keep the real provider discontinuities disconnected, never replace them.
        for family in ('legacy','working'):
            group=0;previous=None;buffer=[]
            def flush():
                if len(buffer)>=2:curves.extend(buffer)
                buffer.clear()
            for r in grids:
                if family=='legacy':mean,sigma=reference.legacy_parameters(index,r);state='legacy'
                else:
                    mean,sigma,state=reference.calibrated_parameters(calibration,index,r)
                    state='legacy' if state.startswith('legacy_') else 'working'
                invalid=not math.isfinite(mean) or not math.isfinite(sigma) or mean<=0 or sigma<=0
                discontinuity=previous is not None and (state!=previous[0] or abs(mean-previous[1])>d['workingDiscontinuityThresholdKeV'])
                if invalid or discontinuity:flush();group+=1
                if invalid:previous=None;continue
                buffer.append([name,family,group,r,mean,sigma]);previous=(state,mean)
            flush()
        total=0
        for group,segment in enumerate(candidate['species'][name]['segments']):
            for interval_index,(rlo,rhi) in enumerate(segment['supportedNativeIntervals']):
                n=max(1,int(math.ceil((rhi-rlo)/step)))
                for j in range(n+1):
                    r=rlo+(rhi-rlo)*j/n;value=evaluate_curve(candidate,name,r)
                    if not value['supported']:raise ValueError('Candidate evaluated outside explicit supported native interval.')
                    curves.append([name,'candidate',str(group)+'_'+str(interval_index),r,value['mean'],value['sigma']]);total+=1
        counts[name]=dict(totalFreePoints=sum(x[0]==name for x in points),acceptedNumericalShapePoints=sum(x[0]==name and x[8] for x in points),candidateDisplayPoints=total,candidateSegments=len(candidate['species'][name]['segments']),identityValidated=False)
    output=a.output_dir.resolve();allowed=PROJECT/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014'
    if output.exists() or allowed not in output.parents or any(folder in output.parents for folder in (inputs,tags,alternative,curve)):raise ValueError('New sibling output required; never overwrite source studies.')
    output.mkdir(parents=True);snapshot=output/'source';snapshot.mkdir()
    for source in own:shutil.copy2(source,snapshot/source.name)
    table(output/'points.tsv',['species','sample','index','lo','hi','r','mean','sigma','acceptedNumericalShape','failureFlags','zeroYield','atBoundary','numericalConverged'],points)
    table(output/'curves.tsv',['species','family','group','r','mean','sigma'],curves)
    (output/'settings.tsv').write_text(''.join('%s\t%.17g\n'%(k,v) for k,v in d.items()))
    save(output/'guide_audit.json',dict(noFit=True,productionEnabled=False,identityValidated=False,species=counts,omitted=omitted,failedPointPolicy='Red unlinked markers with +/-2 projected width bars; never interpolated.',DDisplayPolicy='TAG_D reused for d and He4-compatible pages: same events, not independent/pure He4.',candidatePolicy='Only explicit supportedNativeIntervals, R-unfolded seed widths, NOT validated transfer/intrinsic calibration.',workingPolicy='Exact current knot-provider reference, jumps disconnected; no production changes.',responseModel='gaussian'))
    save(output/'context.json',dict(schemaVersion=1,noFit=True,productionEnabled=False,inputFingerprint=initial,runtimeMainconf=str(mainconf),sourceStudies=[str(x) for x in (inputs,tags,alternative,curve)],protectedPreRun=protected,display=d))
    call=str(snapshot/'PlotNuclearGuides.C')+'++('+','.join(json.dumps(str(x)) for x in (inputs/'cache.root',output/'points.tsv',output/'curves.tsv',output/'settings.tsv',output))+')'
    command=['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',call];save(output/'command.json',command)
    with (output/'root.log').open('w') as log:result=subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=cfg['executionTimeoutSeconds'])
    if result.returncode or 'NUCLEAR_GUIDES_COMPLETE' not in (output/'root.log').read_text():raise RuntimeError('Display failed; log retained.')
    preserve()
    save(output/'completed.json',dict(executionComplete=True,noFit=True,productionEnabled=False,inputFingerprintPreserved=True,protectedPostRunVerified=True,products={x.name:digest(x) for x in output.iterdir() if x.is_file()},sourceSha256={x.name:digest(x) for x in snapshot.iterdir() if x.is_file()}))
    print('NUCLEAR_GUIDES_RENDER_COMPLETE',output,flush=True)

if __name__=='__main__':main()
