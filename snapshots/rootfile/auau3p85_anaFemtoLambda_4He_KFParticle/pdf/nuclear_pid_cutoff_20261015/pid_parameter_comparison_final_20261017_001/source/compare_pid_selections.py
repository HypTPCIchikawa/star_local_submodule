#!/usr/bin/env python3
"""Compare PID response functions and cut windows; never recalibrate or fit."""
import argparse
import datetime
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess
import yaml

HERE=Path(__file__).resolve().parent
PROJECT=HERE.parents[1]
SPECIES=('d','t','he3','he4')

def sha(path):
    value=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):value.update(block)
    return value.hexdigest()

def read(path):return json.loads(path.read_text())
def resolve(path):
    value=Path(path)
    return (value if value.is_absolute() else PROJECT/value).resolve(strict=True)

def write(path,value):
    with path.open('x') as stream:stream.write(json.dumps(value,indent=2,ensure_ascii=False,allow_nan=False)+'\n')

def fingerprint(path):
    v=path.stat()
    return dict(path=str(path),size=v.st_size,mtime_ns=v.st_mtime_ns,device=v.st_dev,inode=v.st_ino)

def tsv(path,keys,rows):
    with path.open('x') as stream:
        stream.write('\t'.join(keys)+'\n')
        for row in rows:
            stream.write('\t'.join('nan' if v is None else format(v,'.17g') if isinstance(v,float) else str(v) for v in row)+'\n')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    args=parser.parse_args();settings=args.settings.resolve(strict=True);cfg=yaml.safe_load(settings.read_text())
    if cfg['schemaVersion']!=1 or not 0<cfg['rigidityMin']<cfg['rigidityMax']:raise ValueError('Invalid geometry.')
    output=args.output_dir.resolve()
    allowed=PROJECT/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015'
    if output.exists() or output.parent!=allowed:raise ValueError('Choose a new direct child of the PID study.')
    mainconf,inputs,original=(resolve(cfg[k]) for k in ('runtimeMainconf','inputStudy','originalRoot'))
    native,done=(read(inputs/n) for n in ('context.json','completed.json'))
    if fingerprint(original)!=native['inputFingerprint']:raise ValueError('Original ROOT fingerprint changed.')
    for name in ('context.json','native_audit.json','cache.root'):
        if sha(inputs/name)!=done['products'][name]:raise ValueError('Changed immutable native input '+name)
    main=yaml.safe_load(mainconf.read_text());info_path=resolve(Path('config')/main['analysis'])
    tags=yaml.safe_load(info_path.read_text())['starTag']
    if tuple(tags[k] for k in ('collision','year','productionTag','libraryTag'))!=('auau3p85','year2021','P24iy','SL24y'):
        raise ValueError('Use 2021 P24iy SL24y only.')
    provider_path=PROJECT/'tools/nuclear_pid_tof/runtime_selection.py'
    spec=importlib.util.spec_from_file_location('pid_compare_runtime',str(provider_path));provider=importlib.util.module_from_spec(spec);spec.loader.exec_module(provider)
    header=resolve(cfg['legacyHeader']);tables=provider.parse_legacy_tables(header)
    profiles={};files=[settings,mainconf,info_path,provider_path,header,Path(__file__).resolve(),HERE/'ComparePIDSelections.C',PROJECT/'script/singularity_root_batch.sh']
    for label,key in (('oldLegacy','legacyNuclearConfig'),('oldCanonical','oldCanonicalNuclearConfig'),('oldTofExtended','oldTofExtendedNuclearConfig'),('new','newNuclearConfig')):
        path=resolve(cfg[key]);cuts=yaml.safe_load(path.read_text());files.append(path);calibration=None
        if cuts.get('dedxCalibration',{}).get('enabled'):
            p=resolve(Path('config')/cuts['dedxCalibration']['calibrationFile']);calibration=yaml.safe_load(p.read_text());files.append(p)
        profiles[label]=dict(tables=tables,calibration=calibration,cuts=cuts)
    if any(p['cuts']['maxNSigmaNuclear']!=cfg['gateNSigma'] for p in profiles.values()):raise ValueError('nSigma windows differ.')
    keys=['species','R']
    for name in ('oldCanonical','oldLegacy','new'):keys.extend([name+'Mean',name+'Sigma',name+'Valid'])
    def values(species,r):
        row=[SPECIES[species],r]
        for label in ('oldCanonical','oldLegacy','new'):
            mean,sigma,_=provider.current_parameters(profiles[label],species,r)
            valid=all(math.isfinite(v) and v>0 for v in (mean,sigma))
            row.extend([mean if valid else None,sigma if valid else None,int(valid)])
        return row
    steps=int(round((cfg['rigidityMax']-cfg['rigidityMin'])/cfg['rigidityStep']))
    rows=[values(s,cfg['rigidityMin']+i*cfg['rigidityStep']) for s in range(4) for i in range(steps+1)]
    reference=[values(s,r) for s in range(4) for r in cfg['referenceRigidity']]
    old=profiles['oldLegacy']['cuts'];final=profiles['new']['cuts'];ext=profiles['oldTofExtended']['cuts'];windows=[]
    for s,(species,key) in enumerate(zip(SPECIES,('Deuteron','Triton','He3','He4'))):
        if s<3:
            mean,sigma=old['nuclearTofMean'+key],old['nuclearTofSigma'+key]
            lo,hi=mean-old['m2SigmaCut']*sigma,mean+old['m2SigmaCut']*sigma
        else:mean,sigma=None,None;lo,hi=ext['he4HighRigidityTof']['mass2Min'],ext['he4HighRigidityTof']['mass2Max']
        current=final['highRigidityTof']['species'][species];suffix=dict(d='D',t='T',he3='HE3',he4='D')[species]
        studied=[native['settings']['tag'+suffix+k] for k in ('Mass2Min','Mass2Max')]
        matches=all(abs(a-b)<1.e-12 for a,b in zip(studied,(current['mass2Min'],current['mass2Max'])))
        windows.append([species,mean,sigma,old['m2SigmaCut'] if s<3 else None,lo,hi,current['mass2Min'],current['mass2Max'],studied[0],studied[1],int(matches)])
    output.mkdir(parents=True);originals={str(p):sha(p) for p in set(files)};source=output/'source';source.mkdir()
    for path in set(files):
        rel=path.relative_to(PROJECT) if PROJECT in path.parents else Path('settings')/path.name
        destination=source/rel;destination.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,destination)
        if sha(destination)!=originals[str(path)]:raise ValueError('Source changed during snapshot.')
    for name in ('context.json','completed.json','native_audit.json'):shutil.copy2(inputs/name,source/('native_'+name))
    tsv(output/'response_curves.tsv',keys,rows);tsv(output/'reference_parameters.tsv',keys,reference)
    tsv(output/'tof_windows.tsv',['species','oldReferenceMean','oldReferenceSigma','oldReferenceNSigma','oldTofExtendedMin','oldTofExtendedMax','currentConfiguredMin','currentConfiguredMax','studiedMin','studiedMax','studyMatchesConfigured'],windows)
    numeric={k:cfg[k] for k in ('rigidityMin','rigidityMax','lowRigidityMax','gateNSigma','tofThreshold','canvasWidth','canvasHeight','dedxDisplayMax','m2DisplayMin','m2DisplayMax','m2DetailMin','m2DetailMax')}
    tsv(output/'plot_settings.tsv',['key','value'],sorted(numeric.items()))
    write(output/'context.json',dict(schemaVersion=1,createdAt=datetime.datetime.now(datetime.timezone.utc).isoformat(),originalFingerprint=native['inputFingerprint'],wholeOriginalRootHashCalculated=False,originalM2Histogram=cfg['originalM2Histogram'],originalOpenedForOneKeyOnly=True,productionEnabled=False,calibrationChanged=False,tofMeanSigmaRefitted=False,oldCanonicalAndNewResponseIdentical=profiles['oldCanonical']['calibration']==profiles['new']['calibration'],responseMeaning='Configured response curves only, not new fitted response or uncertainty.',oldDedicatedProfile='d/t/He3 highpurity entries previously use legacy response, including competing nuclei.',tofWindowsMeaning='Configured gates; He4 studied D window is distinct from production [3,4].',limitations=['No dEdx-m2 joint; M2 map alone cannot estimate final isotope purity.','Inclusive track QA, not pair-conditioned purity or isotope truth.','2021 P24iy partial merge is not full 23B coverage.'],protectedSources=originals))
    call=str(source/'tools/nuclear_pid_cutoff/ComparePIDSelections.C')+'+('+','.join(json.dumps(str(p)) for p in (original,inputs/'cache.root',output,output/'response_curves.tsv',output/'tof_windows.tsv',output/'plot_settings.tsv'))+','+json.dumps(cfg['originalM2Histogram'])+')'
    command=['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',call];write(output/'command.json',command)
    with (output/'root.log').open('x') as log:result=subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=cfg['executionTimeoutSeconds'])
    if result.returncode or 'PID_SELECTION_COMPARISON_PASS' not in (output/'root.log').read_text():raise RuntimeError('See ROOT rendering log '+str(output/'root.log'))
    if fingerprint(original)!=native['inputFingerprint'] or any(sha(Path(p))!=value for p,value in originals.items()):raise ValueError('Source or original changed.')
    products={str(p.relative_to(output)):sha(p) for p in output.rglob('*') if p.is_file() and p.suffix not in ('.so','.d','.pcm')}
    write(output/'completed.json',dict(executionComplete=True,productionEnabled=False,noFit=True,originalFingerprintPreserved=True,protectedSourcesVerified=True,products=products))
    print('PID_SELECTION_COMPARISON_COMPLETE',output)

if __name__=='__main__':main()
