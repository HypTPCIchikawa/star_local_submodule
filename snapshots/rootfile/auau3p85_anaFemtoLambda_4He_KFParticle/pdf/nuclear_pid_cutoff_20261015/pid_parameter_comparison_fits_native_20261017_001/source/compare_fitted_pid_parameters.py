#!/usr/bin/env python3
"""Overlay measured projected fit parameters onto immutable old/new PID functions."""
import argparse
import csv
import importlib.util
import math
from pathlib import Path
import shutil
import statistics
import subprocess
import compare_pid_selections as base
import yaml

HERE=Path(__file__).resolve().parent
PROJECT=base.PROJECT
LABELS=('gaussian_linear_nominal','gaussian_exponential_nominal','lognormal_linear_nominal','lognormal_exponential_nominal')

def checked_artifact(folder,name):
    done=base.read(folder/'completed.json');path=folder/name
    candidates=[]
    for section in ('products','artifacts','inputs'):
        if name in done.get(section,{}):
            entry=done[section][name];candidates.append(entry if isinstance(entry,str) else entry['sha256'])
    if not candidates or len(set(candidates))!=1:raise ValueError('Missing/contradictory completed input hash '+str(path))
    expected=candidates[0]
    if base.sha(path)!=expected:raise ValueError('Changed completed input '+str(path))
    return path

def number(row,key):
    value=row.get(key)
    try:return float(value)
    except (ValueError,TypeError):return float('nan')

def flag(row,key):return number(row,key)>.5

def snapshot_profiles(folder):
    source=folder/'source';provider_path=source/'tools/nuclear_pid_tof/runtime_selection.py'
    rel=str(provider_path.relative_to(folder));checked_artifact(folder,rel)
    spec=importlib.util.spec_from_file_location('frozen_fit_comparison',str(provider_path));provider=importlib.util.module_from_spec(spec);spec.loader.exec_module(provider)
    header=source/'StMaker/common/NuclearIdDeDxVsMom.h';checked_artifact(folder,str(header.relative_to(folder)))
    tables=provider.parse_legacy_tables(header);profiles={}
    for name,filename in (('oldCanonical','nuclearid_femtoLambda_auau3p85_calibrated.yaml'),('oldLegacy','nuclearid_femtoLambda_legacy.yaml'),('new','nuclearid_femtoLambda_auau3p85_purity50_20261015.yaml')):
        path=source/'config/cuts/nuclearid'/filename;checked_artifact(folder,str(path.relative_to(folder)));cuts=yaml.safe_load(path.read_text());calibration=None
        if cuts.get('dedxCalibration',{}).get('enabled'):
            p=source/'config'/cuts['dedxCalibration']['calibrationFile'];checked_artifact(folder,str(p.relative_to(folder)));calibration=yaml.safe_load(p.read_text())
        profiles[name]=dict(tables=tables,calibration=calibration,cuts=cuts)
    return provider,profiles

def enrich(row,region,provider,profiles):
    if row['label'] not in LABELS:return None
    s=base.SPECIES.index(row['species']);r=number(row,'r');m=number(row,'mean');sig=number(row,'sigma');em=number(row,'meanError');es=number(row,'sigmaError')
    current={}
    for name in ('oldCanonical','oldLegacy','new'):
        mu,width,_=provider.current_parameters(profiles[name],s,r);current[name+'Mean']=mu;current[name+'Sigma']=width
    valid=all(math.isfinite(v) for v in (m,sig,em,es,r)) and sig>0 and em>=0 and es>=0
    plotted=valid and flag(row,'numericalValid') and flag(row,'fitResponseQuality')
    pull=(current['newMean']-m)/sig if valid else None
    ratio=current['newSigma']/sig if valid else None
    return dict(region=region,label=row['label'],species=row['species'],rLo=number(row,'rLo'),rHi=number(row,'rHi'),r=r,mean=m if math.isfinite(m) else None,sigma=sig if math.isfinite(sig) else None,meanError=em if math.isfinite(em) else None,sigmaError=es if math.isfinite(es) else None,
        numericalValid=int(flag(row,'numericalValid')),fitResponseQuality=int(flag(row,'fitResponseQuality')),targetFreeSupported=int(flag(row,'targetFreeSupported')),plotted=int(plotted),centerPull=pull,widthRatio=ratio,
        centerPullDiagonalError=math.hypot(em/sig,(current['newMean']-m)*es/(sig*sig)) if valid else None,widthRatioError=current['newSigma']*es/(sig*sig) if valid else None,**current)

def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--settings',required=True,type=Path);parser.add_argument('--output-dir',required=True,type=Path);parser.add_argument('--high-summary',type=Path);parser.add_argument('--coarse-summary',type=Path);args=parser.parse_args()
    settings=args.settings.resolve(strict=True);cfg=yaml.safe_load(settings.read_text())
    if cfg['schemaVersion']!=1:raise ValueError('Unsupported settings.')
    configured=base.resolve(cfg['configuredComparison']);low=base.resolve(cfg['lowViews']);high=base.resolve(args.high_summary or cfg['highSummary'])
    output=args.output_dir.resolve();allowed=PROJECT/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015'
    if output.exists() or output.parent!=allowed:raise ValueError('Require a new direct child of the PID study.')
    provider,profiles=snapshot_profiles(configured)
    paths=[('low_ALL',low,'inputs/native_parameters.tsv'),('high_native_TOF',high,'native_parameters.tsv')]
    if args.coarse_summary:paths.append(('high_coarse_TOF',base.resolve(args.coarse_summary),'native_parameters.tsv'))
    oldpdf=checked_artifact(configured,'pid_parameter_comparison.pdf');curves=checked_artifact(configured,'response_curves.tsv');frozen_macro=checked_artifact(configured,'source/tools/nuclear_pid_cutoff/ComparePIDSelections.C')
    inputs={str(oldpdf):base.sha(oldpdf),str(curves):base.sha(curves),str(configured/'completed.json'):base.sha(configured/'completed.json')};points=[];coverage=[]
    for name in base.read(configured/'completed.json')['products']:
        if name.startswith('source/') and Path(name).suffix in ('.py','.h','.yaml','.C'):
            checked_artifact(configured,name);inputs[str(configured/name)]=base.sha(configured/name)
    for region,folder,name in paths:
        path=checked_artifact(folder,name);inputs[str(path)]=base.sha(path);inputs[str(folder/'completed.json')]=base.sha(folder/'completed.json')
        with path.open() as stream:raw=list(csv.DictReader(stream,delimiter='\t'))
        enriched=[v for v in (enrich(r,region,provider,profiles) for r in raw) if v is not None]
        if not enriched:raise ValueError('No nominal fitted parameters in '+str(path))
        seen=set()
        for row in enriched:
            identity=(row['label'],row['species'],row['rLo'],row['rHi'])
            if identity in seen:raise ValueError('Duplicate fit point in a same-region/model.')
            seen.add(identity)
            if region=='low_ALL' and row['rHi']>2.00000001:raise ValueError('Low ALL sample crossed TOF boundary.')
            if region!='low_ALL' and row['rLo']<1.99999999:raise ValueError('High TOF sample below threshold.')
        for species in base.SPECIES:
            rows=[r for r in enriched if r['species']==species]
            if not rows:continue
            for label in LABELS:
                selected=[r for r in rows if r['label']==label]
                if not selected:raise ValueError('Missing declared nominal model '+species+' '+label)
                good=[r for r in selected if r['plotted']]
                coverage.append(dict(region=region,species=species,label=label,total=len(selected),plotted=len(good),identitySupported=sum(r['targetFreeSupported'] for r in good),medianCenterPull=statistics.median(r['centerPull'] for r in good) if good else None,medianWidthRatio=statistics.median(r['widthRatio'] for r in good) if good else None))
        points.extend(enriched)
    output.mkdir(parents=True);source=output/'source';source.mkdir()
    for path in (Path(__file__).resolve(),Path(base.__file__).resolve(),HERE/'DrawFittedPIDComparison.C',settings):shutil.copy2(path,source/path.name)
    shutil.copy2(configured/'completed.json',source/'configured_comparison_completed.json')
    shutil.copy2(frozen_macro,source/'ComparePIDSelections.C');shutil.copy2(curves,output/'configured_curves.tsv')
    for region,folder,name in paths:
        shutil.copy2(folder/name,source/(region+'_native_parameters.tsv'));shutil.copy2(folder/'completed.json',source/(region+'_completed.json'))
    keys=list(points[0]);base.tsv(output/'fit_points.tsv',keys,[[r[k] for k in keys] for r in points])
    base.tsv(output/'coverage.tsv',list(coverage[0]),[[r[k] for k in coverage[0]] for r in coverage])
    base.tsv(output/'plot_settings.tsv',['key','value'],[(k,cfg[k]) for k in ('canvasWidth','canvasHeight')])
    base.write(output/'context.json',dict(schemaVersion=1,productionEnabled=False,noNewFit=True,oldConfiguredBasis=str(configured),oldBasis='Immutable snapshot only; live mainconf never used as old reference.',fitInputs=inputs,regions=[r[0] for r in paths],lowSample='ALL before TOF selection',highSample='preTPC TOF-tagged histogram, He4 study uses D window',widthMeaning='Gaussian mean/std vs lognormal mode/local-curvature width; projected response, not adopted intrinsic calibration.',plotQuality='Numerical and fit-response precision screen; gray points fail target identity/reference screen. Neither status implies full purity support.',ratioError='Diagonal fit-parameter approximation; center-width covariance and configured-response uncertainty not included.',noOriginalRootRead=True,noM2CacheRead=True,coverage=coverage))
    mainconf=base.resolve(cfg['runtimeMainconf']);call=str(source/'DrawFittedPIDComparison.C')+'+('+','.join(base.json.dumps(str(p)) for p in (output/'fit_points.tsv',output/'configured_curves.tsv',output/'coverage.tsv',output/'plot_settings.tsv',output))+')'
    command=['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',call];base.write(output/'command.json',command)
    with (output/'root.log').open('x') as log:result=subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=cfg['executionTimeoutSeconds'])
    if result.returncode or 'FITTED_PID_COMPARISON_PASS' not in (output/'root.log').read_text():raise RuntimeError('ROOT rendering failed; see log.')
    merge=['gs','-q','-dBATCH','-dNOPAUSE','-sDEVICE=pdfwrite','-sOutputFile='+str(output/'new_old_pid_final.pdf'),str(oldpdf),str(output/'fitted_parameter_comparison.pdf')];base.write(output/'merge_command.json',merge)
    with (output/'merge.log').open('x') as log:subprocess.run(merge,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=cfg['executionTimeoutSeconds'])
    if any(base.sha(Path(p))!=value for p,value in inputs.items()):raise ValueError('Fit/source products changed during rendering.')
    base.write(output/'completed.json',dict(executionComplete=True,noNewFit=True,productionEnabled=False,protectedInputsVerified=True,oldPdfPreserved=base.sha(oldpdf)==inputs[str(oldpdf)],products={str(p.relative_to(output)):base.sha(p) for p in output.rglob('*') if p.is_file() and p.suffix not in ('.so','.d','.pcm')}))
    print('FITTED_PID_COMPARISON_COMPLETE',output)

if __name__=='__main__':main()
