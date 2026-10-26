#!/usr/bin/env python3
"""Run legacy synthetic regression and actual mainconf domains without event jobs."""
import argparse
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import time
import yaml
PROJECT=Path(__file__).resolve().parents[1]
def sha(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1024*1024),b''):h.update(block)
    return h.hexdigest()
def save(path,value):
    with path.open('x') as f:json.dump(value,f,indent=2,ensure_ascii=False,allow_nan=False);f.write('\n')
def resolve(value):
    p=Path(value);return (p if p.is_absolute() else PROJECT/p).resolve(strict=True)
def main():
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument('--mainconf',action='append',required=True,type=Path);parser.add_argument('--ranges',required=True,type=Path);parser.add_argument('--windows',required=True,type=Path);parser.add_argument('--output-dir',required=True,type=Path);args=parser.parse_args()
    mains=list(map(resolve,args.mainconf));ranges,windows=resolve(args.ranges),resolve(args.windows);out=args.output_dir.resolve()
    if out.exists():raise ValueError('Output must be new.')
    for p in (ranges,windows):
        done=json.loads((p.parent/'completed.json').read_text())
        if sha(p)!=done['products'][p.name]:raise ValueError('Expected range/window input changed.')
    sources=[Path(__file__).resolve(),PROJECT/'tests/test_femto_nuclear_pid_domain.C',PROJECT/'tests/run_femto_nuclear_pid_domain_test.C',PROJECT/'tests/test_femto_nuclear_pid_actual_domain.C',PROJECT/'tests/run_femto_nuclear_pid_actual_domain_test.C',PROJECT/'StMaker/common/FemtoLambdaLegacy.cxx',PROJECT/'StMaker/common/FemtoLambdaLegacy.h',PROJECT/'script/singularity_root_batch.sh',ranges,windows]
    for p in mains:
        doc=yaml.safe_load(p.read_text());sources.append(p)
        for ref in doc.values():sources.append(resolve(Path('config')/ref) if not Path(ref).is_absolute() else resolve(ref))
        info=yaml.safe_load(resolve(Path('config')/doc['analysis']).read_text())['starTag']
        if tuple(info[k] for k in ('collision','year','productionTag','libraryTag'))!=('auau3p85','year2021','P24iy','SL24y'):raise ValueError('Wrong production dataset/toolchain.')
        nuclear=yaml.safe_load(resolve(Path('config')/doc['nuclearid']).read_text())
        if nuclear.get('dedxCalibration',{}).get('enabled'):sources.append(resolve(Path('config')/nuclear['dedxCalibration']['calibrationFile']))
    for name in ('libStCommon.so','libStarAnaConfig.so','libStRefMultCorr.so'):sources.append(PROJECT/'lib'/name)
    protected={str(p):sha(p) for p in set(sources)};out.mkdir(parents=True);snap=out/'source';snap.mkdir()
    for p in set(sources):
        if p.suffix=='.so':continue
        dest=snap/p.relative_to(PROJECT);dest.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(p,dest)
        if sha(dest)!=protected[str(p)]:raise ValueError('Source changed during snapshot.')
    save(out/'context.json',dict(productionEnabled=False,noEventLoop=True,mainconfs=list(map(str,mains)),ranges=str(ranges),windows=str(windows),protectedSourcesAndLibraries=protected))
    runs=[('fixture_1359',mains[0],'tests/run_femto_nuclear_pid_domain_test.C',[mains[0],out/'fixture_1359'])]
    runs += [('actual_'+str(i),p,'tests/run_femto_nuclear_pid_actual_domain_test.C',[p,ranges,windows,out/('actual_'+str(i))]) for i,p in enumerate(mains)]
    records=[]
    for label,mainconf,macro,values in runs:
        call=macro+'('+','.join(json.dumps(str(v)) for v in values)+')';command=['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',call];save(out/(label+'_command.json'),command);started=time.monotonic()
        with (out/(label+'.log')).open('x') as log:result=subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=1800)
        text=(out/(label+'.log')).read_text();matches=re.findall(r'PASS (?:actual )?nuclear PID domain: (\d+) checks',text)
        passed=result.returncode==0 and bool(matches) and (label!='fixture_1359' or int(matches[-1])==1359)
        record=dict(label=label,mainconf=str(mainconf),exitCode=result.returncode,passed=passed,checks=int(matches[-1]) if matches else None,runtimeSeconds=time.monotonic()-started);records.append(record);save(out/(label+'_result.json'),record);print('DOMAIN_REGRESSION',label,'PASS' if passed else 'FAIL',record['checks'],flush=True)
        if not passed:raise RuntimeError('Domain regression failed: '+str(out/(label+'.log')))
    if any(sha(Path(p))!=v for p,v in protected.items()):raise ValueError('Source/library changed during regression.')
    save(out/'completed.json',dict(executionComplete=True,productionEnabled=False,sourceLibrariesVerified=True,noEventLoop=True,runs=records,totalChecks=sum(r['checks'] for r in records),products={str(p.relative_to(out)):sha(p) for p in out.rglob('*') if p.is_file()}))
    print('NUCLEAR_PRODUCTION_DOMAIN_REGRESSIONS_PASS',out,flush=True)
if __name__=='__main__':main()
