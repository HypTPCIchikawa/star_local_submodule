#!/usr/bin/env python3
"""Display immutable TAG fits on full physical native/display axes; no fit/PID."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SPECIES = ('d','t','he3','he4')

def digest(path):
    h=hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda:stream.read(1048576),b''):
            h.update(block)
    return h.hexdigest()

def fingerprint(path):
    s=path.stat()
    return dict(path=str(path),size=s.st_size,mtime_ns=s.st_mtime_ns,device=s.st_dev,inode=s.st_ino)

def save(path,value):
    path.write_text(json.dumps(value,indent=2,allow_nan=False)+'\n')

def manifest(data,path):
    if data['schemaVersion']!=1 or data['fitStage']!='independent_tag_shapes' or data['productionEnabled']:
        raise ValueError('Only the explicit study-only TAG result schema is supported.')
    required=('index','sample','lo','hi','r','firstNativeRBin','lastNativeRBin','fitLo','fitHi',
              'firstFitDedxBin','lastFitDedxBin','fitCount','backgroundModel','numericalConverged','status')
    optional=('covarianceStatus','finite','pointAccepted','amplitudeConverged','atBoundary','chi2','conditionalNdf','dataDeviance')
    header=list(required)+list(optional)+['componentNames']
    for name in SPECIES:
        header += [name+'Known',name+'FreeMean',name+'FreeSigma',name+'IdentityValidated',name+'Active',name+'ZeroYield',name+'OutsideDetector',name+'Mean',name+'Sigma']
    seen=set()
    with path.open('w',newline='') as stream:
        writer=csv.writer(stream,delimiter='\t',lineterminator='\n');writer.writerow(header)
        for s in data['slices']:
            if any(key not in s for key in required) or s['index'] in seen:
                raise ValueError('Missing required TAG metadata or duplicate index.')
            seen.add(s['index'])
            empty=s.get('emptyFitWindow',False)
            if not empty and any(key not in s for key in optional):
                raise ValueError('Nonempty fit metadata cannot be approximated.')
            components=s['components']
            if set(components)-set(SPECIES):
                raise ValueError('Unknown saved component name.')
            names=','.join(components) or '-'
            row=[s[key] for key in required]+[s.get(key,-1) for key in optional]+[names]
            for name in SPECIES:
                c=components.get(name)
                row += ([1,c['freeMean'],c['freeSigma'],c['identityValidated'],c['active'],c['zeroYield'],c['centerOutsideDetector'],c['mean'],c['sigma']] if c else [0,-1,-1,-1,-1,-1,-1,-999,-999])
            values=[]
            for value in row:
                if isinstance(value,bool):value=int(value)
                if isinstance(value,float) and not math.isfinite(value):raise ValueError('Nonfinite manifest field.')
                values.append(value)
            writer.writerow(values)
    if not seen:raise ValueError('Empty TAG result.')
    return len(seen)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--study-folder',type=Path)
    p.add_argument('--fit-dir',type=Path)
    p.add_argument('--compile-only',action='store_true')
    for name in ('mainconf','settings','output-dir'):p.add_argument('--'+name,type=Path,required=True)
    args=p.parse_args()
    if args.compile_only and (args.study_folder or args.fit_dir):raise ValueError('Compile-only does not consume study data.')
    if not args.compile_only and not args.study_folder:raise ValueError('--study-folder is required for drawing.')
    settings_path=args.settings.resolve(strict=True);cfg=yaml.safe_load(settings_path.read_text())
    if set(cfg)!={'schemaVersion','executionTimeoutSeconds','display'} or cfg['schemaVersion']!=1:
        raise ValueError('Unsupported display-only settings schema.')
    if not isinstance(cfg['executionTimeoutSeconds'],int) or cfg['executionTimeoutSeconds']<=0:raise ValueError('Invalid explicit timeout.')
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in cfg['display'].values()):raise ValueError('Display settings must be numeric and finite.')
    mainconf=args.mainconf.resolve(strict=True);runtime=yaml.safe_load(mainconf.read_text())
    info_path=(PROJECT/'config'/runtime['analysis']).resolve(strict=True);info=yaml.safe_load(info_path.read_text());tags=info['starTag']
    if (tags['collision'],tags['year'],tags['productionTag'],tags['libraryTag'])!=('auau3p85','year2021','P24iy','SL24y'):
        raise ValueError('Explicit ROOT5 2021P24iy 3p85 runtime required.')
    originals=[Path(__file__).resolve(),HERE/'PlotNuclearStudy.C',settings_path,mainconf,info_path,PROJECT/'script/singularity_root_batch.sh']
    protected={str(x):digest(x) for x in originals};initial=None;study=None;phase=None
    if not args.compile_only:
        study=args.study_folder.resolve(strict=True);phase=(args.fit_dir or study).resolve(strict=True)
        if phase!=study and study not in phase.parents:raise ValueError('Fit phase must belong to the completed study.')
        context_path,completion_path=study/'context.json',study/'completed.json'
        context=json.loads(context_path.read_text());completion=json.loads(completion_path.read_text())
        if context['productionEnabled'] or completion['productionEnabled'] or not completion['executionComplete']:
            raise ValueError('Study must be authoritative complete and unadopted.')
        if mainconf!=Path(context['runtimeMainconf']).resolve(strict=True):raise ValueError('Runtime differs from completed study.')
        initial=context['inputFingerprint']
        if fingerprint(Path(initial['path']))!=initial:raise ValueError('Original ROOT changed.')
        for key in ('protected','snapshots'):
            for encoded,expected in context[key].items():
                if encoded in protected and protected[encoded]!=expected:raise ValueError('Conflicting source protection.')
                protected[encoded]=expected
        for relative,expected in completion['products'].items():
            target=(study/relative).resolve(strict=True)
            if study not in target.parents:raise ValueError('Product escapes study.')
            protected[str(target)]=expected
        protected[str(context_path)]=digest(context_path);protected[str(completion_path)]=digest(completion_path)
        for name in ('tag_fits.root','tag_fits.json','slices.tsv'):
            if str(phase/name) not in protected:raise ValueError('Required fit product absent from authoritative manifest.')
    def preserve():
        if initial and fingerprint(Path(initial['path']))!=initial:raise ValueError('Original ROOT fingerprint changed.')
        for encoded,expected in protected.items():
            if digest(Path(encoded))!=expected:raise ValueError('Protected source/product changed: '+encoded)
    preserve()
    output=args.output_dir.resolve();allowed=PROJECT/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014'
    if output.exists() or allowed not in output.parents or (study and study in output.parents):raise ValueError('Choose a new sibling display output, never overwrite source study.')
    output.mkdir(parents=True);snapshot=output/'source';snapshot.mkdir()
    for source in originals:shutil.copy2(source,snapshot/source.name)
    display=dict(cfg['display']);n=0
    if not args.compile_only:
        data=json.loads((phase/'tag_fits.json').read_text());n=manifest(data,output/'plot_manifest.tsv')
        display['expectedSlices']=n
    (output/'settings.tsv').write_text(''.join('%s\t%.17g\n'%(key,value) for key,value in display.items()))
    save(output/'context.json',dict(schemaVersion=1,noFit=True,productionEnabled=False,compileOnly=args.compile_only,
         sourceStudy=str(study) if study else None,sourcePhase=str(phase) if phase else None,inputFingerprint=initial,
         protectedPreRun=protected,display=display,warning='Nominal isotope components are unvalidated; plots do not adopt calibration.'))
    macro=snapshot/'PlotNuclearStudy.C'
    if args.compile_only:call=str(macro)+'++';marker='NUCLEAR_STUDY_PLOT_COMPILE_OK'
    else:
        call=str(macro)+'++('+','.join(json.dumps(str(x)) for x in (phase/'tag_fits.root',output/'plot_manifest.tsv',output/'settings.tsv',output,'TAG_REFERENCE'))+')'
        marker='NUCLEAR_STUDY_PLOT_COMPLETE'
    command=['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',call];save(output/'command.json',command)
    with (output/'root.log').open('w') as log:result=subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=cfg['executionTimeoutSeconds'])
    if result.returncode or marker not in (output/'root.log').read_text():raise RuntimeError('Plot compile/draw failed; root.log retained.')
    preserve()
    save(output/'completed.json',dict(executionComplete=True,noFit=True,productionEnabled=False,compileOnly=args.compile_only,
        protectedPostRunVerified=True,inputFingerprintPreserved=True,
        products={x.name:digest(x) for x in output.iterdir() if x.is_file()},
        sourceSha256={x.name:digest(x) for x in snapshot.iterdir() if x.is_file()}))
    print('NUCLEAR_STUDY_RENDER_COMPLETE',output,flush=True)

if __name__=='__main__':main()
