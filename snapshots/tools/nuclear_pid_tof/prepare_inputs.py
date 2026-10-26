#!/usr/bin/env python3
"""Prepare immutable 2021 P24iy nuclear-TOF native inputs. No fit or PID writes."""
import argparse
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
KEYS = {'ALL':'true/hDedxP', 'TAG_D':'true/hDedxP_d_m2',
        'TAG_T':'true/hDedxP_t_m2', 'TAG_HE3':'true/hDedxP_3He_m2'}
SETTING_KEYS = set('''expectedXBins expectedYBins expectedXMin expectedXMax expectedYMin expectedYMax
expectedXWidth expectedYWidth axisRelativeTolerance varianceRelativeTolerance displayDedxRebinFactor
overviewRigidityMin overviewRigidityMax atlasRigidityMin atlasRigidityMax coarseRigidityStep
diagnosticRigidityMin diagnosticRigidityMax fullDedxMin fullDedxMax supplementalDedxMin supplementalDedxMax
modeSmoothingHalfBins modeMomentHalfWindow logMinimum logMaximumFactor linearMaximumFactor axisMaxDigits
mapCanvasWidth mapCanvasHeight atlasCanvasWidth atlasCanvasHeight mapMarginLeft mapMarginRight
mapMarginBottom mapMarginTop projectionMarginLeft projectionMarginRight projectionMarginBottom
projectionMarginTop mapAxisTitleSize mapAxisLabelSize projectionAxisTitleSize projectionAxisLabelSize
xTitleOffset yTitleOffset headerTextSize captionTextSize markerSize'''.split())

def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024*1024), b''):
            h.update(chunk)
    return h.hexdigest()

def fingerprint(path):
    s = path.stat()
    return dict(path=str(path),size=s.st_size,mtime_ns=s.st_mtime_ns,device=s.st_dev,inode=s.st_ino)

def save(path, value):
    path.write_text(json.dumps(value,ensure_ascii=False,indent=2,allow_nan=False)+'\n')

def checked_path(value):
    p = Path(value)
    p = (p if p.is_absolute() else PROJECT/p).resolve(strict=True)
    if PROJECT not in p.parents:
        raise ValueError('Input/config must belong to this project: '+str(p))
    return p

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input-root','mainconf','settings','output-dir'):
        parser.add_argument('--'+name,type=Path,required=True)
    args = parser.parse_args()
    settings_path = args.settings.resolve(strict=True)
    cfg = yaml.safe_load(settings_path.read_text())
    if set(cfg) != {'schemaVersion','inputRoot','sourceSummary','tofNuclearConfig','executionTimeoutSeconds','protectedRelativePaths','settings'} or cfg['schemaVersion'] != 1:
        raise ValueError('Unsupported explicit INPUT settings schema.')
    if set(cfg['settings']) != SETTING_KEYS:
        raise ValueError('Missing or unknown explicit numeric INPUT settings.')
    source = args.input_root.resolve(strict=True)
    if source != checked_path(cfg['inputRoot']) or '2021P24iy_partial_merge' not in source.name:
        raise ValueError('Only the fixed 3p85 2021 P24iy source is authorized.')
    initial = fingerprint(source)
    state_path = checked_path(cfg['sourceSummary'])
    state = json.loads(state_path.read_text())
    if state['phase'] != 'complete' or state['provenance']['input'] != initial:
        raise ValueError('Source fingerprint differs from verified ROOT state.')
    mainconf = args.mainconf.resolve(strict=True)
    runtime = yaml.safe_load(mainconf.read_text())
    info_path = checked_path(Path('config')/runtime['analysis'])
    info = yaml.safe_load(info_path.read_text())
    tags = info['starTag']
    if (tags['collision'],tags['year'],tags['productionTag'],tags['libraryTag']) != ('auau3p85','year2021','P24iy','SL24y'):
        raise ValueError('Explicit runtime must be 3p85 / 2021 P24iy / SL24y ROOT5.')
    settings = dict(cfg['settings'],schemaVersion=cfg['schemaVersion'])
    if any(isinstance(v,bool) or not isinstance(v,(int,float)) or not math.isfinite(v) for v in settings.values()):
        raise ValueError('All macro settings must be explicit finite numeric fields.')
    if not isinstance(cfg['executionTimeoutSeconds'],int) or cfg['executionTimeoutSeconds']<=0:
        raise ValueError('Invalid execution timeout.')
    nuclear_path = checked_path(cfg['tofNuclearConfig'])
    nuclear = yaml.safe_load(nuclear_path.read_text())
    for name,key in (('D','Deuteron'),('T','Triton'),('HE3','He3')):
        mean,width = nuclear['nuclearTofMean'+key],nuclear['nuclearTofSigma'+key]
        settings['tag'+name+'Mass2Min'] = mean-nuclear['m2SigmaCut']*width
        settings['tag'+name+'Mass2Max'] = mean+nuclear['m2SigmaCut']*width
    ranges = sorted((settings['tag'+name+'Mass2Min'],settings['tag'+name+'Mass2Max']) for name in ('D','T','HE3'))
    if any(lo>=hi for lo,hi in ranges) or any(ranges[i][1]>=ranges[i+1][0] for i in range(2)):
        raise ValueError('Original TOF tag windows must be disjoint.')
    protected_paths = [checked_path(p) for p in cfg['protectedRelativePaths']]
    protected_paths += [settings_path,mainconf,info_path,state_path,nuclear_path,Path(__file__).resolve(),HERE/'ExtractTofInputs.C']
    protected = {str(p):digest(p) for p in protected_paths}
    # Verify frozen old fits against their completion manifest, not just a fresh hash.
    old = PROJECT/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/he4_purity_20261013/refined_lowrigidity_from040'
    done = json.loads((old/'completed.json').read_text())
    for p in protected_paths:
        if p.parent==old and p.name in done['products'] and digest(p)!=done['products'][p.name]:
            raise ValueError('Old fit/PDF changed before extraction: '+str(p))
    output = args.output_dir.resolve()
    allowed = PROJECT/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014'
    if output.exists() or allowed not in output.parents or source in output.parents:
        raise ValueError('Choose a nonexistent child of the new nuclear_pid_tof output folder.')
    output.mkdir(parents=True)
    snapshot = output/'source'
    snapshot.mkdir()
    for p in (Path(__file__).resolve(),HERE/'ExtractTofInputs.C',settings_path,mainconf,info_path,nuclear_path,PROJECT/'script/singularity_root_batch.sh'):
        shutil.copy2(p,snapshot/p.name)
    save(output/'context.json',dict(schemaVersion=1,inputFingerprint=initial,runtimeMainconf=str(mainconf),
        inputHistogramKeys=KEYS,cacheKeys=list(KEYS)+['REMAINING'],productionEnabled=False,noFit=True,
        originalCountPolicy='Preserve TH2F stored bin contents without integer rounding or Entries renormalization.',
        variancePolicy='Original Sumw2 absent: explicit sqrt(N) Poisson-type errors; correlated REMAINING variance is remaining N, NOT Var(ALL)+sum Var(TAG).',
        protectedPreRun=protected,settings=settings,
        limitations=['TOF tags are conditional samples with efficiency/acceptance/impurity biases.',
          'D tag accepts D/He4-compatible m/q; tags are not pure-isotope truth.',
          'Coarse modes and conditional moments are diagnostics, not Gaussian calibration or purity.',
          'Charge-inclusive, centrality-integrated QA; source partial merge is not full 23B coverage.']))
    (output/'settings.tsv').write_text(''.join('%s\t%.17g\n'%(k,v) for k,v in settings.items()))
    save(output/'effective_settings.json',settings)
    call=str(snapshot/'ExtractTofInputs.C')+'++('+','.join(json.dumps(str(p)) for p in (source,output,output/'settings.tsv'))+')'
    command=['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q',call]
    save(output/'command.json',command)
    with (output/'root.log').open('w') as log:
        result=subprocess.run(command,cwd=PROJECT,stdout=log,stderr=subprocess.STDOUT,timeout=cfg['executionTimeoutSeconds'])
    if result.returncode or 'NUCLEAR_TOF_INPUTS_COMPLETE' not in (output/'root.log').read_text():
        raise RuntimeError('INPUT extraction failed; preserved root.log contains diagnostics.')
    if fingerprint(source)!=initial or any(digest(Path(p))!=expected for p,expected in protected.items()):
        raise ValueError('Original input or protected fit/PDF/source changed.')
    save(output/'completed.json',dict(executionComplete=True,noFit=True,productionEnabled=False,
        inputFingerprintPreserved=True,protectedPostRunVerified=True,
        products={p.name:digest(p) for p in output.iterdir() if p.is_file()},
        sourceSha256={p.name:digest(p) for p in snapshot.iterdir() if p.is_file()}))
    print('NUCLEAR_TOF_PREPARE_COMPLETE',output,flush=True)

if __name__=='__main__':
    main()
