#!/usr/bin/env python3
"""Saved from040 local fits -> soft iterations -> study-only common raw fit."""
import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

from smoothing import normalized_basis, smooth_species, value

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def fingerprint(path):
    s = path.stat()
    return dict(path=str(path), size=s.st_size, mtime_ns=s.st_mtime_ns, device=s.st_dev, inode=s.st_ino)


def save(path, obj):
    path.write_text(json.dumps(obj, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def tsv(path, rows):
    path.write_text(''.join('\t'.join(format(v,'.17g') if isinstance(v,float) else str(v)
                                    for v in row)+'\n' for row in rows))


def probability(lo, hi, mean, sigma):
    a, b = (lo-mean)/(math.sqrt(2)*sigma), (hi-mean)/(math.sqrt(2)*sigma)
    if a >= 0:
        return (math.erfc(a)-math.erfc(b))/2
    if b <= 0:
        return (math.erfc(-b)-math.erfc(-a))/2
    return (math.erf(b)-math.erf(a))/2


def used_window(anchor, sigma, engine, width=.5, upper=120.):
    requested_low = max(0., anchor+engine['fitNSigmaMin']*sigma)
    requested_high = min(upper, anchor+engine['fitNSigmaMax']*sigma)
    first = max(1, int(math.ceil(requested_low/width+.5-1.e-10)))
    last = min(round(upper/width), int(math.ceil(requested_high/width+.5-1.e-10))-1)
    return (first-1)*width, last*width


def old_seeds(audit, engine, indices):
    rows, details = [], []
    width = audit['fitDedxBinWidth']
    for index in indices:
        item = audit['slices'][index]
        eligible = [(j,f) for j,f in enumerate(item['fits'])
                    if len(f.get('parameters',[])) >= 8 and 2 in f.get('neighborSpecies',[]) and
                    all(math.isfinite(v) for v in f['parameters']) and f['parameters'][2] > 0]
        if not eligible:
            raise ValueError('No saved Gaussian seed for slice '+str(index))
        chosen = item['chosen']
        j, fit = next(((j,f) for j,f in eligible if j == chosen),
                      min(eligible,key=lambda pair:pair[1]['deviance']))
        p, anchor, sigma = fit['parameters'], item['calibrationMean'], item['calibrationSigma']
        low, high = used_window(anchor, sigma, engine, width)
        components = {}
        for species, offset in [('he4',0)]+[(('d','t','he3')[s],5+3*k)
                                           for k,s in enumerate(fit['neighborSpecies'])]:
            mu, sd = anchor+sigma*p[offset+1], sigma*p[offset+2]
            full = p[offset]*sd*math.sqrt(2*math.pi)/width
            accepted = probability(low,high,mu,sd)
            components[species] = dict(mean=mu,sigma=sd,windowYield=full*accepted,
                                       oldHeight=p[offset],fullYield=full,windowAcceptance=accepted)
        u0, u1 = (low-anchor)/sigma, (high-anchor)/sigma
        if fit['model'] == 'exponential':
            slope = p[4]*(high-low)/sigma
            bg = p[3]*sigma/width*(u1-u0 if abs(p[4])<1.e-12 else
                                 (math.exp(p[4]*u1)-math.exp(p[4]*u0))/p[4])
        elif fit['model'] == 'linear':
            fit0, fit1 = engine['fitNSigmaMin'], min(engine['fitNSigmaMax'],(120.-anchor)/sigma)
            slope = 0.
            bg = sigma/width*(p[3]*(u1-u0)+(p[4]-p[3])/((fit1-fit0)*2)*
                              ((u1-fit0)**2-(u0-fit0)**2))
        else:
            raise ValueError('Unknown old background')
        t = components.get('t',dict(mean=0.,sigma=1.,windowYield=0.))
        rows.append([index,item['lo'],item['hi'],anchor,sigma,chosen,int(fit['valid']),int(fit['quality']),
                     fit['model'],int('t' in components),t['mean'],t['sigma'],t['windowYield'],
                     components['he3']['mean'],components['he3']['sigma'],components['he3']['windowYield'],
                     components['he4']['mean'],components['he4']['sigma'],components['he4']['windowYield'],
                     max(0.,bg),slope])
        details.append(dict(index=index,chosen=chosen,seedFitIndex=j,oldAccepted=chosen>=0 and fit['valid'],
                            oldQuality=fit['quality'],fitWindow=[low,high],components=components,
                            backgroundModel=fit['model'],backgroundWindowYield=bg,backgroundSlope=slope,
                            interpretation='Saved fit initialization, NOT truth or a fit constraint.'))
    return rows, details


def validate(cfg):
    sample, engine, sm = cfg['sample'],cfg['engine'],cfg['smoothing']
    if cfg['schemaVersion'] != 1 or (sample['collision'],sample['year'],sample['production']) != ('auau3p85',2021,'P24iy'):
        raise ValueError('Only 3.85 GeV / 2021 P24iy')
    if not (engine['rigidityMin'] == .4 and engine['rigidityMax'] == 2.5 and engine['dedxRebinFactor'] == 5):
        raise ValueError('Study must use the frozen native from040 domain and .5 keV/cm bins')
    if engine['rQuadratureOrder'] not in (1,3,5) or not 0 < engine['minimumSigma'] < engine['maximumSigma']:
        raise ValueError('Invalid width or R averaging policy')
    if (sm['weightPolicy'] != 'component_squared_over_total_counts' or not sm['penaltyStrengths'] or
            any(v <= 0 for v in sm['penaltyStrengths']) or sm['iterationsPerStrength'] < 1):
        raise ValueError('Invalid declared iterative policy')
    for block in (engine,sm,cfg['display']):
        for key,v in block.items():
            if isinstance(v,(int,float)) and not isinstance(v,bool) and not math.isfinite(v):
                raise ValueError('Nonfinite setting '+key)


def seed_update(previous, seed_rows):
    result = []
    by_index = {row['index']:row for row in previous['slices']}
    for old in seed_rows:
        row = list(old)
        data = by_index[row[0]]
        for key,offset in [('t',10),('he3',13),('he4',16)]:
            comp = data.get(key,data.get('components',{}).get(key))
            if comp is None:
                if key == 't' and not row[9]:
                    continue
                raise ValueError('Missing saved local seed '+key)
            row[offset:offset+3] = [comp['mean'],comp['sigma'],comp['windowYield']]
        if 'backgroundSlope' in data:
            row[20] = data['backgroundSlope']
        if 'backgroundWindowYield' in data:
            row[19] = data['backgroundWindowYield']
        result.append(row)
    return result


def common(rows, cfg):
    species = {key:smooth_species(rows,key,cfg['smoothing'],
               [cfg['engine']['rigidityMin'],cfg['engine']['rigidityMax']]) for key in ('he3','he4')}
    coefficients = []
    for key in ('he3','he4'):
        coefficients += species[key]['meanCoefficients']+species[key]['sigmaCoefficients']
    return coefficients,species


def penalty(rows, coefficients, fixed):
    result = 0.
    for row in rows:
        for offset,key in [(0,'he3'),(7,'he4')]:
            comp, info = row[key],fixed[row['index']][key]
            x = math.log10(row['r'])
            result += info['effectiveCount']*((comp['mean']-value(coefficients[offset:offset+4],x))**2 +
                       (comp['sigma']-value(coefficients[offset+4:offset+7],x))**2)/info['penaltyScale']**2
    return result


def candidate_yaml(coefficients,cfg,mode):
    return dict(schemaVersion=1,functionalForm='log10_polynomial',rigidityUnit='GeV/c',dedxUnit='keV/cm',
        outsideRangePolicy='legacy',productionEnabled=False,candidateOnly=True,phase=mode,
        widthMeaning='Intrinsic Gaussian width after raw R-bin-averaged joint refit' if mode=='common' else
                     'Folded projected width; initialization only, NOT intrinsic detector width',
        species={key:dict(domain=[cfg['engine']['rigidityMin'],cfg['engine']['rigidityMax']],
                         meanCoefficients=coefficients[offset:offset+4],sigmaCoefficients=coefficients[offset+4:offset+7])
                 for key,offset in [('he3',0),('he4',7)]})


def run_command(command, log, timeout, marker):
    with log.open('w') as handle:
        handle.write(json.dumps(command)+'\n');handle.flush()
        proc = subprocess.run(command,cwd=PROJECT,stdout=handle,stderr=subprocess.STDOUT,timeout=timeout)
    content = log.read_text(errors='replace')
    if proc.returncode or marker not in content or any(s in content for s in
            ['Error in <ACLiC>','segmentation violation','fatal error:','Error: Function']):
        raise RuntimeError('ROOT phase failed; see '+str(log))


def result_summary(result):
    rows=result['slices']
    return dict(dataDeviance=result['dataDeviance'],penaltyLoss=result['penaltyLoss'],
        totalLoss=result['totalLoss'],numericalConverged=result['numericalConverged'],
        nSlices=len(rows),finiteSlices=sum(r['finite'] for r in rows),
        amplitudeConvergedSlices=sum(r['ampConverged'] for r in rows),
        localStatusZeroSlices=sum(r['status']==0 for r in rows),
        numericalValidSlices=sum(r['valid'] for r in rows),
        boundarySlices=sum(r['atBoundary'] for r in rows),
        labelSwappedSlices=sum(r['labelSwapped'] for r in rows),
        common=result['common'],
        interpretation='Numerical convergence only; NOT isotope identification/model adequacy/production validation')


def renderer_manifest(result, path, display):
    header=['index','lo','hi','r','anchorMean','anchorSigma','fitLoDedx','fitHiDedx',
            'firstFitBin','lastFitBin','stageMode','model','valid','quality','status','covarianceStatus',
            'atBoundary','chi2','deviance','penalty','ndf','rawLocalMaximum','oldAccepted',
            'tShapeState','he3ShapeState','he4ShapeState']
    mode={'free':'local_free','regularized':'regularized','common':'shared_fixed'}[result['phase']]
    rows=[header]
    for data in result['slices']:
        states=[('unidentified' if data[key]['shapeState']=='fixed_unidentified' else data[key]['shapeState'])
                for key in ('t','he3','he4')]
        residual_poor=data['ndf']>0 and data['chi2']/data['ndf']>display['maximumReferenceChi2Ndf']
        # The core's quality flag is NUMERICAL; never relabel it scientific quality.
        quality=0 if residual_poor or not data['valid'] or data['atBoundary'] else -1
        peak=data['rawSeparatePeak']
        rows.append([data[k] for k in ('index','lo','hi','r','anchorMean','anchorSigma')]+
            [data['fitLo'],data['fitHi'],data['firstFitBin'],data['lastFitBin'],mode,data['backgroundModel'],
             int(data['valid']),quality,data['status'],data['covarianceStatus'],int(data['atBoundary']),
             data['chi2'],data['dataDeviance'],data['penaltyLoss'],data['ndf'],
             -1 if peak is None else int(peak),int(data['oldChosen']>=0 and data['oldValid'])]+states)
    tsv(path,rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',type=Path,default=HERE/'iterative_settings.yaml')
    parser.add_argument('--mainconf',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    parser.add_argument('--indices',help='Comma-separated cached slice indices; provenance records subset.')
    parser.add_argument('--stop-after-free',action='store_true')
    parser.add_argument('--prepare-only',action='store_true')
    parser.add_argument('--no-render',action='store_true')
    args = parser.parse_args()
    setting_path = args.settings.resolve(strict=True)
    cfg = yaml.safe_load(setting_path.read_text());validate(cfg)
    sample = cfg['sample'];mainconf=args.mainconf.resolve(strict=True)
    if mainconf != (PROJECT/sample['mainconf']).resolve(strict=True):
        raise ValueError('Explicit runtime mainconf differs')
    runtime = yaml.safe_load(mainconf.read_text())
    analysis_path = (PROJECT/'config'/runtime['analysis']).resolve(strict=True)
    analysis = yaml.safe_load(analysis_path.read_text())
    tags = analysis.get('starTag',{})
    if (tags.get('libraryTag'),tags.get('year'),tags.get('productionTag'),tags.get('collision'),
            analysis.get('analysis',{}).get('mode')) != ('SL24y','year2021','P24iy','auau3p85','fxtmult'):
        raise ValueError('Runtime identity is not SL24y/2021P24iy/3p85/FXT')
    old = (PROJECT/sample['oldStudy']).resolve(strict=True)
    audit_path,cache = old/'purity.json',old/'he4_fits.root'
    completed = json.loads((old/'completed.json').read_text())
    if completed.get('phase') != 'complete':
        raise ValueError('Only completed frozen source fits')
    for path in (audit_path,cache,old/'context.json'):
        if digest(path) != completed['products'][path.name]:
            raise ValueError('Old source hash mismatch '+str(path))
    audit=json.loads(audit_path.read_text())
    old_settings=yaml.safe_load((old/'settings.yaml').read_text())
    # The linear endpoints belong to the ORIGINAL fit window.  This study
    # holds that window fixed rather than silently reinterpreting endpoints.
    if any(cfg['engine'][key]!=old_settings[key] for key in ('fitNSigmaMin','fitNSigmaMax')):
        raise ValueError('New fit range differs from frozen source; old linear endpoints require an explicit model transform')
    source=(PROJECT/sample['input']).resolve(strict=True)
    summary=(PROJECT/sample['sourceSummary']).resolve(strict=True)
    initial=fingerprint(source)
    if Path(audit['inputRoot']).resolve() != source or json.loads(summary.read_text())['provenance']['input'] != initial:
        raise ValueError('Source is not verified same 2021P24iy original ROOT')
    if audit['dedxRebinFactor'] != 5 or audit['fitDedxBinWidth'] != .5 or len(audit['slices']) != 210:
        raise ValueError('Frozen source binning differs')
    indices=list(range(210)) if not args.indices else [int(s) for s in args.indices.split(',')]
    if len(indices)!=len(set(indices)) or any(not 0<=v<210 for v in indices) or not indices:
        raise ValueError('Invalid cached indices')
    if not args.stop_after_free and len(indices)<cfg['smoothing']['minimumSmoothingSlices']:
        raise ValueError('Insufficient common-function support')
    output=args.output_dir.resolve();allowed=(PROJECT/sample['outputBase']).resolve()
    if output.exists() or allowed not in output.parents:
        raise ValueError('Choose a nonexistent separate child of '+str(allowed))
    rows,seed_audit=old_seeds(audit,cfg['engine'],indices)
    output.mkdir(parents=True);(output/'logs').mkdir()
    files=[Path(__file__),HERE/'smoothing.py',HERE/'IterativeHeliumFit.C',HERE/'RenderIterativeHelium.C',
           HERE/'DrawIterativeSummary.C',
           HERE.parent/'global_logpoly/GlobalHeliumFit.C',PROJECT/'StMaker/common/NuclearIdDeDxVsMom.h',
           PROJECT/'script/singularity_root_batch.sh',setting_path,mainconf,analysis_path]
    old_figures=[old/'views/he4_low040_components.pdf',old/'views_full_dedx/he4_low040_components_full_dedx.pdf',
                 old.parent/'calibrated_selection/purity_calibrated_windows.pdf',
                 old.parent/'refined_lowrigidity/he4_fits_linear.pdf']
    protected={str(path):digest(path) for path in files+[audit_path,cache,old/'completed.json',old/'context.json',
               old/'settings.yaml',summary]+old_figures}
    snapshots={}
    for path in files:
        target=output/'source'/path.resolve().relative_to(PROJECT)
        target.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(path,target)
        snapshots[str(target)]=digest(target)
    def preserve():
        if fingerprint(source)!=initial:
            raise ValueError('Original ROOT fingerprint changed')
        for p,h in {**protected,**snapshots}.items():
            if digest(Path(p))!=h:
                raise ValueError('Protected source/snapshot changed '+p)
    save(output/'seed_audit.json',seed_audit);tsv(output/'initial_seeds.tsv',rows)
    (output/'settings.yaml').write_text(yaml.safe_dump(cfg,sort_keys=False))
    context=dict(schemaVersion=1,source=initial,cache=str(cache),protected=protected,snapshots=snapshots,
                 settings=cfg,indices=indices,productionEnabled=False,createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat())
    save(output/'context.json',context)
    if args.prepare_only:
        print('ITERATIVE_PREPARED',output,flush=True);return
    macro=output/'source/tools/he4_purity/iterative_logpoly/IterativeHeliumFit.C'
    renderer=output/'source/tools/he4_purity/iterative_logpoly/RenderIterativeHelium.C'
    wrapper=PROJECT/'script/singularity_root_batch.sh'
    history=[];fixed={}
    def execute(name,phase,seed_rows,coefficients=None,strength=0.):
        preserve();folder=output/name;folder.mkdir()
        engine=dict(cfg['engine'],phase=phase,penaltyStrength=strength)
        tsv(folder/'settings.tsv',engine.items());tsv(folder/'seeds.tsv',seed_rows)
        target_rows=([['none']] if coefficients is None else [['coefficients']+coefficients]+
            [['slice',index,fixed[index]['he3']['effectiveCount'],fixed[index]['he3']['penaltyScale'],
              fixed[index]['he3']['penaltyScale'],fixed[index]['he4']['effectiveCount'],
              fixed[index]['he4']['penaltyScale'],fixed[index]['he4']['penaltyScale']] for index in indices])
        tsv(folder/'targets.tsv',target_rows)
        save(folder/'context.json',dict(phase=phase,name=name,inputCache=str(cache),inputCacheSha256=protected[str(cache)],
            settingsSha256=digest(folder/'settings.tsv'),seedsSha256=digest(folder/'seeds.tsv'),
            targetsSha256=digest(folder/'targets.tsv'),parentContext=str(output/'context.json'),
            productionEnabled=False))
        call=str(macro)+'+('+','.join(json.dumps(str(p)) for p in
                                    (cache,folder,folder/'settings.tsv',folder/'seeds.tsv',folder/'targets.tsv'))+')'
        command=['bash',str(wrapper),str(mainconf),'root4star','-l','-b','-q',call]
        print('START',name,flush=True);save(output/'progress.json',dict(phase=name,state='running',history=history))
        run_command(command,folder/'fit.log',cfg['executionTimeoutSeconds'],'ITERATIVE_HELIUM_COMPLETE')
        result=json.loads((folder/'fit_result.json').read_text(),parse_constant=lambda v: (_ for _ in ()).throw(ValueError(v)))
        if [r['index'] for r in result['slices']] != indices:
            raise ValueError('Output slice order changed')
        summary=result_summary(result)
        save(folder/'summary.json',summary)
        if not args.no_render:
            views=folder/'views';views.mkdir()
            renderer_manifest(result,views/'slices.tsv',cfg['display'])
            tsv(views/'settings.tsv',dict(cfg['display'],expectedSlices=len(indices),
                allowSparseSlices=int(indices!=list(range(210)))).items())
            call=str(renderer)+'+('+','.join(json.dumps(str(p)) for p in
                    (folder/'phase.root',views/'slices.tsv',views/'settings.tsv',views,name))+')'
            run_command(['bash',str(wrapper),str(mainconf),'root4star','-l','-b','-q',call],
                        views/'render.log',cfg['executionTimeoutSeconds'],'ITERATIVE_HELIUM_RENDER_COMPLETE')
        preserve()
        history.append(dict(phase=name,mode=phase,penaltyStrength=strength,summary=summary))
        save(output/'history.json',history)
        print('FINISH',name,summary,flush=True)
        return result
    previous=execute('00_local_free','free',rows)
    if not args.stop_after_free:
        fixed={r['index']:{key:dict(r[key],penaltyScale=max(r[key]['sigma'],cfg['smoothing']['minimumPenaltyScale']))
                            for key in ('he3','he4')} for r in previous['slices']}
        influence=[]
        for row in previous['slices']:
            for key in ('he3','he4'):
                item=fixed[row['index']][key]
                eligible=row['finite'] and item['windowYield']>cfg['smoothing']['minimumWindowYield']
                if not eligible:
                    item['effectiveCount']=0.
                influence.append(dict(index=row['index'],r=row['r'],species=key,eligible=eligible,
                    effectiveCount=item['effectiveCount'],penaltyScale=item['penaltyScale'],
                    smoothingWeight=item['effectiveCount']/item['penaltyScale']**2,
                    scope='Fixed free-phase numerical/window support, not nuclear identity or covariance precision'))
        save(output/'smoothing_influence.json',influence)
        def smooth_result(result):
            smoothing_rows=[]
            for row in result['slices']:
                if not row['finite']:
                    raise ValueError('Nonfinite local model; cannot smooth a failed phase')
                copied=dict(row)
                for key in ('he3','he4'):
                    copied[key]=dict(row[key],effectiveCount=fixed[row['index']][key]['effectiveCount'],
                                     penaltyScale=fixed[row['index']][key]['penaltyScale'],
                                     windowYield=fixed[row['index']][key]['windowYield'])
                smoothing_rows.append(copied)
            return common(smoothing_rows,cfg)
        coefficients,fit_audit=smooth_result(previous)
        save(output/'initial_common.json',dict(coefficients=coefficients,audit=fit_audit))
        number=1
        for strength in cfg['smoothing']['penaltyStrengths']:
            for iteration in range(cfg['smoothing']['iterationsPerStrength']):
                name='%02d_regularized_lambda%s_iter%d'%(number,str(strength).replace('.','p'),iteration+1)
                next_result=execute(name,'regularized',seed_update(previous,rows),coefficients,strength)
                proposed,fit_audit=smooth_result(next_result)
                old_penalty=penalty(next_result['slices'],coefficients,fixed)
                adopted=coefficients;fraction=0.
                for trial_fraction in cfg['smoothing']['dampingFractions']:
                    trial=[a+trial_fraction*(b-a) for a,b in zip(coefficients,proposed)]
                    if penalty(next_result['slices'],trial,fixed)<=old_penalty*(1+cfg['smoothing']['objectiveRelativeTolerance']):
                        adopted=trial;fraction=trial_fraction;break
                drift=max(abs(a-b)/max(1.,abs(a)) for a,b in zip(coefficients,adopted))
                save(output/name/'common_update.json',dict(coefficients=adopted,proposed=proposed,dampingFraction=fraction,
                    oldPenalty=old_penalty,newPenalty=penalty(next_result['slices'],adopted,fixed),
                    coefficientRelativeDrift=drift,audit=fit_audit,
                    updateState='accepted' if fraction>0 else 'stalled_or_unchanged',
                    interpretation='QR block update only; not a claim of local/global likelihood convergence'))
                coefficients,previous=adopted,next_result;number+=1
                if drift<cfg['smoothing']['coefficientRelativeTolerance']:
                    break
        final=execute('%02d_common_raw_refit'%number,'common',seed_update(previous,rows),coefficients)
        coefficients=final['sharedCoefficients']
        if len(coefficients)!=14:
            raise ValueError('Final common-function coefficient count differs')
        (output/'analytic_candidate.yaml').write_text(yaml.safe_dump(candidate_yaml(coefficients,cfg,'common'),sort_keys=False))
        save(output/'final_summary.json',dict(phase=history[-1]['phase'],summary=result_summary(final),
            freeSummary=json.loads((output/'00_local_free/summary.json').read_text()),
            physicsValidationComplete=False,productionEnabled=False))
        if not args.no_render:
            regularized=[h for h in history if h['mode']=='regularized']
            selected=[(history[0],'Local_free_projected'),(regularized[0],'Early_regularized_projected'),
                      (regularized[-1],'Last_regularized_projected'),(history[-1],'Common_raw_intrinsic')]
            unique=[]
            for item,label in selected:
                if item['phase'] not in [p['phase'] for p,label in unique]:
                    unique.append((item,label))
            summary_rows=[]
            for ordinal,(item,label) in enumerate(unique):
                data=json.loads((output/item['phase']/'fit_result.json').read_text())
                for row in data['slices']:
                    summary_rows.append([ordinal,label,row['index'],row['r'],
                        row['he3']['mean'],row['he3']['sigma'],row['he4']['mean'],row['he4']['sigma'],
                        row['dataDeviance'],row['chi2'],row['ndf'],int(row['valid'])])
            tsv(output/'summary_curves.tsv',summary_rows)
            tsv(output/'summary_settings.tsv',cfg['display'].items())
            summary_macro=output/'source/tools/he4_purity/iterative_logpoly/DrawIterativeSummary.C'
            call=str(summary_macro)+'+('+','.join(json.dumps(str(p)) for p in
                (output/'summary_curves.tsv',output/'summary_settings.tsv',output/'iterative_calibration_summary.pdf'))+')'
            run_command(['bash',str(wrapper),str(mainconf),'root4star','-l','-b','-q',call],
                        output/'logs/summary.log',cfg['executionTimeoutSeconds'],'ITERATIVE_SUMMARY_COMPLETE')
    preserve()
    products={str(p.relative_to(output)):digest(p) for p in output.rglob('*') if p.is_file() and
              'source' not in p.relative_to(output).parts and p.name not in ('context.json','completed.json','progress.json')}
    save(output/'completed.json',dict(executionComplete=True,physicsValidationComplete=False,
        productionEnabled=False,products=products,completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))
    save(output/'progress.json',dict(state='complete',history=history))
    print('ITERATIVE_RUNNER_COMPLETE',output,flush=True)


if __name__=='__main__':
    main()
