"""Saved-fit comparison and conservative purity summary, without refitting."""
from collections import Counter
import math
from pathlib import Path
import yaml
from prepare_tag_study import SPECIES,read,save,sha,tsv,require

def comparison(fixed_dir,oracle_dir,all_dir,output):
    plan=read(fixed_dir/'plan.json');cal=Path(plan['coefficientDirectory']);local_path=Path(read(cal.parent/'completed.json')['sourceFitPath'])
    local={q['id']:q for q in read(local_path)['results']};oracle={q['id']:q for q in read(oracle_dir/'fits.json')['results']};all_rows={q['id']:q for q in read(all_dir/'fits.json')['results']};fixed=read(fixed_dir/'fits.json')['results'];rows=[]
    protected={str(p):sha(p) for p in (local_path,oracle_dir/'fits.json',all_dir/'fits.json',fixed_dir/'fits.json')}
    for f in fixed:
        a,o,l=all_rows[f['id']],oracle[f['id']],local[f['id']];s=f['target'];li=next(q for q in l['components'] if q['species']==s);fi=next(q for q in f['components'] if q['species']==s)
        delta=o['dataDeviance']-l['dataDeviance']
        require(abs(delta)<1.e-5*max(1.,abs(l['dataDeviance'])),'Saved local response not reproduced for '+f['id'])
        rows.append(dict(id=f['id'],species=s,r=(f['rLo']+f['rHi'])/2.,localChi2PerNdf=l['pearsonPerNdf'],oracleFixedShapeChi2PerNdf=o['pearsonPerNdf'],commonPolynomialTofChi2PerNdf=f['pearsonPerNdf'],commonPolynomialAllChi2PerNdf=a['pearsonPerNdf'],localDeviance=l['dataDeviance'],oracleDeltaDeviance=delta,polynomialExtraDeviance=f['dataDeviance']-l['dataDeviance'],meanShiftInLocalSigma=(fi['mean']-li['mean'])/li['sigma'],sigmaRatio=fi['sigma']/li['sigma'],localWidthRatio=l['responseWidthRatio'],polynomialWidthRatio=sum(q['weight']*q['widthRatio'] for q in f['fixedResponseNodes'][SPECIES.index(s)]['nodes'])))
    save(output/'transfer_comparison.json',dict(status='completed',productionEnabled=False,oracleReproductionPassed=True,rows=rows,inputs=protected,ndfMeaning='Oracle fixes the three shape parameters: NDF differs from local fit. Compare deviance for engine identity.'))
    header=tuple(rows[0]);tsv(output/'transfer_comparison.tsv',header,([q[k] for k in header] for q in rows))
    text=['# 共通多項式への圧縮・ALL転用の切り分け','', '保存済みTOF局所応答を同じデータ・ROIで新エンジンに渡し、devianceの再現を確認した。', '局所fitから3つのshape parameterを固定すると自由度は変わるため、同値性はχ²/ndfではなくdevianceで判定。', '', '| Species | R | Local TOF χ²/ndf | Fixed polynomial TOF χ²/ndf | Fixed polynomial ALL χ²/ndf | Δμ/σlocal | σpoly/σlocal |', '|---|---:|---:|---:|---:|---:|---:|']
    for q in rows:text.append('| {species} | {r:.3f} | {localChi2PerNdf:.2f} | {commonPolynomialTofChi2PerNdf:.2f} | {commonPolynomialAllChi2PerNdf:.2f} | {meanShiftInLocalSigma:.3f} | {sigmaRatio:.3f} |'.format(**q))
    text+=['','多項式固定による不一致はTOF標本だけでも発生する。ALLへの転用だけを原因と断定しない。','この診断の不良fitは新しいcutoffや系統誤差へ採用しない。']
    with (output/'transfer_comparison.md').open('x') as stream:stream.write('\n'.join(text)+'\n')
    return protected

def good_runs(rows,minimum=3):
    runs=[];current=[]
    for q in rows:
        if q['assessed'] and q['lowerModelPurity']>=.5:
            if current and q['r']-current[-1]['r']>.061:
                if len(current)>=minimum:runs.append(current)
                current=[]
            current.append(q)
        else:
            if len(current)>=minimum:runs.append(current)
            current=[]
    if len(current)>=minimum:runs.append(current)
    return runs

def summary(directories,output):
    require(len(directories)>=2,'At least two explicit model variants required')
    plans=[read(p/'plan.json') for p in directories];fits=[read(p/'fits.json') for p in directories]
    for p in directories:require(read(p/'completed.json')['status']=='completed','Unfinished saved study')
    require(all(p['coefficientDirectory']==plans[0]['coefficientDirectory'] and p['cache']==plans[0]['cache'] for p in plans),'Model variants must share curves and raw input')
    require(all([q['id'] for q in f['results']]==[q['id'] for q in fits[0]['results']] for f in fits),'Compared R grids differ')
    protected={str(d/n):sha(d/n) for d in directories for n in ('fits.json','plan.json','completed.json')}
    points=[];fit_quality={}
    for plan,f in zip(plans,fits):
        fit_quality[plan['variant']]=dict(total=f['totalSlices'],executionFailures=f['fitFailures'],accepted=f['qualityAcceptedSlices'],reasons=dict(Counter(r for q in f['results'] for r in q.get('purityQualityReasons',[]))))
    for index,b in enumerate(fits[0]['results']):
        valid=[(p['variant'],f['results'][index]) for p,f in zip(plans,fits) if f['results'][index].get('purityAvailable') is True]
        assessed=b.get('purityAvailable') is True and len(valid)>=2
        envelope=min(q['purity'] for _,q in valid) if assessed else None
        reasons=b.get('purityQualityReasons',[b.get('error','unassessed')]);gate=b.get('winnerGate',{})
        if b.get('purityAvailable') and len(valid)<2:reasons=['fewer_than_two_valid_models']
        points.append(dict(id=b['id'],species=b['target'],r=(b['rLo']+b['rHi'])/2.,rLo=b['rLo'],rHi=b['rHi'],mode=b['mode'],assessed=assessed,purity=b.get('purity') if assessed else None,conditionalError=b.get('purityConditionalError') if assessed else None,lowerModelPurity=envelope,validModels=[name for name,q in valid],reason=','.join(reasons),signal=gate.get('signal'),background=gate.get('background'),targetEffectiveCount=(b.get('effectiveCount') or [0.]*4)[SPECIES.index(b['target'])],chi2PerNdf=b.get('pearsonPerNdf'),modelValues=[dict(variant=name,purity=q['purity'],conditionalError=q['purityConditionalError']) for name,q in valid]))
    # Reassess the count threshold, not merely recount the >=100 accepted points.
    sensitivity={str(n):{s:[] for s in SPECIES} for n in (20,50,100)}
    for index,b in enumerate(fits[0]['results']):
        target=SPECIES.index(b['target'])
        for threshold in (20,50,100):
            candidates=[]
            for f in fits:
                q=f['results'][index];reasons=set(q.get('purityQualityReasons',['execution_failed']))-{'insufficient_target_effective_count'}
                neff=(q.get('effectiveCount') or [0.]*4)[target]
                gate=q.get('winnerGate',{});p=gate.get('candidatePurity')
                ok=not reasons and neff>=threshold and isinstance(p,(float,int)) and 0<=p<=1
                candidates.append((ok,p))
            if candidates[0][0] and sum(q[0] for q in candidates)>=2:
                sensitivity[str(threshold)][b['target']].append(dict(r=(b['rLo']+b['rHi'])/2.,mode=b['mode'],lowerModelPurity=min(q[1] for q in candidates if q[0])))
    cutoffs={};minimum=plans[0]['config']['cutoffEstimation']['minimumConsecutiveSampledPoints']
    for s in SPECIES:
        selected=sorted((q for q in points if q['species']==s),key=lambda q:q['r']);low=[q for q in selected if q['mode']=='ALL'];high=[q for q in selected if q['mode']!='ALL']
        lr,hr=good_runs(low,minimum),good_runs(high,minimum)
        below=[q for q in low if q['assessed'] and q['lowerModelPurity']<.5]
        isolated=[]
        for i,q in enumerate(selected[1:-1],1):
            left,right=selected[i-1],selected[i+1]
            if q['assessed'] and q['lowerModelPurity']<.5 and all(x['assessed'] and x['lowerModelPurity']>=.5 for x in (left,right)):
                isolated.append(dict(id=q['id'],r=q['r'],action='review_neighbor_fits_no_automatic_veto'))
        cutoffs[s]=dict(productionEnabled=False,adoptionAllowed=False,assessedPoints=sum(q['assessed'] for q in selected),assessedLowRPoints=sum(q['assessed'] for q in low),assessedHighRPoints=sum(q['assessed'] for q in high),candidateLowCoverageStart=lr[0][0]['rLo'] if lr else None,candidateHighCoverageEnd=hr[-1][-1]['rHi'] if hr else None,lowPurityCrossingObserved=bool(below),candidateLowRCutoff=None,finalHighRLimit=None,qualityAcceptedSampledRuns=[dict(mode=run[0]['mode'],rMin=run[0]['rLo'],rMax=run[-1]['rHi'],sampledPoints=len(run)) for run in lr+hr],isolatedLowPurityPoints=isolated,reason='No production cutoff is inferred across unassessed fit intervals or unsupported curve domains.',countThresholdSensitivity={str(n):dict(assessedPoints=len(sensitivity[str(n)][s]),highestAssessedR=max((q['r'] for q in sensitivity[str(n)][s]),default=None),lowRPoints=sum(q['mode']=='ALL' for q in sensitivity[str(n)][s])) for n in (20,50,100)})
    save(output/'purity_summary.json',dict(status='completed',productionEnabled=False,cutoffAdoptionAllowed=False,points=points,fitQuality=fit_quality,cutoffCandidates=cutoffs,unsupportedSlices=plans[0]['skippedSlices'],inputs=protected,interpretation=['P=S/(S+B), not S/B.','Bad fits excluded from the model envelope.','Baseline plus at least one other valid model required.','Track-weighted, charge-inclusive, centrality-integrated; not CF-pair-weighted.','Fixed-curve conditional errors do not include PID model uncertainty.','No automatic isolated vetoes; no extrapolation of TOF-supported curves.']))
    header=('species','r','mode','purity','conditionalError','lowerModelPurity','validModels','assessed','reason','signal','targetEffectiveCount','chi2PerNdf')
    def number(q,k):
        v=q.get(k);return v if isinstance(v,(float,int)) and math.isfinite(v) else -1.
    tsv(output/'purity_points.tsv',header,((q['species'],q['r'],q['mode'],number(q,'purity'),number(q,'conditionalError'),number(q,'lowerModelPurity'),len(q['validModels']),int(q['assessed']),q['reason'] or 'assessed',number(q,'signal'),number(q,'targetEffectiveCount'),number(q,'chi2PerNdf')) for q in points))
    with (output/'cutoff_candidates.yaml').open('x') as stream:yaml.safe_dump(dict(schemaVersion=1,productionEnabled=False,adoptionAllowed=False,purityDefinition='S/(S+B)',minimumPurity=.5,species=cutoffs),stream,sort_keys=False,allow_unicode=True)
    lines=['# 共通多項式を用いたALL/TOF純度の評価','', '| Species | Assessed points | R < 2 | R >= 2 | Final cutoff |','|---|---:|---:|---:|---|']
    for s,q in cutoffs.items():lines.append('| %s | %d | %d | %d | Not adopted |'%(s,q['assessedPoints'],q['assessedLowRPoints'],q['assessedHighRPoints']))
    lines+=['','## 採否','', '候補曲線を固定すると、不良fitが多く残る。χ²不良点を純度0や50%未満と置き換えず、', 'モデル差の帯にも含めない。空白を跨ぐ単一の下限・上限や、孤立vetoは本番へ採用しない。','', '図の灰色×は未評価を示す記号で、純度の測定値ではない。誤差は固定曲線の条件付き誤差。', '低RはALL、高Rは同じ核種別TOF窓。4Heはdと完全に同じm²タグを使用。','', '統計のある高R点でも、曲線共通支持範囲の端をデータの統計限界と呼ばない。', '単一多項式への圧縮誤差とTOF→ALLの転用差の切り分けを別途参照。']
    with (output/'purity_summary.md').open('x') as stream:stream.write('\n'.join(lines)+'\n')
    print('PURITY_SUMMARY_COUNTS '+str({s:(q['assessedLowRPoints'],q['assessedHighRPoints']) for s,q in cutoffs.items()}),flush=True)
    return protected
