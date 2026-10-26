#!/usr/bin/env python3
import collections,csv,math,pathlib,statistics,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
sys.path.insert(0,str(HERE.parent));import summarize as shared
P,OUT,ROOT=run.P,run.OUT,run.ROOT;read,save,sha=run.read,run.save,run.sha
old=shared.load_old();metrics=[];detail=[];stats=[];comparison=[];overview=[]
wide=list(csv.DictReader((OUT/'atlas_low_002/candidate_atlas.tsv').open(),delimiter='\t'));allmodels={};cfgs={};audits={}
def finite(x):return isinstance(x,(int,float)) and math.isfinite(x)
def valid(row):return row['numericalValid'] and finite(row['dataDeviance']) and finite(row['pearson']) and row['conditionalNdf']>0
for response in run.cfg['responses']:
 folders={'none':ROOT/('none_'+response+'_none')}
 for bg in run.cfg['backgroundsForControls']:folders[bg]=ROOT/('baseline_'+response+'_'+bg) if 'fitWindowPolicy' in run.cfg else OUT/('candidate/low_he4_'+response+'_'+bg+'_nominal');folders[bg+'_warm']=ROOT/('warm_'+response+'_'+bg)
 data={key:read(p/'cutoff_fits.json')['slices'] for key,p in folders.items()};allmodels[response]=data
 settings=read(folders['linear']/'effective_settings.json');settings['fitWindowMode']='target_local';cfgs[response]=settings;audit={a['index']:a for a in read(folders['linear']/'tag_shape_audit.json')};audits[response]=audit
 group=collections.defaultdict(dict)
 for i,nobg in enumerate(data['none']):
  variants={'none':(nobg,folders['none'],'no_smooth')}
  for bg in run.cfg['backgroundsForControls']:
   ori,warm=data[bg][i],data[bg+'_warm'][i]
   options=[(r,p,label) for r,p,label in [(ori,folders[bg],'original_start'),(warm,folders[bg+'_warm'],'no_smooth_start')] if valid(r)]
   chosen=min(options,key=lambda t:t[0]['dataDeviance']) if options else (ori,folders[bg],'failed')
   variants[bg]=chosen
  view=[r for r in wide if r['species']=='he4' and r['model'].startswith(response+'_') and int(r['index'])==i]
  bounds={k:min(float(r[k]) for r in view) if k.endswith('Lo') else max(float(r[k]) for r in view) for k in ['linearLo','linearHi','logLo','logHi']}
  bounds['linearLo']=min(bounds['linearLo'],nobg['fitLo']);bounds['logLo']=min(bounds['logLo'],nobg['fitLo'])
  for model,(row,folder,start) in variants.items():
   assess=old.aggregate_one([row],'he4',audit,settings);sel=row['components']['he4']['selection'];mean=row['components']['he4']['mean'];sigma=row['components']['he4']['sigma'];chi=row['pearson']/row['conditionalNdf'] if valid(row) else None
   record=dict(response=response,model=model,index=i,r=row['r'],rLo=row['lo'],rHi=row['hi'],chi2=row['pearson'],ndf=row['conditionalNdf'],chi2Ndf=chi,deviance=row['dataDeviance'],numericalValid=valid(row),identityPass=row['targetIdentityScreenPass'],supported=assess['supported'],purity=sel['purity'],gateResidual=sel['localResidualFrac'],mean=mean,sigma=sigma,start=start,folder=str(folder),reasons=assess['reasons'])
   metrics.append(record);group[i//5].setdefault(model,[]).append(row)
   comparison.append(dict(band=response,species='he4',model=model+'_'+start,index=i,rLo=row['lo'],rHi=row['hi'],fitLo=row['fitLo'],fitHi=row['fitHi'],mean=mean,sigma=sigma,root=str(folder/'cutoff_fits.root'),purity=sel['purity'],numericalValid=row['numericalValid'],supported=assess['supported'],**bounds,chi2=row['pearson'],ndf=row['conditionalNdf'],deviance=row['dataDeviance'],identityPass=row['targetIdentityScreenPass']))
  for bg in run.cfg['backgroundsForControls']:
   full=variants[bg][0];original=data[bg][i];warm=data[bg+'_warm'][i];pair=dict(response=response,bg=bg,index=i,r=nobg['r'],noneValid=valid(nobg),bgValid=valid(full),noneIdentity=nobg['targetIdentityScreenPass'],bgIdentity=full['targetIdentityScreenPass'],nestedViolation=valid(nobg) and valid(full) and full['dataDeviance']>nobg['dataDeviance']+run.cfg['nestedDevianceTolerance'],originalDeviance=original['dataDeviance'],warmDeviance=warm['dataDeviance'],bestStart=variants[bg][2])
   if valid(nobg) and valid(full):pair.update(chi2Ratio=nobg['pearson']/full['pearson'] if full['pearson']>0 else None,chi2NdfRatio=(nobg['pearson']/nobg['conditionalNdf'])/(full['pearson']/full['conditionalNdf']) if full['pearson']>0 else None,deltaDeviance=nobg['dataDeviance']-full['dataDeviance'])
   detail.append(pair)
 for g,models in sorted(group.items()):
  for model,rows in models.items():
   v=old.aggregate_one(rows,'he4',audit,settings);overview.append(dict(response=response,model=model,rLo=rows[0]['lo'],rHi=rows[-1]['hi'],**v))
 for bg in run.cfg['backgroundsForControls']:
  paired=[r for r in detail if r['response']==response and r['bg']==bg and r['noneValid'] and r['bgValid']];physical=[r for r in paired if r['noneIdentity'] and r['bgIdentity']]
  stats.append(dict(response=response,bg=bg,pairedNumerical=len(paired),nonePearsonBetter=sum(r['chi2Ratio']<1 for r in paired),noneReducedChiBetter=sum(r['chi2NdfRatio']<1 for r in paired),medianChi2Ratio=statistics.median(r['chi2Ratio'] for r in paired),medianReducedChi2Ratio=statistics.median(r['chi2NdfRatio'] for r in paired),pairedIdentity=len(physical),identityNonePearsonBetter=sum(r['chi2Ratio']<1 for r in physical),identityNoneReducedChiBetter=sum(r['chi2NdfRatio']<1 for r in physical),identityMedianChi2Ratio=statistics.median(r['chi2Ratio'] for r in physical),nestedViolations=sum(r['nestedViolation'] for r in paired)))
save(ROOT/'comparison.json',dict(productionEnabled=False,scope='He4 ALL R<2; background none means four nucleus columns with zero smooth columns',metrics=metrics,pairs=detail,statistics=stats,aggregated=overview))
shared.tsv(ROOT/'metrics.tsv',metrics,'response model index r rLo rHi chi2 ndf chi2Ndf deviance numericalValid identityPass supported purity gateResidual mean sigma start folder'.split())
shared.tsv(ROOT/'atlas.tsv',comparison,'band species model index rLo rHi fitLo fitHi mean sigma root purity numericalValid supported linearLo linearHi logLo logHi chi2 ndf deviance identityPass'.split())
shared.tsv(ROOT/'aggregate.tsv',overview,'response model rLo rHi purity signalSelected totalSelected rawSelected localResidualFrac supported'.split())
print('NO_SMOOTH_COMPARISON',stats,flush=True)
