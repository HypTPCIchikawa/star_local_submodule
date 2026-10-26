#!/usr/bin/env python3
import collections,math,pathlib,shutil,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
sys.path.insert(0,str(HERE.parent));import summarize as shared
R,C,O=run.R,run.C,run.OUT;old=shared.load_old();native=[];bins=[];examples=[];groups=collections.defaultdict(list);raw={}
for row in run.tab(O/'candidate_gate_counts/gate_counts.tsv'):raw[int(row['rBin']),row['species']]={k:float(row[k]) for k in ['rawGateCount','rawGateCountLowerFullBins','rawGateCountUpperTouchingBins']}
wide={}
for name in ['atlas_low_002/candidate_atlas.tsv','atlas_high_002/candidate_atlas.tsv']:
 for a in run.tab(O/name):
  if a['model']=='gaussian_linear_nominal':wide[a['band'],a['species'],round(float(a['rLo']),6)]=a
for a in run.tab(O/'he4_r061_001/atlas.tsv'):
 if a['band']=='gaussian' and a['model'].startswith('linear_'):wide['low','he4',round(float(a['rLo']),6)]=a
for job in run.jobs():
 folder=R/job['name'];assert run.read(folder/'completed.json')['executionComplete'];sp=job['species'];band=job['band'];cfg=run.read(folder/'effective_settings.json');cfg['fitWindowMode']='target_local';audit={r['index']:r for r in run.read(folder/'tag_shape_audit.json')};linearref={e['index']:e for e in run.read(folder/'linear_seed_provenance.json')};cache={}
 for row in run.read(folder/'cutoff_fits.json')['slices']:
  v=old.aggregate_one([row],sp,audit,cfg,raw if band=='low' else None)
  if band=='high':
   drop=row['components'][sp].get('conditionalDropDeltaD')
   if drop is None or drop<cfg['support']['minDropDeltaD']:v['supported']=False;v['reasons']=sorted(set(v['reasons']+['conditional_drop_deviance']))
  c=row['components'][sp];ref=linearref[row['index']];src=pathlib.Path(ref['reference'])
  if str(src) not in cache:cache[str(src)]={r['index']:r for r in run.read(src)['slices']}
  before=cache[str(src)][row['index']];sel=c['selection'];native.append(dict(job=job['name'],band=band,species=sp,index=row['index'],rLo=row['lo'],rHi=row['hi'],mean=c['mean'],sigma=c['sigma'],purity=sel['purity'],conditionalError=sel['conditionalError'],numericalValid=row['numericalValid'],supported=v['supported'],reasons=v['reasons'],chi2=row['pearson'],ndf=row['conditionalNdf'],deviance=row['dataDeviance'],linearChi2=before['pearson'],linearNdf=before['conditionalNdf'],linearPurity=before['components'][sp]['selection']['purity'],linearMean=ref['targetMean'],linearSigma=ref['targetSigma'],meanFixed=ref['meanFixed'],linearSeedAvailable=ref['linearSeedAvailable'],root=str(folder/'cutoff_fits.root')))
  a=wide[band,sp,round(row['lo'],6)];bounds={k:float(a[k]) for k in ['linearLo','linearHi','logLo','logHi']}
  for kind in ['linear','log']:
   bounds[kind+'Lo']=min(bounds[kind+'Lo'],row['fitLo']);bounds[kind+'Hi']=max(bounds[kind+'Hi'],row['fitHi'])
  examples.append(dict(group=band+'_'+sp,species=sp,index=row['index'],rLo=row['lo'],rHi=row['hi'],root=str(folder/'cutoff_fits.root'),fitLo=row['fitLo'],fitHi=row['fitHi'],**bounds,purity=sel['purity'],chi2=row['pearson'],ndf=row['conditionalNdf'],mean=c['mean'],sigma=c['sigma'],he3Mean=row['components']['he3']['mean'],he3Sigma=row['components']['he3']['sigma'],sourceModel='gaussian_none',supported=v['supported'],numericalValid=row['numericalValid'],meanFixed=ref['meanFixed']))
  group=round((row['lo']-(C['aggregateOrigin'] if band=='low' else C['tofThreshold']))/C['aggregateWidth'],7)
  group=math.floor(group+1e-7);groups[band,sp,group].append((row,audit,cfg,v))
for (band,sp,g),items in sorted(groups.items()):
 rows=[a[0] for a in items];audit={k:v for a in items for k,v in a[1].items()};cfg=items[0][2];v=old.aggregate_one(rows,sp,audit,cfg,raw if band=='low' else None)
 if band=='high' and not all(a[3]['supported'] for a in items):v['supported']=False;v['reasons']=sorted(set(v['reasons']+['one_or_more_native_fits_unassessed']))
 expected=4 if sp=='he4' and band=='low' and g==4 else 1 if sp=='he4' and band=='high' else 5
 assert len(rows)==expected,(band,sp,g,len(rows),expected)
 bins.append(dict(band=band,species=sp,rLo=min(r['lo'] for r in rows),rHi=max(r['hi'] for r in rows),purity=v['purity'],statConditionalBound=v['conditionalErrorWeightedBound'],quality=v['supported'],signalSelected=v['signalSelected'],totalSelected=v['totalSelected'],rawSelected=v['rawSelected'],localResidualFrac=v['localResidualFrac'],reasons=v['reasons'],response='gaussian',background='none',variation='nominal',sourceModel='gaussian_none',assessment=v))
assert len(bins)==364
axes=[];stats=[]
for sp in C['species']:
 good=[r for r in bins if r['species']==sp and r['quality']]
 extent=good or [r for r in bins if r['species']==sp];lo,hi=C['fixedDAxis'] if sp=='d' else (min(r['rLo'] for r in extent)-C['axisPadding'],max(r['rHi'] for r in extent)+C['axisPadding']);axes.append(dict(species=sp,lo=lo,hi=hi))
 for band in ['low','high']:
  rows=[r for r in bins if r['species']==sp and r['band']==band];v=[r for r in rows if r['quality']];stats.append(dict(species=sp,band=band,supported=len(v),total=len(rows),rMin=min((r['rLo'] for r in v),default=None),rMax=max((r['rHi'] for r in v),default=None),minPurity=min((r['purity'] for r in v),default=None),below50=[r['rLo'] for r in v if r['purity']<C['purityThreshold']]))
run.save(R/'summary.json',dict(configuration=C,rows=bins,statistics=stats,axes=axes));run.save(R/'native.json',native);run.save(R/'examples.json',examples)
shared.tsv(R/'purity.tsv',bins,'band species rLo rHi purity statConditionalBound quality signalSelected totalSelected rawSelected localResidualFrac response background variation sourceModel'.split());shared.tsv(R/'axes.tsv',axes,['species','lo','hi']);shared.tsv(R/'native.tsv',native,list(native[0]));shared.tsv(R/'examples.tsv',examples,list(examples[0]))
for n in ['coefficients.json','curves.tsv']:shutil.copy2(O/n,R/n)
selection=run.read(O/'he4_r061_001/pid_selection.yaml');selection['purityEstimation']=C;run.save(R/'pid_selection.yaml',selection)
print('GAUSSIAN_ONLY_SUMMARY',len(native),'native fits;',len(bins),'aggregate bins;',stats,flush=True)
