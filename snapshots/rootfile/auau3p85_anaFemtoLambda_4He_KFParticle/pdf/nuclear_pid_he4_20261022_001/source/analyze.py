#!/usr/bin/env python3
import collections,math,pathlib,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
sys.path.insert(0,str(run.P/'tools/nuclear_pid_smooth_20261021'));import summarize as shared
old=shared.load_old();R,C=run.R,run.C;native=[];examples=[];bins=[];before=run.read(run.OLD/'he4_r061_001/atlas.tsv') if False else None
wide={}
for a in run.tab(run.OLD/'he4_r061_001/atlas.tsv'):
 if a['band']=='gaussian' and a['model'].startswith('linear_'):wide['ALL',round(float(a['rLo']),6)]=a
for a in run.tab(run.OLD/'atlas_high_002/candidate_atlas.tsv'):
 if a['species']=='he4' and a['model']=='gaussian_linear_nominal':wide['TAG_D',round(float(a['rLo']),6)]=a
for region in C['regions']:
 folder=R/('candidate_'+region['name']);cfg=run.read(folder/'effective_settings.json');audit={r['index']:r for r in run.read(folder/'tag_shape_audit.json')};group=collections.defaultdict(list)
 for r in run.read(folder/'cutoff_fits.json')['slices']:
  c=r['components']['he4'];v=old.aggregate_one([r],'he4',audit,cfg)
  if region['sample']!='ALL':
   drop=c.get('conditionalDropDeltaD')
   if drop is None or drop<cfg['support']['minDropDeltaD']:v['supported']=False;v['reasons']=sorted(set(v['reasons']+['conditional_drop_deviance']))
  s=c['selection'];native.append(dict(region=region['name'],index=r['index'],rLo=r['lo'],rHi=r['hi'],mean=c['mean'],sigma=c['sigma'],purity=s['purity'],conditionalError=s['conditionalError'],chi2=r['pearson'],ndf=r['conditionalNdf'],deviance=r['dataDeviance'],numericalValid=r['numericalValid'],supported=v['supported'],reasons=v['reasons'],rawSelected=s['rawGateCount'],signalSelected=s['signalSelected'],totalSelected=s['totalSelected'],localResidualFrac=s['localResidualFrac']))
  a=wide[region['sample'],round(r['lo'],6)];bounds={k:float(a[k]) for k in ['linearLo','linearHi','logLo','logHi']}
  for kind in ['linear','log']:bounds[kind+'Lo']=min(bounds[kind+'Lo'],r['fitLo']);bounds[kind+'Hi']=max(bounds[kind+'Hi'],r['fitHi'])
  examples.append(dict(group=region['name'],region=region['name'],species='he4',index=r['index'],rLo=r['lo'],rHi=r['hi'],root=str(folder/'cutoff_fits.root'),fitLo=r['fitLo'],fitHi=r['fitHi'],**bounds,purity=s['purity'],chi2=r['pearson'],ndf=r['conditionalNdf'],mean=c['mean'],sigma=c['sigma'],he3Mean=r['components']['he3']['mean'],he3Sigma=r['components']['he3']['sigma'],sourceModel='gaussian_'+region['background'],background=region['background'],numericalValid=r['numericalValid'],supported=v['supported'],meanFixed=False))
  g=math.floor((r['lo']-C['aggregateOrigin'])/C['aggregateWidth']+1e-7);group[g].append((r,v))
 for g,items in sorted(group.items()):
  rows=[r for r,v in items];v=old.aggregate_one(rows,'he4',audit,cfg)
  if region['sample']!='ALL' and not all(v['supported'] for r,v in items):v['supported']=False;v['reasons']=sorted(set(v['reasons']+['one_or_more_native_fits_unassessed']))
  lo=min(r['lo'] for r in rows);hi=max(r['hi'] for r in rows);expected=4 if lo==C['minimumRExclusive'] else 1 if region['sample']!='ALL' else 5;assert len(rows)==expected
  bins.append(dict(region=region['name'],band='low' if region['sample']=='ALL' else 'high',species='he4',rLo=lo,rHi=hi,purity=v['purity'],statConditionalBound=v['conditionalErrorWeightedBound'],quality=v['supported'],signalSelected=v['signalSelected'],totalSelected=v['totalSelected'],rawSelected=v['rawSelected'],localResidualFrac=v['localResidualFrac'],reasons=v['reasons'],assessment=v))
assert len(native)==199 and len(bins)==88
statistics=[]
for region in C['regions']:
 rows=[r for r in bins if r['region']==region['name']];valid=[r for r in rows if r['quality']];n=[r for r in native if r['region']==region['name']];statistics.append(dict(region=region['name'],total=len(rows),supported=len(valid),native=len(n),numericalValid=sum(r['numericalValid'] for r in n),minPurity=min((r['purity'] for r in valid),default=None),below50=[r for r in valid if r['purity']<C['purityThreshold']]))
run.save(R/'summary.json',dict(configuration=C,rows=bins,statistics=statistics));run.save(R/'native.json',native);run.save(R/'examples.json',examples)
for filename,rows in [('purity.tsv',bins),('native.tsv',native),('examples.tsv',examples)]:shared.tsv(R/filename,rows,[k for k in rows[0] if k!='assessment'])
run.save(R/'pid_selection.yaml',dict(productionEnabled=False,referenceGit=C['referenceGit'],coefficients=run.read(R/'coefficients.json')['candidate'],strictMinimumRigidityGeVcm={'he4':C['minimumRExclusive']},tofMass2={'he4':C['tofMass2']},purityEstimation=C))
print('SUMMARY',statistics,flush=True)
