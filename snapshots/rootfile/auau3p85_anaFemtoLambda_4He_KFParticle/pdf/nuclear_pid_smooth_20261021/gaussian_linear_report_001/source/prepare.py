#!/usr/bin/env python3
import collections,csv,json,math,pathlib,shutil,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent));from study import P,OUT,read,save,sha
from summarize import tsv
C=read(HERE/'settings.yaml');R=OUT/C['outputDirectory'];LATEST=OUT/'he4_r061_001'
def tab(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def main():
 if R.exists():raise ValueError('Never overwrite '+str(R))
 R.mkdir();(R/'source').mkdir();protected={}
 for base,name in [(OUT,'candidate_summary.json'),(OUT,'coefficients.json'),(OUT,'curves.tsv'),(LATEST,'summary.json'),(LATEST,'completed.json')]:protected[str(base/name)]=sha(base/name)
 prior=read(OUT/'candidate_summary.json');latest=read(LATEST/'summary.json');rows=[]
 for r in latest['bins']:
  key='gaussian_linear' if r['species']=='he4' and r['band']=='low' else 'gaussian_linear_nominal';v=r['models'][key]
  rows.append(dict(band=r['band'],species=r['species'],rLo=r['rLo'],rHi=r['rHi'],purity=v['purity'],statConditionalBound=v['conditionalErrorWeightedBound'],quality=v['supported'],signalSelected=v['signalSelected'],totalSelected=v['totalSelected'],rawSelected=v['rawSelected'],localResidualFrac=v['localResidualFrac'],reasons=v['reasons'],response=C['response'],background=C['background'],variation=C['variation'],source=str(LATEST/'summary.json') if r['species']=='he4' and r['band']=='low' else str(OUT/'candidate_summary.json'),sourceModel=key))
 assert len(rows)==len(latest['bins'])
 for r in rows:
  assert r['response']=='gaussian' and r['background']=='linear'
  if r['totalSelected']>0:assert abs(r['purity']-r['signalSelected']/r['totalSelected'])<1e-12
 axes=[];stats=[]
 for s in ['d','t','he3','he4']:
  valid=[r for r in rows if r['species']==s and r['quality']];assert valid
  lo,hi=C['fixedDAxis'] if s=='d' else (min(r['rLo'] for r in valid)-C['axisPadding'],max(r['rHi'] for r in valid)+C['axisPadding']);axes.append(dict(species=s,lo=lo,hi=hi))
  for band in ['low','high']:
   v=[r for r in rows if r['species']==s and r['band']==band and r['quality']];below=[r for r in v if r['purity']<C['purityThreshold']];stats.append(dict(species=s,band=band,supported=len(v),total=sum(r['species']==s and r['band']==band for r in rows),rMin=min(r['rLo'] for r in v) if v else None,rMax=max(r['rHi'] for r in v) if v else None,minPurity=min(r['purity'] for r in v) if v else None,below50=[dict(rLo=r['rLo'],rHi=r['rHi'],purity=r['purity']) for r in below]))
 tsv(R/'purity.tsv',rows,'band species rLo rHi purity statConditionalBound quality signalSelected totalSelected rawSelected localResidualFrac response background variation sourceModel'.split());tsv(R/'axes.tsv',axes,['species','lo','hi']);save(R/'summary.json',dict(configuration=C,rows=rows,statistics=stats,axes=axes))
 low=tab(OUT/'atlas_low_002/candidate_atlas.tsv');high=tab(OUT/'atlas_high_002/candidate_atlas.tsv');new=tab(LATEST/'atlas.tsv');examples=[]
 for group,targets in [('low',C['fitExamplesLow']),('anchors',C['fitExamplesAnchors']),('tof',C['fitExamplesHigh'])]:
  for target in targets:
   species=target['species'];pool=[r for r in (new if species=='he4' and group!='tof' else high if group=='tof' else low) if r['species']==species and (r.get('band')=='gaussian' and r['model'].startswith('linear_') if species=='he4' and group!='tof' else r['model']=='gaussian_linear_nominal')]
   a=min(pool,key=lambda r:abs(.5*(float(r['rLo'])+float(r['rHi']))-target['r']));folder=pathlib.Path(a['root']).parent;protected[str(folder/'cutoff_fits.json')]=sha(folder/'cutoff_fits.json');protected[a['root']]=sha(pathlib.Path(a['root']));data=read(folder/'cutoff_fits.json');r=next(r for r in data['slices'] if r['index']==int(a['index']));v=r['components'][species]['selection']
   examples.append(dict(group=group,species=species,index=r['index'],rLo=r['lo'],rHi=r['hi'],root=a['root'],fitLo=r['fitLo'],fitHi=r['fitHi'],linearLo=min(float(a['linearLo']),r['fitLo']),linearHi=max(float(a['linearHi']),r['fitHi']),logLo=min(float(a['logLo']),r['fitLo']),logHi=max(float(a['logHi']),r['fitHi']),purity=v['purity'],chi2=r['pearson'],ndf=r['conditionalNdf'],mean=r['components'][species]['mean'],sigma=r['components'][species]['sigma'],he3Mean=r['components']['he3']['mean'],he3Sigma=r['components']['he3']['sigma'],sourceModel='gaussian_linear'))
 tsv(R/'examples.tsv',examples,'group species index rLo rHi root fitLo fitHi linearLo linearHi logLo logHi purity chi2 ndf mean sigma he3Mean he3Sigma sourceModel'.split());save(R/'examples.json',examples)
 shutil.copy2(OUT/'curves.tsv',R/'curves.tsv');shutil.copy2(OUT/'coefficients.json',R/'coefficients.json');selection=read(LATEST/'pid_selection.yaml');selection['purityEstimation']=dict(response=C['response'],background=C['background'],variation=C['variation'],errorPolicy=C['errorPolicy'],configuration='settings.yaml');save(R/'pid_selection.yaml',selection);shutil.copy2(HERE/'settings.yaml',R/'settings.yaml');save(R/'protected.json',protected)
 for p in HERE.glob('*'):
  if p.is_file():shutil.copy2(p,R/'source'/p.name)
 print('LINEAR_PREPARE',len(rows),'purity bins;',stats,flush=True)
if __name__=='__main__':main()
