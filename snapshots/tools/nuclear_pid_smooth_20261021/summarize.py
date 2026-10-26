#!/usr/bin/env python3
"""Use the historical aggregation and quality thresholds without a new mask."""
import argparse,collections,csv,importlib.util,json,math,pathlib,shutil,sys
from study import P,OLD,OUT,SP,read,save,sha,sources

def load_old():
 dst=OUT/'summary_source';dst.mkdir(exist_ok=True)
 for name in ['summarize_cutoff.py','run_cutoff.py']:
  src=OLD/'summary_free_sensitivity_001/source'/name
  if not (dst/name).exists():shutil.copy2(src,dst/name)
  if sha(src)!=sha(dst/name):raise ValueError('Historical summary changed')
 sys.path.insert(0,str(dst));spec=importlib.util.spec_from_file_location('historical_summary',dst/'summarize_cutoff.py');m=importlib.util.module_from_spec(spec);spec.loader.exec_module(m);return m

def tsv(p,rows,fields):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=fields,delimiter='\t',extrasaction='ignore');w.writeheader()
  for r in rows:w.writerow({k:('nan' if v is None else int(v) if isinstance(v,bool) else v) for k,v in r.items() if k in fields})

def summarize(mode):
 settings=read(P/'tools/nuclear_pid_smooth_20261021/settings.yaml');width=settings['aggregateWidth'];threshold=settings['purityThreshold']
 old=load_old();groups=collections.defaultdict(dict);native=[];provenance={};expected=collections.defaultdict(set)
 rawpath=OUT/(mode+'_gate_counts/gate_counts.tsv');raw={}
 for r in csv.DictReader(rawpath.open(),delimiter='\t'):
  raw[int(r['rBin']),r['species']]={k:float(r[k]) for k in ['rawGateCount','rawGateCountLowerFullBins','rawGateCountUpperTouchingBins']}
 if len(raw)!=640:raise ValueError('Incomplete independent gate counts')
 jobs=sources()
 if mode=='candidate':
  jobs+= [('extension',p.name[len('extension_'):],p) for p in sorted((OUT/mode).glob('extension_*'))]
 for band,label,src in jobs:
  if mode=='original' and not label.endswith('_nominal'):continue
  folder=src if band=='extension' else OUT/mode/(band+'_'+label)
  if band=='extension':band='high'
  done=read(folder/'completed.json')
  for n,h in done['products'].items():
   if sha(folder/n)!=h:raise ValueError('Completed input changed')
  log=(folder/'fit.log').read_text()
  if 'NUCLEAR_CUTOFF_READBACK_PASS' not in log or not re_complete(log,band):raise ValueError('Missing ROOT completion')
  s=label.split('_')[0];model=label[len(s)+1:];expected[band,s].add(model)
  cfg=read(folder/'effective_settings.json');cfg['fitWindowMode']='target_local'
  audit={r['index']:r for r in read(folder/'tag_shape_audit.json')};bygroup=collections.defaultdict(list)
  for row in read(folder/'cutoff_fits.json')['slices']:
   # round-to-nearest is inappropriate inside the five native bins: floor with edge tolerance.
   group=int(math.floor((row['lo']-(settings['auditRMin'] if band=='low' else settings['tofThreshold']))/width+1e-7))
   bygroup[group].append(row)
   if label.endswith('_nominal'):
    c=row['components'][s];native.append(dict(linearLo=max(row['fitLo'],c['mean']-settings['plots']['linearFocusHalfSigma']*c['sigma']),linearHi=min(row['fitHi'],c['mean']+settings['plots']['linearFocusHalfSigma']*c['sigma']),supported=old.aggregate_one([row],s,audit,cfg,raw if band=='low' else None)['supported'],mode=mode,band=band,species=s,model=model,index=row['index'],rLo=row['lo'],rHi=row['hi'],fitLo=row['fitLo'],fitHi=row['fitHi'],mean=c['mean'],sigma=c['sigma'],root=str(folder/'cutoff_fits.root'),purity=c['selection'].get('purity'),numericalValid=row['numericalValid']))
  for group,rows in bygroup.items():
   value=old.aggregate_one(rows,s,audit,cfg,raw if band=='low' else None)
   if band=='high':
    for row in rows:
     v=old.aggregate_one([row],s,audit,cfg);drop=row['components'][s].get('conditionalDropDeltaD')
     if not v['supported'] or drop is None or drop<cfg['support']['minDropDeltaD']:
      value['reasons']=sorted(set(value['reasons']+['one_or_more_native_fits_unassessed']));value['supported']=False
   count=1 if band=='high' and s=='he4' else 5
   if len(rows)!=count:raise ValueError('Incomplete R aggregation')
   if model in groups[band,s,group]:raise ValueError('Overlapping R/model inputs')
   groups[band,s,group][model]=value
  provenance[str(folder/'completed.json')]=sha(folder/'completed.json')
 # Validate coverage independently of model support or purity values.
 for band in ['low','high']:
  for species in SP:
   maxR=settings['tofThreshold'] if band=='low' else settings['auditRMax'] if mode=='candidate' or species=='he4' else settings['extensionRMin']
   minR=settings['auditRMin'] if band=='low' else settings['tofThreshold']
   expectedGroups=set(range(round((maxR-minR)/width)))
   actualGroups={g for b,s,g in groups if b==band and s==species}
   if actualGroups!=expectedGroups:raise ValueError('Missing/overlapping R coverage '+band+' '+species)
 bins=[]
 for (band,s,g),models in sorted(groups.items()):
  allok=set(models)==expected[band,s] and all(m['supported'] for m in models.values())
  nominal=models['gaussian_linear_nominal'];good=[m['purity'] for m in models.values() if m['supported'] and m['purity'] is not None]
  ngood=[m['purity'] for k,m in models.items() if k.endswith('_nominal') and m['supported'] and m['purity'] is not None]
  nquality=all(m['supported'] for k,m in models.items() if k.endswith('_nominal'))
  lo=(settings['auditRMin'] if band=='low' else settings['tofThreshold'])+g*width
  row=dict(mode=mode,band=band,species=s,rLo=lo,rHi=lo+width,purity=nominal['purity'],statConditionalBound=nominal['conditionalErrorWeightedBound'],purityLow=min(good) if good else None,purityHigh=max(good) if good else None,quality=allok,nominalQuality=nquality,nominalLow=min(ngood) if ngood else None,nominalHigh=max(ngood) if ngood else None,rawSelected=nominal['rawSelected'],modelsSupported=sum(m['supported'] for m in models.values()),modelCount=len(models),status='supported_above_50' if allok and min(good)>=threshold else 'supported_below_50' if allok else 'unassessed',reasons=sorted({v for m in models.values() for v in m['reasons']}),models=models)
  bins.append(row)
 summary=dict(schemaVersion=1,productionEnabled=False,mode=mode,purityDefinition='S/(S+B); sum signal and total before dividing',statErrorMeaning='Historical conditional yield error bound; excludes response/identity uncertainties',modelEnvelope='Supported fits only. quality requires every declared historical model, without new thresholds.',bins=bins,provenance=provenance)
 save(OUT/(mode+'_summary.json'),summary);save(OUT/(mode+'_atlas.json'),native)
 tsv(OUT/(mode+'_summary.tsv'),bins,'mode band species rLo rHi purity statConditionalBound purityLow purityHigh quality nominalQuality nominalLow nominalHigh rawSelected modelsSupported modelCount status'.split())
 tsv(OUT/(mode+'_atlas.tsv'),native,'band species model index rLo rHi fitLo fitHi mean sigma root purity numericalValid supported linearLo linearHi'.split())
 for s in SP:
  for band in ['low','high']:
   rows=[r for r in bins if r['species']==s and r['band']==band];good=[r for r in rows if r['quality']];print(mode,s,band,'supported',len(good),'/',len(rows),'min',min([r['purityLow'] for r in good],default=None),flush=True)
 return summary

def re_complete(log,band):return ('NUCLEAR_CUTOFF_COMPLETE' if band=='low' else 'NUCLEAR_HIGH_R_COMPLETE') in log
if __name__=='__main__':
 parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['candidate','original']);summarize(parser.parse_args().mode)
