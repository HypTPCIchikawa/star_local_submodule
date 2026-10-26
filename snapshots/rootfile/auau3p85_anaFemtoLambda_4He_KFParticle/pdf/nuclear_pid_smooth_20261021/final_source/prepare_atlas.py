#!/usr/bin/env python3
import argparse,csv
from study import P,OUT,SP,read,sha
from summarize import load_old,tsv
p=argparse.ArgumentParser();p.add_argument('band',choices=['low','high']);args=p.parse_args();cfg=read(P/'tools/nuclear_pid_smooth_20261021/settings.yaml');old=load_old();rows=[]
folders=sorted((OUT/'candidate').glob(args.band+'_*_nominal'))
if args.band=='high':folders+=sorted((OUT/'candidate').glob('extension_*_nominal'))
for folder in folders:
 if not (folder/'completed.json').exists():raise ValueError('Incomplete nominal fit '+str(folder))
 context=read(folder/'context.json');s=context['targetSpecies'];label='_'.join(context[k] for k in ['response','background','variation']);settings=read(folder/'effective_settings.json');settings['fitWindowMode']='target_local';audit={r['index']:r for r in read(folder/'tag_shape_audit.json')}
 for r in read(folder/'cutoff_fits.json')['slices']:
  c=r['components'][s];v=old.aggregate_one([r],s,audit,settings)
  rows.append(dict(band=args.band,species=s,model=label,index=r['index'],rLo=r['lo'],rHi=r['hi'],fitLo=r['fitLo'],fitHi=r['fitHi'],mean=c['mean'],sigma=c['sigma'],root=str(folder/'cutoff_fits.root'),purity=c['selection'].get('purity'),numericalValid=r['numericalValid'],supported=v['supported'],linearLo=max(r['fitLo'],c['mean']-cfg['plots']['linearFocusHalfSigma']*c['sigma']),linearHi=min(r['fitHi'],c['mean']+cfg['plots']['linearFocusHalfSigma']*c['sigma'])))
out=OUT/('atlas_'+args.band+'_001');out.mkdir();tsv(out/'candidate_atlas.tsv',rows,'band species model index rLo rHi fitLo fitHi mean sigma root purity numericalValid supported linearLo linearHi'.split());print(out)
