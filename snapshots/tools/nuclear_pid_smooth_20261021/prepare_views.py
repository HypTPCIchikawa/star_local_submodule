#!/usr/bin/env python3
import csv,math
from study import OUT,OLD,SP,read,save,P
from summarize import tsv
out=[]
for band,source in [('low',OLD/'summary_free_sensitivity_001/summary.json'),('high',OLD/'highrigidity_tof_20261017_001/summary_final_001/summary.json')]:
 data=read(source)
 for r in data['aggregateRows' if band=='low' else 'bins']:
  stat=r.get('modelDiagnostics',{}).get('gaussian_linear_nominal',{}).get('conditionalErrorWeightedBound')
  out.append(dict(mode='historical',band=band,species=r['species'],rLo=r['rLo'],rHi=r['rHi'],purity=r['purityNominal'],statConditionalBound=stat,purityLow=r['purityLow'],purityHigh=r['purityHigh'],quality=r.get('modelsSupported',r['quality'])))
fields='mode band species rLo rHi purity statConditionalBound purityLow purityHigh quality'.split();tsv(OUT/'historical_summary.tsv',out,fields)
coeff=read(OUT/'coefficients.json');settings=read(P/'tools/nuclear_pid_smooth_20261021/settings.yaml');rows=[]
for i in range(settings['auditGridCount']):
 r=settings['auditRMin']*(settings['auditRMax']/settings['auditRMin'])**(i/(settings['auditGridCount']-1));x=math.log10(r);row=[r]
 for model in ['original','candidate']:
  for s in range(4):
   for kind in ['mean','sigma']:row.append(sum(v*x**j for j,v in enumerate(coeff[model][kind][s])))
 rows.append(row)
with (OUT/'curves.tsv').open('w') as f:
 f.write('r\t'+'\t'.join(m+'_'+s+'_'+k for m in ['original','candidate'] for s in SP for k in ['mean','sigma'])+'\n')
 for row in rows:f.write('\t'.join(format(v,'.17g') for v in row)+'\n')
print('VIEWS_INPUTS_COMPLETE')
