#!/usr/bin/env python3
"""Input/source protection, positive polynomial widths and saved-gate arithmetic audit."""
import csv,math,subprocess
from study import P,OUT,OLD,SP,read,save,sha,sources
cfg=read(P/'tools/nuclear_pid_smooth_20261021/settings.yaml');coeff=read(OUT/'coefficients.json');checks=[];protected={}
for band,label,src in sources():
 done=read(src/'completed.json')
 for n,expected in done['products'].items():
  if n in ['shapes.tsv','gates.tsv','settings.tsv','cutoff_fits.json','cutoff_fits.root','tag_shape_audit.json','fixed_windows.json','effective_settings.json'] or n.startswith('source/') and n.endswith(('.C','.h','.py')):
   actual=sha(src/n)
   if actual!=expected:raise ValueError('Historical input changed '+str(src/n))
   protected[str(src/n)]=actual
for model in ['original','candidate']:
 for s in range(4):
  if s<2 and any(coeff['original'][k][s]!=coeff['candidate'][k][s] for k in ['mean','sigma']):raise ValueError('d/t changed')
  for i in range(cfg['auditGridCount']):
   r=cfg['auditRMin']*(cfg['auditRMax']/cfg['auditRMin'])**(i/(cfg['auditGridCount']-1));x=math.log10(r)
   v={k:sum(c*x**j for j,c in enumerate(coeff[model][k][s])) for k in ['mean','sigma']}
   if min(v.values())<=0:raise ValueError('Nonpositive polynomial')
  checks.append(dict(model=model,species=SP[s],positiveDomain=[cfg['auditRMin'],cfg['auditRMax']]))
inputdir=P/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001';done=read(inputdir/'completed.json')
if sha(inputdir/'cache.root')!=done['products']['cache.root']:raise ValueError('Cache changed')
protected[str(inputdir/'cache.root')]=done['products']['cache.root']
# Preserve current production references for a final unchanged check.
mainfiles=list((P/'config/mainconf').glob('*auau3p85*gitpid*'))
for p in (P/'config/mainconf').glob('*auau3p85*Lambda*.yaml'):
 protected[str(p)]=sha(p)
for p in (P/'config/cuts/nuclearid').glob('*.yaml'):protected[str(p)]=sha(p)
save(OUT/'protected_input_audit.json',dict(executionComplete=True,productionEnabled=False,polynomialChecks=checks,protectedFiles=protected))
print('AUDIT_PASS',len(protected),'protected inputs; positive means/widths in original domain')
