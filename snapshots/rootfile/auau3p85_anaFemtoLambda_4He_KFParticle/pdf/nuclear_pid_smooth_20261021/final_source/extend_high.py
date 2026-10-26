#!/usr/bin/env python3
"""Extend d/t/He3 study beyond historical R=2.5 with the frozen high-R engine."""
import concurrent.futures,csv,json,math,pathlib,shutil,subprocess,time
from study import P,OUT,OLD,SP,read,save,sha
cfg=read(P/'tools/nuclear_pid_smooth_20261021/settings.yaml')
inputContext=read(P/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001/context.json')
for sp,suffix in [('d','D'),('t','T'),('he3','HE3'),('he4','D')]:
 recorded=[inputContext['settings']['tag'+suffix+'Mass2'+side] for side in ['Min','Max']]
 if any(abs(a-b)>1e-12 for a,b in zip(recorded,cfg['tofMass2'][sp])):raise ValueError('Saved TOF window mismatch')
coeff=read(OUT/'coefficients.json');high=OLD/'highrigidity_tof_20261017_001'
def poly(cs,r):return sum(c*math.log10(r)**j for j,c in enumerate(cs))
def run(job):
 s,response,bg=job;label=s+'_'+response+'_'+bg+'_nominal';src=high/('native_he4_'+response+'_'+bg+'_001');dst=OUT/'candidate'/('extension_'+label)
 if dst.exists():raise ValueError('Output exists')
 dst.mkdir();shutil.copytree(src/'source',dst/'source',ignore=shutil.ignore_patterns('*.so','*.d','*.pcm','__pycache__'))
 hist=read(src/'completed.json')
 for n,h in hist['products'].items():
  if n in ['shapes.tsv','settings.tsv','gates.tsv'] or n.startswith('source/') and n.endswith(('.C','.h')):
   if sha(src/n)!=h:raise ValueError('Historical source changed')
 oldcmd=read(src/'command.json');macro=str(dst/'source'/pathlib.Path(oldcmd[-1].split('+(',1)[0]).relative_to(src/'source'))
 reference=OUT/'reproduce/high_he4_gaussian_linear_nominal/source';relative=pathlib.Path(oldcmd[-1].split('+(',1)[0]).relative_to(src/'source')
 mathfiles=[p.relative_to(dst/'source') for p in (dst/'source').rglob('*') if p.suffix in ['.C','.h']]
 if all((reference/p).exists() and sha(reference/p)==sha(dst/'source'/p) for p in mathfiles):macro=str(reference/relative)
 settings=dict(x.split('\t',1) for x in (src/'settings.tsv').read_text().splitlines());settings['targetSpecies']=s;settings['sampleHistogram']={'d':'TAG_D','t':'TAG_T','he3':'TAG_HE3'}[s];settings['mass2Min'],settings['mass2Max']=map(str,cfg['tofMass2'][s]);settings['rigidityMin']=str(cfg['extensionRMin']);settings['rigidityMax']=str(cfg['auditRMax'])
 (dst/'settings.tsv').write_text('\n'.join(k+'\t'+v for k,v in sorted(settings.items()))+'\n')
 width=float(settings['nativeRigidityBinWidth']);first=round(cfg['extensionRMin']/width)+1;last=round(cfg['auditRMax']/width);iTarget=SP.index(s);half=float(read(src/'effective_settings.json')['allFit']['fitWindowSigmaHalfWidth']);align=read(src/'effective_settings.json')['allFit']['fitWindowAlignmentKeV']
 shapeRows=[];gates=[];audits=[];windows=[]
 for index,rbin in enumerate(range(first,last+1)):
  lo=(rbin-1)*width;hi=rbin*width;mid=(lo+hi)/2;offset=math.sqrt(3/5)*(hi-lo)/2;nodes=[(mid-offset,5/18),(mid,8/18),(mid+offset,5/18)]
  references=[(poly(coeff['original']['mean'][iTarget],mid),poly(coeff['original']['sigma'][iTarget],mid))]
  for r,w in nodes:references.append((poly(coeff['candidate']['mean'][iTarget],r),poly(coeff['candidate']['sigma'][iTarget],r)))
  fitlo=max(0.,math.floor(min(m-half*w for m,w in references)/align)*align);fithi=min(120.,math.ceil(max(m+half*w for m,w in references)/align)*align)
  for node,(r,w) in enumerate(nodes):
   for k,sp in enumerate(SP):
    mean=poly(coeff['original']['mean'][k],mid);sd=poly(coeff['original']['sigma'][k],mid)
    shapeRows.append([index,rbin,lo,hi,node,w,r,sp,mean,sd,int(sp==s),fitlo,fithi]);gates.append([index,node,sp,poly(coeff['candidate']['mean'][k],r),poly(coeff['candidate']['sigma'][k],r),1])
  audits.append(dict(index=index,rBin=rbin,lo=lo,hi=hi,failures={sp:[] if sp==s else ['provider_only_nuisance_response_unvalidated'] for sp in SP},referencePolicy='Git original initializes target only; target is free in data; other provider-only responses are unvalidated.'))
  windows.append(dict(index=index,rBin=rbin,lo=lo,hi=hi,fitLo=fitlo,fitHi=fithi))
 def table(name,header,rows):
  (dst/name).write_text('\t'.join(header.split())+'\n'+'\n'.join('\t'.join(format(v,'.17g') if isinstance(v,float) else str(v) for v in row) for row in rows)+'\n')
 table('shapes.tsv','index rBin rLo rHi node weight rValue species mean sigma supported fitLo fitHi',shapeRows);table('gates.tsv','index node species mean sigma active',gates)
 save(dst/'tag_shape_audit.json',audits);save(dst/'fixed_windows.json',windows);effective=read(src/'effective_settings.json');effective.update(rigidityMin=cfg['extensionRMin'],rigidityMax=cfg['auditRMax']);save(dst/'effective_settings.json',effective)
 context=read(src/'context.json');context.update(targetSpecies=s,sampleHistogram=settings['sampleHistogram'],historicalSource=str(src),stage='extension_native',exactProductionTofWindow=True);save(dst/'context.json',context)
 cache=json.loads(oldcmd[-1].split('+(',1)[1].split(',',1)[0]);command=['bash',str(P/'script/singularity_root_batch.sh'),str(P/'config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml'),'root4star','-l','-b','-q',macro+'+('+','.join(json.dumps(str(p)) for p in [cache,dst,dst/'settings.tsv',dst/'shapes.tsv',dst/'gates.tsv'])+')'];save(dst/'command.json',command);print('EXTENSION_START',label,flush=True);t=time.time()
 with (dst/'fit.log').open('w') as f:status=subprocess.run(command,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
 if status or 'NUCLEAR_HIGH_R_COMPLETE rows='+str(last-first+1) not in (dst/'fit.log').read_text():raise RuntimeError('Extension failed '+str(dst))
 save(dst/'completed.json',dict(executionComplete=True,productionEnabled=False,seconds=time.time()-t,products={n:sha(dst/n) for n in ['cutoff_fits.json','cutoff_fits.root','gates.tsv','shapes.tsv','settings.tsv']}));print('EXTENSION_DONE',label,round(time.time()-t,1),flush=True)
with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:list(pool.map(run,[(s,r,b) for s in SP[:3] for r in ['gaussian','lognormal'] for b in ['linear','exponential']]))
save(OUT/'extension_completed.json',dict(executionComplete=True,productionEnabled=False,runs=12))
