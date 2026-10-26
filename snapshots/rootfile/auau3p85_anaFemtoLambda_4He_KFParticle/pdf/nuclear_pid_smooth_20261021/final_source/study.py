#!/usr/bin/env python3
"""Replay immutable 17-Oct fits, changing only explicitly recorded polynomial gates."""
import argparse, concurrent.futures, csv, hashlib, json, math, pathlib, re, shutil, subprocess, sys, time
P=pathlib.Path(__file__).resolve().parents[2]
CONFIG=json.loads((P/'tools/nuclear_pid_smooth_20261021/settings.yaml').read_text())
OLD=P/CONFIG['referenceBase']
OUT=P/CONFIG['outputBase']
SP=['d','t','he3','he4']
def read(p): return json.loads(pathlib.Path(p).read_text())
def save(p,x): pathlib.Path(p).write_text(json.dumps(x,indent=2,ensure_ascii=False,allow_nan=False)+'\n')
def sha(p): return hashlib.sha256(pathlib.Path(p).read_bytes()).hexdigest()
def sources():
 low=read(OLD/'summary_free_sensitivity_001/summary.json')['inputStudies']
 runs=[('low',k,pathlib.Path(v)) for k,v in low.items()]
 high=OLD/'highrigidity_tof_20261017_001'
 for s in SP:
  for response in ['gaussian','lognormal']:
   for bg in ['linear','exponential']:
    p=high/(('coarse_native_he4_' if s=='he4' else 'native_'+s+'_')+response+'_'+bg+'_001')
    runs.append(('high',s+'_'+response+'_'+bg+'_nominal',p))
 return runs

def execute(job,mode,coeff=None):
 band,label,src=job
 dst=OUT/mode/(band+'_'+label)
 if dst.exists():
  if (dst/'completed.json').exists(): return str(dst)
  raise ValueError('Incomplete existing output '+str(dst))
 dst.mkdir(parents=True)
 shutil.copytree(src/'source',dst/'source',ignore=shutil.ignore_patterns('*.so','*.d','*.pcm','__pycache__'))
 ctx=read(src/'context.json'); cmd=read(src/'command.json')
 oldmacro=cmd[-1].split('+(',1)[0]; macro=str(dst/'source'/pathlib.Path(oldmacro).relative_to(src/'source'))
 # Reuse an already compiled, byte-identical historical snapshot when available.
 reference=OUT/'reproduce'/('low_he4_gaussian_linear_nominal' if band=='low' else 'high_he4_gaussian_linear_nominal')/'source'
 relative=pathlib.Path(oldmacro).relative_to(src/'source')
 if mode!='reproduce' and (reference/relative).exists():
  mathfiles=[p.relative_to(dst/'source') for p in (dst/'source').rglob('*') if p.suffix in ['.C','.h']]
  if all((reference/p).exists() and sha(reference/p)==sha(dst/'source'/p) for p in mathfiles):macro=str(reference/relative)
 # All mathematical sources are copied from the historical snapshot, never live code.
 protected={str(p):sha(p) for p in (src/'source').rglob('*') if p.is_file() and p.suffix in ['.C','.h','.py']}
 for name in ['settings.tsv','shapes.tsv','gates.tsv','tag_shape_audit.json','fixed_windows.json','effective_settings.json']:
  shutil.copy2(src/name,dst/name)
 rows=read(src/'cutoff_fits.json')['slices']
 if mode=='reproduce':
  wanted={rows[0]['index'],rows[len(rows)//2]['index'],rows[-1]['index']}
  for name in ['shapes.tsv','gates.tsv']:
   lines=(dst/name).read_text().splitlines();(dst/name).write_text('\n'.join([lines[0]]+[x for x in lines[1:] if int(x.split('\t')[0]) in wanted])+'\n')
 else:
  shapes=list(csv.DictReader((dst/'shapes.tsv').open(),delimiter='\t'))
  gates=[]
  for row in shapes:
   s=SP.index(row['species']); r=float(row['rValue']); x=math.log10(r)
   def poly(cs): return sum(c*x**i for i,c in enumerate(cs))
   mu,sd=poly(coeff['mean'][s]),poly(coeff['sigma'][s]);lo,hi=coeff['domain'][s]
   active=lo<=r<=hi and mu>0 and sd>0
   if band=='low': active=active and r<2.5
   gates.append([row['index'],row['node'],SP[s],format(mu,'.17g'),format(sd,'.17g'),str(int(active))])
  (dst/'gates.tsv').write_text('index\tnode\tspecies\tmean\tsigma\tactive\n'+'\n'.join('\t'.join(x) for x in gates)+'\n')
 cache=json.loads(cmd[-1].split('+(',1)[1].split(',',1)[0])
 main=P/CONFIG['runtimeMainconf']
 command=['bash',str(P/'script/singularity_root_batch.sh'),str(main),'root4star','-l','-b','-q',macro+'+('+','.join(json.dumps(str(x)) for x in [cache,dst,dst/'settings.tsv',dst/'shapes.tsv',dst/'gates.tsv'])+')']
 save(dst/'command.json',command);save(dst/'context.json',dict(ctx, historicalSource=str(src),studyMode=mode, baselinePolicy='Historical shape seeds, fit windows and math; only gates differ in candidate mode.'))
 print('START',mode,band,label,flush=True);t=time.time()
 with (dst/'fit.log').open('w') as f: status=subprocess.run(command,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
 if status or not (dst/'cutoff_fits.json').exists(): raise RuntimeError('Fit failed '+str(dst))
 data=read(dst/'cutoff_fits.json');checks=[]
 if mode=='reproduce':
  ref={r['index']:r for r in rows}
  for row in data['slices']:
   a=row['components'][label.split('_')[0]]['selection'];b=ref[row['index']]['components'][label.split('_')[0]]['selection']
   for key in ['signalSelected','totalSelected','rawGateCount','purity']:
    av,bv=a[key],b[key];ok=(av==bv) if av is None or bv is None else abs(av-bv)<=1.e-8*max(1.,abs(bv))
    checks.append(dict(index=row['index'],key=key,old=bv,new=av,passCheck=ok))
  if not all(x['passCheck'] for x in checks): raise RuntimeError('Reproduction mismatch '+str(dst))
 for p,h in protected.items():
  if sha(p)!=h: raise RuntimeError('Historical source changed '+p)
 save(dst/'completed.json',dict(executionComplete=True,productionEnabled=False,seconds=time.time()-t,historicalSource=str(src),historicalSourceSha256=protected,reproductionChecks=checks,products={n:sha(dst/n) for n in ['cutoff_fits.json','cutoff_fits.root','gates.tsv','shapes.tsv','settings.tsv']}))
 print('DONE',mode,band,label,round(time.time()-t,1),flush=True)
 return str(dst)

def main():
 parser=argparse.ArgumentParser();parser.add_argument('mode',choices=['reproduce','candidate','original']);parser.add_argument('--workers',type=int,default=2);args=parser.parse_args()
 OUT.mkdir(parents=True,exist_ok=True)
 jobs=sources()
 if args.mode=='original': jobs=[j for j in jobs if j[1].endswith('_nominal')]
 if args.mode=='reproduce': jobs=[j for j in jobs if '_gaussian_linear_nominal' in j[1]]
 coeff=read(OUT/'coefficients.json')['candidate' if args.mode=='candidate' else 'original'] if args.mode!='reproduce' else None
 with concurrent.futures.ThreadPoolExecutor(max_workers=args.workers) as pool:
  results=list(pool.map(lambda j:execute(j,args.mode,coeff),jobs))
 save(OUT/(args.mode+'_completed.json'),dict(executionComplete=True,runs=results,productionEnabled=False))
if __name__=='__main__': main()
