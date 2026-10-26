#!/usr/bin/env python3
import bisect,concurrent.futures,csv,difflib,json,math,pathlib,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent));from study import P,OUT,read,save,sha
CONFIG=HERE/'settings.yaml';cfg=read(CONFIG);ROOT=OUT/cfg['outputDirectory'];SRC=ROOT/'source';PREVIOUS=OUT/'no_smooth_he4_low_he3window_001'
def table(path):return list(csv.DictReader(path.open(),delimiter='\t'))
def write(path,rows):
 with path.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def oldmean(rows,species,r):
 points=sorted((float(a['rValue']),float(a['mean'])) for a in rows if a['species']==species);xs=[a[0] for a in points];j=max(1,min(len(xs)-1,bisect.bisect_left(xs,r)));x,y=points[j-1];xx,yy=points[j];return y+(r-x)*(yy-y)/(xx-x)
def initial_mean(rows,species,r):
 anchors=cfg['userPeakAnchors'];a,b=anchors;end=cfg['lowPeakAdjustmentMaxR'];slope=(b[species]-a[species])/(b['r']-a['r'])
 if r<=b['r']:return a[species]+slope*(r-a['r'])
 if r>=end:return oldmean(rows,species,r)
 # Cubic Hermite join for optimization starts only. PID selection remains cubic in log10 R.
 width=end-b['r'];t=(r-b['r'])/width;delta=cfg['nativeRigidityBinWidth'];y=oldmean(rows,species,end);derivative=(oldmean(rows,species,end+delta/2)-oldmean(rows,species,end-delta/2))/delta
 return (2*t**3-3*t*t+1)*b[species]+(t**3-2*t*t+t)*width*slope+(-2*t**3+3*t*t)*y+(t**3-t*t)*width*derivative
def prepare_windows():
 result=[]
 for w in read(PREVIOUS/'fit_windows.json'):
  if w['r']<=cfg['minimumRExclusive']:continue
  q=dict(w);q['previousFitLo']=w['fitLo'];q['references']=[]
  if w['r']<cfg['lowPeakAdjustmentMaxR']:
   lower=w['oldFitLo']
   for response in cfg['responses']:
    original=table(OUT/('candidate/low_he4_'+response+'_linear_nominal/shapes.tsv'));widths={r['index']:r['components']['he3']['sigma'] for r in read(OUT/('candidate/low_he3_'+response+'_linear_nominal/cutoff_fits.json'))['slices']}
    for row in original:
     if int(row['index'])!=w['index'] or row['species']!='he3':continue
     r=float(row['rValue']);mu=initial_mean(original,'he3',r);sigma=max(widths[w['index']],float(row['sigma']));left=mu-math.sqrt(2*math.log(2))*sigma if response=='gaussian' else mu*math.exp(-math.sqrt(2*math.log(2))*sigma/mu)
     lower=min(lower,math.floor((left-cfg['halfHeightExtraMarginKeVcm'])/cfg['fitWindowAlignmentKeVcm'])*cfg['fitWindowAlignmentKeVcm']);q['references'].append(dict(response=response,node=int(row['node']),mean=mu,sigma=sigma,leftHalfMaximum=left))
   q['fitLo']=lower
  result.append(q)
 save(ROOT/'fit_windows.json',result)
def prepare():
 if ROOT.exists():raise ValueError('Never overwrite '+str(ROOT))
 SRC.mkdir(parents=True)
 for name in ['NuclearFitCommon.h','NuclearLognormalFitCommon.h','FitLocalNuclearCutoff.C']:
  original=(PREVIOUS/'source'/name).read_text();s=original
  if name.endswith('.C'):
   s=s.replace('c.freeMean=freeTarget&&s==target;c.freeSigma=freeTarget&&s==target;', 'bool lowHe3=s==2&&.5*(slice.lo+slice.hi)>Number(cfg,"he3FreeMinR")&&.5*(slice.lo+slice.hi)<Number(cfg,"he3FreeMaxR");c.freeMean=(freeTarget&&s==target)||(lowHe3&&Flag(Integer(cfg,"he3FreeMean")));c.freeSigma=(freeTarget&&s==target)||(lowHe3&&Flag(Integer(cfg,"he3FreeSigma")));')
   assert s!=original
   s=s.replace('Only target response mean/mode and width fitted to ALL; other projected TAG responses fixed; TAG is seed/identity reference, not additional likelihood','Target response plus explicitly configured low-R He3 mean/width fitted to ALL; remaining nuisance responses fixed; anchors are optimization starts and identity reference, not extra likelihood')
  (SRC/name).write_text(s);(SRC/(name+'.patch')).write_text(''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),fromfile='previous/'+name,tofile='r061/'+name)))
 shutil.copy2(CONFIG,SRC/CONFIG.name);shutil.copy2(HERE/'run.py',SRC/'run.py')
 protected=read(PREVIOUS/'protected.json');protected.update({str(PREVIOUS/n):h for n,h in read(PREVIOUS/'completed.json')['productsSha256'].items()})
 for n in ['coefficients.json','candidate_summary.json','candidate_summary.tsv']:protected[str(OUT/n)]=sha(OUT/n)
 save(ROOT/'protected.json',protected);save(ROOT/'source_manifest.json',dict(previous=str(PREVIOUS/'source'),sourceSha256={p.name:sha(p) for p in SRC.iterdir() if p.suffix in ['.h','.C','.py','.yaml']},productionEnabled=False))
 prepare_windows()
def execute(response,bg,stage):
 ref=OUT/('candidate/low_he4_'+response+'_'+('linear' if bg=='none' else bg)+'_nominal');dst=ROOT/(stage+'_'+response+'_'+bg);dst.mkdir()
 for n in ['tag_shape_audit.json','effective_settings.json']:shutil.copy2(ref/n,dst/n)
 original=table(ref/'shapes.tsv');windows={w['index']:w for w in read(PREVIOUS/'fit_windows.json')};allhe3={r['index']:r for r in read(OUT/('candidate/low_he3_'+response+'_linear_nominal/cutoff_fits.json'))['slices']};shapes=[];changes=[]
 for item in original:
  row=dict(item);i=int(row['index']);r=float(row['rValue']);mid=(float(row['rLo'])+float(row['rHi']))/2
  if mid<=cfg['minimumRExclusive']:continue
  if stage.startswith('pilot') and i not in cfg['pilotIndices']:continue
  if r<cfg['lowPeakAdjustmentMaxR'] and row['species'] in ['he3','he4']:
   row['mean']=format(initial_mean(original,row['species'],r),'.17g')
   if row['species']=='he3':row['sigma']=format(allhe3[i]['components']['he3']['sigma'],'.17g')
  shapes.append(row)
 # Both responses use exactly the same corrected half-height window.
 lower={w['index']:w['fitLo'] for w in read(ROOT/'fit_windows.json')}
 for row in shapes:
  i=int(row['index']);row['fitLo']=format(lower[i],'.17g');row['fitHi']=format(windows[i]['fitHi'],'.17g')
  if row['node']=='1' and row['species'] in ['he3','he4']:changes.append(dict(index=i,r=float(row['rValue']),species=row['species'],mean=float(row['mean']),sigma=float(row['sigma']),fitLo=lower[i],fitHi=windows[i]['fitHi']))
 write(dst/'shapes.tsv',shapes);save(dst/'peak_initial_values.json',changes)
 # Reuse exact gate coefficients; strict lower limit is applied at each quadrature node.
 indices={int(row['index']) for row in shapes};gates=[]
 for item in table(ref/'gates.tsv'):
  if int(item['index']) not in indices:continue
  row=dict(item)
  if row['species']=='he4':
   match=next(x for x in shapes if x['index']==row['index'] and x['node']==row['node'] and x['species']=='he4')
   row['active']=str(int(row['active']=='1' and float(match['rValue'])>cfg['minimumRExclusive']))
  gates.append(row)
 write(dst/'gates.tsv',gates)
 values=dict(line.split('\t',1) for line in (ref/'settings.tsv').read_text().splitlines());values.update(backgroundModel=bg,he3FreeMinR=str(cfg['minimumRExclusive']),he3FreeMaxR=str(cfg['lowPeakAdjustmentMaxR']),he3FreeMean=str(int(cfg['nuisanceHe3FreeMean'])),he3FreeSigma=str(int(cfg['nuisanceHe3FreeSigma'])))
 if bg=='none':values['freeBackgroundSlope']='0'
 if stage=='warm':
  seeds=[]
  for r in read(ROOT/('full_'+response+'_none/cutoff_fits.json'))['slices']:
   c=r['components']['he4']
   if r['numericalValid'] and c['mean']>0 and c['sigma']>0:seeds.append(str(r['index'])+'\t'+format(c['mean'],'.17g')+'\t'+format(c['sigma'],'.17g'))
  (dst/'target_seeds.tsv').write_text('\n'.join(seeds)+'\n');values['targetWarmStartFile']=str(dst/'target_seeds.tsv')
 (dst/'settings.tsv').write_text('\n'.join(k+'\t'+v for k,v in sorted(values.items()))+'\n')
 effective=read(dst/'effective_settings.json');effective['r061Study']=cfg;save(dst/'effective_settings.json',effective)
 command0=read(ref/'command.json');cache=json.loads(command0[-1].split('+(',1)[1].split(',',1)[0]);command=['bash',str(P/'script/singularity_root_batch.sh'),str(P/cfg['runtimeMainconf']),'root4star','-l','-b','-q',str(SRC/'FitLocalNuclearCutoff.C')+'+('+','.join(json.dumps(str(p)) for p in [cache,dst,dst/'settings.tsv',dst/'shapes.tsv',dst/'gates.tsv'])+')'];save(dst/'command.json',command)
 print('START',stage,response,bg,flush=True)
 with (dst/'fit.log').open('w') as f:status=subprocess.run(command,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
 log=(dst/'fit.log').read_text();assert not status and 'NUCLEAR_CUTOFF_COMPLETE' in log and 'NUCLEAR_CUTOFF_READBACK_PASS' in log,dst
 result=read(dst/'cutoff_fits.json')['slices'];assert len(result)==len(indices)
 for row in result:
  assert row['r']>cfg['minimumRExclusive'] and row['fitLo']==lower[row['index']] and row['fitHi']==windows[row['index']]['fitHi']
 save(dst/'completed.json',dict(rows=len(result),productionEnabled=False,products={n:sha(dst/n) for n in ['cutoff_fits.json','cutoff_fits.root','shapes.tsv','gates.tsv','settings.tsv']}));print('COMPLETE',stage,response,bg,flush=True)
if __name__=='__main__':
 mode=sys.argv[1]
 if mode=='prepare':prepare()
 elif mode=='pilot-first':execute('gaussian','none','pilot')
 else:
  if mode=='pilot':jobs=[(r,b,'pilot') for r in cfg['responses'] for b in cfg['backgrounds'] if (r,b)!=('gaussian','none')]
  elif mode=='pilot-anchor':jobs=[(r,b,'pilot_anchor') for r in cfg['responses'] for b in cfg['backgrounds']]
  elif mode=='repair-log-full':jobs=[('lognormal',b,'full') for b in cfg['backgrounds']]
  elif mode=='repair-log-warm':jobs=[('lognormal',b,'warm') for b in cfg['backgrounds'] if b!='none']
  elif mode=='full':jobs=[(r,b,'full') for r in cfg['responses'] for b in cfg['backgrounds']]
  elif mode=='warm':jobs=[(r,b,'warm') for r in cfg['responses'] for b in cfg['backgrounds'] if b!='none']
  else:raise ValueError(mode)
  if jobs and not (SRC/'FitLocalNuclearCutoff_C.so').exists():execute(*jobs.pop(0))
  with concurrent.futures.ThreadPoolExecutor(max_workers=cfg['workers']) as pool:list(pool.map(lambda j:execute(*j),jobs))
