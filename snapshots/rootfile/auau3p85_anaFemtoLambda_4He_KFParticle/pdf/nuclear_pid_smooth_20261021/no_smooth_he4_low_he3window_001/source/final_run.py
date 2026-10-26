#!/usr/bin/env python3
"""Isolated none-background fits using copied historical math, never edits production."""
import concurrent.futures,csv,difflib,json,math,pathlib,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent));from study import P,OUT,SP,read,save,sha
CONFIG_PATH=pathlib.Path(sys.argv[2]).resolve() if len(sys.argv)>2 else HERE/'settings.yaml'
cfg=read(CONFIG_PATH);ROOT=OUT/cfg['outputDirectory'];SRC=ROOT/'source'
def replace(text,old,new):
 if old not in text:raise ValueError('Source contract changed: '+old[:80])
 return text.replace(old,new)
def prepare():
 if ROOT.exists():raise ValueError('Never overwrite '+str(ROOT))
 SRC.mkdir(parents=True);ref=OUT/'candidate/low_he4_gaussian_linear_nominal/source'
 changes={}
 for name in ['NuclearFitCommon.h','NuclearLognormalFitCommon.h','FitLocalNuclearCutoff.C']:
  original=(ref/name).read_text();s=original
  if name.endswith('.h'):
   s=replace(s,'bg=request.backgroundModel=="linear"?2:1','bg=request.backgroundModel=="none"?0:(request.backgroundModel=="linear"?2:1)')
   old='initial.push_back(request.backgroundWindowYield/count/(request.backgroundModel=="linear"?2:1));if(request.backgroundModel=="linear")initial.push_back(initial.back());'
   s=replace(s,old,'if(request.backgroundModel!="none"){'+old+'}')
   s=replace(s,'request.backgroundModel=="linear"||request.backgroundModel=="exponential"','request.backgroundModel=="none"||request.backgroundModel=="linear"||request.backgroundModel=="exponential"')
   s=replace(s,'else result.columnNames.push_back("BG_exp");','else if(fitted.backgroundModel=="exponential")result.columnNames.push_back("BG_exp");')
  else:
   s=replace(s,'bg=request.backgroundModel=="linear"?2:1','bg=request.backgroundModel=="none"?0:(request.backgroundModel=="linear"?2:1)')
   s=replace(s,'double per=count/(active+bg);','double per=(active+bg)>0?count/(active+bg):0;')
   old='typename E::Result fit;std::string failure;'
   new='''typename E::Request identityReference=request;
    if(cfg.count("targetWarmStartFile")){std::ifstream seeds(Text(cfg,"targetWarmStartFile").c_str());Require(bool(seeds),"Missing explicit warm starts");int index;double mean,sigma;while(seeds>>index>>mean>>sigma)if(index==slice.index){Require(Finite(mean)&&Finite(sigma)&&mean>0&&sigma>0,"Invalid warm start");typename E::Component& c=request.components[target];double shift=mean-c.mean,ratio=sigma/c.sigma;for(unsigned q=0;q<c.nodes.size();++q){c.nodes[q].mean+=shift;c.nodes[q].sigma*=ratio;}c.mean=mean;c.sigma=sigma;}Require(seeds.eof(),"Malformed warm starts");}
    typename E::Result fit;std::string failure;'''
   s=replace(s,old,new)
   s=replace(s,'DiagnoseTarget<E>(request,fit,slice,target','DiagnoseTarget<E>(identityReference,fit,slice,target')
  (SRC/name).write_text(s);changes[name]=''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),fromfile='historical/'+name,tofile='no_smooth_private/'+name))
  (SRC/(name+'.patch')).write_text(changes[name])
 for p in [CONFIG_PATH,pathlib.Path(__file__).resolve()]:shutil.copy2(p,SRC/p.name)
 save(ROOT/'source_patch_manifest.json',dict(originalSource=str(ref),originalSha256={n:sha(ref/n) for n in changes},patchedSha256={n:sha(SRC/n) for n in changes},policy='none removes background template columns, amplitudes and shape variables; no floor background or tiny nonzero amplitude. Warm starts preserve original identity reference.'))
 protected={str(OUT/n):h for n,h in read(OUT/'completed.json')['productsSha256'].items()}
 for response in cfg['responses']:
  for bg in cfg['backgroundsForControls']:
   folder=OUT/('candidate/low_he4_'+response+'_'+bg+'_nominal')
   for n in ['cutoff_fits.json','cutoff_fits.root','shapes.tsv','gates.tsv','settings.tsv','tag_shape_audit.json','effective_settings.json']:protected[str(folder/n)]=sha(folder/n)
 save(ROOT/'protected.json',protected)
 if 'fitWindowPolicy' in cfg:
  windows=[];factor=math.sqrt(2*math.log(2));policy=cfg['fitWindowPolicy']
  refs={response:read(OUT/('candidate/low_he4_'+response+'_linear_nominal/cutoff_fits.json'))['slices'] for response in cfg['responses']}
  peaks={response:read(OUT/('candidate/low_he3_'+response+'_linear_nominal/cutoff_fits.json'))['slices'] for response in cfg['responses']}
  for i in range(len(refs[cfg['responses'][0]])):
   oldrow=refs[cfg['responses'][0]][i];left=[];evidence=[]
   for response in cfg['responses']:
    for kind,row in [('fixed_he3',refs[response][i]),('all_he3',peaks[response][i])]:
     if kind=='all_he3' and not (row['numericalValid'] and row['targetFreeResponseSupported']):continue
     c=row['components']['he3'];m,w=c['mean'],c['sigma'];half=m-factor*w if response=='gaussian' else m*math.exp(-factor*w/m)
     assert math.isfinite(half) and half>0
     left.append(half);evidence.append(dict(response=response,reference=kind,mean=m,sigma=w,leftHalfMaximum=half))
   step=policy['alignmentKeVcm'];lo=math.floor((min([oldrow['fitLo']]+left)-policy['extraMarginKeVcm'])/step)*step
   lo=max(policy['minimumKeVcm'],lo);assert lo<=min(left) and lo<=oldrow['fitLo']
   windows.append(dict(index=i,r=oldrow['r'],oldFitLo=oldrow['fitLo'],fitLo=lo,fitHi=oldrow['fitHi'],references=evidence))
  save(ROOT/'fit_windows.json',windows)

def execute(response,bg,stage):
 ref=OUT/('candidate/low_he4_'+response+'_'+('linear' if bg=='none' else bg)+'_nominal');dst=ROOT/(stage+'_'+response+'_'+bg);dst.mkdir()
 for n in ['shapes.tsv','gates.tsv','tag_shape_audit.json','effective_settings.json']:shutil.copy2(ref/n,dst/n)
 if 'fitWindowPolicy' in cfg:
  windows={r['index']:r for r in read(ROOT/'fit_windows.json')};lines=list(csv.DictReader((dst/'shapes.tsv').open(),delimiter='\t'));fields=list(lines[0])
  for row in lines:row['fitLo']=format(windows[int(row['index'])]['fitLo'],'.17g');row['fitHi']=format(windows[int(row['index'])]['fitHi'],'.17g')
  with (dst/'shapes.tsv').open('w') as f:
   writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t',lineterminator='\n');writer.writeheader();writer.writerows(lines)
 values=dict(x.split('\t',1) for x in (ref/'settings.tsv').read_text().splitlines());values['backgroundModel']=bg
 if bg=='none':values['freeBackgroundSlope']='0'
 if stage=='control':
  for n in ['shapes.tsv','gates.tsv']:
   lines=(dst/n).read_text().splitlines();(dst/n).write_text('\n'.join([lines[0]]+[r for r in lines[1:] if int(r.split('\t')[0]) in cfg['controlIndices']])+'\n')
 if stage=='warm':
  prior=read(ROOT/('none_'+response+'_none/cutoff_fits.json'))['slices'];seeds=[]
  for r in prior:
   c=r['components']['he4']
   if r['numericalValid'] and c['mean']>0 and c['sigma']>0:seeds.append(str(r['index'])+'\t'+format(c['mean'],'.17g')+'\t'+format(c['sigma'],'.17g'))
  (dst/'target_seeds.tsv').write_text('\n'.join(seeds)+'\n');values['targetWarmStartFile']=str(dst/'target_seeds.tsv')
 (dst/'settings.tsv').write_text('\n'.join(k+'\t'+v for k,v in sorted(values.items()))+'\n')
 oldcmd=read(ref/'command.json');cache=json.loads(oldcmd[-1].split('+(',1)[1].split(',',1)[0])
 command=['bash',str(P/'script/singularity_root_batch.sh'),str(P/cfg['runtimeMainconf']),'root4star','-l','-b','-q',str(SRC/'FitLocalNuclearCutoff.C')+'+('+','.join(json.dumps(str(p)) for p in [cache,dst,dst/'settings.tsv',dst/'shapes.tsv',dst/'gates.tsv'])+')'];save(dst/'command.json',command)
 print('FIT_START',stage,response,bg,flush=True)
 with (dst/'fit.log').open('w') as f:status=subprocess.run(command,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
 log=(dst/'fit.log').read_text()
 if status or 'NUCLEAR_CUTOFF_READBACK_PASS' not in log or 'NUCLEAR_CUTOFF_COMPLETE' not in log:raise RuntimeError('Fit failed '+str(dst))
 data=read(dst/'cutoff_fits.json');rows=data['slices'];original={r['index']:r for r in read(ref/'cutoff_fits.json')['slices']}
 for row in rows:
  expected=windows[row['index']] if 'fitWindowPolicy' in cfg else original[row['index']]
  assert row['fitLo']==expected['fitLo'] and row['fitHi']==expected['fitHi']
  if stage=='control':
   for k in ['dataDeviance','pearson','columnWindowYields','shapePoint']:
    a,b=row[k],original[row['index']][k]
    aa,bb=(a,b) if isinstance(a,list) else ([a],[b])
    assert len(aa)==len(bb) and all(abs(x-y)<=1e-8*max(1.,abs(y)) for x,y in zip(aa,bb)),k
  if bg=='none':assert not row['columnWindowYields'] or len(row['columnWindowYields'])==4
 save(dst/'completed.json',dict(executionComplete=True,productionEnabled=False,response=response,background=bg,stage=stage,rows=len(rows),products={n:sha(dst/n) for n in ['cutoff_fits.json','cutoff_fits.root','shapes.tsv','gates.tsv','settings.tsv']}));print('FIT_COMPLETE',stage,response,bg,flush=True)
if __name__=='__main__':
 mode=sys.argv[1]
 if mode=='prepare':prepare()
 elif mode=='control':execute('gaussian','linear','control')
 elif mode=='control-lognormal':execute('lognormal','linear','control')
 elif mode=='none-gaussian':execute('gaussian','none','none')
 elif mode=='none-lognormal':execute('lognormal','none','none')
 elif mode=='none':
  if not (SRC/'FitLocalNuclearCutoff_C.so').exists():
   execute(cfg['responses'][0],'none','none');remaining=cfg['responses'][1:]
  else:remaining=cfg['responses']
  with concurrent.futures.ThreadPoolExecutor(max_workers=cfg['workers']) as pool:list(pool.map(lambda r:execute(r,'none','none'),remaining))
 elif mode=='baseline':
  jobs=[(r,b) for r in cfg['responses'] for b in cfg['backgroundsForControls']]
  with concurrent.futures.ThreadPoolExecutor(max_workers=cfg['workers']) as pool:list(pool.map(lambda j:execute(j[0],j[1],'baseline'),jobs))
 elif mode=='warm':
  jobs=[(r,b) for r in cfg['responses'] for b in cfg['backgroundsForControls']]
  with concurrent.futures.ThreadPoolExecutor(max_workers=cfg['workers']) as pool:list(pool.map(lambda j:execute(j[0],j[1],'warm'),jobs))
 else:raise ValueError(mode)
