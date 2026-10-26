#!/usr/bin/env python3
"""Isolated Gaussian-only refits, from frozen historical engines and explicit settings."""
import concurrent.futures,csv,difflib,json,pathlib,shutil,subprocess,sys,time
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE.parent))
from study import P,OUT,SP,read,save,sha
C=read(HERE/'settings.yaml');R=OUT/C['outputDirectory'];SRC=R/'source'
def tab(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def replace(s,a,b):
 assert a in s,a[:120]
 return s.replace(a,b)
def jobs():
 result=[]
 for species in C['species']:
  result.append(dict(name='low_'+species,band='low',species=species,reference=str(OUT/('he4_r061_001/full_gaussian_linear' if species=='he4' else 'candidate/low_'+species+'_gaussian_linear_nominal'))))
  for stage in (['high'] if species=='he4' else ['high','extension']):result.append(dict(name=stage+'_'+species,band='high',species=species,reference=str(OUT/('candidate/'+stage+'_'+species+'_gaussian_linear_nominal'))))
 return result
WARM=r'''
template<class E> void ApplyLinearSeeds(typename E::Request& request,const Slice& slice,const Settings& cfg,int target) {
 std::ifstream seeds(Text(cfg,"gaussianLinearSeeds").c_str());Require(bool(seeds),"Missing Gaussian+linear seeds");
 int index,s,count=0;double mean,sigma,yield;
 while(seeds>>index>>s>>mean>>sigma>>yield)if(index==slice.index){
  Require(s>=0&&s<4&&Finite(mean)&&Finite(sigma)&&sigma>0&&Finite(yield)&&yield>=0,"Invalid linear seed");
  typename E::Component& c=request.components[s];double shift=mean-c.mean,ratio=sigma/c.sigma;
  for(unsigned q=0;q<c.nodes.size();++q){c.nodes[q].mean+=shift;c.nodes[q].sigma*=ratio;}
  c.mean=mean;c.sigma=sigma;c.startWindowYield=c.active?yield:0;++count;
 }
 Require(seeds.eof()&&count==4,"Incomplete Gaussian+linear seeds");
 if(Text(cfg,"targetSpecies")==Text(cfg,"fixedTargetMeanSpecies")&&Text(cfg,"sampleHistogram")=="ALL"&&.5*(slice.lo+slice.hi)>Number(cfg,"fixedTargetMeanMinR")&&.5*(slice.lo+slice.hi)<Number(cfg,"fixedTargetMeanMaxR"))request.components[target].freeMean=false;
}
'''
def prepare():
 assert not R.exists(),'Never overwrite existing output'
 SRC.mkdir(parents=True)
 old=OUT/'he4_r061_001/source'
 for n in ['NuclearFitCommon.h','NuclearLognormalFitCommon.h','FitLocalNuclearCutoff.C']:
  original=(old/n).read_text();s=original
  if n.endswith('.C'):
   a=s.index('    if(cfg.count("targetWarmStartFile"))');b=s.index('    typename E::Result fit;',a)
   s=s[:a]+'    ApplyLinearSeeds<E>(request,slice,cfg,target);\n'+s[b:]
   s=replace(s,'template<class E> void Run(',WARM+'\ntemplate<class E> void Run(')
   s=replace(s,'out.meanError=fit.meanError[target];','out.meanError=req.components[target].freeMean?fit.meanError[target]:0;')
   # Diagnose against the original references, but account for a conditional fixed mean.
   s=replace(s,'ApplyLinearSeeds<E>(request,slice,cfg,target);','ApplyLinearSeeds<E>(request,slice,cfg,target);identityReference.components[target].freeMean=request.components[target].freeMean;')
   s=replace(s,'if(fit.meanError.size()>unsigned(s)&&fit.meanError[s]>=0)Scalar(json,fit.meanError[s]);else json<<"null";', 'if(!request.components[s].freeMean)json<<0;else if(fit.meanError.size()>unsigned(s)&&fit.meanError[s]>=0)Scalar(json,fit.meanError[s]);else json<<"null";json<<",\\\"meanFixed\\\":"<<Boolean(!request.components[s].freeMean)<<",\\\"sigmaFixed\\\":"<<Boolean(!request.components[s].freeSigma);')
   s=s.replace('Target response plus explicitly configured low-R He3 mean/width fitted to ALL; remaining nuisance responses fixed; anchors are optimization starts and identity reference, not extra likelihood','Gaussian-only; linear fit seeds; He4 target center fixed below configured R; width and yields free; original nuisance references and identity screens retained')
  (SRC/n).write_text(s);(SRC/(n+'.patch')).write_text(''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True),fromfile=str(old/n),tofile='gaussian_only/'+n)))
 high=OUT/'candidate/high_d_gaussian_linear_nominal/source/tools/nuclear_pid_cutoff/highrigidity/FitHighRigidity.C';original=high.read_text();s=original.replace('#include "../FitLocalNuclearCutoff.C"','#include "FitLocalNuclearCutoff.C"');s=replace(s,'bg=request.backgroundModel=="linear"?2:1','bg=request.backgroundModel=="none"?0:(request.backgroundModel=="linear"?2:1)');s=replace(s,'typename E::Result fit;std::string failure;','typename E::Request identityReference=request;ApplyLinearSeeds<E>(request,slice,cfg,target);\n    typename E::Result fit;std::string failure;');s=replace(s,'DiagnoseTarget<E>(request,fit','DiagnoseTarget<E>(identityReference,fit');(SRC/'FitHighRigidity.C').write_text(s);(SRC/'FitHighRigidity.C.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True))))
 protected={}
 for d in [OUT/'gaussian_linear_report_002',OUT/'he4_r061_001']:
  for p in d.iterdir():
   if p.is_file() and p.suffix in ['.json','.pdf','.tsv','.yaml']:protected[str(p)]=sha(p)
 for job in jobs():
  ref=pathlib.Path(job['reference'])
  for n in ['cutoff_fits.json','cutoff_fits.root','shapes.tsv','gates.tsv','settings.tsv','effective_settings.json','tag_shape_audit.json']:protected[str(ref/n)]=sha(ref/n)
 save(R/'protected.json',protected);save(R/'jobs.json',jobs());shutil.copy2(HERE/'settings.yaml',R/'settings.yaml')
 for p in HERE.iterdir():
  if p.is_file():shutil.copy2(p,SRC/p.name)
 print('PREPARED',R,flush=True)
def execute(job):
 dst=R/job['name'];ref=pathlib.Path(job['reference']);species=job['species']
 if (dst/'completed.json').exists():return
 assert not dst.exists(),str(dst)
 dst.mkdir()
 data={r['index']:r for r in read(ref/'cutoff_fits.json')['slices']};indices=set(data)
 linear={i:(r,ref) for i,r in data.items()}
 if job['name']=='low_he4':
  selected=[r for r in tab(OUT/'he4_r061_001/atlas.tsv') if r['band']=='gaussian' and r['model'].startswith('linear_')]
  cache={}
  for a in selected:
   f=pathlib.Path(a['root']).parent
   if str(f) not in cache:cache[str(f)]={r['index']:r for r in read(f/'cutoff_fits.json')['slices']}
   i=int(a['index']);linear[i]=(cache[str(f)][i],f)
 for n in ['shapes.tsv','gates.tsv','effective_settings.json','tag_shape_audit.json']:shutil.copy2(ref/n,dst/n)
 # Copy only requested He4 bins; all other original coverage remains unchanged.
 if species=='he4' and job['band']=='low':indices={i for i,r in data.items() if r['r']>C['minimumHe4RExclusive']}
 for n in ['shapes.tsv','gates.tsv']:
  rows=[r for r in tab(dst/n) if int(r['index']) in indices];write(dst/n,rows)
 seed=[];evidence=[]
 for i in sorted(indices):
  row,folder=linear[i]
  for k,species0 in enumerate(SP):
   c=row['components'][species0];seed.append(f"{i}\t{k}\t{c['mean']:.17g}\t{c['sigma']:.17g}\t{(row['columnWindowYields'][k] if len(row['columnWindowYields'])==4 else row['fitCount']/4):.17g}")
  evidence.append(dict(index=i,linearSeedAvailable=(len(row['columnWindowYields'])>=4),reference=str(folder/'cutoff_fits.json'),meanFixed=(species==C['fixedTargetMeanSpecies'] and job['band']=='low' and C['fixedTargetMeanMinR']<row['r']<C['fixedTargetMeanMaxR']),targetMean=row['components'][species]['mean'],targetSigma=row['components'][species]['sigma']))
 (dst/'linear_seeds.tsv').write_text('\n'.join(seed)+'\n');save(dst/'linear_seed_provenance.json',evidence)
 values=dict(line.split('\t',1) for line in (ref/'settings.tsv').read_text().splitlines());values.pop('targetWarmStartFile',None)
 values.update(backgroundModel='none',freeBackgroundSlope='0',gaussianLinearSeeds=str(dst/'linear_seeds.tsv'),fixedTargetMeanSpecies=C['fixedTargetMeanSpecies'],fixedTargetMeanMinR=str(C['fixedTargetMeanMinR']),fixedTargetMeanMaxR=str(C['fixedTargetMeanMaxR']))
 if job['name']!='low_he4':values.update(he3FreeMinR='0',he3FreeMaxR='0',he3FreeMean='0',he3FreeSigma='0')
 (dst/'settings.tsv').write_text('\n'.join(k+'\t'+v for k,v in sorted(values.items()))+'\n');effective=read(dst/'effective_settings.json');effective['gaussianOnlyStudy']=C;save(dst/'effective_settings.json',effective)
 previous=read(ref/'command.json');cache=json.loads(previous[-1].split('+(',1)[1].split(',',1)[0]);macro=SRC/('FitLocalNuclearCutoff.C' if job['band']=='low' else 'FitHighRigidity.C')
 command=['bash',str(P/'script/singularity_root_batch.sh'),str(P/C['runtimeMainconf']),'root4star','-l','-b','-q',str(macro)+'+('+','.join(json.dumps(str(p)) for p in [cache,dst,dst/'settings.tsv',dst/'shapes.tsv',dst/'gates.tsv'])+')'];save(dst/'command.json',command)
 print('START',job['name'],len(indices),flush=True);start=time.time()
 with (dst/'fit.log').open('w') as f:status=subprocess.run(command,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
 log=(dst/'fit.log').read_text();assert not status and 'NUCLEAR_CUTOFF_READBACK_PASS' in log and ('NUCLEAR_CUTOFF_COMPLETE' if job['band']=='low' else 'NUCLEAR_HIGH_R_COMPLETE') in log,dst
 new=read(dst/'cutoff_fits.json');assert new['responseModel']=='gaussian' and new['backgroundModel']=='none' and len(new['slices'])==len(indices)
 fixed=0
 for row in new['slices']:
  orig=data[row['index']];assert row['fitLo']==orig['fitLo'] and row['fitHi']==orig['fitHi']
  assert len(row['columnWindowYields']) in [0,4]
  e=next(e for e in evidence if e['index']==row['index'])
  if e['meanFixed']:
   fixed+=1;assert abs(row['components'][species]['mean']-e['targetMean'])<1e-10;assert row['components'][species]['meanFixed']
 save(dst/'completed.json',dict(executionComplete=True,rows=len(indices),fixedMeans=fixed,seconds=time.time()-start,products={p.name:sha(p) for p in dst.iterdir() if p.is_file()}));print('DONE',job['name'],round(time.time()-start,1),'s',flush=True)
if __name__=='__main__':
 if sys.argv[1]=='prepare':prepare()
 elif sys.argv[1]=='fit':
  todo=[j for j in jobs() if not (R/j['name']/'completed.json').exists()]
  # Serialize first compilation of each macro, then run independent studies.
  for band in ['low','high']:
   first=next((j for j in todo if j['band']==band),None)
   if first:execute(first);todo.remove(first)
  with concurrent.futures.ThreadPoolExecutor(max_workers=C['workers']) as pool:list(pool.map(execute,todo))
 else:raise ValueError(sys.argv[1])
