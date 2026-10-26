#!/usr/bin/env python3
import csv,difflib,json,math,pathlib,shutil,subprocess,sys,time
HERE=pathlib.Path(__file__).resolve().parent;P=HERE.parents[1];sys.path.insert(0,str(P/'tools/nuclear_pid_smooth_20261021'));from study import OUT as OLD,read,save,sha,SP
C=read(HERE/'settings.yaml');R=P/C['outputBase'];PREVIOUS=P/C['previousOutputBase'];SRC=R/'source'
def tab(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def replace(s,a,b):
 assert a in s,a[:100]
 return s.replace(a,b)
def poly(cs,r):return sum(c*math.log10(r)**i for i,c in enumerate(cs))
def ref(region):return PREVIOUS/('candidate_low_linear' if region['name']=='low_fixed' else 'candidate_'+region['name'])
WARM=r'''
template<class E> void ApplyFixedMeanSeeds(typename E::Request& request,const Slice& slice,const Settings& cfg,int target) {
 std::ifstream seeds(Text(cfg,"fixedMeanSeedFile").c_str());Require(bool(seeds),"Missing explicit mean and response seeds");int index,s,count=0;double mean,sigma,yield;
 while(seeds>>index>>s>>mean>>sigma>>yield)if(index==slice.index){Require(s>=0&&s<4&&Finite(mean)&&Finite(sigma)&&sigma>0&&Finite(yield)&&yield>=0,"Invalid fitted seed");typename E::Component& c=request.components[s];double shift=mean-c.mean,ratio=sigma/c.sigma;for(unsigned q=0;q<c.nodes.size();++q){c.nodes[q].mean+=shift;c.nodes[q].sigma*=ratio;}c.mean=mean;c.sigma=sigma;c.startWindowYield=c.active?yield:0;++count;}
 Require(seeds.eof()&&count==4,"Incomplete fitted seeds");request.components[target].freeMean=false;
}
'''
def prepare():
 assert not R.exists(),'Never overwrite';SRC.mkdir(parents=True)
 protected=read(PREVIOUS/'protected.json')
 for p in PREVIOUS.rglob('*'):
  if p.is_file() and p.suffix in ['.json','.root','.pdf','.yaml','.tsv','.h','.C','.py'] and 'diagnostic' not in str(p.relative_to(PREVIOUS)):protected[str(p)]=sha(p)
 save(R/'protected.json',protected)
 for name in ['git_baseline.json','git_NuclearIdDeDxVsMom.h','coefficients.json','coefficients.yaml','curves.tsv','input_provenance.json']:shutil.copy2(PREVIOUS/name,R/name)
 shutil.copytree(PREVIOUS/'candidate_pid',R/'candidate_pid')
 for region in C['regions'][1:]:shutil.copytree(ref(region),R/('candidate_'+region['name']))
 for name in ['NuclearFitCommon.h','NuclearLognormalFitCommon.h','FitLocalNuclearCutoff.C']:shutil.copy2(PREVIOUS/'source'/name,SRC/name)
 p=SRC/'FitLocalNuclearCutoff.C';original=p.read_text();s=replace(original,'template<class E> void Run(',WARM+'\ntemplate<class E> void Run(')
 a=s.index('    if(cfg.count("targetWarmStartFile"))');b=s.index('    typename E::Result fit;',a);s=s[:a]+'    ApplyFixedMeanSeeds<E>(request,slice,cfg,target);identityReference.components[target].freeMean=false;\n'+s[b:]
 s=replace(s,'out.meanError=fit.meanError[target];','out.meanError=req.components[target].freeMean?fit.meanError[target]:0;')
 s=replace(s,'if(fit.meanError.size()>unsigned(s)&&fit.meanError[s]>=0)Scalar(json,fit.meanError[s]);else json<<"null";','if(!request.components[s].freeMean)json<<0;else if(fit.meanError.size()>unsigned(s)&&fit.meanError[s]>=0)Scalar(json,fit.meanError[s]);else json<<"null";json<<",\\\"meanFixed\\\":"<<Boolean(!request.components[s].freeMean)<<",\\\"sigmaFixed\\\":"<<Boolean(!request.components[s].freeSigma);')
 s=s.replace('Target response plus explicitly configured low-R He3 mean/width fitted to ALL; remaining nuisance responses fixed; anchors are optimization starts and identity reference, not extra likelihood','He4 mean fixed to preceding Gaussian+linear fit; sigma and yields free; no smooth BG; conditional fixed-center error excludes calibration uncertainty')
 p.write_text(s);(SRC/'FitLocalNuclearCutoff.C.patch').write_text(''.join(difflib.unified_diff(original.splitlines(True),s.splitlines(True))))
 shutil.copy2(HERE/'settings.yaml',R/'settings.yaml')
 for p in HERE.iterdir():
  if p.is_file():shutil.copy2(p,SRC/p.name)
 print('PREPARED',R,flush=True)
def fit():
 region=C['regions'][0];previous=ref(region);folder=R/'candidate_low_fixed';assert not folder.exists();folder.mkdir()
 for n in ['shapes.tsv','gates.tsv','tag_shape_audit.json','effective_settings.json','independent_gate_counts.tsv']:shutil.copy2(previous/n,folder/n)
 data=read(previous/'cutoff_fits.json')['slices'];seeds=[];evidence=[]
 for r in data:
  assert len(r['columnWindowYields'])==6
  for k,sp in enumerate(SP):
   c=r['components'][sp];seeds.append(f"{r['index']}\t{k}\t{c['mean']:.17g}\t{c['sigma']:.17g}\t{r['columnWindowYields'][k]:.17g}")
  evidence.append(dict(index=r['index'],source=str(previous/'cutoff_fits.json'),mean=r['components']['he4']['mean'],sigma=r['components']['he4']['sigma']))
 (folder/'fixed_mean_seeds.tsv').write_text('\n'.join(seeds)+'\n');save(folder/'seed_provenance.json',evidence)
 cfg=dict(line.split('\t',1) for line in (previous/'settings.tsv').read_text().splitlines());cfg.pop('targetWarmStartFile',None);cfg.update(backgroundModel='none',freeBackgroundSlope='0',fixedMeanSeedFile=str(folder/'fixed_mean_seeds.tsv'));(folder/'settings.tsv').write_text('\n'.join(k+'\t'+v for k,v in sorted(cfg.items()))+'\n');effective=read(folder/'effective_settings.json');effective['he4FixedMeanRevision']=C;save(folder/'effective_settings.json',effective)
 cache=json.loads(read(previous/'command.json')[-1].split('+(',1)[1].split(',',1)[0]);cmd=['bash',str(P/'script/singularity_root_batch.sh'),str(P/C['runtimeMainconf']),'root4star','-l','-b','-q',str(SRC/'FitLocalNuclearCutoff.C')+'+('+','.join(json.dumps(str(v)) for v in [cache,folder,folder/'settings.tsv',folder/'shapes.tsv',folder/'gates.tsv'])+')'];save(folder/'command.json',cmd)
 with (folder/'fit.log').open('w') as f:status=subprocess.run(cmd,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
 log=(folder/'fit.log').read_text();assert not status and 'NUCLEAR_CUTOFF_COMPLETE rows=19' in log and 'NUCLEAR_CUTOFF_READBACK_PASS' in log
 result=read(folder/'cutoff_fits.json');assert result['backgroundModel']=='none'
 for r,e in zip(result['slices'],evidence):
  c=r['components']['he4'];assert r['index']==e['index'] and c['mean']==e['mean'] and c['meanFixed'] and not c['sigmaFixed'] and c['meanError']==0
 save(folder/'completed.json',dict(executionComplete=True,rows=19,meanFixed=True,background='none',products={p.name:sha(p) for p in folder.iterdir() if p.is_file()}));print('FIT_COMPLETE 19 fixed-center Gaussian-only intervals',flush=True)
if __name__=='__main__':
 if sys.argv[1]=='prepare':prepare()
 elif sys.argv[1]=='fit':fit()
 else:raise ValueError(sys.argv[1])
