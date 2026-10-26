#!/usr/bin/env python3
import concurrent.futures,csv,difflib,json,math,pathlib,re,shutil,subprocess,sys,time
HERE=pathlib.Path(__file__).resolve().parent;P=HERE.parents[1];sys.path.insert(0,str(P/'tools/nuclear_pid_smooth_20261021'));from study import OUT as OLD,read,save,sha,SP
C=read(HERE/'settings.yaml');R=P/C['outputBase'];SRC=R/'source'
def tab(p):return list(csv.DictReader(p.open(),delimiter='\t'))
def write(p,rows):
 with p.open('w') as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]),delimiter='\t',lineterminator='\n');w.writeheader();w.writerows(rows)
def replace(s,a,b):
 assert a in s,a[:100]
 return s.replace(a,b)
def poly(cs,r):return sum(c*math.log10(r)**i for i,c in enumerate(cs))
def ref(region):return OLD/('candidate/high_he4_gaussian_linear_nominal' if region['sample']!='ALL' else 'he4_r061_001/full_gaussian_linear')
def prepare():
 assert not R.exists(),'Never overwrite previous output'
 SRC.mkdir(parents=True)
 git=subprocess.run(['git','show',C['referenceGit']+':'+C['referenceHeader']],cwd=P,capture_output=True,check=True).stdout;(R/'git_NuclearIdDeDxVsMom.h').write_bytes(git)
 s=git.decode();original={}
 for kind,array in [('mean','kMean'),('sigma','kSigma'),('domain','kPRange')]:
  block=s.split(array+'[kNSpecies] = {',1)[1].split('};',1)[0];original[kind]=[[float(x) for x in row.split(',')] for row in re.findall(r'\{([^{}]+)\}',block)]
  assert len(original[kind])==4
 assert original==read(OLD/'coefficients.json')['original']
 save(R/'git_baseline.json',dict(revision=C['referenceGit'],header=C['referenceHeader'],sha256=sha(R/'git_NuclearIdDeDxVsMom.h'),tables=original))
 frozen=OLD/'he4_r061_001/source'
 for n in ['NuclearFitCommon.h','NuclearLognormalFitCommon.h','FitLocalNuclearCutoff.C']:shutil.copy2(frozen/n,SRC/n)
 high=OLD/'candidate/high_he4_gaussian_linear_nominal/source/tools/nuclear_pid_cutoff/highrigidity/FitHighRigidity.C';s=high.read_text();s=s.replace('#include "../FitLocalNuclearCutoff.C"','#include "FitLocalNuclearCutoff.C"')
 s=replace(s,'c.active=Finite(acceptance[s])&&acceptance[s]>floor;','c.active=s==target&&Finite(acceptance[s])&&acceptance[s]>floor;')
 s=replace(s,'if(c.active) ++active;else {bool ratio;','if(c.active) ++active;else if(s==target){bool ratio;')
 # Preserve the original identity reference when using explicit linear-derived target starts.
 s=replace(s,'typename E::Result fit;std::string failure;', '''typename E::Request identityReference=request;
    if(cfg.count("targetWarmStartFile")){std::ifstream seeds(Text(cfg,"targetWarmStartFile").c_str());Require(bool(seeds),"Missing high-R seeds");int index;double mean,sigma;while(seeds>>index>>mean>>sigma)if(index==slice.index){typename E::Component& c=request.components[target];double shift=mean-c.mean,ratio=sigma/c.sigma;for(unsigned q=0;q<c.nodes.size();++q){c.nodes[q].mean+=shift;c.nodes[q].sigma*=ratio;}c.mean=mean;c.sigma=sigma;}Require(seeds.eof(),"Malformed seeds");}
    typename E::Result fit;std::string failure;''')
 s=replace(s,'DiagnoseTarget<E>(request,fit','DiagnoseTarget<E>(identityReference,fit')
 s=s.replace('Target mean/mode and width fitted to this conditional sample; nuisance responses conditional on their recorded TAG reference; not intrinsic widths','He4 Gaussian mean/width free plus linear only; d/t/He3 amplitudes disabled by explicit user model; zeros are model assumptions, not measured absence')
 (SRC/'FitHighRigidity.C').write_text(s);(SRC/'FitHighRigidity.C.patch').write_text(''.join(difflib.unified_diff(high.read_text().splitlines(True),s.splitlines(True))))
 protected={}
 for region in C['regions']:
  for n in ['cutoff_fits.json','cutoff_fits.root','shapes.tsv','gates.tsv','settings.tsv','tag_shape_audit.json','effective_settings.json']:protected[str(ref(region)/n)]=sha(ref(region)/n)
 for d in [OLD/'gaussian_only_001',OLD/'gaussian_linear_report_002',OLD/'he4_r061_001']:
  for p in d.iterdir():
   if p.is_file() and p.suffix in ['.pdf','.json','.yaml','.tsv']:protected[str(p)]=sha(p)
 save(R/'protected.json',protected);shutil.copy2(HERE/'settings.yaml',R/'settings.yaml')
 for p in HERE.iterdir():
  if p.is_file():shutil.copy2(p,SRC/p.name)
 print('PREPARED',R,flush=True)
def execute(region,stage):
 folder=R/(stage+'_'+region['name'])
 if (folder/'completed.json').exists():return
 assert not folder.exists(),folder;folder.mkdir()
 reference=ref(region);data={r['index']:r for r in read(reference/'cutoff_fits.json')['slices']};indices={i for i,r in data.items() if region['rMin']<=r['r']<region['rMax']};assert indices
 shapes=[a for a in tab(reference/'shapes.tsv') if int(a['index']) in indices];write(folder/'shapes.tsv',shapes)
 coefficient=read(R/'git_baseline.json')['tables'] if stage=='baseline' else read(R/'coefficients.json')['candidate']
 gates=[]
 for a in shapes:
  s=SP.index(a['species']);rr=float(a['rValue']);mean,sigma=poly(coefficient['mean'][s],rr),poly(coefficient['sigma'][s],rr);lo,hi=coefficient['domain'][s];active=lo<=rr<=hi and (s!=3 or rr>C['minimumRExclusive'])
  gates.append(dict(index=a['index'],node=a['node'],species=a['species'],mean=format(mean,'.17g'),sigma=format(sigma,'.17g'),active=int(active)))
 write(folder/'gates.tsv',gates)
 for n in ['effective_settings.json','tag_shape_audit.json']:shutil.copy2(reference/n,folder/n)
 values=dict(line.split('\t',1) for line in (reference/'settings.tsv').read_text().splitlines());values.pop('targetWarmStartFile',None);values.update(backgroundModel=region['background'],freeBackgroundSlope='0')
 seeds=[];selected={}
 if region['sample']=='ALL':
  for a in tab(OLD/'he4_r061_001/atlas.tsv'):
   if a['band']=='gaussian' and a['model'].startswith('linear_') and int(a['index']) in indices:selected[int(a['index'])]=pathlib.Path(a['root']).parent
 cache={};evidence=[]
 for i in sorted(indices):
  origin=selected.get(i,reference)
  if str(origin) not in cache:cache[str(origin)]={r['index']:r for r in read(origin/'cutoff_fits.json')['slices']}
  row=cache[str(origin)][i];comp=row['components']['he4'];seeds.append(f"{i}\t{comp['mean']:.17g}\t{comp['sigma']:.17g}");evidence.append(dict(index=i,source=str(origin/'cutoff_fits.json'),mean=comp['mean'],sigma=comp['sigma']))
 (folder/'target_seeds.tsv').write_text('\n'.join(seeds)+'\n');save(folder/'seed_provenance.json',evidence);values['targetWarmStartFile']=str(folder/'target_seeds.tsv')
 (folder/'settings.tsv').write_text('\n'.join(k+'\t'+v for k,v in sorted(values.items()))+'\n');effective=read(folder/'effective_settings.json');effective['he4OnlyStudy']=C;effective['fitWindowMode']='target_local';save(folder/'effective_settings.json',effective)
 previous=read(reference/'command.json');cachePath=json.loads(previous[-1].split('+(',1)[1].split(',',1)[0]);macro=SRC/('FitLocalNuclearCutoff.C' if region['sample']=='ALL' else 'FitHighRigidity.C')
 cmd=['bash',str(P/'script/singularity_root_batch.sh'),str(P/C['runtimeMainconf']),'root4star','-l','-b','-q',str(macro)+'+('+','.join(json.dumps(str(p)) for p in [cachePath,folder,folder/'settings.tsv',folder/'shapes.tsv',folder/'gates.tsv'])+')'];save(folder/'command.json',cmd);print('START',stage,region['name'],len(indices),flush=True);t=time.time()
 with (folder/'fit.log').open('w') as f:status=subprocess.run(cmd,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
 log=(folder/'fit.log').read_text();assert status==0 and 'NUCLEAR_CUTOFF_READBACK_PASS' in log and ('NUCLEAR_CUTOFF_COMPLETE' if region['sample']=='ALL' else 'NUCLEAR_HIGH_R_COMPLETE') in log,folder
 results=read(folder/'cutoff_fits.json')['slices'];assert len(results)==len(indices)
 for row in results:
  assert row['fitLo']==data[row['index']]['fitLo'] and row['fitHi']==data[row['index']]['fitHi']
  if region['sample']!='ALL':
   for s in SP[:3]:assert not row['components'][s]['active']
   assert not row['columnWindowYields'] or row['columnWindowYields'][:3]==[0,0,0]
 save(folder/'completed.json',dict(executionComplete=True,rows=len(results),region=region,stage=stage,seconds=time.time()-t,products={n:sha(folder/n) for n in ['cutoff_fits.json','cutoff_fits.root','shapes.tsv','gates.tsv','settings.tsv','target_seeds.tsv','seed_provenance.json']}));print('DONE',stage,region['name'],round(time.time()-t,1),flush=True)
if __name__=='__main__':
 mode=sys.argv[1]
 if mode=='prepare':prepare()
 elif mode in ['baseline','candidate']:
  for region in C['regions']:execute(region,mode)
 else:raise ValueError(mode)
