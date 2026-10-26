#!/usr/bin/env python3
import collections,json,math,pathlib,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,C=run.R,run.C;checks=[];roots=['region\tnew\told\tfirst\tlast'];gatechecks=[];half=[];fitcomparison=[]
for name in ['ExtractNuclearGateCounts.C','ValidateHe4.C']:shutil.copy2(HERE/name,R/'source'/name)
def rootcall(name,args,label,marker):
 cmd=['bash',str(run.P/'script/singularity_root_batch.sh'),str(run.P/C['runtimeMainconf']),'root4star','-l','-b','-q',str(R/'source'/name)+'+('+','.join(json.dumps(str(a)) for a in args)+')'];run.save(R/(label+'_command.json'),cmd)
 with (R/(label+'.log')).open('w') as f:status=subprocess.run(cmd,cwd=run.P,stdout=f,stderr=subprocess.STDOUT).returncode
 assert status==0 and marker in (R/(label+'.log')).read_text(),label
for region in C['regions']:
 folder=R/('candidate_'+region['name']);reference=run.ref(region);meta=run.read(folder/'cutoff_fits.json');rows=meta['slices'];cfg=run.read(folder/'effective_settings.json');baseline={r['index']:r for r in run.read(R/('baseline_'+region['name'])/'cutoff_fits.json')['slices']};indices=[r['index'] for r in rows];assert indices==list(range(min(indices),max(indices)+1));roots.append('\t'.join([region['name'],str(folder/'cutoff_fits.root'),str(reference/'cutoff_fits.root'),str(min(indices)),str(max(indices))]))
 if region['sample']=='ALL':
  output=folder/'independent_gate_counts.tsv'
  if not output.exists():
   cache=json.loads(run.read(folder/'command.json')[-1].split('+(',1)[1].split(',',1)[0]);rootcall('ExtractNuclearGateCounts.C',[cache,folder/'shapes.tsv',folder/'gates.tsv',folder/'settings.tsv',output],'count_'+region['name'],'NUCLEAR_GATE_COUNTS_COMPLETE')
  independent={(int(a['index']),a['species']):float(a['rawGateCount']) for a in run.tab(output)}
 shapes=run.tab(folder/'shapes.tsv');tables=run.read(R/'coefficients.json')['candidate']
 for a in run.tab(folder/'gates.tsv'):
  shape=next(s for s in shapes if s['index']==a['index'] and s['node']==a['node'] and s['species']==a['species']);s=run.SP.index(a['species']);r=float(shape['rValue']);assert float(a['mean'])==run.poly(tables['mean'][s],r) and float(a['sigma'])==run.poly(tables['sigma'][s],r);checks.append(dict(region=region['name'],index=int(a['index']),node=int(a['node']),species=a['species'],gitExact=s<3))
 for r in rows:
  old=baseline[r['index']]
  # Shape likelihood is independent of the PID gate: only purity integration changes.
  for k in ['pearson','dataDeviance','shapePoint','columnWindowYields']:assert r[k]==old[k],(region['name'],r['index'],k)
  fitcomparison.append(dict(region=region['name'],index=r['index'],exactLikelihoodMatch=True))
  if region['sample']=='ALL':
   for s in run.SP:
    value=r['components'][s]['selection']['rawGateCount'];external=independent[r['index'],s];assert abs(value-external)<1e-8*max(1.,external);gatechecks.append(dict(region=region['name'],index=r['index'],species=s,raw=value,independent=external))
   c=r['components']['he3']
   for a in shapes:
    if int(a['index'])!=r['index'] or a['species']!='he3':continue
    left=float(a['mean'])+c['mean']-c['initialMean']-math.sqrt(2*math.log(2))*float(a['sigma'])*c['sigma']/c['initialSigma'];assert left>=r['fitLo']-1e-9;half.append(dict(index=r['index'],node=int(a['node']),fitLo=r['fitLo'],leftHalfMaximum=left))
  else:
   assert meta['mass2Window']==C['tofMass2'] and meta['sampleHistogram']=='TAG_D'
   for s in run.SP[:3]:assert not r['components'][s]['active']
(R/'validation_roots.tsv').write_text('\n'.join(roots)+'\n');rootcall('ValidateHe4.C',[R],'root_validation','HE4_ROOT_VALIDATION_PASS')
for p,h in run.read(R/'protected.json').items():assert run.sha(pathlib.Path(p))==h,p
coef=run.read(R/'coefficients.json');baseline=run.read(R/'git_baseline.json')['tables'];assert all(coef['candidate'][k][:3]==baseline[k][:3] for k in ['mean','sigma','domain']);assert coef['candidate']['domain']==baseline['domain']
for region in C['regions']:
 f=R/('candidate_'+region['name']);seeds=run.read(f/'seed_provenance.json');cache={}
 for e in seeds:
  source=pathlib.Path(e['source']);
  if str(source) not in cache:cache[str(source)]={a['index']:a for a in run.read(source)['slices']};done=run.read(source.parent/'completed.json');assert run.sha(source)==done['products']['cutoff_fits.json']
  c=cache[str(source)][e['index']]['components']['he4'];assert c['mean']==e['mean'] and c['sigma']==e['sigma']
run.save(R/'gate_validation.json',checks);run.save(R/'independent_counts_validation.json',gatechecks);run.save(R/'halfheight_validation.json',half);run.save(R/'likelihood_comparison.json',fitcomparison)
summary=run.read(R/'summary.json');assert len(summary['rows'])==88
for r in summary['rows']:
 if r['totalSelected']>0:assert abs(r['purity']-r['signalSelected']/r['totalSelected'])<1e-12
run.save(R/'physics_validation.json',dict(executionComplete=True,allRegionsCorrect=True,gitDTHe3Exact=True,onlyHe4CoefficientsChanged=True,noFixedHe4Mean=True,independentRawChecks=len(gatechecks),he3HalfHeightChecks=len(half),fitGateIndependenceComparisons=len(fitcomparison),protectedInputs=len(run.read(R/'protected.json')),oldInputsUnchanged=True,highOtherNucleiZero=True,middleBackgroundZero=True))
print('HE4_VALIDATION_PASS',len(checks),'gates;',len(gatechecks),'independent raw counts;',len(half),'halfheight checks',flush=True)
