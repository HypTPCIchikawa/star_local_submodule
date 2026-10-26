#!/usr/bin/env python3
import collections,json,math,pathlib,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,C,OLD=run.R,run.C,run.PREVIOUS;means=[];half=[];counts=[];gates=[];roots=['region\tnew\told\tfirst\tlast'];comparisons=[]
for region in C['regions']:
 folder=R/('candidate_'+region['name']);before=run.ref(region);rows=run.read(folder/'cutoff_fits.json')['slices'];prior={a['index']:a for a in run.read(before/'cutoff_fits.json')['slices']};indices=[r['index'] for r in rows];assert indices==list(range(min(indices),max(indices)+1));roots.append('\t'.join([region['name'],str(folder/'cutoff_fits.root'),str(before/'cutoff_fits.root'),str(min(indices)),str(max(indices))]));shape=run.tab(folder/'shapes.tsv');tables=run.read(R/'coefficients.json')['candidate']
 if region['name']!='low_fixed':
  for name,value in run.read(before/'completed.json')['products'].items():assert run.sha(folder/name)==value
  assert run.read(folder/'cutoff_fits.json')==run.read(before/'cutoff_fits.json')
 else:
  for n in ['shapes.tsv','gates.tsv','independent_gate_counts.tsv','tag_shape_audit.json']:assert run.sha(folder/n)==run.sha(before/n)
  for n,h in run.read(folder/'completed.json')['products'].items():assert run.sha(folder/n)==h
 for a in run.tab(folder/'gates.tsv'):
  s=run.SP.index(a['species']);q=next(q for q in shape if q['index']==a['index'] and q['node']==a['node'] and q['species']==a['species']);rr=float(q['rValue']);assert float(a['mean'])==run.poly(tables['mean'][s],rr) and float(a['sigma'])==run.poly(tables['sigma'][s],rr);gates.append(dict(region=region['name'],index=int(a['index']),node=int(a['node']),species=a['species'],gitExact=s<3))
 if region['sample']=='ALL':independent={(int(a['index']),a['species']):float(a['rawGateCount']) for a in run.tab(folder/'independent_gate_counts.tsv')}
 for r in rows:
  if region['name']=='low_fixed':
   c=r['components']['he4'];old=prior[r['index']]['components']['he4'];assert c['meanFixed'] and not c['sigmaFixed'] and c['meanError']==0 and c['mean']==old['mean'];assert len(r['shapePoint'])==len(prior[r['index']]['shapePoint'])-1;assert len(r['columnWindowYields'])==4
   nparams=len(r['shapePoint'])+sum(y>r['fitCount']*1e-12 for y in r['columnWindowYields']);assert r['conditionalNdf']==round((r['fitHi']-r['fitLo'])/.1)-nparams
   means.append(dict(index=r['index'],rLo=r['lo'],rHi=r['hi'],mean=c['mean'],previousLinearMean=old['mean'],sigma=c['sigma'],previousLinearSigma=old['sigma'],purity=c['selection']['purity'],previousPurity=old['selection']['purity'],meanFixed=True,conditionalMeanError=0,ndf=r['conditionalNdf']))
  if region['sample']=='ALL':
   for s in run.SP:
    x=r['components'][s]['selection']['rawGateCount'];external=independent[r['index'],s];assert abs(x-external)<1e-8*max(1.,external);counts.append(dict(region=region['name'],index=r['index'],species=s,raw=x,independent=external))
   c=r['components']['he3']
   for a in shape:
    if int(a['index'])!=r['index'] or a['species']!='he3':continue
    left=float(a['mean'])+c['mean']-c['initialMean']-math.sqrt(2*math.log(2))*float(a['sigma'])*c['sigma']/c['initialSigma'];assert left>=r['fitLo']-1e-9;half.append(dict(index=r['index'],node=int(a['node']),leftHalfMaximum=left,fitLo=r['fitLo']))
for n in ['coefficients.json','coefficients.yaml','curves.tsv','git_baseline.json','git_NuclearIdDeDxVsMom.h']:assert run.sha(R/n)==run.sha(OLD/n)
for n in ['NuclearIdDeDxVsMom.h','he4_coefficients.patch']:assert run.sha(R/'candidate_pid'/n)==run.sha(OLD/'candidate_pid'/n)
current=run.read(R/'summary.json');prior=run.read(OLD/'summary.json');unchanged=[r for r in current['rows'] if r['region']!='low_fixed'];assert unchanged==[r for r in prior['rows'] if r['region']!='low_linear']
for a,b in zip(current['rows'][:4],prior['rows'][:4]):comparisons.append(dict(rLo=a['rLo'],rHi=a['rHi'],previousLinearPurity=b['purity'],fixedMeanGaussianOnlyPurity=a['purity'],quality=a['quality'],reasons=a['reasons'],conditionalError=a['statConditionalBound']))
(R/'validation_roots.tsv').write_text('\n'.join(roots)+'\n');shutil.copy2(HERE/'ValidateHe4.C',R/'source/ValidateHe4.C');cmd=['bash',str(run.P/'script/singularity_root_batch.sh'),str(run.P/C['runtimeMainconf']),'root4star','-l','-b','-q',str(R/'source/ValidateHe4.C')+'+('+json.dumps(str(R))+')'];run.save(R/'root_validation_command.json',cmd)
with (R/'root_validation.log').open('w') as f:status=subprocess.run(cmd,cwd=run.P,stdout=f,stderr=subprocess.STDOUT).returncode
assert not status and 'HE4_ROOT_VALIDATION_PASS' in (R/'root_validation.log').read_text()
for p,h in run.read(R/'protected.json').items():assert run.sha(pathlib.Path(p))==h,p
assert len(means)==19 and len(unchanged)==84
run.save(R/'fixed_center_validation.json',means);run.save(R/'low_comparison.json',comparisons);run.save(R/'halfheight_validation.json',half);run.save(R/'independent_counts_validation.json',counts);run.save(R/'gate_validation.json',gates)
run.save(R/'physics_validation.json',dict(executionComplete=True,fixedMeanIntervals=19,allMeansExactlyPreviousLinear=True,he4SigmaFree=True,zeroConditionalFixedMeanError=True,noBackgroundBelow2=True,ndfMatchesFreeParameters=True,independentRawChecks=len(counts),he3HalfHeightChecks=len(half),unchangedOtherRegionsNativeFits=180,unchangedAggregateBins=84,pidCoefficientsUnchanged=True,oldInputsUnchanged=True,protectedInputs=len(run.read(R/'protected.json')),highOtherNucleiZero=True,rootAll199Readback=True))
print('REVISION_VALIDATION_PASS',len(means),'fixed means;',len(unchanged),'unchanged aggregate bins;',len(counts),'raw checks',flush=True)
