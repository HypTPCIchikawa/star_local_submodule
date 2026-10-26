#!/usr/bin/env python3
import collections,csv,difflib,json,math,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,C=run.R,run.C
summary=run.read(R/'summary.json');native=run.read(R/'native.json');examples=run.read(R/'examples.json');validation=run.read(R/'pdf_validation.json');assert validation['allPagesReadback']
fixed=[];windows=[];seedchecks=[];inputchecks=[];verifiedLinearSources=set()
for job in run.jobs():
 folder=R/job['name'];meta=run.read(folder/'cutoff_fits.json');data=meta['slices'];assert meta['responseModel']=='gaussian' and meta['backgroundModel']=='none'
 if job['band']=='low':assert meta['sampleHistogram']=='ALL'
 else:assert meta['mass2Window']==run.read(HERE.parent/'settings.yaml')['tofMass2'][job['species']]
 settings=dict(line.split('\t',1) for line in (folder/'settings.tsv').read_text().splitlines());assert settings['backgroundModel']=='none' and settings['freeBackgroundSlope']=='0'
 source=pathlib.Path(job['reference']);assert run.sha(folder/'gates.tsv')==run.sha(source/'gates.tsv')
 for n,h in run.read(folder/'completed.json')['products'].items():assert run.sha(folder/n)==h
 provenance={e['index']:e for e in run.read(folder/'linear_seed_provenance.json')};shapes=run.tab(folder/'shapes.tsv');rawseed={}
 for line in (folder/'linear_seeds.tsv').read_text().splitlines():
  i,sp,mean,sigma,yield0=line.split('\t');rawseed[int(i),run.SP[int(sp)]]=(float(mean),float(sigma),float(yield0))
 cache={}
 for row in data:
  e=provenance[row['index']];sourceJson=pathlib.Path(e['reference'])
  if str(sourceJson) not in cache:cache[str(sourceJson)]={a['index']:a for a in run.read(sourceJson)['slices']}
  old=cache[str(sourceJson)][row['index']]
  if str(sourceJson) not in verifiedLinearSources:
   priorDone=run.read(sourceJson.parent/'completed.json')
   assert run.sha(sourceJson)==priorDone['products']['cutoff_fits.json']
   assert run.sha(sourceJson.parent/'cutoff_fits.root')==priorDone['products']['cutoff_fits.root']
   verifiedLinearSources.add(str(sourceJson))
  assert len(row['columnWindowYields']) in [0,4]
  for s in run.SP:
   before=old['components'][s];seed=rawseed[row['index'],s];assert seed[0]==before['mean'] and seed[1]==before['sigma']
   expected=old['columnWindowYields'][run.SP.index(s)] if len(old['columnWindowYields'])==4 else old['fitCount']/4;assert seed[2]==expected
  seedchecks.append(dict(job=job['name'],index=row['index'],exactLinearSeeds=True,linearSeedAvailable=e['linearSeedAvailable']))
  if e['meanFixed']:
   c=row['components'][job['species']];assert c['meanFixed'] and c['meanError']==0 and not c['sigmaFixed'] and c['mean']==e['targetMean'];fixed.append(dict(index=row['index'],rLo=row['lo'],rHi=row['hi'],mean=c['mean'],linearMean=e['targetMean'],sigma=c['sigma']))
  if job['name']=='low_he4':
   c=row['components']['he3']
   for a in shapes:
    if int(a['index'])!=row['index'] or a['species']!='he3':continue
    left=float(a['mean'])+(c['mean']-c['initialMean'])-math.sqrt(2*math.log(2))*float(a['sigma'])*(c['sigma']/c['initialSigma']);assert left>=row['fitLo']-1e-9;windows.append(dict(index=row['index'],node=int(a['node']),leftHalfMaximum=left,fitLo=row['fitLo']))
assert len(fixed)==9 and len(native)==1579 and len(examples)==1579 and len(summary['rows'])==364
for row in summary['rows']:
 v=row['assessment'];assert row['quality']==v['supported'];assert row['purity']==v['purity'];assert row['signalSelected']==v['signalSelected'] and row['totalSelected']==v['totalSelected'];assert row['statConditionalBound']==v['conditionalErrorWeightedBound']
 if row['totalSelected']>0:assert abs(row['purity']-row['signalSelected']/row['totalSelected'])<1e-12
for p,h in run.read(R/'protected.json').items():assert run.sha(pathlib.Path(p))==h,p
run.save(R/'fixed_center_validation.json',fixed);run.save(R/'halfheight_validation.json',windows);run.save(R/'warm_seed_validation.json',seedchecks)
# Refresh source snapshots, including a complete diff against the frozen historical math.
for p in HERE.iterdir():
 if p.is_file():shutil.copy2(p,R/'source'/p.name)
old=run.OUT/'he4_r061_001/source'
for name in ['NuclearFitCommon.h','NuclearLognormalFitCommon.h','FitLocalNuclearCutoff.C']:
 (R/'source'/(name+'.patch')).write_text(''.join(difflib.unified_diff((old/name).read_text().splitlines(True),(R/'source'/name).read_text().splitlines(True),fromfile=str(old/name),tofile='gaussian_only/'+name)))
# Assemble convenient complete and He4-only volumes. Never overwrite an existing merge.
commands=[];products=dict(validation['products']);order=['fit_low_'+sp+'.pdf' for sp in C['species']]+['fit_high_'+sp+'.pdf' for sp in C['species']]
for name,inputs,expected in [('gaussian_all_intervals.pdf',order,1579),('fit_gaussian_he4_all.pdf',['fit_low_he4.pdf','fit_high_he4.pdf'],199)]:
 pdf=R/name
 if not pdf.exists():
  cmd=['gs','-q','-dSAFER','-dNOPAUSE','-dBATCH','-dPDFSTOPONERROR','-dUseCropBox','-sDEVICE=pdfwrite','-sOutputFile='+str(pdf)]+[str(R/n) for n in inputs];commands.append(cmd)
  with (R/(name+'.merge.log')).open('w') as f:status=subprocess.run(cmd,stdout=f,stderr=subprocess.STDOUT).returncode
  assert status==0
 cmd=['gs','-dSAFER','-dNOPAUSE','-dBATCH','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)];commands.append(cmd);p=subprocess.run(cmd,capture_output=True,text=True);log=p.stdout+p.stderr;(R/(name+'.readback.log')).write_text(log);match=re.search(r'Processing pages 1 through (\d+)',log);assert not p.returncode and 'Error:' not in log and match and int(match.group(1))==expected;products[name]=dict(pages=expected,sha256=run.sha(pdf));print('MERGED_PDF_PASS',name,expected,flush=True)
run.save(R/'merge_commands.json',commands)
# Raster representative bins only; vector PDFs contain every requested interval.
previews=[]
for name,pages in [('fit_low_he4.pdf',[1,4,20,60,100]),('fit_low_he3.pdf',[22,41,81,121]),('fit_low_d.pdf',[41]),('fit_low_t.pdf',[41]),('fit_high_he4.pdf',[1,12]),('purity_vs_p_over_q.pdf',[1,2])]:
 for page in pages:
  image=R/(pathlib.Path(name).stem+f'_page_{page:03d}.png');cmd=['gs','-q','-dSAFER','-dNOPAUSE','-dBATCH','-dPDFSTOPONERROR','-dUseCropBox','-sDEVICE=png16m','-r100','-dFirstPage='+str(page),'-dLastPage='+str(page),'-sOutputFile='+str(image),str(R/name)];p=subprocess.run(cmd,capture_output=True,text=True);assert not p.returncode;previews.append(str(image))
run.save(R/'preview_images.json',previews)
run.save(R/'completed.json',dict(executionComplete=True,productionEnabled=False,response='gaussian',background='none',nativeFits=1579,aggregateBins=364,supportedBins=sum(r['quality'] for r in summary['rows']),meanFixedIntervals=9,fixedMeansExactlyEqualLinear=True,meanFixedErrorMeaning='zero conditional error; calibration uncertainty excluded',linearSeedsChecked=len(seedchecks),linearSeedsUnavailable=sum(not e['linearSeedAvailable'] for e in seedchecks),he3HalfHeightChecks=len(windows),rootReadbackZeroSmoothAndRawInvariant='GAUSSIAN_ONLY_ROOT_PASS' in (R/'ValidateGaussianOnly.log').read_text(),allPdfPagesReadback=True,uniqueFitPages=1579,protectedInputsUnchanged=True,protectedFiles=len(run.read(R/'protected.json')),sourceSha256={p.name:run.sha(p) for p in (R/'source').iterdir() if p.is_file() and p.suffix in ['.py','.h','.C','.yaml']},products=products,numericalProductsSha256={p.name:run.sha(p) for p in R.iterdir() if p.is_file() and p.suffix in ['.tsv','.yaml','.json'] and p.name!='completed.json'}))
print('GAUSSIAN_ONLY_COMPLETE',1579,'fits;',sum(r['quality'] for r in summary['rows']),'supported bins;',len(fixed),'fixed centers',flush=True)
