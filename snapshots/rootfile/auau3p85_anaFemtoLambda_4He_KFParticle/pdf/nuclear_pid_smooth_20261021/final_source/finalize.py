#!/usr/bin/env python3
"""Final numerical and provenance checks. Does not alter production files."""
import collections,csv,json,math,pathlib,shutil,subprocess
from study import P,OLD,OUT,SP,CONFIG,read,save,sha,sources
settings=read(P/'tools/nuclear_pid_smooth_20261021/settings.yaml');candidate=read(OUT/'candidate_summary.json');original=read(OUT/'original_summary.json');protected=read(OUT/'protected_input_audit.json')['protectedFiles']
for name,h in protected.items():
 if sha(name)!=h:raise ValueError('Protected input/config changed '+name)
checks=[];nativeCount=0
for band,label,src in sources():
 dst=OUT/'candidate'/(band+'_'+label);a=read(dst/'cutoff_fits.json')['slices'];b=read(src/'cutoff_fits.json')['slices']
 for x,y in zip(a,b):
  for key in ['fitLo','fitHi','columnWindowYields','shapePoint','pearson','dataDeviance','numericalValid']:
   if x[key]!=y[key]:raise ValueError('Unexpected change outside the PID gate '+str(dst))
 checks.append(dict(run=str(dst),nativeRows=len(a),parametersAndYieldsExactlyUnchanged=True));nativeCount+=len(a)
for folder in sorted((OUT/'candidate').glob('extension_*')):
 if not read(folder/'completed.json')['executionComplete']:raise ValueError('Incomplete extension')
 rows=read(folder/'cutoff_fits.json')['slices'];nativeCount+=len(rows)
 if len(rows)!=round((settings['auditRMax']-settings['extensionRMin'])/.01):raise ValueError('Incomplete extension range')
 for i,row in enumerate(rows):
  if abs(row['lo']-(settings['extensionRMin']+.01*i))>1e-8:raise ValueError('Extension coordinate mismatch')
# Every declared model must cover each bin, regardless of support.
for sp in SP:
 for band in ['low','high']:
  rows=[r for r in candidate['bins'] if r['species']==sp and r['band']==band]
  lo=settings['auditRMin'] if band=='low' else settings['tofThreshold'];hi=settings['tofThreshold'] if band=='low' else settings['auditRMax'];n=round((hi-lo)/settings['aggregateWidth'])
  if len(rows)!=n:raise ValueError('Final R coverage mismatch')
  for i,r in enumerate(sorted(rows,key=lambda r:r['rLo'])):
   if abs(r['rLo']-(lo+i*settings['aggregateWidth']))>1e-8 or r['modelCount']!=(20 if band=='low' else 4):raise ValueError('Missing R/model coverage')
   for model,v in r['models'].items():
    if v['totalSelected']>0 and abs(v['purity']-v['signalSelected']/v['totalSelected'])>1e-12:raise ValueError('Purity is not S/total')
    if v['supported'] and v['reasons']:raise ValueError('Inconsistent support flags')
    if v['supported'] and not 0<=v['purity']<=1:raise ValueError('Invalid supported purity')
# Preserve old-vs-new parameter identity and compare purity only where both were assessed.
historical=read(OLD/'summary_free_sensitivity_001/summary.json');oldindex={(r['species'],round(r['rLo'],6)):r for r in historical['aggregateRows']};differences={}
for s in SP:
 rows=[r for r in candidate['bins'] if r['species']==s and r['band']=='low' and r['quality'] and oldindex[s,round(r['rLo'],6)]['modelsSupported']]
 differences[s]=dict(commonAssessedBins=len(rows),maximumAbsoluteNominalPurityDifference=max((abs(r['purity']-oldindex[s,round(r['rLo'],6)]['purityNominal']) for r in rows),default=None))
# Explicit intervals are diagnostics, never automatic production veto/cutoff updates.
intervals={}
for s in SP:
 intervals[s]={}
 for band in ['low','high']:
  out=[]
  for r in sorted([r for r in candidate['bins'] if r['species']==s and r['band']==band],key=lambda r:r['rLo']):
   if not r['quality'] or r['purityLow']<settings['purityThreshold']:continue
   if out and abs(out[-1][1]-r['rLo'])<1e-8:out[-1][1]=r['rHi']
   else:out.append([r['rLo'],r['rHi']])
  intervals[s][band]=out
save(OUT/'study_results.yaml',dict(schemaVersion=1,productionEnabled=False,referenceGit=settings['referenceGit'],coefficients=read(OUT/'coefficients.json'),tofMass2=settings['tofMass2'],purityThreshold=settings['purityThreshold'],assessedIntervals=intervals,productionCutoffUpdate=None,highUpperLimitPolicy='Evaluation covers 2 to 5; supported intervals are not automatically statistical upper limits.'))
save(OUT/'assessed_ranges.json',dict(productionEnabled=False,policy='Supported intervals with model-low >=0.5. Unassessed gaps are not low purity and are never bridged. No production cutoff/veto changed.',intervals=intervals))
# Closed PDFs must have been verified, not merely generated.
for folder,stage in [(OUT,'summary'),(OUT/'atlas_low_001','low'),(OUT/'atlas_high_001','high')]:
 v=read(folder/('pdf_validation_'+stage+'.json'))
 for n,r in v['products'].items():
  if sha(folder/n)!=r['sha256']:raise ValueError('PDF changed after verification')
source=OUT/'final_source';source.mkdir()
for p in (P/'tools/nuclear_pid_smooth_20261021').iterdir():
 if p.suffix in ['.py','.C','.yaml','.md']:shutil.copy2(p,source/p.name)
logs=['render_summary.log','atlas_low_001/render.log','atlas_high_001/render.log']
for n in logs:
 if 'SMOOTH_PLOT_COMPLETE' not in (OUT/n).read_text():raise ValueError('Plot completion marker missing')
context=read(P/'rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001/context.json');fp=context['inputFingerprint'];root=pathlib.Path(fp['path']);st=root.stat()
for k,value in [('size',st.st_size),('mtime_ns',st.st_mtime_ns),('device',st.st_dev),('inode',st.st_ino)]:
 if fp[k]!=value:raise ValueError('Original ROOT changed')
productionCounts={}
for s in SP:
 productionCounts[s]={band:sum(r['quality'] for r in candidate['bins'] if r['species']==s and r['band']==band) for band in ['low','high']}
save(OUT/'completed.json',dict(executionComplete=True,productionEnabled=False,physicsScope='10/17 conditional response/background model estimates; not complete systematic validation',candidateFitSlices=nativeCount,originalNominalFitSlices=3400,reproductionFitSlices=24,reproductionScalarChecks=96,allHistoricalInputAndProductionConfigsPreserved=True,originalRootFingerprintPreserved=True,fitResponseInvariantChecks=checks,lowHistoricalComparison=differences,supportedBins=productionCounts,allRequiredPDFPagesReadback=1610,sourceSha256={p.name:sha(p) for p in source.iterdir()},productsSha256={n:sha(OUT/n) for n in ['coefficients.json','coefficients.yaml','candidate_summary.json','candidate_summary.tsv','original_summary.json','original_summary.tsv','assessed_ranges.json','study_results.yaml','curves.tsv','protected_input_audit.json','git_original_audit.json','purity_vs_p_over_q.pdf','dedx_vs_p_over_q_curves.pdf']}))
print('FINAL_STUDY_COMPLETE',nativeCount,'candidate slices',differences,flush=True)
