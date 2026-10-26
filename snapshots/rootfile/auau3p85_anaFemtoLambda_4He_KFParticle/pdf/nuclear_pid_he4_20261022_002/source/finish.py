#!/usr/bin/env python3
import json,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,C=run.R,run.C;build=run.read(R/'report_build.json');physics=run.read(R/'physics_validation.json');plots=run.read(R/'pdf_validation.json');assert physics['executionComplete'] and physics['allMeansExactlyPreviousLinear'] and physics['unchangedAggregateBins']==84;assert build['missingGlyphs']==False and build['overfullBoxes']==0
for n,v in plots['products'].items():assert run.sha(R/n)==v['sha256']
for region in C['regions']:
 folder=R/('candidate_'+region['name']);done=run.read(folder/'completed.json')
 for n,h in done['products'].items():assert run.sha(folder/n)==h
for p,h in run.read(R/'protected.json').items():assert run.sha(pathlib.Path(p))==h,p
for n in ['coefficients.json','curves.tsv','git_baseline.json']:assert run.sha(R/n)==run.sha(run.PREVIOUS/n)
priorInput=run.read(R/'input_provenance.json');assert run.sha(pathlib.Path(priorInput['cachePath']))==priorInput['cacheSha256'];meta=priorInput['inputFingerprint'];st=pathlib.Path(meta['path']).stat();assert st.st_size==meta['size'] and st.st_mtime_ns==meta['mtime_ns']
header=run.read(R/'coefficients.json');baseline=run.read(R/'git_baseline.json')['tables'];assert all(header['candidate'][k][:3]==baseline[k][:3] for k in ['mean','sigma','domain'])
pdf=R/(build['stem']+'.pdf');cmd=['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)];p=subprocess.run(cmd,capture_output=True,text=True);log=p.stdout+p.stderr;(R/'report_readback.log').write_text(log);m=re.search(r'Processing pages 1 through (\d+)',log);assert not p.returncode and 'Error:' not in log and m and int(m.group(1))==build['expectedPages']
cmd=['gs','-q','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-dUseCropBox','-sDEVICE=png16m','-r90','-sOutputFile='+str(R/'report_page_%02d.png'),str(pdf)];subprocess.run(cmd,check=True,capture_output=True)
for p in HERE.iterdir():
 if p.is_file() and p.suffix in ['.py','.C','.yaml','.md']:shutil.copy2(p,R/'source'/p.name)
products=dict(plots['products']);products[pdf.name]=dict(pages=build['expectedPages'],sha256=run.sha(pdf));summary=run.read(R/'summary.json')
run.save(R/'completed.json',dict(executionComplete=True,productionEnabled=False,revision='low-R fixed mean and Gaussian-only',previousOutput=str(run.PREVIOUS),newFits=19,reusedNativeFits=180,meanFixedIntervals=19,allFixedMeansExactlyPreviousLinear=True,pidCoefficientsUnchanged=True,unchangedAggregateBins=84,totalNativeIntervals=199,aggregateBins=88,supportedBins=sum(a['quality'] for a in summary['rows']),reportPages=build['expectedPages'],totalPdfPages=sum(v['pages'] for v in products.values()),allPdfPagesReadback=True,noMissingGlyphs=True,noOverfullBoxes=True,physicsValidation=physics,oldOutputsUnchanged=True,products=products,sourceSha256={p.name:run.sha(p) for p in (R/'source').iterdir() if p.is_file() and p.suffix in ['.py','.C','.h','.yaml']},numericProductsSha256={p.name:run.sha(p) for p in R.iterdir() if p.is_file() and p.suffix in ['.json','.tsv','.yaml'] and p.name!='completed.json'}))
print('FIXED_MEAN_GAUSSIAN_ONLY_COMPLETE',19,'new fits;',sum(a['quality'] for a in summary['rows']),'supported;',build['expectedPages'],'report pages;',sum(v['pages'] for v in products.values()),'total PDF pages',flush=True)
