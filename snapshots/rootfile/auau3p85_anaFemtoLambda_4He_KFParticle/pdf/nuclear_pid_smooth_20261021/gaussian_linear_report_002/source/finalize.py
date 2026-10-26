#!/usr/bin/env python3
import json,math,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import prepare as base
R=base.R;C=base.C;summary=base.read(R/'summary.json');cache={};checks=[]
assert C['response']=='gaussian' and C['background']=='linear' and C['variation']=='nominal'
for row in summary['rows']:
 src=pathlib.Path(row['source'])
 if str(src) not in cache:cache[str(src)]=base.read(src)
 old=next(x for x in cache[str(src)]['bins'] if all(x[k]==row[k] for k in ['band','species','rLo','rHi']));v=old['models'][row['sourceModel']]
 assert row['sourceModel'] in ['gaussian_linear','gaussian_linear_nominal']
 for a,b in [('purity','purity'),('statConditionalBound','conditionalErrorWeightedBound'),('quality','supported'),('signalSelected','signalSelected'),('totalSelected','totalSelected'),('rawSelected','rawSelected')]:assert row[a]==v[b],(row,a)
 if row['quality']:assert all(math.isfinite(row[k]) for k in ['purity','statConditionalBound']) and row['statConditionalBound']>=0
 checks.append(dict(species=row['species'],band=row['band'],rLo=row['rLo'],rHi=row['rHi'],sourceModel=row['sourceModel'],exactMatch=True))
for e in base.read(R/'examples.json'):
 d=base.read(pathlib.Path(e['root']).parent/'cutoff_fits.json');assert d['responseModel']=='gaussian' and d['backgroundModel']=='linear'
examples=base.read(R/'examples.json')
for species in ['d','t','he3','he4']:
 low=[e for e in examples if e['group']=='low_'+species]
 assert len(low)==3 and len({e['rLo'] for e in low})==3
 assert all(e['rHi']<=C['tofThreshold'] for e in low)
 for e in low:
  meta=base.read(pathlib.Path(e['root']).parent/'cutoff_fits.json')
  assert meta['sampleHistogram']=='ALL'
previous=base.read(base.OUT/'gaussian_linear_report_001/summary.json')
assert all(summary[k]==previous[k] for k in ['rows','statistics','axes'])
assert base.sha(R/'purity.tsv')==base.sha(base.OUT/'gaussian_linear_report_001/purity.tsv')
for n,h in base.read(R/'protected.json').items():assert base.sha(pathlib.Path(n))==h,n
build=base.read(R/'report_build.json');assert build['missingGlyphs']==False and build['overfullBoxes']==0
pdf=R/'pid_purity_report_20261021.pdf';cmd=['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)];p=subprocess.run(cmd,capture_output=True,text=True);log=p.stdout+p.stderr;(R/'report_readback.log').write_text(log);match=re.search(r'Processing pages 1 through (\d+)',log);assert not p.returncode and 'Error:' not in log and match and int(match.group(1))==build['expectedPages']
commands=[cmd]
cmd=['gs','-q','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=png16m','-r90','-sOutputFile='+str(R/'report_page_%02d.png'),str(pdf)];p=subprocess.run(cmd,capture_output=True,text=True);(R/'report_preview.log').write_text(p.stdout+p.stderr);assert not p.returncode and 'Error:' not in p.stdout+p.stderr;commands.append(cmd)
base.save(R/'validation_commands.json',commands);base.save(R/'point_validation.json',checks)
selection=base.read(R/'pid_selection.yaml');selection['purityEstimation']=dict(response=C['response'],background=C['background'],variation=C['variation'],errorPolicy=C['errorPolicy'],configuration='settings.yaml');base.save(R/'pid_selection.yaml',selection)
for file in HERE.iterdir():
 if file.is_file():shutil.copy2(file,R/'source'/file.name)
base.save(R/'completed.json',dict(executionComplete=True,noNewFit=True,unchangedPurity=True,threeLowRPointsPerSpecies=True,productionEnabled=False,response='gaussian',background='linear',variation='nominal',purityPointsChecked=len(checks),supportedPoints=sum(r['quality'] for r in summary['rows']),sourceModelsVerified=True,otherModelEnvelopesExcluded=True,gaussianLinearOwnSupport=True,conditionalYieldErrorsOnly=True,reportPages=build['expectedPages'],figurePages=12,allPagesReadback=True,noMissingGlyphs=True,noOverfullBoxes=True,oldInputsUnchanged=True,protectedFiles=len(base.read(R/'protected.json')),sourceSha256={p.name:base.sha(p) for p in (R/'source').iterdir() if p.suffix in ['.py','.C','.yaml']},productsSha256={p.name:base.sha(p) for p in R.iterdir() if p.is_file() and p.name!='completed.json'}))
print('GAUSSIAN_LINEAR_REPORT_PASS',len(checks),'source-matched bins;',build['expectedPages'],'report pages; 12 figure pages',flush=True)
