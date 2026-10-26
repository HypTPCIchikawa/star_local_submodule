#!/usr/bin/env python3
import collections,json,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,C=run.R,run.C;build=run.read(R/'report_build.json');physics=run.read(R/'physics_validation.json');plots=run.read(R/'pdf_validation.json');assert physics['executionComplete'] and plots['allPagesReadback'] and build['missingGlyphs']==False and build['overfullBoxes']==0
for n,v in plots['products'].items():assert run.sha(R/n)==v['sha256']
for stage in ['baseline','candidate']:
 for region in C['regions']:
  f=R/(stage+'_'+region['name']);done=run.read(f/'completed.json')
  for n,h in done['products'].items():assert run.sha(f/n)==h
for p,h in run.read(R/'protected.json').items():assert run.sha(pathlib.Path(p))==h
header=(R/'candidate_pid/NuclearIdDeDxVsMom.h').read_text();coefs=run.read(R/'coefficients.json');original=(R/'git_NuclearIdDeDxVsMom.h').read_text()
for k,array in [('mean','kMean'),('sigma','kSigma'),('domain','kPRange')]:
 block=header.split(array+'[kNSpecies] = {',1)[1].split('};',1)[0];parsed=[[float(x) for x in row.split(',')] for row in re.findall(r'\{([^{}]+)\}',block)];assert parsed==coefs['candidate'][k];assert parsed[:3]==coefs['original'][k][:3]
cache=json.loads(run.read(R/'candidate_low_linear/command.json')[-1].split('+(',1)[1].split(',',1)[0]);contextPath=pathlib.Path(cache).parent/'context.json';context=run.read(contextPath);inputfile=pathlib.Path(context['inputFingerprint']['path']);st=inputfile.stat();assert st.st_size==context['inputFingerprint']['size'] and st.st_mtime_ns==context['inputFingerprint']['mtime_ns'];assert '2021P24iy_partial_merge' in inputfile.name
run.save(R/'input_provenance.json',dict(contextPath=str(contextPath),contextSha256=run.sha(contextPath),cachePath=cache,cacheSha256=run.sha(pathlib.Path(cache)),inputFingerprint=context['inputFingerprint'],fingerprintVerified=True,statistics='existing partial 2021 P24iy merged sample; not full-statistics assertion',catalogSettings=str(run.P/'config/analysis/analysis_info_auau3p85_anaFemtoLambda_4He_KFParticle.yaml')))
pdf=R/(build['stem']+'.pdf');cmd=['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)];p=subprocess.run(cmd,capture_output=True,text=True);log=p.stdout+p.stderr;(R/'report_readback.log').write_text(log);m=re.search(r'Processing pages 1 through (\d+)',log);assert p.returncode==0 and 'Error:' not in log and m and int(m.group(1))==build['expectedPages']
cmd=['gs','-q','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-dUseCropBox','-sDEVICE=png16m','-r90','-sOutputFile='+str(R/'report_page_%02d.png'),str(pdf)];subprocess.run(cmd,check=True,capture_output=True)
for p in HERE.iterdir():
 if p.is_file() and p.suffix in ['.py','.C','.yaml','.md']:shutil.copy2(p,R/'source'/p.name)
products=dict(plots['products']);products[pdf.name]=dict(pages=build['expectedPages'],sha256=run.sha(pdf));summary=run.read(R/'summary.json')
run.save(R/'completed.json',dict(executionComplete=True,productionEnabled=False,target='he4_only',referenceGit=C['referenceGit'],d_t_he3_pid_exact_git=True,onlyHe4PidCoefficientsChanged=True,baselineFits=199,candidateFits=199,aggregateBins=88,supportedBins=sum(r['quality'] for r in summary['rows']),regionPolicies=C['regions'],allPdfPagesReadback=True,reportPages=build['expectedPages'],totalPdfPages=sum(v['pages'] for v in products.values()),noMissingGlyphs=True,noOverfullBoxes=True,physicsValidation=physics,inputProvenance=run.read(R/'input_provenance.json'),oldOutputsUnchanged=True,products=products,sourceSha256={p.name:run.sha(p) for p in (R/'source').iterdir() if p.is_file() and p.suffix in ['.py','.C','.h','.yaml']},numericProductsSha256={str(p.relative_to(R)):run.sha(p) for p in R.iterdir() if p.is_file() and p.suffix in ['.json','.tsv','.yaml'] and p.name!='completed.json'},candidateHeaderSha256=run.sha(R/'candidate_pid/NuclearIdDeDxVsMom.h')))
print('HE4_STUDY_COMPLETE',len(summary['rows']),'bins;',sum(r['quality'] for r in summary['rows']),'supported;',build['expectedPages'],'report pages;',sum(v['pages'] for v in products.values()),'PDF pages',flush=True)
