#!/usr/bin/env python3
import json,math,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,O,C=run.ROOT,run.OUT,run.cfg
protected=run.read(R/'protected.json')
for name,x in run.read(O/'wide_atlas_completed.json')['books'].items():protected[str(O/name)]=x['sha256']
for name,h in run.read(O/'purity_axis_002/completed.json')['products'].items():protected[str(O/'purity_axis_002'/name)]=h
for name,h in protected.items():assert run.sha(pathlib.Path(name))==h,name
for name in ['NuclearFitCommon.h','NuclearLognormalFitCommon.h']:assert run.sha(R/'source'/name)==run.sha(run.PREVIOUS/'source'/name)
windows={w['index']:w for w in run.read(R/'fit_windows.json')};assert len(windows)==139
regression=[];fittedhalves=[];fitcount=0
for response in C['responses']:
 for stage,bg,oldstage in [('full','none','none'),('full','linear','baseline'),('full','exponential','baseline'),('warm','linear','warm'),('warm','exponential','warm')]:
  d=R/(stage+'_'+response+'_'+bg);rows=run.read(d/'cutoff_fits.json')['slices'];previous={r['index']:r for r in run.read(run.PREVIOUS/(oldstage+'_'+response+'_'+bg+'/cutoff_fits.json'))['slices']};shapes=run.table(d/'shapes.tsv');byindex={i:[s for s in shapes if int(s['index'])==i and s['species']=='he3'] for i in windows};assert len(rows)==139 and {r['index'] for r in rows}==set(windows)
  for r in rows:
   w=windows[r['index']];assert r['r']>C['minimumRExclusive'] and r['fitLo']==w['fitLo'] and r['fitHi']==w['fitHi']
   c=r['components']['he3'];ratio=c['sigma']/c['initialSigma'];assert math.isfinite(ratio) and ratio>0
   center=next(s for s in byindex[r['index']] if s['node']=='1');assert abs(c['mean']-float(center['mean']))<1e-8
   expected=(3 if r['r']<C['lowPeakAdjustmentMaxR'] else 2)+(1 if bg=='exponential' else 0)
   assert len(r['shapePoint'])==expected,(d,r['index'],r['shapePoint'])
   half=[]
   for s in byindex[r['index']]:
    m,sd=float(s['mean']),float(s['sigma'])*ratio;left=m-math.sqrt(2*math.log(2))*sd if response=='gaussian' else m*math.exp(-math.sqrt(2*math.log(2))*sd/m)
    assert r['fitLo']<=left+1e-8,(d,r['index'],r['fitLo'],left);half.append(left)
   fittedhalves.append(dict(folder=d.name,index=r['index'],fitLo=r['fitLo'],minimumFittedHe3LeftHalfMaximum=min(half)))
   if r['r']>=C['lowPeakAdjustmentMaxR']:
    prior=previous[r['index']];largest=0
    for k in ['pearson','dataDeviance','conditionalNdf','shapePoint','columnWindowYields']:
     a,b=r[k],prior[k];a=a if isinstance(a,list) else [a];b=b if isinstance(b,list) else [b];assert len(a)==len(b)
     for aa,bb in zip(a,b):
      delta=abs(aa-bb)/max(1.,abs(bb));largest=max(largest,delta);assert delta<=1e-8,(d,r['index'],k,aa,bb)
    for s in ['d','t','he3','he4']:
     for k in ['mean','sigma']:
      a,b=r['components'][s][k],prior['components'][s][k];assert abs(a-b)<=1e-8*max(1.,abs(b)),(d,r['index'],s,k)
    for k in ['signalSelected','totalSelected','rawGateCount','purity']:
     a,b=r['components']['he4']['selection'][k],prior['components']['he4']['selection'][k];assert (a==b) if a is None or b is None else abs(a-b)<=1e-8*max(1.,abs(b))
    regression.append(dict(folder=d.name,index=r['index'],maximumRelativeDifference=largest))
  fitcount+=len(rows)
summary=run.read(R/'summary.json');prior=run.read(O/'candidate_summary.json');unchanged=[r for r in prior['bins'] if not(r['species']=='he4' and r['band']=='low')];assert [r for r in summary['bins'] if not(r['species']=='he4' and r['band']=='low')]==unchanged
assert not any(r['nestedViolation'] for r in summary['backgroundComparisons'])
assert 'R061_ROOT_VALIDATION_PASS' in (R/'root_validation.log').read_text()
# Restore the original preparation snapshot, while retaining the final reproducible runner/config separately.
shutil.copy2(R/'source/diagnostic_free_mean_settings.yaml',R/'source/settings.yaml')
for n,h in run.read(R/'source_manifest.json')['sourceSha256'].items():assert run.sha(R/'source'/n)==h,n
for name in ['run.py','analyze.py','finish.py','settings.yaml','ValidateR061.C','PlotR061.C','PlotR061Purity.C','PlotR061Detail.C','render.py']:shutil.copy2(HERE/name,R/'source'/('final_'+name))
run.save(R/'regression_checks.json',dict(unchangedAboveR=C['lowPeakAdjustmentMaxR'],comparisons=regression,otherSpeciesAndHighHe4Exact=True,halfMaximumChecks=fittedhalves))
if len(sys.argv)>1 and sys.argv[1]=='numerical':
 print('NUMERICAL_PASS',len(regression),'comparisons; maximum relative difference',max(r['maximumRelativeDifference'] for r in regression));sys.exit(0)
def command(cmd,log):
 p=subprocess.run(cmd,capture_output=True,text=True);(R/log).write_text(p.stdout+p.stderr);assert not p.returncode,(cmd,log);return p.stdout+p.stderr
products={};commands=[]
# Full PDFs, including the finalized purity redraw and the low-R detail.
for name,pages in [('fit_atlas_gaussian_he4.pdf',139),('fit_atlas_lognormal_he4.pdf',139),('purity_vs_p_over_q.pdf',2),('he4_low_r_detail.pdf',2)]:
 cmd=['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(R/name)];text=command(cmd,name+'.final_readback.log');match=re.search(r'Processing pages 1 through (\d+)',text);assert 'Error:' not in text and match and int(match.group(1))==pages,name;products[name]=dict(pages=pages,sha256=run.sha(R/name));commands.append(cmd)
extract=R/'representative_pages';extract.mkdir(exist_ok=True);files=[]
indices=sorted(windows)
for response in C['responses']:
 for index in C['representativeIndices']:
  page=indices.index(index)+1;dst=extract/(response+'_'+str(index)+'.pdf');cmd=['gs','-q','-dNOPAUSE','-dBATCH','-sDEVICE=pdfwrite','-dFirstPage='+str(page),'-dLastPage='+str(page),'-sOutputFile='+str(dst),str(R/('fit_atlas_'+response+'_he4.pdf'))];command(cmd,'extract_'+response+'_'+str(index)+'.log');commands.append(cmd);files.append(str(dst))
cmd=['gs','-q','-dNOPAUSE','-dBATCH','-sDEVICE=pdfwrite','-sOutputFile='+str(R/'representative_fits.pdf')]+files;command(cmd,'representative_merge.log');commands.append(cmd)
cmd=['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(R/'representative_fits.pdf')];text=command(cmd,'representative_readback.log');match=re.search(r'Processing pages 1 through (\d+)',text);assert 'Error:' not in text and match and int(match.group(1))==len(files);products['representative_fits.pdf']=dict(pages=len(files),sha256=run.sha(R/'representative_fits.pdf'));commands.append(cmd)
for pdf,page,name in [('representative_fits.pdf',1,'fit_low_r.png'),('purity_vs_p_over_q.pdf',1,'purity_final.png'),('he4_low_r_detail.pdf',1,'purity_low_r.png'),('he4_low_r_detail.pdf',2,'peaks_low_r.png')]:
 cmd=['gs','-q','-dNOPAUSE','-dBATCH','-sDEVICE=png16m','-r100','-dFirstPage='+str(page),'-dLastPage='+str(page),'-sOutputFile='+str(R/name),str(R/pdf)];command(cmd,name+'.render.log');commands.append(cmd)
run.save(R/'validation_and_extract_commands.json',commands)
run.save(R/'completed.json',dict(executionComplete=True,productionEnabled=False,mainFits=fitcount,pilotFits=48,diagnosticModelWindowFits=695,protectedFilesUnchanged=len(protected),allFittedHe3HalfMaximaCovered=True,sameFitRangeAcrossAllModels=True,backgroundNestingViolations=0,numericalComparisonsUnchangedAboveR=C['lowPeakAdjustmentMaxR'],unchangedComparisons=len(regression),maximumRelativeDifference=max(r['maximumRelativeDifference'] for r in regression),otherSpeciesAndHighHe4Exact=True,allPdfPagesReadback=True,pdfProducts=products,configuration=C,sourceSha256={p.name:run.sha(p) for p in (R/'source').glob('final_*')},productsSha256={p.name:run.sha(p) for p in R.iterdir() if p.is_file() and p.name!='completed.json'}))
print('R061_FINAL_PASS',fitcount,'fits;',len(regression),'unchanged higher-R comparisons;',len(protected),'protected files;',sum(v['pages'] for v in products.values()),'PDF pages',flush=True)
