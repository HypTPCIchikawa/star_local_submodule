#!/usr/bin/env python3
import collections,json,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,C=run.R,run.C
for n in ['PlotGaussianOnly.C','ValidateGaussianOnly.C']:shutil.copy2(HERE/n,R/'source'/n)
rows=['new\told\tfirst\tlast']
for j in run.jobs():
 folder=R/j['name'];indices=[r['index'] for r in run.read(folder/'cutoff_fits.json')['slices']];assert indices==list(range(min(indices),max(indices)+1));rows.append(str(folder/'cutoff_fits.root')+'\t'+str(pathlib.Path(j['reference'])/'cutoff_fits.root')+'\t'+str(min(indices))+'\t'+str(max(indices)))
(R/'validation_roots.tsv').write_text('\n'.join(rows)+'\n')
commands=[]
for name,args,marker in [('ValidateGaussianOnly',[str(R)],'GAUSSIAN_ONLY_ROOT_PASS'),('PlotGaussianOnly',[str(R),''],'LINEAR_REPORT_PLOTS_COMPLETE')]:
 cmd=['bash',str(run.P/'script/singularity_root_batch.sh'),str(run.P/C['runtimeMainconf']),'root4star','-l','-b','-q',str(R/'source'/(name+'.C'))+'+('+','.join(json.dumps(x) for x in args)+')'];commands.append(cmd)
 with (R/(name+'.log')).open('w') as f:status=subprocess.run(cmd,cwd=run.P,stdout=f,stderr=subprocess.STDOUT).returncode
 assert not status and marker in (R/(name+'.log')).read_text(),name
run.save(R/'plot_commands.json',commands)
counts=collections.Counter(e['group'] for e in run.read(R/'examples.json'));counts={('fit_'+k+'.pdf'):v for k,v in counts.items()};counts['purity_vs_p_over_q.pdf']=2;products={}
for name,expected in sorted(counts.items()):
 pdf=R/name;cmd=['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)];p=subprocess.run(cmd,capture_output=True,text=True);log=p.stdout+p.stderr;(R/(name+'.readback.log')).write_text(log);match=re.search(r'Processing pages 1 through (\d+)',log);assert not p.returncode and 'Error:' not in log and match and int(match.group(1))==expected,name;products[name]=dict(pages=expected,sha256=run.sha(pdf));print('PDF_PASS',name,expected,flush=True)
run.save(R/'pdf_validation.json',dict(allPagesReadback=True,products=products))
