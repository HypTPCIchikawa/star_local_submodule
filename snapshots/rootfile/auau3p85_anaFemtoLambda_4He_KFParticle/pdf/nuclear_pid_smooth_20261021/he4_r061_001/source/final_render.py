#!/usr/bin/env python3
import json,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R=run.ROOT
C=run.cfg
(R/'detail_settings.tsv').write_text('\n'.join(k+'\t'+str(v) for k,v in [('minimumR',C['minimumRExclusive']),('detailMaxR',C['detailRMax']),('peakMaxR',C['peakPlotRMax']),('purityThreshold',C['purityThreshold'])])+'\n')
(R/'anchors.tsv').write_text('r\the3\the4\n'+'\n'.join(str(a['r'])+'\t'+str(a['he3'])+'\t'+str(a['he4']) for a in C['userPeakAnchors'])+'\n')
src=R/'plot_source';src.mkdir()
for name in ['PlotR061.C','PlotR061Purity.C','PlotR061Detail.C','render.py','analyze.py','settings.yaml']:shutil.copy2(HERE/name,src/name)
commands=[]
for name,invocation,marker in [('PlotR061',json.dumps(str(R))+',"",1','SMOOTH_PLOT_COMPLETE'),('PlotR061Purity',json.dumps(str(R)),'PURITY_AXES_COMPLETE'),('PlotR061Detail',json.dumps(str(R)),'R061_DETAIL_COMPLETE')]:
 cmd=['bash',str(run.P/'script/singularity_root_batch.sh'),str(run.P/run.cfg['runtimeMainconf']),'root4star','-l','-b','-q',str(src/(name+'.C'))+'+('+invocation+')'];commands.append(cmd)
 with (R/(name+'.log')).open('w') as f:status=subprocess.run(cmd,cwd=run.P,stdout=f,stderr=subprocess.STDOUT).returncode
 assert status==0 and marker in (R/(name+'.log')).read_text(),name
run.save(R/'plot_commands.json',commands)
products={}
for name,expected in [('fit_atlas_gaussian_he4.pdf',139),('fit_atlas_lognormal_he4.pdf',139),('purity_vs_p_over_q.pdf',2),('he4_low_r_detail.pdf',2)]:
 cmd=['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(R/name)];proc=subprocess.run(cmd,capture_output=True,text=True);text=proc.stdout+proc.stderr;(R/(name+'.readback.log')).write_text(text);match=re.search(r'Processing pages 1 through (\d+)',text);assert not proc.returncode and 'Error:' not in text and match and int(match.group(1))==expected,name;products[name]=dict(pages=expected,sha256=run.sha(R/name));print('PDF_PASS',name,expected,flush=True)
run.save(R/'plot_completed.json',dict(executionComplete=True,allPagesReadback=True,products=products,sourceSha256={p.name:run.sha(p) for p in src.iterdir() if p.suffix in ['.C','.py','.yaml']}))
