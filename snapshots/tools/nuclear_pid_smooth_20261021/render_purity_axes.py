#!/usr/bin/env python3
"""Change only first-page species x axes in the saved purity plot."""
import csv,json,re,shutil,subprocess
from pathlib import Path
from study import P,OUT,read,save,sha
HERE=Path(__file__).resolve().parent;cfg=read(HERE/'purity_axis_settings.yaml');dst=OUT/cfg['outputDirectory']
if dst.exists():raise ValueError('Do not overwrite '+str(dst))
base=read(OUT/'completed.json');protected={str(OUT/n):h for n,h in base['productsSha256'].items()};names=['candidate_summary.tsv','historical_summary.tsv','original_summary.tsv'];rows=[]
for name in names[:2]:rows+=list(csv.DictReader((OUT/name).open(),delimiter='\t'))
axes=[]
for sp in ['d','t','he3','he4']:
 kept=[r for r in rows if r['species']==sp and r['quality']=='1'];pointlo=min(float(r['rLo']) for r in kept);pointhi=max(float(r['rHi']) for r in kept)
 lo,hi=cfg['fixedDRange'] if sp=='d' else [pointlo-cfg['paddingGeV'],pointhi+cfg['paddingGeV']]
 assert all(lo<=float(r['rLo'])+1e-9 and hi>=float(r['rHi'])-1e-9 for r in kept)
 axes.append(dict(species=sp,lo=lo,hi=hi,plottedBinEdges=[pointlo,pointhi]))
dst.mkdir();snapshot=dst/'source';snapshot.mkdir()
for p in [HERE/'PlotPurityAxes.C',HERE/'purity_axis_settings.yaml',Path(__file__).resolve()]:shutil.copy2(p,snapshot/p.name)
for name in names:shutil.copy2(OUT/name,dst/name);protected[str(OUT/name)]=sha(OUT/name)
(dst/'axes.tsv').write_text('species\tlo\thi\n'+'\n'.join(a['species']+'\t'+format(a['lo'],'.17g')+'\t'+format(a['hi'],'.17g') for a in axes)+'\n')
cmd=['bash',str(P/'script/singularity_root_batch.sh'),str(P/cfg['runtimeMainconf']),'root4star','-l','-b','-q',str(snapshot/'PlotPurityAxes.C')+'+('+json.dumps(str(dst))+')'];save(dst/'command.json',cmd)
with (dst/'render.log').open('w') as f:status=subprocess.run(cmd,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
if status or 'PURITY_AXES_COMPLETE' not in (dst/'render.log').read_text():raise RuntimeError('Rendering failed')
pdf=dst/'purity_vs_p_over_q.pdf'
with (dst/'pdf_readback.log').open('w') as f:status=subprocess.run(['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)],stdout=f,stderr=subprocess.STDOUT).returncode
text=(dst/'pdf_readback.log').read_text()
if status or not re.search(r'Processing pages 1 through 2\b',text) or 'Error:' in text:raise RuntimeError('PDF validation failed')
# Compare the entire second page rendered with exactly the same Ghostscript options.
for stem,path,page in [('page1',pdf,1),('page2',pdf,2),('original_page2',OUT/'purity_vs_p_over_q.pdf',2)]:
 subprocess.run(['gs','-q','-dSAFER','-dBATCH','-dNOPAUSE','-sDEVICE=png16m','-r90','-dFirstPage='+str(page),'-dLastPage='+str(page),'-sOutputFile='+str(dst/(stem+'.png')),str(path)],check=True)
if sha(dst/'page2.png')!=sha(dst/'original_page2.png'):raise ValueError('Page 2 changed')
for name,h in protected.items():
 if sha(name)!=h:raise ValueError('Numerical result or original product changed '+name)
save(dst/'completed.json',dict(executionComplete=True,displayOnly=True,noRefit=True,numericalResultsAndOriginalProductsUnchanged=True,page2PixelIdentical=True,pages=2,axes=axes,products={n:sha(dst/n) for n in ['purity_vs_p_over_q.pdf','page1.png','page2.png','axes.tsv']},sourceSha256={p.name:sha(p) for p in snapshot.iterdir() if p.suffix in ['.C','.py','.yaml']}))
print('PURITY_AXES_COMPLETE',axes,flush=True)
