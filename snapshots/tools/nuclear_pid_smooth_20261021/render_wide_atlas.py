#!/usr/bin/env python3
"""Display-only expansion of saved atlas x ranges; never runs a fit."""
import argparse,csv,json,re,shutil,subprocess,math
from pathlib import Path
from study import P,OUT,read,save,sha
HERE=Path(__file__).resolve().parent
parser=argparse.ArgumentParser();parser.add_argument('band',choices=['low','high']);parser.add_argument('--preview',action='store_true');args=parser.parse_args()
cfg=read(HERE/'wide_atlas_settings.yaml');src=OUT/('atlas_'+args.band+'_'+cfg['sourceSuffix']);dst=OUT/('atlas_'+args.band+('_wide_preview_001' if args.preview else '_'+cfg['outputSuffix']))
if dst.exists():raise ValueError('Never overwrite '+str(dst))
original=read(OUT/'completed.json');protected={str(OUT/name):value for name,value in original['productsSha256'].items()}
protected[str(src/'candidate_atlas.tsv')]=sha(src/'candidate_atlas.tsv')
for pdf in src.glob('*.pdf'):protected[str(pdf)]=sha(pdf)
rows=list(csv.DictReader((src/'candidate_atlas.tsv').open(),delimiter='\t'));fields=list(rows[0])+['logLo','logHi'];axislo,axishi=cfg['axisRange'];clipped=0;fallbacks=[]
for row in rows:
 protected[row['root']]=sha(row['root']) if row['root'] not in protected else protected[row['root']]
 for low,high in [('linearLo','linearHi'),('logLo','logHi')]:
  a=float(row['linearLo' if low=='linearLo' else 'fitLo']);b=float(row['linearHi' if high=='linearHi' else 'fitHi'])
  if not math.isfinite(a) or not math.isfinite(b) or b<=a or a>=axishi or b<=axislo:
   if cfg['invalidPreviousDisplayPolicy']!='full_native_axis':raise ValueError('Invalid prior display')
   fallbacks.append(dict(species=row['species'],model=row['model'],index=row['index'],panel=low,priorLo=a,priorHi=b));a,b=axislo,axishi
  mid=(a+b)/2;half=(b-a)*cfg['widthScale']/2
  row[low]=format(max(axislo,mid-half),'.17g');row[high]=format(min(axishi,mid+half),'.17g')
  clipped+=int(mid-half<axislo or mid+half>axishi)
 if float(row['linearHi'])<=float(row['linearLo']) or float(row['logHi'])<=float(row['logLo']):raise ValueError('Invalid display range')
if args.preview:rows=[r for r in rows if r['species']=='he4' and int(r['index']) in [40,80,120]]
dst.mkdir();snapshot=dst/'source';snapshot.mkdir()
for file in [Path(__file__).resolve(),HERE/'PlotWideAtlas.C',HERE/'wide_atlas_settings.yaml']:shutil.copy2(file,snapshot/file.name)
with (dst/'candidate_atlas.tsv').open('w') as f:
 writer=csv.DictWriter(f,fieldnames=fields,delimiter='\t');writer.writeheader();writer.writerows(rows)
# Protect original bounds/fit metadata: only the display columns may differ.
old={(r['species'],r['model'],r['root'],r['index']):r for r in csv.DictReader((src/'candidate_atlas.tsv').open(),delimiter='\t')}
for r in rows:
 ref=old[r['species'],r['model'],r['root'],r['index']]
 for k,v in ref.items():
  if k not in ['linearLo','linearHi'] and r[k]!=v:raise ValueError('Changed fit metadata')
cmd=['bash',str(P/'script/singularity_root_batch.sh'),str(P/cfg['runtimeMainconf']),'root4star','-l','-b','-q',str(HERE/'PlotWideAtlas.C')+'+('+json.dumps(str(dst))+',"",1)']
save(dst/'command.json',cmd);save(dst/'context.json',dict(displayOnly=True,widthScale=cfg['widthScale'],clippedPanels=clipped,invalidPriorDisplayFallbacks=fallbacks,sourceManifest=str(src/'candidate_atlas.tsv'),protected=protected))
print('WIDE_RENDER_START',args.band,'preview',args.preview,flush=True)
with (dst/'render.log').open('w') as f:status=subprocess.run(cmd,cwd=P,stdout=f,stderr=subprocess.STDOUT).returncode
if status or 'SMOOTH_PLOT_COMPLETE' not in (dst/'render.log').read_text():raise RuntimeError('Rendering failed; inspect '+str(dst))
products={}
for sp in sorted({r['species'] for r in rows}):
 pdf=dst/('fit_atlas_'+args.band+'_'+sp+'.pdf');pages=len({r['rLo'] for r in rows if r['species']==sp});log=dst/(pdf.name+'.readback.log')
 with log.open('w') as f:status=subprocess.run(['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)],stdout=f,stderr=subprocess.STDOUT).returncode
 text=log.read_text();match=re.search(r'Processing pages 1 through (\d+)',text)
 if status or not match or int(match.group(1))!=pages or 'Error:' in text:raise RuntimeError('PDF validation failed '+str(pdf))
 products[pdf.name]=dict(pages=pages,sha256=sha(pdf));print('WIDE_PDF_PASS',sp,pages,flush=True)
for path,expected in protected.items():
 if sha(path)!=expected:raise ValueError('Saved numerical data or original PDF changed '+path)
save(dst/'completed.json',dict(executionComplete=True,displayOnly=True,productionEnabled=False,noRefit=True,numericalResultsAndOldPdfsUnchanged=True,widthScale=cfg['widthScale'],products=products,sourceSha256={p.name:sha(p) for p in snapshot.iterdir()}))
print('WIDE_ATLAS_COMPLETE',dst,flush=True)
