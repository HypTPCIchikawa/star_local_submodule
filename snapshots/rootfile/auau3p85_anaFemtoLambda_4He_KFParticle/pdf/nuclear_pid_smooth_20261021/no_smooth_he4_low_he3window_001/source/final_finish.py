#!/usr/bin/env python3
import collections,json,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(HERE));import run
r=run.ROOT;o=run.OUT
protected=run.read(r/'protected.json')
for name,entry in run.read(o/'wide_atlas_completed.json')['books'].items():protected[str(o/name)]=entry['sha256']
for name,h in run.read(o/'purity_axis_002/completed.json')['products'].items():protected[str(o/'purity_axis_002'/name)]=h
for name,h in protected.items():assert run.sha(pathlib.Path(name))==h,name
windows=run.read(r/'fit_windows.json');assert len(windows)==160
for w in windows:
 assert w['fitLo']<=w['oldFitLo'] and all(w['fitLo']<=v['leftHalfMaximum'] for v in w['references'])
rows=0
for response in run.cfg['responses']:
 for stage,bg in [('none','none'),('baseline','linear'),('baseline','exponential'),('warm','linear'),('warm','exponential')]:
  d=r/(stage+'_'+response+'_'+bg);assert run.read(d/'completed.json')['rows']==160
  ref=o/('candidate/low_he4_'+response+'_'+('linear' if bg=='none' else bg)+'_nominal')
  assert run.sha(d/'gates.tsv')==run.sha(ref/'gates.tsv')
  for a,w in zip(run.read(d/'cutoff_fits.json')['slices'],windows):assert a['index']==w['index'] and a['fitLo']==w['fitLo'] and a['fitHi']==w['fitHi']
  rows+=160
assert 'NO_SMOOTH_ROOT_VALIDATION_PASS' in (r/'root_validation.log').read_text()
assert run.read(r/'plot_completed.json')['readbackPass']
# A compact vector PDF extracted from the complete atlases, retaining their annotations.
extract=r/'representative_pages';extract.mkdir(exist_ok=True);pages=[];commands=[]
def captured(cmd,name):
 proc=subprocess.run(cmd,capture_output=True,text=True);(r/name).write_text(proc.stdout+proc.stderr);assert proc.returncode==0

for response in run.cfg['responses']:
 for index in run.cfg['controlIndices']:
  dst=extract/(response+'_'+str(index)+'.pdf')
  cmd=['gs','-q','-dNOPAUSE','-dBATCH','-sDEVICE=pdfwrite','-dFirstPage='+str(index+1),'-dLastPage='+str(index+1),'-sOutputFile='+str(dst),str(r/('fit_atlas_'+response+'_he4.pdf'))]
  captured(cmd,'extract_'+response+'_'+str(index)+'.log');commands.append(cmd);pages.append(str(dst))
cmd=['gs','-q','-dNOPAUSE','-dBATCH','-sDEVICE=pdfwrite','-sOutputFile='+str(r/'representative_fits.pdf')]+pages
captured(cmd,'representative_transform_'+str(len(commands))+'.log');commands.append(cmd)
cmd=['gs','-q','-dNOPAUSE','-dBATCH','-sDEVICE=png16m','-r85','-dFirstPage=1','-dLastPage=1','-sOutputFile='+str(r/'representative_gaussian.png'),str(r/'representative_fits.pdf')]
captured(cmd,'representative_transform_'+str(len(commands))+'.log');commands.append(cmd)
run.save(r/'representative_commands.json',commands)
check=subprocess.run(['gs','-dSAFER','-dNOPAUSE','-dBATCH','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(r/'representative_fits.pdf')],capture_output=True,text=True);assert check.returncode==0 and 'Error:' not in check.stdout+check.stderr
match=re.search(r'Processing pages 1 through (\d+)',check.stdout+check.stderr);assert match and int(match.group(1))==6
(r/'representative_readback.log').write_text(check.stdout+check.stderr)
for name in ['run.py','ValidateNoSmooth.C','finish.py','render.py']:shutil.copy2(HERE/name,r/'source'/('final_'+name))
run.save(r/'completed.json',dict(executionComplete=True,productionEnabled=False,realFits=rows,fullAtlasPages=320,representativePages=6,allPdfReadbackPass=True,oldProductsUnchanged=len(protected),allWindowsIncludeHe3LeftHalfMaximum=True,extendedWindows=sum(w['fitLo']<w['oldFitLo'] for w in windows),sameWindowForAllModels=True,zeroSmoothAll320=True,identicalRawSpectra1280Controls=True,productsSha256={p.name:run.sha(p) for p in r.iterdir() if p.is_file() and p.name!='completed.json'}))
print('NO_SMOOTH_FINAL_VALIDATION_PASS',rows,'fits;',len(protected),'protected files')
