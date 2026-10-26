#!/usr/bin/env python3
import json,pathlib,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent
sys.path.insert(0,str(HERE));import run
root=run.ROOT
src=root/'plot_source';src.mkdir()
for name in ['PlotNoSmooth.C','render.py','analyze.py']:shutil.copy2(HERE/name,src/name)
shutil.copy2(run.CONFIG_PATH,src/run.CONFIG_PATH.name)
command=['bash',str(run.P/'script/singularity_root_batch.sh'),str(run.P/run.cfg['runtimeMainconf']),'root4star','-l','-b','-q',str(src/'PlotNoSmooth.C')+'+('+json.dumps(str(root))+',"",1)']
run.save(root/'plot_command.json',command)
with (root/'plot.log').open('w') as f:status=subprocess.run(command,cwd=run.P,stdout=f,stderr=subprocess.STDOUT).returncode
log=(root/'plot.log').read_text()
if status or 'SMOOTH_PLOT_COMPLETE' not in log or log.count('pages=160')!=2:raise RuntimeError('Plot failed, see '+str(root/'plot.log'))
for response in run.cfg['responses']:
 pdf=root/('fit_atlas_'+response+'_he4.pdf')
 check=subprocess.run(['gs','-q','-dNOPAUSE','-dBATCH','-sDEVICE=nullpage',str(pdf)],capture_output=True,text=True)
 (root/(response+'_pdf_readback.log')).write_text(check.stdout+check.stderr)
 assert check.returncode==0
 info=subprocess.check_output(['pdfinfo',str(pdf)],text=True);(root/(response+'_pdfinfo.txt')).write_text(info)
 assert 'Pages:           160' in info
run.save(root/'plot_completed.json',dict(readbackPass=True,pages=320,productionEnabled=False,products={p.name:run.sha(p) for p in root.glob('*.pdf')}))
print('NO_SMOOTH_PLOT_COMPLETE 320 pages')
