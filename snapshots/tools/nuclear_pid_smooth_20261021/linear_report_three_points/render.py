#!/usr/bin/env python3
import json,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import prepare as base
R=base.R
shutil.copy2(HERE/'PlotLinearReport.C',R/'source/PlotLinearReport.C')
prior=base.read(base.OUT/'candidate/low_he4_gaussian_linear_nominal/command.json');cache=json.loads(prior[-1].split('+(',1)[1].split(',',1)[0])
cmd=['bash',str(base.P/'script/singularity_root_batch.sh'),str(base.P/base.C['runtimeMainconf']),'root4star','-l','-b','-q',str(R/'source/PlotLinearReport.C')+'+('+json.dumps(str(R))+','+json.dumps(cache)+')'];base.save(R/'plot_command.json',cmd)
with (R/'plot.log').open('w') as f:status=subprocess.run(cmd,cwd=base.P,stdout=f,stderr=subprocess.STDOUT).returncode
assert not status and 'LINEAR_REPORT_PLOTS_COMPLETE' in (R/'plot.log').read_text()
products={};trimCommands=[]
for name,pages in [('purity_vs_p_over_q',2),('pid_function_updates',2),('dedx_pid_map',1),('fit_low_d',1),('fit_low_t',1),('fit_low_he3',1),('fit_low_he4',1),('fit_anchors',1),('fit_tof',2)]:
 pdf=R/(name+'.pdf');p=subprocess.run(['gs','-dSAFER','-dNOPAUSE','-dBATCH','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)],capture_output=True,text=True);log=p.stdout+p.stderr;(R/(name+'.readback.log')).write_text(log);match=re.search(r'Processing pages 1 through (\d+)',log);assert not p.returncode and 'Error:' not in log and match and int(match.group(1))==pages,name
 p=subprocess.run(['gs','-q','-dNOPAUSE','-dBATCH','-dPDFSTOPONERROR','-sDEVICE=png16m','-r'+str(base.C['figureRasterDpi']),'-sOutputFile='+str(R/(name+'_%02d.png')),str(pdf)],capture_output=True,text=True);(R/(name+'.raster.log')).write_text(p.stdout+p.stderr);assert not p.returncode and 'Error:' not in p.stdout+p.stderr
 for page in range(1,pages+1):
  png=R/(name+'_%02d.png'%page);dst=png.with_name(png.stem+'_trim.png');crop=['convert',str(png),'-trim','+repage','-bordercolor','white','-border','10',str(dst)];subprocess.run(crop,check=True,capture_output=True);trimCommands.append(crop)
 products[pdf.name]=dict(pages=pages,sha256=base.sha(pdf));print('FIGURE_PASS',name,pages,flush=True)
base.save(R/'figure_trim_commands.json',trimCommands)
base.save(R/'figures_completed.json',dict(allPagesReadback=True,products=products))
