#!/usr/bin/env python3
import collections,json,pathlib,re,shutil,subprocess,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,C=run.R,run.C;examples=run.read(R/'examples.json');full=[dict(e,group='all_regions') for e in examples];run.write(R/'examples_all.tsv',full);representatives=[]
for region,values in C['representatives'].items():
 for target in values:
  e=min((e for e in examples if e['region']==region),key=lambda e:abs((e['rLo']+e['rHi'])/2-target));representatives.append(dict(e,group='representative_'+region))
run.write(R/'examples_report.tsv',representatives);run.save(R/'representatives.json',representatives);shutil.copy2(HERE/'PlotHe4.C',R/'source/PlotHe4.C')
cache=json.loads(run.read(R/'candidate_low_linear/command.json')[-1].split('+(',1)[1].split(',',1)[0]);cmd=['bash',str(run.P/'script/singularity_root_batch.sh'),str(run.P/C['runtimeMainconf']),'root4star','-l','-b','-q',str(R/'source/PlotHe4.C')+'+('+json.dumps(str(R))+','+json.dumps(cache)+')'];run.save(R/'plot_command.json',cmd)
with (R/'plot.log').open('w') as f:status=subprocess.run(cmd,cwd=run.P,stdout=f,stderr=subprocess.STDOUT).returncode
assert status==0 and 'HE4_PLOTS_COMPLETE' in (R/'plot.log').read_text()
counts={'purity_vs_p_over_q':2,'pid_functions':1,'dedx_pid_map':1,'fit_all_regions':199}
for region in C['regions']:counts['fit_'+region['name']]=sum(e['region']==region['name'] for e in examples);counts['fit_representative_'+region['name']]=1
products={}
for name,pages in counts.items():
 pdf=R/(name+'.pdf');p=subprocess.run(['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)],capture_output=True,text=True);log=p.stdout+p.stderr;(R/(name+'.readback.log')).write_text(log);match=re.search(r'Processing pages 1 through (\d+)',log);assert not p.returncode and 'Error:' not in log and match and int(match.group(1))==pages;products[pdf.name]=dict(pages=pages,sha256=run.sha(pdf));print('PDF_PASS',name,pages,flush=True)
 if name in ['purity_vs_p_over_q','pid_functions','dedx_pid_map'] or name.startswith('fit_representative'):
  cmd=['gs','-q','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-dUseCropBox','-sDEVICE=png16m','-r'+str(C['figureDpi']),'-sOutputFile='+str(R/(name+'_%02d.png')),str(pdf)];subprocess.run(cmd,check=True,capture_output=True)
  for page in range(1,pages+1):
   png=R/(name+'_%02d.png'%page);subprocess.run(['convert',str(png),'-trim','+repage','-bordercolor','white','-border','10',str(png.with_name(png.stem+'_trim.png'))],check=True,capture_output=True)
run.save(R/'pdf_validation.json',dict(allPagesReadback=True,products=products))
