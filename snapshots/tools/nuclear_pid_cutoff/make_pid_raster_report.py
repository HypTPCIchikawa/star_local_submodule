#!/usr/bin/env python3
"""Embed exact completed PNG chart pixels full-canvas; retain ROOT5 vector PDFs."""
import argparse
from pathlib import Path
import shutil
import subprocess
import compare_pid_selections as base


def main():
    p=argparse.ArgumentParser(description=__doc__)
    for key in ('configured','fitted','output-dir'):p.add_argument('--'+key,type=Path,required=True)
    a=p.parse_args();old=a.configured.resolve(strict=True);new=a.fitted.resolve(strict=True);out=a.output_dir.resolve()
    if out.exists() or out.parent!=new:raise ValueError('Require new child folder of completed fit report.')
    species=('d','t','he3','he4')
    oldnames=['configuration_table.png']+['response_'+s+'_full.png' for s in species]+['response_he3_lowR.png','response_he4_lowR.png','tpc_maps.png','tof_m2_maps.png']
    newnames=['fitted_parameter_definitions.png']+['fitted_'+s+'_'+r+'.png' for r in ('low_ALL','low_adopted_ALL','high_native_TOF') for s in species]+['fitted_he4_high_coarse_TOF.png']
    pages=[];inputs={}
    for folder,names in ((old,oldnames),(new,newnames)):
        manifest=folder/'completed.json';inputs[str(manifest)]=base.sha(manifest);d=base.read(manifest)
        for name in names:
            path=folder/name;sha=base.sha(path)
            if sha!=d['products'][name]:raise ValueError('Changed declared PNG '+str(path))
            pages.append(path);inputs[str(path)]=sha
    out.mkdir();(out/'source').mkdir()
    for path in (Path(__file__).resolve(),Path(base.__file__).resolve()):shutil.copy2(path,out/'source'/path.name)
    pdf=out/'new_old_pid_final.pdf';command=['convert','-units','PixelsPerInch','-density','100']+[str(v) for v in pages]+['-compress','Zip',str(pdf)]
    base.write(out/'command.json',command)
    with (out/'convert.log').open('x') as log:subprocess.run(command,stdout=log,stderr=subprocess.STDOUT,check=True,timeout=180)
    q=subprocess.run(['gs','-q','-dNODISPLAY','-dNOSAFER','-c','('+str(pdf)+') (r) file runpdfbegin pdfpagecount = quit'],stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,check=True,timeout=60)
    if int(q.stdout.strip())!=len(pages) or q.stderr.strip():raise ValueError('PDF page count/read diagnostics: '+q.stdout+q.stderr)
    q=subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-sDEVICE=nullpage',str(pdf)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,text=True,check=True,timeout=60)
    if q.stdout.strip():raise ValueError('PDF diagnostics: '+q.stdout)
    for n in (1,len(pages)):
        subprocess.run(['gs','-q','-dBATCH','-dNOPAUSE','-dFirstPage='+str(n),'-dLastPage='+str(n),'-sDEVICE=png16m','-r90','-sOutputFile='+str(out/('pdf_page_%02d.png'%n)),str(pdf)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,check=True,timeout=60)
    if any(base.sha(Path(path))!=sha for path,sha in inputs.items()):raise ValueError('PNG/manifest changed during assembly.')
    base.write(out/'context.json',dict(pageCount=len(pages),pageOrder=[str(v) for v in pages],inputs=inputs,noNewFit=True,noPidChange=True,rasterBackedReport=True,pixelValuesUnchanged=True,originalVectorPdf=str(new/'new_old_pid_final.pdf'),why='ROOT5 vector-PDF input diagnostics and canvas offsets; exact completed PNG pixels embedded full-canvas. Original PDFs preserved.'))
    base.write(out/'completed.json',dict(executionComplete=True,pageCount=len(pages),allPagesReadWithoutDiagnostics=True,protectedInputsVerified=True,products={str(v.relative_to(out)):base.sha(v) for v in out.rglob('*') if v.is_file()}))
    print('PID_FULL_CANVAS_REPORT_COMPLETE',pdf)


if __name__=='__main__':main()
