#!/usr/bin/env python3
"""Reproducible saved-fit views and response tests in ROOT5; study output only."""
import argparse
from pathlib import Path
import math
import re
import shutil
import subprocess
from prepare_tag_study import PROJECT,owned,read,save,sha,tsv,quote,require

SOURCES=('DrawPidRecalibration.C','DrawSharedPidFits.C','TestAllFixedResponseFit.C','AllFixedResponseFit.h','Float16BinMapping.h','SplitNormalFit.h','PeakResponseFit.h')

def freeze(out):
    require(not out.exists(),'Choose a fresh report output')
    out.mkdir(parents=True)
    paths=[Path(__file__).resolve()]+[Path(__file__).resolve().with_name(n) for n in SOURCES]+[PROJECT/'tools/nuclear_pid_tof/NuclearFitCommon.h']
    hashes={}
    for p in paths:
        q=out/'source'/p.relative_to(PROJECT);q.parent.mkdir(parents=True,exist_ok=True);shutil.copy2(str(p),str(q));hashes[str(q)]=sha(q)
    return hashes

def run_macro(out,main,source,args,marker):
    macro=out/'source/tools/nuclear_pid_restart_20261020'/source
    cmd=['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(main),'root4star','-l','-b','-q',str(macro)+'+('+','.join(quote(p) for p in args)+')']
    save(out/'command.json',dict(argv=cmd,cwd=str(PROJECT)))
    with (out/'run.log').open('x') as stream:
        subprocess.run(cmd,cwd=PROJECT,stdout=stream,stderr=subprocess.STDOUT,check=True,timeout=1800)
    require(marker in (out/'run.log').read_text(),'Missing completed-output marker: '+str(out/'run.log'))

def pdfcheck(out,pdf,expected):
    result=subprocess.run(['gs','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)],stdout=subprocess.PIPE,stderr=subprocess.STDOUT,universal_newlines=True,check=True)
    pages=re.search(r'Processing pages 1 through ([0-9]+)\.',result.stdout)
    require(pages and int(pages.group(1))==expected and 'Error:' not in result.stdout,'PDF page/readback check failed')
    with (out/'pdf_readback.log').open('x') as stream:stream.write(result.stdout)

def execute(action,fitdir,coeffdir,out):
    plan=read(fitdir/'plan.json');cfg=plan['config'];main=owned(cfg['mainconf']);protected=freeze(out)
    protected.update({str(fitdir/n):sha(fitdir/n) for n in ('plan.json','settings.tsv')})
    if action=='tests':
        protected[str(fitdir/'settings.tsv')]=sha(fitdir/'settings.tsv')
        run_macro(out,main,'TestAllFixedResponseFit.C',[fitdir/'settings.tsv'],'ALL_FIXED_RESPONSE_TEST_COMPLETE')
    elif action=='fits':
        result=read(fitdir/'fits.json');records={q['id']:q for q in plan['slices']};rows=[]
        protected.update({str(fitdir/n):sha(fitdir/n) for n in ('fits.json','native_spectra.root')})
        for i,q in enumerate(result['results']):
            if not q['fitCompleted']:continue
            expected=records[q['id']];chi=q.get('pearsonPerNdf');chi=chi if isinstance(chi,(float,int)) and math.isfinite(chi) else -1.
            gate=q.get('winnerGate',{});p=gate.get('candidatePurity');e=gate.get('conditionalError')
            note='P='+('unassessed' if not q['purityAvailable'] else '%.4f +/- %.4f'%(p,e))+'; '+','.join(q.get('purityQualityReasons',[]))
            rows.append((i,q['id'],q['target'],q['mode'],(q['rLo']+q['rHi'])/2.,q['fitLo'],q['fitHi'],expected['mean'],expected['sigma'],chi,int(q['purityAvailable']),note or 'assessed'))
        tsv(out/'metrics.tsv',('index','id','target','mode','r','fitLo','fitHi','mean','sigma','chi2Ndf','numerical','reason'),rows)
        tsv(out/'view_settings.tsv',('key','value'),dict(displayRebin=1,viewPaddingSigma=.05,linearPadFraction=.25,logMin=.5,maxCasesPerPage=2).items())
        pdf=out/'shared_curve_fit_diagnostics.pdf';run_macro(out,main,'DrawSharedPidFits.C',[fitdir/'native_spectra.root',out/'metrics.tsv',out/'view_settings.tsv',pdf],'SHARED_PID_FIT_DRAW_COMPLETE rows='+str(len(rows)));pdfcheck(out,pdf,(len(rows)+1)//2)
        protected[str(pdf)]=sha(pdf)
    else:
        require(coeffdir is not None,'Explicit coefficient directory required')
        old=Path(read(coeffdir.parent/'completed.json')['sourceFitPath']).parent/'original_curves.tsv';new=coeffdir/'curves.tsv'
        for p in (old,new,Path(plan['cache'])):protected[str(p)]=sha(p)
        tsv(out/'overlay_settings.tsv',('key','value'),dict(allKey='ALL',tagDKey='TAG_D',tagTKey='TAG_T',tagHe3Key='TAG_HE3',rLowLo=.4,rLowHi=2.,rHighLo=2.,rHighHi=5.,dedxLowLo=0.,dedxLowHi=120.,dedxHighLo=0.,dedxHighHi=25.,selectionNSigma=2.).items())
        pdf=out/'pid_polynomial_colz_comparison.pdf';run_macro(out,main,'DrawPidRecalibration.C',[plan['cache'],old,new,out/'overlay_settings.tsv',pdf],'PID_RECALIBRATION_DRAW_COMPLETE pages=8');pdfcheck(out,pdf,8);protected[str(pdf)]=sha(pdf)
    for path,digest in protected.items():require(sha(Path(path))==digest,'Protected source/input changed')
    save(out/'completed.json',dict(status='completed',action=action,productionEnabled=False,protectedReadback=True,products=protected))
    print('SHARED_PID_REPORT_COMPLETE '+action+' '+str(out),flush=True)

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('action',choices=('tests','fits','overlay'));p.add_argument('--fit-dir',required=True);p.add_argument('--coefficient-dir');p.add_argument('--output-dir',required=True);a=p.parse_args();execute(a.action,owned(a.fit_dir),owned(a.coefficient_dir) if a.coefficient_dir else None,owned(a.output_dir,False))
if __name__=='__main__':main()
