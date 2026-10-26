#!/usr/bin/env python3
import argparse,re,subprocess
from study import OUT,save,sha
p=argparse.ArgumentParser();p.add_argument('stage',choices=['low','high','summary']);args=p.parse_args();folder=OUT/('atlas_'+args.stage+'_001') if args.stage!='summary' else OUT
expected={('fit_atlas_'+args.stage+'_'+s+'.pdf'):(160 if args.stage=='low' else 300 if s!='he4' else 60) for s in ['d','t','he3','he4']} if args.stage!='summary' else {'purity_vs_p_over_q.pdf':2,'dedx_vs_p_over_q_curves.pdf':8}
verified={}
for name,pages in expected.items():
 pdf=folder/name;log=folder/(name+'.readback.log')
 with log.open('w') as f:status=subprocess.run(['gs','-dSAFER','-dBATCH','-dNOPAUSE','-dPDFSTOPONERROR','-sDEVICE=nullpage',str(pdf)],stdout=f,stderr=subprocess.STDOUT).returncode
 text=log.read_text();m=re.search(r'Processing pages 1 through (\d+)',text)
 if status or not m or int(m.group(1))!=pages or 'Error:' in text:raise RuntimeError('PDF readback failed '+str(pdf))
 verified[name]=dict(pages=pages,bytes=pdf.stat().st_size,sha256=sha(pdf));print('PDF_PASS',name,pages,flush=True)
save(folder/('pdf_validation_'+args.stage+'.json'),dict(executionComplete=True,allPagesReadback=True,products=verified))
