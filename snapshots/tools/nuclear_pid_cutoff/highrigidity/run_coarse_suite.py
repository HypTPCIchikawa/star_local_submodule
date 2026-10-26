#!/usr/bin/env python3
"""Four declared coarse-R models, retaining each independent fit product."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys
import yaml
import run_highrigidity as high


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',required=True,type=Path)
    parser.add_argument('--stage',required=True,choices=('pilot','native'))
    parser.add_argument('--suffix',required=True)
    parser.add_argument('--workers',required=True,type=int)
    parser.add_argument('--variations',required=True,nargs='+')
    args=parser.parse_args()
    settings=high.common.resolve(args.settings);cfg=yaml.safe_load(settings.read_text())
    base=high.checked_settings(high.common.resolve(cfg['baseSettings']))
    if args.workers not in (1,2) or not args.suffix.isalnum() or not set(args.variations)<=set(base['shapeVariations']):
        raise ValueError('Invalid explicit worker/suffix/variation.')
    root=high.common.PROJECT/cfg['outputBase'];jobs=[]
    for response in base['responseModels']:
        for background in base['backgroundModels']:
            for variation in args.variations:
                name='_'.join(('coarse',args.stage,'he4',response,background))
                if variation!='nominal': name+='_'+variation
                folder=root/(name+'_'+args.suffix)
                command=[sys.executable,'-B',str(high.HERE/'run_coarse.py'),'--settings',str(settings),
                    '--stage',args.stage,'--response',response,'--background',background,'--variation',variation,'--output-dir',str(folder)]
                jobs.append((folder,command))
    def execute(job):
        folder,command=job
        if folder.exists(): raise ValueError('Existing coarse fit will not be overwritten: '+str(folder))
        print('COARSE_SUITE_START',folder,flush=True)
        subprocess.run(command,cwd=high.common.PROJECT,check=True)
        return str(folder)
    with ThreadPoolExecutor(max_workers=args.workers) as pool: results=list(pool.map(execute,jobs))
    print('COARSE_HIGH_R_SUITE_COMPLETE cases='+str(len(results)),flush=True)
