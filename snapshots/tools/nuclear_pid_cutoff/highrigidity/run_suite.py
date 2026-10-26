#!/usr/bin/env python3
"""Run explicit high-R cases, keeping each model in its own immutable output."""
import argparse
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
import subprocess
import sys
import run_highrigidity as high


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',required=True,type=Path)
    parser.add_argument('--stage',required=True,choices=('pilot','upper-pilot','native'))
    parser.add_argument('--targets',required=True,nargs='+',choices=tuple(high.MAPPING))
    parser.add_argument('--variations',required=True,nargs='+')
    parser.add_argument('--suffix',required=True)
    parser.add_argument('--workers',required=True,type=int)
    parser.add_argument('--reuse-complete',action='store_true')
    parser.add_argument('--summarize-to',type=Path)
    args=parser.parse_args();cfg=high.checked_settings(high.common.resolve(args.settings))
    if args.workers<1 or args.workers>2 or not args.suffix.isalnum(): raise ValueError('Explicit workers 1/2 and alphanumeric suffix required.')
    if not set(args.variations)<=set(cfg['shapeVariations']): raise ValueError('Undeclared variation.')
    output=high.common.PROJECT/cfg['outputBase']
    jobs=[]
    for target in args.targets:
        for response in cfg['responseModels']:
            for background in cfg['backgroundModels']:
                for variation in args.variations:
                    name='_'.join((args.stage.replace('-','_'),target,response,background))
                    if variation!='nominal': name+='_'+variation
                    folder=output/(name+'_'+args.suffix)
                    command=[sys.executable,'-B',str(high.HERE/'run_highrigidity.py'),'--settings',str(args.settings),
                        '--target',target,'--stage',args.stage,'--response',response,'--background',background,
                        '--variation',variation,'--output-dir',str(folder)]
                    jobs.append((folder,command))
    def execute(job):
        folder,command=job
        if folder.exists():
            done=high.common.read(folder/'completed.json')
            if not args.reuse_complete or not done.get('executionComplete') or not done.get('protectedSourcesVerified'):
                raise ValueError('Existing incomplete/nonreusable output: '+str(folder))
            for name in ('cutoff_fits.json','cutoff_fits.root','context.json','native_summary.json'):
                if high.common.digest(folder/name)!=done['products'][name]: raise ValueError('Changed completed fit: '+str(folder/name))
            print('HIGH_R_SUITE_REUSE',folder,flush=True)
        else:
            print('HIGH_R_SUITE_START',folder,flush=True)
            subprocess.run(command,cwd=high.common.PROJECT,check=True)
        return str(folder)
    with ThreadPoolExecutor(max_workers=args.workers) as pool: results=list(pool.map(execute,jobs))
    if args.summarize_to:
        command=[sys.executable,'-B',str(high.HERE/'summarize_highrigidity.py'),'--settings',str(args.settings),
            '--stage',args.stage,'--output-dir',str(args.summarize_to),'--variations']+args.variations+['--runs']+results
        subprocess.run(command,cwd=high.common.PROJECT,check=True)
    print('HIGH_R_SUITE_COMPLETE cases='+str(len(results)),flush=True)


if __name__=='__main__': main()
