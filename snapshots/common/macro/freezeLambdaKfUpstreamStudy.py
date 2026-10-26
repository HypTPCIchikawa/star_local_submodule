#!/usr/bin/env python3
"""Create new, frozen best50 jobs for 10k replay and independent 100k validation.

Does not run events, modify defaults, or rank on validation data. Both energy
decisions must already be complete. Refuses existing output configuration.
"""
import argparse
import hashlib
import json
import os
from pathlib import Path
import yaml

REPO = Path(__file__).resolve().parents[2]

def new(path, value):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('x') as stream:
        stream.write(value)

def dump(path, obj):
    new(path, yaml.safe_dump(obj, sort_keys=False))

def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def prepare_scores(master):
    for energy in master['energies']:
        out = REPO/'rootfile'/(energy+'_anaLambda_KFParticle')/master['studyTag']
        meta = (out/'exports/default/metadata.txt').read_text()
        reco = int(next(s.split(':',1)[1] for s in meta.splitlines() if s.startswith('reconstructedEvents:')))
        config = dict(schemaVersion=1, baseConfig=master['baseScoreConfig'],
            eventIdsFile=str(out/'provenance/input/scan_event_ids.tsv'), inputEvents=master['scanEvents'],
            reconstructedEvents=reco, grid=master['finalGrid'],
            title=energy+' FXT: upstream/downstream KF Lambda purity exploration')
        for key in ('includeAllPairs','tripleAnchor','retentionFloors','minimumWindowCount',
                    'poissonUpperConfidenceLevel','maxPointsPerDataset','bootstrapReplicates',
                    'bootstrapSeed','backgroundCheck'):
            config[key]=master[key]
        dump(out/'configs/score.yaml',config)
        new(out/'configs/datasets.tsv','id\texport_dir\n'+''.join(
            '{}\t{}\n'.format(variant,out/'exports'/variant) for variant in master['upstreamVariants']))

def freeze(master):
    decisions=[]
    for energy in master['energies']:
        out=REPO/'rootfile'/(energy+'_anaLambda_KFParticle')/master['studyTag']
        status=json.loads((out/'scan/status.json').read_text())
        if not status['complete']: raise RuntimeError('Incomplete discovery '+energy)
        chosen=json.loads((out/'scan/chosen_points.json').read_text())
        best=next(x for x in chosen['choices'] if x['choice_id']=='best50')
        if (out/'configs/frozen').exists(): raise RuntimeError('Refusing existing frozen config')
        decisions.append((energy,out,best))
    for energy,out,best in decisions:
        folder=out/'configs/frozen'
        original=yaml.safe_load((out/'configs/default_kf.yaml').read_text())
        overrides=best['combinedCutOverrides']
        effective=json.loads((out/'scan/sources.json').read_text())['datasets'][0]['effectiveCuts']
        # Optional valid cuts may be absent in YAML but present in effective metadata.
        if not set(overrides).issubset(effective): raise RuntimeError('Unknown cut override')
        dump(folder/'default_kf.yaml',original)
        dump(folder/'best50_kf.yaml',dict(original,**overrides))
        jobs=[]
        for sample,events in (('replay10k',master['scanEvents']),('validation100k',master['validationEvents'])):
            kinds=('best50',) if sample=='replay10k' else ('default','best50')
            input_path=out/'provenance/input'/('scan.list' if sample=='replay10k' else 'validation.list')
            for kind in kinds:
                job_id=kind+'_'+sample
                info=yaml.safe_load((out/'configs/default_analysis.yaml').read_text())
                info['dataset']['allPicoDstList']=str(input_path)
                name=energy+'_KF_'+master['studyTag']+'_'+job_id
                for key in ('anaName','name','jobName','scratchSubdir','outputFileStem'):
                    info['analysis'][key]=name
                info_path=folder/(job_id+'_analysis.yaml')
                dump(info_path,info)
                main=yaml.safe_load((out/'configs/default_main_resolved.yaml').read_text())
                main['kf']=os.path.relpath(folder/(kind+'_kf.yaml'),REPO/'config')
                main['analysis']=os.path.relpath(info_path,REPO/'config')
                main_path=folder/(job_id+'_main.yaml')
                dump(main_path,main)
                jobs.append(dict(id=job_id,energy=energy,mainconf=str(main_path),input=str(input_path),
                    events=events,output=str(out/sample/(kind+'.root')),
                    log=str(out/'provenance/logs'/(job_id+'.log')),
                    export=str(out/sample/'exports'/kind),overrides=overrides if kind=='best50' else {}))
        new(folder/'jobs.json',json.dumps(jobs,indent=2)+'\n')
        manifest=dict(energy=energy,choice=best,selectionEvents=master['scanEvents'],
            validationEvents=master['validationEvents'],validationExaminedBeforeFreeze=False,
            decisionPath=str(out/'scan/chosen_points.json'),decisionSha256=sha(out/'scan/chosen_points.json'),
            selectionRule='best50: maximize worst-sideband Poisson-upper-background purity, min S retention >= 0.50',
            noRetuningOnValidation=True)
        new(folder/'frozen_decision.json',json.dumps(manifest,indent=2)+'\n')
        products=sorted(p for p in folder.iterdir() if p.is_file())
        new(folder/'SHA256SUMS',''.join('{}  {}\n'.format(sha(p),p.name) for p in products))
        print('FROZEN',energy,best['dataset_id'],overrides,flush=True)

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--config',required=True,type=Path)
    parser.add_argument('--prepare-scores',action='store_true')
    args=parser.parse_args()
    master=yaml.safe_load(args.config.read_text())
    if args.prepare_scores: prepare_scores(master)
    else: freeze(master)

if __name__=='__main__': main()
