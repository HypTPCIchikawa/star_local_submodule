#!/usr/bin/env python3
import argparse,fcntl,json,os,sys,types,traceback
from pathlib import Path
S=Path(__file__).resolve().parent
sys.path.insert(0,str(S))
import watch_partial as m

def read_audit(path):
    data={}
    for line in path.read_text().splitlines():
        f=line.split('\t')
        if f[0]=='HIST':data[f[1]]=(f[2],tuple(float(x) for x in f[3:]))
    if not data:raise RuntimeError('Missing histogram audit '+str(path))
    return data

def add_summaries(total,part):
    for name,(schema,values) in part.items():
        if name in total:
            old_schema,old_values=total[name]
            if schema!=old_schema:raise RuntimeError('Schema mismatch '+name)
            values=tuple(a+b for a,b in zip(old_values,values))
        total[name]=(schema,values)

def checkpoint(paths,label,staging,runner):
    state=S/label;state.mkdir(exist_ok=True)
    done=state/'completed.json'
    before=m.fingerprints(paths)
    if done.exists():
        info=json.loads(done.read_text())
        if before!=info['inputs']:raise RuntimeError('Changed inputs '+label)
        m.archive_helper.verify_fingerprint(Path(info['path']),info['fingerprint'])
        return Path(info['path']),read_audit(Path(info['expectedAudit']))
    # Never overwrite an interrupted attempt or assume its output was valid.
    i=0
    while (state/('attempt_%03d'%i)).exists():i+=1
    attempt=state/('attempt_%03d'%i);attempt.mkdir()
    m.notify.save_json(attempt/'started.json',{'inputs':before,'at':m.notify.now()})
    (attempt/'merge_notify_audit.C').write_text(m.ROOT_AUDIT)
    expected=m.audit(paths,'input',attempt,runner)
    target=staging/(label+'_attempt_%03d.root'%i)
    m.run_logged(runner+['hadd',str(target)]+[str(p) for p in paths],attempt/'hadd.log')
    actual=m.audit([target],'output',attempt,runner)
    m.compare_histograms(expected,actual)
    if m.fingerprints(paths)!=before:raise RuntimeError('Input changed during checkpoint')
    m.notify.save_json(done,{'inputs':before,'path':str(target),'fingerprint':m.archive_helper.regular(target),'expectedAudit':str(attempt/'input.tsv'),'at':m.notify.now()})
    return target,expected

def main():
    parser=argparse.ArgumentParser();parser.add_argument('--toy',type=Path);args=parser.parse_args()
    lock=(S/'monitor.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    repo=Path('/gpfs/mnt/gpfs01/star/pwg/yichikawa/git_STAR_3/star-analyzer')
    if args.toy:
        paths=sorted(args.toy.glob('toy_*.root'));indices=list(range(len(paths)))
        dest=args.toy/'toy_merge.root';batch=2
        status={'phase':'starting','expectedSubjobs':len(paths),'selectedIndices':indices,'excludedIndices':[],'partial':False,'stateDirectory':str(S),'outcomes':{'exit_0':len(paths)}}
        meta={};runmeta=None
    else:
        prior=S.parent/'partial_merge_20261023'
        status=json.loads((prior/'status.json').read_text())
        attempt=json.loads((prior/'merge_attempt.json').read_text())
        paths=[Path(p['path']) for p in attempt['inputs']]
        if m.fingerprints(paths)!=attempt['inputs']:raise RuntimeError('Production inputs changed since original attempt')
        indices=status['selectedIndices'];dest=Path(status['plannedOutput']);batch=100
        runmeta=Path(json.loads((prior/'settings.json').read_text())['runmeta']);meta=json.loads(runmeta.read_text())
        status.update(stateDirectory=str(S),previousAttempt=str(prior))
    if (S/'result.json').exists():raise RuntimeError('Result already exists; inspect before resuming')
    if dest.exists():raise RuntimeError('Destination already exists')
    staging=dest.parent/('.merge_staging_'+S.name);staging.mkdir(exist_ok=True)
    runner=['bash',str(repo/'script/singularity_root_batch.sh'),str(repo/'config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml')]
    if not (S/'input_manifest.json').exists():m.notify.save_json(S/'input_manifest.json',{'inputs':m.fingerprints(paths)})
    if not (S/'inputs.list').exists():m.write_list(S/'inputs.list',paths)
    m.notify.save_json(S/'process.json',{'pid':os.getpid(),'host':m.socket.getfqdn(),'startedAt':m.notify.now(),'mode':'checkpointed partial hadd','mailAfterValidation':not bool(args.toy)})
    try:
        total={};parts=[]
        for i,begin in enumerate(range(0,len(paths),batch)):
            status.update(phase='auditing_and_merging_chunks',checkedAt=m.notify.now(),completedChunks=i,totalChunks=(len(paths)+batch-1)//batch)
            m.notify.save_json(S/'status.json',status)
            target,expected=checkpoint(paths[begin:begin+batch],'chunk_%04d'%i,staging,runner)
            parts.append(target);add_summaries(total,expected)
        round_no=0
        while len(parts)>1:
            status.update(phase='merging_intermediates',checkedAt=m.notify.now(),mergeRound=round_no)
            m.notify.save_json(S/'status.json',status)
            parts=[checkpoint(parts[b:b+batch],'round_%02d_%04d'%(round_no,b//batch),staging,runner)[0] for b in range(0,len(parts),batch)]
            round_no+=1
        (S/'merge_notify_audit.C').write_text(m.ROOT_AUDIT)
        actual=m.audit(parts,'merged_audit',S,runner);m.compare_histograms(total,actual)
        if not args.toy:
            comparison=json.loads((S.parent/'partial_merge_20261023/event_comparison_and_authorization.json').read_text())
            counters={'inputEventsFromRefMultEntries':actual['hRefMult;1'][1][0], 'reconstructedEventsFromRefMultCorrEntries':actual['hRefMultCorr;1'][1][0]}
            if counters['inputEventsFromRefMultEntries']!=comparison['currentInputEvents'] or counters['reconstructedEventsFromRefMultCorrEntries']!=comparison['currentReconstructedEvents']:
                raise RuntimeError('Merged event counters differ from accepted comparison')
            m.notify.save_json(S/'event_counter_validation.json',counters)
        # Re-check original input fingerprints before publishing and archiving.
        original=json.loads((S/'input_manifest.json').read_text())['inputs']
        if m.fingerprints(paths)!=original:raise RuntimeError('Original input changed')
        m.publish_no_clobber(parts[0],dest)
        status.update(phase='merged',checkedAt=m.notify.now(),merge={'path':str(dest),'size':dest.stat().st_size,'fingerprint':m.archive_helper.regular(dest),'histogramCount':len(actual),'validation':'ROOT recursive read + per-checkpoint and final additive histogram fingerprints','inputCount':len(paths),'publishedAt':m.notify.now(),'stagingDirectory':str(staging)})
    except Exception as error:
        status.update(phase='failed',checkedAt=m.notify.now(),error=str(error));traceback.print_exc()
    m.notify.save_json(S/'result.json',status);m.notify.save_json(S/'status.json',status)
    if args.toy:return 0 if status['phase']=='merged' else 1
    m.archive_after_merge(status,runmeta,S)
    return m.finish(status,meta,types.SimpleNamespace(email='yudai.ichikawa.d3@tohoku.ac.jp'),S)
if __name__=='__main__':sys.exit(main())
