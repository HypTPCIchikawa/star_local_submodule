#!/usr/bin/env python3
import fcntl,json,os,sys,traceback
from pathlib import Path
S=Path(__file__).resolve().parent
sys.path.insert(0,str(S))
import watch_partial as m
from accepted_compare import compare

def summary(path,files):
    rows=path.read_text().splitlines()
    if not rows or rows[-1]!='DONE\t'+str(files):raise RuntimeError('Incomplete audit '+str(path))
    out={}
    for row in rows:
        f=row.split('\t')
        if f[0]=='HIST':out[f[1]]=(f[2],tuple(float(x) for x in f[3:]))
    if not out:raise RuntimeError('Empty audit')
    return out

def accumulate(total,part):
    for name,(schema,values) in part.items():
        if name in total:
            old_schema,old_values=total[name]
            if schema!=old_schema:raise RuntimeError('Schema changed')
            values=tuple(a+b for a,b in zip(values,old_values))
        total[name]=(schema,values)

def validate(expected,actual,policy,path):
    waived=compare(expected,actual,policy)
    m.notify.save_json(path,{'checkedAt':m.notify.now(),'histograms':len(actual),'acceptedDifferences':waived,'acceptance':str(S/'acceptance.json')})

def merge(paths,label,expected,staging,runner,policy,status,saved=None):
    a=S/label;a.mkdir()
    m.notify.save_json(a/'inputs.json',{'inputs':m.fingerprints(paths)})
    status.update(phase=label,checkedAt=m.notify.now());m.notify.save_json(S/'status.json',status)
    (a/'merge_notify_audit.C').write_text(m.ROOT_AUDIT)
    if saved is None:
        target=staging/(label+'.root')
        m.run_logged(runner+['hadd',str(target)]+[str(p) for p in paths],a/'hadd.log')
    else:target=saved
    actual=m.audit([target],'output',a,runner)
    validate(expected,actual,policy,a/'accepted_comparison.json')
    if m.fingerprints(paths)!=json.loads((a/'inputs.json').read_text())['inputs']:raise RuntimeError('Intermediate inputs changed')
    m.notify.save_json(a/'completed.json',{'path':str(target),'fingerprint':m.archive_helper.regular(target),'actualAudit':str(a/'output.tsv'),'acceptedComparison':str(a/'accepted_comparison.json')})
    return target,actual,a/'output.tsv'

def main():
    lock=(S/'monitor.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    old=S.parent/'partial_merge_checkpoint_20261023';prior=S.parent/'partial_merge_20261023'
    original=json.loads((old/'input_manifest.json').read_text())['inputs'];paths=[Path(x['path']) for x in original]
    policy=json.loads((S/'acceptance.json').read_text());config=json.loads((prior/'settings.json').read_text())
    runmeta=Path(config['runmeta']);meta=json.loads(runmeta.read_text());repo=Path(meta['projectRoot'])
    status=json.loads((old/'result.json').read_text())
    status.pop('error',None);status.pop('mergeRound',None)
    status.update(phase='checking_saved_chunks',checkedAt=m.notify.now(),stateDirectory=str(S),completedChunks=270,userAcceptedDifferences=str(S/'acceptance.json'),continuedFrom=str(old))
    if (S/'result.json').exists() or (S/'delivery_attempt.json').exists():raise RuntimeError('Already attempted; inspect result')
    dest=Path(status['plannedOutput'])
    if dest.exists():raise RuntimeError('Merged output already exists')
    m.notify.save_json(S/'process.json',{'pid':os.getpid(),'host':m.socket.getfqdn(),'startedAt':m.notify.now(),'mode':'continue accepted seven-histogram differences'})
    m.notify.save_json(S/'status.json',status)
    staging=dest.parent/('.merge_staging_'+S.name);staging.mkdir()
    runner=['bash',str(repo/'script/singularity_root_batch.sh'),str(repo/'config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml')]
    try:
        if m.fingerprints(paths)!=original:raise RuntimeError('Original inputs changed')
        m.notify.save_json(S/'input_manifest.json',{'inputs':original});m.write_list(S/'inputs.list',paths)
        parts=[];part_summaries=[];total={}
        for i in range(270):
            c=old/('chunk_%04d'%i);info=json.loads((c/'completed.json').read_text());p=Path(info['path'])
            m.archive_helper.verify_fingerprint(p,info['fingerprint'])
            e=Path(info['expectedAudit']);a=e.parent/'output.tsv'
            actual=summary(a,1);parts.append(p);part_summaries.append(actual)
            accumulate(total,summary(e,len(info['inputs'])))
        outputs=[];sums=[]
        for i,b in enumerate(range(0,len(parts),100)):
            expected={}
            for d in part_summaries[b:b+100]:accumulate(expected,d)
            saved=None
            if i==0:
                prior_attempt=old/'round_00_0000/attempt_000'
                if m.fingerprints(parts[b:b+100])!=json.loads((prior_attempt/'started.json').read_text())['inputs']:raise RuntimeError('Previously merged inputs changed')
                row=(prior_attempt/'output.tsv').read_text().splitlines()[0].split('\t');saved=Path(row[1])
            out,actual,report=merge(parts[b:b+100],'intermediate_%d'%i,expected,staging,runner,policy,status,saved)
            outputs.append(out);sums.append(actual)
        expected={}
        for d in sums:accumulate(expected,d)
        merged,actual,report=merge(outputs,'final',expected,staging,runner,policy,status)
        validate(total,actual,policy,S/'final_vs_original_comparison.json')
        comparison=json.loads((prior/'event_comparison_and_authorization.json').read_text())
        counters={'inputEventsFromRefMultEntries':actual['hRefMult;1'][1][0],'reconstructedEventsFromRefMultCorrEntries':actual['hRefMultCorr;1'][1][0]}
        if counters['inputEventsFromRefMultEntries']!=comparison['currentInputEvents'] or counters['reconstructedEventsFromRefMultCorrEntries']!=comparison['currentReconstructedEvents']:raise RuntimeError('Merged event counters differ')
        m.notify.save_json(S/'event_counter_validation.json',counters)
        if m.fingerprints(paths)!=original:raise RuntimeError('Original input changed before publication')
        (S/'merged_audit.tsv').write_bytes(report.read_bytes())
        m.publish_no_clobber(merged,dest)
        status.update(phase='merged',checkedAt=m.notify.now(),merge={'path':str(dest),'size':dest.stat().st_size,'fingerprint':m.archive_helper.regular(dest),'histogramCount':len(actual),'validation':'ROOT recursive read; exact Entries and accepted event counts; additive summaries with seven histograms explicitly waived by user for this run','inputCount':len(paths),'publishedAt':m.notify.now(),'stagingDirectory':str(staging)})
    except Exception as error:
        status.update(phase='failed',checkedAt=m.notify.now(),error=str(error));traceback.print_exc()
    m.notify.save_json(S/'result.json',status);m.notify.save_json(S/'status.json',status)
    m.archive_after_merge(status,runmeta,S)
    msg=m.make_message(config['email'],meta,status)
    good=status['phase']=='merged' and status.get('archive',{}).get('phase')=='archived'
    msg.replace_header('Subject','[STAR 部分hadd完了・承認済み]' if good else '[STAR hadd終了・要確認]')
    msg.replace_header('Message-ID','<star-'+meta['jobid']+'-accepted-20261025@'+m.socket.getfqdn()+'>')
    msg.set_content(msg.get_content()+'\n今回の7ヒストグラムのbin合計・重み付き合計差はユーザー承認により許容しました。Entries一致と残りヒストの加算検査を維持。承認記録: '+str(S/'acceptance.json')+'\n')
    (S/'message.eml').write_bytes(msg.as_bytes());m.notify.send_once(S,msg,status,config['email'])
    print(m.notify.now(),'Result mail accepted',status['phase'],flush=True)
    return 0 if good else 1
if __name__=='__main__':sys.exit(main())
