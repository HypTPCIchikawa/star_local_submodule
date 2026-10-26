#!/usr/bin/env python3
import concurrent.futures,fcntl,json,os,sys,traceback,hashlib,datetime
from pathlib import Path
S=Path(__file__).resolve().parent
sys.path.insert(0,str(S))
import verification as v
m=v.m

def invoke(mode,inp,out,report,runner,log):
    macro=str(S/'split_sparse.C')+'+('+','.join(json.dumps(str(x)) for x in (mode,inp,out,'hNuclearPidPairJoint',report))+')'
    m.run_logged(runner+['root4star','-b','-q',macro],log)

def main():
    lock=(S/'monitor.lock').open('a');fcntl.flock(lock,fcntl.LOCK_EX|fcntl.LOCK_NB)
    old=S.parent/'partial_merge_native_verified_20261025';prior=S.parent/'partial_merge_20261023'
    status=json.loads((old/'result.json').read_text());status.pop('error',None)
    config=json.loads((prior/'settings.json').read_text());runmeta=Path(config['runmeta']);meta=json.loads(runmeta.read_text());repo=Path(meta['projectRoot'])
    status.update(phase='preparing_split_outputs',checkedAt=m.notify.now(),stateDirectory=str(S),sparseStorage='3 separate ROOT files; original class, axes, contents, Sumw2 preserved',formatApproval=str(S/'format_approval.json'))
    if (S/'result.json').exists() or (S/'delivery_attempt.json').exists():raise RuntimeError('Prior attempt exists; inspect before restarting')
    dest=Path(status['plannedOutput']);directory=dest.parent;stem=meta['anaName']+'_'+meta['jobid']
    part_dest=[directory/(stem+'_pid_pair_joint_part%03d.root'%i) for i in range(3)]
    manifest_dest=directory/(stem+'_outputs.json');list_dest=directory/(stem+'_pid_pair_joint_parts.list');helper_dest=directory/(stem+'_project_joint_parts.C')
    for p in [dest,manifest_dest,list_dest,helper_dest]+part_dest:
        if os.path.lexists(str(p)):raise RuntimeError('Refusing existing destination '+str(p))
    m.notify.save_json(S/'status.json',status);m.notify.save_json(S/'process.json',{'pid':os.getpid(),'host':m.socket.getfqdn(),'startedAt':m.notify.now()})
    staging=directory/('.merge_staging_'+S.name);staging.mkdir();normal=staging/'normal.root'
    parts=[staging/('pid_pair_joint_part%03d.root'%i) for i in range(3)]
    runner=['bash',str(repo/'script/singularity_root_batch.sh'),str(repo/'config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml')]
    try:
        original=json.loads((old/'input_manifest.json').read_text())['inputs'];original_paths=[Path(x['path']) for x in original]
        if m.fingerprints(original_paths)!=original:raise RuntimeError('Original subjob inputs changed')
        m.notify.save_json(S/'input_manifest.json',{'inputs':original});m.write_list(S/'inputs.list',original_paths)
        sources=[];expected={};source_infos=[]
        for i in range(3):
            info=json.loads((old/('intermediate_%d/completed.json'%i)).read_text());source=Path(info['path']);m.archive_helper.verify_fingerprint(source,info['fingerprint']);sources.append(source);source_infos.append(info)
            v.accumulate(expected,v.summary(Path(info['actualAudit']),1))
        sources_before=m.fingerprints(sources);m.notify.save_json(S/'source_intermediates.json',{'inputs':sources_before})
        source_list=S/'intermediates.list';m.write_list(source_list,sources)
        plan={'schemaVersion':1,'jobid':meta['jobid'],'anaName':meta['anaName'],'normalRoot':str(dest),'normalHistogramCount':1786,'separateObject':'hNuclearPidPairJoint','separateClass':'THnSparseD','parts':[],'combination':'Project each sparse part using identical selections, then add projections; never persist the full 8D sum','formatApproval':str(S/'format_approval.json')}
        for i in range(3):plan['parts'].append({'index':i,'file':str(part_dest[i]),'source':str(sources[i]),'subjobCount':[10000,10000,6955][i]})
        m.notify.save_json(S/'output_plan.json',plan)
        status.update(phase='merging_normal_and_extracting_sparse',checkedAt=m.notify.now());m.notify.save_json(S/'status.json',status)
        # Independent outputs; two low-priority ROOT processes at most.
        with concurrent.futures.ThreadPoolExecutor(max_workers=2) as pool:
            futures=[pool.submit(invoke,'merge',source_list,normal,S/'output_plan.json',runner,S/'normal_merge.log')]
            futures += [pool.submit(invoke,'extract',sources[i],parts[i],S/('sparse_part%03d.tsv'%i),runner,S/('sparse_part%03d.log'%i)) for i in range(3)]
            for future in futures:future.result()
        if m.fingerprints(sources)!=sources_before:raise RuntimeError('Intermediate input changed')
        status.update(phase='validating_split_outputs',checkedAt=m.notify.now());m.notify.save_json(S/'status.json',status)
        sparse_records=[]
        for i in range(3):
            report=S/('sparse_part%03d.tsv'%i);rows=report.read_text().splitlines()
            if rows[-1]!='DONE\t1':raise RuntimeError('Incomplete sparse validation')
            values=dict(row.split('\t',1) for row in rows[:-1])
            if values['EXPLICIT_SUMW2']!='1' or values['DIMENSIONS']!='8':raise RuntimeError('Sparse schema/error array not preserved')
            sparse_records.append({'path':str(part_dest[i]),'size':parts[i].stat().st_size,'fingerprint':m.archive_helper.regular(parts[i]),'validation':'Exact bin-by-bin coordinates, contents, errors2, axes, Entries and moments vs source','report':str(report),'entries':float(values['ENTRIES']),'filledBins':int(values['FILLED_BINS'])})
        (S/'merge_notify_audit.C').write_text(m.ROOT_AUDIT);actual=m.audit([normal],'merged_audit',S,runner)
        if len(actual)!=1786:raise RuntimeError('Wrong normal histogram count '+str(len(actual)))
        policy=json.loads((S/'acceptance.json').read_text());v.validate(expected,actual,policy,S/'normal_vs_intermediates.json',sources,normal,runner)
        comparison=json.loads((prior/'event_comparison_and_authorization.json').read_text())
        counters={'inputEventsFromRefMultEntries':actual['hRefMult;1'][1][0],'reconstructedEventsFromRefMultCorrEntries':actual['hRefMultCorr;1'][1][0]}
        if counters['inputEventsFromRefMultEntries']!=comparison['currentInputEvents'] or counters['reconstructedEventsFromRefMultCorrEntries']!=comparison['currentReconstructedEvents']:raise RuntimeError('Merged event counters differ')
        m.notify.save_json(S/'event_counter_validation.json',counters)
        if m.fingerprints(original_paths)!=original:raise RuntimeError('Original input changed before publication')
        for i in range(3):m.publish_no_clobber(parts[i],part_dest[i])
        m.write_list(list_dest,part_dest)
        manifest=dict(plan,sparseValidation=sparse_records,eventCounters=counters,partial=True,inputCount=len(original),expectedSubjobs=status['expectedSubjobs'],excludedIndices=status['excludedIndices'],stateDirectory=str(S))
        with manifest_dest.open('x') as out:json.dump(manifest,out,ensure_ascii=False,indent=2);out.write('\n')
        with helper_dest.open('xb') as out:out.write((S/'project_joint_parts.C').read_bytes())
        m.publish_no_clobber(normal,dest)
        status.update(phase='merged',checkedAt=m.notify.now(),sparseParts=sparse_records,outputManifest=str(manifest_dest),sparsePartsList=str(list_dest),projectionHelper=str(helper_dest),merge={'path':str(dest),'size':dest.stat().st_size,'fingerprint':m.archive_helper.regular(dest),'histogramCount':len(actual),'validation':'ROOT recursive read; exact Entries and event counts; additive validation with seven user-approved exceptions and exact native float replay; sparse separately copied and compared in every bin','inputCount':len(original),'publishedAt':m.notify.now(),'stagingDirectory':str(staging)})
    except Exception as error:
        status.update(phase='failed',checkedAt=m.notify.now(),error=str(error));traceback.print_exc()
    m.notify.save_json(S/'result.json',status);m.notify.save_json(S/'status.json',status)
    m.archive_after_merge(status,runmeta,S)
    msg=m.make_message(config['email'],meta,status)
    good=status['phase']=='merged' and status.get('archive',{}).get('phase')=='archived'
    msg.replace_header('Subject','[STAR 部分hadd完了・Sparse分割保存]' if good else '[STAR hadd終了・要確認]')
    msg.replace_header('Message-ID','<star-'+meta['jobid']+'-sparse-split-20261026@'+m.socket.getfqdn()+'>')
    msg.set_content(msg.get_content()+'\n承認された保存方式: 通常1786ヒストを1本、巨大hNuclearPidPairJointを3本の別ROOTに保存。\nSparseファイル: '+'; '.join(x['path'] for x in status.get('sparseParts',[]))+'\n出力manifest: '+status.get('outputManifest','未作成')+'\n')
    (S/'message.eml').write_bytes(msg.as_bytes());m.notify.send_once(S,msg,status,config['email'])
    print(m.notify.now(),'Result mail accepted',status['phase'],flush=True)
    return 0 if good else 1
if __name__=='__main__':sys.exit(main())
