#!/usr/bin/env python3
"""Sum-S/sum-total four-model high-R purity, with fail-closed native coverage."""
import argparse
from collections import defaultdict
import math
from pathlib import Path
import shutil
import run_highrigidity as high
from summarize_cutoff import fitted_response_failures
common=high.common


def checked_native(rows,target,audit,cfg):
    values=[]
    for row in rows:
        v=high.aggregate_one([row],target,audit,dict(cfg,fitWindowMode='target_local'))
        drop=row['components'][target].get('conditionalDropDeltaD')
        if drop is None or drop<cfg['support']['minDropDeltaD']:
            v['reasons']=sorted(set(v['reasons']+['fitted_drop_delta_d_too_small']));v['supported']=False
        values.append(v)
    return values


def grouped_model(rows,target,audit,cfg):
    value=high.aggregate_one(rows,target,audit,dict(cfg,fitWindowMode='target_local'))
    if any(not v['supported'] for v in checked_native(rows,target,audit,cfg)):
        value['reasons']=sorted(set(value['reasons']+['one_or_more_native_fits_unassessed']));value['supported']=False
    if len(rows)!=round(cfg['aggregateWidth']/cfg['nativeWidth']):
        value['reasons']=sorted(set(value['reasons']+['incomplete_native_coverage']));value['supported']=False
    return value


def summarize(args):
    cfg=high.checked_settings(common.resolve(args.settings))
    expected={(r,b,v) for r in cfg['responseModels'] for b in cfg['backgroundModels'] for v in args.variations}
    if not set(args.variations)<=set(cfg['shapeVariations']): raise ValueError('Undeclared variation.')
    targets=defaultdict(dict);protected={}
    for folder in args.runs:
        folder=common.resolve(folder);done=common.read(folder/'completed.json');context=common.read(folder/'context.json')
        if not done.get('executionComplete') or not done.get('protectedSourcesVerified') or done.get('selfTestOnly'):
            raise ValueError('Incomplete fit: '+str(folder))
        for name in ('cutoff_fits.json','cutoff_fits.root','tag_shape_audit.json','settings.tsv','context.json','gates.tsv',
                     'shapes.tsv','raw_gate_counts.tsv','fixed_windows.json','current_selection.json','effective_settings.json'):
            actual=common.digest(folder/name)
            if actual!=done['products'][name]: raise ValueError('Changed completed fit artifact: '+str(folder/name))
            protected[str(folder/name)]=actual
        protected[str(folder/'completed.json')]=common.digest(folder/'completed.json')
        for name,expected_sha in done['products'].items():
            if name.startswith('source/'):
                actual=common.digest(folder/name)
                if actual!=expected_sha: raise ValueError('Changed completed source snapshot: '+str(folder/name))
                protected[str(folder/name)]=actual
        target=context['targetSpecies'];key=tuple(context[k] for k in ('response','background','variation'))
        if target not in cfg['targets'] or key not in expected or key in targets[target]: raise ValueError('Repeated/unknown model.')
        if context['sampleHistogram']!=cfg['targets'][target]['sample'] or context['exactProductionTofWindow']!=cfg['targets'][target]['exactProductionTofWindow']:
            raise ValueError('Contradictory sample/window provenance.')
        fit=common.read(folder/'cutoff_fits.json');audit={r['index']:r for r in common.read(folder/'tag_shape_audit.json')}
        local=common.read(folder/'effective_settings.json')
        if context['stage']!=args.stage or [r['rBin'] for r in fit['slices']]!=high.selected_bins(local,args.stage):
            raise ValueError('Incomplete/contradictory native coverage.')
        targets[target][key]=(folder,context,fit,audit,local)
    if set(targets)!=set(cfg['targets']) or any(set(m)!=expected for m in targets.values()):
        raise ValueError('All species and declared models/variations are required.')
    summaries=[];details=[];native=[];parameters=[]
    for target,models in targets.items():
        reference=None;groups=defaultdict(dict)
        for key,(folder,context,fit,audit,local) in models.items():
            source_versions={str(Path(p).relative_to(folder/'source')):sha for p,sha in protected.items()
                             if str(p).startswith(str(folder/'source')+'/')}
            state=((folder/'gates.tsv').read_bytes(),(folder/'fixed_windows.json').read_bytes(),
                   (folder/'effective_settings.json').read_bytes(),source_versions,
                   context['mass2Window'],context['sourceOriginalInput'])
            if reference is None: reference=state
            elif reference!=state: raise ValueError('Models must share TPC gate, fit window, TOF window and input.')
            label='_'.join(key);rows_by_group=defaultdict(list)
            for row,value in zip(fit['slices'],checked_native(fit['slices'],target,audit,local)):
                c=row['components'][target]
                native.append(dict(species=target,rLo=row['lo'],rHi=row['hi'],label=label,**value))
                parameters.append([label,target,row['lo'],row['hi'],row['r'],c['mean'],c['sigma'],c['meanError'],c['sigmaError'],
                    c['initialMean'],c['initialSigma'],row['targetCenterShiftSigma'],row['targetWidthRatio'],
                    int(row['targetFreeResponseSupported']),int(row['numericalValid']),int(not fitted_response_failures(row,target,local)),
                    context['sampleHistogram'],int(context['exactProductionTofWindow'])])
                group=round((row['lo']-cfg['rigidityMin'])/cfg['nativeWidth'])//round(cfg['aggregateWidth']/cfg['nativeWidth'])
                rows_by_group[group].append(row)
            for group,rows in rows_by_group.items(): groups[group][key]=grouped_model(rows,target,audit,local)
        for group,results in sorted(groups.items()):
            lo=cfg['rigidityMin']+group*cfg['aggregateWidth'];hi=lo+cfg['aggregateWidth']
            nominal=results.get(('gaussian','linear','nominal'))
            if nominal is None: raise ValueError('Nominal Gaussian+linear baseline missing.')
            numerical=[m['purity'] for m in results.values() if m['purity'] is not None and math.isfinite(m['purity'])]
            quality=set(results)==expected and all(m['supported'] for m in results.values())
            reasons=sorted({r for m in results.values() for r in m['reasons']})
            summaries.append(dict(species=target,rLo=lo,rHi=hi,purityNominal=nominal['purity'],
                purityLow=min(numerical) if numerical else None,purityHigh=max(numerical) if numerical else None,
                rawSelected=nominal['rawSelected'],quality=quality,modelCount=len(results),expectedModels=len(expected),
                status='supported' if quality else 'unassessed',reasons=reasons))
            details.append(dict(species=target,rLo=lo,rHi=hi,models={'_'.join(k):v for k,v in results.items()}))
    output=args.output_dir.resolve();allowed=(common.PROJECT/cfg['outputBase']).resolve()
    if output.exists() or output.parent!=allowed: raise ValueError('Require new summary child of outputBase.')
    output.mkdir();source=output/'source';source.mkdir();shutil.copy2(Path(__file__),source/Path(__file__).name)
    fields='species rLo rHi purityNominal purityLow purityHigh rawSelected quality'
    common.tsv(output/'summary.tsv',fields,[[int(r[k]) if k=='quality' else r[k] if r[k] is not None else 'nan' for k in fields.split()] for r in summaries])
    high.save(output/'summary.json',dict(schemaVersion=1,productionEnabled=False,fullUncertaintyValidated=False,
        stage=args.stage,expectedVariations=args.variations,expectedModelsPerSpecies=len(expected),
        purityDefinition='sum(S)/sum(total), never mean purity; every native fit/model required.',
        he4Study='Same saved D TOF window; current production 3<m2<4 is not evaluated.',bins=summaries,details=details))
    high.save(output/'native_models.json',native)
    common.tsv(output/'native_parameters.tsv',
        'label species rLo rHi r mean sigma meanError sigmaError initialMean initialSigma centerShiftSigma widthRatio targetFreeSupported numericalValid fitResponseQuality sample exactProductionTofWindow',
        [[v if v is not None else 'nan' for v in row] for row in parameters])
    high.save(output/'context.json',dict(schemaVersion=1,protectedInputs=protected,sourceSha256=common.digest(Path(__file__))))
    for path,expected_sha in protected.items():
        if common.digest(path)!=expected_sha: raise ValueError('Input changed during summary: '+path)
    high.save(output/'completed.json',dict(executionComplete=True,productionEnabled=False,protectedInputsVerified=True,
        fullUncertaintyValidated=False,rows=len(summaries),products={str(p.relative_to(output)):common.digest(p) for p in output.rglob('*') if p.is_file()}))
    print('NUCLEAR_HIGH_R_SUMMARY_COMPLETE',output,'supported',sum(r['quality'] for r in summaries),'/',len(summaries),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings',required=True,type=Path)
    parser.add_argument('--output-dir',required=True,type=Path)
    parser.add_argument('--stage',required=True,choices=('pilot','native'))
    parser.add_argument('--variations',required=True,nargs='+')
    parser.add_argument('--runs',required=True,nargs='+',type=Path)
    summarize(parser.parse_args())
