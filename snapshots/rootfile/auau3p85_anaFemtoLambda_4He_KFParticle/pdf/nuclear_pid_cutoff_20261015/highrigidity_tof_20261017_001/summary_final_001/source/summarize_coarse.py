#!/usr/bin/env python3
"""Independent coarse He4 summary; no splicing of successful native points."""
import argparse
import csv
import math
from pathlib import Path
import shutil
import run_highrigidity as high
from summarize_cutoff import fitted_response_failures
common=high.common


def run(args):
    native=common.resolve(args.native_summary);ndone=common.read(native/'completed.json')
    if not ndone.get('executionComplete') or not ndone.get('protectedInputsVerified'): raise ValueError('Incomplete native summary.')
    protected={}
    for name in ('summary.tsv','summary.json','context.json','native_parameters.tsv'):
        expected=ndone['products'][name]
        if common.digest(native/name)!=expected: raise ValueError('Changed native summary.')
        protected[str(native/name)]=expected
    protected[str(native/'completed.json')]=common.digest(native/'completed.json')
    models={};reference=None;parameters=[];inputs=[]
    for folder in args.coarse_runs:
        folder=common.resolve(folder);done=common.read(folder/'completed.json');context=common.read(folder/'context.json')
        if not done.get('executionComplete') or not done.get('protectedSourcesVerified') or context['stage']!='coarse_native':
            raise ValueError('Require completed full-domain coarse fits.')
        for name,expected in done['products'].items():
            if common.digest(folder/name)!=expected: raise ValueError('Changed coarse product.')
            protected[str(folder/name)]=expected
        protected[str(folder/'completed.json')]=common.digest(folder/'completed.json')
        if context['targetSpecies']!='he4' or context['sampleHistogram']!='TAG_D' or context['fitRigidityWidth']!=.05 or context['sourceNativeRigidityWidth']!=.01:
            raise ValueError('Wrong target/TOF/domain provenance.')
        key=tuple(context[k] for k in ('response','background','variation'))
        if key in models: raise ValueError('Repeated model.')
        settings=common.read(folder/'effective_settings.json')
        state=((folder/'gates.tsv').read_bytes(),(folder/'fixed_windows.json').read_bytes(),
            (folder/'effective_settings.json').read_bytes(),context['sourceOriginalInput'],context['mass2Window'],
            {k:v for k,v in done['products'].items() if k.startswith('source/')})
        if reference is None: reference=state
        elif reference!=state: raise ValueError('Coarse models differ in settings/gate/window/input.')
        rows=common.read(folder/'native_summary.json')['nativeRows'];fit=common.read(folder/'cutoff_fits.json')
        audit={r['index']:r for r in common.read(folder/'tag_shape_audit.json')}
        for saved,row in zip(rows,fit['slices']):
            metrics=high.aggregate_one([row],'he4',audit,dict(settings,fitWindowMode='target_local'))
            for k in ('unsupportedNuisanceFraction','omissionBoundOverRawGateCount','localResidualFrac','integralResidualFrac'):
                saved[k]=metrics[k]
        if len(rows)!=60 or [r['rBin'] for r in rows]!=list(range(41,101)): raise ValueError('Coarse domain has a hole.')
        models[key]={r['rBin']:r for r in rows};label='_'.join(key)
        for row in fit['slices']:
            c=row['components']['he4']
            parameters.append([label,'he4',row['lo'],row['hi'],row['r'],c['mean'],c['sigma'],c['meanError'],c['sigmaError'],
                c['initialMean'],c['initialSigma'],row['targetCenterShiftSigma'],row['targetWidthRatio'],
                int(row['targetFreeResponseSupported']),int(row['numericalValid']),int(not fitted_response_failures(row,'he4',settings)),'TAG_D',0,.05])
        inputs.append(str(folder))
    expected={(r,b,'nominal') for r in ('gaussian','lognormal') for b in ('linear','exponential')}
    if set(models)!=expected: raise ValueError('All four coarse models required.')
    coarse=[]
    for rbin in range(41,101):
        rows=[v[rbin] for v in models.values()];nominal=models['gaussian','linear','nominal'][rbin]
        purity=[r['purityNumerical'] for r in rows if r['purityNumerical'] is not None and math.isfinite(r['purityNumerical'])]
        if max(r['rawGateCount'] for r in rows)-min(r['rawGateCount'] for r in rows)>1.e-8*max(1.,nominal['rawGateCount']):
            raise ValueError('Fixed raw count differs among models.')
        coarse.append(dict(species='he4',rLo=nominal['rLo'],rHi=nominal['rHi'],purityNominal=nominal['purityNumerical'],
            purityLow=min(purity) if purity else None,purityHigh=max(purity) if purity else None,
            rawSelected=nominal['rawGateCount'],quality=all(r['quality'] for r in rows),
            status='supported' if all(r['quality'] for r in rows) else 'unassessed',
            reasons=sorted({reason for r in rows for reason in r['reasons']}),method='independent_0p05_R_projected_refit'))
    native_bins=common.read(native/'summary.json')['bins']
    final=[dict(r,method='sum_of_five_independent_0p01_R_fits') for r in native_bins if r['species']!='he4']+coarse
    fields='species rLo rHi purityNominal purityLow purityHigh rawSelected quality'
    parameter_fields='label species rLo rHi r mean sigma meanError sigmaError initialMean initialSigma centerShiftSigma widthRatio targetFreeSupported numericalValid fitResponseQuality sample exactProductionTofWindow fitRigidityWidth'
    with (native/'native_parameters.tsv').open() as stream:
        old_parameters=[r for r in csv.DictReader(stream,delimiter='\t') if r['species']!='he4']
    native_parameters=[[r[k] for k in parameter_fields.split()[:-1]]+[.01] for r in old_parameters]
    for output,rows,label in ((args.coarse_summary.resolve(),coarse,'independent_he4_coarse_refits'),
                              (args.final_summary.resolve(),final,'native_d_t_he3_and_independent_he4_coarse_refits')):
        if output.exists() or output.parent!=native.parent: raise ValueError('Require new direct summary sibling.')
        output.mkdir();source=output/'source';source.mkdir();shutil.copy2(Path(__file__),source/Path(__file__).name)
        common.tsv(output/'summary.tsv',fields,[[int(r[k]) if k=='quality' else r[k] if r[k] is not None else 'nan' for k in fields.split()] for r in rows])
        high.save(output/'summary.json',dict(schemaVersion=1,productionEnabled=False,fullUncertaintyValidated=False,
            stage=label,purityDefinition='Same-bin sum(S)/sum(total); He4 separately refitted on exact .05-R projected counts.',
            he4RigidityRange=[2.,5.],he4FitRWidth=.05,originalNativeRWidth=.01,he4ProjectedWidthIsIntrinsic=False,
            he4TofPolicy='User-requested same saved D window, distinct from pre-existing production 3<m2<4.',
            bins=rows,modelDetails={'_'.join(k):list(v.values()) for k,v in models.items()},
            nativeSummary=str(native),coarseRuns=inputs))
        common.tsv(output/'native_parameters.tsv',
            parameter_fields,[[v if v is not None else 'nan' for v in row] for row in
                (parameters if label=='independent_he4_coarse_refits' else native_parameters+parameters)])
        high.save(output/'context.json',dict(schemaVersion=1,protectedInputs=protected,sourceSha256=common.digest(Path(__file__))))
        for p,expected_sha in protected.items():
            if common.digest(p)!=expected_sha: raise ValueError('Input changed during coarse summary: '+p)
        high.save(output/'completed.json',dict(executionComplete=True,productionEnabled=False,protectedInputsVerified=True,
            fullUncertaintyValidated=False,rows=len(rows),supportedRows=sum(r['quality'] for r in rows),
            products={str(p.relative_to(output)):common.digest(p) for p in output.rglob('*') if p.is_file()}))
        print('COARSE_HIGH_R_SUMMARY_COMPLETE',output,'supported',sum(r['quality'] for r in rows),'/',len(rows),flush=True)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--native-summary',required=True,type=Path)
    parser.add_argument('--coarse-summary',required=True,type=Path)
    parser.add_argument('--final-summary',required=True,type=Path)
    parser.add_argument('--coarse-runs',required=True,nargs='+',type=Path)
    run(parser.parse_args())
