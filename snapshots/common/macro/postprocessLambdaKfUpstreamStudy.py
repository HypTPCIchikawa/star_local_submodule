#!/usr/bin/env python3
"""Independent ROOT selection replay / background checks / frozen final figure.
Creates new products only. Does not choose cuts or reconstruct PicoDst.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import subprocess

REPO=Path(__file__).resolve().parents[2]
WRAPPER=REPO/'rootfile/auau13p5_anaLambda_KFParticle/purity_study_20260927/provenance/helpers/run_root5_readonly.sh'

def run_qa(source, expression, log):
    log.parent.mkdir(parents=True,exist_ok=True)
    call='common/macro/runLambdaKfPurityQa.C({},{})'.format(json.dumps(source),json.dumps(expression))
    with log.open('x') as stream:
        result=subprocess.run(['bash',str(WRAPPER),call],cwd=REPO,stdout=stream,stderr=subprocess.STDOUT,timeout=1800)
    if result.returncode: raise RuntimeError('ROOT QA failed: '+str(log))
    print('PASS',log.name,flush=True)

def expression(function,*args):
    return function+'('+','.join(json.dumps(str(x)) for x in args)+')'

def with_source_baseline(folder, source):
    """Add only a diagnostic baseline to NEW source-specific check manifests.

    The old independent checker requires a selected==1-only selection in each
    invocation. A winning cut from another upstream source does not include
    that source's uncut selected spectrum. Use its own exported native spectrum,
    never the common default denominator, and leave all original scan files
    and their hashes unchanged.
    """
    cuts_path, hist_path=folder/source['cuts'],folder/source['histograms']
    with cuts_path.open(newline='') as stream:
        reader=csv.DictReader(stream,delimiter='\t')
        cut_fields=reader.fieldnames
        cuts=list(reader)
    with hist_path.open(newline='') as stream:
        reader=csv.DictReader(stream)
        hist_fields=reader.fieldnames
        histograms=list(reader)
    if cut_fields!=['selection','key','operator','value'] or hist_fields!=['selection','bin','low','high','count','error']:
        raise RuntimeError('Unexpected independent-check manifest schema')
    by_selection={}
    for row in cuts: by_selection.setdefault(row['selection'],[]).append(row)
    has_baseline=any(len(rows)==1 and rows[0]['key']=='selected' and
        rows[0]['operator'] in ('eq','==') and float(rows[0]['value'])==1.
        for rows in by_selection.values())
    if not has_baseline:
        label='source_baseline'
        if label in by_selection: raise RuntimeError('Diagnostic baseline label collision')
        provenance=json.loads((folder/'sources.json').read_text())
        records=[item for item in provenance['datasets'] if item['id']==source['dataset_id']]
        if len(records)!=1: raise RuntimeError('Missing or duplicate upstream source provenance')
        record=records[0]
        if Path(record['meta']['inputRoot']).resolve()!=Path(source['inputRoot']).resolve():
            raise RuntimeError('Source ROOT provenance mismatch')
        native_record=record['inputs']['histogram']
        native_path=Path(native_record['path'])
        digest=hashlib.sha256(native_path.read_bytes()).hexdigest()
        if digest!=native_record['sha256']: raise RuntimeError('Native source histogram hash mismatch')
        with native_path.open(newline='') as stream:
            reader=csv.DictReader(stream)
            if reader.fieldnames!=['bin','low','high','count','error']:
                raise RuntimeError('Unexpected native source histogram schema')
            native=list(reader)
        if [int(row['bin']) for row in native]!=list(range(202)):
            raise RuntimeError('Native source histogram requires all 202 bins')
        cuts.append(dict(selection=label,key='selected',operator='eq',value='1'))
        histograms.extend(dict(selection=label,**row) for row in native)
    derived_cuts=cuts_path.with_name(cuts_path.stem+'_with_source_baseline.tsv')
    derived_hist=hist_path.with_name(hist_path.stem+'_with_source_baseline.csv')
    if derived_cuts.exists() or derived_hist.exists():
        raise RuntimeError('Refusing existing derived source-baseline manifests')
    with derived_cuts.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=cut_fields,delimiter='\t',lineterminator='\n')
        writer.writeheader();writer.writerows(cuts)
    with derived_hist.open('x',newline='') as stream:
        writer=csv.DictWriter(stream,fieldnames=hist_fields,lineterminator='\n')
        writer.writeheader();writer.writerows(histograms)
    print('SOURCE_BASELINE',source['dataset_id'],'already_present' if has_baseline else 'added_native_source',flush=True)
    return derived_cuts,derived_hist

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--energy',choices=('auau13p5','auau3p85'),required=True)
    p.add_argument('--phase',choices=('scan','validation'),required=True)
    args=p.parse_args()
    out=REPO/'rootfile'/(args.energy+'_anaLambda_KFParticle')/'purity_upstream_20260928'
    logs=out/'provenance/logs'
    figure=REPO/'share/figure'/(args.energy+'_Lambda_KFParticle_purity_upstream_20260928')
    figure.mkdir(parents=True,exist_ok=True)
    suffix='_v2' if args.phase=='scan' else ''
    if args.phase=='scan':
        folder=out/'scan'
        with (folder/'verification_sources.csv').open() as stream:
            sources=list(csv.DictReader(stream))
        for source in sources:
            label=source['dataset_id']
            cuts,histograms=with_source_baseline(folder,source)
            root=out/'scan/verification'/(label+'_independent'+suffix+'.root')
            call=expression('check_lambda_purity_selection',source['inputRoot'],cuts,histograms,root)
            run_qa('tests/check_lambda_purity_selection.C',call,logs/('scan_independent_'+label+suffix+'.log'))
    else:
        folder=out/'validation100k/comparison'
        run_qa('common/macro/drawLambdaKfFrozenComparison.C',expression('drawLambdaKfFrozenComparison',folder,figure/'validation100k_default_vs_best'),logs/'validation100k_figure.log')
    run_qa('common/macro/checkLambdaKfPurityBackground.C',
        expression('checkLambdaKfPurityBackground',folder,figure/(args.phase+'_background_checks'+suffix)),
        logs/(args.phase+'_background_checks'+suffix+'.log'))

if __name__=='__main__': main()
