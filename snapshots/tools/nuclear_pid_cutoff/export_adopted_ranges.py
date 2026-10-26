#!/usr/bin/env python3
"""Export supported low-TPC ranges as a lower cut plus explicit veto intervals."""
import argparse
import csv
from decimal import Decimal
from pathlib import Path
import shutil
import yaml
from run_cutoff import read, save, resolve, digest, tsv, SPECIES


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--summary', required=True, type=Path)
    parser.add_argument('--raw-gates', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    summary_path, raw_path = resolve(args.summary), resolve(args.raw_gates)
    summary, complete = read(summary_path), read(summary_path.parent/'completed.json')
    if not complete['executionComplete'] or digest(summary_path)!=complete['products'][summary_path.name]:
        raise ValueError('Incomplete or changed summary.')
    cfg = yaml.safe_load((summary_path.parent/'settings.yaml').read_text())
    raw_complete = read(raw_path.parent/'completed.json')
    if not raw_complete['executionComplete'] or digest(raw_path)!=raw_complete['products'][raw_path.name]:
        raise ValueError('Independent source counts changed.')
    if not summary['declaredShapeSensitivityComplete']:
        raise ValueError('Final adopted ranges require the complete declared sensitivity comparison.')
    output = args.output_dir.resolve()
    if output.exists() or output.parent!=resolve(cfg['outputBase']):
        raise ValueError('Choose a new result directory below outputBase.')
    with raw_path.open() as stream:
        raw = list(csv.DictReader(stream,delimiter='\t'))
    digits = max(0,-Decimal(str(cfg['nativeWidth'])).as_tuple().exponent)
    edge = lambda x: round(x,digits)
    tolerance = cfg['numeric']['axisTolerance']
    rows, selections, result = [], [], {}
    for species in SPECIES:
        ranges = [[edge(a),edge(b)] for a,b in summary['recommendations'][species]['safeIntervals']]
        veto = []
        minimum = ranges[0][0] if ranges else None
        if ranges:
            cursor = minimum
            for lo,hi in ranges:
                if lo>cursor+tolerance:
                    veto.append(dict(min=edge(cursor),max=edge(lo)))
                cursor = hi
            if cursor<cfg['rigidityMax']-tolerance:
                veto.append(dict(min=edge(cursor),max=edge(cfg['rigidityMax'])))
        inside = lambda lo,hi:any(lo>=a-tolerance and hi<=b+tolerance for a,b in ranges)
        observed = [r for r in raw if r['species']==species]
        denominator = sum(float(r['rawGateCount']) for r in observed)
        numerator = sum(float(r['rawGateCount']) for r in observed if inside(float(r['lo']),float(r['hi'])))
        accepted = [r for r in summary['aggregateRows'] if r['species']==species and inside(r['rLo'],r['rHi'])]
        if any(not r['quality'] for r in accepted):
            raise ValueError('A range bridged a failed-quality point.')
        lower = min((r['purityLow'] for r in accepted),default=None)
        value = dict(adoptedMin=minimum,adoptedVeto=veto,adoptedLowTpcRanges=ranges,
            baselineRawCount=denominator,adoptedRawCount=numerator,
            rawTrackRetention=numerator/denominator if denominator>0 else None,
            minModelPurityIncludingOmittedBound=lower,
            singleLowerCutSufficientWithoutVeto=summary['recommendations'][species]['singleLowerCutSufficient'])
        result[species] = value
        rows.extend([species,a,b] for a,b in ranges)
        selections.append([species,minimum,numerator,denominator,value['rawTrackRetention'],lower,
                           ';'.join('['+str(v['min'])+','+str(v['max'])+')' for v in veto) or 'none'])
    output.mkdir()
    tsv(output/'adopted_ranges.tsv','species min max',rows)
    tsv(output/'adopted_selection.tsv',
        'species adoptedMin adoptedRaw baselineRaw rawRetention minModelPurity adoptedVeto',selections)
    save(output/'adopted_ranges.json',dict(schemaVersion=1,productionEnabled=False,
        studyRecommendationOnly=True,species=result,
        tpcOnlyRange=[cfg['rigidityMin'],cfg['rigidityMax']],
        policy='Min plus [min,max) veto, not a single unconditional lower cut. '
               'R>=TPC endpoint requires separately configured TOF; no high-R range is inferred here. '
               'Raw retention is inclusive track-count retention, not true efficiency or CF pair retention.'))
    snapshot = output/'source'
    snapshot.mkdir()
    for source in (Path(__file__).resolve(),Path(__file__).resolve().with_name('run_cutoff.py'),
                   summary_path.parent/'settings.yaml'):
        shutil.copy2(source,snapshot/source.name)
    save(output/'context.json',dict(sourceSummary=str(summary_path),summarySha256=digest(summary_path),
        rawCounts=str(raw_path),rawCountsSha256=digest(raw_path),
        sourceSummaryCompletionSha256=digest(summary_path.parent/'completed.json')))
    products = {str(p.relative_to(output)):digest(p) for p in output.rglob('*') if p.is_file()}
    save(output/'completed.json',dict(executionComplete=True,productionEnabled=False,
                                    products=products,derivedOnly=True))
    print('NUCLEAR_ADOPTED_RANGES_COMPLETE',output,flush=True)


if __name__=='__main__':
    main()
