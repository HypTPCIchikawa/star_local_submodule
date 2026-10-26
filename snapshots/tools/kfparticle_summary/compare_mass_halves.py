#!/usr/bin/env python3
"""Mass-half CF comparison, using saved pair histograms without event reanalysis."""
import argparse
import datetime
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys

import yaml

from mass_half_statistics import (ErrorModel, normalized_cf, integrated_cf,
                                  purity_corrected_cf, contrast)

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
REGIONS = ('full', 'low', 'high')
METHODS = ('raw', 'sbsub', 'cfsub', 'method3', 'method5')


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def fingerprint(path):
    s = path.stat()
    return dict(path=str(path), size=s.st_size, mtime_ns=s.st_mtime_ns, inode=s.st_ino)


def command(args, logfile):
    print(datetime.datetime.now().isoformat(), logfile.name, flush=True)
    with logfile.open('w') as out:
        out.write('COMMAND ' + json.dumps(args) + '\n')
        out.flush()
        done = subprocess.run(args, cwd=str(PROJECT), stdout=out, stderr=subprocess.STDOUT, timeout=1800)
    if done.returncode:
        raise RuntimeError('Failed: ' + str(logfile))
    content = logfile.read_text(errors='replace')
    for word in ('Error in <ACLiC>', 'fatal error:', 'segmentation violation', 'Error: Function '):
        if word in content:
            raise RuntimeError(word + ': ' + str(logfile))


def selected_bins(edges, limits):
    for edge in limits:
        if not any(abs(edge - b) < 1.e-9 for b in edges):
            raise ValueError('Requested range is not aligned to saved/rebinned edges.')
    return [i for i in range(len(edges)-1) if limits[0] <= (edges[i]+edges[i+1])/2 < limits[1]]


def masked_cf(model, se, me, norm, region):
    if any(se[i] is None or me[i] is None for i in norm):
        return [None]*len(se), 'invalid mass fit in normalization interval: ' + region
    a = [q if q is not None else model.constant(0) for q in se]
    b = [q if q is not None else model.constant(0) for q in me]
    result = normalized_cf(a, b, norm)
    return [q if se[i] is not None and me[i] is not None else None
            for i, q in enumerate(result['bins'])], result['invalid']


def masked_integrated(se, me, norm, selected):
    if any(se[i] is None or me[i] is None for i in set(norm+selected)):
        return None
    # Unselected invalid bins do not enter either integral.
    model = next(q.model for q in se if q is not None)
    return integrated_cf([q if q is not None else model.constant(0) for q in se],
                         [q if q is not None else model.constant(0) for q in me], norm, selected)


def purity_proxies(data):
    inclusive = data['inclusive']
    bounds = data['massPartition']['binned']
    model = ErrorModel()
    q = {r: model.variable(r, inclusive[r]['count'], inclusive[r]['error2'])
         for r in ('low', 'high', 'leftSB', 'rightSB')}
    q['full'] = q['low'] + q['high'] + (inclusive['full']['count']-q['low'].value-q['high'].value)
    side_width = sum(bounds[r]['upper']-bounds[r]['lower'] for r in ('leftSB', 'rightSB'))
    result = {}
    for region in REGIONS:
        alpha = (bounds[region]['upper']-bounds[region]['lower'])/side_width
        p = 1.-alpha*(q['leftSB']+q['rightSB'])/q[region]
        if not 0 < p.value <= 1:
            raise ValueError('Nonphysical inclusive purity proxy: ' + region)
        result[region] = dict(model.estimate(p), width_scale=alpha)
    return result


def raw_inputs(model, group, settings):
    result = {}
    closures = {}
    for sample in ('SE', 'ME'):
        raw = group['raw'][sample]
        halves = {region: model.series(sample+'_'+region, raw[key]['counts'], raw[key]['error2'])
                  for region, key in [('low','signalLow'), ('high','signalHigh')]}
        closures[sample] = []
        full = []
        for i, (low, high, value, variance) in enumerate(zip(halves['low'], halves['high'],
                raw['signal']['counts'], raw['signal']['error2'])):
            delta = value-low.value-high.value
            vdelta = variance-model.variance(low+high)
            tolerance = settings['closure_absolute_tolerance'] + settings['closure_relative_tolerance']*abs(value)
            if abs(delta) > tolerance or abs(vdelta) > tolerance:
                raise ValueError('Non-rounding closure discrepancy: %s %s %d' % (group['id'], sample, i))
            full.append(low+high+delta)
            closures[sample].append(dict(bin=i, content_residual=delta, variance_residual=vdelta))
        result[sample] = dict(halves, full=full)
        left = model.series(sample+'_leftSB', raw['leftSB']['counts'], raw['leftSB']['error2'])
        right = model.series(sample+'_rightSB', raw['rightSB']['counts'], raw['rightSB']['error2'])
        result[sample]['sideband'] = [l+r for l,r in zip(left,right)]
    return result, closures


def coarse_fit_inputs(model, group, settings):
    """Apply shared 50-MeV purity to original 5-MeV counts, then sum to CF bins."""
    output, reasons = {}, []
    cf_rebin, purity_rebin = settings['rebin'], settings['purityRebin']
    for sample in ('SE', 'ME'):
        base = group['method5Base'][sample]
        cells = [[model.variable(('method5',sample,'base',i,r), row['counts'][r], row['error2'][r])
                  for r in range(2)] for i,row in enumerate(base)]
        factors = []
        for i,row in enumerate(group['methods']['method5'][sample]):
            factor = {r: None for r in REGIONS}
            background, covariance = row['background'], row['backgroundCovariance']
            for r,valid in enumerate(row['valid']):
                if not valid:
                    reasons.append(dict(sample=sample,purity_bin=i+1,region=('low','high')[r],
                                        reason=row['regionReason'][r]))
            if background is None or covariance is None or any(v is None for v in background):
                factors.append(factor)
                continue
            try:
                b = model.correlated_block(('method5',sample,'background',i),background,covariance)
            except ValueError as exc:
                reasons.append(dict(sample=sample,purity_bin=i+1,reason='invalid fit covariance: '+str(exc)))
                factors.append(factor)
                continue
            first, last = row['windowOriginalFirstBin']-1, row['windowOriginalLastBin']
            window = [sum((cells[j][r] for j in range(first,last)),model.constant(0)) for r in range(2)]
            for r,name in enumerate(('low','high')):
                if row['valid'][r] and window[r].value>0:
                    factor[name] = 1-b[r]/window[r]
            total = window[0]+window[1]
            if row['full']['valid'] and total.value>0:
                factor['full'] = 1-(b[0]+b[1])/total
            factors.append(factor)
        if len(cells)%cf_rebin or len(cells)%purity_rebin:
            raise ValueError('CF/purity binning does not divide original bins.')
        output[sample] = {r: [] for r in REGIONS}
        for first in range(0,len(cells),cf_rebin):
            for region in REGIONS:
                value = model.constant(0)
                for i in range(first,first+cf_rebin):
                    factor = factors[i//purity_rebin][region]
                    if factor is None:
                        value = None
                        break
                    count = cells[i][0]+cells[i][1] if region=='full' else cells[i][('low','high').index(region)]
                    value = value+count*factor
                output[sample][region].append(value)
    return output, reasons


def fit_inputs(model, group, method, settings):
    if method == 'method5':
        return coarse_fit_inputs(model,group,settings)
    result = {}
    reasons = []
    for sample in ('SE', 'ME'):
        rows = group['methods'][method][sample]
        central = [[model.variable((method,sample,'count',i,r), row['centralNtot'][r], row['centralError2'][r])
                    for r in range(2)] for i,row in enumerate(rows)]
        output = {r: [] for r in REGIONS}
        for i, row in enumerate(rows):
            for r,valid in enumerate(row['valid']):
                if not valid:
                    reasons.append(dict(sample=sample,bin=i,region=('low','high')[r],
                                        reason=row['regionReason'][r]))
            background, covariance = row['background'], row['backgroundCovariance']
            if background is None or covariance is None or any(v is None for v in background):
                for r in REGIONS: output[r].append(None)
                reasons.append(dict(sample=sample, bin=i, reason='fit unavailable'))
                continue
            try:
                b = model.correlated_block((method,sample,'background',i), background, covariance)
            except ValueError as exc:
                for r in REGIONS: output[r].append(None)
                reasons.append(dict(sample=sample, bin=i, reason='invalid fit covariance: '+str(exc)))
                continue
            if method == 'method3':
                fitted = [central[i][r]-b[r] for r in range(2)]
                full = fitted[0]+fitted[1]
            else:
                first, last = row['windowFirstBin']-1, row['windowLastBin']
                window = []
                for r in range(2):
                    w = sum((central[j][r] for j in range(first,last)), model.constant(0))
                    # Preserve the projected value if ROOT's TH2F summation differs slightly.
                    window.append(w + (row['windowNtot'][r]-w.value))
                fitted = [(central[i][r]*(1-b[r]/window[r]) if window[r].value>0 else None) for r in range(2)]
                total_window = window[0]+window[1]
                full = ((central[i][0]+central[i][1])*(1-(b[0]+b[1])/total_window)
                        if total_window.value>0 else None)
            for r, index in [('low',0),('high',1)]:
                output[r].append(fitted[index] if row['valid'][index] else None)
            output['full'].append(full if row['full']['valid'] else None)
        result[sample] = output
    return result, reasons


def analyze(data, settings, base):
    purity = purity_proxies(data)
    output = dict(schemaVersion=1, inputRoot=data['inputRoot'], purity=purity,
                  normalization=base['normalization'], comparison_range=settings['comparison_range'], groups=[])
    rows = ['method\tgroup\tlabel\tregion\tbin\tlow\thigh\tvalue\terror\tratio\tratio_error']
    labels = {g['id']: g['label'] for g in settings['groups']}
    for group in data['groups']:
        edges = group['kstarEdges']
        norm = selected_bins(edges, base['normalization'])
        selected = selected_bins(edges, settings['comparison_range'])
        displayed = selected_bins(edges, [0., settings['plot_x_max']])
        raw_model = ErrorModel()
        raw, closures = raw_inputs(raw_model, group, settings)
        group_result = dict(id=group['id'], label=labels[group['id']], closure=closures, methods={})
        for method in METHODS:
            model = raw_model if method in ('raw','sbsub','cfsub') else ErrorModel()
            reasons = []
            if method in ('method3','method5'):
                source, reasons = fit_inputs(model, group, method, data['settings'])
            elif method == 'sbsub':
                outer = data['settings']['sidebandOuterFactor']
                source = {s: {r: [a-(1 if r=='full' else .5)/(outer-1)*b for a,b in
                                  zip(raw[s][r],raw[s]['sideband'])] for r in REGIONS} for s in ('SE','ME')}
            else:
                source = raw
            cfs, integrated, invalid = {}, {}, {}
            if method == 'cfsub':
                bg_cf, bg_invalid = masked_cf(model, raw['SE']['sideband'], raw['ME']['sideband'], norm, 'background')
                bg_integrated = masked_integrated(raw['SE']['sideband'],raw['ME']['sideband'],norm,selected)
            for region in REGIONS:
                cfs[region], invalid[region] = masked_cf(model,source['SE'][region],source['ME'][region],norm,region)
                integrated[region] = masked_integrated(source['SE'][region],source['ME'][region],norm,selected)
                if method == 'cfsub':
                    invalid[region] = dict(signal=invalid[region], background=bg_invalid)
                    p = purity[region]['value']
                    cfs[region] = [purity_corrected_cf(s,b,p) if s is not None and b is not None else None
                                   for s,b in zip(cfs[region],bg_cf)]
                    integrated[region] = (purity_corrected_cf(integrated[region],bg_integrated,p)
                        if integrated[region] is not None and bg_integrated is not None else None)
            metrics = {}
            for region in REGIONS:
                numerator, full = integrated[region], integrated['full']
                metrics[region] = dict(cf=model.estimate(numerator) if numerator is not None else None,
                                      vsFull=contrast(model,numerator,full) if numerator is not None and full is not None else None)
            metrics['high_vs_low'] = (contrast(model,integrated['high'],integrated['low'])
                if integrated['high'] is not None and integrated['low'] is not None else None)
            serialized = {r: [] for r in REGIONS}
            for region in REGIONS:
                for i, q in enumerate(cfs[region]):
                    full = cfs['full'][i]
                    estimate = model.estimate(q) if q is not None else None
                    comparison = contrast(model,q,full) if q is not None and full is not None else None
                    ratio = comparison['ratio'] if comparison else None
                    serialized[region].append(dict(cf=estimate,vsFull=comparison))
                    values = [estimate['value'] if estimate else float('nan'), estimate['error'] if estimate else float('nan'),
                              ratio['value'] if ratio else float('nan'), ratio['error'] if ratio else float('nan')]
                    rows.append('\t'.join([method,group['id'],labels[group['id']],region,str(i+1),
                        format(edges[i],'.17g'),format(edges[i+1],'.17g')]+[format(v,'.17g') for v in values]))
            group_result['methods'][method] = dict(curves=serialized,integrated=metrics,
                invalid=invalid, fit_problems=reasons, valid_bins={r:sum(q is not None for q in cfs[r]) for r in REGIONS},
                valid_plot_bins={r:sum(cfs[r][i] is not None for i in displayed) for r in REGIONS},
                total_plot_bins=len(displayed))
        output['groups'].append(group_result)
    return output, '\n'.join(rows)+'\n'


def settings_text(maker, settings, base):
    fit = base['mass_fit']
    records = [('gate',[maker[k] for k in ('lambdaSignalMean','lambdaSignalSigma','lambdaSignalNSigma','lambdaSidebandOuterFactor')]),
        ('rebin',[base['rebin_group']]),('slidingHalfWidth',[base['sliding_half_width']]),
        ('purityRebin',[base['purity_rebin']]),
        ('normalization',base['normalization']),('fitRange',fit['range']),('fitMeanLimits',fit['mean_limits']),
        ('fitSigmaLimits',fit['sigma_limits']),('fitSigmaInitial',[fit['sigma_initial']]),
        ('minimumEntries',[fit['minimum_entries']]),('minimumCovarianceStatus',[fit['minimum_covariance_status']])]
    text = ''.join(key+' '+' '.join(map(str,values))+'\n' for key,values in records)
    return text+''.join('group %s %d %d\n' % (g['id'],g['min'],g['max']) for g in settings['groups'])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mainconf',required=True,type=Path)
    parser.add_argument('--input-root',required=True,type=Path)
    parser.add_argument('--energy',required=True,choices=['13p5','3p85'])
    parser.add_argument('--settings',type=Path,default=HERE/'mass_half_settings.yaml')
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--resume',action='store_true')
    args = parser.parse_args()
    source = args.input_root.resolve(strict=True)
    mainconf = args.mainconf.resolve(strict=True)
    config = yaml.safe_load(mainconf.read_text())
    maker_path = (PROJECT/'config'/config['maker']).resolve(strict=True)
    maker = yaml.safe_load(maker_path.read_text())
    settings_path = args.settings.resolve(strict=True)
    settings = yaml.safe_load(settings_path.read_text())
    base_path = (settings_path.parent/settings['summary_settings']).resolve(strict=True)
    base = yaml.safe_load(base_path.read_text())
    out = (args.output_dir or source.parent/'pdf'/'mass_halves').resolve()
    context = dict(input=fingerprint(source),energy=args.energy,config_hashes={str(p):digest(p) for p in
                   (mainconf,maker_path,settings_path,base_path)},mainconf=str(mainconf),settings=settings,base=base)
    context_path = out/'context.json'
    if context_path.exists():
        if not args.resume or json.loads(context_path.read_text())!=context:
            raise ValueError('Existing output/config differs; use a new directory or --resume for same inputs.')
    else:
        if out.exists() and any(out.iterdir()): raise ValueError('Refusing nonempty unowned directory.')
        out.mkdir(parents=True,exist_ok=True)
        save(context_path,context)
    (out/'logs').mkdir(exist_ok=True)
    (out/'macros').mkdir(exist_ok=True)
    explicit_settings = out/'extract_settings.txt'
    explicit_settings.write_text(settings_text(maker,settings,base))
    root = ['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q']
    extract_source = HERE/'ExtractMassHalfCF.C'
    extracted = out/'mass_half_extract.json'
    if not extracted.exists():
        copy = out/'macros'/extract_source.name
        shutil.copy2(str(extract_source),str(copy))
        call = str(copy)+'+('+','.join(json.dumps(str(p)) for p in (source,out,explicit_settings))+')'
        command(root+[call],out/'logs/extract.log')
    data = json.loads(extracted.read_text())
    if data['inputRoot']!=str(source): raise ValueError('Extraction input mismatch.')
    results, rows = analyze(data,settings,base)
    save(out/'comparison.json',results)
    (out/'cf_points.tsv').write_text(rows)
    plot_source = HERE/'DrawMassHalfCF.C'
    if not (out/'mass_halves_method5.pdf').exists():
        copy = out/'macros'/plot_source.name
        shutil.copy2(str(plot_source),str(copy))
        call = str(copy)+'+('+','.join([json.dumps(str(out/'cf_points.tsv')),json.dumps(str(out)),
                    json.dumps('Au+Au '+('13.5' if args.energy=='13p5' else '3.85')+' GeV, KFParticle')]+[
                    repr(v) for v in base['normalization']+[settings['plot_x_max']]])+')'
        command(root+[call],out/'logs/plot.log')
    command([sys.executable,str(HERE/'report_mass_halves.py'),'--output-dir',str(out)],out/'logs/report.log')
    if fingerprint(source)!=context['input']: raise RuntimeError('Input ROOT changed during comparison.')
    save(out/'completed.json',dict(phase='complete',input_unchanged=True,
        sources={str(p):digest(p) for p in (Path(__file__),extract_source,plot_source,HERE/'mass_half_statistics.py',HERE/'report_mass_halves.py')},
        extract_sha256=digest(extracted),comparison_sha256=digest(out/'comparison.json')))
    print('COMPLETE',out,flush=True)


if __name__=='__main__':
    main()
