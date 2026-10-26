#!/usr/bin/env python3
"""Compare fixed TPC gates using sum(S)/sum(N), never averages of purity."""
import argparse
from collections import defaultdict
import csv
import math
from pathlib import Path
import shutil
import yaml
from run_cutoff import PROJECT, SPECIES, read, save, resolve, digest, module, tsv


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def models_supported(models):
    # Compare only this target's models, not all runs of the other three targets.
    return bool(models) and all(finite(m['purity']) and m['supported'] for m in models.values())


def aggregate_one(rows, species, audit, cfg, raw_counts=None):
    selected = [r['components'][species]['selection'] for r in rows]
    reasons, signal, total, raw, absolute, unsupported, omitted = [], 0., 0., 0., 0., 0., 0.
    conditional_errors = []
    for row, item in zip(rows, selected):
        for key in ('numericalValid', 'available', 'selectionCovered', 'axisCovered', 'gateActive', 'identitySupport'):
            if not item.get(key):
                reasons.append(key)
        if item.get('boundaryTarget'):
            reasons.append('target_yield_boundary')
        if cfg['allFit'].get('targetShapeMode') == 'free_target_mean_sigma' and \
                not item.get('targetFreeResponseSupported', False):
            reasons.append('free_target_shape_or_identity_screen')
        if not finite(item.get('purity')):
            reasons.append('purity_unavailable')
        s, n = item.get('signalSelected'), item.get('totalSelected')
        if finite(s) and finite(n):
            signal += s
            total += n
            if finite(item.get('conditionalError')):
                conditional_errors.append(n*item['conditionalError'])
        else:
            reasons.append('selected_counts_unavailable')
        measured = raw_counts.get((row['rBin'],species)) if raw_counts is not None else None
        if measured is not None:
            count = measured['rawGateCount']
            if item.get('gateIntervalsAtNodes') and abs(count-item['rawGateCount']) > \
                    cfg['numeric']['countVarianceRelativeTolerance']*max(1.,count):
                raise ValueError('Model-independent and fit-stage raw gate counts disagree.')
            raw += count
        elif finite(item.get('rawGateCount')):
            raw += item['rawGateCount']
        else:
            reasons.append('raw_gate_unavailable')
        if finite(item.get('gateAbsoluteResidualCount')):
            absolute += item['gateAbsoluteResidualCount']
        if 'omittedLocalYieldUpperBound' in item and finite(item['omittedLocalYieldUpperBound']):
            omitted += item['omittedLocalYieldUpperBound']
        elif cfg.get('fitWindowMode') == 'target_local':
            reasons.append('omitted_nuisance_bound_unavailable')
        counts = item.get('selectedColumnCounts', [])
        for i, other in enumerate(SPECIES):
            if other != species and i < len(counts) and audit[row['index']]['failures'][other]:
                unsupported += counts[i]
    purity = signal/total if total > 0 else None
    local = absolute/raw if raw > 0 else None
    integral = abs(total-raw)/raw if raw > 0 else None
    nuisance = unsupported/total if total > 0 else None
    omission_fraction = omitted/raw if raw > 0 else None
    if local is None or local > cfg['cutoff']['maxGateAbsoluteResidualFraction']:
        reasons.append('local_gate_absolute_residual')
    if integral is None or integral > cfg['cutoff']['maxGateIntegralResidualFraction']:
        reasons.append('gate_integral_residual')
    if nuisance is None or nuisance > cfg['cutoff']['maxUnsupportedContaminationFraction']:
        reasons.append('unsupported_material_nuisance')
    if omission_fraction is None or omission_fraction > cfg['cutoff']['maxUnsupportedContaminationFraction']:
        reasons.append('omitted_nuisance_bound_material')
    return dict(purity=purity, signalSelected=signal, totalSelected=total, rawSelected=raw,
        localResidualFrac=local, integralResidualFrac=integral, unsupportedNuisanceFraction=nuisance,
        omittedLocalYieldUpperBound=omitted, omissionBoundOverRawGateCount=omission_fraction,
        purityLowerIncludingOmittedBound=signal/(total+omitted) if total+omitted > 0 else None,
        conditionalErrorWeightedBound=sum(conditional_errors)/total if total > 0 and len(conditional_errors)==len(rows) else None,
        reasons=sorted(set(reasons)), supported=not reasons, nativeRows=len(rows))


def safe_intervals(rows, cfg):
    result, start, stop, count = [], None, None, 0
    for row in rows+[None]:
        if row is not None and row['quality']:
            if start is None:
                start = row['rLo']
            stop, count = row['rHi'], count+1
        elif start is not None:
            if count >= cfg['cutoff']['minContiguousBins']:
                result.append([start, stop])
            start, stop, count = None, None, 0
    return result


def recommendation(rows, cfg):
    intervals = safe_intervals(rows, cfg)
    suffix = [r for r in intervals if abs(r[1]-cfg['rigidityMax']) <= cfg['numeric']['axisTolerance']]
    cutoff = suffix[0][0] if suffix else None
    return dict(safeIntervals=intervals, lowerCutoff=cutoff,
        singleLowerCutSufficient=cutoff is not None,
        firstPassingEdge=next((r['rLo'] for r in rows if r['quality']), None),
        semantics='Only a supported contiguous suffix reaching the full TPC-only endpoint defines a lower cut.')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings', required=True, type=Path)
    parser.add_argument('--baseline', required=True, type=Path)
    parser.add_argument('--runs', required=True, type=Path, nargs='+')
    parser.add_argument('--raw-gates', type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    args = parser.parse_args()
    cfg = yaml.safe_load(resolve(args.settings).read_text())
    raw_counts, raw_provenance = None, None
    if args.raw_gates:
        raw_path = resolve(args.raw_gates)
        raw_completion = read(raw_path.parent/'completed.json')
        if not raw_completion.get('executionComplete') or \
                digest(raw_path) != raw_completion['products'][raw_path.name]:
            raise ValueError('Independent raw-gate extraction is incomplete or changed.')
        raw_counts = {}
        with raw_path.open() as stream:
            for row in csv.DictReader(stream,delimiter='\t'):
                key = int(row['rBin']),row['species']
                if key in raw_counts:
                    raise ValueError('Duplicate neutral gate count.')
                value = {name:float(row[name]) for name in ('rawGateCount',
                    'rawGateCountLowerFullBins','rawGateCountUpperTouchingBins')}
                if not all(finite(x) and x>=0 for x in value.values()) or not \
                        value['rawGateCountLowerFullBins']<=value['rawGateCount']<=value['rawGateCountUpperTouchingBins']:
                    raise ValueError('Invalid independent raw-gate count/bounds.')
                raw_counts[key] = value
        raw_provenance = dict(path=str(raw_path),sha256=digest(raw_path),
            completionSha256=digest(raw_path.parent/'completed.json'),
            gatesSha256=digest(raw_path.parent/'gates.tsv'))
    elif cfg['allFit'].get('targetShapeMode')=='free_target_mean_sigma':
        raise ValueError('Free-stage summary requires fit-independent --raw-gates counts.')
    baseline_path = resolve(args.baseline)
    baseline = read(baseline_path)
    output = args.output_dir.resolve()
    if output.exists() or output.parent != resolve(cfg['outputBase']):
        raise ValueError('Choose a new summary sibling output directory.')
    products, runs, native_signatures, nominal = {}, [], {}, defaultdict(set)
    provenance = []
    for folder in map(resolve, args.runs):
        complete = read(folder/'completed.json')
        if not complete.get('executionComplete'):
            raise ValueError('Incomplete study cannot enter the envelope.')
        critical = ['cutoff_fits.json','cutoff_fits.root','tag_shape_audit.json','gates.tsv',
                    'fixed_windows.json','shapes.tsv','settings.tsv','context.json',
                    'effective_settings.json','source/NuclearFitCommon.h',
                    'source/NuclearLognormalFitCommon.h']
        critical += [n for n in complete['products'] if n.startswith('source/') and n.endswith('.C')]
        for name in critical:
            if digest(folder/name) != complete['products'][name]:
                raise ValueError('Changed run product: '+str(folder/name))
        context, result = read(folder/'context.json'), read(folder/'cutoff_fits.json')
        if digest(baseline_path) != context['protectedSources'][str(baseline_path)]:
            raise ValueError('Baseline changed after the fit.')
        provenance.append(dict(folder=str(folder), completedSha256=digest(folder/'completed.json'),
            criticalSha256={name:complete['products'][name] for name in critical}))
        target = context.get('targetSpecies') or 'ALL'
        signatures = digest(folder/'gates.tsv'), digest(folder/'fixed_windows.json')
        if raw_provenance and signatures[0] != raw_provenance['gatesSha256']:
            raise ValueError('Independent counter used a different TPC selector.')
        if target in native_signatures and native_signatures[target] != signatures:
            raise ValueError('A same-target comparison changed its fixed gate or fit windows.')
        native_signatures[target] = signatures
        if context['baseline'] != str(baseline_path):
            raise ValueError('Mixed runtime-gate baselines.')
        label = '_'.join(context[k] for k in ('response', 'background', 'variation'))
        product_label = target+'_'+label
        if product_label in products:
            raise ValueError('Duplicate response/background/variation run.')
        products[product_label] = str(folder)
        if context['variation'] == 'nominal':
            nominal[target].add((context['response'], context['background']))
        runs.append((label, context, result, {r['index']:r for r in read(folder/'tag_shape_audit.json')}))
    expected = {(a,b) for a in ('gaussian','lognormal') for b in ('linear','exponential')}
    if any(models != expected for models in nominal.values()) or \
            set(nominal) not in ({'ALL'}, set(SPECIES)):
        raise ValueError('All four nominal comparisons are required.')
    output.mkdir()
    snapshot = output/'source'
    snapshot.mkdir()
    source_signatures = {}
    for path in (Path(__file__).resolve(),Path(__file__).resolve().with_name('run_cutoff.py'),
                 resolve(args.settings),baseline_path):
        signature = digest(path)
        target = snapshot/path.name
        shutil.copy2(path,target)
        if digest(target) != signature:
            raise ValueError('Changed source while saving summary provenance.')
        source_signatures[target.name] = signature
    detailed, combined, table, native_parameters = [], [], [], []
    ratio = round(cfg['aggregateWidth']/cfg['nativeWidth'])
    if abs(ratio*cfg['nativeWidth']-cfg['aggregateWidth']) > cfg['numeric']['axisTolerance']:
        raise ValueError('Aggregation is not an integer native-bin grouping.')
    groups = defaultdict(dict)
    for label, context, result, audit in runs:
        by_group = defaultdict(list)
        for row in result['slices']:
            if raw_counts is not None and any((row['rBin'],s) not in raw_counts for s in SPECIES):
                raise ValueError('Neutral counter is missing a native R/species row.')
            index = round((row['lo']-cfg['rigidityMin'])/cfg['nativeWidth'])//ratio
            by_group[index].append(row)
            targets = (context['targetSpecies'],) if context.get('targetSpecies') else SPECIES
            for species in targets:
                component = row['components'][species]
                native_parameters.append([label,species,row['lo'],row['hi'],row['r'],
                    component['mean'],component['sigma'],component.get('meanError'),
                    component.get('sigmaError'),component.get('initialMean',component['mean']),
                    component.get('initialSigma',component['sigma']),row.get('targetCenterShiftSigma'),
                    row.get('targetWidthRatio'),int(component['selection']['targetFreeResponseSupported'])
                    if 'targetFreeResponseSupported' in component['selection'] else None,
                    int(row['numericalValid'])])
        for group, rows in sorted(by_group.items()):
            if len(rows) != ratio:
                raise ValueError('Incomplete native interval group; no hidden holes.')
            lo = cfg['rigidityMin']+group*cfg['aggregateWidth']
            hi = lo+cfg['aggregateWidth']
            targets = (context['targetSpecies'],) if context.get('targetSpecies') else SPECIES
            for species in targets:
                value = aggregate_one(rows, species, audit, cfg,raw_counts)
                detailed.append(dict(label=label, species=species, rLo=lo, rHi=hi, **value))
                groups[species, group][label] = value
    for (species, group), models in sorted(groups.items(), key=lambda x:(SPECIES.index(x[0][0]),x[0][1])):
        lo = cfg['rigidityMin']+group*cfg['aggregateWidth']
        hi = lo+cfg['aggregateWidth']
        available = [m['purity'] for m in models.values() if finite(m['purity'])]
        lower_values = [m['purityLowerIncludingOmittedBound'] for m in models.values()
                        if finite(m['purityLowerIncludingOmittedBound'])]
        lower = min(lower_values) if lower_values else None
        upper = max(available) if available else None
        numerical_support = models_supported(models)
        quality = bool(numerical_support and lower >= cfg['cutoff']['purityThreshold'])
        status = 'accepted_model_envelope' if quality else ('below_purity_threshold' if numerical_support else 'unassessed')
        ref = models['gaussian_linear_nominal']
        residuals = [m['localResidualFrac'] for m in models.values() if finite(m['localResidualFrac'])]
        combined.append(dict(species=species, rLo=lo, rHi=hi, purityNominal=ref['purity'],
            purityLow=lower, purityHigh=upper, rawSelected=ref['rawSelected'],
            localResidualMax=max(residuals) if residuals else None, quality=quality, status=status,
            modelsSupported=numerical_support, modelDiagnostics=models))
    recommendations = {}
    for species in SPECIES:
        rows = [r for r in combined if r['species']==species]
        recommendations[species] = recommendation(rows, cfg)
        cutoff = recommendations[species]['lowerCutoff']
        for row in rows:
            table.append([species, row['rLo'], row['rHi'], row['purityNominal'], row['purityLow'],
                          row['purityHigh'], row['rawSelected'], row['localResidualMax'], int(row['quality']), cutoff])
    def clean(rows):
        return [[float('nan') if x is None else x for x in row] for row in rows]
    tsv(output/'summary.tsv', 'species rLo rHi purityNominal purityLow purityHigh rawSelected localResidualMax quality cutoff', clean(table))
    runtime = module('summary_baseline_runtime', baseline_path.parent/'source/tools/nuclear_pid_tof/runtime_selection.py')
    curve = []
    for species_index, species in enumerate(SPECIES):
        for i in range(round(cfg['rigidityMin']/cfg['nativeWidth']),round(cfg['rigidityMax']/cfg['nativeWidth'])):
            r = (i+.5)*cfg['nativeWidth']
            mean, sigma, _ = runtime.current_parameters(baseline['current'], species_index, r)
            curve.append([species, r, mean, sigma, cfg['allFit']['maxNSigma']])
    tsv(output/'baseline_curves.tsv', 'species r mean sigma nSigma', curve)
    tsv(output/'native_parameters.tsv',
        'label species rLo rHi r mean sigma meanError sigmaError initialMean initialSigma centerShiftSigma widthRatio targetFreeSupported numericalValid',
        clean(native_parameters))
    save(output/'summary.json', dict(schemaVersion=1, productionEnabled=False,
        fullUncertaintyValidated=False, inputStudies=products, recommendations=recommendations,
        includedVariations=sorted({context['variation'] for label,context,result,audit in runs}),
        declaredShapeSensitivityComplete={context['variation'] for label,context,result,audit in runs}==set(cfg['shapeVariations']),
        aggregateRows=combined, detailedRows=detailed,
        gateAndFixedWindowSha256={target:dict(gate=v[0],fixedWindows=v[1]) for target,v in native_signatures.items()},
        policy='P=sum(S)/sum(N); lower/upper are response/BG/declared transfer-sensitivity extrema, '
               'not statistical confidence bounds. Local residual and identity support screened separately. '
               'Missing/failed support is unassessed, not zero purity.'))
    tsv(output/'model_results.tsv', 'label species rLo rHi purity signal total raw localResidual integralResidual unsupportedNuisance supported',
        clean([[r[k] for k in ('label','species','rLo','rHi','purity','signalSelected','totalSelected',
              'rawSelected','localResidualFrac','integralResidualFrac','unsupportedNuisanceFraction','supported')] for r in detailed]))
    shutil.copy2(args.settings, output/'settings.yaml')
    save(output/'context.json',dict(schemaVersion=1,productionEnabled=False,
        inputRunManifestSignatures=provenance,sourceSha256=source_signatures,
        rawGateProvenance=raw_provenance,
        baseline=str(baseline_path),summaryOnly=True,
        policy='Critical fit/root/tables/settings/source SHA verified; immutable original run snapshots retained.'))
    manifest = {str(p.relative_to(output)):digest(p) for p in output.rglob('*') if p.is_file()}
    save(output/'completed.json',dict(executionComplete=True,productionEnabled=False,
        physicsValidationComplete=False,summaryOnly=True,products=manifest))
    print('NUCLEAR_CUTOFF_SUMMARY_COMPLETE', output, flush=True)


if __name__ == '__main__':
    main()
