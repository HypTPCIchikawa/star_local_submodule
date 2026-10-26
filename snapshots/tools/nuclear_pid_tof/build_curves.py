#!/usr/bin/env python3
"""Compare smooth TAG-shape candidates. Outputs are explicitly not production PID."""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil

import yaml

from nuclear_curves import evaluate_segment, make_segment
from runtime_selection import SPECIES, digest, gl3_nodes

PRIMARY = {'d': 'TAG_D', 't': 'TAG_T', 'he3': 'TAG_HE3', 'he4': 'TAG_D'}


def numeric_failures(row, name, settings):
    failures = []
    c = row.get('components', {}).get(name)
    if not c:
        return ['missing_component']
    if not row['numericalConverged']:
        failures.append('not_numerically_converged')
    if row.get('atBoundary', True):
        failures.append('shape_boundary')
    if row.get('conditionalNdf', 0) <= 0 or row.get('chi2', float('inf'))/max(1, row.get('conditionalNdf', 0)) > settings['maximumPearsonChi2Ndf']:
        failures.append('poor_spectrum_shape')
    for key, gate in (('effectiveCount', 'minimumEffectiveSignalCount'),
                      ('windowYield', 'minimumWindowSignalCount'),
                      ('conditionalDropDeltaD', 'minimumConditionalDropDeltaD'),
                      ('detectorAcceptance', 'minimumPhysicalGaussianCoverage')):
        if not math.isfinite(c[key]) or c[key] < settings[gate]:
            failures.append(key+'_below_gate')
    if c['sigma'] <= 0 or c['mean'] <= 0 or c['centerOutsideDetector']:
        failures.append('invalid_or_truncated_center')
    if c['meanError']/max(c['sigma'], 1.e-99) > settings['maximumMeanErrorOverSigma']:
        failures.append('poor_mean_precision')
    if c['sigmaError']/max(c['sigma'], 1.e-99) > settings['maximumSigmaRelativeError']:
        failures.append('poor_width_precision')
    if c['meanError'] <= 0 or c['sigmaError'] <= 0:
        failures.append('nonpositive_shape_error')
    return failures


def keyed_rows(document):
    return {(row['sample'], round(row['lo'], 9), round(row['hi'], 9)): row
            for row in document['slices']}


def collect_points(linear, exponential, settings):
    rows = [keyed_rows(d) for d in (linear, exponential)]
    if set(rows[0]) != set(rows[1]):
        raise ValueError('Background comparison requires identical native R/sample slices.')
    accepted, audit = {name: [] for name in SPECIES}, []
    for name in SPECIES:
        for key in sorted(k for k in rows[0] if k[0] == PRIMARY[name]):
            a, b = rows[0][key], rows[1][key]
            if abs(a['hi']-a['lo']-.01) > 1.e-9:
                raise ValueError('Primary curve calibration requires native 0.01 R projections.')
            if (a['fitLo'], a['fitHi']) != (b['fitLo'], b['fitHi']):
                raise ValueError('Background comparison must use identical dE/dx windows.')
            failures = [m+':'+f for m, row in zip(('linear', 'exponential'), (a, b))
                        for f in numeric_failures(row, name, settings)]
            ca, cb = (row.get('components', {}).get(name) for row in (a, b))
            if ca and cb and ca['sigma'] > 0 and cb['sigma'] > 0:
                sigma = .5*(ca['sigma']+cb['sigma'])
                shift_mu = abs(ca['mean']-cb['mean'])/sigma
                shift_sigma = abs(ca['sigma']-cb['sigma'])/sigma
                if shift_mu > settings['maximumBackgroundMeanShiftOverSigma']:
                    failures.append('background_model_mean_instability')
                if shift_sigma > settings['maximumBackgroundSigmaRelativeShift']:
                    failures.append('background_model_width_instability')
            else:
                shift_mu = shift_sigma = None
            record = dict(species=name, sample=key[0], lo=a['lo'], hi=a['hi'], r=a['r'],
                          failures=failures, acceptedNumericalShape=not failures,
                          identityValidated=False, modelMeanShiftOverSigma=shift_mu,
                          modelWidthRelativeShift=shift_sigma)
            if ca and cb:
                record.update(mean=.5*(ca['mean']+cb['mean']),
                              projectedSigma=.5*(ca['sigma']+cb['sigma']),
                              meanError=max(ca['meanError'], cb['meanError'], abs(ca['mean']-cb['mean'])/2),
                              sigmaError=max(ca['sigmaError'], cb['sigmaError'], abs(ca['sigma']-cb['sigma'])/2))
            audit.append(record)
            if not failures:
                accepted[name].append(dict(record, intrinsicSigmaSeed=record['projectedSigma']))
    # A TAG_D He4 decomposition cannot validate fixed nuisance shapes by itself.
    if settings['requireIndependentTagNuisanceValidationForHe4']:
        supported = {name: {(p['lo'], p['hi']) for p in accepted[name]} for name in ('t', 'he3')}
        keep = []
        for point in accepted['he4']:
            missing = [name for name in ('t', 'he3') if (point['lo'], point['hi']) not in supported[name]]
            if missing:
                point['failures'].append('nuisance_TAG_shape_not_supported:'+','.join(missing))
                point['acceptedNumericalShape'] = False
                # audit record shares the mutable failures list, but not the bool.
                for record in audit:
                    if record['species'] == 'he4' and record['lo'] == point['lo']:
                        record['acceptedNumericalShape'] = False
            else:
                keep.append(point)
        accepted['he4'] = keep
    return accepted, audit


def groups(points, max_gap):
    result = []
    for point in points:
        if not result or point['r']-result[-1][-1]['r'] > max_gap+1.e-9:
            result.append([])
        result[-1].append(point)
    return result


def build_candidate(accepted, family, settings):
    output = dict(schemaVersion=1, productionEnabled=False, identityValidated=False,
                  family=family, species={}, uncertaintyScope='conditional TAG errors plus background-model half-difference; not total uncertainty',
                  transferValidated=False, widthCalibrationValidated=False,
                  varianceUnfoldPolicy='Uniform-R GL3 subtraction is a seed approximation; no negative-variance clipping',
                  purityMeaning='nominal Gaussian decomposition conditional on tagged shape transfer')
    for name in SPECIES:
        segments, rejected = [], []
        for group in groups(accepted[name], settings['support']['maximumGapRigidity']):
            if len(group) < settings['support']['minimumAcceptedPoints']:
                rejected.append(dict(domain=[group[0]['lo'], group[-1]['hi']], reason='too_few_supported_points', count=len(group)))
                continue
            points = [dict(p) for p in group]
            try:
                segment = make_segment(points, family, settings['curves'])
                for _ in range(settings['curves']['varianceUnfoldIterations']):
                    for point in points:
                        means = [(evaluate_segment(segment, r)['mean'], w) for r, w in gl3_nodes(point['lo'], point['hi'])]
                        center = sum(m*w for m, w in means)
                        variance = sum(w*(m-center)**2 for m, w in means)
                        remainder = point['projectedSigma']**2-variance
                        if remainder <= 0:
                            raise ValueError('Nonpositive intrinsic variance at R='+str(point['r']))
                        point['intrinsicSigmaSeed'] = math.sqrt(remainder)
                        point['rigiditySmearingVariance'] = variance
                    segment = make_segment(points, family, settings['curves'])
                segment['supportedNativeIntervals'] = [[p['lo'], p['hi']] for p in points]
                segment['pointCount'] = len(points)
                segment['calibrationPoints'] = points
                segments.append(segment)
            except ValueError as error:
                rejected.append(dict(domain=[group[0]['lo'], group[-1]['hi']], reason=str(error), count=len(group)))
        output['species'][name] = dict(segments=segments, rejectedGroups=rejected,
                                      lowestNumericallySupportedRigidity=min((s['domain'][0] for s in segments), default=None),
                                      reliableIsotopeLowerLimit=None)
    return output


def main():
    parser = argparse.ArgumentParser(__doc__)
    parser.add_argument('--linear-dir', type=Path, required=True)
    parser.add_argument('--exponential-dir', type=Path, required=True)
    parser.add_argument('--settings', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    settings = yaml.safe_load(args.settings.read_text())
    inputs = []
    for folder in (args.linear_dir, args.exponential_dir):
        folder = folder.resolve(strict=True)
        complete = json.loads((folder/'completed.json').read_text())
        if not complete.get('executionComplete') or complete.get('productionEnabled') is not False:
            raise ValueError('TAG execution has not completed as a study-only output.')
        for relative, expected in complete['products'].items():
            product = (folder/relative).resolve(strict=True)
            if folder not in product.parents or digest(product) != expected:
                raise ValueError('TAG product path/hash mismatch: '+relative)
        if 'tag_fits.json' not in complete['products']:
            raise ValueError('Authoritative TAG fit is missing from completion manifest.')
        fit = json.loads((folder/'tag_fits.json').read_text())
        if (fit.get('fitStage') != 'independent_tag_shapes'
                or fit.get('responseModel', 'gaussian') != 'gaussian'):
            raise ValueError('This curve builder is for Gaussian projected mean/std only, not another response parameterization.')
        inputs.append(fit)
    accepted, audit = collect_points(*inputs, settings['support'])
    args.output_dir.mkdir(parents=True, exist_ok=False)
    context = dict(settingsPath=str(args.settings.resolve()), settingsSha256=digest(args.settings),
                   tagSources={str((p/'tag_fits.json').resolve()): digest(p/'tag_fits.json')
                               for p in (args.linear_dir, args.exponential_dir)},
                   productionEnabled=False)
    (args.output_dir/'context.json').write_text(json.dumps(context, indent=2)+'\n')
    snapshots = args.output_dir/'source'
    snapshots.mkdir()
    for source in (Path(__file__), Path(__file__).parent/'nuclear_curves.py',
                   Path(__file__).parent/'runtime_selection.py', args.settings):
        target = snapshots/source.name
        shutil.copy2(source, target)
        context.setdefault('snapshots', {})[target.name] = digest(target)
    (args.output_dir/'context.json').write_text(json.dumps(context, indent=2)+'\n')
    (args.output_dir/'support_audit.json').write_text(json.dumps(audit, indent=2, allow_nan=False)+'\n')
    with (args.output_dir/'support_audit.tsv').open('w') as stream:
        writer = csv.writer(stream, delimiter='\t')
        writer.writerow(['species', 'Rlo', 'Rhi', 'acceptedNumericalShape', 'failures', 'identityValidated'])
        for p in audit:
            writer.writerow([p['species'], p['lo'], p['hi'], int(p['acceptedNumericalShape']), ';'.join(p['failures']), 0])
    products = []
    for family in ('logpoly', 'pchip_logR'):
        candidate = build_candidate(accepted, family, settings)
        path = args.output_dir/(family+'.json')
        path.write_text(json.dumps(candidate, indent=2, allow_nan=False)+'\n')
        products.append(path)
    (args.output_dir/'completed.json').write_text(json.dumps(dict(
        executionComplete=True, productionEnabled=False,
        products={p.name: digest(p) for p in products},
        acceptedNumericalPoints={name: len(p) for name, p in accepted.items()},
        scientificValidationComplete=False), indent=2)+'\n')
    print('NUCLEAR_STUDY_CURVES_COMPLETE', args.output_dir)


if __name__ == '__main__':
    main()
