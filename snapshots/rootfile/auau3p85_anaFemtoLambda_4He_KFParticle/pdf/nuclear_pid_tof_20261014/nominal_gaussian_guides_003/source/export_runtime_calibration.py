#!/usr/bin/env python3
"""Export the frozen 3p85 component fits and evaluate the exact new TPC gates.

No refitting, ROOT input modification, TOF purity inference, or CF correction.
The exported calibration intentionally contains model-dependent reference knots.
"""
import argparse
import bisect
import csv
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
TOL = 1.e-12
MEAN = ((10.2875, -36.6639, 67.63, -45.2483),
        (19.2177, -70.3198, 113.876, -67.0591),
        (28.8747, -68.8637, 104.455, -56.9342),
        (39.974, -103.365, 139.182, -67.1446))
SIGMA = ((.771189, -2.63609, 2.8044),
         (1.40996, -4.97147, 4.97608),
         (1.636, -4.35588, 4.40594),
         (2.24127, -7.08328, 8.07124))


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b''):
            h.update(chunk)
    return h.hexdigest()


def save(path, data):
    Path(path).write_text(json.dumps(data, indent=2, allow_nan=False) + '\n')


def legacy_parameters(species, rigidity):
    domain = (.2, 6.) if species < 2 else (.4, 5.)
    if not math.isfinite(rigidity) or not domain[0] <= rigidity <= domain[1]:
        return 0., 0.
    log = math.log10(rigidity)
    return tuple(math.fsum(c * log**i for i, c in enumerate(a))
                 for a in (MEAN[species], SIGMA[species]))


def calibrated_parameters(calibration, species, rigidity):
    """Matches the opt-in C++ provider, including supported half-bins in gaps."""
    if species < 2:
        return (*legacy_parameters(species, rigidity), 'legacy_unchanged')
    entry = calibration['species'][('he3', 'he4')[species - 2]]
    lo, hi = entry['domain']
    knots = entry['knots']
    if rigidity < lo - TOL or rigidity > hi + TOL:
        return (*legacy_parameters(species, rigidity), 'legacy_outside_domain')
    positions = [k['rigidity'] for k in knots]
    i = bisect.bisect_left(positions, rigidity)
    if i < len(knots) and abs(positions[i] - rigidity) <= TOL:
        return knots[i]['mean'], knots[i]['sigma'], 'calibrated_knot'
    if i > 0 and abs(positions[i - 1] - rigidity) <= TOL:
        return knots[i - 1]['mean'], knots[i - 1]['sigma'], 'calibrated_knot'
    if i == 0:
        return knots[0]['mean'], knots[0]['sigma'], 'calibrated_edge_clamp'
    if i == len(knots):
        return knots[-1]['mean'], knots[-1]['sigma'], 'calibrated_edge_clamp'
    left, right = knots[i - 1], knots[i]
    gap = right['rigidity'] - left['rigidity']
    if gap > calibration['maxInterpolationGap'] + TOL:
        half = entry['supportBinHalfWidth']
        if rigidity <= left['rigidity'] + half + TOL:
            return left['mean'], left['sigma'], 'calibrated_gap_edge_clamp'
        if rigidity >= right['rigidity'] - half - TOL:
            return right['mean'], right['sigma'], 'calibrated_gap_edge_clamp'
        return (*legacy_parameters(species, rigidity), 'legacy_unsupported_gap')
    fraction = (rigidity - left['rigidity']) / gap
    return tuple(left[key] + fraction * (right[key] - left[key])
                 for key in ('mean', 'sigma')) + ('calibrated_interpolation',)


def gaussian_integral(parameters, index, lo, hi):
    amplitude, mean, sigma = parameters[index:index + 3]
    scale = math.sqrt(2.) * sigma
    return amplitude * sigma * math.sqrt(math.pi / 2.) * (
        math.erf((hi - mean) / scale) - math.erf((lo - mean) / scale))


def model_integrals(fit, fit_lo, fit_hi, intervals):
    signal = background = 0.
    p = fit['parameters']
    for lo, hi in intervals:
        signal += gaussian_integral(p, 0, lo, hi)
        if fit['model'] == 'linear':
            background += p[3] * (hi - lo) + (p[4] - p[3]) / (2. * (fit_hi - fit_lo)) * (
                (hi - fit_lo)**2 - (lo - fit_lo)**2)
        elif abs(p[4]) < 1.e-8:
            background += p[3] * (hi - lo)
        else:
            background += p[3] * (math.exp(p[4] * hi) - math.exp(p[4] * lo)) / p[4]
        for index in range(5, len(p), 3):
            background += gaussian_integral(p, index, lo, hi)
    return signal, background


def best_species_intervals(rigidity, calibration, gate):
    """Partition at all |pull4|=|pull_j| roots; ties go to the earlier species."""
    parameters = [calibrated_parameters(calibration, s, rigidity)[:2] for s in range(4)]
    mean, sigma = parameters[3]
    lo, hi = mean - gate * sigma, mean + gate * sigma
    boundaries = [lo, hi]
    for other_mean, other_sigma in parameters[:3]:
        denominators = (other_sigma - sigma, other_sigma + sigma)
        numerators = (other_sigma * mean - sigma * other_mean,
                      other_sigma * mean + sigma * other_mean)
        for numerator, denominator in zip(numerators, denominators):
            if abs(denominator) <= TOL:
                continue
            root = numerator / denominator
            if lo < root < hi:
                boundaries.append(root)
    boundaries.sort()
    intervals = []
    for left, right in zip(boundaries, boundaries[1:]):
        if right - left <= TOL:
            continue
        center = (left + right) / 2.
        pull = abs((center - mean) / sigma)
        if all(pull < abs((center - m) / s) for m, s in parameters[:3]):
            if intervals and abs(intervals[-1][1] - left) <= TOL:
                intervals[-1] = (intervals[-1][0], right)
            else:
                intervals.append((left, right))
    return intervals


def make_calibration(data, context, fit_dir, settings):
    selected = []
    for index, slice_ in enumerate(data['slices']):
        if slice_['chosen'] < 0:
            continue
        fit = slice_['fits'][slice_['chosen']]
        if fit['valid'] and fit['signalEvidence']:
            selected.append((index, slice_, fit))
    if len(selected) != 205:
        raise ValueError('This profile requires the frozen 205 chosen component-evidence knots.')
    source = dict(collision='auau3p85', collisionEnergyGeV=3.85, year=2021, production='P24iy',
        sourceRoot=data['inputRoot'], sourceRootFingerprint=context['input'],
        fitSummary=str(fit_dir / 'purity.json'), fitSummarySha256=digest(fit_dir / 'purity.json'),
        fitSource=str(fit_dir / 'source/tools/he4_purity/FitHe4Purity.C'),
        fitSourceSha256=digest(fit_dir / 'source/tools/he4_purity/FitHe4Purity.C'),
        settings=str(fit_dir / 'settings.yaml'), settingsSha256=digest(fit_dir / 'settings.yaml'),
        applicationStatus='Explicitly enabled model-dependent working calibration; not high-precision isotope validation.',
        selectionPolicy='Preserve original chosen model; require valid fit and signalEvidence; no alternative model substitution.',
        fitQualityKnots=sum(f['quality'] for _, _, f in selected),
        referenceKnots=sum(not f['quality'] for _, _, f in selected),
        scope='Charge-inclusive, centrality-integrated pre-PID tracks; no joint TOF or pair-conditioned purity.',
        legacyHeader=str(PROJECT / 'StMaker/common/NuclearIdDeDxVsMom.h'),
        legacyHeaderSha256=digest(PROJECT / 'StMaker/common/NuclearIdDeDxVsMom.h'))
    calibration = dict(schemaVersion=1, profile='auau3p85_2021P24iy_he3_he4_from040_20261013',
        rigidityUnit='GeV/c', dedxUnit='keV/cm', interpolation='linear',
        edgePolicy='clampWithinDomain', outsideRangePolicy='legacy', maxInterpolationGap=.015,
        source=source, species={})
    for species, name in ((2, 'he3'), (3, 'he4')):
        knots = []
        for index, slice_, fit in selected:
            parameter = 0 if species == 3 else 5 + 3 * fit['neighborSpecies'].index(2)
            amplitude, norm_mean, norm_sigma = fit['parameters'][parameter:parameter + 3]
            errors = fit['parameterErrors'][parameter:parameter + 3]
            mean = slice_['calibrationMean'] + slice_['calibrationSigma'] * norm_mean
            sigma = slice_['calibrationSigma'] * norm_sigma
            if not (math.isfinite(mean) and math.isfinite(sigma) and mean > 0 and sigma > 0 and
                    errors[0] > 0 and amplitude / errors[0] >= settings['minimumSignalSignificance']):
                raise ValueError('Selected component lacks finite positive parameters or amplitude evidence.')
            if species == 2:
                old_mean, old_sigma = legacy_parameters(2, slice_['r'])
                nm = (old_mean - slice_['calibrationMean']) / slice_['calibrationSigma']
                ns = old_sigma / slice_['calibrationSigma']
                mean_low = nm - settings['neighborMeanShiftMax'] * ns
                signal_low = (settings['signalLowerMeanMin'] if slice_['r'] < .5 else
                              settings['signalLowMeanMin'] if slice_['r'] < .85 else
                              -settings['signalMeanShiftMax'])
                mean_high = min(nm + settings['neighborMeanShiftMax'] * ns,
                                signal_low - settings['minimumSignalNeighborSeparation'])
                sigma_low, sigma_high = settings['neighborSigmaFactorMin'] * ns, settings['neighborSigmaFactorMax'] * ns
                if any(abs(value - bound) < 1.e-4 * (high - low)
                       for value, low, high in ((norm_mean, mean_low, mean_high),
                                                (norm_sigma, sigma_low, sigma_high))
                       for bound in (low, high)):
                    raise ValueError('He3 component at a material mean/width boundary.')
            knots.append(dict(rigidity=slice_['r'], mean=mean, sigma=sigma, provenance=dict(
                sourceIndex=index, model=fit['model'], fitQuality=bool(fit['quality']),
                referenceOnly=not fit['quality'], sourceHe4SignalEvidence=True,
                overallNuisanceBoundary=bool(fit['nuisanceAtBoundary']),
                componentMeanWidthInterior=True, sourceHe4RawLocalMaximum=bool(fit['rawLocalMaximum']),
                componentAmplitudeOverError=amplitude / errors[0],
                meanStatError=errors[1] * slice_['calibrationSigma'],
                sigmaStatError=errors[2] * slice_['calibrationSigma'])))
        calibration['species'][name] = dict(domain=[.41, 2.5], supportBinHalfWidth=.005, knots=knots)
    return calibration


def make_purity(data, calibration, settings):
    rows = []
    for slice_ in data['slices']:
        r = slice_['r']
        mean, sigma, policy = calibrated_parameters(calibration, 3, r)
        selections = {'old_gate': [(slice_['calibrationMean'] - 2 * slice_['calibrationSigma'],
                                    slice_['calibrationMean'] + 2 * slice_['calibrationSigma'])],
                      'new_tpc_gate': [(mean - 2 * sigma, mean + 2 * sigma)],
                      'new_best_species_proxy': best_species_intervals(r, calibration, 2.)}
        fit_lo = max(settings['fitNSigmaMin'], -slice_['calibrationMean'] / slice_['calibrationSigma'])
        fit_hi = min(settings['fitNSigmaMax'], (120. - slice_['calibrationMean']) / slice_['calibrationSigma'])
        for selection, raw_intervals in selections.items():
            intervals = [((lo - slice_['calibrationMean']) / slice_['calibrationSigma'],
                          (hi - slice_['calibrationMean']) / slice_['calibrationSigma'])
                         for lo, hi in raw_intervals]
            support = bool(intervals) and all(lo >= fit_lo - TOL and hi <= fit_hi + TOL for lo, hi in intervals)
            for index, fit in enumerate(slice_['fits']):
                row = dict(rigidity=r, selection=selection, model=fit['model'], chosen=index == slice_['chosen'],
                    fit_valid=fit['valid'], signal_evidence=fit['signalEvidence'], fit_quality=fit['quality'], reference_only=not fit['quality'],
                    calibration_policy=policy, fit_support_complete=support, intervals_dedx=raw_intervals,
                    purity=None, signal_model_count=None, background_model_count=None, signal_relative_old_gate=None)
                if support and fit['valid'] and fit['signalEvidence']:
                    signal, background = model_integrals(fit, fit_lo, fit_hi, intervals)
                    width = data['fitDedxBinWidth'] / slice_['calibrationSigma']
                    if signal >= 0 and background >= 0 and signal + background > 0:
                        row.update(purity=signal / (signal + background), signal_model_count=signal / width,
                            background_model_count=background / width,
                            signal_relative_old_gate=(signal / width / fit['gateSignal']) if fit['gateSignal'] > 0 else None)
                rows.append(row)
    for row in rows:
        alternatives = [other['purity'] for other in rows if other['rigidity'] == row['rigidity'] and
                        other['selection'] == row['selection'] and other['purity'] is not None]
        row['model_min'] = min(alternatives) if alternatives else None
        row['model_max'] = max(alternatives) if alternatives else None
        row['model_spread'] = max(alternatives) - min(alternatives) if len(alternatives) > 1 else None
    return rows


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-dir', type=Path, required=True)
    parser.add_argument('--calibration-output', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--mainconf', type=Path, required=True)
    args = parser.parse_args()
    fit_dir = args.fit_dir.resolve(strict=True)
    output = args.output_dir.resolve()
    calibration_path = args.calibration_output.resolve()
    if calibration_path.exists() or output.exists():
        raise ValueError('Refusing to overwrite a calibration or QA output.')
    completed = json.loads((fit_dir / 'completed.json').read_text())
    context = json.loads((fit_dir / 'context.json').read_text())
    data = json.loads((fit_dir / 'purity.json').read_text())
    settings = yaml.safe_load((fit_dir / 'settings.yaml').read_text())
    if completed['phase'] != 'complete' or '2021P24iy' not in data['inputRoot'] or context['sample'] != '3p85':
        raise ValueError('Only the completed 2021 P24iy 3p85 fit is accepted.')
    for name in ('purity.json', 'settings.yaml', 'source/tools/he4_purity/FitHe4Purity.C'):
        if digest(fit_dir / name) != completed['products'][name]:
            raise ValueError('Frozen fit product hash differs: ' + name)
    source_stat = Path(data['inputRoot']).stat()
    expected = context['input']
    if (source_stat.st_size, source_stat.st_mtime_ns, source_stat.st_ino, source_stat.st_dev) != (
            expected['size'], expected['mtime_ns'], expected['inode'], expected['device']):
        raise ValueError('Original ROOT fingerprint changed.')
    calibration = make_calibration(data, context, fit_dir, settings)
    output.mkdir(parents=True, exist_ok=False)
    (output / 'source').mkdir()
    for filename in ('export_runtime_calibration.py', 'PlotCalibrationWindows.C'):
        shutil.copy2(str(HERE / filename), str(output / 'source' / filename))
    shutil.copy2(str(PROJECT / 'StMaker/common/NuclearIdDeDxVsMom.h'),
                 str(output / 'source/NuclearIdDeDxVsMom.h'))
    calibration_path.parent.mkdir(parents=True, exist_ok=True)
    calibration_path.write_text(yaml.safe_dump(calibration, sort_keys=False, allow_unicode=True))
    rows = make_purity(data, calibration, settings)
    save(output / 'purity_windows.json', dict(calibration=str(calibration_path),
        calibration_sha256=digest(calibration_path), scope=calibration['source']['scope'],
        limitations=['No new-window statistical errors are synthesized from missing full covariance.',
                     'The model spread is a comparison of fixed common selections, not a calibrated systematic error.',
                     'Best-species is a TPC-only d/t/he3/he4 proxy; no TOF or pair conditioning.',
                     'The frozen background functions are not refitted.'], rows=rows))
    fields = [key for key in rows[0] if key != 'intervals_dedx'] + ['intervals_dedx']
    with (output / 'purity_windows.csv').open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for row in rows:
            writer.writerow(dict(row, intervals_dedx=json.dumps(row['intervals_dedx'])))
    selection_ids = {'old_gate': 0, 'new_tpc_gate': 1, 'new_best_species_proxy': 2}
    with (output / 'plot_points.tsv').open('w') as handle:
        for row in rows:
            if row['purity'] is not None:
                handle.write('%.17g\t%d\t%d\t%d\t%d\t%.17g\t%.17g\t%.17g\t%.17g\n' % (
                    row['rigidity'], selection_ids[row['selection']], int(row['model'] == 'exponential'),
                    row['chosen'], row['fit_quality'], row['purity'], row['model_min'], row['model_max'],
                    row['signal_relative_old_gate']))
    macro = output / 'source/PlotCalibrationWindows.C'
    with (output / 'calibration_curves.tsv').open('w') as handle:
        for slice_ in data['slices']:
            for species in (2, 3):
                mean, sigma, policy = calibrated_parameters(calibration, species, slice_['r'])
                old_mean, old_sigma = legacy_parameters(species, slice_['r'])
                handle.write('%.17g\t%d\t%.17g\t%.17g\t%.17g\t%.17g\t%d\n' % (
                    slice_['r'], species, mean, sigma, old_mean, old_sigma, policy.startswith('calibrated_')))
    call = str(macro) + '++(' + ','.join(json.dumps(str(p)) for p in (output / 'plot_points.tsv', output / 'calibration_curves.tsv', output)) + ')'
    with (output / 'plot.log').open('w') as handle:
        result = subprocess.run(['bash', str(PROJECT / 'script/singularity_root_batch.sh'),
            str(args.mainconf.resolve(strict=True)), 'root4star', '-l', '-b', '-q', call], cwd=str(PROJECT),
            stdout=handle, stderr=subprocess.STDOUT, timeout=600)
    if result.returncode or 'CALIBRATION_WINDOWS_COMPLETE' not in (output / 'plot.log').read_text():
        raise RuntimeError('Purity plot failed; data/calibration are retained: ' + str(output / 'plot.log'))
    (output / 'README.md').write_text('# 新しいヘリウム校正のTPC選択窓QA\n\n'
        '[Purity比較PDF](purity_calibrated_windows.pdf)は元校正±2、新校正±2、新校正best-speciesの3選択を比較します。'
        '背景関数は保存済みlinear/exponentialを同じ選択窓で積分しています。\n\n'
        '205のchosen component-evidence knotsを両核種へ適用し、6品質通過／199参考fitを区別します。'
        '範囲端は支持half-bin内のみclamp、未支持gapや範囲外はlegacy校正です。0.405は採用していません。\n\n'
        'charge-inclusive・中央度積分のpre-PID分布のTPCモデル推定で、TOF後やペア別の確定純度ではありません。'
        '保存されていない共分散から新しい統計誤差を作らず、CFへの純度補正も自動適用しません。\n\n'
        '`purity_windows.json/csv`にmodel別の純度、信号数、元gate比、窓の区間・支持範囲を保存します。\n')
    save(output / 'completed.json', dict(phase='complete', calibration=str(calibration_path),
        calibration_sha256=digest(calibration_path), fit_json_sha256=digest(fit_dir / 'purity.json'),
        source_sha256={str(p): digest(p) for p in (Path(__file__).resolve(), macro)},
        products={str(p.relative_to(output)): digest(p) for p in output.rglob('*') if p.is_file()
                  and p.name != 'completed.json' and p.suffix not in ('.so', '.d', '.o')}))
    print('COMPLETE', calibration_path, output, flush=True)


if __name__ == '__main__':
    main()
