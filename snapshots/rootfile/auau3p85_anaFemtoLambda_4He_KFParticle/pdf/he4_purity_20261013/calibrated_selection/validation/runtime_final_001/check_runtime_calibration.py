#!/usr/bin/env python3
"""Independent mask/integration closure and Python/C++ calibration comparison."""
import argparse
import csv
import json
import math
from pathlib import Path
import shutil

import yaml

from export_runtime_calibration import (best_species_intervals, calibrated_parameters,
    digest, model_integrals, save)


def direct_best(dedx, rigidity, calibration):
    winner = -1
    closest = 2.
    for species in range(4):
        mean, sigma, _ = calibrated_parameters(calibration, species, rigidity)
        if sigma <= 0:
            continue
        pull = abs((dedx - mean) / sigma)
        if pull < closest:
            winner, closest = species, pull
    return winner


def mask_check(calibration, rigidity, points=1001):
    intervals = best_species_intervals(rigidity, calibration, 2.)
    mean, sigma, _ = calibrated_parameters(calibration, 3, rigidity)
    lo, hi = mean - 2 * sigma, mean + 2 * sigma
    checked = skipped = 0
    for i in range(points):
        dedx = lo + (hi - lo) * (i + .5) / points
        pulls = sorted(abs((dedx - m) / s) for m, s, _ in
                       (calibrated_parameters(calibration, j, rigidity) for j in range(4)) if s > 0)
        if len(pulls) > 1 and abs(pulls[0] - pulls[1]) < 1.e-10:
            skipped += 1
            continue
        expected = direct_best(dedx, rigidity, calibration) == 3
        actual = any(left < dedx < right for left, right in intervals)
        if actual != expected:
            raise AssertionError(('best-species mask mismatch', rigidity, dedx, intervals))
        checked += 1
    return dict(rigidity=rigidity, checked=checked, skipped_exact_ties=skipped, intervals=intervals)


def density(fit, coordinate, fit_lo, fit_hi):
    p = fit['parameters']
    signal = p[0] * math.exp(-.5 * ((coordinate - p[1]) / p[2])**2)
    if fit['model'] == 'linear':
        background = p[3] + (p[4] - p[3]) * (coordinate - fit_lo) / (fit_hi - fit_lo)
    else:
        background = p[3] * math.exp(p[4] * coordinate)
    for index in range(5, len(p), 3):
        background += p[index] * math.exp(-.5 * ((coordinate - p[index + 1]) / p[index + 2])**2)
    return signal, background


def simpson_raw(fit, slice_, fit_lo, fit_hi, intervals, steps=4000):
    """Independent raw-dEdx integration, including the coordinate Jacobian."""
    total = [0., 0.]
    for lo, hi in intervals:
        step = (hi - lo) / steps
        contributions = [[], []]
        for index in range(steps + 1):
            coordinate = (lo + index * step - slice_['calibrationMean']) / slice_['calibrationSigma']
            pair = density(fit, coordinate, fit_lo, fit_hi)
            weight = 1 if index in (0, steps) else 4 if index % 2 else 2
            for component in range(2):
                contributions[component].append(weight * pair[component])
        for component in range(2):
            total[component] += math.fsum(contributions[component]) * step / (3. * slice_['calibrationSigma'])
    return total


def numerical_closure(calibration, data, settings):
    cases = []
    for rigidity in (.415, .425, .455, .505, .605, 1.605, 2.005, 2.495):
        slice_ = min(data['slices'], key=lambda item: abs(item['r'] - rigidity))
        mean, sigma, _ = calibrated_parameters(calibration, 3, slice_['r'])
        selections = {'new_tpc': [(mean - 2 * sigma, mean + 2 * sigma)],
                      'new_best': best_species_intervals(slice_['r'], calibration, 2.)}
        fit_lo = max(settings['fitNSigmaMin'], -slice_['calibrationMean'] / slice_['calibrationSigma'])
        fit_hi = min(settings['fitNSigmaMax'], (120. - slice_['calibrationMean']) / slice_['calibrationSigma'])
        for fit in slice_['fits']:
            if not fit['valid'] or not fit['signalEvidence']:
                continue
            for name, intervals in selections.items():
                normalized = [((lo - slice_['calibrationMean']) / slice_['calibrationSigma'],
                               (hi - slice_['calibrationMean']) / slice_['calibrationSigma'])
                              for lo, hi in intervals]
                analytical = model_integrals(fit, fit_lo, fit_hi, normalized)
                numerical = simpson_raw(fit, slice_, fit_lo, fit_hi, intervals)
                errors = [abs(a - n) / max(abs(a), 1.) for a, n in zip(analytical, numerical)]
                if max(errors) > 1.e-8:
                    raise AssertionError(('integral closure mismatch', rigidity, fit['model'], name, errors))
                cases.append(dict(rigidity=slice_['r'], model=fit['model'], selection=name,
                    analytical=list(analytical), independent_raw_simpson=numerical, relative_errors=errors))
    return cases


def toy_calibration(he3, he4):
    return dict(maxInterpolationGap=.015, species={name: dict(domain=[.995, 1.005], supportBinHalfWidth=.005,
        knots=[dict(rigidity=1., mean=mean, sigma=sigma)]) for name, (mean, sigma) in zip(('he3', 'he4'), (he3, he4))})


def runtime_compare(path, calibration):
    """Consume the standalone C++ provider TSV, not a Python-generated fixture."""
    with Path(path).open() as handle:
        reader = csv.DictReader(handle, delimiter='\t')
        fields = reader.fieldnames
        rows = list(reader)
    if not rows:
        raise ValueError('Empty runtime audit TSV.')
    required = ('rigidity', 'species', 'mean', 'sigma')
    if not all(name in fields for name in required):
        raise ValueError('Runtime TSV requires fields: ' + ', '.join(required) + '; actual=' + str(fields))
    maximum = [0., 0.]
    for row in rows:
        species = int(row['species'])
        rigidity = float(row['rigidity'])
        mean, sigma, policy = calibrated_parameters(calibration, species, rigidity)
        differences = [abs(float(row['mean']) - mean), abs(float(row['sigma']) - sigma)]
        maximum = [max(a, b) for a, b in zip(maximum, differences)]
        if differences[0] > 1.e-10 * max(abs(mean), 1.) or differences[1] > 1.e-10 * max(abs(sigma), 1.):
            raise AssertionError(('runtime/export curve mismatch', row, mean, sigma, policy))
        flag = 'calibrated' if 'calibrated' in fields else 'custom' if 'custom' in fields else None
        if flag and bool(int(row[flag])) != policy.startswith('calibrated_'):
            raise AssertionError(('runtime/export policy mismatch', row, policy))
    return dict(rows=len(rows), maximum_absolute_difference_mean=maximum[0],
        maximum_absolute_difference_sigma=maximum[1], source=str(Path(path).resolve()), source_sha256=digest(path))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--calibration', type=Path, required=True)
    parser.add_argument('--qa-dir', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--runtime-tsv', type=Path)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError('Use a fresh validation folder to retain historical results.')
    qa = args.qa_dir.resolve(strict=True)
    calibration = yaml.safe_load(args.calibration.read_text())
    completed = json.loads((qa / 'completed.json').read_text())
    if digest(args.calibration) != completed['calibration_sha256']:
        raise ValueError('Calibration differs from the frozen published QA.')
    if any(digest(qa / name) != checksum for name, checksum in completed['products'].items()):
        raise ValueError('A published QA product changed.')
    source = Path(calibration['source']['fitSummary'])
    data = json.loads(source.read_text())
    settings = yaml.safe_load(Path(calibration['source']['settings']).read_text())
    mask_cases = [mask_check(calibration, s['r']) for s in data['slices']]
    twins = toy_calibration((40., 3.), (40., 3.))
    if best_species_intervals(1., twins, 2.):
        raise AssertionError('Identical He3 and He4 pulls must never choose He4 in a tie.')
    if any(direct_best(dedx, 1., twins) != 2 for dedx in (39., 40., 41.)):
        raise AssertionError('The independent direct loop must resolve equal He3/He4 pulls to He3.')
    islands = toy_calibration((40., 1.), (42., 5.))
    if len(best_species_intervals(1., islands, 2.)) != 2:
        raise AssertionError('The unequal-width toy must retain both separated accepted intervals.')
    toy_masks = [mask_check(twins, 1.), mask_check(islands, 1.), mask_check(toy_calibration((38., 3.), (42., 3.)), 1.)]
    integrals = numerical_closure(calibration, data, settings)
    window_data = json.loads((qa / 'purity_windows.json').read_text())
    old_error = 0.
    for row in window_data['rows']:
        if row['selection'] == 'old_gate' and row['purity'] is not None:
            slice_ = min(data['slices'], key=lambda s: abs(s['r'] - row['rigidity']))
            fit = next(f for f in slice_['fits'] if f['model'] == row['model'])
            old_error = max(old_error, abs(row['purity'] - fit['purity']))
    if old_error > 1.e-12:
        raise AssertionError('The original-gate purity is not reproduced.')
    runtime = runtime_compare(args.runtime_tsv, calibration) if args.runtime_tsv else None
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(str(Path(__file__).resolve()), str(output / 'check_runtime_calibration.py'))
    if args.runtime_tsv:
        shutil.copy2(str(args.runtime_tsv), str(output / 'runtime_curve_audit.tsv'))
        runtime['preserved_copy'] = str(output / 'runtime_curve_audit.tsv')
    save(output / 'validation.json', dict(status='PASS', calibration_sha256=digest(args.calibration),
        source_fit_json_sha256=digest(source), published_qa_hashes_unchanged=True,
        mask_cases=mask_cases, mask_sample_count=sum(s['checked'] for s in mask_cases),
        toy_masks=toy_masks, numerical_integral_cases=integrals,
        maximum_original_gate_purity_difference=old_error, runtime_curve_comparison=runtime,
        source_sha256={str(p): digest(p) for p in (Path(__file__).resolve(), Path(__file__).parent / 'export_runtime_calibration.py')}))
    (output / 'README.md').write_text('# 新校正windowの独立検証\n\n'
        '全210点のbest-speciesを、解析的な区間分割とは別の4核種closest-pullループで照合しました。'
        '同一pullのtie、幅差による分離した2区間、同じ幅の例も検証します。\n\n'
        '新TPC窓と新best窓の信号・背景積分をraw dE/dx座標の独立Simpson積分と比較し、'
        '元gateの純度も再現します。全検証は校正を変更しません。\n\n'
        'C++ runtime TSVが指定された場合は、実provider出力のcenter・width・校正適用flagをPythonと照合します。'
        '詳細・source SHAはvalidation.jsonに保存します。\n')
    print('PASS', output, 'mask samples', sum(s['checked'] for s in mask_cases), 'integral cases', len(integrals), 'runtime', runtime, flush=True)


if __name__ == '__main__':
    main()
