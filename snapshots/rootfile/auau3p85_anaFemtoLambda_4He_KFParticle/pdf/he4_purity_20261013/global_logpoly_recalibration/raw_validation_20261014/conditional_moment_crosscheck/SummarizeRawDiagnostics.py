#!/usr/bin/env python3
"""Read-only archive validation and conditional-moment regression, not a PID fit."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import sys

import yaml


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def verify(folder):
    completed = json.loads((folder / 'completed.json').read_text())
    context = json.loads((folder / 'context.json').read_text())
    if completed['phase'] != 'complete' or not completed['inputPreserved']:
        raise ValueError('Unfinished or modified-input study.')
    for name, expected in completed['products'].items():
        if digest(folder / name) != expected:
            raise ValueError('Product SHA differs: ' + str(folder / name))
    for name, expected in completed['sourceSha256'].items():
        if digest(folder / 'source' / name) != expected:
            raise ValueError('Source SHA differs: ' + str(folder / 'source' / name))
    p = Path(context['input']['path'])
    s = p.stat()
    actual = dict(path=str(p), size=s.st_size, mtime_ns=s.st_mtime_ns, device=s.st_dev, inode=s.st_ino)
    if actual != context['input']:
        raise ValueError('Original ROOT fingerprint differs.')
    return dict(inputFingerprint=actual, productsChecked=len(completed['products']),
                sourceSnapshotsChecked=len(completed['sourceSha256']),
                completionSha256=digest(folder / 'completed.json'), allHashesMatch=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('baseline', 'readable', 'settings', 'mainconf', 'output-dir', 'project'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError('Choose a new diagnostics directory; existing outputs are preserved.')
    baseline, readable = args.baseline.resolve(strict=True), args.readable.resolve(strict=True)
    cfg = yaml.safe_load(args.settings.read_text())
    if cfg['schemaVersion'] != 1 or cfg['productionEnabled'] is not False or cfg['candidateOnly'] is not True:
        raise ValueError('Diagnostic settings must explicitly exclude production use.')
    if cfg['weighting'] != 'equal_R_bins' or cfg['meanDegree'] != 3 or cfg['windowRmsDegree'] != 2:
        raise ValueError('Use the explicitly requested analytic dimensions and weighting.')
    context = json.loads((readable / 'context.json').read_text())
    mainconf = args.mainconf.resolve(strict=True)
    if str(mainconf) != context['runtimeMainconf'] or digest(mainconf) != context['runtimeMainconfSha256']:
        raise ValueError('Explicit runtime mainconf differs from the completed drawing study.')
    audit = dict(baseline=verify(baseline), readable=verify(readable),
                 runtimeMainconf=str(mainconf), runtimeMainconfSha256=digest(mainconf))
    unchanged = ('curves.tsv', 'knots.tsv', 'projection_audit.tsv', 'subset_audit.json',
                 'tof3_moment_summary.json', 'current_calibration.yaml', 'tof_nuclear_source.yaml')
    audit['unchangedNumericProducts'] = {name: digest(baseline / name) == digest(readable / name) for name in unchanged}
    if not all(audit['unchangedNumericProducts'].values()):
        raise ValueError('Drawing-only update changed numeric products.')
    audit['pdfPageObjects'] = {name: len(re.findall(rb'/Type\s*/Page\b', (readable / name).read_bytes()))
        for name in ('raw_dedx_calibration_overlay.pdf', 'low_rigidity_raw_projections.pdf', 'tof3_vs_all_shape_overlay.pdf')}
    helper = args.project.resolve(strict=True) / 'tools/he4_purity/calibration_candidates.py'
    sys.path.insert(0, str(helper.parent))
    from calibration_candidates import polynomial_fit, polynomial_value
    rows = [{key: float(value) for key, value in row.items()}
        for row in csv.DictReader((baseline / 'projection_audit.tsv').open(), delimiter='\t')]
    results = []
    for low in cfg['regressionRigidityMinimums']:
        chosen = [row for row in rows if low < row['rigidity'] < cfg['regressionRigidityMaximum']]
        x = [math.log10(row['rigidity']) for row in chosen]
        y = [row['tag_window_mean'] for row in chosen]
        width = [row['tag_window_rms'] for row in chosen]
        mean = polynomial_fit(x, y, cfg['meanDegree'])
        rms = polynomial_fit(x, width, cfg['windowRmsDegree'])
        results.append(dict(rigidityMinimum=low, rigidityMaximum=cfg['regressionRigidityMaximum'],
            nConditionalWindows=len(chosen), meanCoefficients=mean, windowRmsCoefficients=rms,
            meanResidualRms=math.sqrt(math.fsum((polynomial_value(mean, u) - v)**2 for u, v in zip(x, y)) / len(x)),
            windowRmsResidualRms=math.sqrt(math.fsum((polynomial_value(rms, u) - v)**2 for u, v in zip(x, width)) / len(x)),
            evaluations=[dict(rigidity=r, conditionalMeanPrediction=polynomial_value(mean, math.log10(r)),
                              conditionalWindowRmsPrediction=polynomial_value(rms, math.log10(r)),
                              extrapolation=not low <= r <= cfg['regressionRigidityMaximum'])
                         for r in cfg['evaluationRigidities']]))
    result = dict(schemaVersion=1, candidateOnly=True, productionEnabled=False,
        scope=cfg['scope'], audit=audit, regressions=results,
        coordinate='R=|p/q| GeV/c; mu/window RMS in keV/cm',
        limitations=['Conditional window mean/RMS are NOT Gaussian isotope mu/sigma.',
            'TOF3 is a subset of ALL, with TOF acceptance/quality and impurity biases.',
            'No fit likelihood, parameter covariance or isotope-identification claim is made.',
            'Different fitted R ranges change low-R extrapolation markedly.',
            'No production calibration or PID settings are changed.'])
    output.mkdir(parents=True, exist_ok=False)
    shutil.copy2(Path(__file__).resolve(), output / 'SummarizeRawDiagnostics.py')
    shutil.copy2(helper, output / 'calibration_candidates.py')
    shutil.copy2(args.settings, output / 'settings.yaml')
    (output / 'conditional_moment_crosscheck.json').write_text(json.dumps(result, indent=2, allow_nan=False) + '\n')
    completed = dict(phase='complete', productionEnabled=False,
        products={p.name: digest(p) for p in output.iterdir() if p.is_file()})
    (output / 'completed.json').write_text(json.dumps(completed, indent=2, allow_nan=False) + '\n')
    print('CONDITIONAL_MOMENT_CROSSCHECK_COMPLETE', output)


if __name__ == '__main__':
    main()
