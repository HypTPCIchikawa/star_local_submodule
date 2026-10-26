#!/usr/bin/env python3
"""Draw unselected dE/dx and independent TOF-tag QA; never modify PID or refit."""
import argparse
import csv
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def fingerprint(path):
    path = Path(path).resolve(strict=True)
    value = path.stat()
    return dict(path=str(path), size=value.st_size, mtime_ns=value.st_mtime_ns,
                device=value.st_dev, inode=value.st_ino)


def save(path, value):
    path.write_text(json.dumps(value, indent=2, allow_nan=False) + '\n')


def candidate_parameters(candidate, species, r, legacy):
    curve = candidate['species'][('he3', 'he4')[species - 2]]
    lo, hi = curve['domain']
    if not lo <= r <= hi:
        return (*legacy(species, r), False)
    log = math.log10(r)
    mean = math.fsum(c * log**i for i, c in enumerate(curve['meanCoefficients']))
    sigma = math.fsum(c * log**i for i, c in enumerate(curve['sigmaCoefficients']))
    if not all(math.isfinite(x) and x > 0 for x in (mean, sigma)):
        raise ValueError('Candidate has nonpositive/nonfinite mean or sigma.')
    return mean, sigma, True


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input-root', 'mainconf', 'calibration', 'tof-nuclear-config', 'settings', 'output-dir'):
        parser.add_argument('--' + name, type=Path, required=True)
    parser.add_argument('--candidate-calibration', '--analytic-coefficients', type=Path,
        dest='candidate_calibration', help='Overlay an explicit study-only log-polynomial candidate; never adopt it.')
    args = parser.parse_args()
    source = args.input_root.resolve(strict=True)
    mainconf = args.mainconf.resolve(strict=True)
    calibration_path = args.calibration.resolve(strict=True)
    nuclear_path = args.tof_nuclear_config.resolve(strict=True)
    settings_path = args.settings.resolve(strict=True)
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError('Choose a new output directory; existing results are preserved.')
    cfg = yaml.safe_load(settings_path.read_text())
    required_numbers = ('wideRigidityMin', 'wideRigidityMax', 'lowRigidityMin', 'lowRigidityMax',
        'dedxDisplayMin', 'dedxDisplayMax', 'curveRigidityStep', 'displayDedxRebinFactor',
        'projectionDedxRebinFactor', 'projectionRigidityMin', 'projectionRigidityMax',
        'projectionDisplayNSigmaMin', 'projectionDisplayNSigmaMax', 'curveNSigmaBand',
        'tofModeSearchNSigma', 'tofMomentNSigma', 'projectionLogMinimum',
        'tofShapeNormalizeNSigma', 'tofShapeHighMinNSigma', 'tofShapeHighMaxNSigma',
        'projectionLinearMaximumFactor', 'projectionLogMaximumFactor',
        'mapCanvasWidth', 'mapCanvasHeight', 'projectionCanvasWidth', 'projectionCanvasHeight',
        'shapeCanvasWidth', 'shapeCanvasHeight', 'axisFont',
        'mapMarginLeft', 'mapMarginRight', 'mapMarginBottom', 'mapMarginTop',
        'projectionMarginLeft', 'projectionMarginRight', 'projectionMarginBottom', 'projectionMarginTop',
        'mapAxisTitleSize', 'mapAxisLabelSize', 'projectionAxisTitleSize', 'projectionAxisLabelSize',
        'mapXTitleOffset', 'mapYTitleOffset', 'projectionXTitleOffset', 'projectionYTitleOffset',
        'mapLegendXMin', 'mapLegendXMax', 'mapLegendYMin', 'mapLegendYMax', 'mapLegendTextSize',
        'projectionLegendXMin', 'projectionLegendXMax', 'projectionLegendYMin', 'projectionLegendYMax',
        'projectionLegendTextSize', 'headerY', 'mapHeaderTextSize', 'projectionHeaderTextSize',
        'mapCaptionTextSize', 'projectionCaptionTextSize', 'mapCaptionFirstY', 'mapCaptionSecondY',
        'projectionCaptionFirstY', 'projectionCaptionSecondY', 'projectionCaptionThirdY')
    if cfg['schemaVersion'] != 1 or any(type(cfg[k]) not in (float, int) or not math.isfinite(cfg[k])
                                       for k in required_numbers):
        raise ValueError('Invalid explicit drawing settings.')
    for lo, hi in (('wideRigidityMin', 'wideRigidityMax'), ('lowRigidityMin', 'lowRigidityMax'),
                   ('dedxDisplayMin', 'dedxDisplayMax'), ('projectionRigidityMin', 'projectionRigidityMax'),
                   ('projectionDisplayNSigmaMin', 'projectionDisplayNSigmaMax')):
        if not cfg[lo] < cfg[hi]:
            raise ValueError('Reversed drawing interval: ' + lo)
    if not (cfg['wideRigidityMin'] > 0 and cfg['curveRigidityStep'] > 0 and cfg['curveNSigmaBand'] > 0
            and cfg['tofModeSearchNSigma'] > 0 and cfg['tofMomentNSigma'] > 0
            and cfg['tofShapeNormalizeNSigma'] > 0
            and 0 < cfg['tofShapeHighMinNSigma'] < cfg['tofShapeHighMaxNSigma']
            and cfg['projectionLogMinimum'] > 0 and cfg['projectionLinearMaximumFactor'] > 1
            and cfg['projectionLogMaximumFactor'] > 1):
        raise ValueError('Invalid scales or positive settings.')
    for key in ('displayDedxRebinFactor', 'projectionDedxRebinFactor'):
        if type(cfg[key]) is not int or cfg[key] <= 0:
            raise ValueError('Rebin factors must be positive integers.')
    for key in ('mapCanvasWidth', 'mapCanvasHeight', 'projectionCanvasWidth', 'projectionCanvasHeight',
                'shapeCanvasWidth', 'shapeCanvasHeight', 'axisFont'):
        if type(cfg[key]) is not int or cfg[key] <= 0:
            raise ValueError('Canvas dimensions and font IDs must be positive integers.')
    for prefix in ('map', 'projection'):
        for side in ('Left', 'Right', 'Bottom', 'Top'):
            if not 0 < cfg[prefix + 'Margin' + side] < 0.5:
                raise ValueError('Pad margin must be between zero and one half.')
        if cfg[prefix + 'MarginTop'] + cfg[prefix + 'MarginBottom'] >= 0.8:
            raise ValueError('Pad vertical plotting area is too small.')
        for label in ('X', 'Y'):
            if not 0 < cfg[prefix + 'Legend' + label + 'Min'] < cfg[prefix + 'Legend' + label + 'Max'] < 1:
                raise ValueError('Invalid explicit legend interval.')
        for label in ('XTitleOffset', 'YTitleOffset', 'AxisTitleSize', 'AxisLabelSize',
                      'LegendTextSize', 'HeaderTextSize', 'CaptionTextSize'):
            if cfg[prefix + label] <= 0:
                raise ValueError('Display text and offsets must be positive.')
    if not (0 < cfg['projectionCaptionThirdY'] < cfg['projectionCaptionSecondY']
            < cfg['projectionCaptionFirstY'] < cfg['projectionMarginBottom']
            and 0 < cfg['mapCaptionSecondY'] < cfg['mapCaptionFirstY'] < cfg['mapMarginBottom']
            and max(cfg['mapLegendYMax'], cfg['projectionLegendYMax']) < cfg['headerY'] < 1):
        raise ValueError('Invalid caption/header order.')
    if any(type(cfg[k]) is not bool for k in ('requireTofHistogram', 'logZ')):
        raise ValueError('Explicit boolean flags required.')
    main_doc = yaml.safe_load(mainconf.read_text())
    info_ref = main_doc['analysis']
    info_path = Path(info_ref) if Path(info_ref).is_absolute() else PROJECT/'config'/info_ref
    info = yaml.safe_load(info_path.read_text())
    tag = info['starTag']
    if (tag['collision'], tag['year'], tag['productionTag'], tag['libraryTag']) != (
            'auau3p85', 'year2021', 'P24iy', 'SL24y') or '2021P24iy' not in source.name:
        raise ValueError('This diagnostic only accepts explicit 3p85/2021P24iy/SL24y data/runtime.')
    calibration = yaml.safe_load(calibration_path.read_text())
    if calibration['source']['sourceRootFingerprint'] != fingerprint(source):
        raise ValueError('Calibration and raw ROOT fingerprints differ.')
    if digest(PROJECT/'StMaker/common/NuclearIdDeDxVsMom.h') != calibration['source']['legacyHeaderSha256']:
        raise ValueError('Legacy table changed since the curve export.')
    for key, sha in (('fitSummary', 'fitSummarySha256'), ('fitSource', 'fitSourceSha256'), ('settings', 'settingsSha256')):
        if digest(calibration['source'][key]) != calibration['source'][sha]:
            raise ValueError('Frozen calibration source hash differs: ' + key)
    module_path = HERE/'export_runtime_calibration.py'
    spec = importlib.util.spec_from_file_location('overlay_calibration', module_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    nuclear = yaml.safe_load(nuclear_path.read_text())
    cfg['tofMass2Min'] = nuclear['nuclearTofMeanHe3'] - nuclear['m2SigmaCut'] * nuclear['nuclearTofSigmaHe3']
    cfg['tofMass2Max'] = nuclear['nuclearTofMeanHe3'] + nuclear['m2SigmaCut'] * nuclear['nuclearTofSigmaHe3']
    if not (cfg['tofMass2Min'] > 0 and cfg['tofMass2Min'] < cfg['tofMass2Max']):
        raise ValueError('Invalid recorded TOF-tag window.')
    candidate = None
    if args.candidate_calibration:
        candidate_path = args.candidate_calibration.resolve(strict=True)
        candidate = yaml.safe_load(candidate_path.read_text())
        if (candidate['schemaVersion'] != 1 or candidate['functionalForm'] != 'log10_polynomial'
                or candidate['rigidityUnit'] != 'GeV/c' or candidate['dedxUnit'] != 'keV/cm'
                or candidate['outsideRangePolicy'] != 'legacy' or candidate['productionEnabled'] is not False
                or candidate['candidateOnly'] is not True):
            raise ValueError('Only explicit, unadopted logpoly candidates are overlaid.')
        for key in ('he3', 'he4'):
            curve = candidate['species'][key]
            if len(curve['meanCoefficients']) != 4 or len(curve['sigmaCoefficients']) != 3:
                raise ValueError('Candidate coefficient dimensions must be 4 and 3.')
            if not (len(curve['domain']) == 2 and 0 < curve['domain'][0] < curve['domain'][1]):
                raise ValueError('Invalid candidate domain.')
            if any(type(x) not in (float, int) or not math.isfinite(x)
                   for x in curve['meanCoefficients'] + curve['sigmaCoefficients']):
                raise ValueError('Candidate coefficients must be finite numbers.')
    initial = fingerprint(source)
    output.mkdir(parents=True, exist_ok=False)
    snapshot = output/'source'
    snapshot.mkdir()
    for path in (Path(__file__).resolve(), module_path, HERE/'DrawHeliumRawCalibration.C'):
        shutil.copy2(path, snapshot/path.name)
    shutil.copy2(settings_path, output/'settings.yaml')
    shutil.copy2(calibration_path, output/'current_calibration.yaml')
    shutil.copy2(nuclear_path, output/'tof_nuclear_source.yaml')
    if candidate:
        shutil.copy2(candidate_path, output/'candidate_calibration.yaml')
    numeric = {k: int(v) if type(v) is bool else v for k, v in cfg.items() if type(v) in (bool, int, float)}
    (output/'settings.tsv').write_text(''.join('%s\t%.17g\n' % (k, v) for k, v in numeric.items()))
    positions = set()
    count = int(math.ceil((cfg['wideRigidityMax'] - cfg['wideRigidityMin']) / cfg['curveRigidityStep']))
    if count > 100000:
        raise ValueError('Curve sampling grid is unreasonable.')
    for i in range(count + 1):
        positions.add(min(cfg['wideRigidityMax'], cfg['wideRigidityMin'] + i * cfg['curveRigidityStep']))
    for curve in calibration['species'].values():
        boundaries = list(curve['domain'])
        knots = curve['knots']
        for left, right in zip(knots, knots[1:]):
            if right['rigidity'] - left['rigidity'] > calibration['maxInterpolationGap'] + 1.e-12:
                boundaries += [left['rigidity'] + curve['supportBinHalfWidth'],
                               right['rigidity'] - curve['supportBinHalfWidth']]
        for point in boundaries:
            positions.update(x for x in (point - 1.e-8, point, point + 1.e-8)
                             if cfg['wideRigidityMin'] <= x <= cfg['wideRigidityMax'])
    with (output/'curves.tsv').open('w') as handle:
        for mode in range(3 if candidate else 2):
            for species in range(4):
                if mode and species < 2:
                    continue
                for r in sorted(positions):
                    if mode == 0:
                        mean, sigma = module.legacy_parameters(species, r)
                        custom = False
                    elif mode == 1:
                        mean, sigma, policy = module.calibrated_parameters(calibration, species, r)
                        custom = policy.startswith('calibrated')
                    else:
                        mean, sigma, custom = candidate_parameters(candidate, species, r, module.legacy_parameters)
                    if mean > 0 and sigma > 0:
                        handle.write('%d\t%d\t%.17g\t%.17g\t%.17g\t%d\n' %
                                     (mode, species, r, mean, sigma, int(custom)))
    with (output/'knots.tsv').open('w') as handle:
        for species, key in ((2, 'he3'), (3, 'he4')):
            for knot in calibration['species'][key]['knots']:
                handle.write('%d\t%.17g\t%.17g\t%.17g\t%d\n' % (species, knot['rigidity'],
                    knot['mean'], knot['sigma'], int(knot['provenance']['fitQuality'])))
    context = dict(input=initial, runtimeMainconf=str(mainconf), runtimeMainconfSha256=digest(mainconf),
        runtimeAnalysisInfo=str(info_path.resolve()), runtimeAnalysisInfoSha256=digest(info_path),
        currentCalibration=str(calibration_path), currentCalibrationSha256=digest(calibration_path),
        tofNuclearSource=str(nuclear_path), tofNuclearSourceSha256=digest(nuclear_path),
        settings=cfg, candidateCalibration=str(candidate_path) if candidate else None,
        candidateCalibrationSha256=digest(candidate_path) if candidate else None,
        scope='Read-only unselected charge-inclusive tracks; TOF-tag is a subset, not another independent sample.',
        limitations=['No signed-charge or centrality information in these TH2 histograms.',
                     'TOF selection is already applied by the source producer, not reapplied here.',
                     'Projection TOF mode/window moments are diagnostics, not Gaussian isotope calibration.',
                     'Input R within-bin density and censored high-dEdx identities are not recoverable.'])
    save(output/'context.json', context)
    macro = snapshot/'DrawHeliumRawCalibration.C'
    arguments = (source, output, output/'settings.tsv', output/'curves.tsv', output/'knots.tsv',
                 cfg['rawHistogram'], cfg['tofHistogram'])
    call = str(macro) + '++(' + ','.join(json.dumps(str(x)) for x in arguments) + ')'
    with (output/'root.log').open('w') as handle:
        result = subprocess.run(['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
            'root4star', '-l', '-b', '-q', call], cwd=str(PROJECT), stdout=handle, stderr=subprocess.STDOUT, timeout=600)
    if result.returncode or 'HELIUM_RAW_OVERLAY_COMPLETE' not in (output/'root.log').read_text():
        raise RuntimeError('Raw diagnostic failed; preserved logs: ' + str(output/'root.log'))
    if fingerprint(source) != initial:
        raise ValueError('Read-only input ROOT fingerprint changed.')
    with (output/'projection_audit.tsv').open() as handle:
        rows = [{key: float(value) for key, value in row.items()}
                for row in csv.DictReader(handle, delimiter='\t')]
    save(output/'tof3_moment_summary.json', dict(rows=rows,
        scope='TOF-tag conditional window mode/moments, not fitted Gaussian isotope parameters.',
        statisticalErrors='Delta-method counting-error proxies inside the fixed measured window; no TOF/model/selection systematic error.',
        normalization='Scaled tag matches ALL integral in tag-mode +- configured legacy-He3 sigma core; not a TOF efficiency or background coefficient.',
        shapeCaution='TOF sample is an ALL subset and may have dEdx/rigidity/acceptance bias; no equality of shapes or independent ALL+TAG likelihood is assumed.'))
    save(output/'completed.json', dict(phase='complete', inputPreserved=True,
        sourceSha256={p.name: digest(p) for p in snapshot.iterdir()},
        products={p.name: digest(p) for p in output.iterdir() if p.is_file() and p.name != 'completed.json'}))
    print('COMPLETE', output, flush=True)


if __name__ == '__main__':
    main()
