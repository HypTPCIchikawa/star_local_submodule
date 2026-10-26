#!/usr/bin/env python3
"""Read-only D/4He-compatible TOF control. No Gaussian calibration or PID edits."""
import argparse
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]


def digest(path):
    h = hashlib.sha256()
    with path.open('rb') as source:
        for block in iter(lambda: source.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fingerprint(path):
    s = path.stat()
    return dict(path=str(path), size=s.st_size, mtime_ns=s.st_mtime_ns, device=s.st_dev, inode=s.st_ino)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input-root', 'mainconf', 'tof-nuclear-config', 'old-study', 'settings', 'output-dir'):
        parser.add_argument('--' + name, type=Path, required=True)
    args = parser.parse_args()
    source = args.input_root.resolve(strict=True)
    mainconf = args.mainconf.resolve(strict=True)
    nuclear_path = args.tof_nuclear_config.resolve(strict=True)
    old = args.old_study.resolve(strict=True)
    settings = args.settings.resolve(strict=True)
    output = args.output_dir.resolve()
    if output.exists():
        raise ValueError('Choose a new output directory; never overwrite existing controls.')
    cfg = yaml.safe_load(settings.read_text())
    if cfg['schemaVersion'] != 1 or any(isinstance(v, (int, float)) and not math.isfinite(v) for v in cfg.values()):
        raise ValueError('Invalid explicit control settings.')
    for lo, hi in (('rigidityMin', 'rigidityMax'), ('fullDedxMin', 'fullDedxMax'), ('nuclearDedxMin', 'nuclearDedxMax')):
        if not cfg[lo] < cfg[hi]:
            raise ValueError('Reversed diagnostic interval: ' + lo)
    runtime = yaml.safe_load(mainconf.read_text())
    info_path = (PROJECT / 'config' / runtime['analysis']).resolve(strict=True)
    info = yaml.safe_load(info_path.read_text())
    tag = info['starTag']
    if (tag['collision'], tag['year'], tag['productionTag'], tag['libraryTag']) != ('auau3p85', 'year2021', 'P24iy', 'SL24y') or '2021P24iy' not in source.name:
        raise ValueError('This control requires explicit 3p85/2021P24iy/SL24y input/runtime.')
    completed = json.loads((old / 'completed.json').read_text())
    protected = {p: digest(p) for p in (old / 'purity.json', old / 'he4_fits.root', old / 'context.json')}
    for path, sha in protected.items():
        if completed['products'][path.name] != sha:
            raise ValueError('Frozen old fit SHA differs: ' + str(path))
    summary = json.loads((old / 'purity.json').read_text())
    initial = fingerprint(source)
    if Path(summary['inputRoot']).resolve() != source:
        raise ValueError('Reference fit and original ROOT differ.')
    nuclear = yaml.safe_load(nuclear_path.read_text())
    cfg['tag3Mass2Min'] = nuclear['nuclearTofMeanHe3'] - nuclear['m2SigmaCut'] * nuclear['nuclearTofSigmaHe3']
    cfg['tag3Mass2Max'] = nuclear['nuclearTofMeanHe3'] + nuclear['m2SigmaCut'] * nuclear['nuclearTofSigmaHe3']
    cfg['tagDMass2Min'] = nuclear['nuclearTofMeanDeuteron'] - nuclear['m2SigmaCut'] * nuclear['nuclearTofSigmaDeuteron']
    cfg['tagDMass2Max'] = nuclear['nuclearTofMeanDeuteron'] + nuclear['m2SigmaCut'] * nuclear['nuclearTofSigmaDeuteron']
    cfg['he4MassOverChargeSquared'] = (nuclear['nuclearMassHe4'] / 2.)**2
    if not cfg['tagDMass2Min'] < cfg['he4MassOverChargeSquared'] < cfg['tagDMass2Max']:
        raise ValueError('The existing D window is not compatible with expected He4 mass/charge.')
    helper_path = HERE.parent / 'export_runtime_calibration.py'
    spec = importlib.util.spec_from_file_location('tof_control_legacy', helper_path)
    helper = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(helper)
    selected = [(i, s) for i, s in enumerate(summary['slices'])
                if s['lo'] >= cfg['rigidityMin'] - cfg['axisRelativeTolerance'] and s['hi'] <= cfg['rigidityMax'] + cfg['axisRelativeTolerance']]
    if not selected:
        raise ValueError('No requested frozen R bins.')
    references = []
    for i, s in selected:
        fits = [(j, f) for j, f in enumerate(s['fits']) if len(f.get('parameters', [])) >= 11 and all(math.isfinite(p) for p in f['parameters'])]
        if not fits:
            raise ValueError('No stored T/3He/4He reference model for R slice ' + str(i))
        chosen = next((item for item in fits if item[0] == s['chosen']), min(fits, key=lambda item: item[1]['aic']))
        j, f = chosen
        means = {3: s['calibrationMean'] + s['calibrationSigma'] * f['parameters'][1]}
        for k, species in enumerate(f['neighborSpecies']):
            means[species] = s['calibrationMean'] + s['calibrationSigma'] * f['parameters'][6 + 3*k]
        if 1 not in means or 2 not in means:
            raise ValueError('Reference model lacks T or 3He.')
        legacy = [helper.legacy_parameters(species, s['r'])[0] for species in range(4)]
        references.append([i, s['lo'], s['hi'], s['r'], int(j == s['chosen']), int(f['quality'])] + legacy + [means[1], means[2], means[3]])
    output.mkdir(parents=True, exist_ok=False)
    snapshot = output / 'source'
    snapshot.mkdir()
    for path in (Path(__file__).resolve(), HERE / 'DrawTofCompatibleHelium.C', helper_path,
                 PROJECT / 'StMaker/common/FemtoLambdaLegacy.cxx'):
        shutil.copy2(path, snapshot / path.name)
    shutil.copy2(settings, output / 'settings.yaml')
    shutil.copy2(nuclear_path, output / 'tof_nuclear_source.yaml')
    (output / 'settings.tsv').write_text(''.join('%s\t%.17g\n' % (k, v) for k, v in cfg.items() if isinstance(v, (int, float))))
    (output / 'references.tsv').write_text(''.join('\t'.join('%.17g' % v for v in row) + '\n' for row in references))
    context = dict(input=initial, runtimeMainconf=str(mainconf), runtimeMainconfSha256=digest(mainconf),
        runtimeAnalysisInfo=str(info_path), runtimeAnalysisInfoSha256=digest(info_path),
        tofNuclearSource=str(nuclear_path), tofNuclearSourceSha256=digest(nuclear_path), settings=cfg,
        protectedOldFits={str(p): sha for p, sha in protected.items()}, productionChanged=False,
        scope='Unselected charge-inclusive pre-PID tracks. TOF D gate is D/4He-compatible, NOT a pure isotope tag.',
        limitations=['TOF samples are correlated subsets of ALL, with acceptance, efficiency, and impurity biases.',
            'Window mode/mean/RMS are diagnostics, NOT fitted Gaussian calibration parameters.',
            'Reference old free-fit centers are initialization/model labels, NOT isotope truth.',
            'No signed-charge or centrality conditioning is possible with these saved histograms.'])
    (output / 'context.json').write_text(json.dumps(context, indent=2, allow_nan=False) + '\n')
    macro = snapshot / 'DrawTofCompatibleHelium.C'
    params = (source, output, output / 'settings.tsv', output / 'references.tsv', cfg['rawHistogram'], cfg['tag3Histogram'], cfg['tagDHistogram'])
    call = str(macro) + '++(' + ','.join(json.dumps(str(v)) for v in params) + ')'
    with (output / 'root.log').open('w') as log:
        result = subprocess.run(['bash', str(PROJECT / 'script/singularity_root_batch.sh'), str(mainconf), 'root4star', '-l', '-b', '-q', call], cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT, timeout=cfg['executionTimeoutSeconds'])
    if result.returncode or 'TOF_COMPATIBLE_HELIUM_COMPLETE' not in (output / 'root.log').read_text():
        raise RuntimeError('Control drawing failed; see preserved root.log.')
    if fingerprint(source) != initial or any(digest(p) != sha for p, sha in protected.items()):
        raise ValueError('Original input or frozen reference changed.')
    done = dict(phase='complete', inputPreserved=True, productionChanged=False,
        sourceSha256={p.name: digest(p) for p in snapshot.iterdir() if p.is_file()},
        products={p.name: digest(p) for p in output.iterdir() if p.is_file()})
    (output / 'completed.json').write_text(json.dumps(done, indent=2) + '\n')
    print('TOF_COMPATIBLE_CONTROL_COMPLETE', output, flush=True)


if __name__ == '__main__':
    main()
