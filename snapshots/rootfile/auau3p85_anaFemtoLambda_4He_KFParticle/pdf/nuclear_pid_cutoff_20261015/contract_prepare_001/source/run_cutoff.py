#!/usr/bin/env python3
"""Immutable native-R ALL fits for the TPC nuclear-PID cutoff study."""
import argparse
from datetime import datetime, timezone
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
SPECIES = ('d', 't', 'he3', 'he4')
TAG = dict(d='TAG_D', t='TAG_T', he3='TAG_HE3', he4='TAG_D')
RESPONSES = ('gaussian', 'lognormal')
BACKGROUNDS = ('linear', 'exponential')


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024*1024), b''):
            value.update(block)
    return value.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    Path(path).write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def resolve(value):
    path = Path(value)
    path = (path if path.is_absolute() else PROJECT/path).resolve(strict=True)
    if PROJECT not in path.parents:
        raise ValueError('Expected project-owned input: '+str(path))
    return path


def module(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    value = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(value)
    return value


def fingerprint(path):
    stat = path.stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                device=stat.st_dev, inode=stat.st_ino)


def positive(value):
    return isinstance(value, (int, float)) and math.isfinite(value) and value > 0


def tag_index(folder):
    completed = read(folder/'completed.json')
    if not completed.get('executionComplete') or completed.get('productionEnabled') is not False:
        raise ValueError('Incomplete TAG source.')
    signatures = {}
    for name in ('tag_fits.json', 'tag_fits.root'):
        path = folder/name
        if digest(path) != completed['products'][name]:
            raise ValueError('Changed TAG source: '+str(path))
        signatures[str(path)] = completed['products'][name]
    signatures[str(folder/'completed.json')] = digest(folder/'completed.json')
    index = {}
    for row in read(folder/'tag_fits.json')['slices']:
        if row['firstNativeRBin'] != row['lastNativeRBin']:
            raise ValueError('TAG shapes must use one native R bin.')
        key = row['sample'], row['firstNativeRBin']
        if key in index:
            raise ValueError('Duplicate TAG bin.')
        index[key] = row
    return index, signatures


def shape_failures(row, c, cfg):
    failures = []
    if not row.get('numericalConverged'):
        failures.append('tag_numerical_unconverged')
    if row.get('atBoundary'):
        failures.append('tag_shape_boundary')
    if c.get('zeroYield') or not c.get('active'):
        failures.append('tag_zero_or_inactive')
    checks = [('minWindowYield', c['windowYield'], False),
              ('minEffective', c['effectiveCount'], False),
              ('minDropDeltaD', c['conditionalDropDeltaD'], False),
              ('minDetectorAcceptance', c['detectorAcceptance'], False),
              ('maxMeanErrSigma', c['meanError']/c['sigma'], True),
              ('maxSigmaRelErr', c['sigmaError']/c['sigma'], True)]
    for key, value, maximum in checks:
        if not math.isfinite(value) or value < 0 or (value > cfg[key] if maximum else value < cfg[key]):
            failures.append(key)
    # Global chi2 is deliberately reported, not a sole identity veto.
    return failures


def varied(mean, sigma, name, cfg):
    if name not in cfg['shapeVariations']:
        raise ValueError('Undeclared variation.')
    if name == 'nominal':
        return mean, sigma
    if name in ('mean_minus', 'mean_plus'):
        return mean+(-1 if name == 'mean_minus' else 1)*cfg['meanShiftSigma']*sigma, sigma
    if name in ('width_minus', 'width_plus'):
        return mean, sigma*(1+(-1 if name == 'width_minus' else 1)*cfg['widthRelativeShift'])
    raise ValueError('Unknown variation.')


def tsv(path, header, rows):
    with path.open('w') as stream:
        if header:
            stream.write('\t'.join(header.split())+'\n')
        for row in rows:
            stream.write('\t'.join(format(x, '.17g') if isinstance(x, float) else str(x) for x in row)+'\n')


def prepare(args):
    settings = resolve(args.settings)
    cfg = yaml.safe_load(settings.read_text())
    required = ('schemaVersion', 'runtimeMainconf', 'inputStudy', 'outputBase', 'tagStudies',
                'rigidityMin', 'rigidityMax', 'nativeWidth', 'numeric', 'allFit', 'support',
                'executionTimeoutSeconds', 'shapeVariations')
    if cfg.get('schemaVersion') != 1 or any(k not in cfg for k in required):
        raise ValueError('Missing explicit cutoff settings.')
    mainconf = resolve(args.mainconf)
    if mainconf != resolve(cfg['runtimeMainconf']):
        raise ValueError('Runtime/mainconf mismatch.')
    baseline_path = resolve(args.baseline)
    baseline = read(baseline_path)
    if baseline.get('schemaVersion') != 1 or baseline.get('productionEnabled') is not False:
        raise ValueError('Invalid frozen baseline.')
    protected = {str(baseline_path): digest(baseline_path)}
    changes = []
    for rel, expected in baseline['sourceSha256'].items():
        frozen, live = baseline_path.parent/'source'/rel, PROJECT/rel
        if digest(frozen) != expected:
            raise ValueError('Baseline snapshot changed: '+str(frozen))
        protected[str(frozen)] = expected
        actual = digest(live) if live.exists() else None
        if actual != expected:
            changes.append(dict(path=str(live), baselineSha256=expected, currentSha256=actual))
    runtime = module('frozen_tpc_runtime', baseline_path.parent/'source/tools/nuclear_pid_tof/runtime_selection.py')
    current = baseline['current']
    if current['cuts']['m2_selection'] or current['cuts']['maxNSigmaNuclear'] != cfg['allFit']['maxNSigma']:
        raise ValueError('Reference must reproduce the actual TPC-only gate.')
    inputs = resolve(cfg['inputStudy'])
    completed, context = read(inputs/'completed.json'), read(inputs/'context.json')
    if not completed.get('executionComplete') or not completed.get('protectedPostRunVerified'):
        raise ValueError('Incomplete input extraction.')
    cache = inputs/'cache.root'
    if digest(cache) != completed['products']['cache.root']:
        raise ValueError('Changed immutable cache.')
    protected[str(cache)] = completed['products']['cache.root']
    for name in ('completed.json', 'context.json', 'native_audit.json'):
        protected[str(inputs/name)] = digest(inputs/name)
    original = Path(context['inputFingerprint']['path'])
    if '2021P24iy_partial_merge' not in original.name or fingerprint(original) != context['inputFingerprint']:
        raise ValueError('Original 2021 P24iy input changed.')
    audit = read(inputs/'native_audit.json')
    if not all(audit.get(k) for k in ('axesExactlyEqual', 'cacheCountErrorFlowPreserved', 'writtenCacheReadbackExact')):
        raise ValueError('Native cache audit failed.')
    axis, tolerance = audit['nativeAxis'], cfg['numeric']['axisTolerance']
    if abs(axis['rBinWidth']-cfg['nativeWidth']) > tolerance:
        raise ValueError('Native/configured R width differs.')
    start, stop = round(cfg['rigidityMin']/cfg['nativeWidth'])+1, round(cfg['rigidityMax']/cfg['nativeWidth'])
    if abs((start-1)*cfg['nativeWidth']-cfg['rigidityMin']) > tolerance or abs(stop*cfg['nativeWidth']-cfg['rigidityMax']) > tolerance:
        raise ValueError('R edges are not native-aligned.')
    sources = {}
    for response in RESPONSES:
        for background in BACKGROUNDS:
            sources[response, background], signatures = tag_index(resolve(cfg['tagStudies'][response][background]))
            protected.update(signatures)
    output = args.output_dir.resolve()
    if output.exists() or output.parent != resolve(cfg['outputBase']):
        raise ValueError('Require a new sibling output below outputBase.')
    output.mkdir()
    snapshot = output/'source'
    snapshot.mkdir()
    for source in (Path(__file__).resolve(), settings, HERE/'FitNuclearCutoff.C',
                   PROJECT/'tools/nuclear_pid_tof/NuclearFitCommon.h',
                   PROJECT/'tools/nuclear_pid_tof/NuclearLognormalFitCommon.h',
                   PROJECT/'script/singularity_root_batch.sh'):
        expected = digest(source)
        shutil.copy2(source, snapshot/source.name)
        if digest(snapshot/source.name) != expected:
            raise ValueError('Changed source while snapshotting.')
        protected[str(snapshot/source.name)] = expected
    shapes, gates, diagnostics, windows = [], [], [], []
    for index, rbin in enumerate(range(start, stop+1)):
        lo, hi = (rbin-1)*cfg['nativeWidth'], rbin*cfg['nativeWidth']
        references, selected, failures = [], {}, {}
        for species in SPECIES:
            for response in RESPONSES:
                for background in BACKGROUNDS:
                    row = sources[response, background][TAG[species], rbin]
                    c = row['components'][species]
                    if not positive(c['mean']) or not positive(c['sigma']):
                        raise ValueError('No invented fallback for invalid TAG shape.')
                    if (response, background) == (args.response, args.background):
                        selected[species] = c
                        failures[species] = shape_failures(row, c, cfg['support'])
                    references.extend(varied(c['mean'], c['sigma'], v, cfg) for v in cfg['shapeVariations'])
            a = sources[args.response, 'linear'][TAG[species], rbin]['components'][species]
            b = sources[args.response, 'exponential'][TAG[species], rbin]['components'][species]
            scale = max(a['sigma'], b['sigma'])
            if abs(a['mean']-b['mean'])/scale > cfg['support']['maxBackgroundMeanShiftSigma']:
                failures[species].append('background_mean_ambiguity')
            if abs(a['sigma']-b['sigma'])/scale > cfg['support']['maxBackgroundSigmaRelativeShift']:
                failures[species].append('background_width_ambiguity')
        for node, (r, weight) in enumerate(runtime.gl3_nodes(lo, hi)):
            for s, species in enumerate(SPECIES):
                mean, sigma, policy = runtime.current_parameters(current, s, r)
                active = runtime.current_active(current, s, r)
                gates.append([index, node, species, mean, sigma, int(active)])
                if active:
                    references.append((mean, sigma))
        half, alignment = cfg['allFit']['fitWindowSigmaHalfWidth'], cfg['allFit']['fitWindowAlignmentKeV']
        ylo, yhi = axis['dedxRange']
        intervals = [(max(ylo, m-half*w), min(yhi, m+half*w)) for m, w in references]
        intervals = [(a, b) for a, b in intervals if b > a]
        if not intervals:
            raise ValueError('No observable response window.')
        fitlo = max(ylo, math.floor(min(a for a, b in intervals)/alignment)*alignment)
        fithi = min(yhi, math.ceil(max(b for a, b in intervals)/alignment)*alignment)
        windows.append(dict(index=index, rBin=rbin, lo=lo, hi=hi, fitLo=fitlo, fitHi=fithi))
        for node, (r, weight) in enumerate(runtime.gl3_nodes(lo, hi)):
            for species in SPECIES:
                c = selected[species]
                mean, sigma = varied(c['mean'], c['sigma'], args.variation, cfg)
                shapes.append([index, rbin, lo, hi, node, weight, r, species, mean, sigma,
                               int(not failures[species]), fitlo, fithi])
        diagnostics.append(dict(index=index, rBin=rbin, lo=lo, hi=hi, failures=failures, tagComponents=selected))
    return cfg, mainconf, baseline_path, cache, original, context, changes, protected, output, snapshot, shapes, gates, diagnostics, windows


def run(args):
    (cfg, mainconf, baseline_path, cache, original, input_context, changes, protected,
     output, snapshot, shapes, gates, diagnostics, windows) = prepare(args)
    tsv(output/'shapes.tsv', 'index rBin rLo rHi node weight rValue species mean sigma supported fitLo fitHi', shapes)
    tsv(output/'gates.tsv', 'index node species mean sigma active', gates)
    numeric = dict(cfg['numeric'])
    numeric.update(responseModel=args.response, backgroundModel=args.background,
                   sampleHistogram='ALL', maxNSigma=cfg['allFit']['maxNSigma'],
                   minimumPurityDenominator=cfg['allFit']['minimumPurityDenominator'],
                   backgroundStartSlope=cfg['allFit']['backgroundStartSlope'])
    tsv(output/'settings.tsv', '', sorted(numeric.items()))
    save(output/'tag_shape_audit.json', diagnostics)
    save(output/'fixed_windows.json', windows)
    save(output/'effective_settings.json', cfg)
    save(output/'context.json', dict(schemaVersion=1, productionEnabled=False,
        createdUtc=datetime.now(timezone.utc).isoformat(), response=args.response,
        background=args.background, variation=args.variation, baseline=str(baseline_path),
        inputFingerprint=input_context['inputFingerprint'], runtimeMainconf=str(mainconf),
        cache=str(cache), protectedSources=protected, authorizedLiveReferenceDifferences=changes,
        policies=['ALL only; TAG never added to its likelihood.',
                  'Fixed projected native-bin TAG response is not certified intrinsic shape.',
                  'Gaussian mean/std and lognormal mode/curvature retain separate meanings.',
                  'No new curve fit, no Gaussian variance subtraction.',
                  'Frozen baseline is authoritative; subsequent authorized adapter edits are recorded.',
                  'Inclusive track purity is not automatically pair-conditioned CF purity.']))
    expression = str(snapshot/'FitNuclearCutoff.C')+'+('+','.join(json.dumps(str(p)) for p in
        (cache, output, output/'settings.tsv', output/'shapes.tsv', output/'gates.tsv'))+')'
    command = ['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
               'root4star', '-l', '-b', '-q', expression]
    save(output/'command.json', command)
    if args.prepare_only:
        save(output/'prepared.json', dict(prepared=True, executionComplete=False, nativeRows=len(diagnostics)))
        print('NUCLEAR_CUTOFF_PREPARED', output, flush=True)
        return
    print('NUCLEAR_CUTOFF_RUN', output, flush=True)
    try:
        with (output/'fit.log').open('w') as stream:
            process = subprocess.run(command, cwd=PROJECT, stdout=stream, stderr=subprocess.STDOUT,
                                     timeout=cfg['executionTimeoutSeconds'])
        if process.returncode or 'NUCLEAR_CUTOFF_COMPLETE' not in (output/'fit.log').read_text():
            raise RuntimeError('ROOT failed or completion marker absent; see retained fit.log.')
        result = read(output/'cutoff_fits.json')
        if len(result['slices']) != len(diagnostics):
            raise ValueError('Incomplete native R output.')
        for path, expected in protected.items():
            if digest(path) != expected:
                raise ValueError('Protected snapshot/cache/TAG source changed during run: '+path)
        if fingerprint(original) != input_context['inputFingerprint']:
            raise ValueError('Original ROOT fingerprint changed.')
        products = {str(p.relative_to(output)): digest(p) for p in output.rglob('*')
                    if p.is_file() and p.suffix not in ('.so', '.d', '.pcm')}
        save(output/'completed.json', dict(executionComplete=True, productionEnabled=False,
            physicsValidationComplete=False, sourceInputPreserved=True,
            protectedSnapshotsVerified=True, response=args.response, background=args.background,
            variation=args.variation, products=products, nativeRows=len(diagnostics),
            note='Execution completion is not isotope validation or final cutoff adoption.'))
        print('NUCLEAR_CUTOFF_RUN_COMPLETE', output, flush=True)
    except Exception as error:
        save(output/'failed.json', dict(executionComplete=False, productionEnabled=False, error=str(error)))
        raise


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('settings', 'baseline', 'mainconf', 'output-dir'):
        parser.add_argument('--'+name, required=True, type=Path)
    parser.add_argument('--response', required=True, choices=RESPONSES)
    parser.add_argument('--background', required=True, choices=BACKGROUNDS)
    parser.add_argument('--variation', required=True)
    parser.add_argument('--prepare-only', action='store_true')
    args = parser.parse_args()
    run(args)


if __name__ == '__main__':
    main()
