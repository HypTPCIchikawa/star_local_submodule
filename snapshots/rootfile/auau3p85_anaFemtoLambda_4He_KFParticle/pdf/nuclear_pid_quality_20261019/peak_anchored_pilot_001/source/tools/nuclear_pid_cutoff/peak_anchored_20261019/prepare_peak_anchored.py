#!/usr/bin/env python3
"""Prepare/run an isolated fixed-original-PID t+He3+He4 peak-anchored pilot."""
import argparse
import csv
import datetime
import hashlib
import importlib.util
import json
import math
from pathlib import Path
import shutil
import subprocess

import yaml

PROJECT = Path(__file__).resolve().parents[3]
SPECIES = ('d', 't', 'he3', 'he4')
SHAPE_HEADER = ('index', 'rBin', 'rLo', 'rHi', 'node', 'weight', 'rValue',
                'species', 'mean', 'sigma', 'supported', 'fitLo', 'fitHi')
GATE_HEADER = ('index', 'node', 'species', 'mean', 'sigma', 'active')
REGION_KEYS = ('roiLeftSigma', 'roiRightSigma', 'heRegionLeftSigma',
               'heRegionRightSigma', 'tControlSigma', 'fitWindowAlignmentKeV')


def require(condition, message):
    if not condition:
        raise ValueError(message)


def owned(value, exists=True):
    value = Path(value)
    result = (value if value.is_absolute() else PROJECT / value).resolve(strict=exists)
    require(PROJECT in result.parents, 'Require a project-owned study path: ' + str(value))
    return result


def sha(source):
    result = hashlib.sha256()
    with source.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            result.update(block)
    return result.hexdigest()


def read(source):
    return json.loads(source.read_text())


def write(destination, value):
    with destination.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def table(destination, header, values):
    with destination.open('x', newline='') as stream:
        output = csv.writer(stream, delimiter='\t', lineterminator='\n')
        if header:
            output.writerow(header)
        output.writerows(values)


def module(source):
    spec = importlib.util.spec_from_file_location('peak_anchor_original_pid_reader', str(source))
    result = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(result)
    return result


def now():
    return datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat()


def quote(text):
    text = str(text)
    require('"' not in text and '\\' not in text and '\n' not in text, 'Unsupported ROOT path quoting')
    return '"' + text + '"'


def validate_config(cfg):
    require(cfg['schemaVersion'] == 1 and cfg['responseStudyMode'] == 'peak_anchored_fixed_pid',
            'Wrong explicit peak-anchored study schema')
    require(cfg['shapePolicy'] == 'fixed_original_pid_functions' and cfg['freeSpecies'] == []
            and cfg['fitSpecies'] == ['t', 'he3', 'he4'], 'Nuclear shapes must remain fixed original PID functions')
    require(cfg['targetSpecies'] == 'he4' and cfg['sampleHistogram'] == 'ALL', 'Wrong target/cache sample')
    require(cfg['responseModels'] == ['gaussian', 'lognormal']
            and cfg['backgroundModels'] == ['linear', 'exponential'], 'Compare explicit four response/background models')
    require(cfg['adapterFunction'] == 'FitPeakAnchored'
            and Path(cfg['adapter']).name == 'FitPeakAnchored.C', 'Wrong explicit study adapter')
    for key in REGION_KEYS:
        require(math.isfinite(cfg[key]) and cfg[key] > 0, 'Invalid configured region factor: ' + key)
    require(0 < cfg['materialOmissionFraction'] < 1
            and cfg['axisCountBoundMethod'] == 'observed_plus_nsigma_sqrt'
            and math.isfinite(cfg['axisCountBoundNSigma']) and cfg['axisCountBoundNSigma'] > 0,
            'Invalid explicit descriptive omitted-component bound')
    require(cfg['numeric']['displayDedxRebinFactor'] == cfg['plots']['displayDedxRebinFactor'],
            'Fit/display bookkeeping rebin declarations differ')


def native_bin(rigidity, axis, tolerance):
    lo, hi = axis['rRange']
    width = axis['rBinWidth']
    index = int(round((rigidity - lo) / width - .5))
    center = lo + (index + .5) * width
    require(0 <= index < axis['rBins'] and abs(center - rigidity) < tolerance,
            'Requested rigidity is not an actual saved native-bin center')
    return index + 1, lo + index * width, lo + (index + 1) * width


def snapped(lo, hi, alignment, axis, tolerance):
    require(math.isfinite(lo) and math.isfinite(hi) and hi > lo, 'Bad physical region')
    first = math.floor(lo / alignment) * alignment
    last = math.ceil(hi / alignment) * alignment
    low, high = axis['dedxRange']
    require(first >= low - tolerance and last <= high + tolerance, 'Configured peak ROI leaves the saved dEdx axis')
    width = axis['dedxBinWidth']
    require(abs((first - low) / width - round((first - low) / width)) < tolerance
            and abs((last - low) / width - round((last - low) / width)) < tolerance,
            'ROI snap must retain adjacent whole native bins')
    return float(first), float(last)


def make_slices(cfg, runtime, reader, axis):
    tolerance = cfg['numeric']['axisTolerance']
    selected, shapes, gates = [], [], []
    seen = set()
    for index, rigidity in enumerate(cfg['representativeRigidity']):
        rbin, lo, hi = native_bin(rigidity, axis, tolerance)
        require(rbin not in seen, 'Repeated native rigidity bin')
        seen.add(rbin)
        nodes = []
        for q, (r, weight) in enumerate(reader.gl3_nodes(lo, hi)):
            values = [reader.legacy_parameters(runtime['tables'], s, r) for s in range(4)]
            require(all(math.isfinite(value) and value > 0 for pair in values for value in pair),
                    'Original PID function is unsupported/nonpositive at an R node')
            nodes.append(dict(node=q, r=r, weight=weight, values=values))
        def left(s, factor):
            return min(node['values'][s][0] - factor * node['values'][s][1] for node in nodes)
        def right(s, factor):
            return max(node['values'][s][0] + factor * node['values'][s][1] for node in nodes)
        fitlo, fithi = snapped(min(left(1, cfg['roiLeftSigma']), left(2, cfg['roiLeftSigma'])),
                              right(3, cfg['roiRightSigma']), cfg['fitWindowAlignmentKeV'], axis, tolerance)
        helium = snapped(left(2, cfg['heRegionLeftSigma']), right(3, cfg['heRegionRightSigma']),
                         cfg['fitWindowAlignmentKeV'], axis, tolerance)
        control = snapped(left(1, cfg['tControlSigma']), right(1, cfg['tControlSigma']),
                          cfg['fitWindowAlignmentKeV'], axis, tolerance)
        require(fitlo <= min(helium[0], control[0]) and fithi >= max(helium[1], control[1]),
                'Saved shared fit ROI does not cover both diagnostic peak regions')
        selected.append(dict(index=index, rBin=rbin, r=rigidity, lo=lo, hi=hi,
                             fitLo=fitlo, fitHi=fithi,
                             regions={'helium': dict(lo=helium[0], hi=helium[1]),
                                      't_control': dict(lo=control[0], hi=control[1])}))
        for node in nodes:
            for s, (mean, sigma) in enumerate(node['values']):
                # Shape support is a mathematical original-PID reference, not fit validation.
                shapes.append((index, rbin, lo, hi, node['node'], node['weight'], node['r'],
                               SPECIES[s], mean, sigma, 0, fitlo, fithi))
                active = int(reader.current_active(runtime, s, node['r']))
                gates.append((index, node['node'], SPECIES[s], mean, sigma, active))
    require(len(selected) == 7, 'This diagnostic pilot is explicitly limited to seven native R bins')
    return selected, shapes, gates


def inspect_inputs(cfg, settings, main):
    baseline_path = owned(cfg['baseline'])
    baseline = read(baseline_path)
    require(baseline['productionEnabled'] is False
            and Path(baseline['runtimeMainconf']).resolve() == main,
            'Require the frozen original-polynomial reference mainconf')
    protected = {str(baseline_path): sha(baseline_path)}
    for relative, digest in baseline['sourceSha256'].items():
        snapshot = baseline_path.parent / 'source' / relative
        current = owned(relative)
        require(sha(snapshot) == digest and sha(current) == digest,
                'Frozen original-PID source/current reference changed: ' + relative)
        protected[str(snapshot)] = digest
        protected[str(current)] = digest
    reader_path = owned(cfg['runtimeReader'])
    reader = module(reader_path)
    runtime = reader.load_current_selection(PROJECT, main)
    require(runtime['calibration'] is None and baseline['current']['calibration'] is None
            and runtime['tables'] == baseline['current']['tables']
            and runtime['cuts'] == baseline['current']['cuts'],
            'Pilot must use the exact frozen original log10-polynomial PID selector')
    main_doc = yaml.safe_load(main.read_text())
    analysis_path = owned(Path('config') / main_doc['analysis'])
    tag = yaml.safe_load(analysis_path.read_text())['starTag']
    require(tuple(tag[key] for key in ('collision', 'year', 'productionTag', 'libraryTag'))
            == ('auau3p85', 'year2021', 'P24iy', 'SL24y'), 'Require 2021 P24iy / SL24y only')
    inputs = owned(cfg['inputStudy'])
    done, context, audit = (read(inputs / name) for name in ('completed.json', 'context.json', 'native_audit.json'))
    require(done['executionComplete'] and done['protectedPostRunVerified']
            and '2021P24iy_partial_merge' in Path(context['inputFingerprint']['path']).name,
            'Require the completed 2021 P24iy native QA extraction')
    cache = inputs / 'cache.root'
    require(sha(cache) == done['products']['cache.root'] and audit['cacheCountErrorFlowPreserved']
            and audit['writtenCacheReadbackExact'], 'Native QA cache/count/error/flow changed')
    axis = audit['nativeAxis']
    require(axis['rBinWidth'] == cfg['numeric']['nativeRigidityBinWidth']
            and axis['dedxBinWidth'] == cfg['numeric']['nativeDedxBinWidth'], 'Configured native axis differs')
    for name in ('cache.root', 'completed.json', 'context.json', 'native_audit.json'):
        source = inputs / name
        digest = sha(source)
        if name != 'completed.json':
            require(digest == done['products'][name], 'Completed native extraction metadata changed: ' + name)
        protected[str(source)] = digest
    sources = {settings, Path(__file__).resolve(), reader_path, main, analysis_path,
               owned(cfg['legacyHeader']), owned(cfg['runtimeWrapper']), owned(cfg['adapter']),
               owned('script/analysis_info_helper.py')}
    sources.update(owned(value) for value in cfg['adapterDependencies'])
    sources.update(owned(Path('config') / value) for value in main_doc.values())
    return runtime, reader, axis, sources, protected


def prepare(settings, output, explicit_mainconf):
    cfg = yaml.safe_load(settings.read_text())
    validate_config(cfg)
    main = owned(explicit_mainconf)
    require(main == owned(cfg['runtimeMainconf']), 'Explicit runtime mainconf differs from settings')
    require(output.parent == owned(cfg['outputBase']) and not output.exists(), 'Choose a new direct study sibling output')
    runtime, reader, axis, sources, immutable = inspect_inputs(cfg, settings, main)
    slices, shapes, gates = make_slices(cfg, runtime, reader, axis)
    source_meta = {str(source): sha(source) for source in sources}
    output.mkdir()
    source_root = output / 'source'
    for original, digest in source_meta.items():
        original = Path(original)
        destination = source_root / original.relative_to(PROJECT)
        destination.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(original), str(destination))
        require(sha(original) == digest and sha(destination) == digest, 'Source changed during snapshot')
    jobs = []
    for response in cfg['responseModels']:
        for background in cfg['backgroundModels']:
            identity = response + '_' + background
            folder = output / identity
            folder.mkdir()
            numeric = dict(cfg['numeric'])
            numeric.update(responseModel=response, backgroundModel=background,
                           responseStudyMode=cfg['responseStudyMode'], shapePolicy=cfg['shapePolicy'],
                           fitSpecies=','.join(cfg['fitSpecies']), freeSpecies='', targetShapeMode='fixed', fixedShape=1,
                           sampleHistogram='ALL', targetSpecies='he4',
                           maxNSigma=runtime['cuts']['maxNSigmaNuclear'],
                           materialOmissionFraction=cfg['materialOmissionFraction'],
                           axisCountBoundMethod=cfg['axisCountBoundMethod'], axisCountBoundNSigma=cfg['axisCountBoundNSigma'])
            numeric.update({key: cfg[key] for key in REGION_KEYS})
            table(folder / 'settings.tsv', (), numeric.items())
            table(folder / 'shapes.tsv', SHAPE_HEADER, shapes)
            table(folder / 'gates.tsv', GATE_HEADER, gates)
            table(folder / 'regions.tsv', ('index', 'region', 'lo', 'hi'),
                  [(item['index'], name, region['lo'], region['hi'])
                   for item in slices for name, region in item['regions'].items()])
            macro = source_root / cfg['adapter']
            invocation = cfg['adapterFunction'] + '(' + ','.join(quote(value) for value in
                (owned(cfg['inputStudy']) / 'cache.root', folder, folder / 'settings.tsv',
                 folder / 'shapes.tsv', folder / 'gates.tsv')) + ');'
            bootstrap = 'void RunPeakAnchored(){\n'
            bootstrap += 'gSystem->AddIncludePath(' + quote('-I' + str(source_root / 'tools/nuclear_pid_tof')) + ');\n'
            bootstrap += 'Int_t error=0;gROOT->ProcessLine(' + quote('.L ' + str(macro) + '+') + ',&error);if(error){gSystem->Exit(2);return;}\n'
            bootstrap += 'gROOT->ProcessLine(' + json.dumps(invocation) + ',&error);if(error)gSystem->Exit(2);\n}\n'
            with (folder / 'RunPeakAnchored.C').open('x') as stream:
                stream.write(bootstrap)
            for source in folder.iterdir():
                immutable[str(source)] = sha(source)
            jobs.append(dict(id=identity, response=response, background=background, target='he4',
                             rows=len(slices), slices=slices, folder=str(folder)))
    immutable.update({str(source): sha(source) for source in source_root.rglob('*') if source.is_file()})
    write(output / 'plan.json', dict(schemaVersion=1, preparedJST=now(), config=cfg,
         explicitMainconf=str(main), jobs=jobs, immutableInputs=immutable, sourceInputs=source_meta,
         productionEnabled=False, fixedShapes=True, shapePolicy=cfg['shapePolicy'],
         rawGatePolicy='Recomputed original-polynomial best-species gates; diagnostic counts only',
         originalPidBaseline=str(owned(cfg['baseline'])), fitsExecuted=False,
         totalJobs=len(jobs), totalFits=sum(job['rows'] for job in jobs),
         limitations=cfg['studyLimitations']))
    print('PEAK_ANCHORED_PREPARED jobs=4 fits=28', output, flush=True)


def unchanged(plan):
    for source, digest in plan['immutableInputs'].items():
        require(sha(Path(source)) == digest, 'Immutable study input changed: ' + source)


def run(output):
    plan = read(output / 'plan.json')
    unchanged(plan)
    for job in plan['jobs']:
        folder = Path(job['folder'])
        if (folder / 'completed.json').exists():
            done = read(folder / 'completed.json')
            require(done['executionComplete'], 'Invalid prior completion')
            for name, digest in done['products'].items():
                require(sha(folder / name) == digest, 'Completed study output changed')
            continue
        require(not (folder / 'run.log').exists(), 'Attempted/failed job is preserved; choose a new study sibling')
        wrapper, main = owned(plan['config']['runtimeWrapper']), owned(plan['explicitMainconf'])
        require(sha(wrapper) == plan['sourceInputs'][str(wrapper)]
                and sha(main) == plan['sourceInputs'][str(main)], 'Runtime main/wrapper changed')
        command = ['bash', str(wrapper), str(main), 'root4star', '-l', '-b', '-q', str(folder / 'RunPeakAnchored.C')]
        write(folder / 'command.json', dict(cwd=str(PROJECT), argv=command))
        with (folder / 'run.log').open('x') as stream:
            subprocess.run(command, cwd=str(PROJECT), stdout=stream, stderr=subprocess.STDOUT,
                           check=True, timeout=plan['config']['executionTimeoutSeconds'])
        log = (folder / 'run.log').read_text()
        require('PEAK_ANCHORED_COMPLETE rows=' + str(job['rows']) in log
                and 'NUCLEAR_CUTOFF_READBACK_PASS rows=' + str(job['rows']) in log
                and 'PEAK_ANCHORED_DIAGNOSTIC_READBACK_PASS rows=' + str(job['rows']) in log,
                'Missing study/closed-ROOT readback completion')
        result = read(folder / 'cutoff_fits.json')
        require(result['productionEnabled'] is False and len(result['slices']) == job['rows'], 'Wrong pilot output identity')
        unchanged(plan)
        write(folder / 'completed.json', dict(executionComplete=True, rawReadback=True, productionEnabled=False,
             fixedShapes=True, completedJST=now(), protectedInputsUnchanged=True,
             products={source.name: sha(source) for source in folder.iterdir() if source.is_file()}))
        print('PEAK_ANCHORED_JOB_COMPLETE', job['id'], flush=True)


def check(output):
    plan = read(output / 'plan.json')
    unchanged(plan)
    total = 0
    for job in plan['jobs']:
        folder = Path(job['folder'])
        done = read(folder / 'completed.json')
        require(done['executionComplete'] and done['rawReadback'], 'Unfinished peak-anchored pilot')
        for name, digest in done['products'].items():
            require(sha(folder / name) == digest, 'Completed product changed')
        result = read(folder / 'cutoff_fits.json')
        require(result['productionEnabled'] is False and result['responseModel'] == job['response']
                and result['backgroundModel'] == job['background']
                and result['shapePolicy'] == 'fixed_original_pid_functions', 'Wrong result policy/model')
        for row, expected in zip(result['slices'], job['slices']):
            require(all(row[key] == expected[key] for key in ('index', 'rBin', 'lo', 'hi', 'fitLo', 'fitHi')),
                    'Pilot changed the saved native R/shared ROI')
            for component in row['components'].values():
                require(component['mean'] == component['initialMean']
                        and component['sigma'] == component['initialSigma']
                        and component.get('freeMean', False) is False
                        and component.get('freeSigma', False) is False, 'A nuclear shape was floated')
                require(component['selection']['available'] is False
                        and component['selection']['purity'] is None, 'This diagnostic pilot must not adopt purity')
        require(len(result['slices']) == job['rows'], 'Wrong completed row count')
        total += len(result['slices'])
    require(len(plan['jobs']) == 4 and total == 28, 'Wrong pilot scope')
    write(output / 'completed.json', dict(status='completed', completedJST=now(),
         productionEnabled=False, fixedShapes=True, immutableInputsUnchanged=True,
         totalJobs=4, totalFits=28, originalPidBaseline=plan['originalPidBaseline'],
         sourceAndCompiledArtifacts={str(source.relative_to(output)): sha(source)
             for source in (output / 'source').rglob('*') if source.is_file()},
         limitations=plan['limitations']))
    print('PEAK_ANCHORED_ALL_COMPLETE fits=28', flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'run', 'check'))
    parser.add_argument('--settings')
    parser.add_argument('--mainconf')
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    output = owned(args.output_dir, exists=args.action != 'prepare')
    if args.action == 'prepare':
        require(args.settings and args.mainconf, 'Prepare requires explicit settings and mainconf')
        prepare(owned(args.settings), output, args.mainconf)
    elif args.action == 'run':
        run(output)
    else:
        check(output)


if __name__ == '__main__':
    main()
