#!/usr/bin/env python3
"""既存native cacheをRだけ厳密5bin加算する。元ROOTはstatのみ、fitなし。"""
import argparse
import datetime
import hashlib
import json
from pathlib import Path
import shutil
import subprocess

import yaml

PROJECT = Path(__file__).resolve().parents[3]
HERE = Path(__file__).resolve().parent


def digest(path):
    value = hashlib.sha256()
    with path.open('rb') as stream:
        for chunk in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(chunk)
    return value.hexdigest()


def fingerprint(path):
    st = path.stat()
    return dict(path=str(path), size=st.st_size, mtime_ns=st.st_mtime_ns,
                device=st.st_dev, inode=st.st_ino)


def load(path):
    return json.loads(path.read_text())


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def resolve(value):
    path = Path(value)
    return (path if path.is_absolute() else PROJECT / path).resolve()


def require(condition, reason):
    if not condition:
        raise ValueError(reason)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    for name in ('input-study', 'mainconf', 'output-dir'):
        parser.add_argument('--' + name, required=True, type=Path)
    parser.add_argument('--factor', required=True, type=int)
    parser.add_argument('--execution-timeout', required=True, type=int)
    args = parser.parse_args()
    require(args.factor == 5 and args.execution_timeout > 0, 'Declare exact grouping factor 5 and a positive timeout')
    inputs = resolve(args.input_study)
    expected = resolve('rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_tof_20261014/inputs_001')
    require(inputs == expected, 'Only the immutable authorized native cache is supported')
    done, context, audit = (load(inputs / name) for name in ('completed.json', 'context.json', 'native_audit.json'))
    require(done['executionComplete'] and done['protectedPostRunVerified'] and audit['writtenCacheReadbackExact'],
            'Original native extraction/audit is incomplete')
    original = Path(context['inputFingerprint']['path'])
    require('2021P24iy_partial_merge' in original.name and fingerprint(original) == context['inputFingerprint'],
            'Original 2021 P24iy ROOT fingerprint differs; no alternate dataset is used')
    protected = {}
    for name in ('cache.root', 'context.json', 'native_audit.json'):
        path = inputs / name
        require(digest(path) == done['products'][name], 'Changed immutable native artifact: ' + name)
        protected[str(path)] = done['products'][name]
    protected[str(inputs / 'completed.json')] = digest(inputs / 'completed.json')
    axis = audit['nativeAxis']
    require(axis == dict(rBins=600, rRange=[0, 6], rBinWidth=.01,
                         dedxBins=1200, dedxRange=[0, 120], dedxBinWidth=.1),
            'Unexpected original native geometry')
    mainconf = resolve(args.mainconf)
    refs = yaml.safe_load(mainconf.read_text())
    info_path = resolve(Path('config') / refs['analysis'])
    info = yaml.safe_load(info_path.read_text())
    tags = info['starTag']
    require(tuple(tags[k] for k in ('collision', 'year', 'productionTag', 'libraryTag')) ==
            ('auau3p85', 'year2021', 'P24iy', 'SL24y'), 'Explicit 2021 P24iy / SL24y runtime is required')
    output = resolve(args.output_dir)
    allowed = resolve('rootfile/auau3p85_anaFemtoLambda_4He_KFParticle/pdf/nuclear_pid_cutoff_20261015/highrigidity_tof_20261017_001')
    require(not output.exists() and output.parent == allowed, 'Choose a new coarse-input child below the high-R study')
    settings = dict(schemaVersion=1, coarseningFactor=args.factor,
                    inputXBins=axis['rBins'], inputYBins=axis['dedxBins'],
                    inputXMin=axis['rRange'][0], inputXMax=axis['rRange'][1],
                    inputYMin=axis['dedxRange'][0], inputYMax=axis['dedxRange'][1],
                    inputXWidth=axis['rBinWidth'], inputYWidth=axis['dedxBinWidth'],
                    axisTolerance=1e-12, varianceTolerance=1e-12)
    sources = (Path(__file__).resolve(), HERE / 'ExtractCoarseInputs.C', mainconf,
               info_path, PROJECT / 'script/singularity_root_batch.sh')
    for path in sources:
        protected[str(path)] = digest(path)
    output.mkdir()
    (output / 'source').mkdir()
    (output / 'provenance').mkdir()
    for path in sources:
        shutil.copy2(str(path), str(output / 'source' / path.name))
    for name in ('context.json', 'native_audit.json', 'completed.json'):
        shutil.copy2(str(inputs / name), str(output / 'provenance' / ('native_' + name)))
    with (output / 'settings.tsv').open('x') as stream:
        for key, value in settings.items():
            stream.write('%s\t%.17g\n' % (key, value))
    tagged = {k: v for k, v in context['settings'].items() if k.startswith('tag') and 'Mass2' in k}
    save(output / 'context.json', dict(schemaVersion=1, processingStage='exact_r_coarsening',
         inputFingerprint=context['inputFingerprint'], runtimeMainconf=str(mainconf),
         nativeCache=dict(path=str(inputs / 'cache.root'), sha256=protected[str(inputs / 'cache.root')]),
         nativeInputStudy=str(inputs), nativeInputAuditSha256=protected[str(inputs / 'native_audit.json')],
         nativeCompletedSha256=protected[str(inputs / 'completed.json')], cacheKeys=context['cacheKeys'],
         originalInputHistogramKeys=context['inputHistogramKeys'], noFit=True, productionEnabled=False,
         settings=dict(tagged, coarseningFactor=args.factor, inputRigidityBinWidth=.01,
                       outputRigidityBinWidth=.05, sourceNativeRigidityBinWidth=.01,
                       nativeDedxBinWidth=.1, aggregation='exact_nonoverlapping_five_native_R_bins'),
         mathematicalSettings=settings, inputNativeAxis=axis,
         protectedPreRun=protected,
         createdJST=datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).isoformat(),
         limitations=['Rebinned R inputs are 0.05 wide, not native 0.01 measurements.',
                      'This stage only sums original counts and stored variances; no fit or average of purities.',
                      'Source TOF D window is D/He4-compatible, not pure-He4 truth or the wider 3-4 window.',
                      'ALL/TAG are correlated; REMAINING retains its original subset-error policy.',
                      'Original 10GB ROOT is not read or hashed here; only its recorded fingerprint is checked.']))
    call = str(output / 'source/ExtractCoarseInputs.C') + '+(' + ','.join(json.dumps(str(p)) for p in
            (inputs / 'cache.root', output, output / 'settings.tsv')) + ')'
    command = ['bash', str(PROJECT / 'script/singularity_root_batch.sh'), str(mainconf),
               'root4star', '-l', '-b', '-q', call]
    save(output / 'command.json', dict(cwd=str(PROJECT), argv=command))
    with (output / 'root.log').open('x') as log:
        subprocess.run(command, cwd=str(PROJECT), stdout=log, stderr=subprocess.STDOUT,
                       timeout=args.execution_timeout, check=True)
    log_text = (output / 'root.log').read_text()
    require('COARSE_INPUT_SELFTEST_PASS' in log_text and 'COARSE_INPUTS_COMPLETE' in log_text,
            'Grouping/self-test did not complete; inspect preserved root.log')
    require(fingerprint(original) == context['inputFingerprint'] and
            all(digest(Path(path)) == expected for path, expected in protected.items()),
            'Original input/source/provenance changed during extraction')
    result = load(output / 'native_audit.json')
    require(all(result[k] for k in ('exactFiveBinSumPassed', 'countErrorFlowPreservedBySummation',
                                   'writtenCacheReadbackExact', 'storedSumw2ArrayExact')),
            'Coarse cache verification failed')
    save(output / 'completed.json', dict(executionComplete=True, noFit=True, productionEnabled=False,
         protectedPostRunVerified=True, inputFingerprintPreserved=True,
         originalNativeCacheSha256=protected[str(inputs / 'cache.root')],
         products={str(p.relative_to(output)): digest(p) for p in output.rglob('*') if p.is_file()
                   and p.parent != output / 'source'},
         sourceSha256={p.name: digest(p) for p in (output / 'source').iterdir() if p.is_file()},
         mathematicalSelfTestPassed=True, closedReadbackAllCellsPassed=True))
    print('COARSE_INPUT_PREPARE_COMPLETE', output, flush=True)


if __name__ == '__main__':
    main()
