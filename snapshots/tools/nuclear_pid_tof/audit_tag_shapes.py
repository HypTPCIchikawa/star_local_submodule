#!/usr/bin/env python3
"""Read-only residual diagnostics of completed native TAG fits. Never relax cuts."""
import argparse
from datetime import datetime, timedelta, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[1]
SPECIES = ('d', 't', 'he3', 'he4')


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def fingerprint(path):
    stat = Path(path).stat()
    return dict(path=str(path), size=stat.st_size, mtime_ns=stat.st_mtime_ns,
                device=stat.st_dev, inode=stat.st_ino)


def save(path, value):
    path.write_text(json.dumps(value, ensure_ascii=False, indent=2, allow_nan=False)+'\n')


def read_json(path):
    def invalid(value):
        raise ValueError('Invalid JSON constant '+value)
    return json.loads(Path(path).read_text(), parse_constant=invalid)


def source_path(path):
    path = Path(path).resolve(strict=True)
    if PROJECT not in path.parents:
        raise ValueError('Expected project-owned completed data/config/source.')
    return path


def finite_positive(value):
    return not isinstance(value, bool) and isinstance(value, (int, float)) \
        and math.isfinite(value) and value > 0


def ratio(a, b):
    return a/b if b > 0 else None


def aggregate_focus(records, ranges):
    summary = []
    for lo, hi in ranges:
        groups = {}
        for row in records:
            if row['lo'] < lo-1.e-9 or row['hi'] > hi+1.e-9:
                continue
            key = (row['sample'], row['species'], row['nSigma'], row['aggregateKeV'])
            groups.setdefault(key, []).append(row)
        for key, selected in sorted(groups.items()):
            observed = sum(r['observed'] for r in selected)
            predicted = sum(r['predicted'] for r in selected)
            l1 = sum(r['l1'] for r in selected)
            active = sum(r['activeGroups'] for r in selected)
            pearson = sum(r['pearson'] for r in selected) if all(r['pearson'] is not None for r in selected) else None
            deviance = sum(r['deviance'] for r in selected) if all(r['deviance'] is not None for r in selected) else None
            primary = sum(r['primaryModelCount'] for r in selected)
            summary.append(dict(rigidityRange=[lo, hi], sample=key[0], species=key[1],
                 nSigma=key[2], aggregateKeV=key[3], nativeSlices=len(selected), observed=observed,
                 predicted=predicted, observedOverPredicted=ratio(observed, predicted),
                 primaryModelCount=primary, modeledPrimaryFraction=ratio(primary, predicted),
                 l1OverObserved=ratio(l1, observed), pearson=pearson, deviance=deviance,
                 activeGroups=active, pearsonPerActiveGroup=ratio(pearson, active) if pearson is not None else None,
                 sqrtPearsonOverObserved=math.sqrt(pearson/observed) if pearson is not None and observed > 0 else None,
                 minimumPhysicalGaussianCoverage=min(r['physicalGaussianCoverage'] for r in selected),
                 numericalValidSlices=sum(r['numericalValid'] for r in selected),
                 meanMinimum=min(r['mean'] for r in selected), meanMaximum=max(r['mean'] for r in selected),
                 projectionSigmaMinimum=min(r['projectionSigma'] for r in selected),
                 projectionSigmaMaximum=max(r['projectionSigma'] for r in selected)))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study', type=Path, action='append', required=True)
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--roi-nsigma', type=float, action='append', required=True)
    parser.add_argument('--aggregate-kev', type=float, action='append', required=True)
    parser.add_argument('--focus-range', action='append', required=True, help='Explicit lo:hi rigidity domain for summary only.')
    parser.add_argument('--timeout-seconds', type=int, required=True)
    args = parser.parse_args()
    if args.timeout_seconds <= 0 or not all(finite_positive(x) for x in args.roi_nsigma+args.aggregate_kev):
        raise ValueError('Invalid explicitly declared diagnostic parameters.')
    if len(set(args.roi_nsigma)) != len(args.roi_nsigma) or len(set(args.aggregate_kev)) != len(args.aggregate_kev):
        raise ValueError('Duplicate ROI or grouping widths.')
    ranges = []
    for item in args.focus_range:
        lo, hi = map(float, item.split(':'))
        if not finite_positive(lo) or not finite_positive(hi) or hi <= lo:
            raise ValueError('Bad focus rigidity interval.')
        ranges.append([lo, hi])
    mainconf = source_path(args.mainconf)
    main = yaml.safe_load(mainconf.read_text())
    analysis_path = source_path(PROJECT/'config'/main['analysis'])
    star = yaml.safe_load(analysis_path.read_text())['starTag']
    if (star['collision'], star['year'], star['productionTag'], star['libraryTag']) != \
            ('auau3p85', 'year2021', 'P24iy', 'SL24y'):
        raise ValueError('Only 3p85 2021 P24iy SL24y diagnostic runtime is allowed.')
    studies = sorted(set(source_path(p) for p in args.study))
    if len(studies) != len(args.study) or len({p.parent for p in studies}) != 1:
        raise ValueError('Choose distinct sibling native model studies.')
    output = args.output_dir.resolve()
    if output.exists() or output.parent != studies[0].parent:
        raise ValueError('Diagnostic output must be a new sibling directory.')
    protected = {str(p):digest(p) for p in (mainconf, analysis_path, Path(__file__).resolve(),
                 HERE/'AuditTagShapes.C', HERE/'NuclearFitCommon.h', PROJECT/'script/singularity_root_batch.sh')}
    source_info = []
    original_fingerprints = {}
    for study in studies:
        context = read_json(study/'context.json')
        done = read_json(study/'completed.json')
        if not done.get('executionComplete') or done.get('productionEnabled') is not False \
                or not done.get('inputFingerprintPreserved') or not done.get('protectedSourcesPreserved') \
                or context['selectedPartition'] != 'native' or context['runtimeMainconf'] != str(mainconf):
            raise ValueError('Only complete protected native TAG studies can be diagnosed.')
        for name, expected in done['products'].items():
            path = study/name
            if digest(path) != expected:
                raise ValueError('Completed TAG artifact changed: '+str(path))
            protected[str(path)] = expected
        protected[str(study/'completed.json')] = digest(study/'completed.json')
        for key in ('protected', 'snapshots'):
            for path, expected in context[key].items():
                if digest(Path(path)) != expected:
                    raise ValueError('TAG source/input changed: '+path)
                protected[path] = expected
        original = Path(context['inputFingerprint']['path']).resolve(strict=True)
        if '2021P24iy_partial_merge' not in original.name or fingerprint(original) != context['inputFingerprint']:
            raise ValueError('Original 2021 P24iy ROOT fingerprint changed.')
        original_fingerprints[str(original)] = fingerprint(original)
        fit = read_json(study/'tag_fits.json')
        if len(fit['slices']) != done['nSlices']:
            raise ValueError('Native TAG slice coverage incomplete.')
        rows = []
        for row in fit['slices']:
            if row['firstNativeRBin'] != row['lastNativeRBin'] \
                    or abs(row['hi']-row['lo']-context['sourceNativeWidth']) > context['settings']['numeric']['axisTolerance']:
                raise ValueError('Diagnostic source is not exact single native-R slices.')
            common = [row['index'], row['sample'], row['lo'], row['hi']]
            tail = [row['fitLo'], row['fitHi'], int(row['numericalConverged'])]
            rows.append(common+['full_fit', 0., 0.]+tail)
            for species, component in row['components'].items():
                if component['active'] and (component['freeMean'] or component['freeSigma']):
                    if species not in SPECIES or not finite_positive(component['sigma']) or not math.isfinite(component['mean']):
                        raise ValueError('Invalid retained primary shape.')
                    rows.append(common+[species, component['mean'], component['sigma']]+tail)
        source_info.append((study, context, done, rows))
    if len({context['inputStudy'] for study, context, done, rows in source_info}) != 1:
        raise ValueError('Background models do not share immutable input cache.')
    output.mkdir()
    snapshot = output/'source'
    snapshot.mkdir()
    snapshots = {}
    for src in (Path(__file__).resolve(), HERE/'AuditTagShapes.C', HERE/'NuclearFitCommon.h',
                PROJECT/'script/singularity_root_batch.sh', mainconf, analysis_path):
        dest = snapshot/src.name
        shutil.copy2(src, dest)
        if digest(src) != digest(dest):
            raise ValueError('Diagnostic source changed while copying snapshot.')
        snapshots[str(dest)] = digest(dest)
    save(output/'diagnostic_settings.json', dict(roiNSigma=args.roi_nsigma,
         aggregateWidthsKeV=args.aggregate_kev, focusRanges=ranges, timeoutSeconds=args.timeout_seconds,
         noFit=True, noGateChanges=True, projectionWidthOnly=True,
         comparisonPolicy='Per-model own retained mean/sigma ROI; ROI edges can differ between backgrounds.',
         aggregationPolicy='Counts summed only; not new fit ndf or independently fitted coarse widths.'))
    snapshots[str(output/'diagnostic_settings.json')] = digest(output/'diagnostic_settings.json')
    contexts = []
    for study, context, done, rows in source_info:
        folder = output/context['background']
        if folder.exists():
            raise ValueError('Duplicate background diagnostic destination.')
        folder.mkdir()
        settings = dict(roiNSigma=','.join(format(x, '.17g') for x in args.roi_nsigma),
             aggregateWidthsKeV=','.join(format(x, '.17g') for x in args.aggregate_kev),
             axisTolerance=context['settings']['numeric']['axisTolerance'],
             nativeDedxBinWidth=context['settings']['numeric']['nativeDedxBinWidth'])
        (folder/'settings.tsv').write_text(''.join(str(k)+'\t'+str(v)+'\n' for k, v in settings.items()))
        header='index sample rLo rHi species mean sigma fitLo fitHi numericalValid'.split()
        with (folder/'rows.tsv').open('w') as stream:
            stream.write('\t'.join(header)+'\n')
            for row in rows:
                stream.write('\t'.join(format(x, '.17g') if isinstance(x, float) else str(x) for x in row)+'\n')
        for path in (folder/'settings.tsv', folder/'rows.tsv'):
            snapshots[str(path)] = digest(path)
        save(folder/'context.json', dict(sourceStudy=str(study), sourceCompletionSha256=digest(study/'completed.json'),
             sourceFitRootSha256=digest(study/'tag_fits.root'), noFit=True, productionEnabled=False))
        contexts.append((study, context, folder))
    save(output/'context.json', dict(schemaVersion=1, noFit=True, productionEnabled=False,
         sourceStudies=list(map(str, studies)), originalInputFingerprints=original_fingerprints,
         protected=protected, snapshots=snapshots,
         createdJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))

    def preserve():
        for path, expected in {**protected, **snapshots}.items():
            if digest(Path(path)) != expected:
                raise ValueError('Protected original/snapshot changed: '+path)
        for path, expected in original_fingerprints.items():
            if fingerprint(Path(path)) != expected:
                raise ValueError('Original ROOT fingerprint changed.')

    preserve()
    comparison = []
    for study, context, folder in contexts:
        macro = snapshot/'AuditTagShapes.C'
        call = str(macro)+'++('+','.join(json.dumps(str(path)) for path in
                (study/'tag_fits.root', folder, folder/'settings.tsv', folder/'rows.tsv'))+')'
        command = ['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
                   'root4star', '-l', '-b', '-q', call]
        save(folder/'command.json', command)
        print('NUCLEAR_TAG_AUDIT_RUN', context['background'], folder, flush=True)
        try:
            with (folder/'root.log').open('w') as log:
                run = subprocess.run(command, cwd=PROJECT, stdout=log, stderr=subprocess.STDOUT,
                                     timeout=args.timeout_seconds)
        finally:
            preserve()
        if run.returncode or 'NUCLEAR_TAG_RESIDUAL_AUDIT_COMPLETE' not in (folder/'root.log').read_text():
            raise RuntimeError('Retained TAG diagnostic failed: '+str(folder/'root.log'))
        data = read_json(folder/'residual_audit.json')
        if not data['noFit'] or data['productionEnabled'] or data['recordCount'] != len(data['records']):
            raise ValueError('Invalid residual diagnostic output.')
        summary = aggregate_focus(data['records'], ranges)
        save(folder/'focus_summary.json', summary)
        comparison.append(dict(background=context['background'], sourceStudy=str(study),
                                recordCount=data['recordCount'], focusSummary=summary))
    save(output/'comparison.json', dict(schemaVersion=1, noFit=True, projectionWidthOnly=True,
         productionEnabled=False, identityValidated=False,
         warning='Native and aggregate Pearson statistics use retained predictions, not fitted ndf. Own-model ROIs may differ; no gate thresholds changed.',
         models=comparison))
    preserve()
    save(output/'completed.json', dict(executionComplete=True, noFit=True, noGateChanges=True,
         sourceInputPreserved=True, productionEnabled=False, physicsValidationComplete=False,
         products={str(p.relative_to(output)):digest(p) for p in output.rglob('*')
                   if p.is_file() and snapshot not in p.parents and p.name!='completed.json'},
         completedJst=datetime.now(timezone(timedelta(hours=9))).isoformat()))
    print('NUCLEAR_TAG_AUDIT_COMPLETE', output, flush=True)


if __name__ == '__main__':
    main()
