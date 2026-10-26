#!/usr/bin/env python3
"""Audit old fitted gates, not a new calibration or purity certification. Never re-fit."""
import argparse
from collections import Counter, defaultdict
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import shutil
import sys

import yaml

HERE = Path(__file__).resolve().parent
PROJECT = HERE.parents[2]
sys.path.insert(0, str(HERE.parent))
import summarize_cutoff as old_summary


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def digest(path):
    h = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            h.update(block)
    return h.hexdigest()


def read(path):
    return json.loads(Path(path).read_text())


def save(path, value):
    with Path(path).open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def resolve(value):
    path = Path(value)
    path = (path if path.is_absolute() else PROJECT / path).resolve(strict=True)
    if PROJECT not in path.parents:
        raise ValueError('Input must belong to this project: ' + str(path))
    return path


def percentile(values, fraction):
    values = sorted(v for v in values if finite(v))
    if not values:
        return None
    where = fraction * (len(values) - 1)
    lo = int(where)
    hi = min(lo + 1, len(values) - 1)
    return values[lo] + (where - lo) * (values[hi] - values[lo])


def chi(row):
    x, n = row.get('pearson'), row.get('conditionalNdf')
    return x / n if finite(x) and finite(n) and n > 0 and x >= 0 else None


def choose(models, threshold, cfg):
    # Crucial: fit-quality filtering is independent of the fitted purity.
    kept = [m for m in models if m['supportedWithoutChi2'] and
            finite(m['maxChi2Ndf']) and m['maxChi2Ndf'] <= threshold]
    nominal = [m for m in kept if m['variation'] == 'nominal']
    coverage = (len(kept) >= cfg['minimumModels'] and
                len(nominal) >= cfg['minimumNominalModels'] and
                len({m['response'] for m in nominal}) >= cfg['minimumResponseFamilies'] and
                len({m['background'] for m in nominal}) >= cfg['minimumBackgroundFamilies'])
    purities = [m['purity'] for m in kept if finite(m['purity'])]
    low = min(purities) if purities else None
    supported = coverage and len(purities) == len(kept)
    return dict(retainedModels=len(kept), retainedNominalModels=len(nominal),
                responseFamilies=sorted({m['response'] for m in nominal}),
                backgroundFamilies=sorted({m['background'] for m in nominal}),
                purityLow=low, purityHigh=max(purities) if purities else None,
                qualityCoverage=supported,
                status='accepted_reference_only' if supported and low >= cfg['purityThreshold'] else
                       'below_purity_threshold' if supported else 'unassessed',
                retainedLabels=sorted(m['label'] for m in kept))


def write_tsv(path, rows, fields):
    with path.open('x') as stream:
        writer = csv.DictWriter(stream, fieldnames=fields, delimiter='\t', extrasaction='ignore')
        writer.writeheader()
        for row in rows:
            writer.writerow({k: json.dumps(v, ensure_ascii=False) if isinstance(v, (list, dict))
                             else '' if v is None else v for k, v in row.items() if k in fields})


def run(args):
    settings = args.settings.resolve(strict=True)
    cfg = yaml.safe_load(settings.read_text())
    if cfg.get('schemaVersion') != 1 or cfg['purityThreshold'] != .5:
        raise ValueError('Expected explicit purity S/(S+B) >= 0.5.')
    output = args.output_dir.resolve()
    allowed = (PROJECT / cfg['outputBase']).resolve()
    if output.exists() or output.parent != allowed:
        raise ValueError('Require a new direct child of outputBase, never overwrite.')
    protected = {}

    def checked(folder, name):
        path = folder / name
        done = read(folder / 'completed.json')
        if not done.get('executionComplete') or digest(path) != done['products'][name]:
            raise ValueError('Incomplete/changed input: ' + str(path))
        protected[str(folder / 'completed.json')] = digest(folder / 'completed.json')
        protected[str(path)] = digest(path)
        return read(path)

    low_folder = resolve(cfg['lowSummary'])
    high_folder = resolve(cfg['highFinalSummary'])
    low = checked(low_folder, 'summary.json')
    high = checked(high_folder, 'summary.json')
    raw_path = resolve(cfg['rawGateCounts'])
    raw_done_path = raw_path.parent / 'completed.json'
    raw_done = read(raw_done_path)
    if not raw_done.get('executionComplete'):
        raise ValueError('Independent raw-gate extraction is incomplete.')
    for path in (raw_path, raw_path.parent/'gates.tsv'):
        expected = raw_done['products'][path.name]
        if digest(path) != expected:
            raise ValueError('Independent raw-gate product changed: ' + str(path))
        protected[str(path)] = expected
    protected[str(raw_done_path)] = digest(raw_done_path)
    raw_gates_sha = digest(raw_path.parent/'gates.tsv')
    raw_counts = {}
    with raw_path.open() as stream:
        for row in csv.DictReader(stream, delimiter='\t'):
            key = int(row['rBin']), row['species']
            if key in raw_counts:
                raise ValueError('Duplicate independent gate-count row.')
            value = {name:float(row[name]) for name in ('rawGateCount',
                     'rawGateCountLowerFullBins','rawGateCountUpperTouchingBins')}
            if not all(finite(x) and x >= 0 for x in value.values()) or not (
                    value['rawGateCountLowerFullBins'] <= value['rawGateCount'] <=
                    value['rawGateCountUpperTouchingBins']):
                raise ValueError('Invalid independent raw gate count or bounds.')
            raw_counts[key] = value
    if len(raw_counts) != 640:
        raise ValueError('Expected complete 4-species x 160-native-R raw-gate coverage.')
    low_models = {(r['species'], round(r['rLo'], 6), r['label']): r for r in low['detailedRows']}
    native_rows, grouped = [], []
    runs = [('low_ALL', Path(p)) for p in low['inputStudies'].values()]
    high_native = resolve(high['nativeSummary'])
    high_context = checked(high_native, 'context.json')
    high_runs = sorted({Path(p).parent for p in high_context['protectedInputs'] if p.endswith('/cutoff_fits.json')})
    runs += [('high_TOF_native', p) for p in high_runs]
    runs += [('high_TOF_coarse', Path(p)) for p in high['coarseRuns']]
    for stage, folder in runs:
        context = checked(folder, 'context.json')
        fit = checked(folder, 'cutoff_fits.json')
        cfg_fit = checked(folder, 'effective_settings.json')
        audits = {r['index']: r for r in checked(folder, 'tag_shape_audit.json')}
        target = context['targetSpecies']
        if stage == 'low_ALL':
            run_gate = folder/'gates.tsv'
            done = read(folder/'completed.json')
            if digest(run_gate) != done['products']['gates.tsv'] or digest(run_gate) != raw_gates_sha:
                raise ValueError('Low-R fit and independent counter used different TPC gates.')
            protected[str(run_gate)] = raw_gates_sha
        response, background, variation = (context[k] for k in ('response', 'background', 'variation'))
        label = '_'.join((response, background, variation))
        if stage != 'low_ALL':
            saved_native = checked(folder, 'native_summary.json')
            saved_native = {r['rBin']: r for r in saved_native['nativeRows']}
        groups = defaultdict(list)
        for row in fit['slices']:
            value = old_summary.aggregate_one([row], target, audits,
                                             dict(cfg_fit, fitWindowMode='target_local'),
                                             raw_counts=raw_counts if stage=='low_ALL' else None)
            if stage != 'low_ALL':
                saved = saved_native[row['rBin']]
                value['supported'] = saved['quality']
                value['reasons'] = saved['reasons']
            x = chi(row)
            rlo, rhi = row['lo'], row['hi']
            inside = any(rhi > a - cfg['focusPadding'] and rlo < b + cfg['focusPadding']
                         for a, b in cfg['focusIntervals'][target])
            native_rows.append(dict(stage=stage, target=target, response=response, background=background,
                variation=variation, label=label, run=str(folder), rowIndex=row['index'], rBin=row['rBin'],
                rLo=rlo, rHi=rhi, chi2=row.get('pearson'), ndf=row.get('conditionalNdf'), chi2Ndf=x,
                purity=value['purity'], localResidualFrac=value['localResidualFrac'],
                integralResidualFrac=value['integralResidualFrac'], numericalValid=row['numericalValid'],
                supportedWithoutChi2=value['supported'], reasons=value['reasons'], inFocus=inside))
            edge = round(math.floor((rlo + 1.e-7) / cfg['aggregateWidth']) * cfg['aggregateWidth'], 6)
            groups[edge].append(row)
        for edge, rows in sorted(groups.items()):
            value = old_summary.aggregate_one(rows, target, audits,
                                             dict(cfg_fit, fitWindowMode='target_local'),
                                             raw_counts=raw_counts if stage=='low_ALL' else None)
            if stage == 'low_ALL':
                original = low_models[target, edge, label]
                for key in ('signalSelected', 'totalSelected', 'rawSelected', 'purity'):
                    a, b = value[key], original[key]
                    if a is None or b is None:
                        if a != b:
                            raise ValueError('Old aggregation availability differs.')
                    elif abs(a-b) > 1.e-8 * max(1., abs(a), abs(b)):
                        raise ValueError('Old aggregation was not reproduced: ' + key)
                if value['reasons'] != original['reasons']:
                    raise ValueError('Old aggregate support reasons differ.')
                supported, reasons = original['supported'], original['reasons']
            else:
                bad = [saved_native[r['rBin']] for r in rows if not saved_native[r['rBin']]['quality']]
                supported = value['supported'] and not bad
                reasons = sorted(set(value['reasons'] + [s for r in bad for s in r['reasons']]))
            xs = [chi(r) for r in rows]
            grouped.append(dict(stage=stage, target=target, response=response, background=background,
                variation=variation, label=label, run=str(folder), rLo=edge, rHi=edge+cfg['aggregateWidth'],
                purity=value['purity'], supportedWithoutChi2=supported, reasons=reasons,
                maxChi2Ndf=max(xs) if all(finite(x) for x in xs) else None,
                summedChi2Ndf=sum(r['pearson'] for r in rows)/sum(r['conditionalNdf'] for r in rows)
                    if all(finite(x) for x in xs) and sum(r['conditionalNdf'] for r in rows)>0 else None,
                localResidualFrac=value['localResidualFrac'], integralResidualFrac=value['integralResidualFrac'],
                nativeRows=len(rows)))
    distributions = []
    for stage, target, variation in sorted({(r['stage'], r['target'], 'nominal' if r['variation']=='nominal' else 'sensitivity') for r in native_rows}):
        sample = [r for r in native_rows if r['stage']==stage and r['target']==target and
                  ('nominal' if r['variation']=='nominal' else 'sensitivity')==variation]
        xs = [r['chi2Ndf'] for r in sample if finite(r['chi2Ndf'])]
        distributions.append(dict(stage=stage, target=target, variation=variation, rows=len(sample), finiteRows=len(xs),
            p0=percentile(xs,0), p25=percentile(xs,.25), p50=percentile(xs,.5),
            p75=percentile(xs,.75), p90=percentile(xs,.9), p100=percentile(xs,1),
            passCounts={str(t):sum(x<=t for x in xs) for t in cfg['chi2NdfThresholds']}))
    scan = []
    group_map = defaultdict(list)
    for item in grouped:
        group_map[item['stage'],item['target'],round(item['rLo'],6)].append(item)
    old_bins = {(r['species'],round(r['rLo'],6)):r for r in low['aggregateRows']}
    for (stage,target,edge), models in sorted(group_map.items()):
        old_bin = old_bins.get((target,edge)) if stage=='low_ALL' else None
        for threshold in cfg['chi2NdfThresholds']:
            result = choose(models, threshold, cfg)
            scan.append(dict(stage=stage,target=target,rLo=edge,rHi=edge+cfg['aggregateWidth'],
                chi2Threshold=threshold, originalQuality=old_bin['quality'] if old_bin else None,
                originalStatus=old_bin['status'] if old_bin else None, **result))
    sources = [settings, Path(__file__).resolve(), Path(old_summary.__file__).resolve(),
               HERE.parent/'run_cutoff.py']
    protected.update({str(p):digest(p) for p in sources})
    output.mkdir(parents=True)
    (output/'source').mkdir()
    for path in sources:
        shutil.copy2(path, output/'source'/path.name)
    write_tsv(output/'native_quality.tsv', native_rows, list(native_rows[0]))
    write_tsv(output/'focus_native_quality.tsv', [r for r in native_rows if r['inFocus']], list(native_rows[0]))
    write_tsv(output/'grouped_quality.tsv', grouped, list(grouped[0]))
    write_tsv(output/'threshold_scan.tsv', scan, list(scan[0]))
    save(output/'distributions.json', distributions)
    save(output/'context.json', dict(createdUtc=datetime.datetime.now(datetime.timezone.utc).isoformat(),
        protectedInputs=protected, noFit=True, coefficientsChanged=False, productionEnabled=False,
        semantics='Old fitted-gate reference only, NOT a certification of restored polynomial gates.',
        thresholdPolicy='Exclude quality failures before inspecting purity; worst native chi2/ndf screen.',
        coveragePolicy='Nominal response/background family coverage required, variations cannot replace nominal families.'))
    reason_counts = Counter(reason for row in native_rows if row['inFocus'] for reason in row['reasons'])
    text = ['# 保存fitの品質監査（2026-10-19）', '',
        '旧補間gateのfitを読取専用で再集計した診断。元多項式gateのpurity再認証ではない。',
        '元ROOT・fit・PID設定・旧PDFは変更していない。新規fitなし。', '',
        f'対象: {len(runs)}実行、{len(native_rows)} fit行。低R80実行、高Rnative16実行、He4 coarse4実行。',
        'χ²/ndfは有限R投影のtarget-local窓のPearson値。品質screenをPの判定より先に適用。',
        '0.05区間は全native fitの最大χ²/ndfを使う。支持済み条件のPはΣS/Σ(S+B)。',
        '高RのHe4 nativeとcoarseは別stageで比較し、良い点の継ぎ接ぎはしない。', '',
        '## 低R nominalのχ²/ndf分布', '', '| Species | N | Median | 90% quantile | Max |', '|---|---:|---:|---:|---:|']
    for row in distributions:
        if row['stage']=='low_ALL' and row['variation']=='nominal':
            text.append('| {} | {} | {:.3f} | {:.3f} | {:.3f} |'.format(row['target'],row['finiteRows'],row['p50'],row['p90'],row['p100']))
    text += ['', '## 主な注意', '',
        '不良modelを除外しても、応答2family・背景2family・nominal2条件以上が残らなければ未評価。',
        '2条件による両family被覆は2×2全組の保持とは異なる。各行のretainedLabelsを必ず確認する。',
        'χ²cutの採用値をこの監査だけで決定しない。隣接binと除外fitの局所図を確認する。',
        '飛び地はfit未収束・帰属・応答精度・残差等で発生し、真の低purityと区別する。', '',
        '## 出力', '',
        '- native_quality.tsv: 失敗を省かない全fitのχ²・残差・旧支持理由。',
        '  native支持は行単位の診断。旧0.05区間の集約支持はgrouped_quality.tsvを参照。',
        '- focus_native_quality.tsv: 下限、飛び地とその両隣。',
        '- grouped_quality.tsv: 0.05区間ごとの元支持とχ²。低Rの旧集約数値・理由を照合済み。',
        '- threshold_scan.tsv: 1.5/2/3/5/10/20の品質限定envelope（旧gate参考）。',
        '- distributions.json: 低R、高Rnative、高Rcoarseを分けた分位値。', '',
        '注目域の理由集計（同区間・感度条件の重複を含み、独立事象数ではない）:', '',
        '```json', json.dumps(reason_counts,ensure_ascii=False,indent=2), '```', '']
    with (output/'README.md').open('x') as stream:
        stream.write('\n'.join(text))
    if any(digest(Path(p))!=expected for p,expected in protected.items()):
        raise ValueError('Protected input/source changed during audit.')
    save(output/'completed.json', dict(executionComplete=True, noFit=True, productionEnabled=False,
        originalSourcesPreserved=True, oldLowAggregationReproduced=True, fitRows=len(native_rows), runs=len(runs),
        products={str(p.relative_to(output)):digest(p) for p in output.rglob('*') if p.is_file()}))
    print('AUDIT_FIT_QUALITY_COMPLETE', output)
    for row in distributions:
        if row['variation']=='nominal':
            print('DIST', json.dumps(row, ensure_ascii=False))
    for threshold in cfg['chi2NdfThresholds']:
        print('SCAN',threshold,Counter((r['stage'],r['target'],r['status']) for r in scan if r['chi2Threshold']==threshold))


if __name__ == '__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings', required=True, type=Path)
    parser.add_argument('--output-dir', required=True, type=Path)
    run(parser.parse_args())
