#!/usr/bin/env python3
"""Prepare empirical-mode-centred ROI diagnostics, without changing Git PID seeds.

Only a NEW eight-slice smoke output is prepared. Fits are launched separately
through prepare_tag_study.py, which checks the completed plan's frozen hashes.
Three-bin smoothing is used solely to locate an ROI, never as fit data.
"""
import argparse
import copy
import csv
import json
import math
from pathlib import Path
import shutil

import yaml

import prepare_tag_study as base


def native_spectra(path, planned, axis):
    """Validate the exact representative export, including flow and zero bins."""
    expected = {row['id']: row for row in planned}
    result = {identity: {} for identity in expected}
    with path.open() as stream:
        reader = csv.DictReader(stream, delimiter='\t')
        base.require(reader.fieldnames == ['id', 'mode', 'target', 'rBin', 'dedxBin',
                                          'lo', 'hi', 'count', 'error'],
                     'Unexpected native export column contract')
        for row in reader:
            identity = row['id']
            base.require(identity in expected, 'Raw export has an unrequested slice')
            planned_row = expected[identity]
            base.require(row['mode'] == planned_row['mode'] and
                         row['target'] == planned_row['target'] and
                         int(row['rBin']) == planned_row['rBin'],
                         'Raw export slice identity differs from fresh plan')
            index = int(row['dedxBin'])
            base.require(str(index) == row['dedxBin'] and
                         0 <= index <= axis['dedxBins'] + 1 and
                         index not in result[identity],
                         'Invalid/duplicate native bin index')
            values = {key: float(row[key]) for key in ('lo', 'hi', 'count', 'error')}
            base.require(all(math.isfinite(v) for v in values.values()) and
                         values['count'] >= 0 and values['error'] >= 0 and
                         values['count'].is_integer(),
                         'Native export requires finite unweighted integer counts')
            expected_lo = axis['dedxRange'][0] + (index - 1) * axis['dedxBinWidth']
            base.require(abs(values['lo'] - expected_lo) < 1.e-8 and
                         abs(values['hi'] - expected_lo - axis['dedxBinWidth']) < 1.e-8,
                         'Native export axes differ from raw cache')
            base.require(abs(values['error'] ** 2 - values['count']) <=
                         1.e-8 * max(1., values['count']),
                         'Native error/count contract changed')
            values['index'] = index
            values['center'] = (values['lo'] + values['hi']) / 2.
            result[identity][index] = values
    required = set(range(axis['dedxBins'] + 2))
    base.require(all(set(spectrum) == required for spectrum in result.values()),
                 'Full native spectrum including zero/flow bins is required')
    return result


def validate_mode_settings(settings):
    required = {'searchHalfWidthSigma', 'smoothingWindowBins', 'modeTiePolicy',
                'searchCenter', 'roiWidth', 'centerDefinition'}
    base.require(set(settings) == required, 'Mode-finder settings must be explicit')
    base.require(settings['searchCenter'] == 'original_git_mean' and
                 settings['roiWidth'] == 'original_git_sigma' and
                 settings['centerDefinition'] == 'native_bin_center' and
                 settings['modeTiePolicy'] == 'nearest_original_mean_then_lower_bin',
                 'Unsupported mode/search/width policy')
    base.require(type(settings['smoothingWindowBins']) is int and
                 settings['smoothingWindowBins'] == 3 and
                 type(settings['searchHalfWidthSigma']) in (float, int) and
                 math.isfinite(settings['searchHalfWidthSigma']) and
                 settings['searchHalfWidthSigma'] > 0,
                 'This diagnostic requires a positive search width and a three-bin finder')


def locate_mode(spectrum, mean, sigma, mode_settings, axis):
    validate_mode_settings(mode_settings)
    base.require(math.isfinite(mean) and math.isfinite(sigma) and sigma > 0,
                 'Require finite original Git seed geometry')
    raw_lo = mean - mode_settings['searchHalfWidthSigma'] * sigma
    raw_hi = mean + mode_settings['searchHalfWidthSigma'] * sigma
    # Search centers, not flow bins. The immediate raw neighbors are allowed
    # outside the search interval, but must remain regular native bins.
    candidates = [index for index, row in spectrum.items()
                  if 2 <= index <= axis['dedxBins'] - 1 and
                  raw_lo <= row['center'] <= raw_hi]
    base.require(candidates, 'No complete three-bin native search neighborhood')
    scores = {index: sum(spectrum[neighbor]['count']
                        for neighbor in (index - 1, index, index + 1))
              for index in candidates}
    maximum = max(scores.values())
    base.require(maximum > 0, 'Empty search interval: no fallback or zero-purity inference')
    tied = [index for index in candidates if scores[index] == maximum]
    chosen = min(tied, key=lambda index: (abs(spectrum[index]['center'] - mean), index))
    return dict(searchLo=raw_lo, searchHi=raw_hi,
                searchAxisClipped=raw_lo < axis['dedxRange'][0] or
                                  raw_hi > axis['dedxRange'][1],
                searchFirstCandidate=min(candidates), searchLastCandidate=max(candidates),
                modeBin=chosen, modeCenter=spectrum[chosen]['center'],
                modeRawCount=spectrum[chosen]['count'], modeThreeBinSum=maximum,
                tiedMaximumBins=sorted(tied),
                searchBoundaryMode=chosen in (min(candidates), max(candidates)),
                neighborBins=[chosen - 1, chosen, chosen + 1],
                smoothingAppliedToFitCounts=False)


def recentered_roi(center, sigma, half_width, alignment, axis):
    base.require(all(math.isfinite(v) for v in (center, sigma, half_width, alignment)) and
                 sigma > 0 and half_width > 0 and alignment > 0,
                 'Finite positive ROI width/alignment required')
    base.require(abs(alignment / axis['dedxBinWidth'] -
                     round(alignment / axis['dedxBinWidth'])) < 1.e-9,
                 'ROI alignment must preserve whole native bins')
    raw_lo, raw_hi = center - half_width * sigma, center + half_width * sigma
    lo = max(axis['dedxRange'][0], math.floor(raw_lo / alignment) * alignment)
    hi = min(axis['dedxRange'][1], math.ceil(raw_hi / alignment) * alignment)
    base.require(hi > lo, 'No native axis support for recentered ROI')
    return dict(fitLo=lo, fitHi=hi, fitWindowAxisClipped=lo > raw_lo or hi < raw_hi,
                unalignedFitLo=raw_lo, unalignedFitHi=raw_hi)


def validate_export(export_dir, cfg, raw_table, completion, provenance):
    base.require(completion.get('status') == 'completed' and
                 completion.get('nativeReadbackExact') is True and
                 completion.get('inputHashesUnchanged') is True and
                 completion.get('noFit') is True and
                 completion.get('productionEnabled') is False and
                 completion.get('tableSha256') == base.sha(raw_table),
                 'Native export is incomplete, changed or not diagnostic-only')
    source_fit = base.owned(provenance['sourceFit'])
    source_plan = base.read(source_fit / 'plan.json')
    expected_cache = base.owned(cfg['inputDirectory']) / 'cache.root'
    base.require(Path(source_plan['cache']).resolve() == expected_cache and
                 source_plan['inputs'][str(expected_cache)] == base.sha(expected_cache),
                 'Native export must originate from the same immutable raw cache')
    # Raw cache provenance matters here; old live configuration hashes do not
    # become new PID seeds. The original executed source snapshot remains historical.
    for name in ('cache.root', 'context.json', 'native_audit.json', 'completed.json'):
        path = base.owned(cfg['inputDirectory']) / name
        base.require(provenance['protected'].get(str(path)) == base.sha(path),
                     'Export raw-input provenance changed: ' + name)
    base.require(completion['summary']['sliceCount'] == 8 and
                 completion['summary']['nonintegerCountBins'] == 0 and
                 completion['summary']['countErrorConsistent'] is True,
                 'Only the verified eight-slice raw export is in scope')


def freeze_extra(plan, output, source):
    copy_path = output / 'source' / source.relative_to(base.PROJECT)
    copy_path.parent.mkdir(parents=True, exist_ok=True)
    base.require(not copy_path.exists(), 'Fresh recenter source snapshot already exists')
    shutil.copy2(str(source), str(copy_path))
    base.require(base.sha(copy_path) == base.sha(source), 'Recenter source changed during freeze')
    plan['inputs'][str(copy_path)] = base.sha(copy_path)
    plan['sourceSha256'][str(source)] = base.sha(source)


def prepare(settings, output, variant):
    stage_inputs = {str(settings): base.sha(settings),
                    str(Path(__file__).resolve()): base.sha(Path(__file__).resolve())}
    cfg = yaml.safe_load(settings.read_text())
    base.require(cfg['schemaVersion'] == 1 and variant in cfg['variants'],
                 'Select an explicit recenter configuration variant')
    validate_mode_settings(cfg['modeFinder'])
    selection = cfg['variants'][variant]
    base.require(set(selection) == {'roiHalfWidthSigma'} and
                 type(selection['roiHalfWidthSigma']) in (float, int) and
                 math.isfinite(selection['roiHalfWidthSigma']) and
                 selection['roiHalfWidthSigma'] > 0,
                 'Explicit finite positive ROI half-width required')
    base_settings = base.owned(cfg['baseSettings'])
    stage_inputs[str(base_settings)] = base.sha(base_settings)
    base_cfg = yaml.safe_load(base_settings.read_text())
    base.require(Path(base_cfg['fitSource']).name in ('FitStorageSlices.C', 'FitSplitNormalSlices.C') and
                 base_cfg['numeric']['storageMapping'] == 'root_default12',
                 'This comparison must retain the explicit ROOT-default12 response')
    base.require(base_cfg['smokeNativeRBins'] == cfg['smokeNativeRBins'] and
                 len(cfg['smokeNativeRBins']) == 2,
                 'Only the configured two-R-bin/eight-slice smoke study is in scope')
    export_dir = base.owned(cfg['rawExportDirectory'])
    raw_table = export_dir / 'native_counts.tsv'
    for name in ('native_counts.tsv', 'completed.json', 'provenance.json'):
        stage_inputs[str(export_dir / name)] = base.sha(export_dir / name)
    completion, provenance = (base.read(export_dir / name)
                              for name in ('completed.json', 'provenance.json'))
    validate_export(export_dir, base_cfg, raw_table, completion, provenance)
    axis = base.read(base.owned(base_cfg['inputDirectory']) / 'native_audit.json')['nativeAxis']
    header = base.owned(base_cfg['legacyHeader'])
    stage_inputs[str(header)] = base.sha(header)
    tables = base.tables_from_header(header.read_text())
    _, _, initial_records = base.make_slices(base_cfg, tables, axis, True)
    base.require(len(initial_records) == 8, 'Exactly eight native diagnostic slices required')
    spectra = native_spectra(raw_table, initial_records, axis)
    planned = []
    for row in initial_records:
        species_index = base.SPECIES.index(row['target'])
        mean, sigma = base.parameters(tables, species_index, row['r'])
        mode = locate_mode(spectra[row['id']], mean, sigma, cfg['modeFinder'], axis)
        roi = recentered_roi(mode['modeCenter'], sigma, selection['roiHalfWidthSigma'],
                             base_cfg['fitEdgeAlignment'], axis)
        planned.append(dict(id=row['id'], target=row['target'], mode=row['mode'],
                            rBin=row['rBin'], originalGitMean=mean, originalGitSigma=sigma,
                            originalFitLo=row['fitLo'], originalFitHi=row['fitHi'],
                            **mode, **roi))
    # Existing fresh-plan preparation validates Git, raw cache and live mainconf.
    # Only files just created by this call will subsequently be renamed/replaced.
    for name, expected in stage_inputs.items():
        base.require(base.sha(Path(name)) == expected, 'Recenter input changed during planning: ' + name)
    base.prepare(base_settings, output, True)
    plan = base.read(output / 'plan.json')
    base.require(plan['config'] == base_cfg, 'Base configuration changed during preparation')
    base.unchanged(plan)
    node_hash = base.sha(output / 'nodes.tsv')
    curves_hash = base.sha(output / 'original_curves.tsv')
    with (output / 'slices.tsv').open() as stream:
        old_rows = list(csv.reader(stream, delimiter='\t'))
    rows = old_rows[1:]
    base.require(len(rows) == len(planned) and
                 all(row[0] == record['id'] for row, record in zip(rows, planned)),
                 'Fresh prepare/raw export ordering mismatch')
    (output / 'plan.json').rename(output / 'pre_recenter_plan.json')
    (output / 'slices.tsv').rename(output / 'pre_recenter_slices.tsv')
    new_rows = []
    for row, original, change in zip(rows, plan['slices'], planned):
        updated = list(row)
        updated[4], updated[5] = change['fitLo'], change['fitHi']
        new_rows.append(updated)
        original['fitLo'], original['fitHi'] = change['fitLo'], change['fitHi']
        original['fitWindowAxisClipped'] = change['fitWindowAxisClipped']
    base.tsv(output / 'slices.tsv', old_rows[0], new_rows)
    plan['inputs'][str(output / 'slices.tsv')] = base.sha(output / 'slices.tsv')
    for name in ('pre_recenter_plan.json', 'pre_recenter_slices.tsv'):
        plan['inputs'][str(output / name)] = base.sha(output / name)
    protected_export = {}
    for name in ('native_counts.tsv', 'completed.json', 'provenance.json'):
        source = export_dir / name
        target = output / 'recenter_inputs' / name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(source), str(target))
        base.require(base.sha(source) == base.sha(target), 'Raw export freeze mismatch')
        protected_export[str(source)] = base.sha(source)
        plan['inputs'][str(source)] = base.sha(source)
        plan['inputs'][str(target)] = base.sha(target)
    for source in (settings, Path(__file__).resolve()):
        freeze_extra(plan, output, source)
    base.require(base.sha(output / 'nodes.tsv') == node_hash and
                 base.sha(output / 'original_curves.tsv') == curves_hash,
                 'Original Git nodes or seed curves changed during ROI recentering')
    record = dict(schemaVersion=1, variant=variant, modeFinder=copy.deepcopy(cfg['modeFinder']),
                  roiHalfWidthSigma=selection['roiHalfWidthSigma'],
                  sourceRawExport=str(export_dir), rawSha256=base.sha(raw_table),
                  protectedExport=protected_export, slices=planned,
                  fitCountsModified=False, fitBinRebinning=1,
                  originalNodesUnchanged=True, originalSeedCurvesUnchanged=True,
                  originalDiagnosticGateUnchanged=True,
                  nativeNodesSha256=node_hash, originalCurvesSha256=curves_hash,
                  modeIsPhysicalPidMean=False, coefficientsUpdated=False,
                  purityAvailable=False, productionEnabled=False)
    for name, expected in stage_inputs.items():
        base.require(base.sha(Path(name)) == expected, 'Recenter input changed during preparation: ' + name)
    base.save(output / 'recenter_provenance.json', record)
    plan['inputs'][str(output / 'recenter_provenance.json')] = base.sha(output / 'recenter_provenance.json')
    plan['roiRecenter'] = dict(provenance='recenter_provenance.json', variant=variant,
                              originalNodesUnchanged=True, countsModified=False)
    plan['config']['originalPrepareFitHalfWidthSigma'] = plan['config']['fitHalfWidthSigma']
    plan['config']['fitHalfWidthSigma'] = selection['roiHalfWidthSigma']
    plan['config']['roiCenterPolicy'] = 'empirical_mode_diagnostic_only'
    plan['limitations'].extend([
        'ROI is selected using the same raw data; empirical mode is not an adopted PID mean.',
        'Three-bin smoothing locates only the ROI; all fits use unchanged native counts.',
        'Chi-square values with different ROI domains are not likelihood-ratio tests.',
        'Original Git seed curves, node shapes and diagnostic mass gates remain unchanged.'])
    base.save(output / 'plan.json', plan)
    base.unchanged(plan)
    print('RESTART_RECENTER_PREPARED fits=8 variant=' + variant, output, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings', required=True)
    parser.add_argument('--variant', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    prepare(base.owned(args.settings), base.owned(args.output_dir, False), args.variant)
