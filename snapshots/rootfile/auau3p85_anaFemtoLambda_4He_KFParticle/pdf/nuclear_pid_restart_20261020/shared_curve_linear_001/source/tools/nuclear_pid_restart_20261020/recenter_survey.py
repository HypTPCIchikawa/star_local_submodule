#!/usr/bin/env python3
"""Prepare populated native survey ROIs, retaining explicit unassessed records."""
import argparse
import csv
import math
from pathlib import Path
import shutil
import yaml

import prepare_tag_study as base
import recenter_tag_study as local


def selected_tables(slices, nodes, changes):
    identities = {row[0] for row in slices}
    base.require(set(changes) <= identities and changes, 'Unknown or empty selection')
    selected = []
    for row in slices:
        if row[0] in changes:
            new = list(row)
            new[4], new[5] = changes[row[0]]['fitLo'], changes[row[0]]['fitHi']
            selected.append(new)
    chosen_nodes = [row for row in nodes if row[0] in changes]
    base.require({row[0] for row in chosen_nodes} == set(changes), 'Missing selected shape nodes')
    return selected, chosen_nodes


def prepare(settings, output):
    cfg = yaml.safe_load(settings.read_text())
    local.validate_mode_settings(cfg['modeFinder'])
    half = cfg['roiHalfWidthSigma']
    base.require(cfg['schemaVersion'] == 1 and type(half) in (float, int) and
                 math.isfinite(half) and half > 0, 'Require an explicit positive ROI width')
    base_settings = base.owned(cfg['baseSettings'])
    original_cfg = yaml.safe_load(base_settings.read_text())
    base.require(Path(original_cfg['fitSource']).name == 'FitSplitNormalSlices.C' and
                 original_cfg['numeric']['storageMapping'] == 'root_default12', 'Wrong accepted response')
    inputs = base.owned(original_cfg['inputDirectory'])
    axis = base.read(inputs/'native_audit.json')['nativeAxis']
    header = base.owned(original_cfg['legacyHeader'])
    tables = base.tables_from_header(header.read_text())
    _, _, original_records = base.make_slices(original_cfg, tables, axis, False)
    export = base.owned(cfg['rawExportDirectory'])
    raw = export/'native_counts.tsv'
    done, provenance = (base.read(export/name) for name in ('completed.json', 'provenance.json'))
    base.require(done['status'] == 'completed' and done['nativeReadbackExact'] and
                 done['inputHashesUnchanged'] and done['noFit'] and not done['productionEnabled'] and
                 done['tableSha256'] == base.sha(raw), 'Native export integrity failed')
    base.require(done['summary']['sliceCount'] == len(original_records) and
                 done['summary']['countErrorConsistent'] and
                 done['summary']['nonintegerCountBins'] == 0, 'Native export contract differs')
    source_plan = base.read(base.owned(provenance['sourceFit'])/'plan.json')
    cache = inputs/'cache.root'
    base.require(Path(source_plan['cache']).resolve() == cache and
                 source_plan['inputs'][str(cache)] == base.sha(cache) and
                 [r['id'] for r in source_plan['slices']] == [r['id'] for r in original_records],
                 'Native grid or source cache differs')
    for name in ('cache.root', 'context.json', 'native_audit.json', 'completed.json'):
        base.require(provenance['protected'][str(inputs/name)] == base.sha(inputs/name),
                     'Raw provenance changed: '+name)
    sources = (settings, base_settings, Path(__file__).resolve(), Path(local.__file__).resolve(),
               header, raw, export/'completed.json', export/'provenance.json')
    protected = {str(path): base.sha(path) for path in sources}
    spectra = local.native_spectra(raw, original_records, axis)
    attempts, changes = [], {}
    for row in original_records:
        mean, sigma = base.parameters(tables, base.SPECIES.index(row['target']), row['r'])
        record = dict(id=row['id'], target=row['target'], mode=row['mode'], rBin=row['rBin'],
                      r=row['r'], originalGitMean=mean, originalGitSigma=sigma,
                      fitRequested=False, assessment='unassessed', lowPurityInferred=False)
        try:
            mode = local.locate_mode(spectra[row['id']], mean, sigma, cfg['modeFinder'], axis)
        except ValueError as error:
            if not str(error).startswith(('Empty search interval:', 'No complete three-bin')):
                raise
            record['reason'] = str(error)
            attempts.append(record)
            continue
        roi = local.recentered_roi(mode['modeCenter'], sigma, half,
                                   original_cfg['fitEdgeAlignment'], axis)
        known = {name: base.parameters(tables, index, row['r'])[0]
                 for index, name in enumerate(base.SPECIES)}
        record.update(mode, **roi)
        record.update(fitRequested=True, assessment='fit_requested',
            originalForeignCentersInsideRoi=[name for name, m in known.items()
                if name != row['target'] and roi['fitLo'] <= m <= roi['fitHi']],
            nearestOriginalCurveToMode=min(known, key=lambda name:
                (abs(known[name]-mode['modeCenter']), base.SPECIES.index(name))),
            supportReviewRequired=mode['searchBoundaryMode'] or roi['fitWindowAxisClipped'])
        changes[row['id']] = record
        attempts.append(record)
    base.require(changes, 'No populated ROI; no fabricated fit')
    for name, expected in protected.items():
        base.require(base.sha(Path(name)) == expected, 'Planning input changed: '+name)
    base.prepare(base_settings, output, False)
    plan = base.read(output/'plan.json')
    base.unchanged(plan)
    base.require(plan['config'] == original_cfg, 'Base configuration changed')
    tables_saved = {}
    for name in ('slices.tsv', 'nodes.tsv'):
        with (output/name).open() as stream:
            tables_saved[name] = list(csv.reader(stream, delimiter='\t'))
    new_slices, new_nodes = selected_tables(tables_saved['slices.tsv'][1:],
                                           tables_saved['nodes.tsv'][1:], changes)
    curves_hash = base.sha(output/'original_curves.tsv')
    for name in ('plan.json', 'slices.tsv', 'nodes.tsv'):
        (output/name).rename(output/('pre_recenter_'+name))
        plan['inputs'][str(output/('pre_recenter_'+name))] = base.sha(output/('pre_recenter_'+name))
    base.tsv(output/'slices.tsv', tables_saved['slices.tsv'][0], new_slices)
    base.tsv(output/'nodes.tsv', tables_saved['nodes.tsv'][0], new_nodes)
    plan['slices'] = [dict(row, fitLo=changes[row['id']]['fitLo'], fitHi=changes[row['id']]['fitHi'],
                         fitWindowAxisClipped=changes[row['id']]['fitWindowAxisClipped'])
                      for row in plan['slices'] if row['id'] in changes]
    plan['totalFits'] = len(new_slices)
    for name in ('slices.tsv', 'nodes.tsv'):
        plan['inputs'][str(output/name)] = base.sha(output/name)
    for name in ('native_counts.tsv', 'completed.json', 'provenance.json'):
        source, target = export/name, output/'recenter_inputs'/name
        target.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(source), str(target))
        base.require(base.sha(source) == base.sha(target), 'Raw export freeze mismatch')
        plan['inputs'][str(source)], plan['inputs'][str(target)] = base.sha(source), base.sha(target)
    for source in (settings, Path(__file__).resolve(), Path(local.__file__).resolve()):
        local.freeze_extra(plan, output, source)
    base.require(base.sha(output/'original_curves.tsv') == curves_hash, 'Original curves changed')
    record = dict(schemaVersion=1, config=cfg, sourceRawExport=str(export), rawSha256=base.sha(raw),
        requestedSlices=len(original_records), fittedSlices=len(new_slices), slices=attempts,
        unassessedSlices=[row['id'] for row in attempts if not row['fitRequested']],
        fitCountsModified=False, selectedOriginalNodeValuesUnchanged=True,
        originalSeedCurvesUnchanged=True, originalDiagnosticGateUnchanged=True,
        modeIsPhysicalPidMean=False, productionEnabled=False, purityAvailable=False)
    for name, expected in protected.items():
        base.require(base.sha(Path(name)) == expected, 'Planning input changed: '+name)
    base.save(output/'recenter_provenance.json', record)
    plan['inputs'][str(output/'recenter_provenance.json')] = base.sha(output/'recenter_provenance.json')
    plan['recenterSurvey'] = dict(provenance='recenter_provenance.json', requestedSlices=len(original_records),
                                 fittedSlices=len(new_slices), unassessedSlices=record['unassessedSlices'])
    plan['limitations'].extend(['ROI only: original Git seed node values and raw native counts unchanged.',
        'Empty peak search is unassessed, never zero purity or an isolated veto.'])
    base.save(output/'plan.json', plan)
    base.unchanged(plan)
    print('RESTART_RECENTER_SURVEY_PREPARED requested='+str(len(original_records))+
          ' fits='+str(len(new_slices)), output, flush=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--settings', required=True)
    parser.add_argument('--output-dir', required=True)
    args = parser.parse_args()
    prepare(base.owned(args.settings), base.owned(args.output_dir, False))
