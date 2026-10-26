#!/usr/bin/env python3
"""Fresh TOF peak diagnostics seeded only by the explicitly selected Git PID."""
import argparse
import csv
import hashlib
import json
import math
from pathlib import Path
import re
import shutil
import subprocess

import yaml

from audit_restore import polynomial, tables_from_header

PROJECT = Path(__file__).resolve().parents[2]
SPECIES = ('d', 't', 'he3', 'he4')


def require(ok, message):
    if not ok:
        raise ValueError(message)


def owned(value, exists=True):
    value = Path(value)
    result = (value if value.is_absolute() else PROJECT/value).resolve(strict=exists)
    require(PROJECT in result.parents, 'Require an explicit project-owned path: '+str(value))
    return result


def sha(path):
    result = hashlib.sha256()
    with path.open('rb') as stream:
        for block in iter(lambda: stream.read(1048576), b''):
            result.update(block)
    return result.hexdigest()


def read(path):
    return json.loads(path.read_text())


def save(path, value):
    with path.open('x') as stream:
        json.dump(value, stream, indent=2, ensure_ascii=False, allow_nan=False)
        stream.write('\n')


def tsv(path, header, rows):
    with path.open('x', newline='') as stream:
        writer = csv.writer(stream, delimiter='\t', lineterminator='\n')
        writer.writerow(header)
        writer.writerows(rows)


def quote(value):
    value = str(value)
    require(not any(c in value for c in ('"', '\\', '\n')), 'Unsupported ROOT argument quoting')
    return '"'+value+'"'


def parameters(tables, species, rigidity):
    lo, hi = tables['domain'][species]
    require(lo <= rigidity <= hi and rigidity > 0, 'No original-function extrapolation is allowed')
    x = math.log10(rigidity)
    mu = polynomial(tables['mean'][species], x)
    sigma = polynomial(tables['sigma'][species], x)
    require(math.isfinite(mu) and math.isfinite(sigma) and mu > 0 and sigma > 0,
            'Invalid original PID mean/width; no fallback')
    return mu, sigma


def gl3(lo, hi):
    mid, half = (lo+hi)/2., (hi-lo)/2.
    offset = math.sqrt(3./5.)*half
    return ((mid-offset, 5./18.), (mid, 8./18.), (mid+offset, 5./18.))


def make_slices(cfg, tables, axis, smoke):
    require(cfg['quadrature'] == 'gl3', 'Require explicit native-bin quadrature')
    require(cfg['samples'] == dict(d='TAG_D', t='TAG_T', he3='TAG_HE3', he4='TAG_D'),
            'The He4 m2 sample must be exactly the d sample')
    require(cfg['tofWindows']['d'] == cfg['tofWindows']['he4'], 'He4 must share the d m2 window')
    bins = cfg['smokeNativeRBins'] if smoke else sorted(set(range(
        cfg['nativeRBinFirst'], cfg['nativeRBinLast']+1, cfg['nativeRBinStride'])) |
        set(cfg['extraNativeRBins']))
    require(bins and len(set(bins)) == len(bins), 'Nonempty unique native R bins required')
    width, minimum = axis['rBinWidth'], axis['rRange'][0]
    alignment = cfg['fitEdgeAlignment']
    require(alignment > 0 and cfg['fitHalfWidthSigma'] > 0 and cfg['gateNSigma'] > 0,
            'Positive explicit fit/gate geometry required')
    require(abs(alignment/axis['dedxBinWidth']-round(alignment/axis['dedxBinWidth'])) < 1.e-9,
            'Fit alignment must preserve whole native bins')
    slices, nodes, records = [], [], []
    for rbin in bins:
        require(type(rbin) is int and 1 <= rbin <= axis['rBins'], 'Invalid native R bin')
        lo, hi = minimum+(rbin-1)*width, minimum+rbin*width
        shape = {name: [(w, *parameters(tables, s, r)) for r, w in gl3(lo, hi)]
                 for s, name in enumerate(SPECIES)}
        for s, name in enumerate(SPECIES):
            mu, sigma = parameters(tables, s, (lo+hi)/2.)
            rawlo, rawhi = mu-cfg['fitHalfWidthSigma']*sigma, mu+cfg['fitHalfWidthSigma']*sigma
            fitlo = max(axis['dedxRange'][0], math.floor(rawlo/alignment)*alignment)
            fithi = min(axis['dedxRange'][1], math.ceil(rawhi/alignment)*alignment)
            require(fithi > fitlo, 'No saved native axis covers the requested peak window')
            identity = name+'_b'+str(rbin).zfill(4)
            slices.append((identity, cfg['samples'][name], name, rbin, fitlo, fithi,
                           mu-cfg['gateNSigma']*sigma, mu+cfg['gateNSigma']*sigma))
            for other in SPECIES:
                nodes.extend((identity, other, w, m, e) for w, m, e in shape[other])
            records.append(dict(id=identity, target=name, mode=cfg['samples'][name],
                 rBin=rbin, rLo=lo, rHi=hi, r=(lo+hi)/2., fitLo=fitlo, fitHi=fithi,
                 fitWindowAxisClipped=fitlo > rawlo or fithi < rawhi))
    return slices, nodes, records


def select_response_variant(cfg, response_variant):
    if 'responseVariants' in cfg:
        require(response_variant in cfg['responseVariants'], 'Select an explicit response variant from the settings')
        variant = cfg['responseVariants'][response_variant]
        require(set(variant) == {'model', 'alphaStart'}, 'Malformed response variant')
        require(variant['model'] in ('gaussian', 'gaussian_left_exp', 'gaussian_right_exp'),
                'Unknown explicit response family')
        require(type(variant['alphaStart']) in (float, int) and math.isfinite(variant['alphaStart']) and
                cfg['numeric']['responseAlphaMin'] < variant['alphaStart'] < cfg['numeric']['responseAlphaMax'],
                'Response start must be inside the explicitly configured support')
        cfg['numeric'].update(responseModel=variant['model'], responseAlphaStart=variant['alphaStart'])
        cfg['selectedResponseVariant'] = response_variant
    else:
        require(response_variant is None, 'Response variants require an explicit response configuration')


def prepare(settings, output, smoke, response_variant=None):
    cfg = yaml.safe_load(settings.read_text())
    select_response_variant(cfg, response_variant)
    require(cfg['schemaVersion'] == 1 and not output.exists() and output.parent == owned(cfg['outputBase'], False),
            'Require a fresh direct study sibling output')
    header = owned(cfg['legacyHeader'])
    git_header = subprocess.check_output(['git', 'show', cfg['gitRevision']+':'+cfg['legacyHeader']], cwd=PROJECT)
    require(header.read_bytes() == git_header, 'PID header differs from the explicit Git reference')
    tables = tables_from_header(git_header.decode())
    inputs = owned(cfg['inputDirectory'])
    done, context, audit = (read(inputs/name) for name in ('completed.json', 'context.json', 'native_audit.json'))
    require(done['executionComplete'] and done['protectedPostRunVerified'] and
            audit['writtenCacheReadbackExact'] and audit['cacheCountErrorFlowPreserved'], 'Incomplete raw native cache')
    require('2021P24iy_partial_merge' in Path(context['inputFingerprint']['path']).name,
            'Only the original 2021 P24iy input is allowed')
    protected = {str(inputs/name): sha(inputs/name) for name in
                 ('cache.root', 'context.json', 'native_audit.json', 'completed.json')}
    for name in ('cache.root', 'context.json', 'native_audit.json'):
        require(protected[str(inputs/name)] == done['products'][name], 'Raw input hash changed: '+name)
    for name, prefix in (('d', 'D'), ('t', 'T'), ('he3', 'HE3')):
        require(cfg['tofWindows'][name] == [context['settings']['tag'+prefix+'Mass2Min'],
                                          context['settings']['tag'+prefix+'Mass2Max']], 'Saved TOF window differs')
    axis = audit['nativeAxis']
    slices, nodes, records = make_slices(cfg, tables, axis, smoke)
    sources = {settings, Path(__file__).resolve(), Path(__file__).resolve().with_name('audit_restore.py'),
               header, owned(cfg['mainconf']), owned(cfg['wrapper']), owned(cfg['fitSource']),
               owned(cfg['baselineSource']), owned(cfg['fitDependency'])}
    require(type(cfg.get('extraSources', [])) is list, 'extraSources must be an explicit path list')
    sources.update(owned(path) for path in cfg.get('extraSources', []))
    main = yaml.safe_load(owned(cfg['mainconf']).read_text())
    sources.update(owned('config/'+ref) for ref in main.values())
    nuclear = yaml.safe_load(owned('config/'+main['nuclearid']).read_text())
    require(nuclear['dedxCalibration']['enabled'] is False, 'Actual mainconf must use original Git PID')
    output.mkdir(parents=True)
    snapshot = output/'source'
    for original in sources:
        copy = snapshot/original.relative_to(PROJECT)
        copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(str(original), str(copy))
        protected[str(copy)] = sha(copy)
        require(sha(original) == protected[str(copy)], 'Source changed during freeze')
    numeric = dict(cfg['numeric'])
    numeric.update({'hist.'+key: key for key in ('ALL', 'TAG_D', 'TAG_T', 'TAG_HE3')})
    numeric.update(expectedRBins=axis['rBins'], expectedRMin=axis['rRange'][0], expectedRMax=axis['rRange'][1],
                   expectedDedxBins=axis['dedxBins'], expectedDedxMin=axis['dedxRange'][0], expectedDedxMax=axis['dedxRange'][1])
    tsv(output/'settings.tsv', ('key', 'value'), numeric.items())
    tsv(output/'slices.tsv', ('id', 'mode', 'target', 'rBin', 'fitLo', 'fitHi', 'gateLo', 'gateHi'), slices)
    tsv(output/'nodes.tsv', ('id', 'species', 'weight', 'mean', 'sigma'), nodes)
    curves = []
    for s, name in enumerate(SPECIES):
        for k in range(axis['rBins']):
            rigidity = axis['rRange'][0]+(k+.5)*axis['rBinWidth']
            if tables['domain'][s][0] <= rigidity <= tables['domain'][s][1]:
                curves.append((name, rigidity, *parameters(tables, s, rigidity)))
    tsv(output/'original_curves.tsv', ('species', 'r', 'mean', 'sigma'), curves)
    plots = dict(cfg['baselinePlots'])
    for name, prefix in (('d', 'D'), ('t', 'T'), ('he3', 'He3')):
        plots['tag'+prefix+'Mass2Min'], plots['tag'+prefix+'Mass2Max'] = cfg['tofWindows'][name]
    tsv(output/'baseline_settings.tsv', ('key', 'value'), plots.items())
    for source in output.glob('*.tsv'):
        protected[str(source)] = sha(source)
    save(output/'plan.json', dict(schemaVersion=1, config=cfg, smoke=smoke, slices=records,
         totalFits=len(slices), inputs=protected, sourceSha256={str(p): sha(p) for p in sources},
         cache=str(inputs/'cache.root'), inputFingerprint=context['inputFingerprint'],
         productionEnabled=False, coefficientsUpdated=False, purityAvailable=False,
         limitations=['Raw input is charge-inclusive and centrality-integrated.',
                      'TAG_D includes d and He4; TAG and ALL are not independent.',
                      'Initial survey covers R where all four original functions exist; it defines no cutoff.',
                      'Local free-target peak fits are diagnostics, not an adopted response or purity.']))
    print('RESTART_TAG_PREPARED fits='+str(len(slices)), output, flush=True)


def unchanged(plan):
    for name, expected in plan['inputs'].items():
        require(sha(Path(name)) == expected, 'Protected input changed: '+name)
    for name, expected in plan['sourceSha256'].items():
        require(sha(Path(name)) == expected, 'Live source/configuration changed: '+name)


def run(output, baseline):
    plan = read(output/'plan.json')
    unchanged(plan)
    cfg = plan['config']
    wrapper, main = owned(cfg['wrapper']), owned(cfg['mainconf'])
    require(sha(wrapper) == plan['sourceSha256'][str(wrapper)] and
            sha(main) == plan['sourceSha256'][str(main)], 'Actual runtime changed')
    if baseline:
        script = output/'source'/cfg['baselineSource']
        report = output/'git_pid_baseline.pdf'
        invocation = str(script)+'+('+','.join(quote(p) for p in
             (plan['cache'], output/'original_curves.tsv', output/'baseline_settings.tsv', report))+')'
        marker, logname = 'RESTART_BASELINE_PDF_CLOSED pages=8', 'baseline.log'
    else:
        script = output/'source'/cfg['fitSource']
        invocation = str(script)+'+('+','.join(quote(p) for p in
             (plan['cache'], output, output/'settings.tsv', output/'slices.tsv', output/'nodes.tsv'))+')'
        marker, logname = 'RESTART_PEAK_READBACK_OK slices='+str(plan['totalFits']), 'fit.log'
    logpath = output/logname
    require(not logpath.exists(), 'An earlier attempt exists; preserve it and choose a fresh output')
    command = ['bash', str(wrapper), str(main), 'root4star', '-l', '-b', '-q', invocation]
    save(output/(logname+'.command.json'), dict(cwd=str(PROJECT), argv=command))
    with logpath.open('x') as stream:
        subprocess.run(command, cwd=PROJECT, stdout=stream, stderr=subprocess.STDOUT,
                       timeout=cfg['executionTimeoutSeconds'], check=True)
    require(marker in logpath.read_text(), 'Missing closed-output completion marker; inspect '+str(logpath))
    if baseline:
        checked = subprocess.run(['gs', '-dBATCH', '-dNOPAUSE', '-dPDFSTOPONERROR',
                                  '-sDEVICE=nullpage', str(report)], check=True,
                                 stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True)
        pages = re.search(r'Processing pages 1 through ([0-9]+)\.', checked.stdout)
        require(pages and int(pages.group(1)) == 8 and 'Error:' not in checked.stdout,
                'Closed baseline PDF page/readback verification failed')
        with (output/'baseline_readback.log').open('x') as stream:
            stream.write(checked.stdout)
    else:
        result = read(output/'fits.json')
        require(result['rawReadback'] and result['predictionReadback'] and
                result['productionEnabled'] is False and result['purityAvailable'] is False and
                result['totalSlices'] == plan['totalFits'] and len(result['results']) == plan['totalFits'],
                'Incomplete or inappropriate calibration diagnostic output')
        for row, expected in zip(result['results'], plan['slices']):
            require(all(row[key] == expected[key] for key in ('id', 'mode', 'target', 'rBin')),
                    'Fit output changed a requested native R slice')
            require(row['curveAdoptionAllowed'] is False and row['purityAvailable'] is False,
                    'A peak diagnostic must not adopt calibration or purity')
    unchanged(plan)
    save(output/(logname+'.completed.json'), dict(status='completed', inputsUnchanged=True,
         productionEnabled=False, coefficientsUpdated=False, purityAvailable=False))
    print('RESTART_'+('BASELINE' if baseline else 'TAG')+'_COMPLETE', output, flush=True)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('action', choices=('prepare', 'fit', 'baseline'))
    parser.add_argument('--settings')
    parser.add_argument('--output-dir', required=True)
    parser.add_argument('--smoke', action='store_true')
    parser.add_argument('--response-variant', help='Explicit variant key from response settings')
    args = parser.parse_args()
    output = owned(args.output_dir, args.action != 'prepare')
    if args.action == 'prepare':
        require(args.settings, 'Explicit settings required')
        prepare(owned(args.settings), output, args.smoke, args.response_variant)
    else:
        require(args.response_variant is None and args.settings is None and not args.smoke,
                'Run uses the frozen plan; response/settings/smoke are prepare-only options')
        run(output, args.action == 'baseline')


if __name__ == '__main__':
    main()
