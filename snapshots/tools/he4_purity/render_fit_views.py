#!/usr/bin/env python3
"""Re-render verified saved He4 fits on linear and physical dE/dx axes; no refit."""
import argparse
from datetime import datetime, timedelta, timezone
import json
import math
from pathlib import Path
import shutil
import subprocess

import yaml

from run_purity import digest, fingerprint, page_count, save, HERE, PROJECT


def runtime_provenance(context, mainconf, allow_drift=False):
    """Keep strict defaults; allow explicit same-path ROOT-runtime-only drift."""
    saved_path = Path(context['mainconf']).resolve()
    current_sha = digest(mainconf)
    if saved_path != mainconf:
        raise ValueError('Explicit mainconf path differs from the saved study.')
    drift = context['mainconf_sha256'] != current_sha
    if drift and not allow_drift:
        raise ValueError('Mainconf differs from saved study; runtime-only drift requires explicit opt-in.')
    current = yaml.safe_load(mainconf.read_text())
    analysis_ref = current.get('analysis') if isinstance(current, dict) else None
    if not isinstance(analysis_ref, str) or not analysis_ref:
        raise ValueError('Runtime mainconf has no explicit analysis reference.')
    config_root = PROJECT/'config'
    analysis_path = (config_root/analysis_ref).resolve(strict=True)
    if config_root.resolve() not in analysis_path.parents:
        raise ValueError('Runtime analysis_info must reside below project config/.')
    info = yaml.safe_load(analysis_path.read_text())
    tag, analysis = info.get('starTag', {}), info.get('analysis', {})
    identity = {key: tag.get(key) for key in
                ('libraryTag', 'year', 'productionTag', 'collision', 'triggerSets')}
    identity['mode'] = analysis.get('mode')
    if drift:
        saved_input = context.get('input', {}).get('path', '')
        sample_input = context.get('settings', {}).get('samples', {}).get('3p85', {}).get('input', '')
        if (context.get('sample') != '3p85' or '2021P24iy' not in saved_input
                or '2021P24iy' not in sample_input
                or identity != dict(libraryTag='SL24y', year='year2021',
                    productionTag='P24iy', collision='auau3p85',
                    triggerSets='production_3p85GeV_fixedTarget_2021', mode='fxtmult')):
            raise ValueError('Runtime drift requires verified 3p85/2021P24iy/FXT/SL24y identity.')
    return dict(mode='explicit_runtime_only_drift' if drift else 'saved_mainconf_exact',
                allow_runtime_mainconf_drift=allow_drift, drift_detected=drift,
                saved_mainconf=str(saved_path), saved_mainconf_sha256=context['mainconf_sha256'],
                runtime_mainconf=str(mainconf), runtime_mainconf_sha256=current_sha,
                runtime_analysis_info=str(analysis_path), runtime_analysis_info_sha256=digest(analysis_path),
                runtime_identity=identity,
                note='Mainconf selects ROOT/STAR runtime only; current analysis cuts/PID are not loaded or reapplied.')


def verified_previous_views(study):
    """Do not replace or modify the frozen prior view revision."""
    folder = study/'views'
    completed_path = folder/'completed.json'
    if not completed_path.exists():
        return {}
    completed = json.loads(completed_path.read_text())
    if completed.get('phase') != 'complete':
        raise ValueError('Prior views are not completed frozen products.')
    result = {str(completed_path): digest(completed_path)}
    for name, expected in completed['products'].items():
        path = (folder/name).resolve(strict=True)
        if folder.resolve() not in path.parents or digest(path) != expected:
            raise ValueError('Prior frozen view product hash mismatch: '+name)
        result[str(path)] = expected
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study-folder', type=Path, required=True)
    parser.add_argument('--mainconf', type=Path, required=True)
    parser.add_argument('--settings', type=Path, default=HERE/'view_settings.yaml')
    parser.add_argument('--output-dir', type=Path)
    parser.add_argument('--allow-runtime-mainconf-drift', action='store_true',
                        help='Explicit same-path opt-in: saved 3p85/2021P24iy fits use '
                             'current SL24y config for ROOT runtime only, not analysis cuts.')
    args = parser.parse_args()
    study = args.study_folder.resolve(strict=True)
    mainconf = args.mainconf.resolve(strict=True)
    settings = args.settings.resolve(strict=True)
    cfg = yaml.safe_load(settings.read_text())
    for key in ('fullDisplayNSigmaMin', 'fullDisplayNSigmaMax',
                'physicalDisplayNSigmaMin', 'physicalDisplayNSigmaMax',
                'linearMaximumFactor', 'logMaximumFactor', 'logMinimum'):
        if not isinstance(cfg[key], (int, float)) or not math.isfinite(cfg[key]):
            raise ValueError('Invalid view setting: '+key)
    if not (cfg['schemaVersion'] == 1 and cfg['panelsPerPage'] == 4
            and cfg['diagnosticPanelsPerPage'] in (2, 4)
            and isinstance(cfg['showNeighborComponents'], bool)
            and isinstance(cfg.get('requireNeighborCentersVisible', False), bool)
            and cfg['fullDisplayNSigmaMin'] < cfg['fullDisplayNSigmaMax']
            and cfg['physicalDisplayNSigmaMin'] < cfg['physicalDisplayNSigmaMax']
            and cfg['linearMaximumFactor'] > 1 and cfg['logMaximumFactor'] > 1
            and cfg['logMinimum'] > 0):
        raise ValueError('Invalid range or scale settings.')
    diagnostic_name = cfg['diagnosticPdf']
    if (not isinstance(diagnostic_name, str) or Path(diagnostic_name).name != diagnostic_name
            or not diagnostic_name.endswith('.pdf') or diagnostic_name == 'he4_fits_linear.pdf'):
        raise ValueError('A distinct local PDF filename is required for diagnostics.')
    centers = cfg['diagnosticRigidityCenters']
    if not centers or any(not isinstance(x, (float, int)) or not math.isfinite(x)
                          for x in centers) or len(set(centers)) != len(centers):
        raise ValueError('Unique finite diagnostic rigidity centers are required.')
    completed_path, context_path, audit_path = (study/n for n in
                                              ('completed.json', 'context.json', 'purity.json'))
    completed = json.loads(completed_path.read_text())
    if completed['phase'] != 'complete':
        raise ValueError('Only completed saved fit studies may be rendered.')
    verified = {}
    for name in ('context.json', 'purity.json', 'he4_fits.root'):
        path = study/name
        verified[name] = digest(path)
        if completed['products'].get(name) != verified[name]:
            raise ValueError('Saved fit product hash mismatch: '+name)
    context, audit = (json.loads(p.read_text()) for p in (context_path, audit_path))
    runtime = runtime_provenance(context, mainconf, args.allow_runtime_mainconf_drift)
    saved_sources = {}
    for name, checksum in context.get('sources', {}).items():
        path = Path(name).resolve(strict=True)
        if study not in path.parents or digest(path) != checksum:
            raise ValueError('Saved fit source snapshot hash mismatch: '+name)
        saved_sources[str(path)] = checksum
    previous_views = verified_previous_views(study)
    source = Path(audit['inputRoot']).resolve(strict=True)
    initial_source = fingerprint(source)
    if context['input'] != initial_source:
        raise ValueError('Original ROOT differs from the saved study fingerprint.')
    fit_root = study/'he4_fits.root'
    initial_fit = fingerprint(fit_root)
    slices = audit['slices']
    diagnostic_indices = []
    for r in centers:
        candidates = [i for i, s in enumerate(slices) if abs(s['r']-r) < 1.e-7]
        if len(candidates) != 1:
            raise ValueError('Requested center is not a unique saved slice: '+str(r))
        diagnostic_indices.append(candidates[0])
    diagnostic_indices.sort()  # Pages follow increasing saved rigidity.
    output = (args.output_dir or study/'views').resolve()
    if output == study or (output.exists() and any(output.iterdir())):
        raise ValueError('Choose a new empty output directory to preserve saved products.')
    output.mkdir(parents=True, exist_ok=True)
    (output/'logs').mkdir()
    (output/'source').mkdir()
    macro = output/'source/RenderHe4FitViews.C'
    shutil.copy2(HERE/'RenderHe4FitViews.C', macro)
    shutil.copy2(Path(__file__).resolve(), output/'source/render_fit_views.py')
    shutil.copy2(context_path, output/'source/study_context.json')
    shutil.copy2(completed_path, output/'source/study_completed.json')
    (output/'runtime_config').mkdir()
    shutil.copy2(mainconf, output/'runtime_config/mainconf.yaml')
    shutil.copy2(runtime['runtime_analysis_info'], output/'runtime_config/analysis_info.yaml')
    shutil.copy2(settings, output/'settings.yaml')
    numeric = {key: int(value) if isinstance(value, bool) else value for key, value in cfg.items()
               if key not in ('diagnosticRigidityCenters', 'diagnosticPdf')}
    numeric['gateNSigma'] = audit['gateNSigma']
    manifest = output/'settings.tsv'
    manifest.write_text(''.join('%s\t%s\n' % pair for pair in numeric.items()))
    rows = []
    for i, item in enumerate(slices):
        values = [i, item['lo'], item['hi'], item['calibrationMean'],
                  item['calibrationSigma'], item['chosen'], int(i in diagnostic_indices)]
        for model in ('linear', 'exponential'):
            fit = next((f for f in item['fits'] if f['model'] == model), None)
            values += ([1, int(fit['valid']), int(fit['quality']), fit['purity'],
                        fit['purityError'], fit['chi2'], fit['ndf'], fit['status'],
                        fit['covarianceStatus'], fit['aic']] if fit else [0, 0, 0, -1, 0, 0, 0, -1, -1, 0])
            neighbors = fit.get('neighborSpecies', []) if fit else []
            if (any(type(s) is not int or s not in (0, 1, 2) for s in neighbors)
                    or len(neighbors) != len(set(neighbors))):
                raise ValueError('Invalid neighbor species metadata.')
            if (cfg['showNeighborComponents'] and fit
                    and len(fit['parameters']) != 5+3*len(neighbors)):
                raise ValueError('Species metadata does not account for all fitted components.')
            values += [len(neighbors)]+neighbors
        rows.append('\t'.join(str(v) for v in values)+'\n')
    slice_manifest = output/'slices.tsv'
    slice_manifest.write_text(''.join(rows))
    initial = dict(fit_root=initial_fit, original_root=initial_source, study=verified,
                   completed_sha256=digest(completed_path))
    save(output/'context.json', dict(study_folder=str(study), inputs=initial,
        mainconf=str(mainconf), mainconf_sha256=digest(mainconf), settings=cfg,
        runtime_provenance=runtime, saved_source_snapshots=saved_sources,
        previous_view_products=previous_views,
        created_jst=datetime.now(timezone(timedelta(hours=9))).isoformat(),
        note='Stored fits are re-rendered, never fitted again; physical curves preserve stored fit-bin heights.'))
    call = str(macro)+'+('+','.join(json.dumps(str(p)) for p in
          (fit_root, output, manifest, slice_manifest, audit['sampleLabel'], diagnostic_name))+')'
    command = ['bash', str(PROJECT/'script/singularity_root_batch.sh'), str(mainconf),
               'root4star', '-l', '-b', '-q', call]
    log = output/'logs/render.log'
    print('RENDER', output, flush=True)
    with log.open('w') as handle:
        handle.write('COMMAND: '+json.dumps(command)+'\n'); handle.flush()
        result = subprocess.run(command, cwd=str(PROJECT), stdout=handle,
                                stderr=subprocess.STDOUT, timeout=900)
    content = log.read_text(errors='replace')
    if result.returncode or 'HE4_FIT_VIEWS_COMPLETE' not in content or any(t in content for t in
       ('Error in <ACLiC>', 'segmentation violation', 'fatal error:', 'Error: Function')):
        raise RuntimeError('ROOT view rendering failed; see '+str(log))
    render = json.loads((output/'render_audit.json').read_text())
    pages = dict(linear=page_count(output/'he4_fits_linear.pdf'),
                 physical=page_count(output/diagnostic_name))
    if pages != dict(linear=(len(slices)+3)//4, physical=len(centers)):
        raise ValueError('Unexpected PDF page counts.')
    if render['slices'] != len(slices) or render['diagnostic_indices'] != diagnostic_indices:
        raise ValueError('Rendered selection differs from the explicit slice manifest.')
    if render['diagnostic_panels_per_page'] != cfg['diagnosticPanelsPerPage']:
        raise ValueError('Rendered panel layout differs from explicit settings.')
    if (render['individual_neighbor_components'] != cfg['showNeighborComponents']
            or (cfg['showNeighborComponents'] and not render['component_sum_equivalence_checked'])):
        raise ValueError('Component drawing/equivalence checks did not match settings.')
    for item in render['diagnostic_components']:
        stored, physical = item['stored_dedx_range'], item['physical_display_dedx']
        if not stored[0] <= physical[0] < physical[1] <= stored[1]:
            raise ValueError('Physical view extends outside the saved dE/dx axis.')
        same_range = (cfg['fullDisplayNSigmaMin'] == cfg['physicalDisplayNSigmaMin']
                      and cfg['fullDisplayNSigmaMax'] == cfg['physicalDisplayNSigmaMax'])
        if same_range and not (item['same_physical_coverage'] and item['same_visible_bins']
                               and item['same_visible_counts']):
            raise ValueError('Upper/lower panels differ in physical/bin/count coverage.')
        if (cfg.get('requireNeighborCentersVisible', False)
                and not all(neighbor['visible_physical'] for neighbor in item['neighbor_centers'])):
            raise ValueError('A stored neighboring Gaussian center is outside the physical view.')
    if (fingerprint(source) != initial_source or fingerprint(fit_root) != initial_fit
            or digest(completed_path) != initial['completed_sha256']
            or any(digest(study/name) != value for name, value in verified.items())
            or any(digest(Path(name)) != value for name, value in saved_sources.items())
            or any(digest(Path(name)) != value for name, value in previous_views.items())
            or digest(mainconf) != runtime['runtime_mainconf_sha256']
            or digest(Path(runtime['runtime_analysis_info'])) != runtime['runtime_analysis_info_sha256']):
        raise ValueError('A saved input changed during rendering.')
    (output/'README.md').write_text(
        '# 保存済み⁴He fitの追加表示\n\n'
        '[全剛性bin・linear scale](he4_fits_linear.pdf) / '
        '[低剛性・log/linear比較]('+diagnostic_name+')\n\n'
        '- 元study: `'+str(study)+'`。追加のfitは行っていません。\n'
        '- 全binのlinear図は元の校正nSigma軸。低剛性比較は物理dE/dx軸（keV/cm）です。\n'
        + ('- 低剛性PDFは各binを4面で表示：校正nSigma軸のlog/linearと、物理dE/dx軸のlog/linear。\n'
           if cfg['diagnosticPanelsPerPage'] == 4 else '') +
        ('- 上下とも校正nSigma '+str(cfg['fullDisplayNSigmaMin'])+'〜'
         +str(cfg['fullDisplayNSigmaMax'])+' の同じ物理範囲です。保存軸0〜120 keV/cmに達する場合は、'
         '両軸とも同じbinまでclipします（拡大表示ではありません）。\n' if same_range else '') +
        ('- 現在mainconfの差異を明示的に許可したruntime専用実行です。旧/現在SHAと設定全文を保存し、'
         'STAR SL24y・2021 P24iy・3.85 GeV FXTを照合しました。現在PID/cutは元fitへ再適用していません。\n'
         if runtime['drift_detected'] else '') +
        '- 黒点は保存射影、赤/青はlinear/exponential背景を使った全fit曲線、'
        '緑は採用モデルの⁴He Gaussian、紫はその背景です。\n'
        + ('- T・³HeのGaussian背景を個別表示し、D tailが含まれる場合も区別します。'
           '不採用fitしかないbinでは最低AICモデルの成分を診断表示し、`unaccepted`と明記します。\n'
           if cfg['showNeighborComponents'] else '') +
        '- 校正窓±'+str(audit['gateNSigma'])+'σを灰点線、採用Gaussian中心を緑破線で表示。'
        'fit曲線は座標だけを変換し、保存済みfit用binあたりの高さを保っています。'
        'bin幅は各図の縦軸・脚注と`render_audit.json`に記録しています。\n'
        '- `No accepted fit` はfitの品質・境界条件で不採用になったことを意味し、'
        'ピークが存在しないとの判定ではありません。\n'
        '- TPC前PID、正負電荷・全採用中央度の合算。最終⁴He候補やTOF併用後の純度ではありません。\n'
        '- y軸は元fit用binあたりの計数です。座標変換によるJacobian倍率は掛けず、関数値の一致を数値検証しました。\n'
        '- `completed.json`のハッシュと入力fingerprintを照合し、元ROOT・fitソースsnapshot・旧PDFの変更なしを確認しました。\n')
    save(output/'completed.json', dict(phase='complete', no_refit=True,
        inputs_unchanged=True, pages=pages,
        products={str(p.relative_to(output)):digest(p) for p in output.rglob('*')
                  if p.is_file() and p.name != 'completed.json' and p.suffix not in ('.so', '.d', '.o')}))
    print('COMPLETE', output, pages, flush=True)


if __name__ == '__main__':
    main()
