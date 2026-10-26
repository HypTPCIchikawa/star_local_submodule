#!/usr/bin/env python3
"""Produce audited KF-only adaptations; never modify the legacy source macros."""
import argparse
import hashlib
import json
import math
from datetime import datetime, timezone
from pathlib import Path
import re
import yaml

ROOT = Path(__file__).resolve().parents[2]
TASKS = [
    ('Draw_LambdaInvMass.C', 'Lambda_InvMass_legacy.pdf', 1),
    ('draw_LambdaNuclearId_relamom_sum.C', 'kstar_sum_4He.pdf', 12),
    ('draw_LambdaNuclearId_relamom_sum_CFsub.C', 'kstar_sum_CFsub_4He.pdf', 12),
    ('draw_LambdaNuclearId_relamom_sum_CFsub_method3.C', 'kstar_sum_CFsub_method3_4He.pdf', 6),
    ('draw_LambdaNuclearId_relamom_sum_CFsub_method5.C', 'kstar_sum_CFsub_method5_4He.pdf', 6),
    ('draw_LambdaNuclearId_relamom_CFsub_method3.C', 'kstar_cent_CFsub_method3_4He.pdf', 3),
    ('draw_LambdaNuclearId_relamom_CFsub_method5.C', 'kstar_cent_CFsub_method5_4He.pdf', 3),
    ('Fit_LambdaNuclearId_relamom_sum.C', 'Fit_LL_4He.pdf', 2),
    ('Fit_LambdaNuclearId_relamom_sum_CFsub.C', 'Fit_LL_CFsub_4He.pdf', 2),
    ('Fit_LambdaNuclearId_relamom_sum_CFsub_method3.C', 'Fit_LL_method3_4He.pdf', 2),
    ('Fit_LambdaNuclearId_relamom_sum_CFsub_method5.C', 'Fit_LL_method5_4He.pdf', 2),
    ('draw_LambdaNuclearId_relamom.C', 'kstar_4He.pdf', 6),
    ('draw_LambdaNuclearId_relamom_CFsub.C', 'kstar_CFsub_4He.pdf', 6),
]

def sha(data):
    return hashlib.sha256(data).hexdigest()

def cpp(value):
    return json.dumps(str(value), ensure_ascii=True)

def num(value):
    return format(float(value), '.17g')

def stable_manifest_bytes(manifest, existing=None):
    """Preserve equivalent old serialization, including pre-canonical header paths.

    Python 3.8 keeps relative __file__ when invoked with a relative script path;
    absolute and relative invocations must not count as a provenance change.
    Other metadata, hashes and values must still agree exactly.
    """
    if existing is not None:
        old = json.loads(existing.decode())
        source = Path(old['support']['source'])
        if not source.is_absolute():
            source = ROOT/source
        old['support']['source'] = str(source.resolve())
        if old == manifest:
            return existing
    return (json.dumps(manifest,indent=2)+'\n').encode()

def close_brace(text, pos):
    depth = 1
    while depth:
        if pos >= len(text):
            raise ValueError('Unbalanced generated source')
        depth += (text[pos] == '{') - (text[pos] == '}')
        pos += 1
    return pos

def adapt(original, name, input_root, outdir, purity, window, settings):
    text, changes = original, []
    def sub(label, pattern, replacement, expected=None, flags=0):
        nonlocal text
        text, count = re.subn(pattern, replacement, text, flags=flags)
        if expected is not None and count != expected:
            raise ValueError('%s: %s expected %s replacements, got %s' % (name, label, expected, count))
        if count:
            changes.append(dict(label=label, count=count))
    # Legacy commented literals contain semicolons, so use the next statement.
    sub('input_root', r'TString inputFile\s*=.*?(?=\s*TFile\s*\*fIn)',
        lambda _: 'TString inputFile = %s;\n' % cpp(input_root), 1, re.S)
    sub('output_directory', r'TString outDir\s*=.*?(?=\s*gSystem->mkdir)',
        lambda _: 'TString outDir = %s;\n' % cpp(outdir), 1, re.S)
    sub('4He_only', r'const char\s*\*species\[\]\s*=\s*\{[^}]*\};', 'const char *species[] = {"4He"};')
    sub('one_species', r'int nSpecies\s*=\s*4;', 'int nSpecies = 1;')
    sub('KF_true_directory', r'fIn->Get\("nuclearid"\)', 'fIn->Get("true")')
    sub('generic_TH2', r'\bTH2D\b', 'TH2')
    if name == 'Draw_LambdaInvMass.C':
        sub('separate_legacy_mass_pdf', r'Lambda_InvMass\.pdf', 'Lambda_InvMass_legacy.pdf', 1)
    group = re.search(r'if\s*\(groupSet\s*==\s*0\)\s*\{\s*centLabel', text)
    if group:
        end_first = close_brace(text, text.index('{', group.start())+1)
        second = re.match(r'\s*else\s*\{', text[end_first:])
        if not second:
            raise ValueError(name + ': unknown group declaration')
        end = close_brace(text, end_first+second.end())
        blocks = []
        for group_set in settings['group_sets']:
            blocks.append('\n'.join('%s (group == %d) { centLabel = %s; jStart = %d; jEnd = %d; }' %
                ('if' if n == 0 else 'else if', n, cpp(g['label']), g['min'], g['max']) for n, g in enumerate(group_set)))
        text = text[:group.start()] + 'if (groupSet == 0) {\n'+blocks[0]+'\n} else {\n'+blocks[1]+'\n}' + text[end:]
        changes.append(dict(label='configured_group_sets', count=1))
    sub('measured_SN', r'double\s+(sn(?:_global)?)\s*=\s*[^;]+;',
        lambda m: 'double %s = %s;' % (m.group(1), num(purity/(1.0-purity))))
    for variable, key in [('mean','mean'),('meanLambda','mean'),('sigma','sigma'),('sigmaLambda','sigma'),('nsigma','nsigma')]:
        sub('maker_'+variable, r'(double\s+'+variable+r'\s*=\s*)[^;]+;', lambda m,k=key: m.group(1)+num(window[k])+';')
    sub('maker_signal_width', r'3\.0\s*\*\s*sigmaLambda', num(window['nsigma'])+'*sigmaLambda')
    sub('maker_SB_outer', r'4\.0\s*\*\s*nsigma', num(window['sideband_outer_factor'])+'*nsigma')
    sub('maker_SB_scale', r'1\.0\s*/\s*3\.0', num(window['sideband_scale']))
    rebin = settings['rebin_group' if '_sum' in name else 'rebin_fine']
    sub('configured_rebin', r'(else if\s*\(sp\s*==\s*"4He"\)\s*nRebin\s*=\s*)\d+;', lambda m: m.group(1)+str(rebin)+';')
    sub('configured_standalone_rebin', r'int nRebin\s*=\s*10;', 'int nRebin = %d;' % rebin)
    if name != 'Draw_LambdaInvMass.C':
        plot = settings['plot_x_range']
        sub('configured_kstar_display',
            r'GetXaxis\(\)->SetRangeUser\(\s*0\.0\s*,\s*0\.(?:6|15|40)\s*\)',
            'GetXaxis()->SetRangeUser(%s,%s)' % tuple(num(x) for x in plot))
        sub('configured_unity_line', r'new TLine\(0\.0,\s*1\.0,\s*0\.6,\s*1\.0\)',
            'new TLine(%s,1.0,%s,1.0)' % tuple(num(x) for x in plot))
    for old, value in [('0.40001',settings['normalization'][0]+1e-5),('0.59999',settings['normalization'][1]-1e-5)]:
        sub('configured_norm', re.escape(old), num(value))
    sub('configured_sliding_width', r'int W\s*=\s*2;', 'int W = %d;' % settings['sliding_half_width'])
    if '_sum' not in name and '_method' in name:
        method = re.search(r'method[35]', name).group(0)
        sub('fine_cent_pdf_name', r'TString outFile\s*=\s*outDir\s*\+\s*"[^"]+";',
            'TString outFile = outDir + "/kstar_cent_CFsub_%s_4He.pdf";' % method, 1)
    purity_match = re.search(r'void\s+getPurity\s*\([^)]*\)\s*\{', text)
    if purity_match:
        end = close_brace(text, purity_match.end())
        wrapper = '''void getPurity(TH1D* hY, double minM, double maxM, double &P, double &err_P,
  double &N_sig, double &err_N_sig, double &N_tot, double &err_N_tot) {
  KFReport::Purity(hY,minM,maxM,P,err_P,N_sig,err_N_sig,N_tot,err_N_tot);
}'''
        text = text[:purity_match.start()]+wrapper+text[end:]
        changes.append(dict(label='audited_mass_fit', count=1))
        sub('safe_CF_division', r'if\s*\(nm\s*>\s*0\)',
            'if (nt > 0 && nm > 0 && TMath::Finite(nt) && TMath::Finite(nm) && TMath::Finite(et) && TMath::Finite(em))', 1)
        sub('safe_CF_error', r'double err_cf\s*=\s*cf\s*\*\s*TMath::Sqrt\([^;]+;',
            'double err_cf = TMath::Sqrt(TMath::Power(et/nm,2)+TMath::Power(nt*em/(nm*nm),2));', 1)
        sub('mask_missing_points', r'(hCF|hPTrue|hPMix)->Draw\("(E1(?: SAME)?)"\);',
            lambda m: 'KFReport::DrawMeasured(%s,"%s");' % m.groups())
    ll = settings['ll_fit']
    for variable, key in [('a0','fixed_a0'),('reff','fixed_reff'),('fit_max',None)]:
        sub('configured_LL_'+variable, r'(double\s+'+variable+r'\s*=\s*)[0-9.]+;',
            lambda m,k=key: m.group(1)+num(ll[k] if k else ll['range'][1])+';')
    sub('configured_LL_R_bounds', r'SetParLimits\(0,\s*0\.1,\s*6\.\)',
        'SetParLimits(0,%s,%s)' % tuple(num(x) for x in ll['radius_limits']))
    sub('configured_LL_initial_fixed', r'SetParameters\(3\.0,\s*1\.0\)',
        'SetParameters(%s,%s)' % (num(ll['radius_initial']),num(ll['normalization_initial'])))
    sub('configured_LL_initial_free', r'SetParameters\(3\.0,\s*4\.01,\s*1\.84,\s*1\.0\)',
        'SetParameters(%s,%s,%s,%s)' % tuple(num(ll[k]) for k in ['radius_initial','fixed_a0','fixed_reff','normalization_initial']))
    sub('configured_LL_lower_bound', r'(fitLL_(?:fixed|free),\s*)0\.0,', lambda m: m.group(1)+num(ll['range'][0])+',')
    sub('record_LL_status', r'(hCF\w*)->Fit\((fit[12]),\s*"(R\+?)"\);',
        lambda m: 'KFReport::FitLL(%s,%s,"%s");' % m.groups())
    if name.startswith('Fit_'):
        # Replace, rather than overpaint, the long legacy parameter text.
        # Overpainting left trailing digits outside the sub-pad text boxes.
        sub('compact_LL_parameter_panel',
            r'TLatex tex;.*?tex\.DrawLatex\(0\.4,\s*0\.38,.*?;',
            'KFReport::DrawLLWarnings(fit1,fit2);', 1, re.S)
        sub('compact_LL_legend', r'new TLegend\(0\.4,\s*0\.7,\s*0\.9,\s*0\.9\)',
            'new TLegend(0.52,0.76,0.88,0.89)', 1)
        sub('short_LL_fixed_label', r'"Fit \(a0, reff fixed\)"', '"LL: fixed"', 1)
        sub('short_LL_free_label', r'"Fit \(all free\)"', '"LL: free"', 1)
    function = Path(name).stem
    audit = outdir/'logs'/(function+'_fit_audit.tsv')
    mf = settings['mass_fit']
    values = mf['range']+mf['mean_limits']+mf['sigma_limits']+[mf['sigma_initial'],mf['minimum_entries'],mf['minimum_covariance_status']]
    init = '\n const double massFitSettings[9] = {'+','.join(num(x) for x in values)+'};\n KFReport::Init(%s,%s,%s,massFitSettings);' % (cpp(audit),num(window['mean']),num(window['sigma']))
    sub('audit_init', r'(void\s+'+function+r'\s*\(\s*\)\s*\{)', lambda m: m.group(1)+init, 1)
    text = '// Generated adaptation; legacy original preserved.\n#include "KFSummarySupport.h"\n'+text
    if 'anaLambdaNuclearId/' in text or 'nSpecies = 4' in text or 'fIn->Get("nuclearid")' in text:
        raise ValueError(name+': legacy input/species survived adaptation')
    return text, changes, audit

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--energy', choices=['13p5','3p85'], required=True)
    for arg in ['input-root','output-dir','maker-yaml','settings']:
        parser.add_argument('--'+arg, type=Path, required=True)
    parser.add_argument('--purity', type=float, required=True)
    parser.add_argument('--refresh-generated', action='store_true',
        help='Archive and refresh only previously manifested generated files; never replace untracked changes.')
    args = parser.parse_args()
    if not math.isfinite(args.purity) or not 0 < args.purity < 1:
        parser.error('purity must be finite and strictly between zero and one')
    input_root = args.input_root.resolve(strict=True)
    if not input_root.is_file() or input_root.suffix != '.root':
        parser.error('input-root must be an existing ROOT file')
    outdir = args.output_dir.resolve()
    if outdir == input_root.parent or outdir.name != 'pdf':
        parser.error('output-dir must be a dedicated pdf/ subdirectory')
    maker = args.maker_yaml.resolve(strict=True)
    config = yaml.safe_load(maker.read_text())
    window = {label:float(config[key]) for label,key in [('mean','lambdaSignalMean'),('sigma','lambdaSignalSigma'),
        ('nsigma','lambdaSignalNSigma'),('sideband_outer_factor','lambdaSidebandOuterFactor')]}
    if not all(math.isfinite(v) for v in window.values()) or window['sigma'] <= 0 or window['nsigma'] <= 0 or window['sideband_outer_factor'] <= 1:
        parser.error('invalid required maker mass-window values')
    half = window['nsigma']*window['sigma']
    window.update(signal_min=window['mean']-half, signal_max=window['mean']+half,
        sideband_scale=1.0/(window['sideband_outer_factor']-1.0))
    for key,value in [('channel_1_signalMin',window['signal_min']),('channel_1_signalMax',window['signal_max'])]:
        if abs(float(config[key])-value) > 1e-9:
            parser.error('maker signal channel inconsistent with mean/sigma/nsigma: '+key)
    settings_path = args.settings.resolve(strict=True)
    settings = yaml.safe_load(settings_path.read_text())
    if [len(x) for x in settings['group_sets']] != [2,3]:
        parser.error('legacy layout requires two groups, then three groups')
    if not 0 <= settings['normalization'][0] < settings['normalization'][1]:
        parser.error('invalid CF normalization range')
    if (len(settings['plot_x_range']) != 2 or
            not all(math.isfinite(float(x)) for x in settings['plot_x_range']) or
            not 0 <= settings['plot_x_range'][0] < settings['plot_x_range'][1]):
        parser.error('invalid kstar display range')
    for key in ['rebin_group','rebin_fine','sliding_half_width']:
        if not isinstance(settings[key],int) or settings[key] < (0 if key=='sliding_half_width' else 1):
            parser.error('invalid positive integer setting: '+key)
    source_dir = ROOT/'tools'
    if args.energy == '3p85': source_dir /= 'auau3p85GeV'
    outputs, entries = [], []
    for name,pdf,pages in TASKS:
        source = source_dir/name
        original = source.read_bytes()
        text,changes,audit = adapt(original.decode(),name,input_root,outdir,args.purity,window,settings)
        data = text.encode(); target = outdir/'macros'/name
        entries.append(dict(name=name,function=source.stem,source=str(source),path=str(target),generated=str(target),
            source_sha256=sha(original),generated_sha256=sha(data),changes=changes,
            expected_pdfs=[pdf],pdfs=[dict(name=pdf,pages=pages)],audit_file=str(audit)))
        outputs.append((target,data))
    header = Path(__file__).resolve().with_name('KFSummarySupport.h'); header_data = header.read_bytes()
    outputs.append((outdir/'macros'/header.name,header_data))
    manifest = dict(schema_version=1,energy=args.energy,input_root=str(input_root),input_size=input_root.stat().st_size,
        input_mtime_ns=input_root.stat().st_mtime_ns,output_dir=str(outdir),maker_yaml=str(maker),maker_yaml_sha256=sha(maker.read_bytes()),
        settings_path=str(settings_path),settings_sha256=sha(settings_path.read_bytes()),settings=settings,
        purity=args.purity,mass_window=window,macros=entries,support=dict(source=str(header),sha256=sha(header_data)),
        caveats=['LL fits are exploratory: consult status/covariance/boundaries.',
                 'Global Lambda purity is not pair purity resolved in kstar and centrality.',
                 'Legacy CFsub omits global-purity and normalization uncertainty.',
                 'Sliding-window and normalization cross-bin covariance are not propagated.',
                 'M3/M5 whole-bin mass projection approximates continuous Maker gate.'])
    manifest_path = outdir/'macro_manifest.json'
    previous_bytes = manifest_path.read_bytes() if manifest_path.exists() else None
    outputs.append((manifest_path,stable_manifest_bytes(manifest,previous_bytes)))
    differing = [(path,data) for path,data in outputs if path.exists() and path.read_bytes() != data]
    archive = None
    if differing:
        if not args.refresh_generated:
            raise FileExistsError('Refusing to replace different generated file without --refresh-generated: '+str(differing[0][0]))
        old_manifest_path = outdir/'macro_manifest.json'
        if not old_manifest_path.is_file():
            raise FileExistsError('Cannot refresh files without the prior ownership manifest')
        old_manifest = json.loads(old_manifest_path.read_text())
        if old_manifest['input_root'] != str(input_root) or old_manifest['energy'] != args.energy:
            raise ValueError('Refresh may not switch the input dataset in an existing report directory')
        known = {Path(m['generated']):m['generated_sha256'] for m in old_manifest['macros']}
        known[outdir/'macros'/header.name] = old_manifest['support']['sha256']
        for path,data in differing:
            if path == old_manifest_path:
                continue
            if path not in known or sha(path.read_bytes()) != known[path]:
                raise FileExistsError('Refusing to overwrite an unowned or locally edited generated file: '+str(path))
        archive = outdir/'macros'/('archive_'+datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%S%fZ'))
        archive.mkdir(parents=True, exist_ok=False)
        # Preserve all changed originals, including their ownership manifest,
        # before replacing any file. The parent runner controls active processes.
        for path,data in differing:
            (archive/path.name).write_bytes(path.read_bytes())
    (outdir/'macros').mkdir(parents=True,exist_ok=True); (outdir/'logs').mkdir(parents=True,exist_ok=True)
    for path,data in outputs:
        if not path.exists():
            path.write_bytes(data)
        elif path.read_bytes() != data:
            temporary = path.with_name(path.name+'.refresh_tmp')
            with temporary.open('xb') as stream:
                stream.write(data)
            temporary.replace(path)
    print(json.dumps(dict(manifest=str(outdir/'macro_manifest.json'),macros=len(entries),
                          archived_previous=str(archive) if archive else None),indent=2))

if __name__ == '__main__':
    main()
