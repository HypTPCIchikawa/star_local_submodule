#!/usr/bin/env python3
"""Four mT slices with fixed-P CFsub and explicitly shared inclusive Method-5 P(k*)."""
import argparse
import json
from pathlib import Path
import shutil
import sys

import yaml

from compare_mass_halves import (HERE, PROJECT, command, digest, fingerprint, save,
                                 selected_bins, purity_proxies)
from mass_half_statistics import ErrorModel, normalized_cf, purity_corrected_cf
from mt_binning import choose_slices
from report_mass_halves import escape, pdf_pages


def sumq(model, values):
    return sum(values, model.constant(0))


def rebin(model, values, factor):
    if len(values) % factor:
        raise ValueError('Rebin factor must divide the saved k* grid.')
    return [sumq(model,values[i:i+factor]) for i in range(0,len(values),factor)]


def shared_purities(model, source_group, source_settings):
    output = {}
    for sample in ('SE','ME'):
        rows = source_group['method5Base'][sample]
        total = [model.variable(('purity_counts',sample,i),sum(r['counts']),sum(r['error2']))
                 for i,r in enumerate(rows)]
        factors = []
        for i,row in enumerate(source_group['methods']['method5'][sample]):
            if not row['full']['valid']:
                factors.append(None)
                continue
            b = model.correlated_block(('purity_background',sample,i),row['background'],
                                       row['backgroundCovariance'])
            first,last = row['windowOriginalFirstBin']-1,row['windowOriginalLastBin']
            n = sumq(model,total[first:last])
            factors.append(1-sumq(model,b)/n if n.value>0 else None)
        output[sample] = factors
        if len(factors)*source_settings['purityRebin'] != len(total):
            raise ValueError('Purity grid does not cover original counts.')
    return output


def calculate(data, selection, settings, base, purity_source):
    group = next(g for g in purity_source['groups'] if
                 g['centMin']==settings['cent9_min'] and g['centMax']==settings['cent9_max'])
    native_edges = data['channels']['SE']['signal']['kstarEdges']
    if native_edges != group['originalKstarEdges']:
        raise ValueError('mT / mass-fit original k* grids differ.')
    factor = base['rebin_group']
    edges = native_edges[::factor]
    norm = selected_bins(edges,base['normalization'])
    display = selected_bins(edges,base['plot_x_range'])
    model = ErrorModel()
    pconst = purity_proxies(purity_source)['full']
    shared = shared_purities(model,group,purity_source['settings'])
    coarse_factor = purity_source['settings']['purityRebin']
    results = dict(normalization=base['normalization'],plot_x_range=base['plot_x_range'],
                   kstar_edges=edges,constant_purity=pconst,selection=selection,slices=[],
                   method5='mT-integrated P(k*) shared across slices; not mT-resolved fits')
    lines = ['method\tslice\tlabel\tlow\thigh\tvalue\terror']
    all_quantities = {m: [] for m in ('constant_purity','method5_common_purity')}
    for index,part in enumerate(selection['slices']):
        series = {}
        for sample in ('SE','ME'):
            series[sample] = {}
            for region in ('signal','leftSB','rightSB'):
                channel = data['channels'][sample][region]
                vals = [sum(r[part['first']:part['stop']]) for r in channel['counts']]
                err2 = [sum(r[part['first']:part['stop']]) for r in channel['error2']]
                series[sample][region] = model.series((index,sample,region),vals,err2)
            series[sample]['sideband'] = [a+b for a,b in zip(series[sample]['leftSB'],series[sample]['rightSB'])]
        signal = normalized_cf(rebin(model,series['SE']['signal'],factor),
                               rebin(model,series['ME']['signal'],factor),norm)
        background = normalized_cf(rebin(model,series['SE']['sideband'],factor),
                                   rebin(model,series['ME']['sideband'],factor),norm)
        const = [purity_corrected_cf(a,b,pconst['value']) if a is not None and b is not None else None
                 for a,b in zip(signal['bins'],background['bins'])]
        weighted, available = {}, {}
        for sample in ('SE','ME'):
            native = [q*shared[sample][i//coarse_factor] if shared[sample][i//coarse_factor] is not None else None
                      for i,q in enumerate(series[sample]['signal'])]
            available[sample] = [all(q is not None for q in native[i:i+factor])
                                 for i in range(0,len(native),factor)]
            weighted[sample] = rebin(model,[q if q is not None else model.constant(0) for q in native],factor)
        if any(not available[s][i] for s in ('SE','ME') for i in norm):
            method5 = [None]*(len(edges)-1)
        else:
            candidate = normalized_cf(weighted['SE'],weighted['ME'],norm)
            method5 = [q if available['SE'][i] and available['ME'][i] else None
                       for i,q in enumerate(candidate['bins'])]
        record = dict(low=part['low'],high=part['high'],methods={})
        label = '%.2f-%.2f' % (part['low'],part['high'])
        for method,values in [('constant_purity',const),('method5_common_purity',method5)]:
            all_quantities[method].append(values)
            record['methods'][method] = dict(points=[model.estimate(q) if q is not None else None for q in values],
                valid_plot_bins=sum(values[i] is not None for i in display),total_plot_bins=len(display))
            for i,q in enumerate(values):
                estimate = model.estimate(q) if q is not None else None
                lines.append('\t'.join([method,str(index),label,format(edges[i],'.17g'),format(edges[i+1],'.17g'),
                    format(estimate['value'],'.17g') if estimate else 'nan',
                    format(estimate['error'],'.17g') if estimate else 'nan']))
        results['slices'].append(record)
    # Small per-k* slice covariance blocks retain the same shared purity variables.
    results['cross_slice_covariance'] = {m:[
        [[model.covariance(rows[a][i],rows[b][i]) if rows[a][i] is not None and rows[b][i] is not None else None
          for b in range(len(rows))] for a in range(len(rows))] for i in range(len(edges)-1)]
        for m,rows in all_quantities.items()}
    results['shared_purities'] = {s:[model.estimate(p) if p is not None else None for p in shared[s]] for s in shared}
    return results,'\n'.join(lines)+'\n'


def write_plot_inputs(out,data,selection):
    lines = ['sample\tkbin\tmtbin\tklow\tkhigh\tmtlow\tmthigh\tcount']
    for sample in ('SE','ME'):
        h = data['channels'][sample]['signal']
        for x,row in enumerate(h['counts']):
            for y,count in enumerate(row):
                lines.append('\t'.join(map(str,[sample,x,y,h['kstarEdges'][x],h['kstarEdges'][x+1],
                                                h['mtEdges'][y],h['mtEdges'][y+1],count])))
    (out/'mt_grid.tsv').write_text('\n'.join(lines)+'\n')
    boundaries = [s['low'] for s in selection['slices']]+[selection['slices'][-1]['high']]
    (out/'mt_edges.txt').write_text(' '.join(map(str,boundaries))+'\n')


def report(out,context,data,result):
    base = context['base']
    p = result['constant_purity']['value']
    slices = result['selection']['slices']
    k = data['channels']['SE']['signal']['kstarEdges']
    low_bins = selected_bins(k,context['settings']['split_kstar'])
    norm_bins = selected_bins(k,base['normalization'])
    table = []
    for i,part in enumerate(slices):
        counts = []
        for sample,region,xbins in [('SE','signal',low_bins),('ME','signal',low_bins),
                ('SE','signal',norm_bins),('ME','signal',norm_bins),('SE','leftSB',low_bins),('SE','rightSB',low_bins)]:
            h = data['channels'][sample][region]
            counts.append(sum(sum(h['counts'][j][part['first']:part['stop']]) for j in xbins))
        table.append([i+1,'%.2f–%.2f' % (part['low'],part['high'])]+[int(v) for v in counts])
    notes = [
        'Method5は0–10%・m_T積分のP_SE(k*)とP_ME(k*)を4区間で共通使用する参考比較です。m_T別に純度をfitした結果ではありません。',
        '保存済みROOTにk*×m_T×Λ質量の同時分布がないため、純度のm_T依存は検証できません。',
        'constant purityは前回と同じinclusive質量分布から求めた一定値P=%.8fを全m_T・k*に使用します。各m_Tのsideband CFを個別に規格化して補正します。' % p,
        '両方式とも各m_T区間のSE・MEから規格化します。完成したCF同士を加算して中央度を統合しません。',
        'Method5は純度50 MeV/c刻み・5bin（内部250 MeV/c）窓を維持。元5 MeV/cのカウントに共通Pを掛けて20 MeV/cへ合算します。',
        '同一純度Pの再利用と規格化に由来する相関は計算に含めます。異なるfit間の共分散、fit背景とカウントの相関、m_Tカウントとm_T積分純度の相関、イベント・混合ペア再利用の相関は未評価です。',
        'constant purityのPは固定値として誤差計算します。誤差は近似で、最小k*でSEが0または少数の点はPoisson区間ではありません。',
        'Method5純度は質量bin全体、m_T分布は候補への厳密な±3σゲートから作られます。この境界差も残るため、m_T別独立fitの代替や精密な物理結論にはしません。',
        '対象は2021 P24iyの26,956 / 26,994 subjob部分統合。元ROOT、選択cut、mixingは変更していません。']
    displayed = selected_bins(result['kstar_edges'],base['plot_x_range'])
    minimum_se = [min(part['channels']['SE']['signal']['rebinned']['counts'][i] for i in displayed)
                  for part in slices]
    notes.append('表示内の1 binあたり最少signal SEペア数はm_T区間順に%sです。全区間の総数が多くても、最小k*付近は低統計です。' %
                 ' / '.join(str(int(n)) for n in minimum_se))
    if result['selection']['warnings']:
        notes.append('flow監査のwarningは比較・規格化の範囲外のカウントに関するものです。今回のm_T underflow/overflowは0で、k*≥1 GeV/cのoverflowはCFに使用していません。')
    readme = ['# Λ–⁴He CFのm_T 4分割（3.85 GeV、Centrality 0–10%）','',
        '[まとめPDF](cf_mt_cent0_10_constant_vs_method5.pdf)','',
        '- Centrality: cent9=7・8をSE/MEそれぞれ加算。signalは±3σ。',
        '- CF: 20 MeV/c幅、表示0–0.3 GeV/c、規格化%.2f–%.2f GeV/c。' % tuple(base['normalization']),
        '- m_T: sqrt((|pTΛ+pTA|/2)² + ((mΛ+mA)/2)²)。',
        '- 境界はk*＝0–0.3 GeV/cの信号SE数の累積25/50/75%に近い保存bin境界を選択。CF形状は最適化していません。','',
        '| Slice | m_T [GeV/c²] | SE (0–0.3) | ME (0–0.3) | SE norm | ME norm | SE leftSB (0–0.3) | SE rightSB (0–0.3) |',
        '|---|---|---:|---:|---:|---:|---:|---:|']
    readme += ['| '+' | '.join(map(str,row))+' |' for row in table]
    readme += ['', '各区間は[下限,上限)で、同じ境界を両方法・SE/ME・sidebandに使用します。','',
        '## 注意点','']+['- '+s for s in notes]
    readme += ['', '## ファイル','',
        '- `mt_distribution.pdf`: 選択した境界とSE/MEのk*–m_T分布。',
        '- `mt_cf_comparison.pdf`: m_T区間ごとの2方式の比較。',
        '- `mt_cf_overlays.pdf`: 方式ごとの4区間の重ね描き。',
        '- `mt_selection.json`: 区間選択、規格化・sideband統計、ゼロbin・flow検査。',
        '- `mt_comparison.json`: CF・近似誤差・共通純度・同一k*のm_T区間間共分散。',
        '- `mt_extract.json`: 入力ROOTのキー一覧と抽出ヒスト、射影検査。','']
    (out/'README.md').write_text('\n'.join(readme))
    rows = [' & '.join(escape(str(v)) for v in row[:6])+r' \\' for row in table]
    tex = r'''\documentclass[10pt,a4paper]{article}
\usepackage[margin=15mm]{geometry}
\usepackage{amsmath,booktabs,pdfpages,hyperref}
\pagestyle{plain}
\begin{document}
\begin{center}\Large Au+Au 3.85 GeV: $\Lambda$--$^4$He CF in four $m_T$ slices\\
\large Centrality 0--10\%, KFParticle\end{center}
\noindent '''+escape('CF: 20 MeV/c; display 0-0.3 GeV/c; normalization %.2f-%.2f GeV/c.' % tuple(base['normalization']))+r'''

\medskip
\noindent The boundaries are selected at stored 0.05 GeV/$c^2$ edges,
near SE signal quartiles in $0\leq k^*<0.3$ GeV/$c$.
The correlation shapes are not used to optimize the boundaries.
The same slices are used for both methods, SE, ME and sidebands.

\begin{center}\small\begin{tabular}{rrrrrr}\toprule
Slice & $m_T$ [GeV/$c^2$] & SE (0--0.3) & ME (0--0.3) & SE norm & ME norm\\\midrule
'''+ '\n'.join(rows)+r'''
\bottomrule\end{tabular}\end{center}

\noindent\textbf{Method 5 is an approximation, not an $m_T$-resolved mass fit.}
The available histograms have no joint $(k^*,m_T,M_\Lambda)$ distribution.
The same $m_T$-integrated, centrality-0--10\% $P_{SE}(k^*)$ and $P_{ME}(k^*)$
are applied to all four slices. Purity grid: 50 MeV/$c$;
five-bin sliding window: 250 MeV/$c$ in the interior.
Each original 5 MeV/$c$ count is weighted before summing to 20 MeV/$c$.

\medskip
\noindent Constant purity: $P='''+format(p,'.8f')+r'''$ (the same inclusive proxy as before),
with $C_{corr}=[C_{sig}-(1-P)C_{SB}]/P$. Sidebands and normalization
are evaluated separately in each $m_T$ slice.

\medskip
\noindent Approximate errors include shared purity reuse and normalization.
Covariances between mass-fit backgrounds and observed counts, distinct overlapping
fits, $m_T$ counts and inclusive purity, and reused event/mixed pairs are not available.
Constant-$P$ errors are conditional on fixed $P$.
Sparse/zero-SE bins do not have Poisson confidence intervals.
The smallest SE signal bins contain only a few pairs despite the large slice totals.
Exact candidate mass gates and whole mass-bin fit integrals differ slightly.

\medskip
\noindent 2021 P24iy only; 26,956/26,994 subjobs in this partial merge.
The input ROOT and all production cuts are unchanged.
Full provenance, statistics and caveats are in the adjacent README and JSON files.
\clearpage
'''
    figures = ['mt_distribution.pdf','mt_cf_comparison.pdf','mt_cf_overlays.pdf']
    for name in figures:
        pdf_pages(out/name)
        tex += r'\includepdf[pages=-]{'+name+'}\n'
    tex += '\\end{document}\n'
    (out/'cf_mt_cent0_10_constant_vs_method5.tex').write_text(tex)
    import subprocess
    for i in range(2):
        with (out/'logs'/('report_latex_%d.log'%i)).open('w') as log:
            subprocess.run(['pdflatex','-interaction=nonstopmode','-halt-on-error',
                            'cf_mt_cent0_10_constant_vs_method5.tex'],cwd=str(out),stdout=log,stderr=subprocess.STDOUT,check=True)
    return pdf_pages(out/'cf_mt_cent0_10_constant_vs_method5.pdf')


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--mainconf',required=True,type=Path)
    parser.add_argument('--input-root',required=True,type=Path)
    parser.add_argument('--purity-study',required=True,type=Path)
    parser.add_argument('--shared-method5-approved',action='store_true',required=True)
    parser.add_argument('--settings',type=Path,default=HERE/'mt_settings.yaml')
    parser.add_argument('--output-dir',type=Path)
    parser.add_argument('--resume',action='store_true')
    args = parser.parse_args()
    rootfile = args.input_root.resolve(strict=True)
    mainconf = args.mainconf.resolve(strict=True)
    settings_path = args.settings.resolve(strict=True)
    settings = yaml.safe_load(settings_path.read_text())
    base_path = (settings_path.parent/settings['summary_settings']).resolve(strict=True)
    base = yaml.safe_load(base_path.read_text())
    purity_study = args.purity_study.resolve(strict=True)
    old_context = json.loads((purity_study/'context.json').read_text())
    if old_context['input']!=fingerprint(rootfile) or old_context['energy']!='3p85':
        raise ValueError('Purity source must belong to the same unchanged 3.85 GeV ROOT.')
    if old_context['base']!=base or old_context['mainconf']!=str(mainconf):
        raise ValueError('Purity source settings / mainconf differ.')
    purity_path = purity_study/'mass_half_extract.json'
    purity_source = json.loads(purity_path.read_text())
    out = (args.output_dir or rootfile.parent/'pdf'/'mt_cent0_10').resolve()
    out.mkdir(parents=True,exist_ok=True)
    context = dict(input=fingerprint(rootfile),mainconf=str(mainconf),base=base,settings=settings,
        settings_sha256=digest(settings_path),base_sha256=digest(base_path),
        purity_source=str(purity_path),purity_source_sha256=digest(purity_path),
        shared_method5_approved=True,energy='3p85')
    if (out/'context.json').exists():
        if not args.resume or json.loads((out/'context.json').read_text())!=context:
            raise ValueError('Existing context differs; do not overwrite a different study.')
    else:
        allowed={'mt_extract.json','macros','logs'}
        if any(p.name not in allowed for p in out.iterdir()):
            raise ValueError('Unowned existing output directory.')
        save(out/'context.json',context)
    (out/'logs').mkdir(exist_ok=True)
    (out/'macros').mkdir(exist_ok=True)
    root = ['bash',str(PROJECT/'script/singularity_root_batch.sh'),str(mainconf),'root4star','-l','-b','-q']
    if not (out/'mt_extract.json').exists():
        copy = out/'macros'/'ExtractMtCF.C'
        shutil.copy2(str(HERE/copy.name),str(copy))
        call = str(copy)+'+('+','.join([json.dumps(str(rootfile)),json.dumps(str(out)),
            str(settings['cent9_min']),str(settings['cent9_max'])])+')'
        command(root+[call],out/'logs/extract.log')
    data = json.loads((out/'mt_extract.json').read_text())
    if data['inputRoot']!=str(rootfile) or (data['centMin'],data['centMax'])!=(settings['cent9_min'],settings['cent9_max']):
        raise ValueError('mT extraction input / centrality mismatch.')
    selection = choose_slices(data,settings,base)
    result,rows = calculate(data,selection,settings,base,purity_source)
    if (out/'mt_cf_overlays.pdf').exists():
        if (out/'mt_cf_points.tsv').read_text()!=rows or digest(out/'macros'/'DrawMtCF.C')!=digest(HERE/'DrawMtCF.C'):
            raise ValueError('Existing plots differ from new data / renderer; preserve them and use a new output directory.')
    save(out/'mt_selection.json',selection)
    save(out/'mt_comparison.json',result)
    (out/'mt_cf_points.tsv').write_text(rows)
    write_plot_inputs(out,data,selection)
    if not (out/'mt_cf_overlays.pdf').exists():
        copy = out/'macros'/'DrawMtCF.C'
        shutil.copy2(str(HERE/copy.name),str(copy))
        call = str(copy)+'+('+','.join([json.dumps(str(out/'mt_cf_points.tsv')),json.dumps(str(out))]+[
            repr(v) for v in base['normalization']+[base['plot_x_range'][1]]]+[
            json.dumps(str(out/'mt_grid.tsv')),json.dumps(str(out/'mt_edges.txt'))])+')'
        command(root+[call],out/'logs/plot.log')
    pages = report(out,context,data,result)
    if fingerprint(rootfile)!=context['input']:
        raise ValueError('Input ROOT changed.')
    sources = [Path(__file__),HERE/'mt_binning.py',HERE/'ExtractMtCF.C',HERE/'DrawMtCF.C',HERE/'mass_half_statistics.py']
    save(out/'completed.json',dict(phase='complete',input_unchanged=True,pages=pages,
        sources={str(p):digest(p) for p in sources},
        products={p.name:digest(p) for p in out.iterdir() if p.suffix in ('.json','.pdf','.tsv','.txt','.md') and p.name!='completed.json'}))
    print('COMPLETE',out,flush=True)


if __name__=='__main__':
    main()
