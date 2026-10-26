#!/usr/bin/env python3
"""Build a provenance-aware report from already computed mass-half CFs."""
import argparse
import datetime
import hashlib
import json
import math
import os
from pathlib import Path
import shutil
import statistics
import subprocess
import sys

OWNER = 'tools/kfparticle_summary/report_mass_halves.py'
METHODS = ('raw', 'sbsub', 'cfsub', 'method3', 'method5')
REGIONS = ('full', 'low', 'high')
NAMES = dict(raw='Raw', sbsub='SBsub', cfsub='CFsub', method3='Method 3', method5='Method 5')
SAMPLES = {
    '13p5': dict(beam='13.5', year=2020, included=542, expected=545,
                 jobid='1213058CE62BCF113B2E3C75ACBD7053'),
    '3p85': dict(beam='3.85', year=2021, included=26956, expected=26994,
                 jobid='0BEFC64E62F19BBA534A129CABC9B106')}
JP_METHODS = {
    'raw': 'signal SE/MEを独立規格化。背景補正なし。',
    'sbsub': 'SE・MEの各分布から左右SBの合計を引いてから規格化。幅係数はfull=1/3、各half=1/6。平坦背景を仮定。',
    'cfsub': '各窓のCFを (Cwindow − (1−Pwindow) CSB) / Pwindow で補正。Pは現在のinclusive質量分布から窓別に推定。',
    'method3': '各k* binのgaus+pol2 fitでN−Bを求める。low/highの背景積分は同一fitから取り、fullでは合計。',
    'method5': '粗いk* binごとの滑動窓fitでPを推定。元の細粒度カウントへ所属するcoarse binのPを掛け、CF用binへ合算する。窓内総数を中心binの収量にしない。'}
JP_CAVEATS = [
    'fullとhalfの共通カウント、共通sideband、規格化積分による共分散を含む一次誤差伝播です。独立なヒスト同士としてratio誤差を計算していません。',
    'M3/M5では同一fitの低側・高側背景積分の共分散を含めます。M5の共有入力カウントと、同じcoarse binのPを複数のCF binで使う共分散も計算します。',
    '同一イベント・トラックの共有、混合ペア再利用による相関はヒストだけから復元できず未反映です。',
    '観測カウントNとfit背景Bの相関、異なるcoarse binの重なる滑動窓fit間の共分散、CFsubのPの誤差・背景モデル系統は未反映です。完全な誤差として扱わないでください。',
    '必要な規格化binのfitが無効なら、そのCF全体をNAにします。積分比較は対象区間の全binが必要で、使えるbinだけの別区間へ黙って変更しません。失敗fitを0には置き換えません。',
    'M3ではfullのN−Bは両halfの和です。M5のfullはfull専用の窓純度を中心binへ掛けるため、low+highと必ずしも同一ではありません。',
    'CF自体を足してfullにはできません。各窓のSE/MEをそれぞれ独立規格化しています。',
    'gaus+pol2の数値収束はモデル妥当性を保証しません。大きなχ²/ndfや小さいfit確率には注意が必要で、高精度モデルや物理効果の確定を主張しません。',
    '中央度のgroupsは互いに重なります。独立実験としてまとめて有意度を計算できません。',
    'TH2Fの丸めによるfull−low−high残差をcomparison.jsonのclosureに保持しています。元fullを黙ってlow+highに置き換えていません。']


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def save(path, data):
    path.write_text(json.dumps(data, ensure_ascii=False, indent=2, allow_nan=False) + '\n')


def finite(value):
    return isinstance(value, (int, float)) and math.isfinite(value)


def number(value, digits=5):
    return format(value, '.%dg' % digits) if finite(value) else 'NA'


def escape(value):
    chars = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$',
             '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}',
             '~': r'\textasciitilde{}', '^': r'\textasciicircum{}'}
    return ''.join(chars.get(c, c) for c in str(value))


def measure(q, tex=False, percent=False):
    if not q or not finite(q.get('value')) or not finite(q.get('error')):
        return 'NA'
    v, e = q['value'], q['error']
    if percent:
        v, e = 100*(v-1), 100*e
    if tex:
        return '$%s\\pm%s$' % (number(v, 4), number(e, 3))
    return '%s ± %s%s' % (number(v), number(e, 3), '%' if percent else '')


def ratio(metrics, region):
    return (metrics[region].get('vsFull') or {}).get('ratio')


def valid_text(method):
    return '/'.join(str(method['valid_plot_bins'][r]) for r in REGIONS) + ' of ' + str(method['total_plot_bins'])


def metric_cells(method, tex=False):
    m = method['integrated']
    return [measure(m[r]['cf'], tex) for r in REGIONS] + [
        measure(ratio(m, r), tex) for r in ('low', 'high')] + [valid_text(method)]


def reasons(value):
    if isinstance(value, str):
        return [value]
    if isinstance(value, dict):
        return [text for v in value.values() for text in reasons(v)]
    if isinstance(value, list):
        return [text for v in value for text in reasons(v)]
    return []


def invalid_reason(method, region):
    texts = list(dict.fromkeys(reasons(method['invalid'].get(region))))
    if not texts and method['integrated'][region]['cf'] is None:
        texts = ['Invalid fit / denominator in comparison or normalization interval']
    return '; '.join(texts)


def pdf_pages(path):
    with path.open('rb') as stream:
        if stream.read(5) != b'%PDF-':
            raise ValueError('Invalid PDF: ' + str(path))
    gs = shutil.which('gs')
    if not gs:
        raise RuntimeError('Ghostscript is required to validate PDF page counts.')
    result = subprocess.run([gs, '-q', '-dNODISPLAY', '-dSAFER',
                             '--permit-file-read=' + str(path), '-sPDFFile=' + str(path),
                             '-c', 'PDFFile (r) file runpdfbegin pdfpagecount = quit'],
                            stdout=subprocess.PIPE, stderr=subprocess.STDOUT,
                            universal_newlines=True, timeout=60, check=True)
    return int(result.stdout.strip().splitlines()[-1])


def fit_diagnostics(extract, xmax):
    output = []
    for group in extract['groups']:
        for method in ('method3', 'method5'):
            for sample in ('SE', 'ME'):
                fits = [row['fit'] for row in group['methods'][method][sample] if 0 <= row['x'] < xmax]
                chi2 = [f['chi2']/f['ndf'] for f in fits if finite(f.get('chi2')) and f['ndf'] > 0]
                probabilities = [f['probability'] for f in fits if finite(f.get('probability'))]
                output.append(dict(group=group['id'], method=method, sample=sample, total=len(fits),
                                   valid=sum(f['numericallyValid'] for f in fits),
                                   boundary=sum(f['boundary'] for f in fits),
                                   median_chi2_ndf=statistics.median(chi2) if chi2 else None,
                                   max_chi2_ndf=max(chi2) if chi2 else None,
                                   min_probability=min(probabilities) if probabilities else None))
    return output


def gates(extract):
    nominal = extract['massPartition']['nominal']
    return {'full': (nominal['low'], nominal['high']), 'low': (nominal['low'], nominal['mean']),
            'high': (nominal['mean'], nominal['high']), 'leftSB': (nominal['sideLow'], nominal['low']),
            'rightSB': (nominal['high'], nominal['sideHigh'])}


def build_readme(out, context, result, extract, diagnostics):
    sample = SAMPLES[context['energy']]
    root = Path(result['inputRoot'])
    settings, base = context['settings'], context['base']
    edges = extract['groups'][0]['kstarEdges']
    width = 1000*(edges[1]-edges[0])
    native_width = width/base['rebin_group']
    purity_width = native_width*base['purity_rebin']
    labels = {g['id']: g['label'] for g in result['groups']}
    binned = extract['massPartition']['binned']
    all_cent = next(g for g in result['groups'] if g['id'] == 'all_0_80')
    lines = ['<!-- generated by %s -->' % OWNER,
             '# Λ質量窓を低側・高側に分けたCF比較（%s GeV、KFParticle）' % sample['beam'], '',
             '- [比較レポート PDF](mass_halves_comparison.pdf) / [TeX](mass_halves_comparison.tex)',
             '- 入力: [%s](%s)' % (root.name, os.path.relpath(str(root), str(out))),
             '- %d P24iyのみ。%s / %s subjobの部分統合（欠損・除外 %d 件）。全統計と同一視しません。' %
             (sample['year'], format(sample['included'], ','), format(sample['expected'], ','),
              sample['expected']-sample['included']),
             '- SUMS ID: %s。ビームエネルギーのラベルであり、√sNNではありません。' % sample['jobid'],
             '- 保存済みヒストの後処理のみ。元ROOT、Λ選択cut、event mixing設定は変更していません。', '',
             '## 結果の要点（0–80%）', '',
             '%s ≤ k* < %s GeV/c の積分CF同士の比です。CF点の単純平均ではありません。' %
             tuple(number(v) for v in result['comparison_range']), '']
    for key in METHODS:
        metrics = all_cent['methods'][key]['integrated']
        lines.append('- %s: Low/Full = **%s**、High/Full = **%s**。Fullからの相対変化は低側 **%s**、高側 **%s**。' %
                     (NAMES[key], measure(ratio(metrics, 'low')), measure(ratio(metrics, 'high')),
                      measure(ratio(metrics, 'low'), percent=True), measure(ratio(metrics, 'high'), percent=True)))
    lines += ['', 'NAは「差がない」や0ではなく、必要なfit・規格化などが成立せず比較不能であることを示します。以下の近似と未反映の相関があるため、比の誤差を完全な誤差や最終的な有意度にしないでください。', '',
              '## 共通の条件', '',
              '- k* bin幅: **%s MeV/c**。表示: **0–%s GeV/c**。' % (number(width), number(settings['plot_x_max'])),
              '- 規格化: **%s–%s GeV/c**。表示範囲外でも規格化区間の全binを使用します。' %
              tuple(number(v) for v in result['normalization']),
              '- Method 5の純度推定: **%s MeV/c** bins、W=%s（左右各W bin）、内部の滑動窓幅 **%s MeV/c**。端では窓を短縮。' %
              (number(purity_width), base['sliding_half_width'], number((2*base['sliding_half_width']+1)*purity_width)),
              '- Method 5のCF: 元の **%s MeV/c** カウントそれぞれに所属する純度binのPを掛け、**%s MeV/c** へ再集計します。純度fitをCFの細かいbinで再fitすることはしません。' %
              (number(native_width), number(width)), '',
              '| Window | Exact Maker gate [GeV/c²] | Whole-bin fit / purity gate [GeV/c²] |', '|---|---|---|']
    for r, edge in gates(extract).items():
        opening, closing = ('(', ']') if r == 'rightSB' else ('[', ')' if r in ('low', 'leftSB') else ']')
        lines.append('| %s | %s%s, %s%s | [%s, %s) |' %
                     (r, opening, number(edge[0], 7), number(edge[1], 7), closing,
                      number(binned[r]['lower'], 7), number(binned[r]['upper'], 7)))
    lines += ['', 'Raw/SBsub/CFsubのsignalは候補質量によるexact gateで保存されたヒストを使います。M3/M5とinclusive純度proxyではbin全体を使うため、中心だけでなく外端も異なります。bin内の候補を補間・再配分していません。', '',
              '| Method | 内容 |', '|---|---|']
    lines += ['| %s | %s |' % (NAMES[m], JP_METHODS[m]) for m in METHODS]
    lines += ['', '### CFsubのinclusive純度proxy', '',
              '| Window | P ± statistical error | Actual bin-width scale |', '|---|---:|---:|']
    for r in REGIONS:
        lines.append('| %s | %s | %s |' %
                     (r, measure(result['purity'][r]), number(result['purity'][r]['width_scale'], 7)))
    lines += ['', 'Pの誤差は上表のみでCFsubの曲線・比には伝播しません。これはペア・中央度・k*依存の純度を直接測った値ではなく、平坦背景を仮定したinclusive proxyです。左右SBの非対称性とモデル系統を含みません。', '',
              '## 中央度別の積分CFとfullとの比較', '',
              'Valid binsは表示範囲のFull/Low/Highの順。積分NAでも、別k* binの有効点が図に残る場合があります。', '']
    for group in result['groups']:
        lines += ['### Centrality %s' % group['label'], '',
                  '| Method | C full | C low | C high | Low/Full | High/Full | Valid bins |',
                  '|---|---:|---:|---:|---:|---:|---|']
        lines += ['| %s | %s |' % (NAMES[m], ' | '.join(metric_cells(group['methods'][m]))) for m in METHODS]
        lines.append('')
    lines += ['## 欠測・fit診断', '', '| Centrality | Method | Window | NA / invalid-bin reason |', '|---|---|---|---|']
    invalid_rows = []
    for group in result['groups']:
        for m in METHODS:
            for r in REGIONS:
                text = invalid_reason(group['methods'][m], r)
                if text:
                    invalid_rows.append('| %s | %s | %s | %s |' %
                                        (group['label'], NAMES[m], r, text.replace('|', '/')))
    lines += invalid_rows or ['| All | All | All | None |']
    lines += ['', '今回のMethod 3の欠測は低k*のSE質量分布での統計不足、またはfit背景が候補総数を超えることによります。規格化区間が不成立という意味ではありません。数値収束したfitでも0 ≤ B ≤ Nの条件を満たさなければ、その窓は無効です。', '']
    first = extract['groups'][0]['methods']['method3']['SE'][0]
    if first['full']['background'] is not None:
        lines += ['0–80%%・最初のk* bin（0–%s GeV/c）のSEでは、full窓のN=%sに対しfit背景B=%sです。低側のみ有効でも、high/fullが無効ならfullとの積分比較はNAになります。' %
                  (number(edges[1]),number(first['full']['centralNtot']),number(first['full']['background'],7)), '']
    lines += ['', '表示範囲内の質量fit診断。新しいcutは追加していません。ROOTのfit確率0は数値underflowを含み得ます。', '',
              '| Centrality | Method | Sample | Numerically valid / total | Boundary | Median χ²/ndf | Max χ²/ndf | Min fit probability |',
              '|---|---|---|---:|---:|---:|---:|---:|']
    for d in diagnostics:
        lines.append('| %s | %s | %s | %d/%d | %d | %s | %s | %s |' %
                     (labels[d['group']], NAMES[d['method']], d['sample'], d['valid'], d['total'], d['boundary'],
                      number(d['median_chi2_ndf'], 4), number(d['max_chi2_ndf'], 4), number(d['min_probability'], 3)))
    lines += ['', '## 誤差と解釈上の注意', ''] + ['- '+text for text in JP_CAVEATS]
    lines += ['', '## 図と再現記録', '',
              '各PDFは2ページ。1ページ目:0–80%、20–80%、0–20%。2ページ目:0–10%、10–20%、20–60%。上段CF、下段Low/Full・High/Full。', '']
    lines += ['- [%s](mass_halves_%s.pdf) / [PNG preview](mass_halves_%s.png)' % (NAMES[m], m, m) for m in METHODS]
    lines += ['', '- context.json: 入力fingerprint、mainconf・YAML hash、後処理条件。',
              '- mass_half_extract.json: カウントと誤差、質量bin境界、fitと共分散。',
              '- comparison.json / cf_points.tsv: 計算済みCF・比・積分比較と欠測理由。',
              '- plotcanvases.root: 配置修正用の10 canvases。元解析ROOTとは別ファイル。',
              '- report_metadata.json: JSON・図のhash、fit診断、summaryページ数。',
              '- macros/ / logs/: 実行マクロとROOT/LaTeXログ。', '',
              'JSONと図からレポートのみ再生成:', '', '~~~bash',
              'python3 tools/kfparticle_summary/report_mass_halves.py \\',
              '  --output-dir %s' % out, '~~~', '']
    return '\n'.join(lines)


def build_tex(context, result, extract):
    sample = SAMPLES[context['energy']]
    base, settings = context['base'], context['settings']
    edges = extract['groups'][0]['kstarEdges']
    width = 1000*(edges[1]-edges[0])
    native_width = width/base['rebin_group']
    purity_width = native_width*base['purity_rebin']
    all_cent = next(g for g in result['groups'] if g['id'] == 'all_0_80')
    binned = extract['massPartition']['binned']
    date = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).strftime('%Y-%m-%d %H:%M JST')
    parts = [r'''\documentclass[10pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{geometry,amsmath,longtable,graphicx,pdfpages,hyperref}
\geometry{margin=1.7cm}
\hypersetup{colorlinks=true,linkcolor=blue,urlcolor=blue}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.45em}
\title{KFParticle $\Lambda$--$^{4}\mathrm{He}$\\Low / High Invariant-Mass CF Comparison}
\author{Post-processing QA report}
\date{''' + escape(date) + r'''}
\begin{document}
\maketitle
\section{Sample, scope and common choices}
''',
             'Au+Au fixed-target sample, beam-energy label %s GeV; %d P24iy only. The beam label is not $\\sqrt{s_{NN}}$. Partial merge: %s/%s subjobs; %s excluded or missing.\n' %
             (sample['beam'], sample['year'], format(sample['included'], ','), format(sample['expected'], ','), sample['expected']-sample['included']),
             r'Input: {\small\path{' + Path(result['inputRoot']).name + r'}}. No event reanalysis, selection-cut change or original ROOT modification.'+'\n',
             r'$k^*$ bin width: '+number(width)+r' MeV/$c$; display: $0$--$'+number(settings['plot_x_max'])+
             r'$ GeV/$c$; normalization: $'+number(result['normalization'][0])+'$--$'+number(result['normalization'][1])+
             r'$ GeV/$c$. The normalization interval remains in use even outside the displayed range.'+'\n',
             r'The integrated comparison uses $'+number(result['comparison_range'][0])+r'\leq k^*<'+
             number(result['comparison_range'][1])+r'$ GeV/$c$. It is the normalized ratio of summed counts, not an arithmetic average of CF points. Each mass window is independently normalized.'+'\n',
             r'\section{Inclusive 0--80\% comparison}'+'\n',
             'Ratios compare each half to full within the same background method. NA means unavailable, not zero or agreement. Errors are approximate (see limitations below).\n', r'\begin{itemize}'+'\n']
    for m in METHODS:
        metrics = all_cent['methods'][m]['integrated']
        parts.append(r'\item '+NAMES[m]+': low/full = '+measure(ratio(metrics, 'low'), True)+
                     '; high/full = '+measure(ratio(metrics, 'high'), True)+'.\n')
    parts += [r'\end{itemize}\section{Mass windows and methods}'+'\n',
              r'Exact gates use candidate masses. Mass-fit projections and inclusive purity use whole stored mass bins, without split-bin interpolation:'+'\n',
              r'\begin{center}\begin{tabular}{lcc}\hline Window & Exact nominal edges & Whole-bin edges\\\hline'+'\n']
    for r, edge in gates(extract).items():
        parts.append('%s & %s--%s & %s--%s' %
                     (r, number(edge[0], 7), number(edge[1], 7), number(binned[r]['lower'], 7),
                      number(binned[r]['upper'], 7))+r'\\'+'\n')
    parts += [r'\hline\end{tabular}\end{center}'+'\n',
              r'Edges are in GeV/$c^2$. Exact low excludes the mean, exact high includes it; the full outer edges are included. Whole-bin intervals are disjoint and half-open. Raw/SBsub/CFsub signal histograms use exact gates, while Methods 3/5 include this documented mass-binning difference.'+'\n',
              r'\begin{description}'+'\n',
              r'\item[Raw] Normalized signal SE/ME, without background correction.'+'\n',
              r'\item[SBsub] Subtract summed sideband counts from SE and ME before normalization, using $1/3$ for full and $1/6$ for each half. A locally flat background is assumed.'+'\n',
              r'\item[CFsub] $C_r^{\rm corr}=[C_r-(1-P_r)C_{\rm SB}]/P_r$. Each window uses its own current inclusive purity proxy.'+'\n',
              r'\item[Method 3] Fit each SE/ME mass projection with Gaussian plus quadratic, then form the normalized ratio of $N-B$. Low/high background integrals share one fit; full uses their sum.'+'\n',
              r'\item[Method 5] Purity fits use '+number(purity_width)+r' MeV/$c$ coarse bins and $W='+str(base['sliding_half_width'])+
              r'$ neighbors per side: interior window '+number((2*base['sliding_half_width']+1)*purity_width)+
              r' MeV/$c$, truncated at edges. Each original '+number(native_width)+
              r' MeV/$c$ count receives the purity of its containing coarse bin, then contributes to the '+number(width)+
              r' MeV/$c$ CF bin. No fine-bin refit is performed. Full has its own window purity; its corrected yield need not equal the sum of corrected half yields.'+'\n',
              r'\end{description}'+'\n',
              r'\begin{center}\begin{tabular}{lcc}\hline Window & Inclusive purity proxy & Actual width factor\\\hline'+'\n']
    for r in REGIONS:
        parts.append(r+' & '+measure(result['purity'][r], True)+' & '+
                     number(result['purity'][r]['width_scale'], 7)+r'\\'+'\n')
    parts += [r'\hline\end{tabular}\end{center}'+'\n',
              'These purity errors are not propagated into CFsub. These are inclusive proxies, not measured pair-, centrality- or momentum-dependent purities; sideband asymmetry and background-model systematics remain.\n',
              r'\clearpage\section{Integrated CF values and half/full ratios}'+'\n',
              'Valid bins lists full/low/high counts in the displayed interval. A partial curve does not guarantee that the complete integration interval is valid. Invalid bins are not silently removed from the requested integration.\n',
              r'{\scriptsize\setlength{\tabcolsep}{3pt}\begin{longtable}{lrrrrrl}'+'\n',
              r'\hline Method & $C_{\rm full}$ & $C_{\rm low}$ & $C_{\rm high}$ & Low/Full & High/Full & Valid bins\\\hline\endhead'+'\n']
    for group in result['groups']:
        parts.append(r'\multicolumn{7}{l}{\textbf{Centrality '+escape(group['label'])+r'}}\\'+'\n')
        for m in METHODS:
            parts.append(NAMES[m]+' & '+' & '.join(metric_cells(group['methods'][m], True))+r'\\'+'\n')
        parts.append(r'\hline'+'\n')
    parts += [r'\end{longtable}}'+'\n', r'\section{Validity and uncertainty limitations}'+'\n',
              'An invalid fit in a required normalization bin makes its entire CF unavailable. Invalid comparison bins make the integrated contrast NA, not zero. Numerical convergence does not prove mass-model adequacy: large chi-square/ndf or small fit probabilities prohibit precision-model claims. Detailed fit diagnostics and NA reasons are in the Japanese README and extraction JSON.\n',
              r'First-order propagation retains shared full/half counts, sidebands, normalization integrals, the low/high background covariance from each common fit, and reuse of one coarse-bin purity by multiple CF bins. It cannot recover event / track sharing or mixed-pair reuse correlations. Covariance between observed $N$ and fitted background $B$, covariance between overlapping fits in different coarse bins, CFsub purity uncertainty and model systematics are not included. These comparisons are QA contrasts, not final significance tests.'+'\n',
              'Centrality groups overlap and are not independent. Low/high CFs cannot be directly added. Small TH2F rounding residuals are recorded without silently replacing the full histogram.\n',
              r'\section{Figure guide}'+'\n',
              r'The next ten pages show two pages per method: Raw, SBsub, CFsub, Method 3, Method 5. First page: 0--80\%, 20--80\%, 0--20\%; second: 0--10\%, 10--20\%, 20--60\%. Upper panels: full (black), low (red), high (blue). Lower panels: low/full and high/full with supplied covariance errors. Failed bins are omitted, not displayed as measured zero.'+'\n']
    for m in METHODS:
        parts.append(r'\includepdf[pages=-,fitpaper=true,pagecommand={}]{mass_halves_'+m+'.pdf}\n')
    parts.append(r'\end{document}'+'\n')
    return ''.join(parts)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    out = parser.parse_args(argv).output_dir.resolve(strict=True)
    files = [out/name for name in ('context.json', 'comparison.json', 'mass_half_extract.json')]
    context, result, extract = [json.loads(path.read_text()) for path in files]
    sample = SAMPLES[context['energy']]
    if not context['input']['path'] == result['inputRoot'] == extract['inputRoot']:
        raise ValueError('Input ROOT mismatch between context, comparison and extraction.')
    expected = 'auau%s_anaFemtoLambda_4He_KFParticle_%dP24iy_partial_merge.root' % (context['energy'], sample['year'])
    if Path(result['inputRoot']).name != expected:
        raise ValueError('Coverage is unrecorded for this ROOT filename/production.')
    if result['normalization'] != context['base']['normalization'] or result['comparison_range'] != context['settings']['comparison_range']:
        raise ValueError('Ranges differ between result and context.')
    figures = {}
    for m in METHODS:
        path = out/('mass_halves_'+m+'.pdf')
        pages = pdf_pages(path)
        if pages != 2:
            raise ValueError('Expected two pages: '+str(path))
        figures[path.name] = dict(pages=pages, sha256=digest(path), bytes=path.stat().st_size)
    signature = dict(owner=OWNER, inputRoot=result['inputRoot'], energy=context['energy'])
    metadata_path = out/'report_metadata.json'
    if metadata_path.exists():
        old = json.loads(metadata_path.read_text())
        if any(old.get(k) != v for k, v in signature.items()):
            raise ValueError('Existing report belongs to another input or generator.')
    elif any((out/name).exists() for name in ('README.md', 'mass_halves_comparison.tex', 'mass_halves_comparison.pdf')):
        raise ValueError('Refusing to overwrite unowned report products.')
    diagnostics = fit_diagnostics(extract, context['settings']['plot_x_max'])
    metadata = dict(signature, phase='building', sample=sample, figures=figures,
                    sources={path.name: digest(path) for path in files}, fit_diagnostics=diagnostics)
    save(metadata_path, metadata)
    (out/'README.md').write_text(build_readme(out, context, result, extract, diagnostics))
    tex = out/'mass_halves_comparison.tex'
    tex.write_text(build_tex(context, result, extract))
    latex = shutil.which('pdflatex')
    if not latex:
        raise RuntimeError('pdflatex is required.')
    (out/'logs').mkdir(exist_ok=True)
    for iteration in (1, 2):
        with (out/'logs'/('mass_halves_pdflatex_%d.log' % iteration)).open('w') as stream:
            subprocess.run([latex, '-no-shell-escape', '-interaction=nonstopmode', '-halt-on-error',
                            '-file-line-error', tex.name], cwd=str(out), stdout=stream,
                           stderr=subprocess.STDOUT, timeout=120, check=True)
    pdf = out/'mass_halves_comparison.pdf'
    pages = pdf_pages(pdf)
    if pages <= 10:
        raise ValueError('Report is missing front matter or figures.')
    metadata.update(phase='complete', report=dict(name=pdf.name, pages=pages, sha256=digest(pdf),
                                                 bytes=pdf.stat().st_size),
                    generator_sha256=digest(Path(__file__)), tex_sha256=digest(tex))
    save(metadata_path, metadata)
    print('Created %s (%s pages)' % (pdf, pages))


if __name__ == '__main__':
    try:
        main()
    except (ValueError, KeyError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print('report_mass_halves: %s' % exc, file=sys.stderr)
        sys.exit(1)
