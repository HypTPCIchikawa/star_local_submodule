#!/usr/bin/env python3
"""Build a provenance-aware KFParticle Lambda--4He PDF from generated QA plots."""
import argparse
from collections import Counter
import csv
import datetime
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import sys


SAMPLES = {
    '13p5': {'year': 2020, 'included': 542, 'expected': 545,
             'jobid': '1213058CE62BCF113B2E3C75ACBD7053',
             'beam': '13.5', 'stem': 'summary_4He_KFParticle'},
    '3p85': {'year': 2021, 'included': 26956, 'expected': 26994,
             'jobid': '0BEFC64E62F19BBA534A129CABC9B106',
             'beam': '3.85', 'stem': 'summary_4He_KFParticle_3p85'},
}
REQUIRED_PDFS = {
    'Lambda_InvMass.pdf': 1,
    'kstar_sum_4He.pdf': 12,
    'kstar_sum_CFsub_4He.pdf': 12,
    'kstar_sum_CFsub_method3_4He.pdf': 6,
    'kstar_sum_CFsub_method5_4He.pdf': 6,
    'kstar_cent_CFsub_method3_4He.pdf': 3,
    'kstar_cent_CFsub_method5_4He.pdf': 3,
    'Fit_LL_4He.pdf': 2,
    'Fit_LL_CFsub_4He.pdf': 2,
    'Fit_LL_method3_4He.pdf': 2,
    'Fit_LL_method5_4He.pdf': 2,
}
OWNER = 'tools/kfparticle_summary/build_report.py'


def esc(value):
    """Escape text, not mathematics, for LaTeX."""
    chars = {'\\': r'\textbackslash{}', '&': r'\&', '%': r'\%', '$': r'\$',
             '#': r'\#', '_': r'\_', '{': r'\{', '}': r'\}', '~': r'\textasciitilde{}',
             '^': r'\textasciicircum{}'}
    return ''.join(chars.get(c, c) for c in str(value))


def number(value, digits=6):
    if value is None:
        return 'unrecorded'
    if not isinstance(value, (float, int)) or not math.isfinite(value):
        raise ValueError('Nonfinite or nonnumeric measurement: %r' % value)
    return format(value, '.%dg' % digits)


def count(value):
    return 'unrecorded' if value is None else format(value, ',.0f')


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def pdf_pages(path):
    with path.open('rb') as handle:
        if handle.read(5) != b'%PDF-':
            raise ValueError('Invalid or incomplete PDF: %s' % path)
    gs = shutil.which('gs')
    if not gs:
        raise RuntimeError('Ghostscript (gs) is required to verify PDF page counts.')
    result = subprocess.run(
        [gs, '-q', '-dNODISPLAY', '-dSAFER', '--permit-file-read=' + str(path),
         '-sPDFFile=' + str(path), '-c',
         'PDFFile (r) file runpdfbegin pdfpagecount = quit'],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, universal_newlines=True,
        timeout=60, check=True)
    return int(result.stdout.strip().splitlines()[-1])


def audits(out, manifest):
    """Keep final repeated fit result per object; never invent a successful fit."""
    rows = []
    sources = []
    for macro in manifest.get('macros', []):
        name = macro.get('name', Path(macro.get('path', macro.get('generated', 'unknown'))).name)
        audit = macro.get('audit_file')
        if not audit:
            continue
        path = Path(audit)
        if not path.is_absolute():
            path = out / path
        if not path.is_file():
            raise RuntimeError('Missing declared fit audit: %s' % path)
        with path.open() as handle:
            latest = {}
            mass_counts = Counter()
            for row in csv.DictReader(handle, delimiter='\t'):
                if row.get('type') == 'll_fit':
                    latest[row['object']] = row
                elif row.get('type') == 'mass_fit':
                    mass_counts[row['reason']] += 1
        sources.append({'path': str(path), 'sha256': sha256(path), 'macro': name,
                        'mass_fit_counts': dict(mass_counts)})
        for obj, row in sorted(latest.items()):
            row['macro'] = name
            rows.append(row)
    return rows, sources


def figure(name, page, caption, width='0.92'):
    return (r'\begin{figure}[!ht]\centering' + '\n' +
            r'\includegraphics[width=' + width + r'\textwidth,page=' + str(page) +
            ']{' + name + '}\n' + r'\caption{' + caption + '}\n' + r'\end{figure}' + '\n')


def pair(left, right):
    parts = [r'\begin{figure}[!ht]\centering']
    for i, (name, page, caption) in enumerate((left, right)):
        if i:
            parts.append(r'\hfill')
        parts.extend([r'\begin{minipage}{0.49\textwidth}\centering',
                      r'\includegraphics[width=\textwidth,page=' + str(page) + ']{' + name + '}',
                      r'\caption{' + caption + '}', r'\end{minipage}'])
    parts.append(r'\end{figure}')
    return '\n'.join(parts) + '\n'


def method5_grid(settings, output_bin_width):
    """Interpret old saved manifests explicitly; new YAML requires purity_rebin."""
    independent = 'purity_rebin' in settings
    rebin = settings['purity_rebin'] if independent else settings['rebin_group']
    purity_width = output_bin_width/settings['rebin_group']*rebin
    return dict(independent=independent, rebin=rebin, width=purity_width,
                window_width=(2*settings['sliding_half_width']+1)*purity_width)


def build_tex(energy, info, manifest, ll_rows, audit_sources):
    sample = SAMPLES[energy]
    w = info['massWindow']
    p = info['inclusivePurityProxy']
    fit = info['gausPol2Comparison']
    counters = info['eventCounters']
    root_name = Path(info['inputRoot']).name
    purity = number(p['purity'])
    purity_error = format(float(p['purityError']), '.6f')
    left_density = w['left']['count'] / (w['left']['actualUpper'] - w['left']['actualLower'])
    right_density = w['right']['count'] / (w['right']['actualUpper'] - w['right']['actualLower'])
    settings = manifest['settings']
    kstar = next(h for h in info['histograms']
                 if h['name'] == 'true/hKstar_4He_CentBin0' and h['present'])
    bin_width = (kstar['xMax'] - kstar['xMin']) / kstar['xBins'] * settings['rebin_group']
    m5 = method5_grid(settings,bin_width)
    ll = settings['ll_fit']
    date = datetime.datetime.now(datetime.timezone(datetime.timedelta(hours=9))).strftime('%Y-%m-%d %H:%M JST')
    parts = [r'''\documentclass[11pt,a4paper]{article}
\usepackage[utf8]{inputenc}
\usepackage[T1]{fontenc}
\usepackage{graphicx,geometry,amsmath,longtable,hyperref}
\geometry{margin=2.0cm}
\hypersetup{colorlinks=true,urlcolor=blue,linkcolor=blue}
\setlength{\parindent}{0pt}
\setlength{\parskip}{0.5em}
\title{KFParticle $\Lambda$--$^{4}\mathrm{He}$\\Correlation Analysis: QA Summary}
\author{STAR analysis working report}
\date{''' + esc(date) + r'''}
\begin{document}
\maketitle
\section{Sample and scope}
''',
        'Au+Au fixed-target sample, beam-energy label ' + esc(sample['beam']) +
        r' GeV; production ' + str(sample['year']) + r' P24iy. This beam label is \emph{not} '
        r'a claim that $\sqrt{s_{NN}}$ equals this value.' + '\n',
        r'\textbf{Partial merge:} ' + format(sample['included'], ',') + '/' +
        format(sample['expected'], ',') + ' subjobs are included; ' +
        str(sample['expected'] - sample['included']) +
        ' missing or excluded subjobs are not represented. No missing yield is filled in.\n',
        r'\par{\small\textbf{ROOT:} \path{' + root_name + r'}\par' + '\n' +
        r'\textbf{SUMS ID:} \texttt{' + sample['jobid'] + r'}\par}' + '\n',
        'This report reuses the structure of the earlier tools summaries, but reads only the '
        'KFParticle merged ROOT above. Legacy ROOTs, legacy PDFs and selection cuts are unchanged. '
        'The methods below are post-processing cross-checks, not new production selection or '
        'a final physics result.\n',
        r'\begin{center}\begin{tabular}{lr}\hline Measurement & Stored value\\\hline' + '\n',
        'Input-event counter & ' + count(counters.get('inputEvents')) + r'\\' + '\n',
        'Reconstructed-event counter & ' + count(counters.get('reconstructedEvents')) + r'\\' + '\n',
        'Selected-$\Lambda$ counter & ' + count(counters.get('selectedLambdaCandidates')) + r'\\' + '\n',
        'Signal-window candidate counts & ' + count(w['signal']['count']) + r'\\' + '\n',
        'Width-scaled sideband background & ' + count(p['background']) + r'\\' + '\n',
        r'Inclusive purity proxy $P$ & $' + purity + r'\pm' + purity_error + r'$\\\hline' + '\n',
        r'\end{tabular}\end{center}' + '\n',
        figure('Lambda_InvMass.pdf', 1,
               r'Inclusive $p\pi^-$ invariant mass. Vertical markers show the unchanged production '
               r'signal and sideband gates; the fit is a diagnostic only.', '0.58'),
        r'\clearpage\section{Mass gates and purity convention}' + '\n',
        r'The Maker gate parameters are $\mu=' + number(w['mean']) + r'$ GeV/$c^2$, '
        r'$\sigma=' + number(w['sigma']) + r'$ GeV/$c^2$, with signal '
        r'$|M-\mu|<' + number(w['nSigma']) + r'\sigma$. The left/right sidebands extend '
        r'from $' + number(w['nSigma']) + r'\sigma$ to $' +
        number(w['nSigma'] * w['sidebandOuterFactor']) + r'\sigma$. These fixed gate parameters '
        r'are not refitted to redefine the event selection.' + '\n',
        r'\begin{center}\begin{tabular}{lcc}\hline Region & Nominal edges & Integrated bin edges\\\hline' + '\n',
    ]
    bounds = [('Signal', w['nominalSignalLower'], w['nominalSignalUpper'], 'signal'),
              ('Left SB', w['nominalSidebandLower'], w['nominalSignalLower'], 'left'),
              ('Right SB', w['nominalSignalUpper'], w['nominalSidebandUpper'], 'right')]
    for label, low, high, key in bounds:
        parts.append('%s & %s--%s & %s--%s' %
                     (label, number(low), number(high), number(w[key]['actualLower']),
                      number(w[key]['actualUpper'])) + r'\\' + '\n')
    parts.extend([
        r'\hline\end{tabular}\end{center}' + '\n',
        r'All mass edges are in GeV/$c^2$. For the inclusive diagnostic, whole bins with centers '
        r'inside each nominal interval are integrated, without double counting boundary bins. '
        r'The binned diagnostic therefore only approximates the unbinned Maker gate.' + '\n',
        r'\begin{equation}P_{\rm inclusive}=1-\alpha_{\rm bins}\frac{N_L+N_R}{N_W},\qquad '
        r'\alpha_{\rm bins}=\frac{\Delta M_W}{\Delta M_L+\Delta M_R}.\end{equation}' + '\n',
        r'Here $\alpha_{\rm bins}=' + number(p['actualWidthScale']) + r'$, whereas the nominal '
        r'sideband scale used with the stored pair histograms is $\alpha=' +
        number(p['nominalWidthScale']) + r'$. The estimate assumes a locally flat background '
        r'across the sidebands; its displayed error is statistical and excludes model uncertainty.' + '\n',
        r'\textbf{CFsub uses the measured value $P=' + purity + r'$ above, not an inherited fixed '
        r'S/N.} This is an inclusive-candidate proxy, not a measured pair-, centrality- or '
        r'$k^*$-dependent purity. Its uncertainty is not propagated through the legacy constant-purity '
        r'CFsub curves; those errors must not be read as total uncertainties.' + '\n',
        r'The current left and right sidebands contain ' + count(w['left']['count']) +
        ' and ' + count(w['right']['count']) + r' candidates. Their right/left count-density ratio '
        r'is ' + number(right_density / left_density) + r'. This measured asymmetry warns against '
        r'treating the flat-background proxy as a precision purity determination; its small '
        r'statistical error does not cover interpolation or shape uncertainty.' + '\n',
        r'The independent Gaussian-plus-quadratic fit has status ' + str(fit['status']) +
        ', covariance status ' + str(fit['covarianceStatus']) +
        r', $\chi^2/\mathrm{ndf}=' + number(fit['chi2']) + '/' + str(fit['ndf']) +
        r'$, and fit probability ' + number(fit.get('chi2Probability')) +
        '. Its comparison purity is ' + number(fit['purity']) +
        r'. This is a diagnostic value only, not adopted for CF correction. Numerical convergence '
        r'does not establish model adequacy: a very small probability or large $\chi^2$/ndf '
        r'indicates that the simple shape model is inadequate. The same limitation must be '
        r'considered for the projection fits in Methods 3 and 5.' + '\n',
        r'\section{Common post-processing choices}' + '\n',
        r'The comparison macros use $k^*$ rebinning by ' + str(settings['rebin_group']) +
        ' (' + number(1000 * bin_width) + r' MeV/$c$ bins) and normalize SE/ME '
        r'in $' + number(settings['normalization'][0]) + r'<k^*<' +
        number(settings['normalization'][1]) + r'$ GeV/$c$. The displayed $k^*$ range is $' +
        number(settings['plot_x_range'][0]) + r'$--$' + number(settings['plot_x_range'][1]) +
        r'$ GeV/$c$; normalization still uses the full stated interval outside the displayed range. '
        r'These are summary post-processing choices, distinct '
        r'from the Maker CF settings. Raw histograms are summed \emph{before} a grouped CF is formed.' + '\n',
        r'Method 5 estimates purity on ' + number(m5['width']*1000) + r' MeV/$c$ bins with $W=' +
        str(settings['sliding_half_width']) + r'$ neighboring bins on each side: the nominal interior '
        r'fit window spans ' + number(m5['window_width']*1000) + r' MeV/$c$. Edge windows are truncated. ' +
        (r'This purity grid is independent of the CF output grid. Each original saved $k^*$ bin is '
         r'weighted by the purity of its own coarse cell before being summed into an output CF bin. '
         r'Thus an output bin crossing a purity-cell boundary retains both contributions; no '
         r'center-only assignment or interpolation is used.' if m5['independent'] else
         r'This historical manifest uses the original coupled grid: purity and CF bins have the same width.') + '\n',
        r'Two group sets are used consistently: (20--80\%, 0--20\%) and '
        r'(0--10\%, 10--20\%, 20--60\%). In the first group-set plots below, '
        r'the left panel is 20--80\% and the right panel is 0--20\%. Fine-bin plots show '
        r'all nine centralities, from 70--80\% to 0--5\%.' + '\n',
        r'\clearpage\section{Definitions and limitations}' + '\n',
        r'Let $W$ denote the signal mass gate, $L+R$ the summed sidebands, and ${\cal N}$ '
        r'the normalization obtained in the stated $k^*$ interval. For raw sideband subtraction,' + '\n',
        r'\begin{equation}C_{\rm SBsub}(k^*)={\cal N}_{\rm sub}'
        r'\frac{N_{\rm SE}^{W}-\alpha(N_{\rm SE}^{L}+N_{\rm SE}^{R})}'
        r'{N_{\rm ME}^{W}-\alpha(N_{\rm ME}^{L}+N_{\rm ME}^{R})}.\end{equation}' + '\n',
        r'For constant-purity CF subtraction, signal and summed-sideband CFs are normalized separately:' + '\n',
        r'\begin{equation}C_{\rm CFsub}(k^*)=\frac{C_W(k^*)-(1-P)C_{L+R}(k^*)}{P}.\end{equation}' + '\n',
        r'The two sidebands are added at the count level; $C_{L+R}$ is not the simple average '
        r'of independently normalized left and right CFs. The methods assume that the sideband '
        r'background represents the signal-region background, in count shape or correlation '
        r'shape respectively. Neither assumption is validated by plotting alone.' + '\n',
        r'Method 3 fits each mass projection with a Gaussian plus quadratic background, '
        r'separately for SE and ME, and constructs a normalized ratio of fitted signal yields. '
        r'Method 5 estimates the purity using overlapping neighboring $k^*$ windows and applies '
        r'it to the central bin. Overlap introduces covariance between adjacent estimates; '
        r'the legacy diagonal-error treatment does not include that covariance.' + '\n',
        r'\begin{itemize}' + '\n',
        r'\item Statistical errors do not include CF-normalization uncertainty, background-model '
        r'systematics, or the constant inclusive-purity uncertainty.' + '\n',
        r'\item Direct subtraction can produce small or nonpositive denominators. Sparse or failed '
        r'mass fits do not constitute zero physical signal. Review fit audit logs and invalid-bin '
        r'diagnostics before using any plotted point quantitatively.' + '\n',
        r'\item Whole mass bins used in projection fits need not match the exact stored signal-gate '
        r'pair histograms at the boundary. A method difference can therefore include bin-edge effects.' + '\n',
        r'\item Missing subjobs and possible coverage biases are not corrected. This report does '
        r'not compare efficiencies or infer purity for the missing sample.' + '\n',
        r'\item No method is declared superior from smoothness alone, and no source radius is '
        r'established as a final physical result.' + '\n',
        r'\end{itemize}' + '\n',
        r'\clearpage\section{Raw sideband subtraction: grouped spectra}' + '\n',
        pair(('kstar_sum_4He.pdf', 3, 'SE signal-gate and sideband distributions.'),
             ('kstar_sum_4He.pdf', 4, 'ME signal-gate and sideband distributions.')),
        pair(('kstar_sum_4He.pdf', 5, 'Background-subtracted SE and ME spectra.'),
             ('kstar_sum_4He.pdf', 6, 'Normalized CF after raw sideband subtraction.')),
        r'\clearpage\section{Constant-purity CF subtraction}' + '\n',
        pair(('kstar_sum_CFsub_4He.pdf', 1, 'SE and ME in the signal mass gate.'),
             ('kstar_sum_CFsub_4He.pdf', 2, 'Uncorrected signal-gate CF.')),
        pair(('kstar_sum_CFsub_4He.pdf', 5, 'Signal-gate CF and summed-sideband CF.'),
             ('kstar_sum_CFsub_4He.pdf', 6, 'CFsub using the current inclusive-purity proxy.')),
        r'\clearpage\section{Method 3: mass-fit signal yields}' + '\n',
        pair(('kstar_sum_CFsub_method3_4He.pdf', 1, 'Mass-fit SE and ME signal yields.'),
             ('kstar_sum_CFsub_method3_4He.pdf', 2, r'SE and ME fitted purity versus $k^*$.')),
        figure('kstar_sum_CFsub_method3_4He.pdf', 3,
               'Normalized fitted-signal-yield ratio. Fit quality must be checked independently.', '0.85'),
        r'\clearpage\section{Method 5: sliding-window comparison}' + '\n',
        r'Both fine-centrality figures below show all nine centralities. The 0--5\% bin '
        r'is the bottom-right panel; page 2 is a purity summary, not a single-centrality page.' + '\n',
        pair(('kstar_cent_CFsub_method3_4He.pdf', 2,
              'Method 3 purity in nine centrality bins.'),
             ('kstar_cent_CFsub_method5_4He.pdf', 2,
              'Method 5 sliding-window purity in the same bins.')),
        pair(('kstar_sum_CFsub_method3_4He.pdf', 3, 'Method 3 grouped CF.'),
             ('kstar_sum_CFsub_method5_4He.pdf', 3, 'Method 5 grouped CF; adjacent-bin covariance is omitted.')),
        r'\clearpage\section{Exploratory Lednicky--Lyuboshits fits}' + '\n',
        r'The blue fits hold $a_0=' + number(ll['fixed_a0']) + r'$ fm and $r_{\rm eff}=' +
        number(ll['fixed_reff']) + r'$ fm fixed and fit '
        r'$R$ and a normalization. The red fits also vary $a_0$ and $r_{\rm eff}$. '
        r'These reproduce the existing model/sign convention, not an independently validated '
        r'interaction extraction. The fit interval is $' + number(ll['range'][0]) + r'<k^*<' +
        number(ll['range'][1]) + r'$ GeV/$c$: only about ' +
        number((ll['range'][1] - ll['range'][0]) / bin_width) + ' bins of width ' +
        number(bin_width * 1000) + r' MeV/$c$ are available before invalid-bin exclusions. Four free parameters '
        r'can therefore be weakly constrained. Parameter bounds, covariance quality and fit status '
        r'matter in addition to the displayed $\chi^2$/ndf.' + '\n',
        pair(('Fit_LL_4He.pdf', 1, 'LL fits to raw-SB-subtracted CFs.'),
             ('Fit_LL_CFsub_4He.pdf', 1, 'LL fits to constant-purity CFsub CFs.')),
        pair(('Fit_LL_method3_4He.pdf', 1, 'LL fits to Method 3 CFs.'),
             ('Fit_LL_method5_4He.pdf', 1, 'LL fits to Method 5 CFs; window covariance is not fitted.')),
        r'\clearpage\section{Additional centrality groups}' + '\n',
        r'The following pages use the same groups for every method: 0--10\%, 10--20\% and '
        r'20--60\%. This corrects the inconsistent third-group definitions in some legacy macros.' + '\n',
        pair(('kstar_sum_4He.pdf', 12, 'SBsub, second group set.'),
             ('kstar_sum_CFsub_4He.pdf', 12, 'Constant-purity CFsub, second group set.')),
        pair(('kstar_sum_CFsub_method3_4He.pdf', 6, 'Method 3, second group set.'),
             ('kstar_sum_CFsub_method5_4He.pdf', 6, 'Method 5, second group set.')),
        r'\clearpage\section{LL fit audit and reproducibility}' + '\n',
        r'Fit status 0 alone does not establish model adequacy. Covariance status and boundary '
        r'flags are retained separately below. Repeated optimizer calls are reduced to '
        r'the last result for each fit object. Full parameters and diagnostics remain in the TSV files.' + '\n',
    ])
    mass_sources = [item for item in audit_sources if item['mass_fit_counts']]
    if mass_sources:
        parts.append(r'{\small\begin{longtable}{lrrrr}\hline Mass-fit macro & Total & Accepted & Boundary & Rejected\\\hline\endhead' + '\n')
        for item in mass_sources:
            counts = item['mass_fit_counts']
            total = sum(counts.values())
            valid = counts.get('valid', 0) + counts.get('valid_boundary_warning', 0)
            name = item['macro'].replace('LambdaNuclearId_relamom_', '').replace('.C', '')
            parts.append('%s & %d & %d & %d & %d' %
                         (esc(name), total, valid, counts.get('valid_boundary_warning', 0),
                          total - valid) + r'\\' + '\n')
        parts.append(r'\hline\end{longtable}}' + '\n' +
                     r'Accepted here means only the implemented numerical checks passed; '
                     r'it is not a goodness-of-fit or physics validation. Boundary cases are '
                     r'included in accepted counts and are shown separately as warnings.' + '\n')
    if ll_rows:
        parts.append(r'For object \texttt{fit1\_s\_g}, $a_0$ and $r_{\rm eff}$ are fixed; '
                     r'\texttt{fit2\_s\_g} is the all-free fit. The zero-based group-set index '
                     r'$s=0$ means (20--80\%, 0--20\%) and $s=1$ means '
                     r'(0--10\%, 10--20\%, 20--60\%); $g$ is the position in that list.' + '\n')
        parts.append(r'{\scriptsize\begin{longtable}{llrrrr}\hline Macro / fit & Object & Status & Cov. & Bound & $R$ (fm)\\\hline\endhead' + '\n')
        for row in ll_rows:
            short = row['macro'].replace('Fit_LambdaNuclearId_relamom_sum', '').replace('.C', '') or 'SBsub'
            if short.startswith('_'):
                short = short[1:]
            radius = number(float(row['value']), 4)
            if row.get('reason') == 'insufficient_valid_points':
                radius = 'not fitted'
            elif row.get('reason') != 'valid_exploratory':
                radius += '*'
            parts.append('%s & %s & %s & %s & %s & %s' %
                         (esc(short), esc(row['object']), esc(row.get('status', '?')),
                          esc(row.get('covariance', '?')), esc(row.get('boundary', '?')),
                          esc(radius)) + r'\\' + '\n')
        parts.append(r'\hline\end{longtable}}' + '\n' +
                     r'An asterisk marks an unreliable optimizer result. Not fitted suppresses '
                     r'initial parameters left behind when too few valid points were available. '
                     r'Even unmarked fits remain exploratory and omit the covariance sources noted above.' + '\n')
    else:
        parts.append(r'\textbf{No machine-readable LL fit audit was available.} The fit pages are '
                     r'illustrative only; numerical interpretation requires inspection of the ROOT logs.' + '\n')
    parts.extend([
        r'{\small\textbf{Reproduction records:} Measurements: \path{inspect.json}; '
        r'macros/settings: \path{macro_manifest.json}; input/configuration: \path{run_state.json}; '
        r'fit diagnostics: \path{logs/}; figure hashes/pages: \path{report_metadata.json}. '
        r'Commands and full figure index: \path{README.md}.}' + '\n',
        r'\end{document}' + '\n',
    ])
    return ''.join(parts)


def build_readme(energy, info, manifest, pdfs):
    sample = SAMPLES[energy]
    purity = info['inclusivePurityProxy']
    fit = info['gausPol2Comparison']
    root = Path(info['inputRoot'])
    settings = manifest['settings']
    kstar = next(h for h in info['histograms']
                 if h['name'] == 'true/hKstar_4He_CentBin0' and h['present'])
    bin_width = (kstar['xMax'] - kstar['xMin']) / kstar['xBins'] * settings['rebin_group']
    m5 = method5_grid(settings,bin_width)
    fit_bins = (settings['ll_fit']['range'][1] - settings['ll_fit']['range'][0]) / bin_width
    lines = ['<!-- Generated by %s -->' % OWNER,
             '# KFParticle Λ–⁴He summary（%s GeV）' % sample['beam'], '',
             '- [まとめレポート](%s.pdf) / [TeX](%s.tex)' % (sample['stem'], sample['stem']),
             '- 入力 ROOT: [`%s`](../%s)' % (root.name, root.name),
             '- %d P24iy の部分統合: %s / %s subjob。欠損 %s 件を含みません。' %
             (sample['year'], count(sample['included']), count(sample['expected']),
              count(sample['expected'] - sample['included'])),
             '- 表記は fixed-target のビームエネルギーのラベルです。√sNN と同一視しません。', '',
             '## 図の一覧', '',
             '| ファイル | 内容 | ページ数 |', '|---|---|---:|']
    study=root.parent/'pdf'/'mass_halves'
    if (study/'README.md').is_file() and (study/'mass_halves_comparison.pdf').is_file():
        lines[7:7]=['', '- 質量窓の追加比較: [low/high と ±3σ の説明](mass_halves/README.md) / '
                    '[比較PDF](mass_halves/mass_halves_comparison.pdf)。BG subtraction方式ごとの結果をまとめています。', '']
    descriptions = {
        'Lambda_InvMass.pdf': 'inclusive Λ invariant mass、質量ゲート、純度の診断',
        'kstar_sum_4He.pdf': '中央度を合算した SE/ME、sideband、raw SB subtraction CF',
        'kstar_sum_CFsub_4He.pdf': '中央度を合算した SE/ME と一定純度の CF subtraction',
        'kstar_sum_CFsub_method3_4He.pdf': 'Method 3：k* bin ごとの質量 fit による信号数・純度・CF',
        'kstar_sum_CFsub_method5_4He.pdf': 'Method 5：滑動窓の質量 fit による信号数・純度・CF',
        'kstar_cent_CFsub_method3_4He.pdf': 'Method 3：9 中央度の信号数（p1）・純度（p2）・CF（p3）',
        'kstar_cent_CFsub_method5_4He.pdf': 'Method 5：9 中央度の信号数（p1）・純度（p2）・CF（p3）',
        'Fit_LL_4He.pdf': 'raw SB subtraction CF に対する探索的 LL fit',
        'Fit_LL_CFsub_4He.pdf': '一定純度 CF subtraction に対する探索的 LL fit',
        'Fit_LL_method3_4He.pdf': 'Method 3 CF に対する探索的 LL fit',
        'Fit_LL_method5_4He.pdf': 'Method 5 CF に対する探索的 LL fit',
        'kstar_4He.pdf': '9 中央度を別々に描いた raw SB subtraction の各段階',
        'kstar_CFsub_4He.pdf': '9 中央度を別々に描いた一定純度 CF subtraction の各段階',
    }
    for name, record in sorted(pdfs.items()):
        lines.append('| [%s](%s) | %s | %s |' %
                     (name, name, descriptions.get(name, 'マクロ出力'), record['pages']))
    lines.extend(['', '## 計算条件と注意点', '',
        '- production の解析や Λ 選択 cut は変更せず、統合 ROOT を読み取り専用で後処理しています。旧 tools/summary*.pdf は上書きしません。',
        '- Maker YAML の ±3σ signal gate、左右 3–12σ sideband を継承。raw pair sideband 合計の幅補正は 1/3 です。',
        '- inclusive Λ 質量ヒストから求めた純度 proxy は **%s ± %s**（統計誤差のみ）。旧解析の固定 S/N は使用しません。' %
        (number(purity['purity']), number(purity['purityError'])),
        '- inclusive 純度は bin 全体を積分するため、Maker の未ビン化ゲートを近似します。実際の bin 幅比 %s と全境界は inspect.json に保存しています。' % number(purity['actualWidthScale']),
        '- 左右 sideband はそれぞれ %s / %s candidates で非対称です。単純な平坦背景 proxy を精密な純度と同一視しないでください。' %
        (count(info['massWindow']['left']['count']), count(info['massWindow']['right']['count'])),
        '- inclusive gaus+pol2 の χ²/NDF=%s/%s、fit probability=%s。数値収束してもモデル適合性は保証されず、このfit値はCF補正に採用しません。' %
        (number(fit['chi2']), fit['ndf'], number(fit.get('chi2Probability'))),
        '- 純度 proxy は平坦な背景補間を仮定し、ペア・k*・中央度ごとの純度ではありません。この proxy の誤差、規格化の誤差、背景形状の系統誤差は一定純度 CFsub の誤差には入っていません。',
        '- k* rebin=%d（%s MeV/c）。表示範囲は %s–%s GeV/c、SE/ME 規格化範囲は %s–%s GeV/c のままです。表示範囲外の規格化 bin も計算に使用します。Maker の CF 設定とは区別してください。' %
        (settings['rebin_group'], number(bin_width*1000), number(settings['plot_x_range'][0]),
         number(settings['plot_x_range'][1]), number(settings['normalization'][0]), number(settings['normalization'][1])),
        '- 合算中央度を全方式で 20–80% / 0–20%、および 0–10% / 10–20% / 20–60% に統一。raw histogram を合算してから CF を計算しています。',
        '- Method 5 の純度推定は %s MeV/c 刻み、W=%d（左右各 W bin）、内部の fit 窓幅 %s MeV/c。端の窓は短くなります。%s 滑動窓は隣接点に相関を作りますが、従来の独立 bin 誤差では隣接fit間の共分散を含みません。' %
        (number(m5['width']*1000), settings['sliding_half_width'], number(m5['window_width']*1000),
         'CF出力binとは独立です。元の5 MeV binに所属する50 MeV cellの純度を掛けてから20 MeVへ合算します。50/20 MeVの境界が揃わなくても両側の寄与を保持し、bin中心への一括割当や補間は行いません。同一CF binで共有する純度誤差はまとめて伝播します。' if m5['independent'] else
         'この保存版manifestは、純度とCFが同じbin幅の旧方式です。'),
        '- LL fit は a0=%s fm・reff=%s fm 固定（R/Norm free）と全パラメータ free の探索的比較です。fit 範囲に約 %s bin（無効 bin の除外前）があること、失敗・境界到達・共分散の状態に注意し、R や方式の優劣を確定しないでください。' %
        (number(settings['ll_fit']['fixed_a0']), number(settings['ll_fit']['fixed_reff']), number(fit_bins)),
        '- 3.85 GeV では 2021 P24iy のみを使用し、2019 年や pilot ROOT を混ぜていません。', '',
        '## 再現・記録', '',
        'ROOT 図を含めて作り直す場合は新しい出力ディレクトリを指定し、リポジトリ直下から実行します。主設定は対応する mainconf を明示してください。', '',
        '```bash',
        'python3 tools/kfparticle_summary/run_summary.py \\',
        '  --energy %s \\' % energy,
        '  --mainconf config/mainconf/main_auau%s_anaFemtoLambda_4He_KFParticle.yaml \\' % energy,
        '  --input-root %s \\' % str(root),
        '  --output-dir /path/to/new_run/pdf',
        '```', '',
        '完成済みの図からレポートだけを再生成する場合:', '',
        '```bash',
        'python3 tools/kfparticle_summary/build_report.py \\',
        '  --energy %s --output-dir %s' % (energy, str(root.parent / 'pdf')),
        '```', '',
        '- `inspect.json`: 入力 ROOT、event counter、ヒストの型・軸、実測純度、質量 fit の診断。',
        '- `macro_manifest.json`: 使用マクロ、元ファイル・生成ファイルの hash、互換性変更、後処理条件。',
        '- `run_state.json`: mainconf と関連 YAML の hash、入力 fingerprint、実行状態。',
        '- `macros/`、`logs/`: 実行したマクロと ROOT/LaTeX ログ。宣言された TSV に fit 状態と診断を保存。',
        '- `report_metadata.json`: このレポートに使った PDF のページ数・hash、JSON/監査 TSV の hash。',
        '- 部分統合・欠損 subjob の経緯は [親 README](../README.md) を参照してください。', ''])
    return '\n'.join(lines)


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output-dir', required=True, type=Path)
    parser.add_argument('--energy', required=True, choices=sorted(SAMPLES))
    args = parser.parse_args(argv)
    out = args.output_dir.resolve(strict=True)
    info_path, manifest_path = out / 'inspect.json', out / 'macro_manifest.json'
    info = json.loads(info_path.read_text())
    manifest = json.loads(manifest_path.read_text())
    if manifest['energy'] != args.energy:
        raise ValueError('Manifest energy does not match --energy.')
    if Path(manifest['input_root']).resolve() != Path(info['inputRoot']).resolve():
        raise ValueError('Inspection and plotting used different ROOT inputs.')
    expected_token = '%dP24iy_partial_merge.root' % SAMPLES[args.energy]['year']
    if not Path(info['inputRoot']).name.endswith(expected_token):
        raise ValueError('This production-specific report requires the documented partial-merge input.')
    if not info['inclusivePurityProxy']['usableForCFsub'] or not 0 < info['purity_for_cfsub'] < 1:
        raise ValueError('Inclusive purity was not accepted for CFsub.')
    if abs(float(manifest['purity']) - info['purity_for_cfsub']) > 1e-12:
        raise ValueError('Macro purity differs from the inspected current-input purity.')
    pdfs = dict(REQUIRED_PDFS)
    for macro in manifest.get('macros', []):
        for item in macro.get('pdfs', []):
            pdfs[item['name']] = item['pages']
    checked = {}
    for name, pages in sorted(pdfs.items()):
        path = out / name
        found = pdf_pages(path)
        if found != pages:
            raise ValueError('%s: expected %s pages, found %s' % (path, pages, found))
        checked[name] = {'pages': found, 'sha256': sha256(path), 'bytes': path.stat().st_size}
    ll_rows, audit_sources = audits(out, manifest)
    sample = SAMPLES[args.energy]
    metadata_path = out / 'report_metadata.json'
    signature = {'owner': OWNER, 'energy': args.energy, 'inputRoot': info['inputRoot']}
    if metadata_path.exists():
        old = json.loads(metadata_path.read_text())
        if any(old.get(key) != value for key, value in signature.items()):
            raise RuntimeError('Existing report belongs to a different input or generator.')
    else:
        conflicts = [p.name for p in (out / (sample['stem'] + '.tex'),
                                     out / (sample['stem'] + '.pdf'), out / 'README.md') if p.exists()]
        if conflicts:
            raise RuntimeError('Refusing to overwrite unowned report products: ' + ', '.join(conflicts))
    metadata = dict(signature, schemaVersion=1, phase='building', sample=sample,
                    inspect_sha256=sha256(info_path), manifest_sha256=sha256(manifest_path),
                    figures=checked, fit_audits=audit_sources)
    metadata_path.write_text(json.dumps(metadata, indent=2) + '\n')
    tex = out / (sample['stem'] + '.tex')
    tex.write_text(build_tex(args.energy, info, manifest, ll_rows, audit_sources))
    (out / 'README.md').write_text(build_readme(args.energy, info, manifest, checked))
    latex = shutil.which('pdflatex')
    if not latex:
        raise RuntimeError('pdflatex is required.')
    (out / 'logs').mkdir(exist_ok=True)
    for iteration in (1, 2):
        log = out / 'logs' / ('report_pdflatex_%d.log' % iteration)
        with log.open('w') as handle:
            subprocess.run([latex, '-no-shell-escape', '-interaction=nonstopmode',
                            '-halt-on-error', '-file-line-error', tex.name], cwd=str(out),
                           stdout=handle, stderr=subprocess.STDOUT, timeout=120, check=True)
    pdf = out / (sample['stem'] + '.pdf')
    metadata.update(phase='complete', report={'name': pdf.name, 'pages': pdf_pages(pdf),
                    'sha256': sha256(pdf), 'bytes': pdf.stat().st_size},
                    tex_sha256=sha256(tex), generator_sha256=sha256(Path(__file__)))
    metadata_path.write_text(json.dumps(metadata, indent=2) + '\n')
    print('Created %s (%s pages)' % (pdf, metadata['report']['pages']))


if __name__ == '__main__':
    try:
        main()
    except (KeyError, ValueError, RuntimeError, OSError, subprocess.SubprocessError) as exc:
        print('build_report: %s' % exc, file=sys.stderr)
        sys.exit(1)
