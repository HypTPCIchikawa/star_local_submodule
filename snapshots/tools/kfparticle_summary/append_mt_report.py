#!/usr/bin/env python3
"""Append a verified fixed-mT report to an archived original, without reanalysis.

The two component reports remain unmodified. A retained build directory and a
small publication journal make retries safe, including an interrupted publish.
Only the top-level summary PDF, its TeX/README, and completion records change.
"""
import argparse
import copy
from datetime import datetime, timezone
import fcntl
import hashlib
import json
import math
import os
from pathlib import Path
import re
import shutil
import subprocess
import tempfile


DEFAULT_REPORT_STEM = 'cf_mt_cent0_10_constant_vs_method5'
PDF = DEFAULT_REPORT_STEM + '.pdf'
TEX = DEFAULT_REPORT_STEM + '.tex'
RECORD = 'combined_report.json'
OWNER = 'append_mt_report.py/v1'
MARKER = '<!-- BEGIN fixed-mt-appendix -->'
SPLIT_KEYS = {'split_strategy', 'mt_extent', 'fixed_mt_edges'}


def report_names(stem):
    if not isinstance(stem, str) or not re.fullmatch(r'[A-Za-z0-9_-]+', stem):
        raise ValueError('Report stem must be a simple filename without an extension.')
    return stem + '.pdf', stem + '.tex'


def digest(path):
    value = hashlib.sha256()
    with Path(path).open('rb') as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b''):
            value.update(block)
    return value.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text(encoding='utf-8'))


def json_bytes(value):
    return (json.dumps(value, indent=2, ensure_ascii=False, sort_keys=True) + '\n').encode('utf-8')


def atomic_write(path, content):
    """Replace one regular output atomically, using a temporary sibling."""
    path = Path(path)
    if path.is_symlink():
        raise ValueError('Refusing a symbolic-link output: %s' % path)
    fd, temporary = tempfile.mkstemp(prefix='.%s.' % path.name, dir=str(path.parent))
    try:
        with os.fdopen(fd, 'wb') as stream:
            stream.write(content)
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, str(path))
    finally:
        if os.path.exists(temporary):
            os.unlink(temporary)


def regular(path):
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError('Expected a regular file: %s' % path)
    return path


def product_path(directory, name):
    name = Path(name)
    if name.is_absolute() or '..' in name.parts:
        raise ValueError('Unsafe product name: %s' % name)
    candidate = directory / name
    regular(candidate)
    if directory not in candidate.resolve().parents:
        raise ValueError('Product escapes the study directory: %s' % candidate)
    return candidate


def check_products(directory, completed, names=None):
    products = completed.get('products', {})
    for name in names if names is not None else products:
        expected = products.get(name)
        if not expected or digest(product_path(directory, name)) != expected:
            raise ValueError('Product checksum mismatch: %s/%s' % (directory, name))


def pdf_pages(path):
    """Count actual PDF pages; do not trust completion metadata alone."""
    path = Path(path).resolve(strict=True)
    if shutil.which('pdfinfo'):
        output = subprocess.check_output(['pdfinfo', str(path)], text=True, timeout=60)
        found = re.search(r'^Pages:\s+(\d+)\s*$', output, re.M)
        if found:
            return int(found.group(1))
    gs = shutil.which('gs')
    if not gs:
        raise RuntimeError('PDF page verification requires pdfinfo or Ghostscript.')
    result = subprocess.run([
        gs, '-q', '-dNODISPLAY', '-dSAFER', '--permit-file-read=' + str(path),
        '-sPDFFile=' + str(path), '-c',
        'PDFFile (r) file runpdfbegin pdfpagecount = quit'],
        stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True,
        timeout=60, check=True)
    lines = [line.strip() for line in result.stdout.splitlines() if line.strip().isdigit()]
    if len(lines) != 1 or int(lines[0]) < 1:
        raise ValueError('Could not determine PDF page count: %s' % path)
    return int(lines[0])


def load_component(directory, archived=False, report_stem=DEFAULT_REPORT_STEM):
    pdf, tex = report_names(report_stem)
    completed_path = regular(directory / 'completed.json')
    completed = read_json(completed_path)
    if completed.get('phase') != 'complete' or completed.get('input_unchanged') is not True:
        raise ValueError('Component is not a completed, unchanged-input study: %s' % directory)
    # The archive intentionally holds report/provenance files, not the large
    # extraction and all plots. Its missing products remain in the main study.
    required = [pdf, 'README.md', 'context.json']
    check_products(directory, completed, required if archived else None)
    for name in required:
        if name not in completed.get('products', {}):
            raise ValueError('Missing component product: %s' % name)
    regular(directory / tex)
    pages = pdf_pages(directory / pdf)
    if pages != completed.get('pages'):
        raise ValueError('Component page count differs from completed.json: %s' % directory)
    return {
        'directory': str(directory), 'pdf': pdf,
        'pdf_sha256': digest(directory / pdf), 'pages': pages,
        'completed_sha256': digest(completed_path),
        'context_sha256': digest(directory / 'context.json'),
        'readme_sha256': digest(directory / 'README.md'),
        'tex_sha256': digest(directory / tex),
    }, completed, read_json(directory / 'context.json')


def check_compatible(original, appendix, original_done, appendix_done):
    # A reused extraction adds only a provenance record, not a physics setting.
    extra = {'settings', 'settings_sha256', 'extraction_cache'}
    original_common = {k: v for k, v in original.items() if k not in extra}
    appendix_common = {k: v for k, v in appendix.items() if k not in extra}
    if original_common != appendix_common:
        raise ValueError('Input, purity source, or common analysis settings differ.')
    if original.get('energy') != '3p85' or original.get('shared_method5_approved') is not True:
        raise ValueError('Expected the approved 3.85 GeV shared-purity study.')
    left, right = original['settings'], appendix['settings']
    if {k: v for k, v in left.items() if k not in SPLIT_KEYS} != {
            k: v for k, v in right.items() if k not in SPLIT_KEYS}:
        raise ValueError('Non-slicing settings differ between reports.')
    if (left.get('split_strategy') != 'nearest_se_signal_quantiles_at_saved_mt_edges'
            or right.get('split_strategy') != 'fixed_saved_mt_edges'):
        raise ValueError('Expected original quartiles followed by fixed saved mT edges.')
    edges = right.get('fixed_mt_edges', [])
    if (len(edges) != right['n_slices'] + 1 or
            not all(isinstance(x, (int, float)) and math.isfinite(x) for x in edges) or
            not all(a < b for a, b in zip(edges, edges[1:]))):
        raise ValueError('Invalid fixed mT edges.')
    extraction_hash = original_done['products'].get('mt_extract.json')
    if not extraction_hash or appendix_done['products'].get('mt_extract.json') != extraction_hash:
        raise ValueError('The two reports must use the same verified histogram extraction.')
    cache = appendix.get('extraction_cache')
    if cache and cache.get('sha256') != extraction_hash:
        raise ValueError('Appendix extraction-cache provenance differs.')


def tex_document(original_path, appendix_path):
    # \detokenize supports underscores/spaces in the established report paths.
    paths = [str(original_path), str(appendix_path)]
    if any(any(char in value for char in '{}%\n\r') for value in paths):
        raise ValueError('Unsupported TeX path characters.')
    return ('\\documentclass{article}\n\\usepackage{pdfpages}\n'
            '\\pagestyle{empty}\n\\begin{document}\n' +
            ''.join('\\includepdf[pages=-,fitpaper=true,pagecommand={}]{\\detokenize{%s}}\n' %
                    value for value in paths) +
            '\\end{document}\n')


def build_pdf(build, original, appendix):
    """Build from byte-verified copies, with the original component first."""
    for name, component in [('original.pdf', original), ('appendix.pdf', appendix)]:
        target = build / name
        shutil.copy2(str(Path(component['directory']) / component['pdf']), str(target))
        if digest(target) != component['pdf_sha256']:
            raise ValueError('Component changed while being copied.')
    (build / 'combined.tex').write_text(tex_document('original.pdf', 'appendix.pdf'),
                                        encoding='utf-8')
    pdfunite = shutil.which('pdfunite')
    if pdfunite:
        command = [pdfunite, 'original.pdf', 'appendix.pdf', 'combined.pdf']
        method = 'pdfunite (component order: original, appendix)'
    else:
        latex = shutil.which('pdflatex')
        if not latex:
            raise RuntimeError('PDF concatenation requires pdfunite or pdflatex with pdfpages.')
        command = [latex, '-no-shell-escape', '-interaction=nonstopmode',
                   '-halt-on-error', 'combined.tex']
        method = 'pdflatex/pdfpages, pages=-, fitpaper=true (original, appendix)'
    with (build / 'build.log').open('w', encoding='utf-8') as log:
        subprocess.run(command, cwd=str(build), stdout=log, stderr=subprocess.STDOUT,
                       timeout=120, check=True)
    output = regular(build / 'combined.pdf')
    if pdf_pages(output) != original['pages'] + appendix['pages']:
        raise ValueError('Combined PDF does not contain the sum of component pages.')
    return output, method


def relative_link(path, study):
    return os.path.relpath(str(path), str(study)).replace(os.sep, '/')


def appended_readme(study, original, appendix, original_pages, appendix_pages, pdf=PDF):
    text = (original / 'README.md').read_text(encoding='utf-8')
    if MARKER in text:
        raise ValueError('Archived original README already contains an appendix marker.')
    last = original_pages + appendix_pages
    note = [
        '', MARKER, '## 固定幅 m_T 分割の追補', '',
        '上記の表・個別図は元の統計四分位分割を保持しています。'
        '同名の[まとめPDF](%s)には、元の%dページを変更せず先頭に残し、'
        '固定幅分割の%dページを追補しました。' % (pdf, original_pages, appendix_pages),
        '',
        '- 1–%dページ: 元の統計四分位分割。' % original_pages,
        '- %d–%dページ: [固定幅分割の比較](%s)。'
        '設定、各区間の統計、範囲外の除外数は[追補の説明](%s)を参照してください。' %
        (original_pages + 1, last, relative_link(appendix / pdf, study),
         relative_link(appendix / 'README.md', study)),
        '- [追補前のPDF](%s)と[追補前の説明](%s)を別途保存しています。' %
        (relative_link(original / pdf, study), relative_link(original / 'README.md', study)),
        '- 入力ROOT、元の抽出データ、元の個別図、元の解析設定は変更していません。'
        '固定幅分割と元の分割では m_T の選択範囲・区間内統計が異なります。',
        '- 結合の構成順、チェックサム、ページ数は[結合記録](combined_report.json)に記録しています。',
        '<!-- END fixed-mt-appendix -->', '',
    ]
    return text.rstrip() + '\n' + '\n'.join(note)


def same_component(recorded, current):
    return recorded == current


def verify_unchanged_products(study, original_done, changed):
    names = [name for name in original_done['products'] if name not in changed]
    check_products(study, original_done, names)


def finish_publication(study, record, original_done):
    """Resume only this journal's old-or-planned bytes; reject foreign changes."""
    pdf = record['original']['pdf']
    tex = str(Path(pdf).with_suffix('.tex'))
    changed = {pdf, tex, 'README.md'}
    verify_unchanged_products(study, original_done, changed)
    for name, item in record['outputs'].items():
        current = digest(regular(study / name))
        if current not in (item['before_sha256'], item['sha256']):
            raise ValueError('Output changed outside this publication: %s' % name)
        staged = product_path(study, item['staged'])
        if digest(staged) != item['sha256']:
            raise ValueError('Staged output checksum mismatch: %s' % staged)
    # Do not leave a complete flag over an in-progress, mixed product set.
    pending = copy.deepcopy(original_done)
    pending.update(phase='report_publication_pending', combined_report=RECORD)
    atomic_write(study / 'completed.json', json_bytes(pending))
    for name, item in record['outputs'].items():
        if digest(study / name) != item['sha256']:
            atomic_write(study / name, (study / item['staged']).read_bytes())
    if pdf_pages(study / pdf) != record['pages']:
        raise ValueError('Published PDF page-count verification failed.')
    verify_unchanged_products(study, original_done, changed)
    record['phase'] = 'complete'
    atomic_write(study / RECORD, json_bytes(record))
    completed = copy.deepcopy(original_done)
    completed.update(phase='complete', pages=record['pages'], combined_report=RECORD)
    completed['sources'][record['source']['path']] = record['source']['sha256']
    for name, item in record['outputs'].items():
        completed['products'][name] = item['sha256']
    completed['products'][RECORD] = digest(study / RECORD)
    atomic_write(study / 'completed.json', json_bytes(completed))
    check_products(study, completed)
    return record


def append_report(study_dir, original_dir, appendix_dir, report_stem=DEFAULT_REPORT_STEM):
    pdf, tex = report_names(report_stem)
    study, original, appendix = [Path(p).resolve(strict=True) for p in
                                  (study_dir, original_dir, appendix_dir)]
    if len({study, original, appendix}) != 3 or not all(p.is_dir() for p in (study, original, appendix)):
        raise ValueError('Use three distinct existing report directories.')
    if study not in original.parents or study not in appendix.parents:
        raise ValueError('Original archive and appendix must be below the top-level study.')
    with (study / '.append_mt_report.lock').open('a') as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        orig, original_done, original_context = load_component(
            original, archived=True, report_stem=report_stem)
        added, appendix_done, appendix_context = load_component(appendix, report_stem=report_stem)
        check_compatible(original_context, appendix_context, original_done, appendix_done)
        source = {'path': str(Path(__file__).resolve()), 'sha256': digest(__file__)}
        journal = study / RECORD
        if journal.exists():
            record = read_json(regular(journal))
            if (record.get('owner') != OWNER or record.get('source') != source or
                    record.get('study_directory') != str(study) or
                    not same_component(record.get('original'), orig) or
                    not same_component(record.get('appendix'), added) or
                    record.get('phase') not in ('prepared', 'complete')):
                raise ValueError('Existing combination belongs to different components or source.')
            if record['phase'] == 'complete':
                done = read_json(study / 'completed.json')
                if done.get('phase') == 'complete' and done.get('combined_report') == RECORD:
                    check_products(study, done)
                    if (done['products'].get(pdf) != record['outputs'][pdf]['sha256'] or
                            done['products'].get(RECORD) != digest(journal) or
                            done.get('pages') != record['pages']):
                        raise ValueError('Combined completion record differs.')
                    verify_unchanged_products(study, original_done, {pdf, tex, 'README.md'})
                    if (done['products'].get('mt_extract.json') != original_done['products']['mt_extract.json']
                            or any(done['sources'].get(k) != v for k, v in original_done['sources'].items())):
                        raise ValueError('Original extraction / source provenance changed.')
                    return record
            return finish_publication(study, record, original_done)
        # First publication is allowed only over the exact verified original.
        if read_json(regular(study / 'completed.json')) != original_done:
            raise ValueError('Top completion record is not the archived original.')
        check_products(study, original_done)
        if digest(regular(study / tex)) != orig['tex_sha256']:
            raise ValueError('Top TeX differs from the archived original.')
        build_parent = study / 'combined_report_build'
        build_parent.mkdir(exist_ok=True)
        build = Path(tempfile.mkdtemp(prefix='append_', dir=str(build_parent)))
        combined, method = build_pdf(build, orig, added)
        # Retain self-contained component PDFs/build source. Top TeX instead
        # references the preserved public components and can be rebuilt there.
        (build / 'published.tex').write_text(tex_document(
            relative_link(original / pdf, study), relative_link(appendix / pdf, study)),
            encoding='utf-8')
        (build / 'published_README.md').write_text(
            appended_readme(study, original, appendix, orig['pages'], added['pages'], pdf=pdf),
            encoding='utf-8')
        staged = {pdf: combined, tex: build / 'published.tex', 'README.md': build / 'published_README.md'}
        record = dict(
            owner=OWNER, phase='prepared', study_directory=str(study),
            created_utc=datetime.now(timezone.utc).isoformat(), source=source,
            original=orig, appendix=added, pages=orig['pages'] + added['pages'],
            page_order=['original', 'appendix'], concatenation=method,
            input=original_context['input'], common_settings_verified=True,
            original_context_unchanged=True, original_plots_unchanged=True,
            build_directory=relative_link(build, study),
            outputs={name: dict(staged=relative_link(path, study), sha256=digest(path),
                                before_sha256=digest(study / name)) for name, path in staged.items()})
        # Recheck the original PDF after the potentially slow external build.
        check_products(study, original_done)
        if digest(study / tex) != orig['tex_sha256']:
            raise ValueError('Original TeX changed while combining.')
        atomic_write(journal, json_bytes(record))
        return finish_publication(study, record, original_done)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--study-dir', type=Path, required=True)
    parser.add_argument('--original-dir', type=Path, required=True)
    parser.add_argument('--appendix-dir', type=Path, required=True)
    parser.add_argument('--report-stem', default=DEFAULT_REPORT_STEM,
                        help='Summary filename without extension (default: %(default)s).')
    args = parser.parse_args()
    result = append_report(args.study_dir, args.original_dir, args.appendix_dir, args.report_stem)
    print('COMPLETE: %d pages, original followed by fixed-mT appendix: %s' %
          (result['pages'], args.study_dir / result['original']['pdf']))


if __name__ == '__main__':
    main()
