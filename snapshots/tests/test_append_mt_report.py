import copy
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest import mock

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'tools/kfparticle_summary'))
import append_mt_report as report


class AppendMtReportTests(unittest.TestCase):
    def setUp(self):
        self.temporary = tempfile.TemporaryDirectory()
        self.addCleanup(self.temporary.cleanup)
        self.study = Path(self.temporary.name)
        self.original = self.study / 'revisions/original'
        self.appendix = self.study / 'fixed_02'
        self.original.mkdir(parents=True)
        self.appendix.mkdir()
        context = dict(input={'path': 'unchanged.root', 'size': 10},
                       base={'normalization': [.3, .4]}, energy='3p85',
                       shared_method5_approved=True,
                       settings=dict(n_slices=4, cent9_min=7, cent9_max=8,
                           split_strategy='nearest_se_signal_quantiles_at_saved_mt_edges',
                           mt_extent='occupied_regular_bins_all_channels'))
        self.source = '/saved/ExtractMtCF.C'
        for directory, label in [(self.study, 'original'), (self.original, 'original'),
                                 (self.appendix, 'appendix')]:
            ctx = copy.deepcopy(context)
            if label == 'appendix':
                ctx['settings'].update(split_strategy='fixed_saved_mt_edges',
                                       mt_extent='configured_range',
                                       fixed_mt_edges=[2.4, 2.6, 2.8, 3., 3.2])
            for name, content in {
                report.PDF: label, report.TEX: label + ' tex',
                'README.md': label + ' description', 'mt_extract.json': '{"same": true}',
                'context.json': json.dumps(ctx), 'mt_distribution.pdf': 'old plot',
            }.items():
                (directory / name).write_text(content)
            done = dict(phase='complete', input_unchanged=True, pages=4,
                        sources={self.source: 'extractor-sha'},
                        products={p.name: report.digest(p) for p in directory.iterdir()
                                  if p.is_file() and p.name != report.TEX})
            (directory / 'completed.json').write_bytes(report.json_bytes(done))
        self.before = report.read_json(self.study / 'completed.json')
        self.build_calls = 0
        self.pages = mock.patch.object(report, 'pdf_pages', side_effect=self.fake_pages)
        self.pages.start()
        self.addCleanup(self.pages.stop)
        self.build = mock.patch.object(report, 'build_pdf', side_effect=self.fake_build)
        self.build.start()
        self.addCleanup(self.build.stop)

    def fake_pages(self, path):
        return 8 if Path(path).read_text() == 'original+appendix' else 4

    def fake_build(self, build, original, appendix):
        self.build_calls += 1
        self.assertEqual(Path(original['directory']), self.original)
        self.assertEqual(Path(appendix['directory']), self.appendix)
        out = build / 'combined.pdf'
        out.write_text('original+appendix')
        return out, 'fixture concat'

    def run_append(self):
        return report.append_report(self.study, self.original, self.appendix)

    def test_publish_is_idempotent_and_preserves_extraction_and_plots(self):
        first = self.run_append()
        second = self.run_append()
        self.assertEqual(first, second)
        self.assertEqual(self.build_calls, 1)
        self.assertEqual(first['pages'], 8)
        self.assertEqual((self.study / report.PDF).read_text(), 'original+appendix')
        self.assertEqual((self.original / report.PDF).read_text(), 'original')
        self.assertEqual((self.study / 'README.md').read_text().count(report.MARKER), 1)
        done = report.read_json(self.study / 'completed.json')
        self.assertEqual(done['sources'][self.source], 'extractor-sha')
        for name in ('mt_extract.json', 'context.json', 'mt_distribution.pdf'):
            self.assertEqual(done['products'][name], self.before['products'][name])
        report.check_products(self.study, done)

    def test_rejects_modified_original_target(self):
        (self.study / report.PDF).write_text('unrelated report')
        with self.assertRaisesRegex(ValueError, 'checksum'):
            self.run_append()
        self.assertFalse((self.study / report.RECORD).exists())

    def test_alternative_report_stem_and_idempotent_retry(self):
        stem = 'cf_mt_cent0_20_constant_vs_method5'
        pdf, tex = report.report_names(stem)
        for directory in (self.study, self.original, self.appendix):
            (directory / report.PDF).rename(directory / pdf)
            (directory / report.TEX).rename(directory / tex)
            context_path = directory / 'context.json'
            context = report.read_json(context_path)
            context['settings']['cent9_min'] = 6
            context_path.write_bytes(report.json_bytes(context))
            done = report.read_json(directory / 'completed.json')
            done['products'][pdf] = done['products'].pop(report.PDF)
            done['products']['context.json'] = report.digest(context_path)
            (directory / 'completed.json').write_bytes(report.json_bytes(done))
        first = report.append_report(self.study, self.original, self.appendix, report_stem=stem)
        second = report.append_report(self.study, self.original, self.appendix, report_stem=stem)
        self.assertEqual(first, second)
        self.assertEqual(first['original']['pdf'], pdf)
        self.assertEqual(self.build_calls, 1)
        done = report.read_json(self.study / 'completed.json')
        report.check_products(self.study, done)
        self.assertIn(pdf, done['products'])
        self.assertIn(tex, done['products'])
        self.assertNotIn(report.PDF, done['products'])
        self.assertFalse((self.study / report.PDF).exists())
        self.assertIn('fixed_02/' + pdf, (self.study / 'README.md').read_text())
        self.assertIn('fixed_02/' + pdf, (self.study / tex).read_text())

    def test_rejects_unsafe_report_stem(self):
        for stem in ('../report', '/report', 'report.pdf', ''):
            with self.assertRaisesRegex(ValueError, 'Report stem'):
                report.append_report(self.study, self.original, self.appendix, report_stem=stem)
        self.assertFalse((self.study / report.RECORD).exists())

    def test_rejects_different_normalization(self):
        context_path = self.appendix / 'context.json'
        context = report.read_json(context_path)
        context['base']['normalization'] = [.4, .6]
        context_path.write_bytes(report.json_bytes(context))
        done = report.read_json(self.appendix / 'completed.json')
        done['products']['context.json'] = report.digest(context_path)
        (self.appendix / 'completed.json').write_bytes(report.json_bytes(done))
        with self.assertRaisesRegex(ValueError, 'common analysis settings'):
            self.run_append()

    def test_interrupted_publication_resumes_without_duplicate_appendix(self):
        original_write = report.atomic_write
        interrupted = [False]

        def flaky(path, content):
            if Path(path) == self.study / report.TEX and not interrupted[0]:
                interrupted[0] = True
                raise OSError('simulated disconnect')
            original_write(path, content)

        with mock.patch.object(report, 'atomic_write', side_effect=flaky):
            with self.assertRaisesRegex(OSError, 'simulated'):
                self.run_append()
        self.assertEqual(report.read_json(self.study / 'completed.json')['phase'],
                         'report_publication_pending')
        self.run_append()
        self.assertEqual(self.build_calls, 1)
        self.assertEqual((self.study / 'README.md').read_text().count(report.MARKER), 1)
        self.assertEqual(report.read_json(self.study / 'completed.json')['phase'], 'complete')


if __name__ == '__main__':
    unittest.main()
