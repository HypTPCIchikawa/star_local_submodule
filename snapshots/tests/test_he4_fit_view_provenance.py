"""The stored-fit renderer must never silently accept changed runtime config."""
import importlib.util
import json
from pathlib import Path
import sys
import tempfile
import unittest
from unittest.mock import patch

import yaml

PROJECT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT/'tools/he4_purity'))
spec = importlib.util.spec_from_file_location('he4_fit_views', PROJECT/'tools/he4_purity/render_fit_views.py')
views = importlib.util.module_from_spec(spec)
spec.loader.exec_module(views)


class ViewProvenanceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory(prefix='he4-view-provenance-')
        self.project = Path(self.temp.name)
        (self.project/'config/mainconf').mkdir(parents=True)
        (self.project/'config/analysis').mkdir()
        self.main = self.project/'config/mainconf/current.yaml'
        self.info = self.project/'config/analysis/runtime.yaml'
        self.main.write_text('analysis: analysis/runtime.yaml\n')
        self.identity = dict(starTag=dict(libraryTag='SL24y', year='year2021',
            productionTag='P24iy', collision='auau3p85',
            triggerSets='production_3p85GeV_fixedTarget_2021'), analysis=dict(mode='fxtmult'))
        self.info.write_text(yaml.safe_dump(self.identity))
        self.context = dict(mainconf=str(self.main), mainconf_sha256=views.digest(self.main),
            sample='3p85', input=dict(path='saved_2021P24iy.root'),
            settings=dict(samples={'3p85':dict(input='saved_2021P24iy.root')}))
        self.project_patch = patch.object(views, 'PROJECT', self.project)
        self.project_patch.start()

    def tearDown(self):
        self.project_patch.stop()
        self.temp.cleanup()

    def drift(self):
        self.main.write_text(self.main.read_text()+'# Current PID changed, runtime only.\n')

    def test_exact_mainconf_remains_default(self):
        result = views.runtime_provenance(self.context, self.main)
        self.assertEqual(result['mode'], 'saved_mainconf_exact')
        self.assertFalse(result['drift_detected'])

    def test_drift_requires_optin(self):
        self.drift()
        with self.assertRaisesRegex(ValueError, 'explicit opt-in'):
            views.runtime_provenance(self.context, self.main)

    def test_explicit_runtime_only_has_old_and_current_sha(self):
        self.drift()
        result = views.runtime_provenance(self.context, self.main, True)
        self.assertEqual(result['mode'], 'explicit_runtime_only_drift')
        self.assertNotEqual(result['saved_mainconf_sha256'], result['runtime_mainconf_sha256'])
        self.assertEqual(result['runtime_identity']['libraryTag'], 'SL24y')

    def test_optin_does_not_allow_another_mainconf_path(self):
        another = self.main.with_name('another.yaml')
        another.write_text(self.main.read_text())
        with self.assertRaisesRegex(ValueError, 'path differs'):
            views.runtime_provenance(self.context, another, True)

    def test_optin_rejects_non2021_or_wrong_star(self):
        self.drift()
        for key, value in [('year','year2019'),('productionTag','P23ia'),('libraryTag','SL25a')]:
            wrong = json.loads(json.dumps(self.identity))
            wrong['starTag'][key] = value
            self.info.write_text(yaml.safe_dump(wrong))
            with self.assertRaisesRegex(ValueError, 'verified 3p85'):
                views.runtime_provenance(self.context, self.main, True)

    def test_optin_rejects_other_saved_dataset(self):
        self.drift()
        self.context['sample'] = '13p5'
        with self.assertRaisesRegex(ValueError, 'verified 3p85'):
            views.runtime_provenance(self.context, self.main, True)

    def test_previous_view_hashes_are_mandatory(self):
        study = self.project/'study'
        old = study/'views'
        old.mkdir(parents=True)
        pdf = old/'old.pdf'
        pdf.write_bytes(b'frozen PDF fixture')
        (old/'completed.json').write_text(json.dumps(dict(phase='complete',
            products={'old.pdf':views.digest(pdf)})))
        self.assertEqual(views.verified_previous_views(study)[str(pdf)], views.digest(pdf))
        pdf.write_bytes(b'changed')
        with self.assertRaisesRegex(ValueError, 'hash mismatch'):
            views.verified_previous_views(study)


if __name__ == '__main__':
    unittest.main()
