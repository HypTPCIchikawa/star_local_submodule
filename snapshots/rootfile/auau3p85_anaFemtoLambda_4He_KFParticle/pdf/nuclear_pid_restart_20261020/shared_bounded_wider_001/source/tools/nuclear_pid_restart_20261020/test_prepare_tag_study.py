import copy
from pathlib import Path
import unittest

import yaml

import prepare_tag_study as study
from audit_restore import tables_from_header


class FreshTagContracts(unittest.TestCase):
    def test_fit_completion_marker_is_specific_to_the_engine(self):
        self.assertEqual(study.fit_completion_marker('dir/FitStorageSlices.C'), 'RESTART_STORAGE_READBACK_OK')
        self.assertEqual(study.fit_completion_marker('dir/FitSplitNormalSlices.C'), 'RESTART_SPLIT_NORMAL_READBACK_OK')
        self.assertEqual(study.fit_completion_marker('dir/FitRestartSlices.C'), 'RESTART_PEAK_READBACK_OK')

    def setUp(self):
        self.cfg = yaml.safe_load(Path(study.__file__).with_name('tag_study.yaml').read_text())
        self.tables = tables_from_header((study.PROJECT/self.cfg['legacyHeader']).read_text())
        self.axis = dict(rBinWidth=.01, rRange=[0, 6], rBins=600,
                         dedxBinWidth=.1, dedxRange=[0, 120])

    def test_original_polynomial_and_positive_width(self):
        mu, sigma = study.parameters(self.tables, 3, .875)
        self.assertAlmostEqual(mu, 46.44951110046477)
        self.assertAlmostEqual(sigma, 2.679187310317944)

    def test_no_extrapolation_or_constant_width_fallback(self):
        with self.assertRaises(ValueError):
            study.parameters(self.tables, 3, .3)
        broken = copy.deepcopy(self.tables)
        broken['sigma'][3] = [-1, 0, 0, 0]
        with self.assertRaises(ValueError):
            study.parameters(broken, 3, .875)

    def test_he4_and_d_share_same_m2_histogram(self):
        slices, nodes, records = study.make_slices(self.cfg, self.tables, self.axis, True)
        self.assertEqual(len(slices), 8)
        self.assertEqual({row[1] for row in slices if row[2] in ('d', 'he4')}, {'TAG_D'})
        self.assertEqual(len(nodes), len(slices)*12)
        self.assertEqual(len(records), len(slices))

    def test_distinct_he4_tof_window_is_rejected(self):
        self.cfg['tofWindows']['he4'] = [3, 4]
        with self.assertRaises(ValueError):
            study.make_slices(self.cfg, self.tables, self.axis, True)

    def test_wrong_tag_or_all_fit_cannot_be_silent(self):
        self.cfg['samples']['he4'] = 'ALL'
        with self.assertRaises(ValueError):
            study.make_slices(self.cfg, self.tables, self.axis, True)

    def test_exact_native_bins_and_gl3_weights(self):
        slices, nodes, records = study.make_slices(self.cfg, self.tables, self.axis, True)
        self.assertEqual({r['rBin'] for r in records}, {88, 178})
        for row in slices:
            for species in study.SPECIES:
                chosen = [node for node in nodes if node[0] == row[0] and node[1] == species]
                self.assertAlmostEqual(sum(node[2] for node in chosen), 1.)
                self.assertTrue(all(node[3] > 0 and node[4] > 0 for node in chosen))

    def test_non_native_fit_alignment_rejected(self):
        self.cfg['fitEdgeAlignment'] = .13
        with self.assertRaises(ValueError):
            study.make_slices(self.cfg, self.tables, self.axis, True)

    def test_bad_r_bin_or_missing_original_support_rejected(self):
        for bins in ([0], [601], [30], [88, 88]):
            self.cfg['smokeNativeRBins'] = bins
            with self.subTest(bins=bins), self.assertRaises(ValueError):
                study.make_slices(self.cfg, self.tables, self.axis, True)

    def test_no_purity_or_cutoff_is_inferred_from_grid(self):
        self.assertNotIn('cutoff', self.cfg)
        self.assertNotIn('veto', self.cfg)
        self.assertEqual(self.cfg['numeric']['backgroundModel'], 'linear')
        self.assertEqual(study.parameters(self.tables, 2, .875),
                         (33.23063300832934, 1.9034234286354532))


if __name__ == '__main__':
    unittest.main()
