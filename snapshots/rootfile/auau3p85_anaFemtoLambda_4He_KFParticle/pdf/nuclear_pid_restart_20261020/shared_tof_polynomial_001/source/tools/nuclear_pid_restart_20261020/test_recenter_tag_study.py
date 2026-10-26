import copy
import math
import tempfile
from pathlib import Path
import unittest

import recenter_tag_study as study


class RecenterContracts(unittest.TestCase):
    def setUp(self):
        self.axis = dict(dedxBins=20, dedxBinWidth=.1, dedxRange=[0, 2])
        self.settings = dict(searchCenter='original_git_mean', searchHalfWidthSigma=2.5,
                             smoothingWindowBins=3,
                             modeTiePolicy='nearest_original_mean_then_lower_bin',
                             roiWidth='original_git_sigma', centerDefinition='native_bin_center')
        self.spectrum = {i: dict(index=i, lo=(i-1)*.1, hi=i*.1,
                                center=(i-.5)*.1, count=0., error=0.)
                         for i in range(22)}

    def test_mode_finder_uses_three_raw_bins_and_changes_no_counts(self):
        for index, count in ((8, 10.), (9, 20.), (10, 30.), (11, 25.), (12, 5.)):
            self.spectrum[index]['count'] = count
        before = copy.deepcopy(self.spectrum)
        result = study.locate_mode(self.spectrum, 1., .2, self.settings, self.axis)
        self.assertEqual(result['modeBin'], 10)
        self.assertEqual(result['modeThreeBinSum'], 75.)
        self.assertFalse(result['smoothingAppliedToFitCounts'])
        self.assertEqual(before, self.spectrum)

    def test_ties_choose_nearest_original_mean_then_lower_bin(self):
        for index in range(22):
            self.spectrum[index]['count'] = 10.
        result = study.locate_mode(self.spectrum, 1., .2, self.settings, self.axis)
        self.assertEqual(result['modeBin'], 10)
        self.assertGreater(len(result['tiedMaximumBins']), 1)

    def test_mode_search_stays_in_original_window(self):
        self.spectrum[10]['count'] = 10.
        self.spectrum[19]['count'] = 1000000.
        result = study.locate_mode(self.spectrum, 1., .1, self.settings, self.axis)
        self.assertLess(result['modeCenter'], 1.25)
        self.assertLess(result['modeThreeBinSum'], 100.)

    def test_zero_population_cannot_become_fake_mode_or_zero_purity(self):
        with self.assertRaisesRegex(ValueError, 'Empty search'):
            study.locate_mode(self.spectrum, 1., .1, self.settings, self.axis)

    def test_boundary_peak_is_flagged_not_silently_vetoed(self):
        self.spectrum[7]['count'] = 100.
        result = study.locate_mode(self.spectrum, 1., .1, self.settings, self.axis)
        self.assertTrue(result['searchBoundaryMode'])

    def test_only_explicit_three_bin_raw_mode_policy_is_supported(self):
        for key, value in (('smoothingWindowBins', 5), ('searchCenter', 'fitted_mean'),
                           ('roiWidth', 'fitted_sigma'), ('modeTiePolicy', 'first')):
            broken = copy.deepcopy(self.settings)
            broken[key] = value
            with self.subTest(key=key), self.assertRaises(ValueError):
                study.validate_mode_settings(broken)

    def test_roi_uses_original_sigma_and_native_edges(self):
        result = study.recentered_roi(1.05, .12, 2.5, .1, self.axis)
        self.assertAlmostEqual(result['fitLo'], .7)
        self.assertAlmostEqual(result['fitHi'], 1.4)
        narrow = study.recentered_roi(1.05, .12, 1.5, .1, self.axis)
        self.assertAlmostEqual(narrow['fitLo'], .8)
        self.assertAlmostEqual(narrow['fitHi'], 1.3)

    def test_roi_clipping_is_explicit(self):
        result = study.recentered_roi(.15, .2, 2.5, .1, self.axis)
        self.assertEqual(result['fitLo'], 0.)
        self.assertTrue(result['fitWindowAxisClipped'])

    def test_fractional_native_bin_rebinning_cannot_be_silent(self):
        with self.assertRaises(ValueError):
            study.recentered_roi(1., .2, 2.5, .13, self.axis)

    def export(self, directory, omit=None, corrupt=None):
        table = directory / 'native.tsv'
        rows = []
        for i, values in self.spectrum.items():
            if i == omit:
                continue
            count = 4.
            rows.append(['he4_b0088', 'TAG_D', 'he4', 88, i,
                         values['lo'], values['hi'], count, math.sqrt(count)])
        if corrupt:
            corrupt(rows)
        study.base.tsv(table, ['id', 'mode', 'target', 'rBin', 'dedxBin',
                               'lo', 'hi', 'count', 'error'], rows)
        return table

    def test_native_reader_preserves_full_zero_flow_contract(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = self.export(Path(tmp))
            result = study.native_spectra(table, [dict(id='he4_b0088', mode='TAG_D',
                                                     target='he4', rBin=88)], self.axis)
            self.assertEqual(len(result['he4_b0088']), 22)

    def test_missing_flow_bin_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            table = self.export(Path(tmp), omit=0)
            with self.assertRaisesRegex(ValueError, 'Full native spectrum'):
                study.native_spectra(table, [dict(id='he4_b0088', mode='TAG_D',
                                                target='he4', rBin=88)], self.axis)

    def test_wrong_m2_sample_cannot_be_silent(self):
        def corrupt(rows):
            rows[5][1] = 'TAG_HE3'
        with tempfile.TemporaryDirectory() as tmp:
            table = self.export(Path(tmp), corrupt=corrupt)
            with self.assertRaisesRegex(ValueError, 'identity differs'):
                study.native_spectra(table, [dict(id='he4_b0088', mode='TAG_D',
                                                target='he4', rBin=88)], self.axis)


if __name__ == '__main__':
    unittest.main()
