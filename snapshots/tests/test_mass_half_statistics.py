#!/usr/bin/env python3
"""Shared-bin first-order error contracts; standard library and no ROOT I/O."""

import importlib.util
import math
from pathlib import Path
import unittest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    'mass_half_statistics', ROOT / 'tools/kfparticle_summary/mass_half_statistics.py')
STATS = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(STATS)


class MassHalfStatistics(unittest.TestCase):
    def setUp(self):
        self.model = STATS.ErrorModel()
        self.al = self.model.series('AL', [30., 100., 200.], [30., 110., 230.])
        self.ah = self.model.series('AH', [50., 120., 210.], [50., 130., 240.])
        self.bl = self.model.series('BL', [100., 1000., 2000.], [120., 1100., 2300.])
        self.bh = self.model.series('BH', [120., 1200., 2100.], [140., 1300., 2400.])
        self.af = STATS.series_add(self.al, self.ah)
        self.bf = STATS.series_add(self.bl, self.bh)

    def test_identical_quantity_ratio_has_zero_error(self):
        quantity = self.al[0]
        self.assertEqual(self.model.estimate(quantity / quantity),
                         dict(value=1.0, variance=0.0, error=0.0))
        self.assertEqual(self.model.variance(quantity - quantity), 0.0)

    def test_normalization_bin_is_not_independent_duplicate(self):
        # One bin normalizes itself; the resulting CF is exactly one.
        cf = STATS.normalized_cf(self.al, self.bl, [0])['bins'][0]
        self.assertAlmostEqual(cf.value, 1.0)
        self.assertLess(self.model.variance(cf), 1e-28)

    def test_finite_difference_derivatives_with_shared_full_window(self):
        def compute(model_inputs):
            al, ah, bl, bh = model_inputs
            norm = lambda x: x[1] + x[2]
            cf = lambda a, b: norm(b) / norm(a) * a[1] / b[1]
            af = [a + b for a, b in zip(al, ah)]
            bf = [a + b for a, b in zip(bl, bh)]
            return cf(al, bl) / cf(af, bf)
        low = STATS.normalized_cf(self.al, self.bl, [1, 2])['bins'][1]
        full = STATS.normalized_cf(self.af, self.bf, [1, 2])['bins'][1]
        ratio = low / full
        inputs = [[q.value for q in series] for series in (self.al, self.ah, self.bl, self.bh)]
        self.assertAlmostEqual(ratio.value, compute(inputs))
        finite_variance = 0.0
        for series_index, prefix in enumerate(('AL', 'AH', 'BL', 'BH')):
            for bin_index in range(3):
                step = inputs[series_index][bin_index] * 1e-5
                inputs[series_index][bin_index] += step
                plus = compute(inputs)
                inputs[series_index][bin_index] -= 2.0 * step
                minus = compute(inputs)
                inputs[series_index][bin_index] += step
                derivative = (plus - minus) / (2.0 * step)
                key = (prefix, bin_index)
                self.assertAlmostEqual(ratio.gradient.get(key, 0.0), derivative, places=10)
                finite_variance += derivative ** 2 * self.model.variances[key]
        self.assertAlmostEqual(self.model.variance(ratio), finite_variance, places=11)

    def test_low_full_errors_include_covariance(self):
        low = STATS.normalized_cf(self.al, self.bl, [1, 2])['bins'][0]
        full = STATS.normalized_cf(self.af, self.bf, [1, 2])['bins'][0]
        covariance = self.model.covariance(low, full)
        self.assertGreater(covariance, 0.0)
        difference = STATS.contrast(self.model, low, full)
        self.assertAlmostEqual(difference['difference']['variance'],
                               self.model.variance(low) + self.model.variance(full) - 2.0 * covariance)
        self.assertLess(difference['difference']['variance'],
                        self.model.variance(low) + self.model.variance(full))

    def test_shared_background_subtraction_keeps_half_full_closure(self):
        background = self.model.series('SB', [15., 90., 180.], [20., 100., 200.])
        low = STATS.series_subtract(self.al, background, 1. / 6.)
        high = STATS.series_subtract(self.ah, background, 1. / 6.)
        full = STATS.series_subtract(self.af, background, 1. / 3.)
        for index in range(3):
            closure = low[index] + high[index] - full[index]
            self.assertAlmostEqual(closure.value, 0.0)
            self.assertLess(self.model.variance(closure), 1e-26)
            self.assertAlmostEqual(self.model.covariance(low[index], high[index]),
                                   background[index].model.variance(background[index]) / 36.)

    def test_fixed_purity_correction_has_shared_background_covariance(self):
        low = STATS.normalized_cf(self.al, self.bl, [1, 2])['bins'][0]
        high = STATS.normalized_cf(self.ah, self.bh, [1, 2])['bins'][0]
        sideband = self.model.variable('CSB', 1.1, .03)
        corrected_low = STATS.purity_corrected_cf(low, sideband, .9)
        corrected_high = STATS.purity_corrected_cf(high, sideband, .8)
        expected_covariance = ((1.0 - .9) / .9) * ((1.0 - .8) / .8) * .03
        self.assertAlmostEqual(self.model.covariance(corrected_low, corrected_high), expected_covariance)
        self.assertAlmostEqual(STATS.purity_corrected_cf(low, sideband, 1.0).value, low.value)
        for purity in (0.0, -1.0, 1.1, float('nan')):
            with self.assertRaises(ValueError):
                STATS.purity_corrected_cf(low, sideband, purity)

    def test_zero_and_negative_bins_are_not_fabricated_positive(self):
        a = self.model.series('ZERO_A', [0., -2., 10.], [0., 3., 10.])
        b = self.model.series('ZERO_B', [5., 0., 10.], [5., 0., 10.])
        output = STATS.normalized_cf(a, b, [2])
        self.assertEqual(output['bins'][0].value, 0.0)
        self.assertIn(0, output['warnings'])
        self.assertIsNone(output['bins'][1])
        self.assertIn(1, output['invalid'])
        other_b = self.model.series('POS_B', [5., 5., 10.], [5., 5., 10.])
        output = STATS.normalized_cf(a, other_b, [2])
        self.assertLess(output['bins'][1].value, 0.0)
        self.assertIn(1, output['warnings'])
        invalid_norm = STATS.normalized_cf(a, other_b, [0, 1])
        self.assertTrue(all(q is None for q in invalid_norm['bins']))

    def test_integrated_cf_matches_ratio_of_integrals_and_retains_norm_covariance(self):
        quantity = STATS.integrated_cf(self.al, self.bl, [1, 2], [0, 1])
        expected = ((3000. / 300.) * (130. / 1100.))
        self.assertAlmostEqual(quantity.value, expected)
        unit = STATS.integrated_cf(self.al, self.bl, [1, 2], [1, 2])
        self.assertAlmostEqual(unit.value, 1.0)
        self.assertLess(self.model.variance(unit), 1e-28)

    def test_correlated_background_fit_block(self):
        low, high = self.model.correlated_block('fit', [10., 12.], [[4., 3.], [3., 9.]])
        self.assertEqual(self.model.variance(low + high), 19.)
        self.assertEqual(self.model.variance(low - high), 7.)
        self.assertEqual(self.model.covariance(low, high), 3.)
        self.assertEqual(self.model.covariance_matrix([low, high]), [[4., 3.], [3., 9.]])

    def test_sliding_window_contains_central_bin(self):
        central = self.model.variable('central', 20., 25.)
        neighbors = self.model.variable('neighbors', 80., 100.)
        window = central + neighbors
        background = self.model.variable('window_background', 10., 4.)
        corrected = central * (1.0 - background / window)
        self.assertEqual(self.model.covariance(central, window), 25.)
        self.assertAlmostEqual(corrected.value, 18.)
        # d/d central = 1-B/W+central*B/W^2 = .92, not .9.
        self.assertAlmostEqual(corrected.gradient['central'], .92)

    def test_bad_inputs_fail_loudly(self):
        with self.assertRaises(ValueError):
            self.model.series('bad', [1.], [-1.])
        with self.assertRaises(ValueError):
            self.model.series('bad', [1.], [1., 2.])
        with self.assertRaises(ValueError):
            STATS.integral(self.al, [])
        with self.assertRaises(ValueError):
            STATS.integral(self.al, [0, 0])
        with self.assertRaises(IndexError):
            STATS.integral(self.al, [-1])
        with self.assertRaises(ValueError):
            self.model.variable(('AL', 0), 1., 1.)
        with self.assertRaises(ValueError):
            self.al[0] + STATS.ErrorModel().constant(1.)
        with self.assertRaises(ValueError):
            self.model.correlated_block('bad', [1., 1.], [[1., 2.], [2., 1.]])
        # Pairwise-valid but globally indefinite covariance must also fail.
        with self.assertRaises(ValueError):
            self.model.correlated_block('bad', [1., 1., 1.],
                                        [[1., -.9, -.9], [-.9, 1., -.9], [-.9, -.9, 1.]])
        self.assertNotIn(('bad', 0), self.model.variances)

    def test_singular_correlated_block_is_allowed(self):
        first, second = self.model.correlated_block('singular', [1., 1.], [[4., 4.], [4., 4.]])
        self.assertEqual(self.model.variance(first - second), 0.)


if __name__ == '__main__':
    unittest.main()
