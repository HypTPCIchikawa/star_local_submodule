import math
import unittest

from nuclear_curves import (evaluate_curve, local_smooth, pchip_slopes, pchip_value,
                            weighted_polynomial)


class CurveTests(unittest.TestCase):
    def grid(self, lo, hi, n):
        return [lo+(hi-lo)*i/(n-1) for i in range(n)]

    def close(self, left, right):
        for a, b in zip(left, right):
            self.assertAlmostEqual(a, b, places=10)

    def test_exact_weighted_cubic(self):
        x = self.grid(-.4, .4, 20)
        coeff = weighted_polynomial(x, [2+3*v-4*v*v+5*v**3 for v in x], [1.]*20, 3)
        self.close(coeff, [2, 3, -4, 5])

    def test_pchip_monotone_no_overshoot(self):
        x, y = [0., .1, .3, .8, 1.], [0., .01, .2, .21, 1.]
        slopes = pchip_slopes(x, y)
        values = [pchip_value(x, y, slopes, v) for v in self.grid(0, 1, 1001)]
        self.assertTrue(all(b >= a-1.e-14 for a, b in zip(values, values[1:])))
        for xx, yy in zip(x, y):
            self.assertAlmostEqual(pchip_value(x, y, slopes, xx), yy)

    def test_two_points_and_ends(self):
        self.assertEqual(pchip_slopes([0, 1], [2, 4]), [2, 2])
        self.assertAlmostEqual(pchip_value([0, 1], [2, 4], [2, 2], -.1), 1.8)

    def test_smoothing_exact_quadratic(self):
        x = self.grid(-.4, .4, 20)
        y = [2+v+v*v for v in x]
        self.close(local_smooth(x, y, [1.]*20, 9), y)

    def test_no_gap_or_domain_extrapolation(self):
        segment = dict(family='logpoly', domain=[.4, .5], coeffMean=[20], coeffSigma=[2],
                       uncertaintyKnots=[dict(r=.45, meanError=.1, sigmaError=.1)])
        c = dict(schemaVersion=1, productionEnabled=False, species={'he4': {'segments': [segment]}})
        self.assertTrue(evaluate_curve(c, 'he4', .45)['supported'])
        self.assertFalse(evaluate_curve(c, 'he4', .39)['supported'])
        self.assertFalse(evaluate_curve(c, 'he4', .51)['supported'])
        self.assertFalse(evaluate_curve(c, 'd', .45)['supported'])
        with self.assertRaises(ValueError):
            evaluate_curve(dict(c, productionEnabled=True), 'he4', .45)

    def test_positive_log_width(self):
        knots = [dict(r=.4, mean=100, logSigma=math.log(3)),
                 dict(r=.6, mean=50, logSigma=math.log(2))]
        xx = [math.log10(k['r']) for k in knots]
        segment = dict(family='pchip_logR', domain=[.4, .6], knots=knots,
                       meanSlopes=pchip_slopes(xx, [100, 50]),
                       logSigmaSlopes=pchip_slopes(xx, [math.log(3), math.log(2)]),
                       uncertaintyKnots=[dict(r=.5, meanError=1, sigmaError=.1)])
        c = dict(schemaVersion=1, productionEnabled=False, species={'he4': {'segments': [segment]}})
        for r in self.grid(.4, .6, 100):
            self.assertGreater(evaluate_curve(c, 'he4', float(r))['sigma'], 0.)


if __name__ == '__main__':
    unittest.main()
