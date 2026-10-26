"""Analytic tests for the study-only common-function initialization."""
import math
import random
import unittest

from smoothing import positive_quadratic, quadratic_minimum, to_original, value, weighted_qr


class SmoothingTest(unittest.TestCase):
    def test_exact_weighted_cubic(self):
        x = [-1.+2*i/40 for i in range(41)]
        expected = [3., -7., 1.2, -.03]
        self.assertLess(max(abs(a-b) for a,b in zip(expected, weighted_qr(
            x, [value(expected,t) for t in x], [1.+i*i for i in range(41)], 3))), 1.e-10)

    def test_original_basis(self):
        c, center, scale = [3., -7., 1.2, -.03], -.005, .398
        converted = to_original(c, center, scale)
        for i in range(100):
            t = -.5+i/100
            self.assertAlmostEqual(value(converted,t), value(c,(t-center)/scale), places=11)

    def test_positive_not_bernstein_family(self):
        # Positive polynomial with a NEGATIVE central Bernstein coefficient.
        # On normalized [-1,1] this is t^2/4 + .01.
        x = [-1.+i/10 for i in range(21)]
        c, audit = positive_quadratic(x,[.25*t*t+.01 for t in x],[1.]*21,.005)
        self.assertFalse(audit['constrained'])
        self.assertAlmostEqual(c[2],.25,places=12)
        self.assertAlmostEqual(c[0]-.25,-.24,places=12)

    def test_negative_interior(self):
        x = [-1.+i/10 for i in range(21)]
        c, audit = positive_quadratic(x,[t*t-.2 for t in x],[1.]*21,.05)
        self.assertTrue(audit['constrained'])
        self.assertGreaterEqual(quadratic_minimum(c)[1],.05-1.e-9)

    def test_offcenter_vertex(self):
        x = [-1.+i/15 for i in range(31)]
        c, audit = positive_quadratic(x,[(t-.37)**2-.2 for t in x],[1.]*31,.03)
        self.assertGreaterEqual(quadratic_minimum(c)[1],.03-1.e-9)
        self.assertTrue(audit['constrained'])

    def test_roundtrip_random(self):
        rng = random.Random(4030)
        for _ in range(20):
            c = [rng.uniform(-10,10) for _ in range(4)]
            converted = to_original(c,.004,.39)
            for t in [-.39,0.,.39]:
                self.assertAlmostEqual(value(c,(t-.004)/.39),value(converted,t),places=10)

    def test_bad_weights(self):
        with self.assertRaises(ValueError):
            weighted_qr([0,1,2],[0,1,2],[1,-1,1],1)
        with self.assertRaises(ValueError):
            weighted_qr([1,1,1],[0,1,2],[1,1,1],1)


if __name__ == '__main__':
    unittest.main()
