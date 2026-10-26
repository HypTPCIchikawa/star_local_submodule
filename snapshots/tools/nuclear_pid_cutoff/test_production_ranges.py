import unittest
from prepare_production_ranges import merge, complement


class ProductionRangesTest(unittest.TestCase):
    def test_decimal_roundoff(self):
        self.assertEqual(merge([(2.0, 2.05), (2.0500000000000003, 2.1)]), [[2.0, 2.1]])

    def test_hole_not_bridged(self):
        self.assertEqual(merge([(.7, .85), (.95, 1.75), (1.8, 2), (2, 3.2)]), [[.7, .85], [.95, 1.75], [1.8, 3.2]])

    def test_veto(self):
        self.assertEqual(complement([[.75, 2.5]], 5), [dict(min=2.5, max=5)])

    def test_tail_and_holes(self):
        self.assertEqual(complement([[.7, .85], [.95, 1.75], [1.8, 3.2]], 5), [dict(min=.85, max=.95), dict(min=1.75, max=1.8), dict(min=3.2, max=5)])

    def test_overlap_fails(self):
        with self.assertRaises(ValueError):
            merge([(1, 2), (1.9, 3)])


if __name__ == '__main__':
    unittest.main()
