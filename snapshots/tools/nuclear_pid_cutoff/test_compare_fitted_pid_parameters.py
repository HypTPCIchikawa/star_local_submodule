import math
import unittest
import compare_fitted_pid_parameters as comparison

class Provider:
    @staticmethod
    def current_parameters(runtime,species,r):return 10.,2.,'fixture'

class FittedComparisonTest(unittest.TestCase):
    def row(self):
        return dict(label='gaussian_linear_nominal',species='d',r='2.005',rLo='2',rHi='2.01',mean='9',sigma='1',meanError='.2',sigmaError='.1',numericalValid='1',fitResponseQuality='1',targetFreeSupported='0')

    def test_fit_and_selection_are_separate(self):
        row=comparison.enrich(self.row(),'high_native_TOF',Provider,{k:{} for k in ('oldCanonical','oldLegacy','new')})
        self.assertEqual(row['centerPull'],1)
        self.assertEqual(row['widthRatio'],2)
        self.assertAlmostEqual(row['centerPullDiagonalError'],math.sqrt(.05))
        self.assertAlmostEqual(row['widthRatioError'],.2)
        self.assertEqual(row['plotted'],1)
        self.assertEqual(row['targetFreeSupported'],0)

    def test_failed_point_is_not_a_zero_fit(self):
        row=self.row();row['mean']='nan';row['numericalValid']='0'
        value=comparison.enrich(row,'low_ALL',Provider,{k:{} for k in ('oldCanonical','oldLegacy','new')})
        self.assertEqual(value['plotted'],0)
        self.assertIsNone(value['mean'])
        self.assertIsNone(value['centerPull'])

    def test_only_nominal_models_are_displayed(self):
        row=self.row();row['label']='gaussian_linear_width_plus'
        self.assertIsNone(comparison.enrich(row,'low_ALL',Provider,{}))

if __name__=='__main__':unittest.main()
