import math
import unittest

import run_cutoff as run
from freeze_baseline import validate_reference_cuts


class CutoffRunnerTests(unittest.TestCase):
    def test_nominal_and_sensitivity_keep_width_positive(self):
        cfg = dict(shapeVariations=['nominal', 'mean_minus', 'mean_plus', 'width_minus', 'width_plus'],
                   meanShiftSigma=.1, widthRelativeShift=.1)
        self.assertEqual(run.varied(10., 2., 'nominal', cfg), (10., 2.))
        self.assertEqual(run.varied(10., 2., 'mean_minus', cfg), (9.8, 2.))
        self.assertEqual(run.varied(10., 2., 'mean_plus', cfg), (10.2, 2.))
        self.assertEqual(run.varied(10., 2., 'width_minus', cfg), (10., 1.8))
        self.assertEqual(run.varied(10., 2., 'width_plus', cfg), (10., 2.2))

    def test_failed_tag_center_is_not_purity_zero_or_valid_support(self):
        support = dict(minWindowYield=50., minEffective=25., minDropDeltaD=25.,
                       minDetectorAcceptance=.95, maxMeanErrSigma=.25, maxSigmaRelErr=.25)
        component = dict(mean=-120., sigma=45., windowYield=500., effectiveCount=100.,
                         conditionalDropDeltaD=100., detectorAcceptance=.01,
                         meanError=1., sigmaError=1., zeroYield=False, active=True)
        failures = run.shape_failures(dict(numericalConverged=True, atBoundary=True), component, support)
        self.assertIn('tag_center_nonphysical_or_failed', failures)
        self.assertIn('minDetectorAcceptance', failures)
        self.assertIn('tag_shape_boundary', failures)

    def test_poor_global_chi2_is_not_sole_support_veto(self):
        support = dict(minWindowYield=50., minEffective=25., minDropDeltaD=25.,
                       minDetectorAcceptance=.95, maxMeanErrSigma=.25, maxSigmaRelErr=.25)
        component = dict(mean=10., sigma=1., windowYield=500., effectiveCount=100.,
                         conditionalDropDeltaD=100., detectorAcceptance=1.,
                         meanError=.01, sigmaError=.01, zeroYield=False, active=True)
        self.assertEqual(run.shape_failures(dict(numericalConverged=True, atBoundary=False,
                                                chi2=1.e8, conditionalNdf=100), component, support), [])

    def test_unknown_variation_is_rejected(self):
        with self.assertRaises(ValueError):
            run.varied(1., 1., 'silent_fallback', dict(shapeVariations=['nominal']))

    def test_named_target_boundaries_ignore_other_species_boundary(self):
        row = dict(shapeParameterNames=['d_mean','d_logSigma','he4_mean','he4_logSigma'],
                   shapePoint=[10., math.log(1.), -120., math.log(45.)],
                   _sourceNumeric=dict(meanAxisPadding=120.,minimumSigma=.05,
                       maximumSigma=120.,boundaryTolerance=1.e-6))
        d = run.named_component_boundaries(row, 'd', 'gaussian', dict(dedxRange=[0.,120.]))
        h = run.named_component_boundaries(row, 'he4', 'gaussian', dict(dedxRange=[0.,120.]))
        self.assertFalse(any(x['atBoundary'] for x in d))
        self.assertTrue(any(x['atBoundary'] for x in h))

    def test_free_target_sensitivity_changes_nuisance_not_target(self):
        cfg = dict(allFit=dict(targetShapeMode='free_target_mean_sigma'),
                   shapeVariations=['nominal','width_plus'],widthRelativeShift=.1)
        c = dict(mean=10.,sigma=2.)
        self.assertEqual(run.varied_response(c,'width_plus',cfg,'he4','he4'), (10.,2.))
        self.assertEqual(run.varied_response(c,'width_plus',cfg,'he3','he4'), (10.,2.2))
        cfg['allFit']['targetShapeMode'] = 'fixed'
        self.assertEqual(run.varied_response(c,'width_plus',cfg,'he4','he4'), (10.,2.2))

    def test_reference_freeze_rejects_new_unmodeled_production_domain(self):
        validate_reference_cuts(dict(highRigidityTof=dict(enabled=False)))
        for cuts in (dict(highRigidityTof=dict(enabled=True)),
                     dict(speciesMinPOverQ=dict(he4=.65)),dict(speciesPOverQVeto=dict(d=[]))):
            with self.assertRaises(ValueError):
                validate_reference_cuts(cuts)


if __name__ == '__main__':
    unittest.main()
