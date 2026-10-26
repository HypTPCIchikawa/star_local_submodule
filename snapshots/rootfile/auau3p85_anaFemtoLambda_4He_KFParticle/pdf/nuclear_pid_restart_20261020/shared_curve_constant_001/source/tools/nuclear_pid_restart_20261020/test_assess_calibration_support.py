import copy
import unittest

import assess_calibration_support as study


class CalibrationSupportContracts(unittest.TestCase):
    def row(self, species='he4', identity='he4_b0051'):
        row = dict(id=identity, target=species, mode='TAG_D', rBin=51,
            fitCompleted=True, status=0, finite=True, pointAccepted=True,
            amplitudeConverged=True, yieldCovarianceValid=True, shapeCovarianceValid=True,
            atBoundary=False, foreignPeakInsideFit=False, curveAdoptionAllowed=False,
            purityAvailable=False, pearsonPerNdf=1., columnNames=[species, 'BG_left', 'BG_right'],
            columnWindowYields=[100., 10., 20.], effectiveCount=[90.],
            components=[dict(species=species, initialMean=84., initialSigma=5.,
                mean=84.5, sigma=5.1, meanError=.1, sigmaError=.2)])
        planned = dict(id=identity, target=species, mode='TAG_D', rBin=51,
                       r=.505, fitLo=72., fitHi=97., fitWindowAxisClipped=False)
        return row, planned

    def test_real_he4_zero_yield_pattern_cannot_identify_a_peak(self):
        row, planned = self.row()
        row['columnWindowYields'] = [0., 123.96212953468303, 2.03786996201282]
        row['effectiveCount'] = [0.]
        row['components'][0].update(meanError=814750.6366594653, sigmaError=2457223.7796063176)
        result = study.analyze_row(row, planned)
        self.assertTrue(result['numericalFlagsOK'])
        self.assertTrue(result['zeroTargetYield'])
        self.assertFalse(result['supportEligible'])

    def test_real_triton_zero_yield_pattern_cannot_identify_a_peak(self):
        row, planned = self.row('t', 't_b0061')
        row['columnWindowYields'] = [0., 0., 1166.9999944844699]
        row['effectiveCount'] = [0.]
        row['components'][0].update(meanError=132639.52992414622, sigmaError=271169.5436440366)
        result = study.analyze_row(row, planned)
        self.assertFalse(result['supportEligible'])
        self.assertTrue(result['zeroTargetEffectiveCount'])

    def test_valid_case_and_source_are_not_mutated(self):
        row, planned = self.row()
        before = copy.deepcopy(row)
        result = study.analyze_row(row, planned)
        self.assertTrue(result['supportEligible'])
        self.assertEqual(row, before)
        self.assertFalse(result['curveAdoptionAllowed'])

    def test_foreign_peak_is_not_an_independent_target_fit(self):
        row, planned = self.row()
        row['foreignPeakInsideFit'] = True
        self.assertFalse(study.analyze_row(row, planned)['supportEligible'])

    def test_failed_actual_migrad_cannot_be_erased_by_successful_hesse(self):
        row, planned = self.row()
        row.update(responseOptimizerDiagnosticsAvailable=True, responseMinimizeReturned=False,
            responseMinimizeConverged=False, responseHesseAttempted=True,
            responseHesseSucceeded=True, responseMigradStatus=4,
            responseHesseStatus=0, responseMigradEdm=.1)
        result = study.analyze_row(row, planned)
        self.assertFalse(result['phaseFlagsOK'])
        self.assertFalse(result['supportEligible'])

    def test_bad_errors_or_effective_counts_are_rejected_without_physics_threshold(self):
        for field, value in (('meanError', 0.), ('sigmaError', float('nan'))):
            row, planned = self.row()
            row['components'][0][field] = value
            self.assertFalse(study.analyze_row(row, planned)['supportEligible'])
        row, planned = self.row()
        row['effectiveCount'] = [0.]
        self.assertFalse(study.analyze_row(row, planned)['supportEligible'])

    def test_large_error_and_partial_roi_are_reported_without_an_implicit_cut(self):
        row, planned = self.row()
        row['components'][0].update(mean=96., meanError=20., sigmaError=10.)
        result = study.analyze_row(row, planned)
        self.assertTrue(result['supportEligible'])
        self.assertGreater(result['meanErrorOverSigma'], 1.)
        self.assertLess(result['rightCoverageSigma'], 1.)

    def test_split_widths_define_coverage_but_not_new_gaussian_sigma(self):
        row, planned = self.row()
        row.update(responseSigmaLeft=10., responseSigmaRight=2.)
        result = study.analyze_row(row, planned)
        self.assertAlmostEqual(result['leftCoverageSigma'], 1.25)
        self.assertAlmostEqual(result['rightCoverageSigma'], 6.25)
        self.assertEqual(result['sigma'], 5.1)

    def test_sampled_runs_break_on_an_unassessed_neighbor(self):
        rows = [dict(id=str(i), rBin=i, r=float(i)) for i in range(4)]
        runs = study.contiguous_runs(rows, {'0', '1', '3'})
        self.assertEqual([r['sampledPointCount'] for r in runs], [2, 1])
        self.assertTrue(all(not r['unsampledInteriorsAreValidated'] for r in runs))

    def test_target_mode_outside_or_on_fit_edge_cannot_calibrate_peak(self):
        for value in (71., 72., 97., 100.):
            row, planned = self.row()
            row['components'][0]['mean'] = value
            result = study.analyze_row(row, planned)
            self.assertFalse(result['targetModeInsideFit'])
            self.assertFalse(result['supportEligible'])

    def test_chi2_reporting_thresholds_require_explicit_positive_unique_values(self):
        for values in ([], [2., 2.], [0.], [float('nan')]):
            with self.assertRaises(ValueError):
                study.summarize([], values)


if __name__ == '__main__':
    unittest.main()
