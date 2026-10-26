import unittest
from unittest.mock import patch
from copy import deepcopy
from pathlib import Path
from summarize_cutoff import aggregate_one, safe_intervals, models_supported, recommendation, numerical_replacement, fitted_response_failures


class CutoffSummaryTests(unittest.TestCase):
    def settings(self):
        return dict(allFit={}, support=dict(maxMeanErrSigma=.25,maxSigmaRelErr=.25,
            minWindowYield=50.,minEffective=25.),cutoff=dict(maxGateAbsoluteResidualFraction=.15,
            maxGateIntegralResidualFraction=.1, maxUnsupportedContaminationFraction=.01,
            minContiguousBins=3))

    def row(self, index, signal, total, supported=True):
        value = dict(signalSelected=signal, totalSelected=total, rawGateCount=total,
            gateAbsoluteResidualCount=0., conditionalError=.01, purity=signal/total,
            numericalValid=True, available=True, selectionCovered=True, axisCovered=True,
            gateActive=True, identitySupport=supported, boundaryTarget=False,
            selectedColumnCounts=[signal, 0., 0., 0., total-signal])
        return dict(index=index,columnWindowYields=[1000.,0.,0.,0.],
            components=dict(d=dict(selection=value,mean=10.,sigma=1.,meanError=.01,sigmaError=.01)))

    def test_sum_signal_over_sum_total_not_mean_purity(self):
        rows = [self.row(0, 900., 1000.), self.row(1, 10., 100.)]
        audit = {i:dict(failures={s:[] for s in ('d','t','he3','he4')}) for i in (0,1)}
        value = aggregate_one(rows, 'd', audit, self.settings())
        self.assertAlmostEqual(value['purity'], 910./1100.)
        self.assertNotAlmostEqual(value['purity'], .5)
        self.assertTrue(value['supported'])

    def test_unsupported_does_not_turn_into_zero_purity(self):
        row = self.row(0, 900., 1000., supported=False)
        audit = {0:dict(failures={s:[] for s in ('d','t','he3','he4')})}
        value = aggregate_one([row], 'd', audit, self.settings())
        self.assertFalse(value['supported'])
        self.assertEqual(value['purity'], .9)
        self.assertIn('identitySupport', value['reasons'])

    def test_holes_are_not_bridged(self):
        cfg = self.settings()
        rows = [dict(rLo=.4+i*.05, rHi=.45+i*.05, quality=i!=3) for i in range(7)]
        self.assertEqual(len(safe_intervals(rows,cfg)), 2)
        self.assertAlmostEqual(safe_intervals(rows,cfg)[0][1], .55)
        self.assertAlmostEqual(safe_intervals(rows,cfg)[1][0], .6)

    def test_each_target_support_uses_its_own_four_models(self):
        models = {str(i):dict(purity=.9,supported=True) for i in range(4)}
        self.assertTrue(models_supported(models))
        models['3']['supported'] = False
        self.assertFalse(models_supported(models))

    def test_four_targets_times_four_models_do_not_poison_other_targets(self):
        targets = {s:{str(i):dict(purity=.9,supported=True) for i in range(4)}
                   for s in ('d','t','he3','he4')}
        self.assertEqual(sum(len(v) for v in targets.values()), 16)
        self.assertTrue(all(models_supported(v) for v in targets.values()))
        targets['he4']['1']['supported'] = False
        self.assertFalse(models_supported(targets['he4']))
        self.assertTrue(all(models_supported(targets[s]) for s in ('d','t','he3')))

    def test_a_high_r_hole_prevents_single_lower_cut_certification(self):
        cfg = self.settings()
        cfg.update(rigidityMax=2., numeric=dict(axisTolerance=1.e-9))
        rows = [dict(rLo=1.5+i*.05,rHi=1.55+i*.05,quality=True) for i in range(10)]
        self.assertEqual(recommendation(rows,cfg)['lowerCutoff'], 1.5)
        rows[8]['quality'] = False
        value = recommendation(rows,cfg)
        self.assertIsNone(value['lowerCutoff'])
        self.assertFalse(value['singleLowerCutSufficient'])

    def test_free_target_requires_independent_covariance_and_identity_screen(self):
        cfg = self.settings()
        cfg['allFit']['targetShapeMode'] = 'free_target_mean_sigma'
        row = self.row(0,900.,1000.)
        audit = {0:dict(failures={s:[] for s in ('d','t','he3','he4')})}
        first = aggregate_one([row],'d',audit,cfg)
        self.assertFalse(first['supported'])
        self.assertIn('free_target_shape_or_identity_screen',first['reasons'])
        row['components']['d']['selection']['targetFreeResponseSupported'] = True
        self.assertTrue(aggregate_one([row],'d',audit,cfg)['supported'])

    def test_failed_fit_does_not_erase_independently_measured_track_count(self):
        row = self.row(0,900.,1000.)
        row['rBin'] = 41
        row['components']['d']['selection'].update(signalSelected=0.,totalSelected=0.,
            purity=None,numericalValid=False,available=False,rawGateCount=0.)
        audit = {0:dict(failures={s:[] for s in ('d','t','he3','he4')})}
        neutral = {(41,'d'):dict(rawGateCount=1000.)}
        value = aggregate_one([row],'d',audit,self.settings(),neutral)
        self.assertEqual(value['rawSelected'],1000.)
        self.assertIsNone(value['purity'])
        self.assertFalse(value['supported'])

    def test_full_run_numerical_replacement_cannot_change_physics(self):
        cfg = dict(numeric=dict(amplitudeMaxIterations=1000,
            amplitudeGradientTolerance=1.e-6,shapeMaxCalls=6000),
            cutoff=dict(purityThreshold=.5))
        modified = deepcopy(cfg)
        modified['numeric']['amplitudeGradientTolerance'] = 1.e-8
        context = dict(targetSpecies='he4',response='gaussian',background='exponential',variation='nominal')
        def reader(path):
            if path.name=='effective_settings.json':
                return deepcopy(cfg if path.parent.name=='old' else modified)
            return context.copy()
        with patch('summarize_cutoff.read',side_effect=reader), patch('summarize_cutoff.digest',return_value='sha'):
            result = numerical_replacement(Path('/tmp/old'),Path('/tmp/new'))
            self.assertIn('amplitudeGradientTolerance',result['numericDifferences'])
            modified['cutoff']['purityThreshold'] = .9
            with self.assertRaises(ValueError):
                numerical_replacement(Path('/tmp/old'),Path('/tmp/new'))

    def test_positive_hessian_but_weak_all_response_is_not_supported(self):
        cfg = self.settings()
        cfg['allFit']['targetShapeMode'] = 'free_target_mean_sigma'
        row = self.row(0,900.,1000.)
        row['components']['d'].update(meanError=3.,sigmaError=4.)
        failures = fitted_response_failures(row,'d',cfg)
        self.assertIn('fitted_meanError_over_sigma',failures)
        self.assertIn('fitted_sigmaError_over_sigma',failures)


if __name__ == '__main__':
    unittest.main()
