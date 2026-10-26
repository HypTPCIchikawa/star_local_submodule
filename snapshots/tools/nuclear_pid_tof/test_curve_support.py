import copy
from pathlib import Path
import unittest

import yaml

from build_curves import collect_points


class SupportTests(unittest.TestCase):
    def documents(self):
        settings = yaml.safe_load((Path(__file__).parent/'tag_study_settings.yaml').read_text())['support']
        document = {'slices': []}
        for sample, names in (('TAG_D', ('d', 't', 'he3', 'he4')),
                              ('TAG_T', ('t',)), ('TAG_HE3', ('he3',))):
            components = {name: dict(mean=20+10*i, sigma=2, meanError=.1, sigmaError=.1,
                                    effectiveCount=1000, windowYield=1000,
                                    conditionalDropDeltaD=100, detectorAcceptance=1,
                                    centerOutsideDetector=False) for i, name in enumerate(names)}
            document['slices'].append(dict(sample=sample, lo=.8, hi=.81, r=.805,
                                          fitLo=5, fitHi=75, numericalConverged=True,
                                          atBoundary=False, conditionalNdf=50, chi2=50,
                                          components=components))
        return document, copy.deepcopy(document), settings

    def test_support_is_not_identity(self):
        a, b, s = self.documents()
        accepted, audit = collect_points(a, b, s)
        self.assertEqual({n: len(p) for n, p in accepted.items()}, dict(d=1, t=1, he3=1, he4=1))
        self.assertTrue(all(not row['identityValidated'] for row in audit))

    def test_background_shift_rejected(self):
        a, b, s = self.documents()
        b['slices'][0]['components']['d']['mean'] += 10
        accepted, audit = collect_points(a, b, s)
        self.assertEqual(len(accepted['d']), 0)
        self.assertIn('background_model_mean_instability', audit[0]['failures'])

    def test_he4_missing_nuisance_reference(self):
        a, b, s = self.documents()
        a['slices'][1]['numericalConverged'] = False
        accepted, audit = collect_points(a, b, s)
        self.assertEqual(len(accepted['he4']), 0)
        row = next(row for row in audit if row['species'] == 'he4')
        self.assertFalse(row['acceptedNumericalShape'])
        self.assertTrue(any(f.startswith('nuisance_TAG_shape_not_supported:') for f in row['failures']))

    def test_coarse_cannot_be_intrinsic_calibration(self):
        a, b, s = self.documents()
        for d in (a, b):
            for row in d['slices']:
                row['hi'] = .85
        with self.assertRaises(ValueError):
            collect_points(a, b, s)

    def test_bad_shape_not_hidden_by_high_yield(self):
        a, b, s = self.documents()
        a['slices'][0]['chi2'] = 1.e12
        accepted, audit = collect_points(a, b, s)
        self.assertEqual(len(accepted['d']), 0)
        self.assertIn('linear:poor_spectrum_shape', audit[0]['failures'])


if __name__ == '__main__':
    unittest.main()
