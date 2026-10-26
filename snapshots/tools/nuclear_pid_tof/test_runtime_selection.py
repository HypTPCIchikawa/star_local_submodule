import importlib.util
import math
from pathlib import Path
import unittest

from runtime_selection import (SPECIES, current_active, current_parameters, gl3_nodes,
                               legacy_parameters, load_current_selection, parse_legacy_tables)

PROJECT = Path(__file__).resolve().parents[2]


class RuntimeSelectorTest(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.runtime = load_current_selection(PROJECT, PROJECT/'config/mainconf/main_auau3p85_anaFemtoLambda_4He_KFParticle.yaml')
        spec = importlib.util.spec_from_file_location('frozen_helium_export', PROJECT/'tools/he4_purity/export_runtime_calibration.py')
        cls.reference = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(cls.reference)

    def test_actual_header_tables(self):
        tables = self.runtime['tables']
        self.assertEqual(tables['mean'], [list(x) for x in self.reference.MEAN])
        self.assertEqual([r[:3] for r in tables['sigma']], [list(x) for x in self.reference.SIGMA])

    def test_legacy_guards_and_seed_extrapolation(self):
        t = self.runtime['tables']
        self.assertEqual(legacy_parameters(t, 2, .35), (0., 0.))
        self.assertGreater(legacy_parameters(t, 2, .35, True)[0], 0)
        self.assertEqual(legacy_parameters(t, 0, 0), (0., 0.))

    def test_matches_existing_optin_export_at_6000_points(self):
        for species in range(4):
            for i in range(1, 1501):
                r = .1+i*.002
                got = current_parameters(self.runtime, species, r)
                want = self.reference.calibrated_parameters(self.runtime['calibration'], species, r)
                for a, b in zip(got[:2], want[:2]):
                    self.assertAlmostEqual(a, b, delta=1.e-10*max(1., abs(b)))

    def test_all_knots_edges_and_gaps(self):
        for species in (2, 3):
            item = self.runtime['calibration']['species'][SPECIES[species]]
            values = item['domain']+[(a['rigidity']+b['rigidity'])/2. for a,b in zip(item['knots'],item['knots'][1:])]
            for knot in item['knots']:
                values.extend([knot['rigidity'],knot['rigidity']-1.e-13,knot['rigidity']+1.e-13])
            for r in values:
                got = current_parameters(self.runtime, species, r)
                want = self.reference.calibrated_parameters(self.runtime['calibration'], species, r)
                for a,b in zip(got[:2],want[:2]):
                    self.assertAlmostEqual(a,b,delta=1.e-10*max(1.,abs(b)))

    def test_strict_rigidity_upper_cut(self):
        self.assertTrue(current_active(self.runtime, 0, 2.49999))
        self.assertFalse(current_active(self.runtime, 0, 2.5))
        self.assertFalse(current_active(self.runtime, 3, .35))

    def test_gauss_legendre_average(self):
        nodes=gl3_nodes(.40,.41)
        self.assertAlmostEqual(sum(w for r,w in nodes),1.)
        self.assertAlmostEqual(sum(r*w for r,w in nodes),.405)
        self.assertAlmostEqual(sum((r-.405)**2*w for r,w in nodes),.01**2/12.)


if __name__ == '__main__':
    unittest.main()
