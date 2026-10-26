"""Small independent checks for the shared-polynomial launcher and seed policy."""
import copy
import math
import unittest
from pathlib import Path
from unittest.mock import patch

import yaml
import run_global_fit as fit


class GlobalLauncherTests(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.cfg = yaml.safe_load((fit.HERE/'global_logpoly_settings.yaml').read_text())
        cls.legacy = fit.legacy_coefficients(fit.PROJECT/'StMaker/common/NuclearIdDeDxVsMom.h')

    def test_legacy_same_formula(self):
        self.assertEqual(len(self.legacy), 14)
        self.assertEqual(self.legacy[:4], [28.8747,-68.8637,104.455,-56.9342])
        self.assertEqual(self.legacy[7:11], [39.974,-103.365,139.182,-67.1446])
        self.assertTrue(fit.feasible(self.legacy,self.legacy,self.cfg['engine']))

    def test_extrema_include_interior(self):
        self.assertEqual(fit.extrema([2.,-2.,1.],0.,2.),(1.,2.))
        for coefficients in ([4.,1.,-3.,1.],[2.,-3.,2.],[1.,1.]):
            exact = fit.extrema(coefficients,-.4,.4)
            grid = [fit.polynomial_value(coefficients,-.4+i*.0001) for i in range(8001)]
            self.assertLessEqual(exact[0],min(grid)+1.e-12)
            self.assertGreaterEqual(exact[1],max(grid)-1.e-12)
            self.assertAlmostEqual(exact[0],min(grid),places=7)
            self.assertAlmostEqual(exact[1],max(grid),places=7)

    def test_positive_width_and_order(self):
        for indices in ((4,5,6),(11,12,13)):
            coefficients = self.legacy[:]
            for index in indices:
                coefficients[index]=0.
            self.assertFalse(fit.feasible(coefficients,self.legacy,self.cfg['engine']))
        reversed_species = self.legacy[7:]+self.legacy[:7]
        self.assertFalse(fit.feasible(reversed_species,self.legacy,self.cfg['engine']))
        broken = self.legacy[:];broken[0]=float('nan')
        self.assertFalse(fit.feasible(broken,self.legacy,self.cfg['engine']))

    def test_blend_is_feasible_and_bounded(self):
        candidate = self.legacy[:];candidate[7]-=100.
        blended, fraction = fit.blend_candidate(candidate,self.legacy,self.cfg['engine'],self.cfg['seedGeneration'])
        self.assertGreater(fraction,0.);self.assertLess(fraction,1.)
        self.assertTrue(fit.feasible(blended,self.legacy,self.cfg['engine']))
        for old,target,actual in zip(self.legacy,candidate,blended):
            self.assertAlmostEqual(actual,old+fraction*(target-old))

    def test_common_low_seed_not_per_slice(self):
        generation = copy.deepcopy(self.cfg['seedGeneration'])
        generation['modes']=['legacy','legacy_low_shift']
        seeds,audit = fit.make_seeds(self.legacy,self.cfg['engine'],generation,
                                    fit.PROJECT/self.cfg['sample']['legacyKnots'])
        self.assertEqual(len(seeds),2)
        self.assertEqual(seeds[0][1],self.legacy)
        self.assertTrue(fit.feasible(seeds[1][1],self.legacy,self.cfg['engine']))
        self.assertEqual(len(seeds[1][1]),14)
        self.assertEqual(seeds[1][1][:7],self.legacy[:7])
        self.assertIn('not truth',audit[1]['interpretation'])

    def test_native_alignment_and_dataset(self):
        fit.validate(self.cfg)
        changes = [('year',2019),('production','P19ik')]
        for key,value in changes:
            cfg=copy.deepcopy(self.cfg);cfg['sample'][key]=value
            with self.assertRaises(ValueError):fit.validate(cfg)
        for key,value in [('sliceIndices',[0,0]),('sliceIndices',[True]),
                          ('rigidityStep',.015),('rigidityMin',.415),
                          ('dedxRebinFactor',1),('rQuadratureOrder',2),('localTolerance',float('nan'))]:
            cfg=copy.deepcopy(self.cfg);cfg['engine'][key]=value
            with self.assertRaises(ValueError):fit.validate(cfg)

    def test_original_summary_and_consumed_snapshot_guard(self):
        source,original,snapshot,summary=map(Path,('data.root','engine.C','snapshot.C','summary.json'))
        hashes={original:'original_sha',snapshot:'snapshot_sha',summary:'summary_sha'}
        initial={'size':10,'mtime_ns':123,'path':'data.root'}
        with patch.object(fit,'fingerprint',return_value=initial),patch.object(fit,'digest',side_effect=lambda p:hashes[Path(p)]):
            fit.check_preserved(source,initial,{str(original):'original_sha'},
                                {str(snapshot):'snapshot_sha'},summary,'summary_sha')
            hashes[summary]='changed_summary'
            with self.assertRaisesRegex(ValueError,'summary changed'):
                fit.check_preserved(source,initial,{str(original):'original_sha'},
                                    {str(snapshot):'snapshot_sha'},summary,'summary_sha')
            hashes[summary]='summary_sha';hashes[snapshot]='changed_compiled_source'
            with self.assertRaisesRegex(ValueError,'snapshot changed'):
                fit.check_preserved(source,initial,{str(original):'original_sha'},
                                    {str(snapshot):'snapshot_sha'},summary,'summary_sha')

    def test_input_fingerprint_guard(self):
        with patch.object(fit,'fingerprint',return_value={'size':11}):
            with self.assertRaisesRegex(ValueError,'ROOT changed'):
                fit.check_preserved(Path('data.root'),{'size':10},{},{},Path('summary.json'),'sha')


if __name__=='__main__':
    unittest.main()
