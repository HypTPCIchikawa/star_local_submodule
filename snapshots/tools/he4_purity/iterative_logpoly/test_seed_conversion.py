import json
import math
from pathlib import Path
import unittest

import yaml

from run_iterative_fit import HERE, PROJECT, old_seeds, probability, validate


class SeedConversionTest(unittest.TestCase):
    def setUp(self):
        self.cfg=yaml.safe_load((HERE/'iterative_settings.yaml').read_text())
        self.audit=json.loads((PROJECT/self.cfg['sample']['oldStudy']/'purity.json').read_text())

    def test_settings(self):
        validate(self.cfg)

    def test_all_seeds_preserve_old_model(self):
        rows,audits=old_seeds(self.audit,self.cfg['engine'],list(range(210)))
        self.assertEqual(len(rows),210)
        for row,detail in zip(rows,audits):
            self.assertEqual(len(row),21)
            self.assertTrue(all(math.isfinite(v) for v in row if isinstance(v,(int,float))))
            self.assertGreater(row[14],0);self.assertGreater(row[17],0)
            self.assertGreaterEqual(row[15],0);self.assertGreaterEqual(row[18],0)
            old=self.audit['slices'][row[0]]
            if old['chosen']>=0:
                self.assertEqual(row[8],old['fits'][old['chosen']]['model'])
            self.assertTrue(0<=detail['fitWindow'][0]<detail['fitWindow'][1]<=120)

    def test_height_binwidth_yield(self):
        _,details=old_seeds(self.audit,self.cfg['engine'],list(range(210)))
        for item in details:
            source=self.audit['slices'][item['index']]
            old=source['fits'][item['seedFitIndex']]
            a=source['calibrationMean'];s=source['calibrationSigma'];p=old['parameters']
            component=item['components']['he4']
            # Independently recover the saved +/-2 legacy gate count.
            got=component['fullYield']*probability(a-2*s,a+2*s,component['mean'],component['sigma'])
            self.assertAlmostEqual(got/max(1.,old['gateSignal']),old['gateSignal']/max(1.,old['gateSignal']),places=9)
            self.assertAlmostEqual(component['mean'],a+s*p[1],places=12)
            self.assertAlmostEqual(component['sigma'],s*p[2],places=12)

    def test_small_tail_probability(self):
        self.assertGreater(probability(9,10,0,1),0)
        self.assertAlmostEqual(probability(-10,-9,0,1),probability(9,10,0,1),places=30)


if __name__=='__main__':
    unittest.main()
