import copy
from pathlib import Path
import sys
import unittest

sys.path.insert(0, str(Path(__file__).resolve().parents[1]/'tools/kfparticle_summary'))
from compare_mass_halves import coarse_fit_inputs
from mass_half_statistics import ErrorModel


class CoarsePurityTests(unittest.TestCase):
    def fixture(self):
        rows = []
        for i,bg in enumerate(([10.,40.],[20.,20.])):
            rows.append(dict(background=bg,backgroundCovariance=[[25.,2.],[2.,36.]],
                valid=[True,True],regionReason=['',''],full={'valid':True},
                windowOriginalFirstBin=1+10*i,windowOriginalLastBin=10+10*i))
        return dict(method5Base={s:[dict(counts=[10.,20.],error2=[10.,20.]) for _ in range(20)]
                                for s in ('SE','ME')},
                    methods={'method5':{s:copy.deepcopy(rows) for s in ('SE','ME')}})

    def test_boundary_uses_both_purities_before_rebin(self):
        model = ErrorModel()
        data,problems = coarse_fit_inputs(model,self.fixture(),dict(rebin=4,purityRebin=10))
        self.assertEqual(problems,[])
        self.assertEqual(len(data['SE']['low']),5)
        self.assertAlmostEqual(data['SE']['low'][0].value,36.)
        self.assertAlmostEqual(data['SE']['low'][2].value,34.)
        self.assertAlmostEqual(data['SE']['low'][4].value,32.)
        self.assertAlmostEqual(data['SE']['high'][2].value,68.)
        self.assertAlmostEqual(data['SE']['full'][2].value,102.)
        self.assertNotEqual(model.covariance(data['SE']['low'][0],data['SE']['low'][1]),0.)

    def test_failed_coarse_region_masks_only_touching_fine_bins(self):
        group = self.fixture()
        group['methods']['method5']['SE'][0]['valid'][0] = False
        group['methods']['method5']['SE'][0]['regionReason'][0] = 'invalid background'
        data,problems = coarse_fit_inputs(ErrorModel(),group,dict(rebin=4,purityRebin=10))
        self.assertEqual([q is None for q in data['SE']['low']],[True,True,True,False,False])
        self.assertTrue(all(q is not None for q in data['SE']['high']))
        self.assertEqual(len(problems),1)


if __name__ == '__main__':
    unittest.main()
