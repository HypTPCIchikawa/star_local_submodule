import copy
import math
import unittest
from pathlib import Path
import yaml
import shared_curve_study as s
from shared_curve_analysis import good_runs
from polynomial_calibration import merge_recenter_audit
from test_polynomial_calibration import settings,saved
from polynomial_calibration import select_point

class Tests(unittest.TestCase):
    def test_no_extrapolation(self):
        c=dict(rMin=.5,rMax=5.,mean0=10.,mean1=-1.,mean2=.1,mean3=0.,sigma0=1.,sigma1=0.,sigma2=0.,logRatio0=0.,logRatio1=0.,logRatio2=0.)
        self.assertEqual(s.evaluate(c,1.),(10.,1.,1.))
        for r in (.499,5.001):
            with self.assertRaises(ValueError):s.evaluate(c,r)
    def test_invalid_width(self):
        c=dict(rMin=.5,rMax=5.,mean0=10.,mean1=0.,mean2=0.,mean3=0.,sigma0=-1.,sigma1=0.,sigma2=0.,logRatio0=0.,logRatio1=0.,logRatio2=0.)
        with self.assertRaises(ValueError):s.evaluate(c,1.)
    def test_he3_peak_and_left_shoulder_in_he4_roi(self):
        cfg=yaml.safe_load(Path(s.__file__).with_name('shared_curve_study.yaml').read_text())
        sh={name:[(1.,(m,1.,1.))] for name,m in zip(s.SPECIES,(10.,20.,35.,48.))}
        lo,hi,active,center,clip=s.select_roi(sh,'he4',cfg,cfg['variants']['linear'],dict(dedxRange=[0.,120.],dedxBinWidth=.1))
        self.assertIn('he3',active);self.assertLessEqual(lo,33.5);self.assertGreaterEqual(hi,51.);self.assertFalse(clip)
        self.assertNotIn('t',active)
    def test_runs_do_not_bridge_bad_fits(self):
        rows=[dict(r=.6+.05*i,assessed=i!=3,lowerModelPurity=.9) for i in range(7)]
        self.assertEqual([len(q) for q in good_runs(rows)],[3,3])
    def test_no_single_isolated_veto(self):
        rows=[dict(r=.6+.05*i,assessed=True,lowerModelPurity=.4 if i==2 else .9) for i in range(5)]
        self.assertEqual(good_runs(rows),[])
    def test_foreign_old_curve_is_review_not_rejection(self):
        p=select_point(saved(),settings()['quality'],25.,False)
        a=dict(searchBoundaryMode=False,searchAxisClipped=False,fitWindowAxisClipped=False,supportReviewRequired=True,nearestOriginalCurveToMode='he3')
        self.assertTrue(merge_recenter_audit(p,a,settings()['quality'])['selected'])
    def test_search_boundary_requires_review(self):
        p=select_point(saved(),settings()['quality'],25.,False)
        a=dict(searchBoundaryMode=True,searchAxisClipped=False,fitWindowAxisClipped=False,supportReviewRequired=True,nearestOriginalCurveToMode='he4')
        self.assertFalse(merge_recenter_audit(p,a,settings()['quality'])['selected'])
    def test_roi_clip_rejected(self):
        p=select_point(saved(),settings()['quality'],25.,False)
        a=dict(searchBoundaryMode=False,searchAxisClipped=False,fitWindowAxisClipped=True,supportReviewRequired=False,nearestOriginalCurveToMode='he4')
        self.assertFalse(merge_recenter_audit(p,a,settings()['quality'])['selected'])
if __name__=='__main__':unittest.main()
