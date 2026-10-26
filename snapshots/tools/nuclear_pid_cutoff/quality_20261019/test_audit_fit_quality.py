import unittest
from audit_fit_quality import choose, percentile


class AuditQualityTests(unittest.TestCase):
    cfg=dict(minimumModels=2,minimumNominalModels=2,minimumResponseFamilies=2,
             minimumBackgroundFamilies=2,purityThreshold=.5)

    def model(self,response='gaussian',background='linear',purity=.8,chi=1.,supported=True):
        return dict(label=response+'_'+background+'_nominal',response=response,background=background,
                    variation='nominal',purity=purity,maxChi2Ndf=chi,supportedWithoutChi2=supported)

    def test_quality_before_purity(self):
        models=[self.model(),self.model('lognormal','exponential'),self.model(purity=.1,chi=99)]
        result=choose(models,3,self.cfg)
        self.assertEqual(result['purityLow'],.8)
        self.assertEqual(result['retainedModels'],2)

    def test_no_cherrypick_good_fit_low_purity(self):
        result=choose([self.model(purity=.4),self.model('lognormal','exponential')],3,self.cfg)
        self.assertEqual(result['status'],'below_purity_threshold')

    def test_missing_family_is_unassessed(self):
        result=choose([self.model(),self.model('gaussian','exponential')],3,self.cfg)
        self.assertEqual(result['status'],'unassessed')

    def test_all_failed_is_unassessed_not_zero(self):
        result=choose([self.model(chi=50)],3,self.cfg)
        self.assertEqual(result['status'],'unassessed')
        self.assertIsNone(result['purityLow'])

    def test_existing_bad_support_excluded(self):
        result=choose([self.model(supported=False),self.model('lognormal','exponential')],3,self.cfg)
        self.assertEqual(result['status'],'unassessed')

    def test_quantile(self):
        self.assertEqual(percentile([1,2,3,4],.5),2.5)
        self.assertIsNone(percentile([], .5))


if __name__=='__main__':
    unittest.main()
