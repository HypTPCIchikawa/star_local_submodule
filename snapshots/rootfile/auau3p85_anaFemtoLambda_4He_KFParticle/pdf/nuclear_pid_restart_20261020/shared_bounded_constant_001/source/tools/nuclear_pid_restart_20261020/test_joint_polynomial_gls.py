import copy
import unittest
import joint_polynomial_gls as g
import polynomial_calibration as b
from test_polynomial_calibration import saved,settings

def support():
    result=[]
    for i in range(10):
        q=b.select_point(saved(index=i,lo=.5+.3*i),settings()['quality'],25.,False)
        q['shapeNames']=['he4_mean','he4_logSigma','split_logRatio']
        q['shapeCovariance']=[[.0004,.00002,-.0001],[.00002,.0001,.00002],[-.0001,.00002,.0009]]
        result.append(q)
    return result

class Tests(unittest.TestCase):
    def test_cholesky(self):
        m=[[4.,1.,.5],[1.,2.,.25],[.5,.25,1.]];l=g.cholesky(m,1.e-10)
        for i in range(3):
            for j in range(3):self.assertAlmostEqual(sum(l[i][k]*l[j][k] for k in range(3)),m[i][j],places=12)
    def test_recovery(self):
        q=g.fit_joint_polynomials(support(),settings());targets=dict(mean=[20.,-5.,1.,-.1],sigma=[1.,-.1,.05],logRatio=[.1,.04,.02])
        for key,cs in targets.items():
            for a,z in zip(q['models'][key]['coefficients'],cs):self.assertAlmostEqual(a,z,places=8)
        self.assertLess(q['jointRegression']['chi2'],1.e-15)
    def test_cross_covariance(self):
        c=g.fit_joint_polynomials(support(),settings())['jointRegression']['conditionalCoefficientCovariance']
        self.assertTrue(any(abs(c[i][j])>1.e-12 for i in range(4) for j in range(7,10)))
    def test_birge(self):
        a=support()
        for q in a:q['chi2PerNdf']=1.
        c=copy.deepcopy(a)
        for q in c:q['chi2PerNdf']=4.
        x=g.fit_joint_polynomials(a,settings());y=g.fit_joint_polynomials(c,settings())
        for name in ('mean','sigma','logRatio'):
            for u,v in zip(x['models'][name]['coefficients'],y['models'][name]['coefficients']):self.assertAlmostEqual(u,v,places=9)
        for i in range(10):self.assertAlmostEqual(y['jointRegression']['conditionalCoefficientCovariance'][i][i],4*x['jointRegression']['conditionalCoefficientCovariance'][i][i],places=10)
    def test_invalid(self):
        a=support();a[0]['shapeNames'].reverse()
        with self.assertRaises(ValueError):g.fit_joint_polynomials(a,settings())
        a=support();a[0]['shapeCovariance'][0][1]=a[0]['shapeCovariance'][1][0]=100.
        with self.assertRaises(ValueError):g.fit_joint_polynomials(a,settings())
    def test_no_hidden_floor(self):
        c=settings();c['jointGls']['modelFloorsEnabled']=True
        with self.assertRaises(ValueError):g.fit_joint_polynomials(support(),c)
    def project(self,ratio=0.,sigma=1.,maximum=None):
        c=settings()
        if maximum is not None:c['jointGlsBounds']['maximumIterations']=maximum
        cs=[20.,0.,0.,0.,sigma,0.,0.,ratio,0.,0.]
        cov=[[float(i==j) for j in range(10)] for i in range(10)]
        return g.constrain_joint(cs,cov,support(),c)
    def test_feasible_projection_identity(self):
        q=self.project();self.assertTrue(q['diagnostics']['converged']);self.assertEqual(q['coefficients'],[20.,0.,0.,0.,1.,0.,0.,0.,0.,0.])
    def test_ratio_bound_projection(self):
        q=self.project(ratio=2.);self.assertTrue(q['diagnostics']['converged'])
        for i in range(101):
            x=__import__('math').log10(.5+(3.21-.5)*i/100)
            self.assertLessEqual(b.polynomial(q['coefficients'][7:],x),__import__('math').log(4.)+1.e-8)
        self.assertGreaterEqual(q['diagnostics']['mahalanobisPenalty'],0.)
    def test_sigma_bound_projection(self):
        q=self.project(sigma=.01);self.assertTrue(q['diagnostics']['converged'])
        for i in range(101):
            x=__import__('math').log10(.5+(3.21-.5)*i/100)
            self.assertGreaterEqual(b.polynomial(q['coefficients'][4:7],x),.05-1.e-8)
    def test_infeasible_projection_must_not_silently_pass(self):
        self.assertFalse(self.project(ratio=2.,maximum=0)['diagnostics']['converged'])
if __name__=='__main__':unittest.main()
