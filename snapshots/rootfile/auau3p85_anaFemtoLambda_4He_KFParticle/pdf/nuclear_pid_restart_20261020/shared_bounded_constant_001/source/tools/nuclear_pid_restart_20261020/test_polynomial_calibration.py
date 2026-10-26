import json
import math
from pathlib import Path
import tempfile
import unittest
import yaml
import polynomial_calibration as c

def settings():
    return yaml.safe_load(Path(__file__).with_name('polynomial_calibration.yaml').read_text())

def saved(species='he4',index=0,lo=.8):
    ns=c.gl3_nodes(lo,lo+.01)
    mean=sum(w*c.polynomial([20.,-5.,1.,-.1],math.log10(r)) for r,w in ns)
    sigma=sum(w*c.polynomial([1.,-.1,.05],math.log10(r)) for r,w in ns)
    ratio=math.exp(sum(w*c.polynomial([.1,.04,.02],math.log10(r)) for r,w in ns))
    ix=c.SPECIES.index(species)
    y,e,a=[0.]*6,[0.]*4,[False]*6
    y[ix],e[ix],a[ix]=1000.,500.,True
    row=dict(id=species+'_'+str(index),target=species,mode=c.EXPECTED_MODES[species],rLo=lo,rHi=lo+.01,fitLo=mean-3*sigma,fitHi=mean+3*sigma,status=0,responseMigradStatus=0,responseHesseStatus=0,covarianceStatus=3,atBoundary=False,foreignPeakInsideFit=False,columnWindowYields=y,effectiveCount=e,yieldCovarianceColumnActive=a,components=[dict(species=species,mean=mean,sigma=sigma,meanError=.02,sigmaError=.01)],responseWidthRatio=ratio,responseWidthRatioError=.03,pearsonPerNdf=2.)
    for key in ('fitCompleted','finite','pointAccepted','amplitudeConverged','yieldCovarianceValid','shapeCovarianceValid','responseOptimizerDiagnosticsAvailable','responseMinimizeReturned','responseMinimizeConverged','responseHesseAttempted','responseHesseSucceeded'):
        row[key]=True
    return row

class Tests(unittest.TestCase):
    def select(self,row,limit=25.,foreign=False):
        return c.select_point(row,settings()['quality'],limit,foreign)
    def test_gl3(self):
        with self.assertRaises(ValueError): c.gl3_nodes(0.,1.)
        self.assertAlmostEqual(sum(w for r,w in c.gl3_nodes(1.,2.)),1.)
    def test_degrees(self):
        self.assertEqual(settings()['polynomialDegrees'],dict(mean=3,sigma=2,logRatio=2))
        self.assertFalse(settings()['productionEnabled'])
    def test_good(self): self.assertTrue(self.select(saved())['selected'])
    def test_he4_tag(self):
        q=saved();q['mode']='TAG_HE4'
        with self.assertRaises(ValueError): self.select(q)
    def test_foreign(self):
        q=saved();q['foreignPeakInsideFit']=True
        self.assertTrue(self.select(q)['selected'])
        self.assertFalse(self.select(q,foreign=True)['selected'])
    def test_bad_not_purity(self):
        q=saved();q['responseHesseSucceeded']=False
        self.assertFalse(self.select(q)['selected'])
        self.assertNotIn('purity',self.select(q))
    def test_zero(self):
        q=saved();q['columnWindowYields'][3]=0.
        self.assertFalse(self.select(q)['selected'])
    def test_ratio_boundary(self):
        q=saved();q['atBoundary']=True;q['responseWidthRatio']=4.
        self.assertIn('response_ratio_boundary',self.select(q)['rejectionReasons'])
    def test_bad_error(self):
        q=saved();q['components'][0]['meanError']=-1
        self.assertFalse(self.select(q)['selected'])
    def test_chi_sensitivity(self):
        q=saved();q['pearsonPerNdf']=12.
        self.assertTrue(self.select(q)['selected'])
        self.assertFalse(self.select(q,10.)['selected'])
    def test_recovery(self):
        points=[self.select(saved(index=i,lo=.5+.3*i)) for i in range(10)]
        q=c.calibrate_species(points,settings())
        self.assertTrue(q['curveCandidateValid'])
        for actual,expected in zip(q['models']['mean']['coefficients'],[20.,-5.,1.,-.1]): self.assertAlmostEqual(actual,expected,places=8)
        for actual,expected in zip(q['models']['sigma']['coefficients'],[1.,-.1,.05]): self.assertAlmostEqual(actual,expected,places=8)
        self.assertFalse(q['extrapolationAllowed'])
        self.assertTrue(q['supportGapWarnings'])
    def test_rank(self):
        with self.assertRaises(ValueError): c.calibrate_species([self.select(saved(index=i)) for i in range(8)],settings())
    def test_turning(self): self.assertEqual(c.turning_points([1.,-2.,1.],1.,100.),[10.])
    def test_cli_preserves(self):
        with tempfile.TemporaryDirectory() as tmp:
            base=Path(tmp);source=base/'in';source.mkdir()
            rows=[saved(s,i,.5+.3*i) for s in c.SPECIES for i in range(8)]
            f=source/'fits.json'
            f.write_text(json.dumps(dict(status='completed',model='split_normal',rawCountsModified=False,rawReadback=True,predictionReadback=True,totalSlices=len(rows),results=rows)))
            before=f.read_bytes();out=base/'out'
            result=c.run(source,Path(__file__).with_name('polynomial_calibration.yaml'),out)
            self.assertEqual(before,f.read_bytes())
            self.assertFalse(result['productionEnabled'])
            self.assertEqual(len((out/'baseline/coefficients.tsv').read_text().splitlines()),5)
            with self.assertRaises(FileExistsError): c.run(source,Path(__file__).with_name('polynomial_calibration.yaml'),out)

if __name__=='__main__': unittest.main()
