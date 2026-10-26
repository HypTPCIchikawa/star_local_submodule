import json,pathlib,sys,numpy as np
P=pathlib.Path(__file__).resolve().parents[2];sys.path.insert(0,str(P/'tools/nuclear_pid_smooth_20261021'));import study
cfg=json.loads((P/'tools/nuclear_pid_smooth_20261021/settings.yaml').read_text())
b=json.loads((study.OLD/'baseline_001/baseline.json').read_text())['current']['tables'];original={k:b[k] for k in ['mean','sigma','domain']};candidate=json.loads(json.dumps(original));evidence={}
for s in ['he3','he4']:
 i=study.SP.index(s);points=[]
 low=study.OLD/('all_'+s+'_gaussian_linear_nominal_free001')
 high=study.OLD/'highrigidity_tof_20261017_001'/('native_'+s+'_gaussian_linear_001')
 for folder in [low,high]:
  for r in json.loads((folder/'cutoff_fits.json').read_text())['slices']:
   c=r['components'][s];sel=c['selection']
   if cfg['calibrationRMin']<=r['r']<=cfg['calibrationRMax'] and r.get('targetFreeResponseSupported') and sel.get('purity') is not None and sel['purity']>cfg['calibrationMinimumPurity'] and sel['localResidualFrac']<=cfg['calibrationMaxLocalResidual']:
    points.append((r['r'],c['mean'],c['sigma'],str(folder)))
 evidence[s]={'points':points,'fits':{}}
 rr=np.array([x[0] for x in points]);xx=np.log10(rr)
 for kind,degree,limit in [('mean',3,cfg['meanMaximumShiftOverOriginalSigma']),('sigma',2,cfg['sigmaMaximumRelativeChange'])]:
  old=np.array(original[kind][i][:degree+1]);A=np.vander(xx,degree+1,increasing=True);target=np.array([p[1 if kind=='mean' else 2] for p in points]);scale=np.polynomial.polynomial.polyval(xx,original['sigma'][i]);weights=1/scale
  # Equal leverage in relative detector widths; no high-statistics point dominates.
  delta=np.linalg.lstsq(A*weights[:,None],(target-A@old)*weights,rcond=None)[0]
  grid=np.geomspace(cfg['auditRMin'],cfg['auditRMax'],cfg['auditGridCount']);gx=np.log10(grid);gs=np.polynomial.polynomial.polyval(gx,original['sigma'][i]);pred=np.polynomial.polynomial.polyval(gx,delta)
  factor=min(1.,limit/np.max(np.abs(pred)/gs));new=old+factor*delta
  candidate[kind][i]=list(new)+([0.] if degree==2 else [])
  evidence[s]['fits'][kind]={'unscaledDelta':delta.tolist(),'appliedScale':float(factor),'maximumChangeOverOriginalSigma':float(np.max(np.abs(factor*pred)/gs)),'normalizedRmsBefore':float(np.sqrt(np.mean(((target-A@old)/scale)**2))),'normalizedRmsAfter':float(np.sqrt(np.mean(((target-A@new)/scale)**2)))}
settings={'schemaVersion':1,'productionEnabled':False,'sourceGitRevision':cfg['referenceGit'],'original':original,'candidate':candidate,'adjustmentPolicy':{'d_t':'unchanged','heMeanMaximumShiftOriginalSigma':cfg['meanMaximumShiftOverOriginalSigma'],'heSigmaMaximumRelativeChange':cfg['sigmaMaximumRelativeChange'],'coefficientFitR':[cfg['calibrationRMin'],cfg['calibrationRMax']],'auditDomainR':[cfg['auditRMin'],cfg['auditRMax']],'weight':'equal residual in original detector sigma; Gaussian ALL and matching TOF free-target estimates, not a fixed response model'},'evidence':evidence}
study.OUT.mkdir(parents=True,exist_ok=True);study.save(study.OUT/'coefficients.json',settings);study.save(study.OUT/'coefficients.yaml',settings)
print(json.dumps({s:evidence[s]['fits'] for s in evidence},indent=2));print('candidate',candidate)
