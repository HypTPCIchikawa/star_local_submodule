#!/usr/bin/env python3
"""Weighted polynomial least squares via reorthogonalized QR, no host BLAS dependency."""
import copy,math,pathlib,sys
HERE=pathlib.Path(__file__).resolve().parent;sys.path.insert(0,str(HERE));import run
R,C=run.R,run.C;original=run.read(R/'git_baseline.json')['tables'];candidate=copy.deepcopy(original);points=[]
def dot(a,b):return math.fsum(x*y for x,y in zip(a,b))
def qrsolve(columns,y):
 n=len(columns);q=[];tri=[[0.]*n for _ in range(n)]
 for j,col in enumerate(columns):
  v=list(col)
  for again in range(2):
   for k in range(j):
    a=dot(q[k],v);tri[k][j]+=a;v=[x-a*z for x,z in zip(v,q[k])]
  norm=math.sqrt(dot(v,v));assert norm>1e-12;tri[j][j]=norm;q.append([x/norm for x in v])
 rhs=[dot(z,y) for z in q];x=[0.]*n
 for j in reversed(range(n)):x[j]=(rhs[j]-math.fsum(tri[j][k]*x[k] for k in range(j+1,n)))/tri[j][j]
 residual=[y[i]-math.fsum(columns[j][i]*x[j] for j in range(n)) for i in range(len(y))]
 assert max(abs(dot(col,residual)) for col in columns)<1e-8*max(1.,math.sqrt(dot(y,y)))
 return x
for region in C['regions']:
 for row in run.read(R/('baseline_'+region['name'])/'cutoff_fits.json')['slices']:
  c=row['components']['he4'];s=c['selection'];reason=[]
  if not C['calibrationRMin']<=row['r']<=C['calibrationRMax']:reason.append('outside_calibration_R')
  if not row['numericalValid'] or not row['targetFreeResponseSupported']:reason.append('response_not_supported')
  if s['localResidualFrac'] is None or s['localResidualFrac']>C['calibrationMaxLocalResidual']:reason.append('local_residual')
  points.append(dict(region=region['name'],r=row['r'],mean=c['mean'],sigma=c['sigma'],used=not reason,reasons=reason))
good=[p for p in points if p['used']];assert len(good)>=10
rr=[p['r'] for p in good];xx=[math.log10(r) for r in rr];grid=[C['auditRMin']*(C['auditRMax']/C['auditRMin'])**(i/(C['auditGridCount']-1)) for i in range(C['auditGridCount'])];gs=[run.poly(original['sigma'][3],r) for r in grid];scale=[run.poly(original['sigma'][3],r) for r in rr];fits={}
for kind,degree,limit in [('mean',3,C['meanMaximumShiftOverOriginalSigma']),('sigma',2,C['sigmaMaximumRelativeChange'])]:
 old=original[kind][3][:degree+1];target=[p[kind] for p in good];y=[(t-run.poly(old,r))/s for t,r,s in zip(target,rr,scale)];A=[[x**j/s for x,s in zip(xx,scale)] for j in range(degree+1)];delta=qrsolve(A,y);change=[run.poly(delta,r) for r in grid];factor=min(1.,limit/max(abs(a)/s for a,s in zip(change,gs)));new=[a+factor*b for a,b in zip(old,delta)];candidate[kind][3]=new+([0.] if kind=='sigma' else [])
 def rms(cs):return math.sqrt(math.fsum(((t-run.poly(cs,r))/s)**2 for t,r,s in zip(target,rr,scale))/len(rr))
 assert rms(new)<=rms(old)+1e-12
 fits[kind]=dict(rawDelta=delta,scaleApplied=factor,maximumChangeOverOriginalSigma=max(abs(factor*a)/s for a,s in zip(change,gs)),normalizedRmsBefore=rms(old),normalizedRmsAfter=rms(new))
for s in range(3):assert all(candidate[k][s]==original[k][s] for k in ['mean','sigma','domain'])
for kind in ['mean','sigma']:assert min(run.poly(candidate[kind][3],r) for r in grid)>0
result=dict(sourceGitRevision=C['referenceGit'],original=original,candidate=candidate,onlyHe4Changed=True,configuration=C,evidence=points,coefficientFits=fits);run.save(R/'coefficients.json',result);run.save(R/'coefficients.yaml',result)
rows=['r\t'+'\t'.join(f'{mode}_{s}_{kind}' for mode in ['original','candidate'] for s in run.SP for kind in ['mean','sigma'])]
for r in grid:rows.append(format(r,'.17g')+'\t'+'\t'.join(format(run.poly(result[mode][kind][s],r),'.17g') for mode in ['original','candidate'] for s in range(4) for kind in ['mean','sigma']))
(R/'curves.tsv').write_text('\n'.join(rows)+'\n');print('CALIBRATION',len(good),'/',len(points),'points;',fits,'He4:',candidate['mean'][3],candidate['sigma'][3],flush=True)

# Export the complete original-function header with only the two He4 rows changed.
import re,difflib
source=(R/'git_NuclearIdDeDxVsMom.h').read_text();header=source
for kind,array in [('mean','kMean'),('sigma','kSigma')]:
 block=header.split(array+'[kNSpecies] = {',1)[1].split('};',1)[0];row=re.findall(r'\{([^{}]+)\}',block)[3];before='{'+row+'}';after='{'+', '.join(format(v,'.17g') for v in candidate[kind][3])+'}';assert header.count(before)==1;header=header.replace(before,after,1)
folder=R/'candidate_pid';folder.mkdir(exist_ok=True);(folder/'NuclearIdDeDxVsMom.h').write_text(header);(folder/'he4_coefficients.patch').write_text(''.join(difflib.unified_diff(source.splitlines(True),header.splitlines(True),fromfile='Git8f7d539/NuclearIdDeDxVsMom.h',tofile='candidate/NuclearIdDeDxVsMom.h')))
