"""Same ten polynomial coefficients, using the full saved shape covariance."""
import math
import polynomial_calibration as base

def lower_solve(matrix,rhs):
    result=[0.]*len(rhs)
    for i in range(len(rhs)):result[i]=(rhs[i]-sum(matrix[i][j]*result[j] for j in range(i)))/matrix[i][i]
    return result

def cholesky(matrix,tolerance):
    n=len(matrix);base.require(n>0 and all(len(q)==n for q in matrix),'Covariance must be square')
    base.require(all(base.finite(v) for q in matrix for v in q),'Nonfinite covariance')
    for i in range(n):
        base.require(matrix[i][i]>0,'Nonpositive covariance diagonal')
        for j in range(i):base.require(abs(matrix[i][j]-matrix[j][i])<=tolerance*max(1.,abs(matrix[i][j]),abs(matrix[j][i])),'Asymmetric covariance')
    result=[[0.]*n for _ in range(n)]
    for i in range(n):
        for j in range(i+1):
            v=(matrix[i][j]+matrix[j][i])/2.-sum(result[i][k]*result[j][k] for k in range(j))
            if i==j:base.require(v>0,'Covariance not positive definite');result[i][j]=math.sqrt(v)
            else:result[i][j]=v/result[j][j]
    return result

def qr_solve(design,observations,tolerance):
    n,p=len(observations),len(design[0]);base.require(n>p and len(design)==n and all(len(q)==p for q in design),'Invalid joint design')
    columns=[[q[j] for q in design] for j in range(p)];largest=max(math.sqrt(base.dot(q,q)) for q in columns)
    upper,orthogonal=[[0.]*p for _ in range(p)],[]
    for j,column in enumerate(columns):
        v=list(column)
        for _ in range(2):
            for i in range(j):
                z=base.dot(orthogonal[i],v);upper[i][j]+=z;v=[a-z*b for a,b in zip(v,orthogonal[i])]
        norm=math.sqrt(base.dot(v,v));base.require(norm>tolerance*largest,'Rank-deficient joint GLS');upper[j][j]=norm;orthogonal.append([a/norm for a in v])
    coefficients=base.upper_solve(upper,[base.dot(q,observations) for q in orthogonal]);inverse=[[0.]*p for _ in range(p)]
    for j in range(p):
        column=base.upper_solve(upper,[float(i==j) for i in range(p)])
        for i in range(p):inverse[i][j]=column[i]
    return coefficients,[[base.dot(inverse[i],inverse[j]) for j in range(p)] for i in range(p)]

def constrain_joint(coefficients,covariance,points,configuration):
    controls=configuration['jointGlsBounds'];original=list(coefficients)
    result=dict(coefficients=original,diagnostics=dict(converged=False))
    try:
        base.require(len(original)==10,'Ten coefficients required');cholesky(covariance,controls['symmetryTolerance'])
        rmin=min(q['rLo'] for q in points);rmax=max(q['rHi'] for q in points);margin=controls['logRatioInteriorMargin']
        bounds=dict(mean_min=controls['meanMin'],sigma_min=controls['sigmaMin'],logRatio_min=math.log(configuration['quality']['ratioMin'])+margin,logRatio_max=math.log(configuration['quality']['ratioMax'])-margin)
        constraints,keys=[],set()
        def make(kind,r):
            x=math.log10(r);a=[0.]*10
            if kind=='mean_min':start,degree,sign=0,3,1.
            elif kind=='sigma_min':start,degree,sign=4,2,1.
            else:start,degree,sign=7,2,-1. if kind=='logRatio_max' else 1.
            for k in range(degree+1):a[start+k]=sign*x**k
            ca=[base.dot(q,a) for q in covariance];variance=base.dot(a,ca);base.require(base.finite(variance) and variance>0,'Invalid constraint metric')
            scale=math.sqrt(variance)
            return dict(kind=kind,r=r,a=[v/scale for v in a],ca=[v/scale for v in ca],b=sign*bounds[kind]/scale,physicalScale=scale)
        def add(kind,r):
            key=(kind,format(r,'.17g'))
            if key not in keys:constraints.append(make(kind,r));keys.add(key)
        count=controls['denseGridPoints'];base.require(count>=3,'Dense grid too short')
        for i in range(count):
            for kind in bounds:add(kind,rmin+(rmax-rmin)*i/(count-1))
        if rmin<=1.<=rmax:
            for kind in bounds:add(kind,1.)
        def dual_solve():
            passive,multipliers=[],{};steps=0
            while steps<controls['maximumIterations']:
                current=[original[j]+sum(constraints[i]['ca'][j]*v for i,v in multipliers.items()) for j in range(10)]
                inactive=[i for i in range(len(constraints)) if i not in multipliers]
                if not inactive:return current,multipliers,steps
                worst=max(inactive,key=lambda i:constraints[i]['b']-base.dot(constraints[i]['a'],current))
                violation=constraints[worst]['b']-base.dot(constraints[worst]['a'],current)
                if violation<=controls['normalizedFeasibilityTolerance']:return current,multipliers,steps
                passive.append(worst);multipliers[worst]=0.
                while steps<controls['maximumIterations']:
                    steps+=1;gram=[[base.dot(constraints[i]['a'],constraints[j]['ca']) for j in passive] for i in passive]
                    lower=cholesky(gram,controls['symmetryTolerance']);base.require(min(lower[i][i]**2 for i in range(len(passive)))>controls['gramRankTolerance'],'Active Gram numerically singular')
                    rhs=[constraints[i]['b']-base.dot(constraints[i]['a'],original) for i in passive]
                    upper=[[lower[j][i] if j>=i else 0. for j in range(len(passive))] for i in range(len(passive))]
                    trial=dict(zip(passive,base.upper_solve(upper,lower_solve(lower,rhs))))
                    negative=[i for i in passive if trial[i]<-controls['dualTolerance']]
                    if not negative:
                        multipliers={i:max(0.,trial[i]) for i in passive if trial[i]>controls['dualTolerance']};passive=[i for i in passive if i in multipliers];break
                    alpha=min(multipliers[i]/(multipliers[i]-trial[i]) for i in negative)
                    for i in passive:multipliers[i]+=alpha*(trial[i]-multipliers[i])
                    passive=[i for i in passive if multipliers[i]>controls['dualTolerance']];multipliers={i:multipliers[i] for i in passive}
                else:raise ValueError('Dual active-set iteration limit')
            raise ValueError('Dual active-set iteration limit')
        total=0
        for iteration in range(controls['maximumCutPasses']):
            projected,multipliers,steps=dual_solve();total+=steps;analytic=[]
            for start,degree,kinds in ((0,3,('mean_min',)),(4,2,('sigma_min',)),(7,2,('logRatio_min','logRatio_max'))):
                for r in [rmin,rmax]+base.turning_points(projected[start:start+degree+1],rmin,rmax):
                    for kind in kinds:
                        c=make(kind,r);v=c['b']-base.dot(c['a'],projected);analytic.append(dict(kind=kind,r=r,physicalViolation=v*c['physicalScale'],normalizedViolation=v))
            bad=[q for q in analytic if q['physicalViolation']>controls['physicalFeasibilityTolerance']]
            grid=max(q['b']-base.dot(q['a'],projected) for q in constraints)
            if not bad and grid<=controls['normalizedFeasibilityTolerance']:
                delta=[a-b for a,b in zip(projected,original)]
                return dict(coefficients=projected,diagnostics=dict(converged=True,method='dual_active_set_metric_projection_analytic_cuts',bounds=bounds,rMin=rmin,rMax=rmax,cutPasses=iteration+1,activeSetIterations=total,constraintCount=len(constraints),maximumNormalizedGridViolation=max(0.,grid),analyticExtrema=analytic,mahalanobisPenalty=sum(v*base.dot(constraints[i]['a'],delta) for i,v in multipliers.items()),activeConstraints=[dict(kind=constraints[i]['kind'],r=constraints[i]['r'],multiplier=v) for i,v in multipliers.items()],originalCoefficients=original,coefficientDelta=delta,covarianceMeaning='Unconstrained GLS covariance is projection metric, NOT constrained-fit uncertainty',productionEnabled=False))
            for q in bad:add(q['kind'],q['r'])
        raise ValueError('Analytic-extremum cutting-plane limit')
    except (ValueError,ZeroDivisionError,OverflowError) as error:
        result['diagnostics'].update(error=str(error),productionEnabled=False);return result

def fit_joint_polynomials(points,configuration):
    controls=configuration['jointGls'];base.require(controls['birgeScale'] is True and controls['modelFloorsEnabled'] is False,'Explicit Birge/no-floor scenario required')
    degrees=configuration['polynomialDegrees'];base.require(degrees==dict(mean=3,sigma=2,logRatio=2),'Unchanged function degrees required')
    base.require(len(points)>=configuration['minimumSupportPoints'],'Insufficient joint support')
    names=('mean','sigma','logRatio');offsets=dict(mean=0,sigma=4,logRatio=7);p=10
    xs=[math.log10((q['rLo']+q['rHi'])/2.) for q in points];center,scale=(min(xs)+max(xs))/2.,(max(xs)-min(xs))/2.;base.require(scale>0,'Rank-deficient R support')
    rows,observations,support=[],[],[]
    for q in points:
        base.require(q['shapeNames']==[q['species']+'_mean',q['species']+'_logSigma','split_logRatio'],'Unexpected covariance ordering')
        cov=q['shapeCovariance'];base.require(len(cov)==3 and all(len(a)==3 for a in cov),'Missing shape covariance')
        jac=[1.,q['sigma'],1.];inflation=max(1.,q['chi2PerNdf'])
        transformed=[[inflation*jac[i]*cov[i][j]*jac[j] for j in range(3)] for i in range(3)];lower=cholesky(transformed,controls['symmetryTolerance'])
        basis=[sum(w*((math.log10(r)-center)/scale)**k for r,w in base.gl3_nodes(q['rLo'],q['rHi'])) for k in range(4)];block=[[0.]*p for _ in range(3)]
        for i,name in enumerate(names):
            for k in range(degrees[name]+1):block[i][offsets[name]+k]=basis[k]
        columns=[lower_solve(lower,[block[i][j] for i in range(3)]) for j in range(p)]
        rows.extend([[columns[j][i] for j in range(p)] for i in range(3)]);observed=[q[name] for name in names];observations.extend(lower_solve(lower,observed));support.append(dict(point=q,lower=lower,covariance=transformed,observed=observed))
    scaled,cov_scaled=qr_solve(rows,observations,configuration['rankTolerance']);transform=[[0.]*p for _ in range(p)]
    for name in names:
        start,degree=offsets[name],degrees[name]
        for i in range(degree+1):
            for j in range(i,degree+1):transform[start+i][start+j]=math.comb(j,i)*(-center)**(j-i)/scale**j
    coefficients=[base.dot(q,scaled) for q in transform]
    covariance=[[sum(transform[i][a]*cov_scaled[a][b]*transform[j][b] for a in range(p) for b in range(p)) for j in range(p)] for i in range(p)]
    bound_audit=None
    if configuration.get('boundedJointGls'):
        constrained=constrain_joint(coefficients,covariance,points,configuration);bound_audit=constrained['diagnostics']
        base.require(bound_audit['converged'],'Bounded GLS numerical failure: '+bound_audit.get('error','unknown'))
        coefficients=constrained['coefficients']
    models={}
    for name in names:
        start,degree=offsets[name],degrees[name]
        models[name]=dict(degree=degree,coefficients=coefficients[start:start+degree+1],regularizedConditionalCovariance=[q[start:start+degree+1] for q in covariance[start:start+degree+1]],rank=degree+1,scaledLogCenter=center,scaledLogScale=scale,residuals=[],weightMeaning='Marginal display of joint shape-covariance GLS, Birge inflation, no floors',chi2Meaning='Marginal residual sum is NOT joint objective')
    residuals=[]
    for saved in support:
        q=saved['point'];pred=[sum(w*base.polynomial(models[name]['coefficients'],math.log10(r)) for r,w in base.gl3_nodes(q['rLo'],q['rHi'])) for name in names];delta=[a-b for a,b in zip(saved['observed'],pred)];white=lower_solve(saved['lower'],delta)
        residuals.append(dict(id=q['id'],observed=saved['observed'],predicted=pred,residual=delta,whitenedResidual=white,jointChi2Contribution=base.dot(white,white),transformedInflatedCovariance=saved['covariance']))
        for i,name in enumerate(names):
            error=math.sqrt(saved['covariance'][i][i]);models[name]['residuals'].append(dict(id=q['id'],observed=saved['observed'][i],predicted=pred[i],residual=delta[i],regularizedError=error,pull=delta[i]/error))
    for name in names:
        m=models[name];m['regressionChi2']=sum(q['pull']**2 for q in m['residuals']);m['ndf']=len(points)-degrees[name]-1;m['regressionChi2PerNdf']=m['regressionChi2']/m['ndf']
    chi2=sum(q['jointChi2Contribution'] for q in residuals);ndf=3*len(points)-p
    return dict(models=models,jointRegression=dict(method='joint_shape_covariance_GLS',covarianceDeltaTransform='diag(1,sigma,1)',originalCovarianceOrder=['mean','logSigma','logRatio'],transformedCovarianceOrder=list(names),birgeScale=True,modelFloorsEnabled=False,chi2=chi2,ndf=ndf,chi2PerNdf=chi2/ndf,rank=p,conditionalCoefficientCovariance=covariance,boundedProjection=bound_audit,residuals=residuals,qualityMeaning='Gaussian approximation of saved local-shape likelihood; NOT forward Poisson validation',productionEnabled=False))
