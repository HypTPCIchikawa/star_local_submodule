#!/usr/bin/env python3
"""Candidate-only original-degree PID polynomial calibration; no production writes."""
import argparse
import csv
import hashlib
import json
import math
import copy
from pathlib import Path
import shutil
from collections import Counter
import yaml

SPECIES = ('d', 't', 'he3', 'he4')
EXPECTED_MODES = dict(d='TAG_D', t='TAG_T', he3='TAG_HE3', he4='TAG_D')
COEFFICIENT_HEADER = ['species', 'rMin', 'rMax'] + ['mean%d' % i for i in range(4)] + ['sigma%d' % i for i in range(3)] + ['logRatio%d' % i for i in range(3)]

def finite(v):
    return isinstance(v, (int, float)) and not isinstance(v, bool) and math.isfinite(v)

def require(ok, message):
    if not ok:
        raise ValueError(message)

def polynomial(coefficients, x):
    answer = 0.
    for coefficient in reversed(coefficients):
        answer = answer*x+coefficient
    return answer

def gl3_nodes(lo, hi):
    require(finite(lo) and finite(hi) and 0 < lo < hi, 'invalid positive R interval')
    mid, half = (lo+hi)/2., (hi-lo)/2.
    offset = math.sqrt(3./5.)*half
    return [(mid-offset, 5./18.), (mid, 4./9.), (mid+offset, 5./18.)]

def dot(a, b):
    return sum(x*y for x,y in zip(a,b))

def upper_solve(matrix, rhs):
    n = len(rhs)
    answer = [0.]*n
    for i in reversed(range(n)):
        answer[i] = (rhs[i]-sum(matrix[i][j]*answer[j] for j in range(i+1,n)))/matrix[i][i]
    return answer

def regress(points, variable, degree, weights, rank_tolerance):
    """Twice-reorthogonalized QR of the GL3-moment log10(R) design."""
    n, p = len(points), degree+1
    require(n > p, 'insufficient support points for regression and diagnostics')
    xs = [math.log10((q['rLo']+q['rHi'])/2.) for q in points]
    center, scale = (min(xs)+max(xs))/2., (max(xs)-min(xs))/2.
    require(scale > 0, 'rank-deficient R support')
    errors = []
    for point in points:
        factor = math.sqrt(max(1.,point['chi2PerNdf'])) if weights['birgeScale'] else 1.
        floor = weights['floors'][variable]
        floor_abs = floor if variable == 'logRatio' else floor*point['sigma']
        errors.append(math.hypot(factor*point[variable+'Error'], floor_abs))
    design = [[sum(w*((math.log10(r)-center)/scale)**k for r,w in gl3_nodes(q['rLo'],q['rHi'])) for k in range(p)] for q in points]
    columns = [[design[i][j]/errors[i] for i in range(n)] for j in range(p)]
    largest = max(math.sqrt(dot(c,c)) for c in columns)
    orthogonal, upper = [], [[0.]*p for _ in range(p)]
    for j,column in enumerate(columns):
        vector = list(column)
        for _ in range(2):
            for i in range(j):
                projection = dot(orthogonal[i],vector)
                upper[i][j] += projection
                vector = [a-projection*b for a,b in zip(vector,orthogonal[i])]
        norm = math.sqrt(dot(vector,vector))
        require(norm > rank_tolerance*largest, 'rank-deficient polynomial design')
        upper[j][j] = norm
        orthogonal.append([v/norm for v in vector])
    values = [q[variable]/e for q,e in zip(points,errors)]
    scaled = upper_solve(upper,[dot(q,values) for q in orthogonal])
    def choose(n,k):
        return math.factorial(n)//(math.factorial(k)*math.factorial(n-k))
    transform = [[choose(j,i)*(-center)**(j-i)/scale**j if j>=i else 0. for j in range(p)] for i in range(p)]
    coefficients = [dot(row,scaled) for row in transform]
    inverse = [[0.]*p for _ in range(p)]
    for j in range(p):
        column = upper_solve(upper,[float(i==j) for i in range(p)])
        for i in range(p):
            inverse[i][j] = column[i]
    cov_scaled = [[dot(inverse[i],inverse[j]) for j in range(p)] for i in range(p)]
    covariance = [[sum(transform[i][a]*cov_scaled[a][b]*transform[j][b] for a in range(p) for b in range(p)) for j in range(p)] for i in range(p)]
    residuals = []
    for point,error in zip(points,errors):
        predicted = sum(w*polynomial(coefficients,math.log10(r)) for r,w in gl3_nodes(point['rLo'],point['rHi']))
        residuals.append(dict(id=point['id'],observed=point[variable],predicted=predicted,residual=point[variable]-predicted,regularizedError=error,pull=(point[variable]-predicted)/error))
    chi2 = sum(q['pull']**2 for q in residuals)
    return dict(degree=degree,coefficients=coefficients,regularizedConditionalCovariance=covariance,regressionChi2=chi2,ndf=n-p,regressionChi2PerNdf=chi2/(n-p),rank=p,scaledLogCenter=center,scaledLogScale=scale,residuals=residuals,weightMeaning='Birge-scaled error plus model floor; not a statistical confidence interval')

def select_point(row, configuration, chi2_limit, reject_foreign):
    species = row.get('target')
    require(species in SPECIES, 'unknown saved target')
    require(row.get('mode') == EXPECTED_MODES[species], 'wrong m2 tag')
    reasons, flags = [], []
    checks = ('fitCompleted','finite','pointAccepted','amplitudeConverged','yieldCovarianceValid','shapeCovarianceValid','responseOptimizerDiagnosticsAvailable','responseMinimizeReturned','responseMinimizeConverged','responseHesseAttempted','responseHesseSucceeded')
    reasons.extend(k+'_false_or_missing' for k in checks if row.get(k) is not True)
    reasons.extend(k+'_nonzero_or_missing' for k in ('status','responseMigradStatus','responseHesseStatus') if row.get(k) != 0)
    if row.get('covarianceStatus',-1) < configuration['minCovarianceStatus']:
        reasons.append('covariance_status')
    if row.get('atBoundary') is not False:
        reasons.append('optimizer_boundary')
    if row.get('foreignPeakInsideFit'):
        flags.append('foreign_original_Git_curve_inside_ROI')
        if reject_foreign:
            reasons.append('foreign_original_Git_curve_inside_ROI')
    component = next((q for q in row.get('components',[]) if q.get('species')==species),{})
    index = SPECIES.index(species)
    columns, neffs, active = (row.get(k,[]) for k in ('columnWindowYields','effectiveCount','yieldCovarianceColumnActive'))
    target_yield = columns[index] if len(columns)>index else None
    neff = neffs[index] if len(neffs)>index else None
    if not (len(active)>index and active[index] is True and finite(target_yield) and target_yield>configuration['minTargetYield']):
        reasons.append('unidentified_or_zero_target_yield')
    if not (finite(neff) and neff>configuration['minEffectiveCount']):
        reasons.append('insufficient_effective_count')
    mean,sigma,mean_error,sigma_error = (component.get(k) for k in ('mean','sigma','meanError','sigmaError'))
    ratio,ratio_error = row.get('responseWidthRatio'),row.get('responseWidthRatioError')
    for key,value in [('mean',mean),('sigma',sigma),('meanError',mean_error),('sigmaError',sigma_error),('ratio',ratio),('ratioError',ratio_error)]:
        if not (finite(value) and value>0):
            reasons.append('invalid_'+key)
    lo,hi = row.get('fitLo'),row.get('fitHi')
    if not (finite(lo) and finite(hi) and finite(mean) and lo<mean<hi):
        reasons.append('mode_outside_ROI')
    if not (finite(lo) and finite(hi) and configuration['nativeDedxMin']<lo<hi<configuration['nativeDedxMax']):
        reasons.append('native_axis_clipping')
    if finite(ratio) and not configuration['ratioMin']<ratio<configuration['ratioMax']:
        reasons.append('response_ratio_boundary')
    chi2 = row.get('pearsonPerNdf')
    if not (finite(chi2) and 0<=chi2<=chi2_limit):
        reasons.append('local_chi2_quality_threshold')
    point = dict(id=row['id'],species=species,mode=row['mode'],rLo=row['rLo'],rHi=row['rHi'],r=(row['rLo']+row['rHi'])/2.,mean=mean,sigma=sigma,meanError=mean_error,sigmaError=sigma_error,ratio=ratio,ratioError=ratio_error,chi2PerNdf=chi2,targetYield=target_yield,effectiveCount=neff,flags=flags,rejectionReasons=sorted(set(reasons)),selected=not reasons,qualityLabel='TOF_TAG_RESPONSE_CANDIDATE_NOT_PURITY_VALIDATION')
    point.update(shapeCovariance=row.get('shapeCovariance',[]),shapeNames=row.get('shapeNames',[]))
    if finite(ratio) and ratio>0 and finite(ratio_error):
        point.update(logRatio=math.log(ratio),logRatioError=ratio_error/ratio)
    return point

def turning_points(coefficients, r_min, r_max):
    derivative = [i*c for i,c in enumerate(coefficients)][1:]
    roots = []
    if len(derivative)==2 and derivative[1] != 0:
        roots = [-derivative[0]/derivative[1]]
    elif len(derivative)==3:
        a,b,c = derivative[2],derivative[1],derivative[0]
        if abs(a)<1e-20:
            if b != 0:
                roots = [-c/b]
        elif b*b-4*a*c>=0:
            disc = math.sqrt(b*b-4*a*c)
            roots = [(-b-disc)/(2*a),(-b+disc)/(2*a)]
    return sorted(10**v for v in roots if math.log10(r_min)<v<math.log10(r_max))

def merge_recenter_audit(point,audit,quality):
    keys=('searchBoundaryMode','searchAxisClipped','fitWindowAxisClipped','supportReviewRequired','nearestOriginalCurveToMode')
    require(all(k in audit for k in keys),'Missing empirical search audit: '+point['id'])
    point['recenterAudit']={k:audit[k] for k in keys}
    for key in keys[:4]:
        require(type(audit[key]) is bool,'Nonboolean search audit')
        if audit[key]:point['flags'].append('recenter_'+key)
    if audit['fitWindowAxisClipped']:point['rejectionReasons'].append('recenter_fit_window_axis_clipped')
    reviewed=quality['reviewedSearchBoundarySupport']
    if point['id'] in reviewed:
        require(isinstance(reviewed[point['id']],str) and reviewed[point['id']].strip(),'Explicit support review reason required')
        point['searchBoundaryReviewReason']=reviewed[point['id']]
    if audit['searchBoundaryMode'] and quality['rejectUnreviewedSearchBoundaryMode'] and point['id'] not in reviewed:
        point['rejectionReasons'].append('recenter_search_boundary_requires_review')
    point['flags']=sorted(set(point['flags']));point['rejectionReasons']=sorted(set(point['rejectionReasons']));point['selected']=not point['rejectionReasons']
    return point

def calibrate_species(points, configuration):
    degrees = configuration['polynomialDegrees']
    require(len(points)>=configuration['minimumSupportPoints'],'too few eligible support points')
    r_min,r_max = min(q['rLo'] for q in points),max(q['rHi'] for q in points)
    joint=None
    if configuration.get('regressionMethod')=='joint_gls':
        from joint_polynomial_gls import fit_joint_polynomials
        joint=fit_joint_polynomials(points,configuration);models=joint['models']
        held_models=[fit_joint_polynomials(points[:i]+points[i+1:],configuration)['models'] for i in range(len(points))]
    else:
        models = {name:regress(points,name,degree,configuration['regressionWeights'],configuration['rankTolerance']) for name,degree in degrees.items()}
    diagnostics,invalid = {},[]
    count = configuration['denseGridPoints']
    for name,model in models.items():
        samples = [polynomial(model['coefficients'],math.log10(r_min+(r_max-r_min)*i/(count-1))) for i in range(count)]
        if not all(finite(v) for v in samples):
            invalid.append(name+'_nonfinite_dense_grid')
        elif name!='logRatio' and min(samples)<=0:
            invalid.append(name+'_nonpositive_dense_grid')
        elif name=='logRatio' and (min(samples)<math.log(configuration['quality']['ratioMin']) or max(samples)>math.log(configuration['quality']['ratioMax'])):
            invalid.append('response_ratio_outside_declared_bounds')
        loo = []
        for i,point in enumerate(points):
            held = held_models[i][name] if joint else regress(points[:i]+points[i+1:],name,degrees[name],configuration['regressionWeights'],configuration['rankTolerance'])
            predicted = sum(w*polynomial(held['coefficients'],math.log10(r)) for r,w in gl3_nodes(point['rLo'],point['rHi']))
            error = model['residuals'][i]['regularizedError']
            loo.append(dict(id=point['id'],prediction=predicted,residual=point[name]-predicted,regularizedPull=(point[name]-predicted)/error))
        diagnostics[name] = dict(turningPointsR=turning_points(model['coefficients'],r_min,r_max),minDenseValue=min(samples),maxDenseValue=max(samples),leaveOneOut=loo,maximumAbsoluteLOOPull=max(abs(q['regularizedPull']) for q in loo),maximumAbsoluteSupportPull=max(abs(q['pull']) for q in model['residuals']))
    ordered = sorted(points,key=lambda q:q['r'])
    gaps = [dict(leftId=a['id'],rightId=b['id'],leftRHi=a['rHi'],rightRLo=b['rLo']) for a,b in zip(ordered,ordered[1:]) if b['rLo']-a['rHi']>configuration['gapWarningWidth']]
    result=dict(curveCandidateValid=not invalid,invalidReasons=invalid,rMin=r_min,rMax=r_max,supportCount=len(points),supportIds=[q['id'] for q in points],models=models,diagnostics=diagnostics,supportGapWarnings=gaps,gapMeaning='Unmeasured interpolation intervals, NOT production vetoes',extrapolationAllowed=False,qualityLabel='CANDIDATE_POLYNOMIAL_NOT_PRODUCTION_VALIDATED')
    if joint:result['jointRegression']=joint['jointRegression']
    return result

def save_json(path,value):
    with path.open('x') as stream:
        json.dump(value,stream,indent=2,ensure_ascii=False,allow_nan=False)
        stream.write('\n')

def write_table(path,header,rows):
    with path.open('x',newline='') as stream:
        writer=csv.writer(stream,delimiter='\t',lineterminator='\n')
        writer.writerow(header)
        for row in rows:
            writer.writerow([format(v,'.17g') if isinstance(v,float) else v for v in row])

def run(fit_dir,settings_path,output_dir):
    configuration=yaml.safe_load(settings_path.read_text())
    require(configuration['schemaVersion']==1,'unsupported settings')
    require(configuration['polynomialDegrees']==dict(mean=3,sigma=2,logRatio=2),'unchanged PID degrees required')
    require(configuration['productionEnabled'] is False,'production writes disabled')
    require(configuration['denseGridPoints']>=3,'dense grid too short')
    source_path=fit_dir/'fits.json'
    source_bytes=source_path.read_bytes()
    source=json.loads(source_bytes)
    require(source['status']=='completed' and source['model']=='split_normal','completed split-normal fits required')
    require(source['rawCountsModified'] is False and source['rawReadback'] is True and source['predictionReadback'] is True,'unverified saved fits')
    require(source['totalSlices']==len(source['results']),'incomplete fits')
    require(len({q['id'] for q in source['results']})==len(source['results']),'duplicate fit ids')
    audit_path=fit_dir/'recenter_provenance.json'
    audit_bytes=audit_path.read_bytes() if audit_path.exists() else None
    audits={q['id']:q for q in json.loads(audit_bytes)['slices']} if audit_bytes else {}
    if audits:require(set(audits)=={q['id'] for q in source['results']},'Search audit IDs differ')
    output_dir.mkdir(parents=True,exist_ok=False)
    snapshot=output_dir/'sources'
    snapshot.mkdir()
    for path in (settings_path,Path(__file__)):
        shutil.copy2(str(path),str(snapshot/path.name))
    scenarios=[dict(name='baseline',chi2PerNdfMax=configuration['quality']['chi2PerNdfMax'],rejectForeignCurveInsideRoi=configuration['quality']['rejectForeignCurveInsideRoi'])]+configuration['sensitivityScenarios']
    summaries=[]
    for scenario in scenarios:
        folder=output_dir/scenario['name']
        folder.mkdir()
        points=[select_point(q,configuration['quality'],scenario['chi2PerNdfMax'],scenario['rejectForeignCurveInsideRoi']) for q in source['results']]
        if audits:points=[merge_recenter_audit(q,audits[q['id']],configuration['quality']) for q in points]
        results={}
        scenario_config=copy.deepcopy(configuration)
        scenario_config['regressionMethod']=scenario.get('regressionMethod','independent_wls')
        for species in SPECIES:
            selected=sorted((q for q in points if q['species']==species and q['selected']),key=lambda q:q['r'])
            try:
                results[species]=calibrate_species(selected,scenario_config)
            except (ValueError,ZeroDivisionError,OverflowError) as error:
                results[species]=dict(curveCandidateValid=False,invalidReasons=[str(error)],supportCount=len(selected),supportIds=[q['id'] for q in selected],extrapolationAllowed=False)
        candidate_rows,valid_rows,curves,responses=[],[],[],[]
        for species in SPECIES:
            record=results[species]
            if 'models' not in record:
                continue
            row=[species,record['rMin'],record['rMax']]+sum((record['models'][name]['coefficients'] for name in ('mean','sigma','logRatio')),[])
            candidate_rows.append(row)
            if record['curveCandidateValid']:
                valid_rows.append(row)
                grid=configuration['denseGridPoints']
                for i in range(grid):
                    r=record['rMin']+(record['rMax']-record['rMin'])*i/(grid-1)
                    x=math.log10(r)
                    curves.append([species,r,polynomial(record['models']['mean']['coefficients'],x),polynomial(record['models']['sigma']['coefficients'],x)])
                    responses.append([species,r,math.exp(polynomial(record['models']['logRatio']['coefficients'],x))])
        write_table(folder/'coefficients.tsv',COEFFICIENT_HEADER,valid_rows)
        write_table(folder/'coefficients_all_candidates.tsv',COEFFICIENT_HEADER,candidate_rows)
        write_table(folder/'curves.tsv',['species','r','mean','sigma'],curves)
        write_table(folder/'response.tsv',['species','r','ratio'],responses)
        save_json(folder/'support_and_rejections.json',dict(scenario=scenario,points=points,rejectionCounts=dict(Counter(r for q in points for r in q['rejectionReasons']))))
        save_json(folder/'calibration.json',dict(schemaVersion=1,scenario=scenario,productionEnabled=False,purityAvailable=False,species=results,coefficientMeaning=dict(mean='mode; cubic log10(R)',sigma='geometric mean sqrt(left*right), not RMS; quadratic log10(R)',logRatio='log(right/left); quadratic response nuisance')))
        yaml_results={}
        for species,record in results.items():
            simple={k:v for k,v in record.items() if k!='models'}
            if 'models' in record:
                simple['coefficients']={name:model['coefficients'] for name,model in record['models'].items()}
            yaml_results[species]=simple
        with (folder/'coefficients.yaml').open('x') as stream:
            yaml.safe_dump(dict(schemaVersion=1,productionEnabled=False,scenario=scenario,basis='powers_of_log10_R',species=yaml_results),stream,sort_keys=False,allow_unicode=True)
        simple_results={s:{k:r[k] for k in ('curveCandidateValid','supportCount','supportIds','invalidReasons','rMin','rMax') if k in r} for s,r in results.items()}
        summaries.append(dict(name=scenario['name'],species=simple_results))
    require(source_path.read_bytes()==source_bytes,'source fits changed')
    if audit_bytes:require(audit_path.read_bytes()==audit_bytes,'Search provenance changed')
    result=dict(schemaVersion=1,status='completed',sourceFitPath=str(source_path.resolve()),sourceFitSha256=hashlib.sha256(source_bytes).hexdigest(),settingsSha256=hashlib.sha256(settings_path.read_bytes()).hexdigest(),sourceSha256=hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),productionEnabled=False,purityAvailable=False,rawInputModified=False,storageMappingIsHypothesis=source.get('storageMappingIsHypothesis'),onFileEncodingAsserted=source.get('onFileEncodingAsserted'),meanFunctionUnchanged=True,sigmaFunctionUnchanged=True,extrapolationAllowed=False,scenarios=summaries)
    result['searchProvenanceSha256']=hashlib.sha256(audit_bytes).hexdigest() if audit_bytes else None
    save_json(output_dir/'completed.json',result)
    print('PID_POLYNOMIAL_CALIBRATION_COMPLETE '+str(output_dir),flush=True)
    for item in summaries:
        print(item['name'],{s:(q['supportCount'],q['curveCandidateValid'],q['invalidReasons']) for s,q in item['species'].items()},flush=True)
    return result

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--fit-dir',type=Path,required=True)
    parser.add_argument('--settings',type=Path,required=True)
    parser.add_argument('--output-dir',type=Path,required=True)
    args=parser.parse_args()
    run(args.fit_dir.resolve(),args.settings.resolve(),args.output_dir.resolve())

if __name__=='__main__':
    main()
