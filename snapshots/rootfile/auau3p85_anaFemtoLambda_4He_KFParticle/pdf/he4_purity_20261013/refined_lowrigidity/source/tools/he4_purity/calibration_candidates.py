"""Publish study-only fitted He4 locus parameters without changing production PID."""
import csv
import json
import math
from pathlib import Path

import yaml


def polynomial_value(coefficients, x):
    value=0.
    for coefficient in reversed(coefficients):
        value=value*x+coefficient
    return value


def polynomial_fit(x, y, degree):
    """Small unweighted QR fit; avoids a host-dependent NumPy/Fortran runtime."""
    if isinstance(degree, bool) or not isinstance(degree, int) or degree<0 or len(x)!=len(y) or len(x)<=degree:
        raise ValueError('Invalid polynomial degree or candidate knots.')
    columns=[]
    upper=[[0.]*(degree+1) for _ in range(degree+1)]
    for k in range(degree+1):
        values=[float(v)**k for v in x]
        # Reorthogonalize to keep the small log10(R) polynomial system stable.
        for _ in range(2):
            for j, column in enumerate(columns):
                projection=math.fsum(a*b for a,b in zip(column, values))
                upper[j][k]+=projection
                values=[a-projection*b for a,b in zip(values, column)]
        norm=math.sqrt(math.fsum(v*v for v in values))
        if norm<=1.e-14*math.sqrt(len(x)):
            raise ValueError('Candidate polynomial fit is rank deficient.')
        upper[k][k]=norm
        columns.append([v/norm for v in values])
    rhs=[math.fsum(a*b for a,b in zip(column, y)) for column in columns]
    coefficients=[0.]*(degree+1)
    for k in range(degree, -1, -1):
        coefficients[k]=(rhs[k]-math.fsum(upper[k][j]*coefficients[j] for j in range(k+1, degree+1)))/upper[k][k]
    return coefficients


def write_candidate_tables(output, settings, sample):
    output = Path(output)
    data = json.loads((output/'purity.json').read_text())
    destination = output/'he4_locus_fit_candidates.csv'
    fields = ['rigidity_low', 'rigidity_high', 'rigidity_center', 'model', 'chosen',
              'fit_valid', 'fit_quality', 'signal_boundary', 'nuisance_boundary',
              'raw_peak_search_applied', 'raw_local_maximum', 'raw_peak_prominence',
              'signal_component_evidence', 'signal_amplitude_over_error', 'null_delta_aic',
              'null_fit_valid', 'initial_mean_nsigma', 'fit_mean_nsigma', 'fit_sigma_nsigma',
              'original_mean_dedx', 'original_sigma_dedx', 'fitted_mean_dedx', 'fitted_sigma_dedx',
              'fitted_mean_stat_error', 'fitted_sigma_stat_error',
              'original_gate_low_dedx', 'original_gate_high_dedx',
              'proposed_gate_low_dedx', 'proposed_gate_high_dedx',
              'original_gate_purity', 'original_gate_purity_stat_error',
              'proposed_gate_within_fit_range', 'proposed_gate_purity', 'proposed_gate_purity_stat_error',
              'original_gate_model_signal', 'proposed_gate_model_signal',
              'original_to_proposed_model_signal_ratio', 'chi2_ndf']
    knots = []
    gate = float(data['gateNSigma'])
    with destination.open('w', newline='') as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for s in data['slices']:
            for index, f in enumerate(s['fits']):
                p = f['parameters']
                errors = f['parameterErrors']
                if len(p)<3 or len(errors)<3:
                    writer.writerow(dict(rigidity_low=s['lo'], rigidity_high=s['hi'], rigidity_center=s['r'],
                        model=f['model'], chosen=False, fit_valid=False, fit_quality=False))
                    continue
                proposed_signal = (p[0]*p[2]*math.sqrt(2*math.pi)*math.erf(gate/math.sqrt(2)) /
                                   (data['fitDedxBinWidth']/s['calibrationSigma']))
                within = f['fittedGateWithinFitRange']
                row = dict(rigidity_low=s['lo'], rigidity_high=s['hi'], rigidity_center=s['r'],
                    model=f['model'], chosen=index==s['chosen'], fit_valid=f['valid'], fit_quality=f['quality'],
                    signal_boundary=f.get('signalAtBoundary', f['atBoundary']),
                    nuisance_boundary=f.get('nuisanceAtBoundary', f['atBoundary']),
                    raw_peak_search_applied=f['searchApplied'], raw_local_maximum=f['rawLocalMaximum'],
                    raw_peak_prominence=f['rawPeakProminence'], signal_component_evidence=f['signalEvidence'],
                    signal_amplitude_over_error=f['signalAmplitudeOverError'], null_delta_aic=f['nullDeltaAic'],
                    null_fit_valid=f['nullFitValid'], initial_mean_nsigma=f['initialMean'],
                    fit_mean_nsigma=p[1], fit_sigma_nsigma=p[2], original_mean_dedx=s['calibrationMean'],
                    original_sigma_dedx=s['calibrationSigma'], fitted_mean_dedx=f['rawSignalMean'],
                    fitted_sigma_dedx=f['rawSignalSigma'], fitted_mean_stat_error=errors[1]*s['calibrationSigma'],
                    fitted_sigma_stat_error=errors[2]*s['calibrationSigma'],
                    original_gate_low_dedx=s['calibrationMean']-gate*s['calibrationSigma'],
                    original_gate_high_dedx=s['calibrationMean']+gate*s['calibrationSigma'],
                    proposed_gate_low_dedx=f['rawSignalMean']-gate*f['rawSignalSigma'],
                    proposed_gate_high_dedx=f['rawSignalMean']+gate*f['rawSignalSigma'],
                    original_gate_purity=f['purity'], original_gate_purity_stat_error=f['purityError'],
                    proposed_gate_within_fit_range=within,
                    proposed_gate_purity=f['purityFittedMeanSigmaGate'] if within else '',
                    proposed_gate_purity_stat_error=f['purityFittedGateError'] if within else '',
                    original_gate_model_signal=f['gateSignal'],
                    proposed_gate_model_signal=proposed_signal if within else '',
                    original_to_proposed_model_signal_ratio=f['gateSignal']/proposed_signal if within and proposed_signal>0 else '',
                    chi2_ndf=f['chi2']/f['ndf'] if f['ndf']>0 else '')
                writer.writerow(row)
                if row['chosen'] and f['valid'] and f['signalEvidence']:
                    knots.append(dict(rigidity=s['r'], mean_dedx=f['rawSignalMean'], sigma_dedx=f['rawSignalSigma'],
                        background_model=f['model'], fit_quality=f['quality'],
                        raw_peak_searched=f['searchApplied'], raw_local_maximum=f['rawLocalMaximum'],
                        nuisance_boundary=f.get('nuisanceAtBoundary', f['atBoundary']),
                        mean_stat_error=errors[1]*s['calibrationSigma'], sigma_stat_error=errors[2]*s['calibrationSigma']))
    low, high = settings['calibrationCandidateMinRigidity'], settings['calibrationCandidateMaxRigidity']
    chosen = [k for k in knots if low <= k['rigidity'] < high]
    candidate = dict(schemaVersion=1, sample=sample, inputRoot=data['inputRoot'],
        coordinates=dict(rigidity='|p/q| in GeV/c', mean='raw dE/dx in keV/cm', sigma='Gaussian raw dE/dx width in keV/cm'),
        candidateOnly=True, productionEnabled=False, originalGateUnchanged=True,
        useDomain=[low, high],
        interpretation='Empirical component-locus proposal, not a validated change to the production PID.',
        limitations=['Background-model dependence and nuisance boundaries remain.',
                     'Conditional amplitude/error and null-model AIC are model evidence, not a formal discovery significance.',
                     'A shoulder can have a fitted component without a separate raw peak.',
                     'Statistical fit errors do not include model uncertainty or calibration uncertainty.'],
        knots=knots)
    candidate['observedKnotDomain'] = [min(k['rigidity'] for k in chosen), max(k['rigidity'] for k in chosen)] if chosen else []
    candidate['extrapolationPolicy'] = 'Requested useDomain does not establish calibration outside observedKnotDomain; no production extrapolation is approved.'
    if len(chosen) >= settings['calibrationCandidateMinimumKnots']:
        x = [math.log10(k['rigidity']) for k in chosen]
        for coordinate, degree in [('mean', settings['calibrationCandidateMeanDegree']),
                                   ('sigma', settings['calibrationCandidateSigmaDegree'])]:
            y = [k[coordinate+'_dedx'] for k in chosen]
            coefficients = polynomial_fit(x, y, degree)
            residual = [v-polynomial_value(coefficients,u) for u,v in zip(x,y)]
            candidate[coordinate+'LogPolyCandidate'] = dict(
                coefficients=[float(v) for v in coefficients],
                formula='sum(c[i] * log10(R)^i)', fitPolicy='unweighted provisional fit to numerically valid component-evidence knots',
                nKnots=len(chosen), rms_dedx=math.sqrt(math.fsum(v*v for v in residual)/len(residual)))
    else:
        candidate['polynomialStatus'] = 'insufficient component-evidence knots'
    (output/('he4_pid_curve_candidate_'+sample+'.yaml')).write_text(
        yaml.safe_dump(candidate, sort_keys=False, allow_unicode=True))
    return dict(csv=destination.name, yaml='he4_pid_curve_candidate_'+sample+'.yaml',
                componentKnots=len(knots), polynomialKnots=len(chosen), productionEnabled=False)
