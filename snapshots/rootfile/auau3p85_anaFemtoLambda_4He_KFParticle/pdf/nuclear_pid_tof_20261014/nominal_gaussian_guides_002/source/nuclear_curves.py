"""Study-only smooth curves. No runtime PID files are modified by this module."""
import bisect
import math

from runtime_selection import polynomial


def weighted_polynomial(x, y, error, degree):
    x, y, error = (list(map(float, a)) for a in (x, y, error))
    if (len(x) != len(y) or len(x) != len(error) or len(x) < degree+1
            or not all(math.isfinite(v) for a in (x, y, error) for v in a)
            or any(e <= 0 for e in error)):
        raise ValueError('Insufficient finite points for weighted polynomial.')
    # Small-degree QR with reorthogonalisation and scaled x; no BLAS/runtime dependency.
    center = (min(x)+max(x))/2
    scale = max(abs(v-center) for v in x)
    if scale == 0:
        if degree:
            raise ValueError('Rank-deficient polynomial.')
        scale = 1.
    columns = [[((v-center)/scale)**k/e for v, e in zip(x, error)] for k in range(degree+1)]
    q, rr = [], [[0.]*(degree+1) for _ in range(degree+1)]
    for j, column in enumerate(columns):
        v = list(column)
        for _ in range(2):
            for i in range(j):
                projection = math.fsum(a*b for a, b in zip(q[i], v))
                rr[i][j] += projection
                v = [a-projection*b for a, b in zip(v, q[i])]
        norm = math.sqrt(math.fsum(a*a for a in v))
        if norm <= 1.e-12*math.sqrt(math.fsum(a*a for a in column)):
            raise ValueError('Rank-deficient polynomial.')
        rr[j][j] = norm
        q.append([a/norm for a in v])
    rhs = [math.fsum(a*b/e for a, b, e in zip(column, y, error)) for column in q]
    scaled_coeff = [0.]*(degree+1)
    for j in reversed(range(degree+1)):
        scaled_coeff[j] = (rhs[j]-math.fsum(rr[j][k]*scaled_coeff[k]
                                          for k in range(j+1, degree+1)))/rr[j][j]
    coefficients = [0.]*(degree+1)
    for k, value in enumerate(scaled_coeff):
        for j in range(k+1):
            coefficients[j] += value*math.comb(k, j)*(-center)**(k-j)/scale**k
    return coefficients


def pchip_slopes(x, y):
    """Fritsch-Carlson shape-preserving slopes, including limited end slopes."""
    if len(x) != len(y) or len(x) < 2 or any(b <= a for a, b in zip(x, x[1:])):
        raise ValueError('PCHIP requires at least two strictly ordered knots.')
    h = [b-a for a, b in zip(x, x[1:])]
    delta = [(b-a)/width for a, b, width in zip(y, y[1:], h)]
    slopes = [0.]*len(x)
    if len(x) == 2:
        return [float(delta[0])]*2
    for i in range(1, len(x)-1):
        if delta[i-1]*delta[i] > 0:
            w1, w2 = 2*h[i]+h[i-1], h[i]+2*h[i-1]
            slopes[i] = (w1+w2)/(w1/delta[i-1]+w2/delta[i])
    for target, h0, h1, d0, d1 in ((0, h[0], h[1], delta[0], delta[1]),
                                   (-1, h[-1], h[-2], delta[-1], delta[-2])):
        value = ((2*h0+h1)*d0-h0*d1)/(h0+h1)
        if value*d0 <= 0:
            value = 0.
        elif d0*d1 < 0 and abs(value) > 3*abs(d0):
            value = 3*d0
        slopes[target] = value
    return slopes


def pchip_value(x, y, slopes, at):
    # Only the explicit native half-bin endpoint extension is used by a segment.
    if at <= x[0]:
        return y[0]+slopes[0]*(at-x[0])
    if at >= x[-1]:
        return y[-1]+slopes[-1]*(at-x[-1])
    i = bisect.bisect_right(x, at)-1
    width = x[i+1]-x[i]
    t = (at-x[i])/width
    return ((2*t**3-3*t*t+1)*y[i]+(t**3-2*t*t+t)*width*slopes[i]
            +(-2*t**3+3*t*t)*y[i+1]+(t**3-t*t)*width*slopes[i+1])


def local_smooth(x, y, errors, neighbours):
    result = []
    for at in x:
        indices = sorted(range(len(x)), key=lambda i: abs(x[i]-at))[:min(neighbours, len(x))]
        xx = [x[i]-at for i in indices]
        coeff = weighted_polynomial(xx, [y[i] for i in indices],
                                    [errors[i] for i in indices], min(2, len(indices)-1))
        result.append(coeff[0])
    return result


def evaluate_segment(segment, rigidity):
    x = math.log10(rigidity)
    if segment['family'] == 'logpoly':
        mean = polynomial(segment['coeffMean'], x)
        sigma = polynomial(segment['coeffSigma'], x)
    elif segment['family'] == 'pchip_logR':
        knots = segment['knots']
        positions = [math.log10(k['r']) for k in knots]
        mean = pchip_value(positions, [k['mean'] for k in knots], segment['meanSlopes'], x)
        sigma = math.exp(pchip_value(positions, [k['logSigma'] for k in knots],
                                    segment['logSigmaSlopes'], x))
    else:
        raise ValueError('Unknown study curve family.')
    uncertainties = segment['uncertaintyKnots']
    positions = [k['r'] for k in uncertainties]
    def interpolate(field):
        i = bisect.bisect_right(positions, rigidity)
        if i == 0:
            return uncertainties[0][field]
        if i == len(positions):
            return uncertainties[-1][field]
        f = (rigidity-positions[i-1])/(positions[i]-positions[i-1])
        return (1-f)*uncertainties[i-1][field]+f*uncertainties[i][field]
    mean_error, sigma_error = interpolate('meanError'), interpolate('sigmaError')
    return dict(mean=float(mean), sigma=float(sigma), meanError=mean_error,
                sigmaError=sigma_error, supported=bool(mean > 0 and sigma > 0
                and math.isfinite(mean) and math.isfinite(sigma)), identityValidated=False)


def evaluate_curve(candidate, species, rigidity):
    if candidate.get('schemaVersion') != 1 or candidate.get('productionEnabled') is not False:
        raise ValueError('Only explicit study-only schemaVersion=1 is accepted.')
    if not math.isfinite(rigidity) or rigidity <= 0:
        return dict(mean=0., sigma=0., meanError=0., sigmaError=0., supported=False,
                    identityValidated=False)
    for segment in candidate['species'].get(species, {}).get('segments', []):
        if segment['domain'][0] <= rigidity <= segment['domain'][1]:
            result = evaluate_segment(segment, rigidity)
            if 'supportedNativeIntervals' in segment:
                result['supported'] = result['supported'] and any(
                    lo <= rigidity <= hi for lo, hi in segment['supportedNativeIntervals'])
            return result
    return dict(mean=0., sigma=0., meanError=0., sigmaError=0., supported=False,
                identityValidated=False)


def make_segment(points, family, settings):
    x = [math.log10(p['r']) for p in points]
    means = [p['mean'] for p in points]
    sigmas = [p['intrinsicSigmaSeed'] for p in points]
    mean_errors = [max(p['meanError'], settings['relativeMeanWeightFloorOverSigma']*s)
                   for p, s in zip(points, sigmas)]
    sigma_errors = [max(p['sigmaError'], settings['relativeSigmaWeightFloor']*s)
                    for p, s in zip(points, sigmas)]
    segment = dict(family=family, domain=[points[0]['lo'], points[-1]['hi']],
                   uncertaintyKnots=[dict(r=p['r'], meanError=e, sigmaError=f)
                   for p, e, f in zip(points, mean_errors, sigma_errors)],
                   endpointExtension='native half-bin only; tangent continuation for PCHIP',
                   widthMeaning='R-variance-subtracted seed; requires R-averaged spectrum validation')
    if family == 'logpoly':
        segment['coeffMean'] = weighted_polynomial(x, means, mean_errors, settings['meanDegree'])
        segment['coeffSigma'] = weighted_polynomial(x, sigmas, sigma_errors, settings['sigmaDegree'])
    elif family == 'pchip_logR':
        smooth_mean = local_smooth(x, means, mean_errors, settings['smoothingNeighbours'])
        logs = [math.log(s) for s in sigmas]
        smooth_log_sigma = local_smooth(x, logs, [e/s for e, s in zip(sigma_errors, sigmas)],
                                       settings['smoothingNeighbours'])
        segment['knots'] = [dict(r=p['r'], mean=m, logSigma=ls)
                            for p, m, ls in zip(points, smooth_mean, smooth_log_sigma)]
        segment['meanSlopes'] = pchip_slopes(x, smooth_mean)
        segment['logSigmaSlopes'] = pchip_slopes(x, smooth_log_sigma)
    else:
        raise ValueError('Unknown family.')
    lo, hi = segment['domain']
    grid = [lo+(hi-lo)*i/1000 for i in range(1001)]
    if not all(evaluate_segment(segment, float(r))['supported'] for r in grid):
        raise ValueError('Candidate mean or width is nonpositive within support.')
    residuals = [evaluate_segment(segment, p['r']) for p in points]
    segment['residualMeanRmsOverSigma'] = math.sqrt(math.fsum(
        ((v['mean']-p['mean'])/p['intrinsicSigmaSeed'])**2 for p, v in zip(points, residuals))/len(points))
    segment['residualSigmaRelativeRms'] = math.sqrt(math.fsum(
        ((v['sigma']-p['intrinsicSigmaSeed'])/p['intrinsicSigmaSeed'])**2
        for p, v in zip(points, residuals))/len(points))
    return segment
