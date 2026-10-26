"""Small weighted-QR polynomials; no NumPy or change to the PID formula family."""
import itertools
import math


def value(c, x):
    result = 0.
    for a in reversed(c):
        result = result*x+a
    return result


def solve(matrix, rhs, tolerance=1.e-12):
    n = len(rhs)
    a = [list(row)+[b] for row, b in zip(matrix, rhs)]
    if any(len(row) != n+1 for row in a):
        raise ValueError('Nonsquare system')
    for j in range(n):
        pivot = max(range(j, n), key=lambda k: abs(a[k][j]))
        if abs(a[pivot][j]) < tolerance:
            raise ValueError('Singular constrained polynomial system')
        a[j], a[pivot] = a[pivot], a[j]
        scale = a[j][j]
        a[j] = [v/scale for v in a[j]]
        for k in range(n):
            if k != j:
                scale = a[k][j]
                a[k] = [b-scale*v for b, v in zip(a[k], a[j])]
    return [row[-1] for row in a]


def weighted_qr(x, y, weights, degree, rank_tolerance=1.e-12):
    if len(x) != len(y) or len(x) != len(weights) or len(x) <= degree:
        raise ValueError('Insufficient matching polynomial inputs')
    if any(not math.isfinite(v) for row in (x, y, weights) for v in row):
        raise ValueError('Nonfinite polynomial input')
    if any(w < 0 for w in weights) or not max(weights) > 0:
        raise ValueError('Weights must be nonnegative and nonempty')
    # Weight normalization does not change the minimizer.
    roots = [math.sqrt(w/max(weights)) for w in weights]
    q, upper = [], [[0.]*(degree+1) for _ in range(degree+1)]
    for j in range(degree+1):
        column = [root*v**j for root, v in zip(roots, x)]
        original_norm = math.sqrt(math.fsum(v*v for v in column))
        for _ in range(2):
            for k, basis in enumerate(q):
                dot = math.fsum(a*b for a, b in zip(basis, column))
                upper[k][j] += dot
                column = [a-dot*b for a, b in zip(column, basis)]
        norm = math.sqrt(math.fsum(v*v for v in column))
        if norm <= rank_tolerance*max(1., original_norm):
            raise ValueError('Rank-deficient weighted polynomial')
        upper[j][j] = norm
        q.append([v/norm for v in column])
    rhs = [math.fsum(a*b*r for a, b, r in zip(column, y, roots)) for column in q]
    c = [0.]*(degree+1)
    for j in range(degree, -1, -1):
        c[j] = (rhs[j]-math.fsum(upper[j][k]*c[k] for k in range(j+1, degree+1)))/upper[j][j]
    return c


def normalized_basis(lo, hi):
    if not 0 < lo < hi:
        raise ValueError('Invalid rigidity interval')
    a, b = math.log10(lo), math.log10(hi)
    return (a+b)/2., (b-a)/2.


def to_original(c, center, scale):
    """p((log10(R)-center)/scale) -> original ascending log10(R) coefficients."""
    if not scale > 0:
        raise ValueError('Nonpositive basis scale')
    output = [0.]*len(c)
    for j, coefficient in enumerate(c):
        for k in range(j+1):
            output[k] += coefficient*math.comb(j, k)*(-center)**(j-k)/scale**j
    return output


def quadratic_minimum(c):
    points = [-1., 1.]
    if c[2] > 0 and -1 < -c[1]/(2*c[2]) < 1:
        points.append(-c[1]/(2*c[2]))
    point = min(points, key=lambda t: value(c, t))
    return point, value(c, point)


def positive_quadratic(x, y, weights, minimum, tolerance=1.e-9, maximum_iterations=50):
    """Convex WLS with exact endpoint/vertex positivity, not Bernstein NNLS.

    A cutting-plane active-set QP adds the current violating quadratic vertex.
    Every retained constraint is a true pointwise positivity constraint.  Final
    acceptance uses the analytic global minimum, so sampling cannot hide a dip.
    """
    if not minimum > 0 or tolerance <= 0:
        raise ValueError('Invalid positive-width policy')
    free = weighted_qr(x, y, weights, 2)
    point, low = quadratic_minimum(free)
    if low >= minimum:
        return free, dict(constrained=False, minimum=low, minimumAt=point, iterations=0)
    weight_scale = max(weights)
    w = [v/weight_scale for v in weights]
    h = [[math.fsum(a*t**(j+k) for a, t in zip(w, x)) for k in range(3)] for j in range(3)]
    f = [math.fsum(a*b*t**j for a, b, t in zip(w, y, x)) for j in range(3)]
    constraints = [-1., 1.]
    for iteration in range(maximum_iterations):
        best, best_loss = None, float('inf')
        for count in range(4):
            for active in itertools.combinations(constraints, count):
                rows = [[1., t, t*t] for t in active]
                # H*c - A^T*lambda=f, A*c=minimum; lambda>=0.
                matrix = [h[j]+[-row[j] for row in rows] for j in range(3)]
                matrix += [row+[0.]*count for row in rows]
                try:
                    fitted = solve(matrix, f+[minimum]*count)
                except ValueError:
                    continue
                c, multipliers = fitted[:3], fitted[3:]
                if any(v < -tolerance for v in multipliers):
                    continue
                if any(value(c, t) < minimum-tolerance for t in constraints):
                    continue
                loss = math.fsum(a*(value(c, t)-b)**2 for a, b, t in zip(w, y, x))
                if loss < best_loss:
                    best, best_loss = c, loss
        if best is None:
            raise ValueError('No feasible positive-width QP solution')
        point, low = quadratic_minimum(best)
        if low >= minimum-tolerance:
            # Numerical roundoff correction only, bounded by tolerance.
            correction = max(0., minimum-low)
            best[0] += correction
            return best, dict(constrained=True, minimum=low+correction, minimumAt=point,
                              iterations=iteration+1, roundoffCorrection=correction)
        if any(abs(point-t) <= tolerance for t in constraints):
            raise ValueError('Positive-width QP stagnated at a violating vertex')
        constraints.append(point)
    raise ValueError('Positive-width polynomial did not converge')


def smooth_species(rows, species, cfg, domain):
    center, scale = normalized_basis(*domain)
    x, mean, sigma, weights = [], [], [], []
    for row in rows:
        item = row[species]
        weight = item['effectiveCount']/max(item['penaltyScale'], cfg['minimumPenaltyScale'])**2
        if not row['finite'] or item['windowYield'] <= cfg['minimumWindowYield'] or weight <= 0:
            continue
        x.append((math.log10(row['r'])-center)/scale)
        mean.append(item['mean']); sigma.append(item['sigma']); weights.append(weight)
    if len(x) < cfg['minimumSmoothingSlices']:
        raise ValueError('Too few finite supported slices for '+species)
    cm = weighted_qr(x, mean, weights, 3, cfg['rankTolerance'])
    cs, audit = positive_quadratic(x, sigma, weights, cfg['minimumSigma'],
                                  cfg['positivityTolerance'], cfg['positivityMaximumIterations'])
    return dict(meanCoefficients=to_original(cm, center, scale),
                sigmaCoefficients=to_original(cs, center, scale),
                normalizedMeanCoefficients=cm, normalizedSigmaCoefficients=cs,
                basisCenter=center, basisScale=scale, widthPositivity=audit,
                nSlices=len(x), weightedMeanResidual=math.fsum(w*(value(cm,t)-y)**2 for w,t,y in zip(weights,x,mean)),
                weightedSigmaResidual=math.fsum(w*(value(cs,t)-y)**2 for w,t,y in zip(weights,x,sigma)))
