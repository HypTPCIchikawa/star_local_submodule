"""First-order, shared-input errors for mass-window correlation functions.

Each input histogram bin is a random variable with the supplied Sumw2 variance.
Different input keys are independent unless a covariance is explicitly supplied.
Build the full window by adding its disjoint halves, and reuse the same sideband
objects, rather than registering their sums as independent input histograms.

This model preserves correlations caused by common input bins and normalization
integrals. It cannot recover event-level or mixed-pair correlations not stored in
the histograms. Its errors and contrasts are approximate, not discovery tests.
There are no analysis cut, binning, or normalization-range defaults in this file.
"""

import math


class Quantity:
    """A scalar value and its sparse derivative with respect to input bins."""

    def __init__(self, model, value, gradient=None):
        self.model = model
        self.value = float(value)
        if not math.isfinite(self.value):
            raise ValueError("non-finite quantity")
        self.gradient = {key: float(derivative) for key, derivative in
                         (gradient or {}).items() if derivative != 0.0}
        if any(not math.isfinite(v) for v in self.gradient.values()):
            raise ValueError("non-finite derivative")

    def _coerce(self, other):
        if not isinstance(other, Quantity):
            return self.model.constant(other)
        if other.model is not self.model:
            raise ValueError("cannot combine quantities from different models")
        return other

    def _combine(self, other, value, left_factor, right_factor):
        gradient = {key: left_factor * derivative
                    for key, derivative in self.gradient.items()}
        for key, derivative in other.gradient.items():
            gradient[key] = gradient.get(key, 0.0) + right_factor * derivative
        return Quantity(self.model, value, gradient)

    def __add__(self, other):
        other = self._coerce(other)
        return self._combine(other, self.value + other.value, 1.0, 1.0)

    __radd__ = __add__

    def __neg__(self):
        return Quantity(self.model, -self.value,
                        {key: -v for key, v in self.gradient.items()})

    def __sub__(self, other):
        return self + -self._coerce(other)

    def __rsub__(self, other):
        return self._coerce(other) - self

    def __mul__(self, other):
        other = self._coerce(other)
        return self._combine(other, self.value * other.value,
                             other.value, self.value)

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = self._coerce(other)
        if other.value == 0.0:
            raise ZeroDivisionError("zero-valued denominator")
        return self._combine(other, self.value / other.value,
                             1.0 / other.value,
                             -self.value / (other.value * other.value))

    def __rtruediv__(self, other):
        return self._coerce(other) / self


class ErrorModel:
    """Named random inputs and their variances / optional covariances."""

    def __init__(self):
        self.variances = {}
        self.cross_covariances = {}

    def constant(self, value):
        return Quantity(self, value)

    def variable(self, key, value, variance):
        if key in self.variances:
            raise ValueError("duplicate random-variable key: %r" % (key,))
        variance = float(variance)
        if not math.isfinite(variance) or variance < 0.0:
            raise ValueError("input variance must be finite and nonnegative")
        quantity = Quantity(self, value, {key: 1.0})
        self.variances[key] = variance
        return quantity

    def series(self, prefix, contents, sumw2):
        if len(contents) != len(sumw2):
            raise ValueError("contents and Sumw2 sizes differ")
        return [self.variable((prefix, index), value, variance)
                for index, (value, variance) in enumerate(zip(contents, sumw2))]

    def correlated_block(self, prefix, values, covariance):
        """Register a block with an externally supplied covariance matrix.

        For example, use a 2x2 shared-polynomial background covariance for the
        low/high mass integrals. This does not create missing fit/data or
        neighboring sliding-window covariances. Such approximations must be
        recorded by the caller.
        """
        values = list(values)
        covariance = [list(row) for row in covariance]
        size = len(values)
        if not size or len(covariance) != size or any(len(row) != size for row in covariance):
            raise ValueError("covariance block shape differs from values")
        matrix = [[float(value) for value in row] for row in covariance]
        if any(not math.isfinite(value) for row in matrix for value in row):
            raise ValueError("non-finite covariance block")
        if any(matrix[index][index] < 0.0 for index in range(size)):
            raise ValueError("negative block variance")
        tolerance = 1e-12 * max((abs(value) for row in matrix for value in row), default=0.0)
        if any(abs(matrix[i][j] - matrix[j][i]) > tolerance
               for i in range(size) for j in range(size)):
            raise ValueError("covariance matrix is not symmetric")
        # Semidefinite LDL decomposition, also allowing exactly singular blocks.
        lower = [[0.0] * size for _ in range(size)]
        diagonal = [0.0] * size
        for i in range(size):
            lower[i][i] = 1.0
            for j in range(i):
                residual = matrix[i][j] - math.fsum(
                    lower[i][k] * lower[j][k] * diagonal[k] for k in range(j))
                if abs(diagonal[j]) <= tolerance:
                    if abs(residual) > tolerance:
                        raise ValueError("covariance matrix is not positive semidefinite")
                else:
                    lower[i][j] = residual / diagonal[j]
            diagonal[i] = matrix[i][i] - math.fsum(
                lower[i][k] ** 2 * diagonal[k] for k in range(i))
            if diagonal[i] < -tolerance:
                raise ValueError("covariance matrix is not positive semidefinite")
        if any((prefix, index) in self.variances for index in range(size)):
            raise ValueError("correlated block keys already exist")
        if any(not math.isfinite(float(value)) for value in values):
            raise ValueError("non-finite block value")
        output = self.series(prefix, values, [matrix[i][i] for i in range(size)])
        for i in range(size):
            for j in range(i):
                self.set_covariance((prefix, i), (prefix, j), matrix[i][j])
        return output

    def set_covariance(self, first_key, second_key, covariance):
        """Register an externally known covariance; does not infer fit errors.

        This checks pairwise Cauchy--Schwarz consistency, not positive
        semidefiniteness of an arbitrary multi-input covariance matrix.
        """
        if first_key not in self.variances or second_key not in self.variances:
            raise KeyError("covariance references an unregistered variable")
        covariance = float(covariance)
        if not math.isfinite(covariance):
            raise ValueError("non-finite covariance")
        if first_key == second_key:
            if covariance != self.variances[first_key]:
                raise ValueError("diagonal covariance must equal input variance")
            return
        limit = math.sqrt(self.variances[first_key] * self.variances[second_key])
        if abs(covariance) > limit * (1.0 + 1e-12):
            raise ValueError("covariance violates Cauchy--Schwarz bound")
        key = frozenset((first_key, second_key))
        if covariance == 0.0:
            self.cross_covariances.pop(key, None)
        else:
            self.cross_covariances[key] = covariance

    def covariance(self, first, second):
        if first.model is not self or second.model is not self:
            raise ValueError("quantity belongs to another error model")
        terms = [derivative * second.gradient.get(key, 0.0) * self.variances[key]
                 for key, derivative in first.gradient.items()]
        for keys, covariance in self.cross_covariances.items():
            key1, key2 = tuple(keys)
            terms.append(covariance * (
                first.gradient.get(key1, 0.0) * second.gradient.get(key2, 0.0) +
                first.gradient.get(key2, 0.0) * second.gradient.get(key1, 0.0)))
        return math.fsum(terms)

    def variance(self, quantity):
        variance = self.covariance(quantity, quantity)
        if variance < 0.0:
            scale = math.fsum(derivative * derivative * self.variances[key]
                              for key, derivative in quantity.gradient.items())
            if variance < -1e-12 * scale:
                raise ValueError("negative propagated variance; invalid covariance matrix")
            return 0.0
        return variance

    def estimate(self, quantity):
        variance = self.variance(quantity)
        return dict(value=quantity.value, variance=variance,
                    error=math.sqrt(variance))

    def covariance_matrix(self, quantities):
        return [[self.covariance(left, right) for right in quantities]
                for left in quantities]


def _same_length(first, second):
    if len(first) != len(second):
        raise ValueError("series sizes differ")


def series_add(first, second):
    _same_length(first, second)
    return [left + right for left, right in zip(first, second)]


def series_subtract(signal, background, scale):
    """Subtract a scaled common background without losing its shared errors."""
    _same_length(signal, background)
    return [observed - scale * sideband
            for observed, sideband in zip(signal, background)]


def integral(series, indices):
    """Sum explicitly selected bins, including their derivative information."""
    indices = list(indices)
    if not indices:
        raise ValueError("integral requires at least one bin")
    if len(set(indices)) != len(indices):
        raise ValueError("duplicate bin indices")
    if any(not isinstance(index, int) or index < 0 or index >= len(series)
           for index in indices):
        raise IndexError("bin index out of range")
    return sum((series[index] for index in indices),
               series[indices[0]].model.constant(0.0))


def normalized_cf(numerator, denominator, norm_indices):
    """C_i = (B_norm/A_norm) A_i/B_i; nonpositive ME bins are masked.

    A zero SE bin is retained (finite derivative, no division by SE). Negative
    subtracted SE values are also retained and flagged, not silently clipped.
    Delta-method errors of sparse bins are not Poisson confidence intervals.
    Both normalization integrals must be strictly positive. Shared bins in the
    plotted and normalization intervals are automatically counted only once.
    """
    _same_length(numerator, denominator)
    norm_indices = list(norm_indices)
    anorm = integral(numerator, norm_indices)
    bnorm = integral(denominator, norm_indices)
    output = dict(bins=[None] * len(numerator), normalization=None,
                  numerator_norm=anorm, denominator_norm=bnorm,
                  invalid={}, warnings={})
    if anorm.value <= 0.0 or bnorm.value <= 0.0:
        output['invalid'] = {index: 'nonpositive normalization integral'
                             for index in range(len(numerator))}
        return output
    factor = bnorm / anorm
    output['normalization'] = factor
    for index, (same, mixed) in enumerate(zip(numerator, denominator)):
        if mixed.value <= 0.0:
            output['invalid'][index] = 'nonpositive mixed-event denominator'
            continue
        output['bins'][index] = factor * same / mixed
        if same.value <= 0.0:
            output['warnings'][index] = 'nonpositive SE; Gaussian delta error is not a sparse-count interval'
    return output


def purity_corrected_cf(signal_cf, background_cf, purity):
    """(C_signal - (1-P) C_background)/P, conditional on fixed purity P.

    Quantity inputs preserve common normalization / sideband bin covariance.
    The uncertainty of the inclusive purity and its correlation with pair
    histograms are intentionally not added by this conditional estimator.
    """
    purity = float(purity)
    if not math.isfinite(purity) or not 0.0 < purity <= 1.0:
        raise ValueError("purity must be finite and in (0, 1]")
    return (signal_cf - (1.0 - purity) * background_cf) / purity


def contrast(model, first, second):
    """Difference and ratio with shared-input covariance, no global p-value."""
    difference = first - second
    result = dict(difference=model.estimate(difference),
                  covariance=model.covariance(first, second), ratio=None,
                  approximate_difference_over_error=None)
    if second.value != 0.0:
        result['ratio'] = model.estimate(first / second)
    error = result['difference']['error']
    if error > 0.0:
        result['approximate_difference_over_error'] = difference.value / error
    return result


def integrated_cf(numerator, denominator, norm_indices, selected_indices):
    """CF of selected integrated counts, not an arithmetic mean of bin CFs.

    This is a low-dimensional contrast useful when plotted bins share a common
    normalization. None marks nonpositive ME or normalization integrals.
    """
    _same_length(numerator, denominator)
    norm_indices = list(norm_indices)
    selected_indices = list(selected_indices)
    anorm = integral(numerator, norm_indices)
    bnorm = integral(denominator, norm_indices)
    selected_a = integral(numerator, selected_indices)
    selected_b = integral(denominator, selected_indices)
    if anorm.value <= 0.0 or bnorm.value <= 0.0 or selected_b.value <= 0.0:
        return None
    return (bnorm / anorm) * (selected_a / selected_b)
