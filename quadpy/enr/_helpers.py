import math

import numpy
import sympy


class EnrScheme:
    def __init__(self, name, dim, weights, points, degree, citation):
        self.name = name
        self.dim = dim
        self.degree = degree
        self.citation = citation

        if weights.dtype == numpy.float64:
            self.weights = weights
        else:
            assert weights.dtype in [numpy.dtype("O"), numpy.int_]
            self.weights = weights.astype(numpy.float64)
            self.weights_symbolic = weights

        if points.dtype == numpy.float64:
            self.points = points
        else:
            assert points.dtype in [numpy.dtype("O"), numpy.int_]
            self.points = points.astype(numpy.float64)
            self.points_symbolic = points
        return

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        return dot(f(flt(self.points).T), flt(self.weights))


# The closed formula is
#
#   2
#   * math.factorial(sum(alpha) + n - 1)
#   * prod([math.gamma((k + 1) / 2.0) for k in alpha])
#   / math.gamma((sum(alpha) + n) / 2).
#
# Care must be taken when evaluating this expression as numerator or denominator will
# quickly overflow. A better representation is via a recurrence. This numerically stable
# and can easily be used for symbolic computation.
def integrate_monomial_over_enr(k, symbolic=False):
    n = len(k)
    if any(a % 2 == 1 for a in k):
        return 0
    if all(a == 0 for a in k):
        return enr_volume(n, symbolic)

    # find first nonzero
    idx = next((i for i, j in enumerate(k) if j > 0), None)
    alpha = (k[idx] - 1) * (sum(k) + n - 1)
    k2 = k.copy()
    k2[idx] -= 2
    return integrate_monomial_over_enr(k2, symbolic) * alpha


# 2 * sqrt(pi) ** n * gamma(n) / gamma(frac(n, 2))
def enr_volume(n, symbolic=False):
    pi = sympy.pi if symbolic else math.pi
    if n == 1:
        return 2
    elif n == 2:
        return 2 * pi
    return 2 * pi * (n - 1) * enr_volume(n - 2)
    # Then n-sphere has
    # return 2 * pi / (n - 2) * surface_hypersphere(n - 2)
    # here.
