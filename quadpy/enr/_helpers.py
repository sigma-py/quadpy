import math

import numpy
import sympy

from ..helpers import QuadratureScheme


class EnrScheme(QuadratureScheme):
    def __init__(self, name, dim, weights, points, degree, source, tol=1.0e-14):
        super().__init__(name, weights, points, degree, source, tol)
        self.domain = f"Enr (n={dim})"
        self.dim = dim

    def integrate(self, f, dot=numpy.dot):
        flt = numpy.vectorize(float)
        ref_vol = enr_volume(self.dim)
        return ref_vol * dot(f(flt(self.points).T), flt(self.weights))


# The closed formula is
#
#   2
#   * math.factorial(sum(alpha) + n - 1)
#   * prod([math.gamma((k + 1) / 2.0) for k in alpha])
#   / math.gamma((sum(alpha) + n) / 2).
#
# Care must be taken when evaluating this expression as numerator or denominator will
# quickly overflow. A better representation is via a recurrence. This is numerically
# stable and can easily be used for symbolic computation.
def integrate_monomial_over_enr(k, symbolic=False):
    n = len(k)
    if any(a % 2 == 1 for a in k):
        return 0
    if all(a == 0 for a in k):
        return enr_volume(n, symbolic)

    # find first nonzero
    idx = next(i for i, j in enumerate(k) if j > 0)
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
