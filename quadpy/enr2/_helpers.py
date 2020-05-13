import math

import numpy
import sympy


class Enr2Scheme:
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
        ref_vol = volume_enr2(self.dim)
        return ref_vol * dot(f(flt(self.points).T), flt(self.weights))


# One could simply write
#
#     sqrt(pi) ** n
#
# but all other n-dimensional volumes are recurrences as well, so keep the tradition
# alive. :)
def volume_enr2(n, symbolic=False):
    sqrt = sympy.sqrt if symbolic else math.sqrt
    pi = sympy.pi if symbolic else math.pi
    if n == 0:
        return 1
    elif n == 1:
        return sqrt(pi)
    return volume_enr2(n - 2, symbolic) * pi


#  numpy.prod([math.gamma(0.5 * (k + 1)) for k in alpha])
def integrate_monomial_over_enr2(k, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    if any(a % 2 == 1 for a in k):
        return 0

    n = len(k)
    if all(a == 0 for a in k):
        return volume_enr2(n, symbolic)

    # find first nonzero
    idx = next(i for i, j in enumerate(k) if j > 0)
    alpha = frac(k[idx] - 1, 2)
    k2 = k.copy()
    k2[idx] -= 2
    return integrate_monomial_over_enr2(k2, symbolic) * alpha
