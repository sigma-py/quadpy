import math

import numpy
import sympy


class UnScheme:
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

    def integrate(self, f, center, radius, dot=numpy.dot):
        center = numpy.array(center)
        rr = numpy.multiply.outer(radius, self.points)
        rr = numpy.swapaxes(rr, 0, -2)
        ff = numpy.array(f((rr + center).T))
        return numpy.array(radius) ** (self.dim - 1) * dot(ff, self.weights)


# The article
#
#     Gerald B. Folland,
#     How to Integrate a Polynomial over a Sphere,
#     The American Mathematical Monthly,
#     Vol. 108, No. 5 (May, 2001), pp. 446-448,
#     <https://doi.org/10.2307/2695802>
#
# gives the formula
#
#     2 * (
#         prod([gamma(Rational(a + 1, 2)) for a in alpha])
#         / gamma(sum([Rational(a + 1, 2) for a in alpha]))
#     )
#
# which is unsuitable for numerical calculations because of quick overflows in numerator
# and denominator. This can be saved by the of exp-lgamma, but a more reasonable
# approach is to use recurrence.
def integrate_monomial_over_unit_nsphere(k, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    if any(a % 2 == 1 for a in k):
        return 0

    n = len(k)
    if all(a == 0 for a in k):
        return sphere_volume(n - 1, symbolic, r=1)

    # find first nonzero
    idx = next(i for i, j in enumerate(k) if j > 0)
    alpha = frac(k[idx] - 1, sum(k) + n - 2)
    k2 = k.copy()
    k2[idx] -= 2
    return integrate_monomial_over_unit_nsphere(k2, symbolic) * alpha


# n sqrt(pi) ** 2 / gamma(n/2 + 1)
def sphere_volume(n, symbolic=False, r=1):
    pi = sympy.pi if symbolic else math.pi
    if n == 0:
        return 2
    elif n == 1:
        return 2 * pi * r
    return (2 * pi) / (n - 1) * r ** 2 * sphere_volume(n - 2, symbolic, r=r)
