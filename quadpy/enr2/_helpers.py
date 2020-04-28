from functools import reduce
import operator
import math

import numpy
from sympy import sqrt, pi


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
        return dot(f(flt(self.points).T), flt(self.weights))


def integrate_monomial_over_enr2(alpha, symbolic=False):
    if any(k % 2 == 1 for k in alpha):
        return 0

    if symbolic:

        def prod(factors):
            return reduce(operator.mul, factors, 1)

        k2 = [kk // 2 for kk in alpha]
        return prod([sqrt(pi) * prod(range(kk + 1, 2 * kk + 1)) / 4 ** kk for kk in k2])

    # return numpy.prod([math.gamma(0.5 * (k + 1)) for k in alpha])
    k2 = [kk // 2 for kk in alpha]
    return numpy.prod(
        [
            numpy.sqrt(numpy.pi) * numpy.prod(numpy.arange(kk + 1, 2 * kk + 1) / 4)
            for kk in k2
        ]
    )
