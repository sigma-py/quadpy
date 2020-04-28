import math
import operator
from functools import reduce

import numpy
from sympy import pi, sqrt


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


def integrate_monomial_over_enr(alpha, symbolic=False):
    if any(k % 2 == 1 for k in alpha):
        return 0

    if symbolic:

        def prod(factors):
            return reduce(operator.mul, factors, 1)

        def fact(k):
            return prod(range(1, k + 1))

        n = len(alpha)

        alpha2 = [k // 2 for k in alpha]

        # TODO find a nicer expression for this whole thing
        # Check out <https://tauday.com/tau-manifesto>
        if n % 2 == 0:
            b = fact(sum(alpha2) + n // 2 - 1)
        else:
            # b = gamma(sum(alpha2) + n // 2 + 0.5)
            k = sum(alpha2) + n // 2
            # b = sqrt(pi) * fact(2 * sum(alpha2) + n) / fact(k) / 4 ** k
            b = sqrt(pi) * fact(2 * k) / fact(k) / 4 ** k

        out = (
            2
            * fact(sum(alpha) + n - 1)
            * numpy.prod(
                [sqrt(pi) * prod(range(k + 1, 2 * k + 1)) / 4 ** k for k in alpha2]
            )
            / b
        )
        return out

    n = len(alpha)
    return (
        2
        * math.factorial(sum(alpha) + n - 1)
        * numpy.prod([math.gamma((k + 1) / 2.0) for k in alpha])
        / math.gamma((sum(alpha) + n) / 2)
    )
