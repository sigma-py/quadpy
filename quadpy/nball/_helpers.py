import numpy

from ..nsphere._helpers import integrate_monomial_over_unit_nsphere


class NBallScheme:
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
        return numpy.array(radius) ** self.dim * dot(ff, self.weights)


def volume_unit_ball(n, symbolic=False):
    return integrate_monomial_over_unit_nball(n * [0], symbolic=symbolic)


def integrate_monomial_over_unit_nball(exponents, symbolic=False):
    """
    Gerald B. Folland,
    How to Integrate a Polynomial over a Sphere,
    The American Mathematical Monthly,
    Vol. 108, No. 5 (May, 2001), pp. 446-448,
    <https://doi.org/10.2307/2695802>.
    """
    radius = 1
    n = len(exponents)
    alpha = n + sum(exponents)
    return (
        radius ** alpha
        * integrate_monomial_over_unit_nsphere(exponents, symbolic=symbolic)
        / alpha
    )
