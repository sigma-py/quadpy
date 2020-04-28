import numpy
from sympy import sqrt, Rational as frac, pi

from ..helpers import article
from ._helpers import NSphereScheme

citation = article(
    authors=["I.P. Mysovskikh"],
    title="The approximation of multiple integrals by using interpolatory cubature formulae",
    journal="Proceedings of a Symposium on Quantitative Approximation Held in Bonn, West Germany, August 20–24, 1979",
    year="1980",
    pages="217-243",
    url="https://doi.org/10.1016/B978-0-12-213650-4.50025-8",
)


# Compute the area of a hypersphere, <https://tauday.com/tau-manifesto>.
# No roots, no gamma functions. Such a nice recursion!
def surface_hypersphere(n):
    assert n > 1
    if n == 2:
        return 2 * pi
    elif n == 3:
        return 4 * pi
    return 2 * pi / (n - 2) * surface_hypersphere(n - 2)


def get_nsimplex_points(n):
    # vertices of the n-simplex
    points = []
    for r in range(n + 1):
        point = [-sqrt(frac(n + 1, n * (n - i + 1) * (n - i))) for i in range(r)]
        if r < n:
            point += [sqrt(frac((n + 1) * (n - r), n * (n - r + 1)))]
            point += [0] * (n - r - 1)
        points.append(point)
    return numpy.array(points)


def mysovskikh_1(n):
    points = get_nsimplex_points(n)
    weights = numpy.full(n + 1, surface_hypersphere(n) / (n + 1))
    return NSphereScheme("Mysovskikh 1", n, weights, points, 2, citation)


def mysovskikh_2(n):
    A = frac((7 - n) * n, 2 * (n + 1) ** 2 * (n + 2)) * surface_hypersphere(n)
    B = frac(2 * (n - 1) ** 2, n * (n + 1) ** 2 * (n + 2)) * surface_hypersphere(n)

    a = get_nsimplex_points(n)
    b = numpy.array([
        sqrt(frac(n, 2 * (n - 1))) * (a[k] + a[l])
        for k in range(n + 1)
        for l in range(k)
    ])

    points = numpy.concatenate([a, -a, b, -b])
    weights = numpy.concatenate([
        numpy.full(n + 1, A),
        numpy.full(n + 1, A),
        numpy.full((n * (n + 1)) // 2, B),
        numpy.full((n * (n + 1)) // 2, B),
    ])

    return NSphereScheme("Mysovskikh 2", n, weights, points, 5, citation)
