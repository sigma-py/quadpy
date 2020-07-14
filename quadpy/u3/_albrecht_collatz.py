import numpy
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, pm, pm_roll, untangle
from ._helpers import U3Scheme, cartesian_to_spherical_sympy

source = article(
    authors=["J. Albrecht", "L. Collatz"],
    title="Zur numerischen Auswertung mehrdimensionaler Integrale",
    journal="ZAMM",
    volume="38",
    number="1-2",
    year="1958",
    pages="1–15",
    url="https://doi.org/10.1002/zamm.19580380102",
)


def albrecht_collatz_1():
    r, s = [sqrt((5 + t * sqrt(5)) / 10) for t in [+1, -1]]
    data = [(frac(1, 12), pm_roll([r, s, 0]))]

    points, weights = untangle(data)
    points = numpy.ascontiguousarray(points.T)
    theta_phi = cartesian_to_spherical_sympy(points)
    return U3Scheme("Albrecht-Collatz 1", weights, points, theta_phi, 5, source)


def albrecht_collatz_2():
    r = 1
    s = sqrt(frac(1, 3))
    data = [(frac(8, 120), fsd(3, (r, 1))), (frac(9, 120), pm([s, s, s]))]

    points, weights = untangle(data)
    points = numpy.ascontiguousarray(points.T)
    theta_phi = cartesian_to_spherical_sympy(points)
    return U3Scheme("Albrecht-Collatz 2", weights, points, theta_phi, 5, source)


def albrecht_collatz_3():
    r = 1
    s = sqrt(frac(1, 2))
    data = [(frac(1, 30), fsd(3, (r, 1))), (frac(2, 30), fsd(3, (s, 2)))]

    points, weights = untangle(data)
    points = numpy.ascontiguousarray(points.T)
    theta_phi = cartesian_to_spherical_sympy(points)
    return U3Scheme("Albrecht-Collatz 3", weights, points, theta_phi, 5, source)


def albrecht_collatz_4():
    r, s = [sqrt((3 + t * sqrt(5)) / 6) for t in [+1, -1]]
    t = sqrt(frac(1, 3))
    data = [
        (frac(1, 20), pm_roll([r, s, 0])),
        (frac(1, 20), pm([t, t, t])),
    ]

    points, weights = untangle(data)
    points = numpy.ascontiguousarray(points.T)
    theta_phi = cartesian_to_spherical_sympy(points)
    return U3Scheme("Albrecht-Collatz 4", weights, points, theta_phi, 5, source)


def albrecht_collatz_5():
    r = 1
    s = sqrt(frac(1, 2))
    t = sqrt(frac(1, 3))

    data = [
        (frac(40, 840), fsd(3, (r, 1))),
        (frac(32, 840), fsd(3, (s, 2))),
        (frac(27, 840), pm([t, t, t])),
    ]

    points, weights = untangle(data)
    points = numpy.ascontiguousarray(points.T)
    theta_phi = cartesian_to_spherical_sympy(points)
    return U3Scheme("Albrecht-Collatz 5", weights, points, theta_phi, 7, source)
