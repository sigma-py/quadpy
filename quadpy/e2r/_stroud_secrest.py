import numpy
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, pm0, untangle
from ._helpers import E2rScheme

_source = article(
    authors=["A.H. Stroud", "D. Secrest"],
    title="Approximate integration formulas for certain spherically symmetric regions",
    journal="Math. Comp.",
    volume="17",
    year="1963",
    pages="105-135",
    url="https://doi.org/10.1090/S0025-5718-1963-0161473-0",
)


def stroud_secrest_5():
    nu = 2 * sqrt(5)
    xi = sqrt(5)
    eta = sqrt(15)

    data = [
        (frac(7, 10), numpy.array([[0, 0]])),
        (frac(1, 20), numpy.array([[+nu, 0], [-nu, 0]])),
        (frac(1, 20), pm0([xi, eta])),
    ]
    points, weights = untangle(data)
    return E2rScheme("Stroud-Secrest V", weights, points, 5, _source)


def stroud_secrest_6():
    sqrt74255 = sqrt(74255)

    nu = sqrt(42)
    xi, eta = [sqrt((6615 - p_m * 21 * sqrt74255) / 454) for p_m in [+1, -1]]
    A = frac(5, 588)
    B, C = [(5272105 + p_m * 18733 * sqrt74255) / 43661940 for p_m in [+1, -1]]

    data = [(A, fsd(2, (nu, 1))), (B, pm0([xi, xi])), (C, pm0([eta, eta]))]

    points, weights = untangle(data)
    return E2rScheme("Stroud-Secrest VI", weights, points, 7, _source)
