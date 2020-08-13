from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import E2r2Scheme, expand_symmetries

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
    nu = sqrt(2)
    xi = nu / 2
    eta = sqrt(6) / 2
    A = frac(1, 2)
    B = frac(1, 12)

    d = {
        "zero": [[A]],
        "pmx": [[B], [nu]],
        "pm2": [[B], [xi], [eta]]
    }
    points, weights = expand_symmetries(d)
    return E2r2Scheme("Stroud-Secrest V", weights, points, 5, _source)


def stroud_secrest_6():
    sqrt5 = sqrt(5)
    nu = sqrt(3)
    xi, eta = [sqrt((9 - p_m * 3 * sqrt5) / 8) for p_m in [+1, -1]]
    A = frac(1, 36)
    B, C = [(5 + p_m * 2 * sqrt5) / 45 for p_m in [+1, -1]]

    d = {
        "s40": [[A], [nu]],
        "s4": [[B, C], [xi, eta]]
    }
    points, weights = expand_symmetries(d)
    return E2r2Scheme("Stroud-Secrest VI", weights, points, 7, _source)
