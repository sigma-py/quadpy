from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import E2rScheme, register

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

    d = {
        "zero2": [[frac(7, 10)]],
        "c2_a0": [[frac(1, 20)], [nu]],
        "sxy": [[frac(1, 20)], [xi], [eta]],
    }
    return E2rScheme("Stroud-Secrest V", d, 5, _source)


def stroud_secrest_6():
    sqrt74255 = sqrt(74255)

    nu = sqrt(42)
    xi, eta = (sqrt((6615 - p_m * 21 * sqrt74255) / 454) for p_m in [+1, -1])
    A = frac(5, 588)
    B, C = ((5272105 + p_m * 18733 * sqrt74255) / 43661940 for p_m in [+1, -1])

    d = {
        "d4_a0": [[A], [nu]],
        "d4_aa": [[B, C], [xi, eta]],
    }
    return E2rScheme("Stroud-Secrest VI", d, 7, _source)


register([stroud_secrest_5, stroud_secrest_6])
