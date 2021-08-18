from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import E2r2Scheme, register

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

    d = {"zero2": [[A]], "c2_a0": [[B], [nu]], "sxy": [[B], [xi], [eta]]}
    return E2r2Scheme("Stroud-Secrest V", d, 5, _source, 4.360e-16)


def stroud_secrest_6():
    sqrt5 = sqrt(5)
    nu = sqrt(3)
    xi, eta = (sqrt((9 - p_m * 3 * sqrt5) / 8) for p_m in [+1, -1])
    A = frac(1, 36)
    B, C = ((5 + p_m * 2 * sqrt5) / 45 for p_m in [+1, -1])

    d = {"d4_a0": [[A], [nu]], "d4_aa": [[B, C], [xi, eta]]}
    return E2r2Scheme("Stroud-Secrest VI", d, 7, _source, 6.661e-16)


register([stroud_secrest_5, stroud_secrest_6])
