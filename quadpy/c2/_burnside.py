from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, expand_symmetries

source = article(
    authors=["W. Burnside"],
    title="An approximate quadrature formula",
    journal="Messenger of Math.",
    volume="37",
    year="1908",
    pages="166-167",
)


def burnside():
    d = {
        "symm_r0": [[frac(10, 49)], [sqrt(frac(7, 15))]],
        "symm_s": [[frac(9, 196)], [sqrt(frac(7, 9))]],
    }
    points, weights = expand_symmetries(d)
    return C2Scheme("Burnside", weights, points, 5, source)
