from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, register

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
        "d4_a0": [[frac(10, 49)], [sqrt(frac(7, 15))]],
        "d4_aa": [[frac(9, 196)], [sqrt(frac(7, 9))]],
    }
    return C2Scheme("Burnside", d, 5, source)


register([burnside])
