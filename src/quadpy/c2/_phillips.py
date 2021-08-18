from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, register

source = article(
    authors=["G.M. Phillips"],
    title="Numerical integration in two and three dimensions",
    journal="Comput J",
    year="1967",
    volume="10",
    number="2",
    pages="202-204",
    url="https://doi.org/10.1093/comjnl/10.2.202",
)
# TODO add scheme for hex


def phillips():
    c = 3 * sqrt(385)
    r, s = (sqrt((105 + i * c) / 140) for i in [+1, -1])
    t = sqrt(frac(3, 5))

    B1, B2 = ((77 - i * c) / 891 for i in [+1, -1])
    B3 = frac(25, 324)

    d = {
        "d4_a0": [[B1, B2], [r, s]],
        "d4_aa": [[B3], [t]],
    }
    return C2Scheme("Phillips", d, 7, source)


register([phillips])
