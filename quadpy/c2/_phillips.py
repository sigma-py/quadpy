from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, concat, pm2, symm_r0

citation = article(
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
    r, s = [sqrt((105 + i * c) / 140) for i in [+1, -1]]
    t = sqrt(frac(3, 5))

    B1, B2 = [(77 - i * c) / 891 for i in [+1, -1]]
    B3 = frac(25, 324)

    weights, points = concat(symm_r0([B1, r], [B2, s]), pm2([B3, t, t]))
    return C2Scheme("Phillips", weights, points, 7, citation)
