from sympy import Rational as frac

from ..helpers import article
from ._helpers import C2Scheme, expand_symmetries

source = article(
    authors=["Bernd Meister"],
    title="On a Family of Cubature Formulae",
    journal="Comput J",
    year="1966",
    volume="8",
    number="4",
    pages="368-371",
    url="https://doi.org/10.1093/comjnl/8.4.368",
)


def meister():
    r = frac(2, 3)
    s = frac(1, 3)

    d = {
        "zero": [[frac(1024, 6720)]],
        "symm_s": [[frac(576, 6720), -frac(9, 6720), frac(47, 6720)], [r, s, 1]],
        "symm_r0": [[frac(576, 6720)], [r]],
        "symm_s_t": [[frac(117, 6720)], [1], [s]]
    }
    points, weights = expand_symmetries(d)
    return C2Scheme("Meister", weights, points, 7, source)
