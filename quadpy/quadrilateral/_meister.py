from sympy import Rational as frac

from ..helpers import article
from ._helpers import QuadrilateralScheme, concat, symm_r0, symm_s, symm_s_t, zero

citation = article(
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

    weights, points = concat(
        zero(frac(1024, 6720)),
        symm_s([frac(576, 6720), r], [-frac(9, 6720), s], [frac(47, 6720), 1]),
        symm_r0([frac(576, 6720), r]),
        symm_s_t([frac(117, 6720), 1, s]),
    )

    weights *= 4
    return QuadrilateralScheme("Meister", weights, points, 7, citation)
