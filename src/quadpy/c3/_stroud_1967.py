from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C3Scheme, register

source = article(
    authors=["A.H. Stroud"],
    title="Some fifth degree integration formulas for symmetric regions II",
    journal="Numerische Mathematik",
    volume="9",
    number="5",
    month="apr",
    year="1967",
    pages="460-468",
    url="https://doi.org/10.1007/BF02162160",
)


def stroud_1967():
    # ENH Analytic expression for all quantities are given in
    #
    # J.W. Peterson,
    # Analytical Formulae for Two of A. H. Stroud's Quadrature Rules,
    # Sep. 2009,
    # <https://arxiv.org/pdf/0909.5106.pdf>.
    sqrt19 = sqrt(19)
    t = sqrt(71440 + 6802 * sqrt19)

    lmbd, gmma = (sqrt((1919 - 148 * sqrt19 + i * 4 * t) / 3285) for i in [+1, -1])
    xi, mu = (-sqrt((1121 + 74 * sqrt19 - i * 2 * t) / 3285) for i in [+1, -1])
    mu *= -1

    B, C = (
        133225 / (260072 - 1520 * sqrt19 + i * (133 - 37 * sqrt19) * t)
        for i in [+1, -1]
    )

    d = {
        "zero3": [[frac(32, 19) / 8]],
        "symm_rss_pm": [[B / 8, C / 8], [lmbd, gmma], [xi, mu]],
    }
    return C3Scheme("Stroud 1967", d, 5, source)


register([stroud_1967])
