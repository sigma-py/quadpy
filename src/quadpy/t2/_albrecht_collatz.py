from sympy import Rational as frac

from ..helpers import article
from ._helpers import T2Scheme, register

source = article(
    authors=["J. Albrecht", "L. Collatz"],
    title="Zur numerischen Auswertung mehrdimensionaler Integrale",
    journal="ZAMM",
    volume="38",
    number="1-2",
    year="1958",
    pages="1–15",
    url="https://doi.org/10.1002/zamm.19580380102",
)


def albrecht_collatz():
    d = {"d3_aa": [[frac(1, 30), frac(9, 30)], [frac(1, 2), frac(1, 6)]]}
    return T2Scheme("Albrecht-Collatz", d, 3, source, tol=2.776e-16)


register([albrecht_collatz])
