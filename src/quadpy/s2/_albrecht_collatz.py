import sympy

from ..helpers import article
from ._helpers import S2Scheme, register

_source = article(
    authors=["J. Albrecht", "L. Collatz"],
    title="Zur numerischen Auswertung mehrdimensionaler Integrale",
    journal="ZAMM",
    volume="38",
    number="1-2",
    year="1958",
    month="jan",
    pages="1–15",
    url="https://doi.org/10.1002/zamm.19580380102",
)

frac = sympy.Rational


def albrecht_collatz():
    # ERR Incorrectly stated in Stroud as sqrt(1/2) instead of 1/2
    d = {"d4_aa": [[frac(1, 4)], [frac(1, 2)]]}
    return S2Scheme("Albrecht-Collatz", d, 3, _source)


register([albrecht_collatz])
