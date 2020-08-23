from sympy import Rational as frac

from ..helpers import article
from ._helpers import C3Scheme, register

_source = article(
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
    d = {
        "zero3": [[frac(1, 4)]],
        "symm_r00": [[frac(1, 12)], [1]],
        "symm_rr0": [[frac(1, 48)], [1]],
    }
    return C3Scheme("Albrecht-Collatz", d, 3, _source)


register([albrecht_collatz])
