from sympy import Rational as frac

from ..helpers import article, untangle
from ._helpers import HexahedronScheme, fs_r00, fs_rr0, z

_citation = article(
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
    data = [(frac(1, 4), z()), (frac(1, 12), fs_r00(1)), (frac(1, 48), fs_rr0(1))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Albrecht-Collatz", weights, points, 3, _citation)
