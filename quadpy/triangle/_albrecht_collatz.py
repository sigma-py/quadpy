# -*- coding: utf-8 -*-
#

import sympy

from ..helpers import article
from ._helpers import TriangleScheme, s2

citation = article(
    authors=["J. Albrecht", "L. Collatz"],
    title="Zur numerischen Auswertung mehrdimensionaler Integrale",
    journal="ZAMM",
    volume="38",
    number="1-2",
    year="1958",
    pages="1–15",
    url="https://doi.org/10.1002/zamm.19580380102",
)


def albrecht_collatz(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = s2([frac(2, 30), frac(1, 2)], [frac(9, 15), frac(1, 6)])
    weights /= 2
    return TriangleScheme("Albrecht-Collatz", weights, points, 3, citation)
