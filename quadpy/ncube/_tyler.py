# -*- coding: utf-8 -*-
#

import sympy

from ..helpers import untangle, fsd, z, article
from ._helpers import NCubeScheme

_citation = article(
    authors=["G.W. Tyler"],
    title="Numerical integration of functions of several variables",
    journal="Canad. J. Math.",
    volume="5",
    year="1953",
    pages="393-412",
    url="https://doi.org/10.4153/CJM-1953-044-1",
)


def tyler(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [(frac(3 - n, 3), z(n)), (frac(1, 6), fsd(n, (1, 1)))]

    points, weights = untangle(data)
    reference_volume = 2 ** n
    weights *= reference_volume
    return NCubeScheme("Tyler", n, weights, points, 3, _citation)
