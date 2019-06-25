# -*- coding: utf-8 -*-
#

import sympy

from ..helpers import article, pm, untangle, z
from ._helpers import NCubeScheme

_citation = article(
    authors=["G.M. Ewing"],
    title="On Approximate Cubature",
    journal="The American Mathematical Monthly",
    volume="48",
    number="2",
    month="feb",
    year="1941",
    pages="134-136",
    url="https://doi.org/10.2307/2303604",
)


def ewing(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [(frac(2, 3), z(n)), (frac(1, 3 * 2 ** n), pm(n, 1))]

    points, weights = untangle(data)
    reference_volume = 2 ** n
    weights *= reference_volume
    return NCubeScheme("Ewing", n, weights, points, 3, _citation)
