# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ..helpers import untangle, fsd, z, pm, article
from ._helpers import NCubeScheme

_citation = article(
    authors=["D. Mustard", "J.N. Lyness", "J.M. Blatt"],
    title="Numerical quadrature in n dimensions",
    journal="Comput J",
    year="1963",
    volume="6",
    number="1",
    pages="75-87",
    url="https://doi.org/10.1093/comjnl/6.1.75",
)


def mustard_lyness_blatt(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = sqrt(frac(2, 5))
    data = [
        (frac(8 - 5 * n, 9), z(n)),
        (frac(5, 18), fsd(n, (r, 1))),
        (frac(1, 9 * 2 ** n), pm(n, 1)),
    ]

    points, weights = untangle(data)
    reference_volume = 2 ** n
    weights *= reference_volume
    return NCubeScheme("Mustard-Lyness-Blatt", n, weights, points, 5, _citation)
