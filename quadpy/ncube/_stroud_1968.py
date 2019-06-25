# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ._helpers import _s2, _s11, NCubeScheme
from ..helpers import untangle, fsd, z, article

_citation = article(
    authors=["A.H. Stroud"],
    title="Extensions of Symmetric Integration Formulas",
    journal="Mathematics of Computation",
    volume="22",
    number="102",
    month="apr",
    year="1968",
    pages="271-274",
    url="https://doi.org/10.2307/2004655",
)


def stroud_1968(n, symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    r = sqrt(frac(7, 15))
    s, t = [sqrt((7 + i * sqrt(24)) / 15) for i in [+1, -1]]
    data = [
        (frac(5 * n ** 2 - 15 * n + 14, 14), z(n)),
        (frac(25, 168), _s2(n, +r)),
        (frac(25, 168), _s2(n, -r)),
        (frac(-25 * (n - 2), 168), fsd(n, (r, 1))),
        (frac(5, 48), _s11(n, +s, -t)),
        (frac(5, 48), _s11(n, -s, +t)),
        (frac(-5 * (n - 2), 48), fsd(n, (s, 1))),
        (frac(-5 * (n - 2), 48), fsd(n, (t, 1))),
    ]

    points, weights = untangle(data)
    weights *= 2 ** n
    return NCubeScheme("Stroud 1968", n, weights, points, 5, _citation)
