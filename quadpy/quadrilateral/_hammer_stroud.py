# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ..helpers import article, fsd, pm, untangle, z
from ._helpers import QuadrilateralScheme

citation = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    volume="12",
    year="1958",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)


def hammer_stroud_1_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    data = [(1, fsd(2, (sqrt(frac(2, 3)), 1)))]
    points, weights = untangle(data)
    return QuadrilateralScheme("Hammer-Stroud 1-2", weights, points, 3, citation)


def hammer_stroud_2_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    alpha = sqrt(frac(3, 5))
    data = [
        (frac(64, 81), z(2)),
        (frac(40, 81), fsd(2, (alpha, 1))),
        (frac(25, 81), pm(2, alpha)),
    ]
    points, weights = untangle(data)
    return QuadrilateralScheme("Hammer-Stroud 2-2", weights, points, 5, citation)


def hammer_stroud_3_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    xi1, xi2 = [sqrt(frac(3, 287) * (38 - i * sqrt(583))) for i in [+1, -1]]
    data = [
        (frac(98, 405), fsd(2, (sqrt(frac(6, 7)), 1))),
        (0.5205929166673945, pm(2, xi1)),
        (0.2374317746906302, pm(2, xi2)),
    ]
    points, weights = untangle(data)
    return QuadrilateralScheme("Hammer-Stroud 3-2", weights, points, 7, citation)
