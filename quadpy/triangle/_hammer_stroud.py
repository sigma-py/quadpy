# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from ._helpers import s2, TriangleScheme, s3, concat
from ..helpers import article

citation = article(
    authors=["Preston C. Hammer", "Arthur H. Stroud"],
    title="Numerical Evaluation of Multiple Integrals II",
    journal="Math. Comp.",
    volume="12",
    year="1958",
    pages="272-280",
    url="https://doi.org/10.1090/S0025-5718-1958-0102176-6",
)


def hammer_stroud_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = s2([frac(1, 3), frac(1, 6)])
    return TriangleScheme("Hammer-Stroud 2", weights, points, 2, citation)


def hammer_stroud_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = concat(s3(-frac(27, 48)), s2([+frac(25, 48), frac(1, 5)]))
    return TriangleScheme("Hammer-Stroud 3", weights, points, 3, citation)
