# -*- coding: utf-8 -*-
#
"""
Preston C. Hammer and Arthur H. Stroud,
Numerical Evaluation of Multiple Integrals II,
Math. Comp. 12 (1958), 272-280,
<https://doi.org/10.1090/S0025-5718-1958-0102176-6>.
"""
from __future__ import division

import sympy

from .helpers import s2, TriangleScheme, s3, concat


def hammer_stroud_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = s2([frac(1, 3), frac(1, 6)])
    return TriangleScheme("Hammer-Stroud 2", 2, weights, bary)


def hammer_stroud_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, bary = concat(s3(-frac(27, 48)), s2([+frac(25, 48), frac(1, 5)]))
    return TriangleScheme("Hammer-Stroud 3", 3, weights, bary)


HammerStroud = {2: hammer_stroud_2, 3: hammer_stroud_3}
