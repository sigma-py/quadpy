# -*- coding: utf-8 -*-
#
"""
Joseph Oscar Irwin,
On quadrature and cubature,
Cambridge University Press, 1923,
<https://books.google.de/books/about/On_quadrature_and_cubature.html?id=SuruAAAAMAAJ&redir_esc=y>
"""
from __future__ import division

import sympy

from .helpers import concat, symm_s, symm_s_t, QuadrilateralScheme


def irwin_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = concat(symm_s([frac(14, 48), 1]), symm_s_t([-frac(1, 48), 3, 1]))
    weights *= 4
    return QuadrilateralScheme("Irwin 1", 3, weights, points)


def irwin_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = concat(
        symm_s([frac(889, 2880), 1], [frac(5, 2880), 3]),
        symm_s_t([-frac(98, 2880), 3, 1], [frac(11, 2880), 5, 1]),
    )
    weights *= 4
    return QuadrilateralScheme("Irwin 2", 5, weights, points)


Irwin = {1: irwin_1, 2: irwin_2}
