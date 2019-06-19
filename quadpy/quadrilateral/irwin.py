# -*- coding: utf-8 -*-
#
from __future__ import division

import sympy

from .helpers import concat, symm_s, symm_s_t, QuadrilateralScheme
from ..helpers import book

citation = book(
    authors=["Joseph Oscar Irwin"],
    title="On quadrature and cubature",
    publisher="Cambridge University Press",
    year="1923",
    url="https://books.google.de/books/about/On_quadrature_and_cubature.html?id=SuruAAAAMAAJ",
)


def irwin_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = concat(symm_s([frac(14, 48), 1]), symm_s_t([-frac(1, 48), 3, 1]))
    weights *= 4
    return QuadrilateralScheme("Irwin 1", weights, points, 3, citation)


def irwin_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    weights, points = concat(
        symm_s([frac(889, 2880), 1], [frac(5, 2880), 3]),
        symm_s_t([-frac(98, 2880), 3, 1], [frac(11, 2880), 5, 1]),
    )
    weights *= 4
    return QuadrilateralScheme("Irwin 2", weights, points, 5, citation)
