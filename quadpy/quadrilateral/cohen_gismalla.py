# -*- coding: utf-8 -*-
#
"""
A. M. Cohen, D. A. Gismalla,
Some integration formulae for symmetric functions of two variables,
International Journal of Computer Mathematics, 1986, 19:1, 57-68,
<https://doi.org/10.1080/00207168608803504>.
"""
from __future__ import division

import warnings

import sympy

from .helpers import concat, zero, pm, QuadrilateralScheme


def cohen_gismalla_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    # TODO improve precision
    warnings.warn("The Cohen-Gismalla schemes are only given in single-precision.")

    u = 0.84623312
    v = 0.46607171
    weights, points = concat(
        zero(frac(8, 7)), pm([frac(5, 7), u, -v], [frac(5, 7), v, u])
    )
    # This scheme is of order 5 for symmetric integrands
    return QuadrilateralScheme("CohenGismalla 1", 3, weights, points)


def cohen_gismalla_2():
    # TODO improve precision
    warnings.warn("The Cohen-Gismalla schemes are only given in single-precision.")

    r = 0.5878606
    s = 0.9353943
    u = 0.6105540
    v = 0.1109710
    A = 0.1856914
    B = 0.5951448
    C = 0.3584324
    weights, points = concat(zero(A), pm([B, u, -v], [B, v, u], [C, r, -s], [C, r, s]))
    # ERR this scheme only has order 1
    # According to the article, it has order 7 for symmetric integrands.
    # Something is fishy...
    return QuadrilateralScheme("CohenGismalla 2", 1, weights, points)


CohenGismalla = {1: cohen_gismalla_1, 2: cohen_gismalla_2}
