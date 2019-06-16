# -*- coding: utf-8 -*-
#
"""
D. Mustard, J.N. Lyness, J.M. Blatt,
Numerical quadrature in n dimensions,
Comput J (1963) 6 (1): 75-87,
<https://doi.org/10.1093/comjnl/6.1.75>.
"""
from __future__ import division

import sympy

from .helpers import fs_rr0, fs_r00, pm_rrr, z, HexahedronScheme
from ..helpers import untangle


def mustard_lyness_blatt_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    data = [(frac(1, 2), z()), (frac(1, 24), fs_rr0(1))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 1", 3, weights, points)


def mustard_lyness_blatt_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [(frac(2, 9), z()), (frac(1, 9), fs_r00(1)), (frac(1, 72), pm_rrr(1))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 2", 3, weights, points)


def mustard_lyness_blatt_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [(+frac(1, 6), fs_rr0(1)), (-frac(1, 8), pm_rrr(1))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 3", 3, weights, points)


def mustard_lyness_blatt_4(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [
        (-frac(2, 45), z()),
        (+frac(2, 45), fs_r00(1)),
        (+frac(4, 45), pm_rrr(frac(1, 2))),
        (frac(1, 120), pm_rrr(1)),
    ]

    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 4", 5, weights, points)


def mustard_lyness_blatt_5(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [
        (-frac(19, 15), z()),
        (+frac(16, 45), fs_r00(frac(1, 2))),
        (-frac(1, 30), fs_r00(1)),
        (+frac(1, 36), fs_rr0(1)),
    ]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 5", 5, weights, points)


def mustard_lyness_blatt_6(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [
        (-frac(4, 3), z()),
        (+frac(16, 45), fs_r00(frac(1, 2))),
        (frac(1, 90), fs_rr0(1)),
        (frac(1, 120), pm_rrr(1)),
    ]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 6", 5, weights, points)


def mustard_lyness_blatt_7(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    data = [
        (frac(2, 45), z()),
        (frac(1, 45), fs_rr0(1)),
        (frac(4, 45), pm_rrr(frac(1, 2))),
        (frac(-1, 360), pm_rrr(1)),
    ]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 7", 5, weights, points)


MustardLynessBlatt = {
    1: mustard_lyness_blatt_1,
    2: mustard_lyness_blatt_2,
    3: mustard_lyness_blatt_3,
    4: mustard_lyness_blatt_4,
    5: mustard_lyness_blatt_5,
    6: mustard_lyness_blatt_6,
    7: mustard_lyness_blatt_7,
}
