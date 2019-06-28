# -*- coding: utf-8 -*-
#
from sympy import Rational as frac

from ..helpers import article, untangle
from ._helpers import HexahedronScheme, fs_r00, fs_rr0, pm_rrr, z

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


def mustard_lyness_blatt_1():
    data = [(frac(1, 2), z()), (frac(1, 24), fs_rr0(1))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 1", weights, points, 3, _citation)


def mustard_lyness_blatt_2():
    data = [(frac(2, 9), z()), (frac(1, 9), fs_r00(1)), (frac(1, 72), pm_rrr(1))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 2", weights, points, 3, _citation)


def mustard_lyness_blatt_3():
    data = [(+frac(1, 6), fs_rr0(1)), (-frac(1, 8), pm_rrr(1))]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 3", weights, points, 3, _citation)


def mustard_lyness_blatt_4():
    data = [
        (-frac(2, 45), z()),
        (+frac(2, 45), fs_r00(1)),
        (+frac(4, 45), pm_rrr(frac(1, 2))),
        (frac(1, 120), pm_rrr(1)),
    ]

    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 4", weights, points, 5, _citation)


def mustard_lyness_blatt_5():
    data = [
        (-frac(19, 15), z()),
        (+frac(16, 45), fs_r00(frac(1, 2))),
        (-frac(1, 30), fs_r00(1)),
        (+frac(1, 36), fs_rr0(1)),
    ]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 5", weights, points, 5, _citation)


def mustard_lyness_blatt_6():
    data = [
        (-frac(4, 3), z()),
        (+frac(16, 45), fs_r00(frac(1, 2))),
        (frac(1, 90), fs_rr0(1)),
        (frac(1, 120), pm_rrr(1)),
    ]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 6", weights, points, 5, _citation)


def mustard_lyness_blatt_7():
    data = [
        (frac(2, 45), z()),
        (frac(1, 45), fs_rr0(1)),
        (frac(4, 45), pm_rrr(frac(1, 2))),
        (frac(-1, 360), pm_rrr(1)),
    ]
    points, weights = untangle(data)
    weights *= 8
    return HexahedronScheme("Mustard-Lyness-Blatt 7", weights, points, 5, _citation)
