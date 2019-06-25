# -*- coding: utf-8 -*-
#
from __future__ import division

import numpy
import sympy

from ._helpers import _s4, _s4_0, PyramidScheme
from ..helpers import untangle, article

citation = article(
    authors=["Carlos Felippa"],
    title="A compendium of FEM integration formulas for symbolic work",
    journal="Engineering Computation",
    volume="21",
    number="8",
    year="2004",
    pages="867-890",
    url="https://doi.org/10.1108/02644400410554362",
)


def felippa_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 1
    data = [(frac(128, 27), numpy.array([[0, 0, -frac(1, 2)]]))]
    points, weights = untangle(data)
    return PyramidScheme("Felippa 1", weights, points, degree, citation)


def felippa_2(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 2
    data = [
        (frac(81, 100), _s4(8 * sqrt(frac(2, 15)) / 5, -frac(2, 3))),
        (frac(125, 27), numpy.array([[0, 0, frac(2, 5)]])),
    ]
    points, weights = untangle(data)
    return PyramidScheme("Felippa 2", weights, points, degree, citation)


def felippa_3(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 2
    data = [
        (frac(504, 625), _s4(sqrt(frac(12, 35)), -frac(2, 3))),
        (frac(576, 625), numpy.array([[0, 0, frac(1, 6)]])),
        (frac(64, 15), numpy.array([[0, 0, frac(1, 2)]])),
    ]
    points, weights = untangle(data)
    return PyramidScheme("Felippa 3", weights, points, degree, citation)


def felippa_4(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 3
    w1 = 5 * (68 + 5 * sqrt(10)) / 432
    w2 = frac(85, 54) - w1
    g1 = sqrt(frac(1, 3))
    g2 = (2 * sqrt(10) - 5) / 15
    data = [(w1, _s4(g1, g2)), (w2, _s4(g1, -frac(2, 3) - g2))]
    points, weights = untangle(data)
    return PyramidScheme("Felippa 4", weights, points, degree, citation)


def felippa_5(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 2
    w1 = (11764 - 461 * sqrt(51)) / 15300
    w2 = frac(346, 225) - w1
    g1, g2 = [sqrt(frac(2, 15) * (573 - i * 2 * sqrt(51))) / 15 for i in [+1, -1]]
    g3, g4 = [-i * (2 * sqrt(51) + i * 13) / 35 for i in [+1, -1]]
    data = [(w1, _s4(g1, g3)), (w2, _s4(g2, g4))]
    points, weights = untangle(data)
    return PyramidScheme("Felippa 5", weights, points, degree, citation)


def felippa_6(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 2
    w1 = 7 * (11472415 - 70057 * sqrt(2865)) / 130739500
    w2 = frac(84091, 68450) - w1

    g1 = 8 * sqrt((573 + 5 * sqrt(2865)) / (109825 + 969 * sqrt(2865)))
    g2 = sqrt(2 * (8025 + sqrt(2865)) / 35) / 37
    g3, g4 = [-i * (+i * 87 + sqrt(2865)) / 168 for i in [+1, -1]]

    data = [
        (w1, _s4(g1, g3)),
        (w2, _s4(g2, g4)),
        (frac(18, 5), numpy.array([[0, 0, frac(2, 3)]])),
    ]
    points, weights = untangle(data)
    return PyramidScheme("Felippa 6", weights, points, degree, citation)


def felippa_7(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    degree = 2
    w1 = frac(170569, 331200)
    w2 = frac(276710106577408, 1075923777052725)
    w3 = frac(12827693806929, 30577384040000)
    w4 = frac(10663383340655070643544192, 4310170528879365193704375)
    g1 = 7 * sqrt(frac(35, 59)) / 8
    g2 = 224 * sqrt(frac(336633710, 33088740423)) / 37
    g3 = sqrt(frac(37043, 35)) / 56
    g4 = -frac(127, 153)
    g5 = frac(1490761, 2842826)
    data = [
        (w1, _s4(g1, -frac(1, 7))),
        (w2, _s4_0(g2, -frac(9, 28))),
        (w3, _s4(g3, g4)),
        (w4, numpy.array([[0, 0, g5]])),
    ]
    points, weights = untangle(data)
    return PyramidScheme("Felippa 7", weights, points, degree, citation)


def felippa_8(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    wg9 = numpy.array([frac(64, 81), frac(40, 81), frac(25, 81)])

    degree = 3
    w1 = 5 * (68 + 5 * sqrt(10)) / 432
    w2 = frac(85, 54) - w1
    g1 = sqrt(frac(3, 5))
    g2 = 1 - 2 * (10 - sqrt(10)) / 15
    g3 = -frac(2, 3) - g2
    data = [
        (w1 * wg9[2], _s4(g1, g2)),
        (w1 * wg9[1], _s4_0(g1, g2)),
        (w1 * wg9[0], numpy.array([[0, 0, g2]])),
        (w2 * wg9[2], _s4(g1, g3)),
        (w2 * wg9[1], _s4_0(g1, g3)),
        (w2 * wg9[0], numpy.array([[0, 0, g3]])),
    ]
    points, weights = untangle(data)
    return PyramidScheme("Felippa 8", weights, points, degree, citation)


def felippa_9(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = numpy.vectorize(sympy.sqrt) if symbolic else numpy.sqrt

    wg9 = numpy.array([frac(64, 81), frac(40, 81), frac(25, 81)])

    degree = 5
    g1 = sqrt(frac(3, 5))
    g3 = -0.854011951853700535688324041975993416
    g4 = -0.305992467923296230556472913192103090
    g5 = +0.410004419776996766244796955168096505
    w1 = (
        frac(4, 15)
        * (4 + 5 * (g4 + g5) + 10 * g4 * g5)
        / ((g3 - g4) * (g3 - g5) * (1 - g3) ** 2)
    )
    w2 = (
        frac(4, 15)
        * (4 + 5 * (g3 + g5) + 10 * g3 * g5)
        / ((g3 - g4) * (g5 - g4) * (1 - g4) ** 2)
    )
    w3 = (
        frac(4, 15)
        * (4 + 5 * (g3 + g4) + 10 * g3 * g4)
        / ((g3 - g5) * (g4 - g5) * (1 - g5) ** 2)
    )
    data = [
        (w1 * wg9[2], _s4(g1, g3)),
        (w1 * wg9[1], _s4_0(g1, g3)),
        (w1 * wg9[0], numpy.array([[0, 0, g3]])),
        (w2 * wg9[2], _s4(g1, g4)),
        (w2 * wg9[1], _s4_0(g1, g4)),
        (w2 * wg9[0], numpy.array([[0, 0, g4]])),
        (w3 * wg9[2], _s4(g1, g5)),
        (w3 * wg9[1], _s4_0(g1, g5)),
        (w3 * wg9[0], numpy.array([[0, 0, g5]])),
    ]
    points, weights = untangle(data)
    return PyramidScheme("Felippa 9", weights, points, degree, citation)
