# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ..helpers import article
from ._helpers import TetrahedronScheme, _s4

citation = article(
    authors=["Y. Liu", "M. Vinokur"],
    title="Exact Integrations of Polynomials and Symmetric Quadrature Formulas over Arbitrary Polyhedral Grids",
    journal="Journal of Computational Physics",
    volume="140",
    pages="122–147",
    year="1998",
    url="https://doi.org/10.1006/jcph.1998.5884",
)
# TODO update weight/points specification


def liu_vinokur_01(symbolic=False):
    weights = numpy.concatenate([numpy.full(1, 1)])
    points = numpy.concatenate([_s4()])
    degree = 1
    return TetrahedronScheme("Liu-Vinokur 1", weights, points, degree, citation)


def liu_vinokur_02(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights = numpy.concatenate([numpy.full(4, frac(1, 4))])
    points = numpy.concatenate([_r_alpha(1.0)])
    degree = 1
    return TetrahedronScheme("Liu-Vinokur 2", weights, points, degree, citation)


def liu_vinokur_03(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    weights = numpy.concatenate([numpy.full(4, frac(1, 4))])
    points = numpy.concatenate([_r_alpha(1 / sqrt(5))])
    degree = 2
    return TetrahedronScheme("Liu-Vinokur 3", weights, points, degree, citation)


def liu_vinokur_04(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights = numpy.concatenate([numpy.full(1, frac(4, 5)), numpy.full(4, frac(1, 20))])
    points = numpy.concatenate([_s4(), _r_alpha(1)])
    degree = 2
    return TetrahedronScheme("Liu-Vinokur 4", weights, points, degree, citation)


def liu_vinokur_05(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights = numpy.concatenate(
        [numpy.full(1, -frac(4, 5)), numpy.full(4, frac(9, 20))]
    )
    points = numpy.concatenate([_s4(), _r_alpha(frac(1, 3))])
    degree = 3
    return TetrahedronScheme("Liu-Vinokur 5", weights, points, degree, citation)


def liu_vinokur_06(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights = numpy.concatenate(
        [numpy.full(4, frac(1, 40)), numpy.full(4, frac(9, 40))]
    )
    points = numpy.concatenate([_r_alpha(1), _r_alpha(-frac(1, 3))])
    degree = 3
    return TetrahedronScheme("Liu-Vinokur 6", weights, points, degree, citation)


def liu_vinokur_07(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    weights = numpy.concatenate(
        [
            numpy.full(1, -frac(148, 1875)),
            numpy.full(4, frac(343, 7500)),
            numpy.full(6, frac(56, 375)),
        ]
    )
    points = numpy.concatenate([_s4(), _r_alpha(frac(5, 7)), _r_beta(sqrt(70) / 28)])
    degree = 4
    return TetrahedronScheme("Liu-Vinokur 7", weights, points, degree, citation)


def liu_vinokur_08(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    alpha1 = (+sqrt(65944 - 19446 * sqrt(11)) + 51 * sqrt(11) - 154) / 89
    alpha2 = (-sqrt(65944 - 19446 * sqrt(11)) + 51 * sqrt(11) - 154) / 89
    weights = numpy.concatenate(
        [
            numpy.full(4, (17 * alpha2 - 7) / (420 * alpha1 ** 2 * (alpha2 - alpha1))),
            numpy.full(4, (17 * alpha1 - 7) / (420 * alpha2 ** 2 * (alpha1 - alpha2))),
            numpy.full(6, frac(2, 105)),
        ]
    )
    points = numpy.concatenate(
        [_r_alpha(alpha1), _r_alpha(alpha2), _r_beta(frac(1, 2))]
    )
    degree = 4
    return TetrahedronScheme("Liu-Vinokur 8", weights, points, degree, citation)


def liu_vinokur_09(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights = numpy.concatenate(
        [
            numpy.full(1, -frac(32, 15)),
            numpy.full(4, frac(3, 280)),
            numpy.full(4, frac(125, 168)),
            numpy.full(6, frac(2, 105)),
        ]
    )
    points = numpy.concatenate(
        [_s4(), _r_alpha(1), _r_alpha(frac(1, 5)), _r_beta(frac(1, 2))]
    )
    degree = 4
    return TetrahedronScheme("Liu-Vinokur 9", weights, points, degree, citation)


def liu_vinokur_10(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    weights = numpy.concatenate(
        [
            numpy.full(1, frac(32, 105)),
            numpy.full(4, -frac(31, 840)),
            numpy.full(4, frac(27, 280)),
            numpy.full(12, frac(4, 105)),
        ]
    )
    points = numpy.concatenate(
        [
            _s4(),
            _r_alpha(1),
            _r_alpha(-frac(1, 3)),
            _r_gamma_delta((2 + sqrt(2)) / 4, (2 - sqrt(2)) / 4),
        ]
    )
    degree = 4
    return TetrahedronScheme("Liu-Vinokur 10", weights, points, degree, citation)


def liu_vinokur_11(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    weights = numpy.concatenate(
        [
            (11 - 4 * sqrt(2)) / numpy.full(4, 840),
            (243 - 108 * sqrt(2)) / numpy.full(4, 1960),
            (62 + 44 * sqrt(2)) / numpy.full(4, 735),
            numpy.full(6, frac(2, 105)),
        ]
    )
    points = numpy.concatenate(
        [_r_alpha(1), _r_alpha(-frac(1, 3)), _r_alpha(sqrt(2) - 1), _r_beta(frac(1, 2))]
    )
    degree = 4
    return TetrahedronScheme("Liu-Vinokur 11", weights, points, degree, citation)


def liu_vinokur_12(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    cos = sympy.cos if symbolic else numpy.cos
    acos = sympy.acos if symbolic else numpy.arccos

    lmbda = frac(4, 27) * (
        4 * sqrt(79) * cos((acos(67 * sqrt(79) / 24964) + 2 * numpy.pi) / 3) + 71
    )
    alpha1 = (+sqrt(9 * lmbda ** 2 - 248 * lmbda + 1680) + 28 - 3 * lmbda) / (
        112 - 10 * lmbda
    )
    alpha2 = (-sqrt(9 * lmbda ** 2 - 248 * lmbda + 1680) + 28 - 3 * lmbda) / (
        112 - 10 * lmbda
    )
    w1 = ((21 - lmbda) * alpha2 - 7) / (420 * alpha1 ** 2 * (alpha2 - alpha1))
    w2 = ((21 - lmbda) * alpha1 - 7) / (420 * alpha2 ** 2 * (alpha1 - alpha2))
    weights = numpy.concatenate(
        [numpy.full(4, w1), numpy.full(4, w2), numpy.full(6, lmbda ** 2 / 840)]
    )
    points = numpy.concatenate(
        [_r_alpha(alpha1), _r_alpha(alpha2), _r_beta(1 / sqrt(lmbda))]
    )
    degree = 5

    return TetrahedronScheme("Liu-Vinokur 12", weights, points, degree, citation)


def liu_vinokur_13(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y
    sqrt = sympy.sqrt if symbolic else numpy.sqrt

    weights = numpy.concatenate(
        [
            numpy.full(1, -frac(16, 21)),
            numpy.full(4, (2249 - 391 * sqrt(13)) / 10920),
            numpy.full(4, (2249 + 391 * sqrt(13)) / 10920),
            numpy.full(6, frac(2, 105)),
        ]
    )
    points = numpy.concatenate(
        [
            _s4(),
            _r_alpha((2 + sqrt(13)) / 9),
            _r_alpha((2 - sqrt(13)) / 9),
            _r_beta(frac(1, 2)),
        ]
    )
    degree = 5

    return TetrahedronScheme("Liu-Vinokur 13", weights, points, degree, citation)


def liu_vinokur_14(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights = numpy.concatenate(
        [
            numpy.full(1, frac(16, 105)),
            numpy.full(4, frac(1, 280)),
            numpy.full(4, frac(81, 1400)),
            numpy.full(4, frac(64, 525)),
            numpy.full(6, frac(2, 105)),
        ]
    )
    points = numpy.concatenate(
        [
            _s4(),
            _r_alpha(1),
            _r_alpha(-frac(1, 3)),
            _r_alpha(frac(1, 2)),
            _r_beta(frac(1, 2)),
        ]
    )
    degree = 5

    return TetrahedronScheme("Liu-Vinokur 14", weights, points, degree, citation)


def _r_alpha(alpha):
    """From the article:

    mu_i = (1 + (n-1) alpha) / n,
    mu_j = (1 - alpha) / n    for j!=i,

    where n is the number of vertices
    """
    a = (1 + 3 * alpha) / 4
    b = (1 - alpha) / 4
    return numpy.array([[a, b, b, b], [b, a, b, b], [b, b, a, b], [b, b, b, a]])


def _r_beta(beta):
    """From the article:

    mu_i = (1+(n-2)*beta) / n,
    mu_j = mu_i,
    mu_k = (1 - 2*beta) / n    for k!=i, k!=j,

    where n is the number of vertices.
    """
    a = (1 + 2 * beta) / 4
    b = (1 - 2 * beta) / 4
    return numpy.array(
        [
            [a, a, b, b],
            [a, b, a, b],
            [b, a, a, b],
            [a, b, b, a],
            [b, a, b, a],
            [b, b, a, a],
        ]
    )


def _r_gamma_delta(gamma, delta):
    """From the article:

    mu_i = (1 + (n-1) gamma - delta) / n,
    mu_j = (1 + (n-1) delta - gamma) / n,
    mu_k = (1 - gamma - delta) / n    for k!=i, k!=j,

    where n is the number of vertices
    """
    b = (1 + 3 * gamma - delta) / 4
    c = (1 + 3 * delta - gamma) / 4
    a = (1 - gamma - delta) / 4
    return numpy.array(
        [
            [a, a, b, c],
            [a, b, a, c],
            [b, a, a, c],
            [a, b, c, a],
            [b, a, c, a],
            [b, c, a, a],
            [a, a, c, b],
            [a, c, a, b],
            [c, a, a, b],
            [a, c, b, a],
            [c, a, b, a],
            [c, b, a, a],
        ]
    )
