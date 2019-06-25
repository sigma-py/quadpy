# -*- coding: utf-8 -*-
#

import sympy

from ._helpers import untangle2, TetrahedronScheme
from ..helpers import article


citation = article(
    authors=["P. Keast"],
    title="Moderate degree tetrahedral quadrature formulas",
    journal="Computer Methods in Applied Mechanics and Engineering",
    volume="55",
    number="3",
    pages="339-348",
    month="may",
    year="1986",
    url="https://doi.org/10.1016/0045-7825(86)90059-9",
)

# https://people.sc.fsu.edu/~jburkardt/datasets/quadrature_rules_tet/quadrature_rules_tet.html


def keast_0(symbolic=False):
    # Does not appear in Keast's article.
    degree = 1
    data = {"s4": [[1.0]]}
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 0", weights, points, degree, citation)


def keast_1(symbolic=False):
    # Does not appear in Keast's article.
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 2
    data = {"s31": [[frac(1, 4), 0.1381966011250105]]}
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 1", weights, points, degree, citation)


def keast_2(symbolic=False):
    # Does not appear in Keast's article.
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 3
    data = {"s4": [[-frac(4, 5)]], "s31": [[frac(9, 20), frac(1, 6)]]}
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 2", weights, points, degree, citation)


def keast_3(symbolic=False):
    # Does not appear in Keast's article.
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 3
    data = {
        "s31": [[0.2177650698804054, 0.1438564719343852]],
        "s22": [[0.0214899534130631, frac(1, 2)]],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 3", weights, points, degree, citation)


def keast_4(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 4
    data = {
        "s4": [[-frac(148, 1875)]],
        "s31": [[+frac(343, 7500), frac(1, 14)]],
        "s22": [[+frac(56, 375), 0.3994035761667992]],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 4", weights, points, degree, citation)


def keast_5(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 4
    data = {
        "s22": [[frac(2, 105), frac(1, 2)]],
        "s31": [
            [0.0885898247429807, 0.1005267652252045],
            [0.1328387466855907, 0.3143728734931922],
        ],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 5", weights, points, degree, citation)


def keast_6(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 5
    data = {
        "s4": [[frac(6544, 36015)]],
        "s31": [[frac(81, 2240), frac(1, 3)], [frac(161051, 2304960), frac(1, 11)]],
        "s22": [[frac(338, 5145), 0.0665501535736643]],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 6", weights, points, degree, citation)


def keast_7(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 6
    data = {
        "s31": [
            [0.0399227502581679, 0.2146028712591517],
            [0.0100772110553207, 0.0406739585346113],
            [0.0553571815436544, 0.3223378901422757],
        ],
        "s211": [[frac(27, 560), 0.0636610018750175, 0.2696723314583159]],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 7", weights, points, degree, citation)


def keast_8(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 7
    data = {
        "s4": [[+0.1095853407966528]],
        "s31": [
            [+0.0635996491464850, 0.0782131923303186],
            [-0.3751064406859797, 0.1218432166639044],
            [+0.0293485515784412, 0.3325391644464206],
        ],
        "s22": [[+0.0058201058201058, frac(1, 2)]],
        "s211": [[+0.1653439153439105, frac(1, 10), frac(1, 5)]],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 8", weights, points, degree, citation)


def keast_9():
    degree = 8
    data = {
        "s4": [[-0.2359620398477557]],
        "s31": [
            [+0.0244878963560562, 0.1274709365666390],
            [+0.0039485206398261, 0.0320788303926323],
        ],
        "s22": [
            [+0.0263055529507371, 0.0497770956432810],
            [+0.0829803830550589, 0.1837304473985499],
        ],
        "s211": [
            [+0.0254426245481023, 0.2319010893971509, 0.5132800333608811],
            [+0.0134324384376852, 0.0379700484718286, 0.1937464752488044],
        ],
    }
    points, weights = untangle2(data)
    return TetrahedronScheme("Keast 9", weights, points, degree, citation)
