from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import T3Scheme, concat, s4, s22, s31, s211

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


def keast_0():
    # Does not appear in Keast's article. TODO remove
    degree = 1
    weights, points = s4(1)
    return T3Scheme("Keast 0", weights, points, degree, citation)


def keast_1():
    # Does not appear in Keast's article.
    degree = 2
    weights, points = s31([frac(1, 4), frac(1, 4) - sqrt(5) / 20])
    return T3Scheme("Keast 1", weights, points, degree, citation)


def keast_2():
    # Does not appear in Keast's article.
    degree = 3
    weights, points = concat(s4(-frac(4, 5)), s31(frac(9, 20), frac(1, 6)))
    return T3Scheme("Keast 2", weights, points, degree, citation)


def keast_3():
    # Does not appear in Keast's article.
    degree = 3
    weights, points = concat(
        s31([0.2177650698804054, 0.1438564719343852]), s22([0.0214899534130631, 0.5])
    )
    return T3Scheme("Keast 3", weights, points, degree, citation)


def keast_4():
    degree = 4
    weights, points = concat(
        s4(-frac(148, 1875)),
        s31([+frac(343, 7500), frac(1, 14)]),
        s22([+frac(56, 375), frac(1, 4) + sqrt(frac(5, 14)) / 4]),
    )
    return T3Scheme("Keast 4", weights, points, degree, citation)


def keast_5():
    degree = 4
    weights, points = concat(
        s22([2.0 / 105.0, 0.5]),
        s31(
            [0.0885898247429807, 0.1005267652252045],
            [0.1328387466855907, 0.3143728734931922],
        ),
    )
    return T3Scheme("Keast 5", weights, points, degree, citation)


def keast_6():
    degree = 5
    weights, points = concat(
        s4(frac(6544, 36015)),
        s31([frac(81, 2240), frac(1, 3)], [frac(161051, 2304960), frac(1, 11)]),
        s22([frac(338, 5145), frac(1, 4) - sqrt(91) / 52]),
    )
    return T3Scheme("Keast 6", weights, points, degree, citation)


def keast_7():
    degree = 6
    weights, points = concat(
        s31(
            [0.0399227502581679, 0.2146028712591517],
            [0.0100772110553207, 0.0406739585346113],
            [0.0553571815436544, 0.3223378901422757],
        ),
        s211([27.0 / 560.0, 0.0636610018750175, 0.2696723314583159]),
    )
    return T3Scheme("Keast 7", weights, points, degree, citation)


def keast_8():
    degree = 7
    weights, points = concat(
        s4(+0.1095853407966528),
        s31(
            [+0.0635996491464850, 0.0782131923303186],
            [-0.3751064406859797, 0.1218432166639044],
            [+0.0293485515784412, 0.3325391644464206],
        ),
        s22([+0.0058201058201058, 0.5]),
        s211([+0.1653439153439105, 0.1, 0.2]),
    )
    return T3Scheme("Keast 8", weights, points, degree, citation)


def keast_9():
    degree = 8
    weights, points = concat(
        s4(-0.2359620398477557),
        s31(
            [+0.0244878963560562, 0.1274709365666390],
            [+0.0039485206398261, 0.0320788303926323],
        ),
        s22(
            [+0.0263055529507371, 0.0497770956432810],
            [+0.0829803830550589, 0.1837304473985499],
        ),
        s211(
            [+0.0254426245481023, 0.2319010893971509, 0.5132800333608811],
            [+0.0134324384376852, 0.0379700484718286, 0.1937464752488044],
        ),
    )
    return T3Scheme("Keast 9", weights, points, degree, citation)
