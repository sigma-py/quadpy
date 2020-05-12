from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import TetrahedronScheme, concat, s4, s31, s211

citation = article(
    authors=["Yu Jinyun"],
    title="Symmetyric Gaussian quadrature formulae for tetrahedronal regions",
    journal="Computer Methods in Applied Mechanics and Engineering",
    volume="43",
    year="1984",
    pages="349-353",
    url="https://doi.org/10.1016/0045-7825(84)90072-0",
)


def yu_1():
    degree = 2
    weights, points = s31([frac(1, 4), frac(1, 4) - sqrt(5) * frac(1, 20)])
    return TetrahedronScheme("Yu 1", weights, points, degree, citation)


def yu_2():
    degree = 3
    weights, points = concat(s4(-frac(4, 5)), s31([frac(9, 20), frac(1, 6)]))
    return TetrahedronScheme("Yu 2", weights, points, degree, citation)


def yu_3():
    degree = 4
    weights, points = concat(
        s31([0.5037379410012282e-01, 0.7611903264425430e-01]),
        s211([0.6654206863329239e-01, 0.4042339134672644, 0.1197005277978019]),
    )
    return TetrahedronScheme("Yu 3", weights, points, degree, citation)


def yu_4():
    degree = 5
    weights, points = concat(
        s4(0.1884185567365411),
        s31([0.6703858372604275e-01, 0.8945436401412733e-01]),
        s211([0.4528559236327399e-01, 0.4214394310662522, 0.1325810999384657]),
    )
    return TetrahedronScheme("Yu 4", weights, points, degree, citation)


def yu_5():
    degree = 6
    weights, points = concat(
        s4(0.9040129046014750e-01),
        s31([0.1911983427899124e-01, 0.5742691731735682e-01]),
        s211(
            [0.4361493840666568e-01, 0.2312985436519147, 0.5135188412556341e-01],
            [0.2581167596199161e-01, 0.4756909881472290e-01, 0.2967538129690260],
        ),
    )
    return TetrahedronScheme("Yu 5", weights, points, degree, citation)
