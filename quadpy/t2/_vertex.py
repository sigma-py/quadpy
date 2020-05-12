from sympy import Rational as frac

from ._helpers import T2Scheme, s2


def vertex():
    weights, points = s2([frac(1, 3), 0])
    return T2Scheme("Vertex scheme", weights, points, 1)
