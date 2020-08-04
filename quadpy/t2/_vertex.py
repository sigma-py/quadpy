from sympy import Rational as frac

from ._helpers import T2Scheme, expand_symmetries


def vertex():
    d = {"s2": [[frac(1, 3)], [0]]}
    points, weights = expand_symmetries(d)
    return T2Scheme("Vertex scheme", weights, points, 1)
