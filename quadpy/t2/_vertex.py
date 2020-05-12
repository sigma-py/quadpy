from sympy import Rational as frac

from ._helpers import TriangleScheme, s2


def vertex():
    weights, points = s2([frac(1, 3), 0])
    return TriangleScheme("Vertex scheme", weights, points, 1)
