from sympy import Rational as frac

from ._helpers import T2Scheme


def vertex():
    d = {"s2": [[frac(1, 3)], [0]]}
    return T2Scheme("Vertex scheme", d, 1)
