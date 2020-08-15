from sympy import Rational as frac

from ._helpers import T2Scheme, register


def vertex():
    d = {"vertex": [[frac(1, 3)]]}
    return T2Scheme("Vertex scheme", d, 1)


register([vertex])
