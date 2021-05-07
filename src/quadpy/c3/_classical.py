from sympy import Rational as frac

from ._helpers import C3Scheme, register


def midpoint():
    d = {"zero3": [[1]]}
    return C3Scheme("Midpoint", d, 1)


def vertex():
    d = {"symm_rrr": [[frac(1, 8)], [1]]}
    return C3Scheme("Vertex", d, 1)


def face_midpoint():
    d = {"symm_r00": [[frac(1, 6)], [1]]}
    return C3Scheme("Face-Midpoint", d, 3)


register([midpoint, face_midpoint, vertex])
