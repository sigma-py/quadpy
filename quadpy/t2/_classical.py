from sympy import Rational as frac

from ._helpers import T2Scheme, register


def centroid():
    d = {"centroid": [[1]]}
    return T2Scheme("Centroid rule", d, 1, tol=7.850e-17)


def vertex():
    d = {"vertex": [[frac(1, 3)]]}
    return T2Scheme("Vertex scheme", d, 1)


def seven_point():
    d = {
        "centroid": [[frac(9, 20)]],
        "vertex": [[frac(1, 20)]],
        "d3_aa": [[frac(2, 15)], [frac(1, 2)]],
    }
    return T2Scheme("Seven-point scheme", d, 3)


register([centroid, vertex, seven_point])
