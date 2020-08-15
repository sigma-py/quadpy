from sympy import Rational as frac

from ._helpers import T2Scheme, register


def seven_point():
    d = {"centroid": [[frac(9, 20)]], "d3_aa": [[frac(1, 20), frac(2, 15)], [0, frac(1, 2)]]}
    return T2Scheme("Seven-point scheme", d, 3)


register([seven_point])
