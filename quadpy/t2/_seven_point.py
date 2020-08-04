from sympy import Rational as frac

from ._helpers import T2Scheme, expand_symmetries


def seven_point():
    d = {"s3": [[frac(9, 20)]], "s2": [[frac(1, 20), frac(2, 15)], [0, frac(1, 2)]]}
    points, weights = expand_symmetries(d)
    return T2Scheme("Seven-point scheme", weights, points, 3)
