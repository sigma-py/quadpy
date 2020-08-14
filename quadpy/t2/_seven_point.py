from sympy import Rational as frac

from ._helpers import T2Scheme


def seven_point():
    d = {"s3": [[frac(9, 20)]], "s2": [[frac(1, 20), frac(2, 15)], [0, frac(1, 2)]]}
    return T2Scheme("Seven-point scheme", d, 3)
