from sympy import Rational as frac

from ._helpers import TriangleScheme, concat, s2, s3


def seven_point():
    weights, points = concat(
        s3(frac(9, 20)), s2([frac(1, 20), 0], [frac(2, 15), frac(1, 2)])
    )
    return TriangleScheme("Seven-point scheme", weights, points, 3)
