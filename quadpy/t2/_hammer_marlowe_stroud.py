"""
Two of the schemes also appear in

P.C. Hammer, Arthur H. Stroud,
Numerical Evaluation of Multiple Integrals II,
Mathematical Tables and Other Aids to Computation.
Vol. 12, No. 64 (Oct., 1958), pp. 272-280,
<https://www.jstor.org/stable/2002370>
"""
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import T2Scheme, concat, r, s3

source = article(
    authors=["P.C. Hammer", "O.J. Marlowe", "A.H. Stroud"],
    title="Numerical Integration Over Simplexes and Cones",
    journal="Mathematical Tables and Other Aids to Computation",
    volume="10",
    number="55",
    month="jul",
    year="1956",
    pages="130-137",
    url="https://doi.org/10.1090/S0025-5718-1956-0086389-6",
)


def hammer_marlowe_stroud_1():
    weights, points = s3(1)
    return T2Scheme("Hammer-Marlowe-Stroud 1", weights, points, 1, source)


def hammer_marlowe_stroud_2():
    weights, points = r([frac(1, 3), frac(1, 2)])
    return T2Scheme("Hammer-Marlowe-Stroud 2", weights, points, 2, source)


def hammer_marlowe_stroud_3():
    weights, points = r([frac(1, 3), -frac(1, 2)])
    return T2Scheme("Hammer-Marlowe-Stroud 3", weights, points, 2, source)


def hammer_marlowe_stroud_4():
    weights, points = concat(s3(-frac(9, 16)), r([frac(25, 48), frac(2, 5)]))
    return T2Scheme("Hammer-Marlowe-Stroud 4", weights, points, 3, source)


def hammer_marlowe_stroud_5():
    w1, w2 = [(155 - i * sqrt(15)) / 1200 for i in [+1, -1]]
    x1, x2 = [(1 + i * sqrt(15)) / 7 for i in [+1, -1]]
    weights, points = concat(s3(frac(9, 40)), r([w1, x1], [w2, x2]))
    return T2Scheme("Hammer-Marlowe-Stroud 5", weights, points, 5, source)
