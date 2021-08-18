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
from ._helpers import T2Scheme, register

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
    d = {"centroid": [[1]]}
    return T2Scheme("Hammer-Marlowe-Stroud 1", d, 1, source, 7.850e-17)


def hammer_marlowe_stroud_2():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 6)]]}
    return T2Scheme("Hammer-Marlowe-Stroud 2", d, 2, source, 2.220e-16)


def hammer_marlowe_stroud_3():
    d = {"d3_aa": [[frac(1, 3)], [frac(1, 2)]]}
    return T2Scheme("Hammer-Marlowe-Stroud 3", d, 2, source, 3.074e-16)


def hammer_marlowe_stroud_4():
    d = {"centroid": [[-frac(9, 16)]], "d3_aa": [[frac(25, 48)], [frac(1, 5)]]}
    return T2Scheme("Hammer-Marlowe-Stroud 4", d, 3, source, 6.661e-16)


def hammer_marlowe_stroud_5():
    w1, w2 = ((155 - i * sqrt(15)) / 1200 for i in [+1, -1])
    x1, x2 = ((1 + i * sqrt(15)) / 7 for i in [+1, -1])
    b1 = (1 - x1) / 3
    b2 = (1 - x2) / 3
    d = {"centroid": [[frac(9, 40)]], "d3_aa": [[w1, w2], [b1, b2]]}
    return T2Scheme("Hammer-Marlowe-Stroud 5", d, 5, source, 2.776e-16)


register(
    [
        hammer_marlowe_stroud_1,
        hammer_marlowe_stroud_2,
        hammer_marlowe_stroud_3,
        hammer_marlowe_stroud_4,
        hammer_marlowe_stroud_5,
    ]
)
