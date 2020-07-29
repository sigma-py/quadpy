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
from ._helpers import T3Scheme, concat, s31, s4

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
# source = book(
#     authors=["Olgierd Zienkiewicz"],
#     title="The Finite Element Method, Sixth Edition",
#     publisher="Butterworth-Heinemann",
#     year="2005",
#     isbn="0750663200",
#     url="http://www.sciencedirect.com/science/book/9780750664318",
# )


# Used in Zienkiewicz 4
def hammer_marlowe_stroud_1():
    degree = 2
    r = 1 / sqrt(5)
    a = (1 - r) / 4
    weights, points = s31([frac(1, 4), a])
    return T3Scheme("Hammer-Marlowe-Stroud 1", weights, points, degree, source)


def hammer_marlowe_stroud_2():
    degree = 2
    r = - 1 / sqrt(5)
    a = (1 - r) / 4
    weights, points = s31([frac(1, 4), a])
    return T3Scheme("Hammer-Marlowe-Stroud 2", weights, points, degree, source)


# Used in Zienkiewicz 5
def hammer_marlowe_stroud_3():
    degree = 3
    r = frac(1, 3)
    a = (1 - r) / 4
    weights, points = concat(s4(-frac(4, 5)), s31([frac(9, 20), a]))
    return T3Scheme("Hammer-Marlowe-Stroud 3", weights, points, degree, source)
