from sympy import Rational as frac
from sympy import sqrt

from ..helpers import techreport
from ._helpers import T2Scheme, expand_symmetries

source = techreport(
    authors=["Noel J. Walkington"],
    title="Quadrature on simplices of arbitrary dimension",
    institution="CMU",
    year="2000",
    url="https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/",
)


def walkington_p5():
    a1, a2 = [(155 + i * sqrt(15)) / 1200 for i in [+1, -1]]
    x1, x2 = [(6 + i * sqrt(15)) / 21 for i in [+1, -1]]
    d = {"s3": [[frac(9, 40)]], "s2": [[a1, a2], [x1, x2]]}
    points, weights = expand_symmetries(d)
    return T2Scheme("Walkington p5", weights, points, 5, source)
