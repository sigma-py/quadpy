from sympy import Rational as frac
from sympy import sqrt

from ..helpers import techreport
from ._helpers import T2Scheme, register

source = techreport(
    authors=["Noel J. Walkington"],
    title="Quadrature on simplices of arbitrary dimension",
    institution="CMU",
    year="2000",
    url="https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/",
)


def walkington_p5():
    a1, a2 = ((155 + i * sqrt(15)) / 1200 for i in [+1, -1])
    x1, x2 = ((6 + i * sqrt(15)) / 21 for i in [+1, -1])
    d = {"centroid": [[frac(9, 40)]], "d3_aa": [[a1, a2], [x1, x2]]}
    return T2Scheme("Walkington p5", d, 5, source, 2.637e-16)


register([walkington_p5])
