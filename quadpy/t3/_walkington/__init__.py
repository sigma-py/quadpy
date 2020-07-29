from ...helpers import techreport
from .._helpers import T3Scheme, s31, s22, concat

source = techreport(
    authors=["Noel J. Walkington"],
    title="Quadrature on simplices of arbitrary dimension",
    institution="CMU",
    year="2000",
    url="https://www.math.cmu.edu/~nw0z/publications/00-CNA-023/023abs/",
)


def walkington_p5():
    degree = 5
    weights, points = concat(
        s31(0.018781320953002641800, 0.31088591926330060980),
        s31(0.012248840519393658257, 0.092735250310891226402),
        s22(0.0070910034628469110730, 0.045503704125649649492),
    )
    weights *= 6.0
    return T3Scheme("Walkington p5", weights, points, degree, source)
