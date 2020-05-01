import numpy
from sympy import sqrt, pi, Rational as frac
from ..helpers import article
from ..nsphere._mysovskikh import get_nsimplex_points
from ._helpers import Enr2Scheme
from ._stroud_secrest import stroud_secrest_4, _nsimplex
from ._stroud import stroud_enr2_5_1a, stroud_enr2_5_1b
from ._phillips import phillips

citation = article(
    authors=["James Lu", "David L. Darmofal"],
    title="Higher-Dimensional Integration with Gaussian Weight for Applications in Probabilistic Design",
    journal="SIAM J. Sci. Comput.",
    volume="26",
    number="2",
    year="2004",
    pages="613–624",
    url="https://doi.org/10.1137/S1064827503426863",
)


def lu_darmofal_1(n):
    assert n >= 4
    # a = get_nsimplex_points(n)
    a = _nsimplex(n)
    b = numpy.array([
        sqrt(frac(n, 2 * (n - 1))) * (a[k] + a[l])
        for k in range(n + 1)
        for l in range(k)
    ])
    points = numpy.concatenate([
        [[0] * n],
        +sqrt(frac(n, 2) - 1) * a,
        -sqrt(frac(n, 2) - 1) * a,
        +sqrt(frac(n, 2) - 1) * b,
        -sqrt(frac(n, 2) - 1) * b,
    ])

    p = frac(2, n + 2)
    A = frac(n ** 2 * (7 - n), 2 * (n + 1) ** 2 * (n + 2) ** 2)
    B = frac(2 * (n - 1) ** 2, (n + 1) ** 2 * (n + 2) ** 2)
    weights = numpy.concatenate([
        [p],
        numpy.full(len(a), A),
        numpy.full(len(a), A),
        numpy.full(len(b), B),
        numpy.full(len(b), B),
    ])
    weights *= sqrt(pi) ** n

    # ERR the article lists degree 5
    return Enr2Scheme("Lu-Darmofal I", n, weights, points, 1, citation)


lu_darmofal_2 = stroud_secrest_4
lu_darmofal_3 = phillips
lu_darmofal_4a = stroud_enr2_5_1a
lu_darmofal_4b = stroud_enr2_5_1b
