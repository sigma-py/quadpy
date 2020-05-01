import numpy
from sympy import Rational as frac
from sympy import pi, sqrt

from ..helpers import article
from ..nsphere._mysovskikh import get_nsimplex_points
from ._helpers import Enr2Scheme
from ._phillips import phillips as lu_darmofal_3
from ._stroud import stroud_enr2_5_1a as lu_darmofal_4a
from ._stroud import stroud_enr2_5_1b as lu_darmofal_4b
from ._stroud_secrest import stroud_secrest_4 as lu_darmofal_2

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
    # ENH The article says n>=4, but the scheme also works for 2, 3
    assert n >= 2
    a = get_nsimplex_points(n)
    b = numpy.array(
        [
            sqrt(frac(n, 2 * (n - 1))) * (a[k] + a[l])
            for k in range(len(a))
            for l in range(k)
        ]
    )
    points = numpy.concatenate(
        [
            [[0] * n],
            +sqrt(frac(n, 2) + 1) * a,
            -sqrt(frac(n, 2) + 1) * a,
            +sqrt(frac(n, 2) + 1) * b,
            -sqrt(frac(n, 2) + 1) * b,
        ]
    )

    p = frac(2, n + 2)
    A = frac(n ** 2 * (7 - n), 2 * (n + 1) ** 2 * (n + 2) ** 2)
    B = frac(2 * (n - 1) ** 2, (n + 1) ** 2 * (n + 2) ** 2)
    weights = numpy.concatenate(
        [
            [p],
            numpy.full(len(a), A),
            numpy.full(len(a), A),
            numpy.full(len(b), B),
            numpy.full(len(b), B),
        ]
    )
    weights *= sqrt(pi) ** n

    return Enr2Scheme("Lu-Darmofal I", n, weights, points, 5, citation)


__all__ = [
    "lu_darmofal_1",
    "lu_darmofal_2",
    "lu_darmofal_3",
    "lu_darmofal_4a",
    "lu_darmofal_4b",
]
