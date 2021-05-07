import numpy as np
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ..un._mysovskikh import get_nsimplex_points
from ._helpers import Enr2Scheme
from ._phillips import phillips as lu_darmofal_3
from ._stroud import stroud_enr2_5_1a as lu_darmofal_4a
from ._stroud import stroud_enr2_5_1b as lu_darmofal_4b
from ._stroud_secrest import stroud_secrest_4 as lu_darmofal_2

source = article(
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
    a = get_nsimplex_points(n, sqrt, frac)
    b = np.array(
        [
            sqrt(frac(n, 2 * (n - 1))) * (a[k] + a[l])
            for k in range(len(a))
            for l in range(k)
        ]
    )
    points = np.concatenate(
        [
            [[0] * n],
            +sqrt(frac(n, 2) + 1) * a,
            -sqrt(frac(n, 2) + 1) * a,
            +sqrt(frac(n, 2) + 1) * b,
            -sqrt(frac(n, 2) + 1) * b,
        ]
    )
    points = np.ascontiguousarray(points.T)

    p = frac(2, n + 2)
    A = frac(n ** 2 * (7 - n), 2 * (n + 1) ** 2 * (n + 2) ** 2)
    B = frac(2 * (n - 1) ** 2, (n + 1) ** 2 * (n + 2) ** 2)
    weights = np.concatenate(
        [
            [p],
            np.full(len(a), A),
            np.full(len(a), A),
            np.full(len(b), B),
            np.full(len(b), B),
        ]
    )
    return Enr2Scheme("Lu-Darmofal I", n, weights, points, 5, source)


__all__ = [
    "lu_darmofal_1",
    "lu_darmofal_2",
    "lu_darmofal_3",
    "lu_darmofal_4a",
    "lu_darmofal_4b",
]
