# -*- coding: utf-8 -*-
#
"""
"""

import sympy

from ._helpers import TriangleScheme, s1, s2, concat
from ..helpers import article


citation = article(
    authors=["Mark A. Taylor", "Beth A. Wingate", "Len P. Bos"],
    title="Several new quadrature formulas for polynomial integration in the triangle",
    journal="arXiv Mathematics e-prints",
    year="2005",
    month="jan",
    url="https://arxiv.org/abs/math/0501496",
)


# TODO missing Taylor-Wingate-Bos schemes


def taylor_wingate_bos_1(symbolic=False):
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    weights, points = s2([frac(2, 3), frac(1, 6)])
    weights /= 2
    return TriangleScheme("Taylor-Wingate-Bos 1", weights, points, 2, citation)


def taylor_wingate_bos_2():
    weights, points = s2(
        [0.2199034873106, 0.0915762135098], [0.4467631793560, 0.4459484909160]
    )
    weights /= 2
    return TriangleScheme("Taylor-Wingate-Bos 2", weights, points, 4, citation)


def taylor_wingate_bos_4():
    weights, points = concat(
        s2(
            [0.0102558174092, 0],
            [0.1679775595335, 0.4743880861752],
            [0.2652238803946, 0.2385615300181],
        ),
        s1([0.1116047046647, 0.7839656651012, 0.0421382841642]),
    )
    weights /= 2
    return TriangleScheme("Taylor-Wingate-Bos 4", weights, points, 7, citation)


def taylor_wingate_bos_5():
    weights, points = concat(
        s2(
            [0.0519871420646, 0.0451890097844],
            [0.1032344051380, 0.4815198347833],
            [0.1881601469167, 0.4036039798179],
        ),
        s1(
            [0.0707034101784, 0.7475124727339, 0.0304243617288],
            [0.0909390760952, 0.1369912012649, 0.2182900709714],
        ),
    )
    weights /= 2
    return TriangleScheme("Taylor-Wingate-Bos 5", weights, points, 9, citation)


def taylor_wingate_bos_8(symbolic=False):
    weights, points = concat(
        s2(
            [0.0010616711990, 0],
            [0.0349317947036, 0.4903668903754],
            [0.0383664533945, 0.0875134669581],
            [0.0897856524107, 0.2217145894873],
            [0.1034544533617, 0.3860471669296],
        ),
        s1(
            [0.0131460236101, 0.0573330873026, 0.0151382269814],
            [0.0242881926949, 0.8159625040711, 0.1659719969565],
            [0.0316799866332, 0.3165475556378, 0.0186886898773],
            [0.0578369491210, 0.0935526036219, 0.2079865423167],
            [0.0725821687394, 0.0974892983467, 0.5380088595149],
        ),
    )
    weights /= 2
    return TriangleScheme("Taylor-Wingate-Bos 8", weights, points, 14, citation)
