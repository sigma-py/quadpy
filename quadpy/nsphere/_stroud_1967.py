# -*- coding: utf-8 -*-
#

import numpy
import sympy

from ..helpers import untangle, fsd, pm, article
from ._helpers import integrate_monomial_over_unit_nsphere, NSphereScheme

citation = article(
    authors=["A.H. Stroud"],
    title="Some Seventh Degree Integration Formulas for Symmetric Regions",
    journal="SIAM J. Numer. Anal.",
    volume="4",
    number="1",
    pages="37–44",
    url="https://doi.org/10.1137/0704004",
)


def stroud_1967(n, symbolic=False):
    sqrt = sympy.sqrt if symbolic else numpy.sqrt
    frac = sympy.Rational if symbolic else lambda x, y: x / y

    degree = 7

    r = 1
    s = sqrt(frac(1, n))
    t = sqrt(frac(1, 2))
    B = frac(8 - n, n * (n + 2) * (n + 4))
    C = frac(n ** 3, 2 ** n * n * (n + 2) * (n + 4))
    D = frac(4, n * (n + 2) * (n + 4))

    data = [
        (B, fsd(n, (r, 1))),
        (C, pm(n, s)),
        # ERR Stroud's book wrongly states (t, t,..., t)_FS instead of
        # (t, t, 0, ..., 0)_FS.
        (D, fsd(n, (t, 2))),
    ]

    points, weights = untangle(data)
    weights *= integrate_monomial_over_unit_nsphere(n * [0], symbolic=symbolic)
    return NSphereScheme("Stroud 1967", n, weights, points, degree, citation)
