import math

import numpy as np
import sympy

from .. import un
from ..helpers import article, fsd, pm, untangle
from ._helpers import SnScheme

source = article(
    authors=["A.H. Stroud"],
    title="Some Seventh Degree Integration Formulas for Symmetric Regions",
    journal="SIAM J. Numer. Anal.",
    volume="4",
    number="1",
    pages="37–44",
    url="https://doi.org/10.1137/0704004",
)


def _stroud_1967_7_ab(n, variant_a, symbolic):
    sqrt = sympy.sqrt if symbolic else math.sqrt

    if variant_a:
        assert 3 <= n <= 7
        t = 1
    else:
        # ERR Stroud mentions nothing of variant b being only valid up to dimension 6,
        # but that's the way it is.
        assert 3 <= n <= 6
        t = -1

    alpha = sqrt(6 * (n + 6) * (8 - n))

    r2 = (3 * (n + 6) * (8 - n) - t * (n - 2) * alpha) / ((n + 6) * (34 - 5 * n))
    s2 = (3 * n * (n + 6) - t * 2 * alpha) / ((n + 6) * (3 * n ** 2 + 6 * n - 16))
    t2 = (6 * (n + 6) + t * alpha) / (14 * (n + 6))

    B = (8 - n) / (n + 2) / (n + 4) / (n + 6) / r2 ** 3
    C = 1 / (n + 2) / (n + 4) / (n + 6) / s2 ** 3 / 2 ** n
    D = 1 / (n + 2) / (n + 4) / (n + 6) / t2 ** 3 / 2
    A = 1 - 2 * n * B - 2 ** n * C - 2 * n * (n - 1) * D

    r = sqrt(r2)
    s = sqrt(s2)
    t = sqrt(t2)

    data = [(A, [n * [0]]), (B, fsd(n, (r, 1))), (C, pm(n * [s])), (D, fsd(n, (t, 2)))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)

    variant = "a" if variant_a else "b"
    return SnScheme(f"Stroud 1967-7{variant}", n, weights, points, 7, source)


def stroud_1967_7_a(n, symbolic=False):
    return _stroud_1967_7_ab(n, variant_a=True, symbolic=symbolic)


def stroud_1967_7_b(n, symbolic=False):
    return _stroud_1967_7_ab(n, variant_a=False, symbolic=symbolic)


def stroud_1967_7_c(n, symbolic=False):
    sqrt = sympy.sqrt if symbolic else math.sqrt
    assert n >= 3

    alpha = sqrt(2 * (n + 2) * (n + 4))

    r1, r2 = (
        sqrt(((n + 2) * (n + 4) + i * 2 * alpha) / (n + 4) / (n + 6)) for i in [+1, -1]
    )
    A1, A2 = (
        (2 * (n + 2) ** 2 + i * (n - 2) * alpha) / (4 * n * (n + 2) ** 2)
        for i in [+1, -1]
    )

    s = un.stroud_1967(n)

    points = np.concatenate([r1 * s.points, r2 * s.points])
    points = np.ascontiguousarray(points.T)
    weights = np.concatenate([A1 * s.weights, A2 * s.weights])

    # weights *= un.volume_nsphere(n - 1, symbolic) / volume_nball(n, symbolic)
    weights *= n

    return SnScheme("Stroud 1967-7 c", n, weights, points, 7, source)
