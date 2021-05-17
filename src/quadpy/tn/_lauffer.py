import numpy as np
from sympy import Rational as frac

from ..helpers import article, rd, untangle
from ._helpers import TnScheme

source = article(
    authors=["R. Lauffer"],
    title="Interpolation mehfacher Integrale",
    journal="Arch. Math.",
    volume="6",
    year="1955",
    pages="159-164",
    url="https://doi.org/10.1007/BF01900222",
)


def lauffer_1(n: int):
    data = [(frac(1, n + 1), rd(n + 1, [(1, 1)]))]
    points, weights = untangle(data)
    return TnScheme("Lauffer 1", n, weights, points, 1, source)


def lauffer_2(n: int):
    B = frac(2 - n, (n + 1) * (n + 2))
    C = frac(4, (n + 1) * (n + 2))
    data = [(B, rd(n + 1, [(1, 1)])), (C, rd(n + 1, [(frac(1, 2), 2)]))]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return TnScheme("Lauffer 2", n, weights, points, 2, source)


def lauffer_3(n: int):
    B = frac(n ** 2 - 4 * n + 6, (n + 1) * (n + 2) * (n + 3))
    C = frac(27 - 9 * n, 2 * (n + 1) * (n + 2) * (n + 3))
    D = frac(27, (n + 1) * (n + 2) * (n + 3))

    r = frac(1, 3)
    s = frac(2, 3)

    data = [
        (B, rd(n + 1, [(1, 1)])),
        (C, rd(n + 1, [(r, 1), (s, 1)])),
        (D, rd(n + 1, [(r, 3)])),
    ]
    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return TnScheme("Lauffer 3", n, weights, points, 3, source)


def lauffer_4(n: int):
    assert n >= 3

    nprod = (n + 1) * (n + 2) * (n + 3) * (n + 4)
    B1 = frac(-3 * n ** 3 + 17 * n ** 2 - 58 * n + 72, 3 * nprod)
    B2 = frac(16 * (n ** 2 - 5 * n + 12), 3 * nprod)
    B3 = frac(4 * (n ** 2 - 9 * n + 12), nprod)
    B4 = frac(64 * (4 - n), 2 * nprod)
    B5 = frac(256, nprod)

    r = frac(1, 4)
    s = frac(3, 4)
    t = frac(1, 2)

    data = [
        (B1, rd(n + 1, [(1, 1)])),
        (B2, rd(n + 1, [(r, 1), (s, 1)])),
        (B3, rd(n + 1, [(t, 2)])),
        (B4, rd(n + 1, [(r, 2), (t, 1)])),
        (B5, rd(n + 1, [(r, 4)])),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return TnScheme("Lauffer 4", n, weights, points, 4, source)


def lauffer_5(n: int):
    assert n >= 4

    nprod = (n + 1) * (n + 2) * (n + 3) * (n + 4) * (n + 5)

    # ERR Stroud is missing the factor 1/12 in B1.
    B1 = frac((12 * n ** 4 - 82 * n ** 3 + 477 * n ** 2 - 1277 * n + 1440), 12 * nprod)
    B2 = frac(5 ** 2 * (-3 * n ** 3 + 19 * n ** 2 - 96 * n + 170), 12 * nprod)
    B3 = frac(5 ** 2 * (-(n ** 3) + 13 * n ** 2 - 47 * n + 65), 6 * nprod)
    B4 = frac(5 ** 3 * (n ** 2 - 6 * n + 20), 3 * nprod)
    B5 = frac(5 ** 3 * (n ** 2 - 11 * n + 20), 4 * nprod)
    B6 = frac(5 ** 4 * (5 - n), 2 * nprod)
    B7 = frac(5 ** 5, nprod)

    r = frac(1, 5)
    s = frac(4, 5)
    u = frac(2, 5)
    v = frac(3, 5)

    data = [
        (B1, rd(n + 1, [(1, 1)])),
        (B2, rd(n + 1, [(r, 1), (s, 1)])),
        (B3, rd(n + 1, [(u, 1), (v, 1)])),
        (B4, rd(n + 1, [(r, 2), (v, 1)])),
        (B5, rd(n + 1, [(r, 1), (u, 2)])),
        (B6, rd(n + 1, [(r, 3), (u, 1)])),
        (B7, rd(n + 1, [(r, 5)])),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return TnScheme("Lauffer 5", n, weights, points, 5, source)
