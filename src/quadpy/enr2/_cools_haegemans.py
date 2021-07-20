import math

import numpy as np
import sympy

from ..cn._cools_haegemans import _gener
from ..helpers import article, fsd, pm, untangle, z
from ._helpers import Enr2Scheme

_source = article(
    authors=["Ronald Cools", "Ann Haegemans"],
    title="An imbedded family of cubature formulae for n-dimensional product regions",
    journal="Journal of Computational and Applied Mathematics",
    volume="51",
    year="1994",
    pages="251-262",
    url="https://doi.org/10.1016/0377-0427%2892%2900007-V",
)


def _mu(j, symbolic):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    # 1/sqrt(pi) ** n int int ... int exp(-r ** 2) dr
    if j == 0:
        return 1
    elif j == 1:
        return 0
    return frac(j - 1, 2) * _mu(j - 2, symbolic)


def cools_haegemans_1(n, delta2=1, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt
    assert frac(1, 2) <= delta2
    m = 1

    w0 = frac(2 * delta2 - 1, 2 * delta2)
    w = frac(_mu(2, symbolic) ** m * _mu(0, symbolic) ** (n - m), 2 ** n * delta2 ** m)

    data = [
        (w0, z(n)),
        (w, pm(n * [sqrt(delta2)])),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Cools-Haegemans 1", n, weights, points, 3, _source)


def cools_haegemans_2(n, delta2=1, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt
    assert n >= 1
    assert frac(1, 2) <= delta2
    if n > 2:
        assert delta2 <= frac(n + 2, 2 * n - 4)
    m = 2

    lmbdas2 = _gener(delta2, 1, _mu, symbolic)

    w0 = frac(
        -(2 * delta2 - 1) * (2 * delta2 * n - 4 * delta2 - n - 2), 8 * delta2 ** 2
    )
    w1 = frac((2 * delta2 - 1) ** 2, 16 * delta2 ** 2)
    w = frac(_mu(2, symbolic) ** m * _mu(0, symbolic) ** (n - m), 2 ** n * delta2 ** m)

    lmbdas = [sqrt(lmbda2) for lmbda2 in lmbdas2]

    data = [
        (w0, z(n)),
        (w1, fsd(n, (lmbdas[0], 1))),
        (w, pm(n * [sqrt(delta2)])),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Cools-Haegemans 2", n, weights, points, 5, _source)


def cools_haegemans_3(n, delta2=(2, 3), symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt

    if isinstance(delta2, tuple):
        delta2 = frac(*delta2)

    assert n >= 2
    m = 3

    lmbdas2 = _gener(delta2, 2, _mu, symbolic)

    delta4 = delta2 ** 2
    delta6 = delta2 ** 3

    w0 = frac(
        (
            delta6 * (16 * n ** 2 - 72 * n + 128)
            + delta4 * (-28 * n ** 2 + 48 * n - 32)
            + delta2 * (16 * n ** 2 + 30 * n - 16)
            + (-3 * n ** 2 - 18 * n - 24)
        )
        * (2 * delta2 - 1),
        64 * (4 * delta2 - 3) * delta6,
    )
    w1 = -frac(
        (12 * delta4 * n - 32 * delta4 - 12 * delta2 * n + 3 * n + 12)
        * (2 * delta2 - 1),
        192 * delta6,
    )
    w2 = frac((2 * delta2 - 1) * (2 * delta2 - 3) ** 3, 384 * (4 * delta2 - 3) * delta6)
    w11 = frac((2 * delta2 - 1) ** 3, 128 * delta6)

    w = frac(_mu(2, symbolic) ** m * _mu(0, symbolic) ** (n - m), 2 ** n * delta2 ** m)

    delta = sqrt(delta2)
    lmbdas = [sqrt(lmbda2) for lmbda2 in lmbdas2]

    data = [
        (w0, z(n)),
        (w1, fsd(n, (lmbdas[0], 1))),
        (w2, fsd(n, (lmbdas[1], 1))),
        (w11, fsd(n, (lmbdas[0], 2))),
        (w, pm(n * [delta])),
    ]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return Enr2Scheme("Cools-Haegemans 3", n, weights, points, 7, _source)


# ERR There is a mistake here somewhere in the weights, but it's unclear where.
# TODO fix this
# def cools_haegemans_4(n, delta2=frac(2, 3)):
#     assert n >= 2
#     m = 4
#
#     lmbdas2 = _gener(delta2, 3, _mu)
#
#     delta4 = delta2 ** 2
#     delta6 = delta2 ** 3
#     delta8 = delta2 ** 4
#
#     w2 = -frac(
#         (
#             delta8 * (48 * n - 160)
#             + delta6 * (192 * n - 480)
#             + delta4 * (-408 * n + 96)
#             + delta2 * (240 * n + 600)
#             + (-45 * n - 270)
#         )
#         * (2 * delta2 - 3) ** 3,
#         4608 * (4 * delta4 + 20 * delta2 - 15) * (4 * delta2 - 3) * delta8,
#     )
#     w11 = -frac(
#         (12 * delta4 * n - 40 * delta4 - 12 * delta2 * n + 3 * n + 18)
#         * (2 * delta2 - 1) ** 2,
#         1536 * delta8,
#     )
#     w3 = frac(
#         (4 * delta4 - 24 * delta2 + 15) ** 4 * (2 * delta2 - 3) ** 3,
#         4608
#         * (8 * delta6 - 324 * delta4 + 414 * delta2 - 135)
#         * (4 * delta4 + 20 * delta2 - 15)
#         * (22 * delta2 - 15)
#         * delta8,
#     )
#     w21 = frac(
#         (2 * delta2 - 1) ** 2 * (2 * delta2 - 3) ** 3, 3072 * (4 * delta2 - 3) * delta8
#     )
#     w111 = frac((2 * delta2 - 1) ** 4, 1024 * delta8)
#
#     w = frac(_mu(2) ** m * _mu(0) ** (n - m), 2 ** n * delta2 ** m)
#
#     delta = sqrt(delta2)
#     lmbdas = [sqrt(lmbda2) for lmbda2 in lmbdas2]
#
#     data = [
#         (w2, fsd(n, (lmbdas[1], 1))),
#         (w11, fsd(n, (lmbdas[0], 2))),
#         (w3, fsd(n, (lmbdas[2], 1))),
#         (w21, fsd(n, (lmbdas[1], 1), (lmbdas[0], 1))),
#         (w111, fsd(n, (lmbdas[0], 3))),
#         (w, pm(n, delta)),
#     ]
#
#     # This is an attempt to find the correct symmetries, but to no avail. Something
#     # appears to be wrong with the weights.
#     print(n)
#     maxk = 30
#     k = 2 ** n
#     print("try")
#     for k2 in range(0, maxk, 2):
#         for k11 in range(0, maxk, 2):
#             for k3 in range(0, maxk, 2):
#                 for k21 in range(0, maxk, 2):
#                     for k111 in range(0, maxk, 2):
#                         val = (
#                             k2 * w2
#                             + k11 * w11
#                             + k3 * w3
#                             + k21 * w21
#                             + k111 * w111
#                             + k * w
#                         )
#                         if val == 1:
#                             print("SUCCESS", k2, k11, k3, k21, k111, k)
#     exit(1)
#
#     points, weights = untangle(data)
#
#     print(points)
#     print(weights)
#     print(sum(weights))
#     exit(1)
#     return Enr2Scheme("Cools-Haegemans 4", n, weights, points, 9, _source)
