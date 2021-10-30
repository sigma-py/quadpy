import math

import numpy as np
import sympy

from ..helpers import article, expand_symmetries, prod
from ._helpers import CnScheme

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
    # 1/2 int_{-1}^1 x^j
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    return frac(1, j + 1) if j % 2 == 0 else 0


def cools_haegemans_1(n, delta2=1, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt
    assert frac(1, 3) <= delta2
    m = 1

    w0 = frac(3 * delta2 - 1, 3 * delta2)
    w = frac(1, 3) ** m / (2 ** n * delta2 ** m)

    d = {"0": [[w0]], "a": [[w], [sqrt(delta2)]]}
    points, weights = expand_symmetries(d, n)
    return CnScheme(f"Cools-Haegemans 1 (dim={n})", n, weights, points, 3, _source)


def cools_haegemans_2(n, delta2=1, symbolic=False):
    frac = sympy.Rational if symbolic else lambda a, b: a / b
    sqrt = sympy.sqrt if symbolic else math.sqrt
    # assert frac(5, 11) <= delta2 <= frac(5 * n + 4, 3 * (5 * n - 4))
    m = 2

    lmbdas2 = _gener(delta2, 1, _mu, symbolic)
    w0 = frac(
        -(15 * delta2 * n - 12 * delta2 - 5 * n - 4) * (3 * delta2 - 1),
        36 * delta2 ** 2,
    )
    w1 = frac(5 * (3 * delta2 - 1) ** 2, 72 * delta2)
    w = frac(1, 3) ** m / (2 ** n * delta2 ** m)

    lmbdas = [sqrt(lmbda2) for lmbda2 in lmbdas2]

    d = {"0": [[w0]], "a0": [[w1], [lmbdas[0]]], "a": [[w], [sqrt(delta2)]]}
    points, weights = expand_symmetries(d, n)
    return CnScheme("Cools-Haegemans 2", n, weights, points, 5, _source, 6.312e-14)


# ERR There is a mistake here somewhere in the weights, but it's unclear where.
# TODO fix this
# def cools_haegemans_3(n, delta2=1):
#     assert n >= 2
#     m = 3
#
#     lmbdas2 = _gener(delta2, 2, _mu)
#
#     delta4 = delta2 ** 2
#     delta6 = delta2 ** 3
#     w1 = frac(
#         -5
#         * (3 * delta2 - 1)
#         * (
#             135 * delta4 * n
#             - 68 * delta4
#             - 90 * delta2 * n
#             - 120 * delta2
#             + 15 * n
#             + 60
#         ),
#         2592 * delta6,
#     )
#     w2 = frac(
#         245 * (5 * delta2 - 3) ** 3 * (3 * delta2 - 1),
#         5184 * (31 * delta2 - 15) * delta6,
#     )
#     w11 = frac(25 * (3 * delta2 - 1) ** 3, 1728 * delta6)
#     w0 = frac(
#         -27 * (2 * n * (w1 + w2 - w11) + 2 * n ** 2 * w11 - 1) * delta6 + 1,
#         27 * delta6,
#     )
#     w = frac(_mu(2) ** m * _mu(0) ** (n - m), 2 ** n * delta2 ** m)
#
#     delta = sqrt(delta2)
#     lmbdas = [sqrt(lmbda2) for lmbda2 in lmbdas2]
#
#     # It's not entirely clear which points are meant to be part of the scheme. The
#     # weights add up to more than 1.
#     data = [
#         (w0, z(n)),
#         (w1, fsd(n, (lmbdas[0], 1))),
#         (w2, fsd(n, (lmbdas[1], 1))),
#         (w11, fsd(n, (lmbdas[0], 2))),
#         (w, pm(n, delta)),
#     ]
#
#     # This is an attempt to find the correct symmetries, but to no avail. Something
#     # appears to be wrong with the weights.
#     # print(n)
#     # maxk = 20
#     # k0 = 1
#     # k = 2 ** n
#     # print("try")
#     # for k1 in range(0, maxk, 2):
#     #     for k2 in range(0, maxk, 2):
#     #         for k11 in range(0, maxk, 2):
#     #             val = k0 * w0 + k1 * w1 + k2 * w2 + k11 * w11 + k * w
#     #             if val == 1:
#     #                 print("SUCCESS", k0, k1, k2, k11, k)
#     # exit(1)
#
#     points, weights = untangle(data)
#
#     print(points)
#     print(weights)
#     print(sum(weights))
#     exit(1)
#     return CnScheme("Cools-Haegemans 3", n, weights, points, 7, _source)


def _gener(delta2, m, mu, symbolic):
    # Computes the lambda_q from the article, eq. (9).
    lmbdas2 = []
    for q in range(1, m + 1):
        if not lmbdas2:
            # https://github.com/numpy/numpy/issues/16152
            coeffs = [1]
        else:
            coeffs = np.poly(lmbdas2)

        a0 = [c * mu(2 * (q - k) + 2, symbolic) for k, c in enumerate(coeffs)]
        a1 = [c * mu(2 * (q - k), symbolic) for k, c in enumerate(coeffs)]
        prod_ = prod([1 - lmbda2 / delta2 for lmbda2 in lmbdas2])
        a = sum(a0) - mu(2, symbolic) ** (q + 1) / mu(0, symbolic) ** q * prod_
        b = sum(a1) - mu(2, symbolic) ** (q + 1) / mu(0, symbolic) ** q * prod_ / delta2
        lmbdas2.append(a / b)
    return lmbdas2
