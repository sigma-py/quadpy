from sympy import sqrt, Rational as frac
import numpy

from ..helpers import article, fsd, untangle, z, pm
from ._helpers import NCubeScheme

_citation = article(
    authors=["Ronald Cools", "Ann Haegemans"],
    title="An imbedded family of cubature formulae for n-dimensional product regions",
    journal="Journal of Computational and Applied Mathematics",
    volume="51",
    year="1994",
    pages="251-262",
    url="https://doi.org/10.1016/0377-0427(92)00007-V",
)


def _mu(j):
    # 1/2 int_{-1}^1 x^j
    return frac(1, j + 1) if j % 2 == 0 else 0


def cools_haegemans_1(n, delta2=1):
    m = 1

    w0 = frac(3 * delta2 - 1, 3 * delta2)
    w = frac(_mu(2) ** m * _mu(0) ** (n - m), 2 ** n * delta2 ** m)

    data = [
        (w0, z(n)),
        (w, pm(n, sqrt(delta2))),
    ]

    points, weights = untangle(data)
    weights *= 2 ** n
    return NCubeScheme("Cools-Haegemans 1", n, weights, points, 3, _citation)


def cools_haegemans_2(n, delta2=1):
    assert n >= 1
    assert (
        frac(5, 11) <= delta2 <= 1
    ), "delta ** 2 must be in [5/11, 1] for all lmbda to be inside the domain"
    m = 2

    lmbdas2 = _gener(delta2, 1, _mu)
    w0 = frac(
        -(15 * delta2 * n - 12 * delta2 - 5 * n - 4) * (3 * delta2 - 1),
        36 * delta2 ** 2,
    )
    w1 = frac(5 * (3 * delta2 - 1) ** 2, 72 * delta2)
    w = frac(_mu(2) ** m * _mu(0) ** (n - m), 2 ** n * delta2 ** m)

    lmbdas = [sqrt(lmbda2) for lmbda2 in lmbdas2]

    data = [
        (w0, z(n)),
        (w1, fsd(n, (lmbdas[0], 1))),
        (w, pm(n, sqrt(delta2))),
    ]

    points, weights = untangle(data)
    weights *= 2 ** n
    return NCubeScheme("Cools-Haegemans 2", n, weights, points, 5, _citation)


def cools_haegemans_3(n, delta2=1):
    # some restrictions on delta2 apply here
    assert n >= 2
    m = 3

    lmbdas2 = _gener(delta2, 2, _mu)
    w1 = frac(
        -5
        * (3 * delta2 - 1)
        * (
            135 * delta2 ** 2 * n
            - 68 * delta2 ** 2
            - 90 * delta2 * n
            - 120 * delta2
            + 15 * n
            + 60
        ),
        2592 * delta2 ** 3,
    )
    w2 = frac(
        245 * (5 * delta2 - 3) ** 3 * (3 * delta2 - 1),
        5184 * (31 * delta2 - 15) * delta2 ** 3,
    )
    w11 = frac(25 * (3 * delta2 - 1) ** 3, 1728 * delta2 ** 3)
    w0 = frac(
        -27 * (2 * n * (w1 + w2 - w11) + 2 * n ** 2 * w11 - 1) * delta2 ** 3 + 1,
        27 * delta2 ** 3,
    )
    w = frac(_mu(2) ** m * _mu(0) ** (n - m), 2 ** n * delta2 ** m)

    lmbdas = [sqrt(lmbda2) for lmbda2 in lmbdas2]

    data = [
        (w0, z(n)),
        (w1, fsd(n, (lmbdas[0], 1))),
        (w1, fsd(n, (lmbdas[1], 1))),
        (w2, fsd(n, (lmbdas[0], 2))),
        (w2, fsd(n, (lmbdas[1], 2))),
        (w11, fsd(n, (lmbdas[0], 1), (lmbdas[1], 1))),
        (w, pm(n, sqrt(delta2))),
    ]

    points, weights = untangle(data)

    print(points)
    print(weights)
    print(sum(weights))
    exit(1)
    weights *= 2 ** n
    return NCubeScheme("Cools-Haegemans 3", n, weights, points, 7, _citation)


def _prod(iterable):
    from functools import reduce
    import operator

    # https://stackoverflow.com/a/7948307/353337
    return reduce(operator.mul, iterable, 1)


def _gener(delta2, m, mu):
    # Computes the lambda_q from the article, eq. (9).
    lmbdas2 = []
    for q in range(1, m + 1):
        if not lmbdas2:
            # https://github.com/numpy/numpy/issues/16152
            coeffs = [1]
        else:
            coeffs = numpy.poly(lmbdas2)

        a0 = [c * mu(2 * (q - k) + 2) for k, c in enumerate(coeffs)]
        a1 = [c * mu(2 * (q - k)) for k, c in enumerate(coeffs)]
        prod = _prod([1 - lmbda2 / delta2 for lmbda2 in lmbdas2])
        a = sum(a0) - mu(2) ** (q + 1) / mu(0) ** q * prod
        b = sum(a1) - mu(2) ** (q + 1) / mu(0) ** q * prod / delta2
        lmbdas2.append(a / b)
    return lmbdas2
