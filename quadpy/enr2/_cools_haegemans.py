from sympy import sqrt, Rational as frac, pi

from ..helpers import article, fsd, untangle, z, pm
from ._helpers import Enr2Scheme
from ..ncube._cools_haegemans import _gener

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
    # 1/sqrt(pi) ** n int int ... int exp(-r ** 2) dr
    if j == 0:
        return 1
    elif j == 1:
        return 0
    return frac(j - 1, 2) * _mu(j - 2)


def cools_haegemans_1(n, delta2=1):
    assert frac(1, 2) <= delta2
    m = 1

    w0 = frac(2 * delta2 - 1, 2 * delta2)
    w = frac(_mu(2) ** m * _mu(0) ** (n - m), 2 ** n * delta2 ** m)

    data = [
        (w0, z(n)),
        (w, pm(n, sqrt(delta2))),
    ]

    points, weights = untangle(data)
    weights *= sqrt(pi) ** n
    return Enr2Scheme("Cools-Haegemans 1", n, weights, points, 3, _citation)


def cools_haegemans_2(n, delta2=1):
    assert n >= 1
    assert frac(1, 2) <= delta2
    if n > 2:
        assert delta2 <= frac(n + 2, 2 * n - 4)
    m = 2

    lmbdas2 = _gener(delta2, 1, _mu)

    w0 = frac(
        -(2 * delta2 - 1) * (2 * delta2 * n - 4 * delta2 - n - 2), 8 * delta2 ** 2
    )
    w1 = frac((2 * delta2 - 1) ** 2, 16 * delta2 ** 2)
    w = frac(_mu(2) ** m * _mu(0) ** (n - m), 2 ** n * delta2 ** m)

    lmbdas = [sqrt(lmbda2) for lmbda2 in lmbdas2]

    data = [
        (w0, z(n)),
        (w1, fsd(n, (lmbdas[0], 1))),
        (w, pm(n, sqrt(delta2))),
    ]

    points, weights = untangle(data)
    weights *= sqrt(pi) ** n
    return Enr2Scheme("Cools-Haegemans 2", n, weights, points, 5, _citation)


def cools_haegemans_3(n, delta2=frac(2, 3)):
    assert n >= 2
    m = 3

    lmbdas2 = _gener(delta2, 2, _mu)

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

    w = frac(_mu(2) ** m * _mu(0) ** (n - m), 2 ** n * delta2 ** m)

    delta = sqrt(delta2)
    lmbdas = [sqrt(lmbda2) for lmbda2 in lmbdas2]

    data = [
        (w0, z(n)),
        (w1, fsd(n, (lmbdas[0], 1))),
        (w2, fsd(n, (lmbdas[1], 1))),
        (w11, fsd(n, (lmbdas[0], 2))),
        (w, pm(n, delta)),
    ]

    points, weights = untangle(data)
    weights *= sqrt(pi) ** n
    return Enr2Scheme("Cools-Haegemans 3", n, weights, points, 7, _citation)
