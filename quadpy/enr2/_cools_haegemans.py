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
