import math

from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article
from ._helpers import C2Scheme, expand_symmetries, register

source = article(
    authors=["A.M. Cohen", "D.A. Gismalla"],
    title="Some integration formulae for symmetric functions of two variables",
    journal="International Journal of Computer Mathematics",
    year="1986",
    volume="19",
    number="1",
    pages="57-68",
    url="https://doi.org/10.1080/00207168608803504",
)


def cohen_gismalla_1():
    B = frac(5, 7)
    u, v = [sqrt((frac(1, 3) + i * sqrt(frac(2, 63))) / B) for i in [+1, -1]]
    d = {"zero": [[frac(8, 7)]], "c2": [[B, B], [u, v], [-v, u]]}
    points, weights = expand_symmetries(d)
    weights /= 4
    # This scheme is of order 5 for symmetric integrands
    return C2Scheme("Cohen-Gismalla 1", weights, points, 3, source, 4.996e-16)


def cohen_gismalla_2():
    # TODO improve precision
    X = 0.3850907
    Y = (24 / 35 - 28 / 45 * X) / (28 / 45 - 2 / 3 * X)
    alpha = (2 / 3 * Y - 28 / 45) / (Y - X)
    beta = (28 / 45 - 2 / 3 * X) / (Y - X)

    B = alpha / X
    C = beta / Y
    A = 4 * (1 - B - C)

    g1 = alpha / B / 2
    g2 = X * (Y / 9 - 2 / 15) / (Y - X) / alpha
    u = math.sqrt(g1 + math.sqrt(g1 ** 2 - g2))
    v = math.sqrt(g1 - math.sqrt(g1 ** 2 - g2))

    h1 = beta / C / 2
    h2 = Y * (2 / 15 - X / 9) / (Y - X) / beta
    r = math.sqrt(h1 - math.sqrt(h1 ** 2 - h2))
    s = math.sqrt(h1 + math.sqrt(h1 ** 2 - h2))

    d = {"zero": [[A]], "c2": [[B, B, C, C], [u, v, r, r], [-v, u, -s, s]]}
    points, weights = expand_symmetries(d)
    weights /= 4
    # ERR this scheme only has order 1
    # According to the article, it has order 7 for symmetric integrands.
    # Something is fishy...
    return C2Scheme("Cohen-Gismalla 2", weights, points, 1, source, 4.441e-16)


register([cohen_gismalla_1, cohen_gismalla_2])
