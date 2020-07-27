from ..helpers import article
from ._helpers import C2Scheme, concat, pm, zero

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
    # TODO improve precision
    u = 0.84623312
    v = 0.46607171
    weights, points = concat(zero(8.0 / 7.0), pm([5.0 / 7.0, u, -v], [5.0 / 7.0, v, u]))
    weights /= 4
    # This scheme is of order 5 for symmetric integrands
    return C2Scheme("Cohen-Gismalla 1", weights, points, 3, source, 2.603e-09)


def cohen_gismalla_2():
    # TODO improve precision
    r = 0.5878606
    s = 0.9353943
    u = 0.6105540
    v = 0.1109710
    A = 0.1856914
    B = 0.5951448
    C = 0.3584324
    weights, points = concat(zero(A), pm([B, u, -v], [B, v, u], [C, r, -s], [C, r, s]))
    weights /= 4
    # ERR this scheme only has order 1
    # According to the article, it has order 7 for symmetric integrands.
    # Something is fishy...
    return C2Scheme("Cohen-Gismalla 2", weights, points, 1, source, 1.001e-07)
