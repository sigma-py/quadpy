import numpy as np

from ..helpers import article, compute_dobrodeev, fsd, untangle
from ._helpers import CnScheme, integrate_monomial_over_ncube

_source = article(
    authors=["L.N. Dobrodeev"],
    title="Cubature rules with equal coefficients for integrating functions with respect to symmetric domains",
    journal="USSR Computational Mathematics and Mathematical Physics",
    volume="18",
    number="4",
    year="1978",
    pages="27-34",
    url="https://doi.org/10.1016/0041-5553%2878%2990064-2",
)


def dobrodeev_1978(n):
    assert 2 <= n <= 20

    dim_config = {
        2: ("II", 1, 1, 1),
        3: ("I", 3, 1, 1),
        4: ("I", 4, 1, 1),
        5: ("I", 5, 1, 1),
        6: ("I", 5, 1, 1),
        7: ("I", 6, 1, 2),
        8: ("I", 7, 1, 2),
        9: ("I", 6, 1, 2),
        10: ("I", 6, 1, 9),
        11: ("I", 6, 1, 10),
        12: ("I", 7, 1, 11),
        13: ("I", 8, 2, 2),
        14: ("I", 8, 1, 13),
        15: ("I", 9, 1, 14),
        16: ("I", 10, 2, 3),
        17: ("I", 10, 2, 3),
        18: ("I", 15, 3, 3),
        19: ("II", 9, 1, 18),
        20: ("I", 9, 1, 18),
    }

    ncube_limits = [[-1.0, 1.0]] * n
    I0 = integrate_monomial_over_ncube(ncube_limits, n * [0])
    I2 = integrate_monomial_over_ncube(ncube_limits, [2] + (n - 1) * [0])
    I22 = integrate_monomial_over_ncube(ncube_limits, [2, 2] + (n - 2) * [0])
    I4 = integrate_monomial_over_ncube(ncube_limits, [4] + (n - 1) * [0])

    pm_type, i, j, k = dim_config[n]

    G, a, b, c = compute_dobrodeev(n, I0, I2, I22, I4, pm_type, i, j, k, symbolic=True)

    data = [(G, fsd(n, (a, i))), (G, fsd(n, (b, j), (c, k)))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return CnScheme("Dobrodeev 1978", n, weights, points, 5, _source, 8.350e-14)
