import ndim
import numpy as np

from ..helpers import article, compute_dobrodeev, fsd, untangle
from ._helpers import SnScheme

source = article(
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
        3: ("II", 1, 1, 1),
        4: ("II", 4, 1, 1),
        5: ("I", 3, 1, 1),
        6: ("I", 3, 1, 1),
        7: ("I", 7, 1, 1),
        8: ("I", 8, 1, 1),
        9: ("I", 9, 1, 1),
        10: ("I", 9, 1, 2),
        11: ("I", 5, 1, 2),
        12: ("I", 12, 1, 2),
        13: ("I", 13, 1, 2),
        14: ("I", 14, 1, 2),
        15: ("I", 15, 1, 2),
        16: ("II", 15, 1, 4),
        17: ("I", 17, 1, 3),
        18: ("I", 18, 1, 3),
        19: ("I", 19, 1, 3),
        20: ("I", 20, 1, 3),
    }

    I0 = ndim.nball.integrate_monomial(n * [0])
    I2 = ndim.nball.integrate_monomial([2] + (n - 1) * [0])
    I22 = ndim.nball.integrate_monomial([2, 2] + (n - 2) * [0])
    I4 = ndim.nball.integrate_monomial([4] + (n - 1) * [0])
    pm_type, i, j, k = dim_config[n]
    G, a, b, c = compute_dobrodeev(n, I0, I2, I22, I4, pm_type, i, j, k)

    data = [(G, fsd(n, (a, i))), (G, fsd(n, (b, j), (c, k)))]

    points, weights = untangle(data)
    points = np.ascontiguousarray(points.T)
    return SnScheme("Dobrodeev 1978", n, weights, points, 5, source)
