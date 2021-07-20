from math import factorial as fact

import ndim
from sympy import Rational as frac
from sympy import sqrt

from ..helpers import article, fsd, untangle
from ._helpers import UnScheme

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


def dobrodeev_1978(n: int):
    # from sympy import sqrt, factorial as fact, Rational as fr

    assert 2 <= n <= 20

    degree = 5

    dim_config = {
        2: ("I", None, 1, 1),
        3: ("I", None, 1, 1),
        4: (None, 2, None, None),
        5: ("I", 5, 1, 1),
        6: ("I", 5, 1, 1),
        7: (None, 3, None, None),
        8: ("I", 8, 1, 1),
        9: ("I", 9, 1, 1),
        10: ("I", 10, 1, 1),
        11: ("I", 11, 1, 1),
        12: ("I", 5, 1, 1),
        13: ("I", 13, 1, 2),
        14: ("I", 14, 1, 2),
        15: ("I", 15, 1, 2),
        16: ("I", 16, 1, 2),
        17: ("I", 17, 1, 2),
        18: ("I", 18, 1, 3),
        19: ("I", 19, 1, 3),
        20: ("I", 20, 1, 3),
    }

    pm_type, i, j, k = dim_config[n]

    I0 = ndim.nsphere.integrate_monomial(n * [0], symbolic=True)

    if i is None:
        G, b, c = _generate_jk(n, pm_type, j, k)
        data = [(G, fsd(n, (b, j), (c, k)))]
    elif j is None:
        assert k is None
        assert pm_type is None
        G, a = _generate_i(n, i)
        data = [(G, fsd(n, (a, i)))]
    else:
        I2 = ndim.nsphere.integrate_monomial([2] + (n - 1) * [0], symbolic=True)
        I22 = ndim.nsphere.integrate_monomial([2, 2] + (n - 2) * [0], symbolic=True)
        I4 = ndim.nsphere.integrate_monomial([4] + (n - 1) * [0], symbolic=True)
        G, a, b, c = _compute_dobrodeev(n, I0, I2, I22, I4, pm_type, i, j, k)

        data = [(G, fsd(n, (a, i))), (G, fsd(n, (b, j), (c, k)))]

    points, weights = untangle(data)
    return UnScheme("Dobrodeev 1978", n, weights, points, degree, source)


def _generate_i(n: int, i: int):
    L = fact(n) // fact(i) // fact(n - i) * 2 ** i
    G = frac(1, L)
    a = sqrt(frac(3, n + 2))
    return G, a


def _generate_jk(n, pm_type, j, k):
    M = fact(n) // fact(j) // fact(k) // fact(n - j - k) * 2 ** (j + k)
    G = frac(1, M)

    t = 1 if pm_type == "I" else -1
    b = sqrt(frac(1, j + k) * (1 + t * (k / j * sqrt(3 * (j + k) / (n + 2) - 1))))
    c = sqrt(frac(1, j + k) * (1 - t * (j / k * sqrt(3 * (j + k) / (n + 2) - 1))))
    return G, b, c


def _compute_dobrodeev(n, I0, I2, I22, I4, pm_type, i, j, k):
    """Same as the helper function in ..helpers, making use of the fact that
    `F == 0` for the sphere
    """
    # TODO prove F==0 analytically
    t = 1 if pm_type == "I" else -1

    L = fact(n) // (fact(i) * fact(n - i)) * 2 ** i
    M = fact(n) // (fact(j) * fact(k) * fact(n - j - k)) * 2 ** (j + k)
    N = L + M
    R = (
        -frac(j + k - i, i) * I2 ** 2 / I0 ** 2
        + frac(j + k - 1, n) * I4 / I0
        - frac(n - 1, n) * I22 / I0
    )
    H = frac(1, i) * (
        +(j + k - i) * I2 ** 2 / I0 ** 2
        + frac(j + k, n) * ((i - 1) * I4 / I0 - (n - 1) * I22 / I0)
    )
    Q = L / M * R + H

    G = frac(1, N)
    a = sqrt(frac(n, i) * I2 / I0)
    b = sqrt(frac(n, j + k) * (I2 / I0 + t * sqrt(k / j * Q)))
    c = sqrt(frac(n, j + k) * (I2 / I0 - t * sqrt(j / k * Q)))
    return G, a, b, c
