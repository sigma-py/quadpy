import numpy as np

from ..helpers import article
from ._helpers import C1Scheme

source = article(
    authors=["J. Waldvogel"],
    title="Fast Construction of the Fejér and Clenshaw–Curtis Quadrature Rules",
    journal="BIT Numerical Mathematics",
    month="mar",
    year="2006",
    volume="46",
    number="1",
    pages="195–202",
    url="https://doi.org/10.1007/s10543-006-0045-4",
)


def clenshaw_curtis(n):
    degree = n

    points = -np.cos((np.pi * np.arange(n)) / (n - 1))

    if n == 2:
        weights = np.array([1.0, 1.0])
        return C1Scheme("Clenshaw-Curtis", degree, weights, points)

    n -= 1
    N = np.arange(1, n, 2)
    length = len(N)
    m = n - length
    v0 = np.concatenate([2.0 / N / (N - 2), np.array([1.0 / N[-1]]), np.zeros(m)])
    v2 = -v0[:-1] - v0[:0:-1]
    g0 = -np.ones(n)
    g0[length] += n
    g0[m] += n
    g = g0 / (n ** 2 - 1 + (n % 2))

    w = np.fft.ihfft(v2 + g)
    assert max(w.imag) < 1.0e-15
    w = w.real

    if n % 2 == 1:
        weights = np.concatenate([w, w[::-1]])
    else:
        weights = np.concatenate([w, w[len(w) - 2 :: -1]])

    return C1Scheme("Clenshaw-Curtis", degree, weights, points, source)
