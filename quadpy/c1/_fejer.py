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


def fejer_1(n):
    degree = n

    points = -np.cos(np.pi * (np.arange(n) + 0.5) / n)

    # n -= 1
    N = np.arange(1, n, 2)
    length = len(N)
    m = n - length
    K = np.arange(m)

    v0 = np.concatenate(
        [
            2 * np.exp(1j * np.pi * K / n) / (1 - 4 * K ** 2),
            np.zeros(length + 1),
        ]
    )
    v1 = v0[:-1] + np.conjugate(v0[:0:-1])

    w = np.fft.ifft(v1)
    assert max(w.imag) < 1.0e-15
    weights = w.real

    return C1Scheme("Fejér 1", degree, weights, points, source)


def fejer_2(n):
    degree = n

    points = -np.cos((np.pi * np.arange(1, n + 1)) / (n + 1))

    n += 1
    N = np.arange(1, n, 2)
    length = len(N)
    m = n - length
    v0 = np.concatenate([2.0 / N / (N - 2), np.array([1.0 / N[-1]]), np.zeros(m)])
    v2 = -v0[:-1] - v0[:0:-1]

    w = np.fft.ihfft(v2)
    assert max(w.imag) < 1.0e-15
    w = w.real

    if n % 2 == 1:
        weights = np.concatenate([w, w[::-1]])
    else:
        weights = np.concatenate([w, w[len(w) - 2 :: -1]])

    # cut off first and last
    weights = weights[1:-1]
    return C1Scheme("Fejér 2", degree, weights, points, source)
