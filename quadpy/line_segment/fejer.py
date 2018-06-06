# -*- coding: utf-8 -*-
#
import numpy


class Fejer1(object):
    """
    Fejér-type-1 quadrature.

    Weights are constructed after

    J. Waldvogel,
    Fast Construction of the Fejér and Clenshaw–Curtis Quadrature Rules,
    BIT Numerical Mathematics, March 2006, Volume 46, Issue 1, pp 195–202,
    DOI: 10.1007/s10543-006-0045-4,
    <https://doi.org/10.1007/s10543-006-0045-4>.
    """

    def __init__(self, n):
        self.degree = n

        self.points = -numpy.cos(numpy.pi * (numpy.arange(n) + 0.5) / n)

        # n -= 1
        N = numpy.arange(1, n, 2)
        length = len(N)
        m = n - length
        K = numpy.arange(m)

        v0 = numpy.concatenate(
            [
                2 * numpy.exp(1j * numpy.pi * K / n) / (1 - 4 * K ** 2),
                numpy.zeros(length + 1),
            ]
        )
        v1 = v0[:-1] + numpy.conjugate(v0[:0:-1])

        w = numpy.fft.ifft(v1)
        assert max(w.imag) < 1.0e-15
        self.weights = w.real

        return


class Fejer2(object):
    """
    Fejér-type-2 quadrature.

    Weights are constructed after

    J. Waldvogel,
    Fast Construction of the Fejér and Clenshaw–Curtis Quadrature Rules,
    BIT Numerical Mathematics, March 2006, Volume 46, Issue 1, pp 195–202,
    DOI: 10.1007/s10543-006-0045-4,
    <https://doi.org/10.1007/s10543-006-0045-4>.
    """

    def __init__(self, n):
        self.degree = n

        self.points = -numpy.cos((numpy.pi * numpy.arange(1, n + 1)) / (n + 1))

        n += 1
        N = numpy.arange(1, n, 2)
        length = len(N)
        m = n - length
        v0 = numpy.concatenate(
            [2.0 / N / (N - 2), numpy.array([1.0 / N[-1]]), numpy.zeros(m)]
        )
        v2 = -v0[:-1] - v0[:0:-1]

        w = numpy.fft.ihfft(v2)
        assert max(w.imag) < 1.0e-15
        w = w.real

        if n % 2 == 1:
            self.weights = numpy.concatenate([w, w[::-1]])
        else:
            self.weights = numpy.concatenate([w, w[len(w) - 2 :: -1]])

        # cut off first and last
        self.weights = self.weights[1:-1]
        return
