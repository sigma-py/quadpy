# -*- coding: utf-8 -*-
#
import numpy

# [1] Remarks on the Disposition of Points in Numerical Integration
#     Formulas,
#     A. H. Stroud,
#     Mathematical Tables and Other Aids to Computation,
#     Vol. 11, No. 60 (Oct., 1957), pp. 257-261,
#     <https://dx.doi.org/10.2307/2001945>.


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    def __init__(self, n, index):
        self.name = 'Stroud({})'.format(index)
        reference_volume = 2.0**n
        if index == 'Cn 1-1':
            # centroid formula
            self.degree = 1
            self.weights = numpy.array([reference_volume])
            self.points = numpy.array([
                numpy.full(n, 0.0)
                ])
        elif index == 'Cn 1-2':
            # product trapezoidal formula
            self.degree = 1
            self.weights = numpy.full(2**n, 1.0)
            self.points = _pm(n, 1.0)
        elif index == 'Cn 2-1':
            # [1]
            self.degree = 2
            self.weights = numpy.full(n+1, reference_volume / (n+1))
            i = numpy.arange(n+1)
            n2 = n / 2 if n % 2 == 0 else (n-1)/2
            pts = [[
                numpy.sqrt(2.0/3.0) * numpy.cos(2*i*k*numpy.pi / (n+1)),
                numpy.sqrt(2.0/3.0) * numpy.sin(2*i*k*numpy.pi / (n+1))
                ] for k in range(1, n2+1)][0]
            if n % 2 == 1:
                sqrt3pm = numpy.full(n+1, 1.0 / numpy.sqrt(3.0))
                sqrt3pm[1::2] *= -1
                pts.append(sqrt3pm)

            self.points = numpy.vstack(pts).T
        elif index == 'Cn 2-2':
            # Thacher,
            # An efficient composite formula for multidimensional quadrature,
            # Comm. A.C.M. 7, 1964, 23-25.
            self.degree = 2
            r = numpy.sqrt(3.0) / 6.0
            self.weights = numpy.concatenate([
                numpy.full(1, reference_volume),
                numpy.full(n, r*reference_volume),
                numpy.full(n, -r*reference_volume),
                ])
            self.points = numpy.concatenate([
                numpy.array([numpy.full(n, 2*r)]),
                _s(n, -1.0, r),
                _s(n, +1.0, r),
                ])
        elif index == 'Cn 3-1':
            # [1]
            self.degree = 3
            self.weights = numpy.full(2*n, reference_volume / (2*n))
            i = numpy.arange(1, 2*n+1)
            n2 = n / 2 if n % 2 == 0 else (n-1)/2
            pts = [[
                numpy.sqrt(2.0/3.0) * numpy.cos((2*k-1)*i*numpy.pi / n),
                numpy.sqrt(2.0/3.0) * numpy.sin((2*k-1)*i*numpy.pi / n),
                ] for k in range(1, n2+1)][0]
            if n % 2 == 1:
                sqrt3pm = numpy.full(2*n, 1.0 / numpy.sqrt(3.0))
                sqrt3pm[1::2] *= -1
                pts.append(sqrt3pm)

            self.points = numpy.vstack(pts).T
        elif index == 'Cn 3-2':
            self.degree = 3
            self.weights = numpy.full(2*n, reference_volume / (2*n))
            r = numpy.sqrt(n / 3.0)
            self.points = _fs(n, r)
        elif index == 'Cn 3-3':
            # G.W. Tyler,
            # Numerical integration of functions of several variables,
            # Canad. J. Math. 5(1953), 393-412,
            # <https://dx.doi.org/10.4153/CJM-1953-044-1>.
            self.degree = 3
            r = numpy.sqrt(n / 3.0)
            self.weights = numpy.concatenate([
                numpy.full(1, (3.0 - n)/3.0 * reference_volume),
                numpy.full(2*n, reference_volume/6.0),
                ])
            self.points = numpy.concatenate([
                _z(n),
                _fs(n, 1.0)
                ])
        elif index == 'Cn 3-4':
            # product Gauss formula
            self.degree = 3
            self.weights = numpy.concatenate([
                numpy.full(2**n, reference_volume / 2**n),
                ])
            r = numpy.sqrt(3.0) / 3.0
            self.points = _pm(n, r)
        elif index == 'Cn 3-5':
            # G.M. Ewing,
            # On Approximate Cubature,
            # The American Mathematical Monthly,
            # Vol. 48, No. 2 (Feb., 1941), pp. 134-136,
            # DOI: 10.2307/2303604.
            self.degree = 3
            self.weights = numpy.concatenate([
                numpy.full(1, 2.0/3.0 * reference_volume),
                numpy.full(2**n, 1.0/3.0 / 2**n * reference_volume),
                ])
            r = numpy.sqrt(3.0) / 3.0
            self.points = numpy.concatenate([
                _z(n),
                _pm(n, 1.0),
                ])
        elif index == 'Cn 3-6':
            # product Simpson's formula
            from itertools import product
            self.degree = 3
            pts = product([-1, 0, 1], repeat=n)
            A = {-1: 1.0/3.0, 0: 4.0/3.0, 1: 1.0/3.0}
            self.weights = numpy.array([
                numpy.product([A[p] for p in pt])
                for pt in pts
                ])
            pts = product([-1, 0, 1], repeat=n)
            self.points = numpy.array(list(pts))
        else:
            assert False

        return


def _z(n):
    return numpy.zeros((1, n))


def _fs(n, r):
    if n == 2:
        return numpy.array([
            [+r, 0.0],
            [0.0, +r],
            [-r, 0.0],
            [0.0, -r],
            ])
    assert n == 3
    return numpy.array([
        [+r, 0.0, 0.0],
        [0.0, +r, 0.0],
        [0.0, 0.0, +r],
        [-r, 0.0, 0.0],
        [0.0, -r, 0.0],
        [0.0, 0.0, -r],
        ])


def _s(n, a, b):
    if n == 2:
        return numpy.array([
            [a, b],
            [b, a],
            ])
    assert n == 3
    return numpy.array([
        [a, b, b],
        [b, a, b],
        [b, b, a],
        ])


def _pm(n, a):
    if n == 2:
        return numpy.array([
            [+a, +a],
            [+a, -a],
            [-a, +a],
            [-a, -a],
            ])
    assert n == 3
    return numpy.array([
        [+a, +a, +a],
        [+a, +a, -a],
        [+a, -a, +a],
        [+a, -a, -a],
        [-a, +a, +a],
        [-a, +a, -a],
        [-a, -a, +a],
        [-a, -a, -a],
        ])
