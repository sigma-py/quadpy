# -*- coding: utf-8 -*-
#
import itertools
import numpy

# [1] Remarks on the Disposition of Points in Numerical Integration
#     Formulas,
#     A. H. Stroud,
#     Mathematical Tables and Other Aids to Computation,
#     Vol. 11, No. 60 (Oct., 1957), pp. 257-261,
#     <https://dx.doi.org/10.2307/2001945>.
# Stroud [6] A.H. Stroud,
#            Some Fifth Degree Integration Formulas for Symmetric Regions,
#            Mathematics of Computation,
#            Vol. 20, No. 93 (Jan., 1966), pp. 90-97,
#            Published by: American Mathematical Society,
#            DOI: 10.2307/2004272.


class Stroud(object):
    '''
    Arthur Stroud,
    Approximate Calculation of Multiple Integrals,
    Prentice Hall, 1971.
    '''
    # pylint: disable=too-many-locals
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
            n2 = n // 2 if n % 2 == 0 else (n-1)//2
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
            n2 = n // 2 if n % 2 == 0 else (n-1)//2
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
            self.points = _fs1(n, r)
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
                _fs1(n, 1.0)
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
            self.degree = 3
            lst = n * [[1.0/3.0, 4.0/3.0, 1.0/3.0]]
            self.weights = numpy.product(
                numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n),
                axis=-1
                )
            lst = n * [[-1.0, 0.0, 1.0]]
            self.points = numpy.array(numpy.meshgrid(*lst)).T.reshape(-1, n)
        # elif index == 'Cn 5-1':
        # Cn 5-1 is not implemented because it's based on explicit values only
        # given for n=4,5,6.
        elif index == 'Cn 5-2':
            # Preston C. Hammer and Arthur H. Stroud,
            # Numerical Evaluation of Multiple Integrals II.
            self.degree = 5
            r = numpy.sqrt(3.0 / 5.0)
            self.points = numpy.concatenate([
                _z(n),
                _fs1(n, r),
                _fs2(n, r),
                ])
            self.weights = numpy.concatenate([
                numpy.full(
                    1, (25*n**2 - 115*n + 162)/162.0 * reference_volume
                    ),
                numpy.full(2*n, (70 - 25*n)/162.0 * reference_volume),
                numpy.full(2*n*(n-1), 25.0/324.0 * reference_volume),
                ])
        elif index == 'Cn 5-3':
            # A. H. Stroud,
            # Extensions of Symmetric Integration Formulas,
            # Mathematics of Computation,
            # Vol. 22, No. 102 (Apr., 1968), pp. 271-274,
            # Published by: American Mathematical Society,
            # DOI: 10.2307/2004655.
            self.degree = 5
            r = numpy.sqrt(7.0 / 15.0)
            s = numpy.sqrt((7.0 + numpy.sqrt(24.0)) / 15.0)
            t = numpy.sqrt((7.0 - numpy.sqrt(24.0)) / 15.0)
            self.points = numpy.concatenate([
                _z(n),
                _s2(n, r),
                _s2(n, -r),
                _fs1(n, r),
                _s11(n, +s, -t),
                _s11(n, -s, +t),
                _fs1(n, s),
                _fs1(n, t),
                ])
            self.weights = numpy.concatenate([
                numpy.full(1, (5*n**2 - 15*n+14)/14.0 * reference_volume),
                numpy.full(n*(n-1), 25.0/168.0 * reference_volume),
                numpy.full(2*n, -25*(n-2)/168.0 * reference_volume),
                numpy.full(2*n*(n-1), 5.0/48.0 * reference_volume),
                numpy.full(4*n, -5*(n-2)/48.0 * reference_volume),
                ])
        elif index == 'Cn 5-4':
            # Stroud [6]
            self.degree = 5
            r = numpy.sqrt((5*n + 4) / 30.0)
            s = numpy.sqrt((5*n + 4.0) / (15*n - 12.0))
            self.points = numpy.concatenate([
                _fs1(n, r),
                _pm(n, s),
                ])
            self.weights = numpy.concatenate([
                numpy.full(2*n, 40.0 / (5*n+4)**2 * reference_volume),
                numpy.full(
                    2**n,
                    ((5*n - 4.0) / (5*n + 4))**2 / 2**n * reference_volume
                    ),
                ])
        elif index == 'Cn 5-5':
            # Mustard, Lyness, Blatt,
            # Numerical quadrature in N dimensions.
            self.degree = 5
            r = numpy.sqrt(2.0 / 5.0)
            self.points = numpy.concatenate([
                _z(n),
                _fs1(n, r),
                _pm(n, 1.0),
                ])
            self.weights = numpy.concatenate([
                numpy.full(1, (8 - 5*n)/9.0 * reference_volume),
                numpy.full(2*n, 5.0/18.0 * reference_volume),
                numpy.full(2**n, 1.0/9.0 / 2**n * reference_volume),
                ])
        elif index == 'Cn 5-6':
            # Stroud [6]
            self.degree = 5
            s = numpy.sqrt(1.0 / 3.0)

            pts = [_z(n)]
            wts = [[4.0 / (5*n + 4) * reference_volume]]
            for k in range(1, n+1):
                r = numpy.sqrt((5*k + 4) / 15.0)
                arr = numpy.zeros((2**(n-k+1), n))
                arr[:, k-1:] = _pm(n-k+1, 1.0)
                arr[:, k-1] *= r
                arr[:, k:] *= s
                pts.append(arr)
                num_pts = len(pts[-1])
                b = 5.0 * 2.0**(k-n+1) / (5.0*k-1.0) / (5.0*k+4.0) \
                    * reference_volume
                wts.append(num_pts * [b])

            self.points = numpy.vstack(pts)
            self.weights = numpy.concatenate(wts)
        elif index == 'Cn 5-7':
            # Stroud [6]
            self.degree = 5
            r = numpy.sqrt((5*n + 4 + 2*(n-1)*numpy.sqrt(5*n+4)) / (15.0*n))
            s = numpy.sqrt((5*n + 4 - 2*numpy.sqrt(5*n+4)) / (15.0*n))
            self.points = numpy.concatenate([
                _z(n),
                _fs11(n, r, s),
                ])
            self.weights = numpy.concatenate([
                numpy.full(1, 4.0/(5*n+4) * reference_volume),
                numpy.full(n * 2**n, 5.0/(5*n+4) / 2**n * reference_volume),
                ])
        elif index == 'Cn 5-8':
            # Stroud [6]
            assert n >= 3
            self.degree = 5
            r = numpy.sqrt(
                (5*n - 2*numpy.sqrt(5.0) + 2*(n-1)*numpy.sqrt(5*n+5))
                / (15.0*n)
                )
            # This sqrt() is imaginary for negative for n=2.
            s = numpy.sqrt(
                (5*n - 2*numpy.sqrt(5.0) - 2*numpy.sqrt(5*n+5)) / (15.0*n)
                )
            t = numpy.sqrt((5.0 + 2*numpy.sqrt(5)) / 15.0)
            self.points = numpy.concatenate([
                _fs11(n, r, s),
                _pm(n, t)
                ])
            self.weights = \
                numpy.full((n+1) * 2**n, reference_volume / 2**n / (n+1))
        else:
            assert False

        return


def _z(n):
    return numpy.zeros((1, n))


def _fs1(n, r):
    return _combine([[+r, -r]] + [[0.0]] * (n-1))


def _fs11(n, r, s):
    '''Get all permutations of [+-r, +-s, ..., +-s] of length n.
    len(out) == n * 2**n.
    '''
    # <https://stackoverflow.com/a/45321972/353337>
    k1 = 1
    k2 = n-1
    idx = itertools.combinations(range(k1 + k2), k1)
    vs = ((s if j not in i else r for j in range(k1 + k2)) for i in idx)
    return numpy.array(list(itertools.chain.from_iterable(
        itertools.product(*((+vij, -vij) for vij in vi)) for vi in vs
        )))


def _fs2(n, r):
    '''Get all permutations of [+-r, +-r, 0, ..., 0] of length n.
    len(out) == 2 * n * (n-1).
    '''
    return _combine([[+r, -r]] * 2 + [[0.0]] * (n-2))


def _combine(pools):
    '''Given an input array with lists of options, e.g.,

    [[a, b], [c], [d]],

    this methods returns all combinations with one element from each
    subset, e.g.,

    [a, c, d], [a, d, c], [c, d, a], ...
    [b, c, d], [b, d, c], [c, d, b], ...
    '''
    # https://stackoverflow.com/a/45322199/353337
    return numpy.array(list(set(itertools.chain.from_iterable([
        itertools.permutations(x) for x in itertools.product(*pools)
        ]))))


def _s(n, a, b):
    '''Get all permutations of [a, b, ..., b] of length n.
    len(out) == n.
    '''
    out = numpy.full((n, n), b)
    numpy.fill_diagonal(out, a)
    return out


def _s2(n, a):
    '''Get all permutations of [a, a, 0, 0, ..., 0] of length n.
    len(out) == (n-1)*n / 2.
    '''
    return _s11(n, a, a)


def _s11(n, a, b):
    '''Get all permutations of [a, b, 0, 0, ..., 0] of length n.
    len(out) == (n-1)*n.
    '''
    s = [a, b] + (n-2) * [0]
    # First permutations, then set can be really inefficient if items are
    # repeated. Check out <https://stackoverflow.com/q/6284396/353337> for
    # improvements.
    return numpy.array(list(set(itertools.permutations(s, n))))


def _pm(n, a):
    '''Return all combinations of [+a, -a] with length n (with repetition).
    len(out) == 2**n.
    '''
    return numpy.array(list(itertools.product([+a, -a], repeat=n)))
