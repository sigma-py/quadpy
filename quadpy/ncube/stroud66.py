# -*- coding: utf-8 -*-
#
import numpy

from .helpers import _fsd, _pm, _z, _fs11


class Stroud66(object):
    '''
    A.H. Stroud,
    Some Fifth Degree Integration Formulas for Symmetric Regions,
    Mathematics of Computation,
    Vol. 20, No. 93 (Jan., 1966), pp. 90-97,
    Published by: American Mathematical Society,
    <https://doi.org/10.1090/S0025-5718-1966-0191094-8>.
    '''
    def __init__(self, n, variant):
        self.name = 'Stroud57{}'.format(variant)
        reference_volume = 2.0**n
        self.degree = 5
        if variant == 'a':
            r = numpy.sqrt((5*n + 4) / 30.0)
            s = numpy.sqrt((5*n + 4.0) / (15*n - 12.0))
            self.points = numpy.concatenate([
                _fsd(n, r, 1),
                _pm(n, s),
                ])
            self.weights = numpy.concatenate([
                numpy.full(2*n, 40.0 / (5*n+4)**2 * reference_volume),
                numpy.full(
                    2**n,
                    ((5*n - 4.0) / (5*n + 4))**2 / 2**n * reference_volume
                    ),
                ])
        elif variant == 'b':
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
        elif variant == 'c':
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
        else:
            assert variant == 'd'
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

        return
